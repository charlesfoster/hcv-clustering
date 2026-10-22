#!/usr/bin/env python3
"""Persistent cluster labels for hcv-cluster runs.

Each run re-clusters from scratch, so run-local IDs such as ``1a_C0001`` change whenever
samples are added. A *registry* (a directory of plain CSV files) remembers which label
each sample belonged to, and every new run maps its connected components back onto
those labels by sample overlap:

* a component whose members carry one old label keeps it (the cluster grew);
* a component claiming several old labels merges them; the **oldest** label survives and
  the others are retired with ``merged_into`` recorded;
* an old label whose members now sit in several components stays with the component
  holding most of them; the other components get new labels with ``split_from`` recorded;
* a component of two or more samples with no old label is named
  ``<prefix>-<YYMM>-<NNN>``, where YYMM is the month the label was created and NNN counts
  up within that month;
* singletons never receive a label.

Labels are tracked per *track* (one per clustering algorithm, e.g. ``tn93``). Each track
records the definition it was created under (metric, region, threshold) and refuses to
be updated under a different one, because a label is only meaningful for the clustering
that produced it. Label numbers are shared across tracks, so a label never names two
different clusters. Any algorithm that can produce a sample -> component partition can
be tracked by calling :func:`assign_labels` with its own track name.

Registry files: ``tracks.csv``, ``labels.csv``, ``membership.csv``, ``events.csv`` and
``runs.csv``. ``events.csv`` is append-only history of every create/join/merge/split.
"""

from __future__ import annotations

import contextlib
import copy
import csv
import datetime
import fcntl
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Iterator, Mapping

DEFAULT_PREFIX = "HCV"
UNASSIGNED = "UNASSIGNED"
SINGLETON_NOTE = "singleton"

TRACKS_FILE = "tracks.csv"
LABELS_FILE = "labels.csv"
MEMBERSHIP_FILE = "membership.csv"
EVENTS_FILE = "events.csv"
RUNS_FILE = "runs.csv"
LOCK_FILE = ".registry.lock"

# (genotype, base sample id): stable across runs, unlike run-local sequence IDs
# (a mixed-infection variant may be v1 in one run and v2 in the next).
Member = tuple[str, str]


class RegistryError(ValueError):
    """The registry cannot be used safely for this run."""


@dataclass
class LabelRecord:
    label: str
    track: str
    genotype: str
    created_seq: int
    created_run: str
    created_date: str
    status: str = "active"  # active | merged
    merged_into: str = ""
    split_from: str = ""


@dataclass
class Registry:
    tracks: dict[str, str] = field(default_factory=dict)  # track -> definition JSON
    labels: dict[str, LabelRecord] = field(default_factory=dict)
    membership: dict[tuple[str, str, str], str] = field(default_factory=dict)  # (track, gt, sample) -> label
    events: list[dict[str, str]] = field(default_factory=list)
    runs: list[dict[str, str]] = field(default_factory=list)

    def resolve(self, label: str) -> str:
        """Follow ``merged_into`` pointers to the label currently in force."""
        seen = set()
        while label in self.labels and self.labels[label].status == "merged" and label not in seen:
            seen.add(label)
            label = self.labels[label].merged_into
        return label

    def check_definition(self, track: str, definition: Mapping[str, object]) -> None:
        """Pin a new track to ``definition``; reject a different one for a known track."""
        encoded = json.dumps(dict(definition), sort_keys=True)
        known = self.tracks.get(track)
        if known is None:
            self.tracks[track] = encoded
        elif known != encoded:
            raise RegistryError(
                f"Registry track '{track}' was created with {known}, but this run uses {encoded}. "
                "Cluster labels are only valid for the clustering definition that created them; "
                "use a different registry directory (or track) for a different definition."
            )

    def _next_seq(self) -> int:
        return max((record.created_seq for record in self.labels.values()), default=0) + 1

    def _mint(self, prefix: str, yymm: str) -> str:
        pattern = re.compile(rf"^{re.escape(prefix)}-{yymm}-(\d+)$")
        used = [int(m.group(1)) for label in self.labels if (m := pattern.match(label))]
        return f"{prefix}-{yymm}-{max(used, default=0) + 1:03d}"


@dataclass
class LabelResult:
    """Outcome of one labelling pass."""

    labels: dict[Member, str]  # "" for singletons
    events: list[dict[str, str]]
    new_labels: set[str]


def assign_labels(
    registry: Registry,
    track: str,
    definition: Mapping[str, object],
    partition: Mapping[Member, str],
    *,
    run_id: str,
    run_date: str,
    prefix: str = DEFAULT_PREFIX,
) -> LabelResult:
    """Label the components of ``partition`` (member -> run-local component ID).

    Mutates ``registry``; callers wanting a dry run pass a copy.
    """
    registry.check_definition(track, definition)
    parsed_date = datetime.date.fromisoformat(run_date)
    run_date = parsed_date.isoformat()  # fromisoformat also accepts 20261023 / 2026-W43-5
    yymm = parsed_date.strftime("%y%m")

    components: dict[str, list[Member]] = defaultdict(list)
    for member, component in partition.items():
        components[component].append(member)
    labelled = {c: sorted(ms) for c, ms in components.items() if len(ms) >= 2}

    # Snapshot of who held which label before this run changes anything, so provenance
    # (split_from, "joined" events) does not depend on the order components are processed.
    prior: dict[Member, str] = {
        member: registry.resolve(old)
        for member in partition
        if (old := registry.membership.get((track, *member)))
    }

    # Where do the members of each previously known label sit now?
    spread: dict[str, Counter[str]] = defaultdict(Counter)
    for member, old in prior.items():
        spread[old][partition[member]] += 1

    def home_rank(component: str, count: int) -> tuple[int, int, Member]:
        return (-count, -len(components[component]), min(components[component]))

    home = {
        label: min(counts.items(), key=lambda item: home_rank(item[0], item[1]))[0]
        for label, counts in spread.items()
    }
    claims: dict[str, list[str]] = defaultdict(list)
    for label, component in home.items():
        claims[component].append(label)

    events: list[dict[str, str]] = []
    new_labels: set[str] = set()

    def log(event: str, label: str, related: str = "", detail: str = "") -> None:
        events.append(
            {
                "run_id": run_id,
                "date": run_date,
                "track": track,
                "event": event,
                "label": label,
                "related_label": related,
                "detail": detail,
            }
        )

    assigned: dict[Member, str] = {}
    for component in sorted(labelled, key=lambda c: labelled[c][0]):
        members = labelled[component]
        genotype = members[0][0]
        claimed = sorted(claims.get(component, ()), key=lambda l: registry.labels[l].created_seq)
        if claimed:
            label = claimed[0]
            for other in claimed[1:]:
                registry.labels[other].status = "merged"
                registry.labels[other].merged_into = label
                for key, held in registry.membership.items():
                    if held == other:
                        registry.membership[key] = label
                log("merged", other, label, f"merged into older label {label}")
        else:
            # Members displaced from an old label whose home is elsewhere mark a split.
            origin = Counter(prior[m] for m in members if m in prior)
            split_from = origin.most_common(1)[0][0] if origin else ""
            label = registry._mint(prefix, yymm)
            registry.labels[label] = LabelRecord(
                label=label,
                track=track,
                genotype=genotype,
                created_seq=registry._next_seq(),
                created_run=run_id,
                created_date=run_date,
                split_from=split_from,
            )
            new_labels.add(label)
            log(
                "split" if split_from else "created",
                label,
                split_from,
                f"{len(members)} members" + (f"; split from {split_from}" if split_from else ""),
            )
        for member in members:
            key = (track, *member)
            previous = prior.get(member, "")
            if registry.membership.get(key) != label:
                registry.membership[key] = label
            if previous != label and (label not in new_labels or previous):
                log("joined", label, previous, member[1])
            assigned[member] = label

    registry.events.extend(events)
    result = {member: assigned.get(member, "") for member in partition}
    return LabelResult(result, events, new_labels)


def label_table(registry: Registry, track: str, present: Mapping[Member, str]) -> list[dict[str, str | int]]:
    """Describe every active label of ``track`` that has members in this run."""
    totals = Counter(label for (t, _g, _s), label in registry.membership.items() if t == track)
    now = Counter(label for label in present.values() if label)
    rows = []
    for label in sorted(now):
        record = registry.labels.get(label)
        if record is None:  # a new label in a dry run that was never registered
            rows.append({"label": label, "track": track, "n_members_this_run": now[label]})
            continue
        aliases = sorted(
            other for other in registry.labels if other != label and registry.resolve(other) == label
        )
        rows.append(
            {
                "label": label,
                "track": track,
                "genotype": record.genotype,
                "created_date": record.created_date,
                "split_from": record.split_from,
                "absorbed_labels": ";".join(aliases),
                "n_members_this_run": now[label],
                "n_members_total": totals[label],
            }
        )
    return rows


# --- persistence -----------------------------------------------------------------------

_LABEL_FIELDS = [f.name for f in fields(LabelRecord)]
_EVENT_FIELDS = ["run_id", "date", "track", "event", "label", "related_label", "detail"]
_RUN_FIELDS = ["run_id", "date", "track", "n_samples", "n_labelled", "committed", "outdir"]


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv_atomic(path: Path, fieldnames: list[str], rows: list[Mapping[str, object]]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def load_registry(directory: Path) -> Registry:
    registry = Registry()
    registry.tracks = {row["track"]: row["definition"] for row in _read_csv(directory / TRACKS_FILE)}
    for row in _read_csv(directory / LABELS_FILE):
        row["created_seq"] = int(row["created_seq"])
        registry.labels[row["label"]] = LabelRecord(**row)
    for row in _read_csv(directory / MEMBERSHIP_FILE):
        registry.membership[(row["track"], row["genotype"], row["sample_id"])] = row["label"]
    registry.events = _read_csv(directory / EVENTS_FILE)
    registry.runs = _read_csv(directory / RUNS_FILE)
    return registry


def save_registry(registry: Registry, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    _write_csv_atomic(
        directory / TRACKS_FILE,
        ["track", "definition"],
        [{"track": t, "definition": d} for t, d in sorted(registry.tracks.items())],
    )
    _write_csv_atomic(
        directory / LABELS_FILE,
        _LABEL_FIELDS,
        [asdict(r) for r in sorted(registry.labels.values(), key=lambda r: r.created_seq)],
    )
    _write_csv_atomic(
        directory / MEMBERSHIP_FILE,
        ["track", "genotype", "sample_id", "label"],
        [
            {"track": t, "genotype": g, "sample_id": s, "label": label}
            for (t, g, s), label in sorted(registry.membership.items())
        ],
    )
    _write_csv_atomic(directory / EVENTS_FILE, _EVENT_FIELDS, registry.events)
    _write_csv_atomic(directory / RUNS_FILE, _RUN_FIELDS, registry.runs)


@contextlib.contextmanager
def locked(directory: Path) -> Iterator[None]:
    """Serialize registry updates (e.g. two scheduled runs on one server)."""
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / LOCK_FILE).open("w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def label_run(
    directory: Path,
    track: str,
    definition: Mapping[str, object],
    partition: Mapping[Member, str],
    *,
    run_id: str,
    run_date: str,
    outdir: str,
    commit: bool,
    prefix: str = DEFAULT_PREFIX,
) -> tuple[LabelResult, list[dict[str, str | int]], Registry]:
    """Label a partition against the registry in ``directory``; save only if ``commit``.

    Without ``commit`` the registry on disk is untouched and clusters that would receive a
    brand-new label are reported as ``UNASSIGNED`` instead of a name that was never saved.
    """
    with locked(directory) if commit else contextlib.nullcontext():
        registry = load_registry(directory)
        working = registry if commit else copy.deepcopy(registry)
        result = assign_labels(
            working, track, definition, partition, run_id=run_id, run_date=run_date, prefix=prefix
        )
        if not commit:
            def mask(label: str) -> str:
                return UNASSIGNED if label in result.new_labels else label

            result.labels = {m: mask(label) for m, label in result.labels.items()}
            result.events = [
                {**e, "label": mask(e["label"]), "related_label": mask(e["related_label"])}
                for e in result.events
            ]
        working.runs.append(
            {
                "run_id": run_id,
                "date": run_date,
                "track": track,
                "n_samples": str(len(partition)),
                "n_labelled": str(sum(1 for label in result.labels.values() if label and label != UNASSIGNED)),
                "committed": "true" if commit else "false",
                "outdir": outdir,
            }
        )
        table = label_table(working, track, result.labels)
        if commit:
            save_registry(working, directory)
    return result, table, working
