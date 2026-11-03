#!/usr/bin/env python3
"""IAP-only clusters, and how they compare with the main clusters.

The IAP clusters use the relatedness index alone: every pair of dated sequences within a
subtype with an index of at least ``alpha`` is linked, whatever its genetic distance, and
the clusters are the connected components of those links (the same code as the main
clustering). The main clusters (TN93, or SNP) are untouched, so the two sets are
independent and can be compared:

* ``links.iap_information.csv``: every main link with its index and whether the model
  supports it (a pair the model cannot judge is not supported, and ``note`` says why);
* ``links.iap.csv``: every IAP link, and whether the main clustering also linked it;
* ``clusters.comparison.csv``: for each cluster of two or more sequences in either set,
  how the other set sees its members.
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Sequence

from . import iap, workflow

TRACK = "iap"  # the label registry track of the IAP clusters
DEFAULT_LABEL_PREFIX = "HCV_IAP"  # its own prefix, so its label numbers never depend on the main labels
LINKS_FILENAME = "links.iap_information.csv"
EDGES_FILENAME = "links.iap.csv"
CLUSTERS_FILENAME = "clusters.iap.csv"
COMPARISON_FILENAME = "clusters.comparison.csv"
LABELS_FILENAME = "cluster_labels.iap.csv"
CHANGES_FILENAME = "label_changes.iap.csv"

LINK_FIELDS = ("raw_dist", "time_diff_years", "relatedness_index", "supported", "note", "flags")
EDGE_FIELDS = ("source", "target", "distance", "raw_dist", "time_diff_years", "relatedness_index", "in_main", "flags")
CLUSTER_FIELDS = ("sample_id", "genotype", "cluster_id", "source_cluster_id", "cluster_size")
NOTE_NO_DATE = "no_collection_date"  # a main link whose sequence is not in iap_pairs.csv

IDENTICAL, CONTAINED, SPLIT, ABSENT = "identical", "contained", "split", "not_clustered_in_other"


@dataclass
class SupportSummary:
    links: int = 0
    supported_links: int = 0
    iap_only_links: int = 0
    flagged_iap_links: int = 0  # IAP links whose pair is flagged low_coverage or noisy
    main: Counter = field(default_factory=Counter)  # relationship of each main cluster -> count
    iap: Counter = field(default_factory=Counter)


def label_definition(alpha: float) -> dict[str, object]:
    """What the IAP labels are valid for: a different alpha or model fit changes who clusters."""
    return {
        "metric": TRACK,
        "region": iap.CALIBRATED_REGION,
        "region_strategy": "fixed",
        "alpha": alpha,
        "model": [iap.MODEL_A, iap.MODEL_B, iap.MODEL_C],
    }


def _read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or ()), list(reader)


def _write(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_link_tables(
    links: Path, pairs: Path, alpha: float
) -> tuple[list[str], list[dict[str, str]], list[dict[str, str]]]:
    """Annotate the main links with the index, and list every IAP link, in one pass over the pairs."""
    link_fields, main_links = _read(links)
    wanted = {frozenset((link["source"], link["target"])): link for link in main_links}
    found: dict[frozenset[str], dict[str, str]] = {}
    edges: list[dict[str, str]] = []
    with pairs.open(newline="", encoding="utf-8") as handle:
        for pair in csv.DictReader(handle):
            key = frozenset((pair["seqid_1"], pair["seqid_2"]))
            main = wanted.get(key)
            if main is not None:
                found[key] = pair
            index = pair["relatedness_index"]
            if index and float(index) >= alpha:
                edges.append(
                    {
                        "source": pair["seqid_1"],
                        "target": pair["seqid_2"],
                        "distance": main["distance"] if main else "",
                        "raw_dist": pair["raw_dist"],
                        "time_diff_years": pair["time_diff_years"],
                        "relatedness_index": index,
                        "in_main": "TRUE" if main else "FALSE",
                        "flags": pair.get("flags", ""),
                    }
                )
    annotated = []
    for link in main_links:
        pair = found.get(frozenset((link["source"], link["target"])))
        index = pair["relatedness_index"] if pair else ""
        annotated.append(
            {
                **link,
                "raw_dist": pair["raw_dist"] if pair else "",
                "time_diff_years": pair["time_diff_years"] if pair else "",
                "relatedness_index": index,
                "supported": "TRUE" if index and float(index) >= alpha else "FALSE",
                "note": pair["note"] if pair else NOTE_NO_DATE,
                "flags": pair.get("flags", "") if pair else "",
            }
        )
    return [*link_fields, *LINK_FIELDS], annotated, edges


def recluster(
    clusters: Sequence[Mapping[str, str]], edges: Sequence[Mapping[str, str]]
) -> list[dict[str, str | int]]:
    """Connected components of the IAP links, one genotype at a time (as in the main run).

    Every sequence of the main cluster table is a node, so one with no collection date (or
    no IAP link) is a singleton.
    """
    samples: dict[str, list[str]] = defaultdict(list)
    genotype_of: dict[str, str] = {}
    for row in clusters:
        samples[row["genotype"]].append(row["sample_id"])
        genotype_of[row["sample_id"]] = row["genotype"]
    links: dict[str, list[workflow.LinkRow]] = defaultdict(list)
    for edge in edges:
        genotype = genotype_of.get(edge["source"])
        if genotype is not None and genotype_of.get(edge["target"]) == genotype:
            links[genotype].append(workflow.LinkRow(edge["source"], edge["target"]))

    rows: list[dict[str, str | int]] = []
    for genotype in sorted(samples):
        assignments = workflow.connected_components(links[genotype], samples[genotype])
        sizes = Counter(assignments.values())
        for sample_id in sorted(assignments):
            # IAPnnnn, not Cnnnn: an ID then says which set of clusters it belongs to.
            source_cluster_id = assignments[sample_id].replace("C", "IAP", 1)
            rows.append(
                {
                    "sample_id": sample_id,
                    "genotype": genotype,
                    "cluster_id": f"{workflow.safe_genotype_name(genotype)}_{source_cluster_id}",
                    "source_cluster_id": source_cluster_id,
                    "cluster_size": sizes[assignments[sample_id]],
                }
            )
    return rows


def _relationships(
    rows: Sequence[Mapping[str, str]], other: Sequence[Mapping[str, str]], name: str
) -> list[dict[str, str | int]]:
    """How the other cluster set sees each cluster of two or more sequences in ``rows``."""
    members: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        members[row["cluster_id"]].append(row["sample_id"])
    cluster_of = {row["sample_id"]: row["cluster_id"] for row in other}
    other_sizes = Counter(cluster_of.values())
    labelled = any("cluster_label" in row for row in rows)
    label_of = {row["cluster_id"]: row.get("cluster_label", "") for row in rows}
    genotype_of = {row["cluster_id"]: row["genotype"] for row in rows}

    result = []
    for cluster_id in sorted(members):
        if len(members[cluster_id]) < 2:
            continue
        seen = Counter(cluster_of.get(sample, "") for sample in members[cluster_id])
        clustered = {c: n for c, n in seen.items() if c and other_sizes[c] > 1}
        if not clustered:
            relationship = ABSENT
        elif len(seen) == 1:
            relationship = IDENTICAL if other_sizes[next(iter(seen))] == len(members[cluster_id]) else CONTAINED
        else:
            relationship = SPLIT
        result.append(
            {
                "set": name,
                "cluster_id": cluster_id,
                "cluster_label": label_of[cluster_id] if labelled else "",
                "genotype": genotype_of[cluster_id],
                "sequences": len(members[cluster_id]),
                "concordant": "TRUE" if relationship == IDENTICAL else "FALSE",
                "relationship": relationship,
                "counterparts": ";".join(f"{c}({n})" for c, n in sorted(clustered.items(), key=lambda x: (-x[1], x[0]))),
                "unclustered_in_other": sum(n for c, n in seen.items() if not c or other_sizes[c] == 1),
            }
        )
    return result


COMPARISON_FIELDS = ("set", "cluster_id", "cluster_label", "genotype", "sequences", "concordant", "relationship",
                     "counterparts", "unclustered_in_other")


def write_links_and_clusters(
    *, links: Path, pairs: Path, clusters: Path, alpha: float,
    out_links: Path, out_edges: Path, out_clusters: Path,
) -> SupportSummary:
    """Write the annotated main links, the IAP links and the (unlabelled) IAP clusters."""
    link_fields, link_rows, edges = build_link_tables(links, pairs, alpha)
    _write(out_links, link_fields, link_rows)
    _write(out_edges, EDGE_FIELDS, edges)
    _write(out_clusters, CLUSTER_FIELDS, recluster(_read(clusters)[1], edges))
    return SupportSummary(
        links=len(link_rows),
        supported_links=sum(link["supported"] == "TRUE" for link in link_rows),
        iap_only_links=sum(edge["in_main"] == "FALSE" for edge in edges),
        flagged_iap_links=sum(bool(edge["flags"]) for edge in edges),
    )


def write_genotype_tables(
    clusters_csv: Path, edges_csv: Path, path_for: Callable[[str, str], Path]
) -> list[Path]:
    """Split the IAP clusters and links per genotype, in the form of a per-genotype main run.

    ``path_for(genotype, kind)`` names each file (kind is ``clusters`` or ``links``); the
    cluster ids are the genotype-local ones (``C0001``). Returns the clusters tables written.
    """
    node_rows = _read(clusters_csv)[1]
    genotype_of = {row["sample_id"]: row["genotype"] for row in node_rows}
    edges_by_genotype: dict[str, list[dict[str, str]]] = defaultdict(list)
    for edge in _read(edges_csv)[1]:
        genotype = genotype_of.get(edge["source"])
        if genotype is not None:
            edges_by_genotype[genotype].append(edge)
    nodes_by_genotype: dict[str, list[dict[str, str | int]]] = defaultdict(list)
    for row in node_rows:
        nodes_by_genotype[row["genotype"]].append(
            {"sample_id": row["sample_id"], "cluster_id": row["source_cluster_id"], "cluster_size": row["cluster_size"]}
        )
    written = []
    for genotype in sorted(nodes_by_genotype):
        written.append(path_for(genotype, "clusters"))
        _write(written[-1], ("sample_id", "cluster_id", "cluster_size"), nodes_by_genotype[genotype])
        _write(path_for(genotype, "links"), EDGE_FIELDS, edges_by_genotype[genotype])
    return written


def write_comparison(
    main_clusters: Path, iap_clusters: Path, out: Path, summary: SupportSummary
) -> SupportSummary:
    """Compare the two cluster tables (after labelling, so the labels are shown) and count the outcomes."""
    main_rows, iap_rows = _read(main_clusters)[1], _read(iap_clusters)[1]
    rows = _relationships(main_rows, iap_rows, "main") + _relationships(iap_rows, main_rows, "iap")
    _write(out, COMPARISON_FIELDS, rows)
    summary.main = Counter(row["relationship"] for row in rows if row["set"] == "main")
    summary.iap = Counter(row["relationship"] for row in rows if row["set"] == "iap")
    return summary


def describe(summary: SupportSummary, alpha: float) -> str:
    def clusters(counts: Counter) -> str:
        return f"{sum(counts.values())} ({counts[IDENTICAL]} identical in the other set)"

    return (
        f"IAP model (alpha {alpha:g}): {summary.supported_links} of {summary.links} main link(s) supported, "
        f"{summary.iap_only_links} link(s) only the model makes; clusters: main {clusters(summary.main)}, "
        f"IAP {clusters(summary.iap)}"
        + (f"; {summary.flagged_iap_links} of the model's links involve a flagged pair (low coverage or ambiguous)"
           if summary.flagged_iap_links else "")
    )
