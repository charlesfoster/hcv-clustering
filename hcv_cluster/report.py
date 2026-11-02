#!/usr/bin/env python3
"""A one-page progress report (HTML, A4), run on demand with ``hcv-report RUN_DIR``.

It is a standalone command, not a step of the pipeline: point it at a finished run's output
folder when a report is wanted. The page shows the sequence funnel (received, passed QC,
in a cluster), cluster statistics, who is in the data (age, sex and so on, from the
metadata), clusters by location, collection months, a map of the largest clusters and what
changed in the label registry. Sections whose data are missing (no metadata, no registry,
no dates) are left out. It shows counts only: no sample identifiers.

The page is one self-contained HTML file (inline SVG, no scripts or network access), so
it renders and prints (A4, "Save as PDF") in any browser.
"""

from __future__ import annotations

import argparse
import base64
import calendar
import csv
import datetime
import html
import json
import mimetypes
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Collection, Mapping, Sequence

from . import __version__, methods

REPORT_FILENAME = "report.html"
DEFAULT_TITLE = "HCV sequencing & transmission-cluster update"
DATE_FORMATS = ("%d/%m/%Y", "%Y-%m-%d", "%d-%b-%y", "%d-%b-%Y", "%d-%B-%Y", "%d %b %Y", "%d %B %Y")  # d/m/yyyy, yyyy-mm-dd, 5-Mar-24
MAX_LOCATIONS = 20  # locations named in the location table (all of them, unless there are more than fit); the rest are combined
MAP_CLUSTERS = 5  # clusters drawn on the map, keyed a, b, c… in the largest-clusters table
# Metadata columns recognised without being named (compared lower-case, spaces as underscores).
COLUMN_NAMES = {
    "date": ("sample_collection_date", "collection_date", "doc"),  # DOC: date of collection, as the sample is taken at the interview
    "location": ("location", "site", "setting"),
    "age": ("age",),
    "sex": ("sex", "gender"),
    "diagnosis": ("how_diagnosed", "reason_for_test"),
    "injecting": ("injecting", "injecting_status", "injecting_drug_use", "idu"),
    "indigenous": ("indigenous_status", "indigenous", "aboriginal_status"),
}
INJECTING_BINS = ("Last month", "Last 6 months", "Last 2 years", "More than 2 years ago", "Ever (recency unknown)", "Never")
MISSING = {"", "n/a", "na", "unknown", "(missing)", "none", "pending"}

BLUE, RED, INK, MUTE, LINE, GREY = "#1d4f9f", "#d6292f", "#1b2430", "#5d6877", "#d9dee6", "#9aa3b2"
CATEGORY_COLOURS = (BLUE, RED, "#2e9b57", "#e0a526", "#7b4fb3", "#17a2b8", "#e377a2")

_CSS = """
@page { size: A4; margin: 0 }
* { box-sizing: border-box }
body { margin: 0; background: #eef0f4; font: 10.5px/1.4 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: @INK@ }
.page { width: 210mm; min-height: 296mm; margin: 0 auto; background: #fff; padding: 9mm 12mm 6mm; display: flex; flex-direction: column; gap: 5px }
header { display: flex; align-items: center; gap: 14px; border-bottom: 3px solid; border-image: linear-gradient(90deg,@RED@,@BLUE@) 1; padding-bottom: 8px }
header img { height: 62px }
header h1 { margin: 0; font-size: 17px; line-height: 1.15 }
header p { margin: 2px 0 0; color: @MUTE@; font-size: 11px }
header p.range { margin-top: 1px; color: @INK@; font-size: 10.5px; font-weight: 600 }
header .meta { margin-left: auto; text-align: right; font-size: 10.5px; color: @MUTE@; white-space: nowrap }
header .meta b { display: block; font-size: 13px; color: @INK@ }
.kpis { display: flex; gap: 8px }
.kpi { flex: 1; border: 1px solid @LINE@; border-top: 3px solid @BLUE@; border-radius: 6px; padding: 5px 9px }
.kpi.red { border-top-color: @RED@ }
.kpi b { display: block; font-size: 19px; line-height: 1.1 }
.kpi span { color: @MUTE@; font-size: 9.5px }
.funnel { display: flex; align-items: stretch }
.funnel .step { flex: 1; border: 1px solid @LINE@; border-top: 3px solid @BLUE@; border-radius: 6px; padding: 5px 10px }
.funnel .step.red { border-top-color: @RED@ }
.funnel .step b { display: block; font-size: 25px; line-height: 1.1 }
.funnel .step span { color: @MUTE@; font-size: 9.5px }
.funnel .arrow { display: flex; flex-direction: column; justify-content: center; align-items: center; width: 44px; color: @MUTE@ }
.funnel .arrow span { font-size: 17px; line-height: 1 } .funnel .arrow small { font-size: 9px }
.band { display: flex; align-items: center; gap: 8px; margin: 2px 0 6px; font-size: 11.5px; font-weight: 700; letter-spacing: .08em; text-transform: uppercase; color: @BLUE@ }
.band::after { content: ""; flex: 1; height: 1px; background: @LINE@ }
section > .funnel, section > .kpis { margin-bottom: 9px }
.cols { display: grid; gap: 18px; align-items: start } .stack { display: flex; flex-direction: column; gap: 10px }
h2 { margin: 0 0 5px; font-size: 11px; font-weight: 700; color: @INK@ }
th small { font-weight: 400; color: @MUTE@; font-size: 8.5px; margin-left: 3px }
table { width: 100%; border-collapse: collapse; font-size: 10px }
th { text-align: left; font-weight: 600; color: @MUTE@; border-bottom: 1.5px solid @INK@; padding: 3px 4px; font-size: 9px; white-space: nowrap }
th.num, td.num { text-align: right; font-variant-numeric: tabular-nums; padding-left: 8px }
td { padding: 1px 4px; border-bottom: 1px solid @LINE@; white-space: nowrap }
td.last { white-space: nowrap } i.dot { display: inline-block; width: 6px; height: 6px; border-radius: 50%; margin-right: 4px; vertical-align: 1px }
i.dot.on { background: #2e9b57 } i.dot.off { border: 1.2px solid @GREY@; width: 4.6px; height: 4.6px }
td.gt { font-weight: 700 } td.key { font-weight: 700; color: @BLUE@; width: 14px }
tr.total td { font-weight: 700; border-top: 1.5px solid @INK@; border-bottom: 0 }
td.barcell { width: 24% }
td.barcell span { display: block; height: 8px; background: @BLUE@; border-radius: 2px }
table.demo { table-layout: fixed; margin-bottom: 9px } table.demo td.lab { overflow: hidden; text-overflow: ellipsis; padding-left: 4px }
table.demo th.title { font-size: 10.5px; color: @INK@; padding-left: 4px; overflow: visible }
.demos { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; align-items: start }
.funnel .step b em { font-style: normal; font-size: 11px; font-weight: 600; color: #1c6b3a }
b.new { color: @RED@; font-size: 8.5px } td.up { color: #1c6b3a; font-weight: 600 }
svg { width: 100%; height: auto; display: block }
.legend { font-size: 9px; color: @MUTE@; display: flex; gap: 3px 8px; flex-wrap: wrap; margin-top: 2px }
.legend i { display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 3px }
.note { font-size: 9.5px; color: @MUTE@; margin: 2px 0 0 }
.callout { border-left: 3px solid @RED@; background: #fdf3f3; padding: 6px 9px; border-radius: 0 5px 5px 0; font-size: 10px }
footer { margin-top: auto; border-top: 1px solid @LINE@; padding-top: 6px; color: @MUTE@; font-size: 9px; line-height: 1.45 }
@media print { body { background: #fff } .page { margin: 0 } }
"""


def _css() -> str:
    out = _CSS
    for name, value in (("INK", INK), ("RED", RED), ("BLUE", BLUE), ("MUTE", MUTE), ("LINE", LINE), ("GREY", GREY)):
        out = out.replace(f"@{name}@", value)
    return out


def _e(value: object) -> str:
    return html.escape(str(value))


def _read(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _long_date(d: datetime.date) -> str:
    return f"{d.day} {d:%B %Y}"


def _parse_date(value: str | None) -> datetime.date | None:
    text = (value or "").strip()
    for fmt in DATE_FORMATS:
        try:
            return datetime.datetime.strptime(text, fmt).date()
        except ValueError:
            continue
    return None


def _pick(fieldnames: Sequence[str], explicit: str | None, kind: str) -> str | None:
    """The metadata column for ``kind`` (named by the user, else recognised by name), or None."""
    if explicit is not None:
        if explicit not in fieldnames:
            raise ValueError(f"'{explicit}' is not a column of the sample table (available: {', '.join(fieldnames)})")
        return explicit
    return next((f for name in COLUMN_NAMES[kind] for f in fieldnames if f.lower().replace(" ", "_") == name), None)  # earlier names win


def _bars(values: Sequence[int], labels: Sequence[str], colour: str, every: int = 1,
          width: int = 330, height: int = 100) -> str:
    left, bottom, top, right = 26, 20, 8, 14  # the right margin keeps the last label inside the picture
    step = (width - left - right) / len(values)
    peak = max(values) or 1
    out = [f'<svg viewBox="0 0 {width} {height}" role="img">']
    for frac in (0.5, 1.0):
        y = height - bottom - (height - bottom - top) * frac
        out.append(f'<line x1="{left}" x2="{width - right}" y1="{y:.1f}" y2="{y:.1f}" stroke="{LINE}" stroke-width="1"/>')
        out.append(f'<text x="{left - 4}" y="{y + 3:.1f}" text-anchor="end" font-size="8.5" fill="{MUTE}">{peak * frac:.0f}</text>')
    for i, value in enumerate(values):
        h = (height - bottom - top) * value / peak
        x = left + i * step + step * 0.15
        out.append(f'<rect x="{x:.1f}" y="{height - bottom - h:.1f}" width="{step * 0.7:.1f}" height="{h:.1f}" rx="1.5" fill="{colour}"/>')
        if i % every == 0:
            out.append(f'<text x="{x + step * 0.35:.1f}" y="{height - 7}" text-anchor="middle" font-size="8.5" fill="{MUTE}">{_e(labels[i])}</text>')
    out.append("</svg>")
    return "".join(out)


def _timelines(rows: Sequence[tuple[str, Sequence[tuple[datetime.date, str, bool]]]], span: Sequence[tuple[int, int]],
               labels: Sequence[str], every: int = 4, width: int = 330, height: int = 125) -> str:
    """One row per cluster: a dot for each sequence at its collection date, on the month axis of ``_bars``.

    ``rows`` is (key, [(date, colour, new since the last report)]). A line joins the first and last sequence.
    """
    left, bottom, top, right = 26, 20, 8, 14  # the same margins as _bars, so the two charts share an axis
    step = (width - left - right) / len(span)
    rowh = (height - bottom - top) / max(len(rows), 1)

    def x_of(d: datetime.date) -> float:  # inside the month's bar: left edge to right edge over the days of the month
        month = (d.year - span[0][0]) * 12 + d.month - span[0][1]
        return left + (month + 0.15 + 0.7 * (d.day - 1) / calendar.monthrange(d.year, d.month)[1]) * step

    def month_of(d: datetime.date) -> int:
        return (d.year - span[0][0]) * 12 + d.month - span[0][1]

    out = [f'<svg viewBox="0 0 {width} {height}" role="img">']
    for j, (key, dots) in enumerate(rows):
        y = top + rowh * (j + 0.5)
        out.append(f'<line x1="{left}" x2="{width - right}" y1="{y:.1f}" y2="{y:.1f}" stroke="{LINE}" stroke-width="0.6"/>')
        out.append(f'<text x="{left - 11}" y="{y + 3.5:.1f}" text-anchor="middle" font-size="10" font-weight="700" fill="{BLUE}">{_e(key)}</text>')
        inside = [dot for dot in dots if 0 <= month_of(dot[0]) < len(span)]
        earlier = any(month_of(d) < 0 for d, _c, _n in dots)  # sequences before the chart: the line runs off its left edge
        xs = [x_of(d) for d, _c, _n in inside]
        start = left + 2 if earlier or not xs else min(xs)
        end = max(xs) if xs else left + 8
        out.append(f'<line x1="{start:.1f}" x2="{end:.1f}" y1="{y:.1f}" y2="{y:.1f}" stroke="{GREY}" stroke-width="1.4"/>')
        if earlier:
            out.append(f'<polygon points="{left - 1},{y:.1f} {left + 5},{y - 3:.1f} {left + 5},{y + 3:.1f}" fill="{GREY}"/>')
        for (d, colour, new), x in sorted(zip(inside, xs), key=lambda pair: pair[1]):
            if new:
                out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.2" fill="none" stroke="{INK}" stroke-width="1.2"/>')
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{colour}" stroke="#fff" stroke-width="0.7"/>')
    for i, label in enumerate(labels):
        if i % every == 0:
            out.append(f'<text x="{left + (i + 0.5) * step:.1f}" y="{height - 7}" text-anchor="middle" font-size="8.5" fill="{MUTE}">{_e(label)}</text>')
    out.append("</svg>")
    return "".join(out)


def _cluster_map(
    members: Mapping[str, list[dict[str, str]]],
    links: Sequence[Mapping[str, str]],
    category_of: Mapping[str, str],
    names: Mapping[str, str],
    *,
    top: int = 12,
    labelled: int = 3,
    label_size: float = 7,
    added: Collection[str] = (),
    width: int = 330,
    height: int = 128,
) -> tuple[str, str] | None:
    """Inline-SVG map of the largest clusters (singletons are not drawn); None if it cannot be drawn."""
    try:
        from . import viz
    except ImportError:  # the plotting extras (networkx, plotly) are optional
        return None
    chosen = list(names)[:top]  # the clusters the caller chose to show, in order
    nodes = sorted(
        ({"sample_id": r["sample_id"], "cluster_id": c} for c in chosen for r in members[c]),
        key=lambda node: node["sample_id"],
    )
    ids = {n["sample_id"] for n in nodes}
    # TN93 writes its links in a different order each run; sort them so the page is reproducible.
    edges = sorted(
        (dict(link) for link in links if link["source"] in ids and link["target"] in ids),
        key=lambda link: tuple(sorted((link["source"], link["target"]))),
    )
    positions = viz.compute_network_layout(nodes, edges, mode="components", node_spacing="compact")
    if not positions:
        return None
    xs, ys = [p[0] for p in positions.values()], [p[1] for p in positions.values()]
    pad = 12
    scale = min((width - 2 * pad) / ((max(xs) - min(xs)) or 1), (height - 2 * pad) / ((max(ys) - min(ys)) or 1))

    def at(sample: str) -> tuple[float, float]:
        x, y = positions[sample]
        return pad + (x - min(xs)) * scale + ((width - 2 * pad) - (max(xs) - min(xs)) * scale) / 2, \
            pad + (y - min(ys)) * scale + ((height - 2 * pad) - (max(ys) - min(ys)) * scale) / 2

    out = [f'<svg viewBox="0 0 {width} {height}" role="img">']
    for edge in edges:
        first, second = sorted((edge["source"], edge["target"]))  # TN93 may report a pair either way round
        (x1, y1), (x2, y2) = at(first), at(second)
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{GREY}" stroke-width="0.8"/>')
    for node in nodes:
        x, y = at(node["sample_id"])
        if node["sample_id"] in added:  # a dark halo, so it shows against every location colour
            out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.2" fill="none" stroke="{INK}" stroke-width="1.2"/>')
        out.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{category_of[node["sample_id"]]}" stroke="#fff" stroke-width="0.7"/>')
    for cluster in chosen[:labelled]:  # name only the largest few; the rest would crowd the map
        points = [at(r["sample_id"]) for r in members[cluster]]
        x = sum(p[0] for p in points) / len(points)
        y = min(p[1] for p in points) - 6
        out.append(
            f'<text x="{x:.1f}" y="{max(y, 7):.1f}" text-anchor="middle" font-size="{label_size}" font-weight="600" fill="{INK}" '
            f'stroke="#fff" stroke-width="2.5" paint-order="stroke">{_e(names[cluster])}</text>'
        )
    out.append("</svg>")
    return "".join(out), f"{len(chosen)} largest of {len(members)} clusters"


def _method_text(params: Mapping[str, object]) -> str:
    distance = params.get("distance") or "tn93"
    primary = "snp" if distance == "snp" else "tn93"
    region = str(params.get("region") or "e1-e2")
    raw = params.get("threshold")
    try:
        threshold = methods.resolve_threshold(region, float(raw) if raw not in (None, "", "None", "null") else None, distance=primary)
        rule = f"pairs with {'SNP p' if primary == 'snp' else 'TN93 '}distance ≤ {threshold:g} are linked"
    except (ValueError, TypeError):
        rule = "pairs closer than the clustering threshold are linked"
    text = (
        f"Sequences are assigned a subtype, aligned to the {region.upper()} region and clustered within subtype: "
        f"{rule} and clusters are the connected groups. "
    )
    return text


def _logo_html(logo: Path | None) -> str:
    if logo is None:
        return ""
    mime = mimetypes.guess_type(logo.name)[0] or "image/png"
    return f'<img src="data:{mime};base64,{base64.b64encode(logo.read_bytes()).decode()}" alt="logo">'


WINDOW_CAP_MONTHS = 36  # the time charts show at most this many months unless told otherwise


def _months_before(d: datetime.date, n: int) -> datetime.date:
    year, month = divmod(d.year * 12 + d.month - 1 - n, 12)
    return datetime.date(year, month + 1, min(d.day, calendar.monthrange(year, month + 1)[1]))


def _window_start(first: tuple[int, int], latest: tuple[int, int], chart_from: datetime.date | None,
                  chart_months: int | None) -> tuple[int, int]:
    """First month of the time charts: a given date, or ``chart_months`` back from the latest month, else the whole
    study up to ``WINDOW_CAP_MONTHS``. Never before the first sample, never after the latest."""
    if chart_from is not None:
        start = (chart_from.year, chart_from.month)
    else:
        index = latest[0] * 12 + latest[1] - 1 - ((chart_months or WINDOW_CAP_MONTHS) - 1)
        start = (index // 12, index % 12 + 1)
    return min(max(start, first), latest)


def _delta(now: int, before: int | None) -> str:
    return f" <em>+{now - before}</em>" if before is not None and now > before else ""


def _pct(a: int, b: int) -> str:
    return f"{100 * a / b:.0f}%" if b else "–"


def _band(age: int) -> str:
    return "<25" if age < 25 else "25–34" if age < 35 else "35–44" if age < 45 else "45–54" if age < 55 else "55+"


def _present(value: str | None) -> str:
    value = (value or "").strip()
    return "Not recorded" if value.lower() in MISSING else value


def _breakdown(title: str, keyed: Sequence[tuple[str, bool]], min_cell: int, order: Sequence[str] | None = None,
               note: str = "") -> str:
    """A small table of (label, sequences, % in a cluster). Tiny categories are pooled; missing is shown, in grey."""
    n_by: Counter[str] = Counter()
    c_by: Counter[str] = Counter()
    for label, in_cluster in keyed:
        n_by[label] += 1
        c_by[label] += in_cluster
    named = {k: v for k, v in n_by.items() if k != "Not recorded"}
    if order is not None:
        order = [*order, *(k for k in sorted(named) if k not in order)]  # a value nobody anticipated is still shown
    else:
        order = [k for k, _ in sorted(named.items(), key=lambda kv: (-kv[1], kv[0]))]
        small = [k for k in order if named[k] < min_cell]
        if len(small) > 1:
            for k in small:
                n_by["Other"] += n_by.pop(k)
                c_by["Other"] += c_by.pop(k)
            order = [k for k in order if k not in small] + ["Other"]
    labels = [k for k in order if n_by.get(k)] + (["Not recorded"] if n_by.get("Not recorded") else [])
    peak = max((n_by[k] for k in labels), default=1)
    body = ""
    for k in labels:
        missing = k == "Not recorded"
        share = 100 * n_by[k] / peak
        shade = "#e6e9ef" if missing else "#d3e0f3"  # the bar is the label cell's background, so every panel lines up
        style = f"background:linear-gradient(90deg,{shade} {share:.0f}%,transparent {share:.0f}%)" + (f";color:{MUTE}" if missing else "")
        body += (f"<tr><td class='lab' style='{style}'>{_e(k)}</td>"
                 f"<td class='num'>{n_by[k]}</td><td class='num'>{_pct(c_by[k], n_by[k])}</td></tr>")
    small = f" <small>{_e(note)}</small>" if note else ""
    return (f"<table class='demo'><colgroup><col style='width:60%'><col style='width:15%'><col style='width:25%'></colgroup>"
            f"<tr class='hd'><th class='title'>{_e(title)}{small}</th><th class='num'>Seqs</th><th class='num'>In cluster</th></tr>"
            f"{body}</table>")


def _compare(since: Path, rows: Sequence[Mapping[str, str]], members: Mapping[str, list[dict[str, str]]],
             name_of: Mapping[str, str]) -> dict:
    """What changed since the earlier run in ``since``, matching clusters by their persistent label."""
    earlier = _read(since / "metadata.surviving_samples.csv") or _read(since / "clusters.csv")
    if not earlier:
        raise FileNotFoundError(f"{since} has no metadata.surviving_samples.csv: is it a finished run?")
    if not any(r.get("cluster_label") not in (None, "", "UNASSIGNED") for r in rows):
        raise ValueError("--since needs persistent cluster labels: run both analyses with the same --registry (and --update-registry)")
    was = Counter(r["cluster_label"] for r in earlier if int(r["cluster_size"]) > 1 and r.get("cluster_label") not in (None, "", "UNASSIGNED"))
    now = {c: len(m) for c, m in members.items()}
    change = {c: now[c] - was.get(name_of[c], 0) for c in members}
    info_file = since / "run_info.json"
    created = json.loads(info_file.read_text()).get("created", "") if info_file.exists() else ""
    return {
        "date": datetime.date.fromisoformat(created[:10]) if created else None,
        "received": len(_read(since / "sample_audit.csv")) or None,
        "passed": len(earlier),
        "clustered": sum(int(r["cluster_size"]) > 1 for r in earlier),
        "new_ids": {r["sample_id"] for r in rows} - {r["sample_id"] for r in earlier},
        "new": {c for c in members if name_of[c] not in was},
        "grown": {c for c, d in change.items() if d > 0 and name_of[c] in was},
        "change": change,
    }


def build_report(
    run_dir: Path,
    out: Path | None = None,
    *,
    title: str = DEFAULT_TITLE,
    logo: Path | None = None,
    location_column: str | None = None,
    date_column: str | None = None,
    age_column: str | None = None,
    sex_column: str | None = None,
    diagnosis_column: str | None = None,
    injecting_column: str | None = None,
    indigenous_column: str | None = None,
    min_cell: int = 5,
    since: Path | None = None,
    chart_from: datetime.date | None = None,
    chart_months: int | None = None,
    active_months: int = 6,
    active_only: bool = False,
    report_date: datetime.date | None = None,
) -> Path:
    """Write the one-page report for the run in ``run_dir`` (default ``run_dir/report.html``).

    Reads ``sample_audit.csv`` and ``metadata.surviving_samples.csv`` (or ``clusters.csv``)
    and, where present, ``links.csv``, ``label_changes.csv`` and ``run_info.json``. Metadata
    columns are recognised by name (see ``COLUMN_NAMES``) unless named; a column the run's
    table lacks leaves its section out. ``since`` is the output folder of an earlier run
    (made with the same label registry): clusters are then flagged new or grown since it.
    The two time charts start at ``chart_from``, or ``chart_months`` back from the latest sample;
    by default they show the whole study up to 36 months. A cluster is active if it has a sample
    in the last ``active_months`` months of the data; ``active_only`` limits the map, timelines and
    largest-clusters table to those (they otherwise show the largest clusters of all time).
    """
    out = out or run_dir / REPORT_FILENAME
    audit = _read(run_dir / "sample_audit.csv")
    rows = _read(run_dir / "metadata.surviving_samples.csv") or _read(run_dir / "clusters.csv")
    if not audit or not rows:
        raise FileNotFoundError(f"{run_dir} has no sample_audit.csv and metadata.surviving_samples.csv: is it a finished run?")
    links = _read(run_dir / "links.csv")
    events_file = run_dir / "label_changes.csv"
    events = _read(events_file)
    info_file = run_dir / "run_info.json"
    info = json.loads(info_file.read_text()) if info_file.exists() else {}
    params = info.get("parameters", {})
    report_date = report_date or datetime.date.today()
    created = info.get("created", "")
    cutoff = datetime.date.fromisoformat(created[:10]) if created else report_date

    fields = list(rows[0])
    loc = _pick(fields, location_column, "location")
    age_col, sex_col = _pick(fields, age_column, "age"), _pick(fields, sex_column, "sex")
    diag_col, inj_col = _pick(fields, diagnosis_column, "diagnosis"), _pick(fields, injecting_column, "injecting")
    ind_col = _pick(fields, indigenous_column, "indigenous")
    date_col = _pick(fields, date_column, "date")

    received = len(audit)
    passed = sum(r["overall_status"] == "pass" for r in audit)
    in_cluster = [int(r["cluster_size"]) > 1 for r in rows]
    clustered = [r for r, c in zip(rows, in_cluster) if c]
    members: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in clustered:
        members[r["cluster_id"]].append(r)
    spanning = {c for c, m in members.items() if len({r[loc] for r in m}) > 1} if loc else set()
    in_spanning = sum(len(members[c]) for c in spanning)
    sizes = sorted((len(m) for m in members.values()), reverse=True)
    name_of = {c: (m[0].get("cluster_label") or c) if m[0].get("cluster_label") != "UNASSIGNED" else c for c, m in members.items()}
    growth = _compare(since, rows, members, name_of) if since is not None else None

    # the two rows of headline numbers
    collected = received  # placeholder: the number of samples collected is not in the run's files yet
    before = growth and growth["received"]
    funnel = [(collected, "samples collected", None, before),
              (received, "sequenced", _pct(received, collected), before),
              (passed, "passed QC and genotyped", _pct(passed, received), growth and growth["passed"]),
              (len(clustered), "in a cluster", _pct(len(clustered), passed), growth and growth["clustered"])]
    funnel_html = "".join(
        (f"<div class='arrow'><span>→</span><small>{share}</small></div>" if share else "")
        + f"<div class='step{' red' if i == len(funnel) - 1 else ''}'><b>{n}{_delta(n, was)}</b><span>{_e(label)}</span></div>"
        for i, (n, label, share, was) in enumerate(funnel)
    )
    kinds = Counter(ev["event"] for ev in events)
    joined = len({ev["label"] for ev in events if ev["event"] == "joined"})  # existing clusters this run added sequences to
    labelled = len({r["cluster_label"] for r in clustered if r.get("cluster_label") not in (None, "", "UNASSIGNED")})
    first_run = bool(events) and kinds["created"] >= labelled > 0 and not (joined or kinds["merged"] or kinds["split"])
    # the same tiles in every report; the last two are "–" when there is nothing earlier to compare with
    if growth:
        new_clusters, grown_clusters = len(growth["new"]), len(growth["grown"])
    elif events and not first_run:
        new_clusters, grown_clusters = kinds["created"], joined
    else:
        new_clusters = grown_clusters = "–"
    tiles = [(len(members), "clusters (2+ sequences)"), (sum(s >= 3 for s in sizes), "clusters of 3+ sequences"),
             (sizes[0] if sizes else 0, "sequences in largest cluster")]
    if loc:
        tiles.append((len(spanning), f"across 2+ locations ({_pct(len(spanning), len(members))})"))
    tiles += [(new_clusters, "new clusters since last report"), (grown_clusters, "existing clusters that grew")]
    tiles_html = "".join(f"<div class='kpi red'><b>{n}</b><span>{_e(label)}</span></div>" for n, label in tiles)

    # who is in the data
    demo: list[str] = []
    if age_col:
        ages = [int(r[age_col]) if r[age_col].strip().isdigit() else None for r in rows]
        known = sorted(a for a in ages if a is not None)
        note = ""
        if len(known) >= 4:
            q1, _, q3 = statistics.quantiles(known, n=4)
            note = f"median {statistics.median(known):g}, IQR {q1:g}–{q3:g}"
        demo.append(_breakdown("Age", [(_band(a) if a is not None else "Not recorded", c) for a, c in zip(ages, in_cluster)],
                               min_cell, order=["<25", "25–34", "35–44", "45–54", "55+"], note=note))
    sex_names = {"M": "Male", "F": "Female"}
    indigenous = lambda v: "Not recorded" if v == "Not recorded" else ("Neither" if v.lower() == "neither" else "Aboriginal and/or TSI")
    for column, heading, names, order in (
        (sex_col, "Sex", sex_names, None),
        (ind_col, "Indigenous status", {}, ["Aboriginal and/or TSI", "Neither"]),
        (inj_col, "Injecting drug use", {}, INJECTING_BINS),
        (diag_col, "How diagnosed", {}, None),
    ):
        if column:
            keyed = [(_present(r[column]), c) for r, c in zip(rows, in_cluster)]
            keyed = [((indigenous(v) if column == ind_col else names.get(v, v)), c) for v, c in keyed]
            demo.append(_breakdown(heading, keyed, min_cell, order=order))
    first_column: list[str] = []  # fill the left column to about half the table rows, so the two end level
    for panel in demo:
        if sum(p.count("<tr>") for p in first_column) < sum(p.count("<tr>") for p in demo) / 2:
            first_column.append(panel)
    second_column = [p for p in demo if p not in first_column]
    demo_html = (f"<div><h2>Who is in the data</h2><div class='demos'><div>{''.join(first_column)}</div><div>{''.join(second_column)}</div></div>"
                 f"<p class='note'>Passing QC. Bars: relative size. In cluster: share of the group in a cluster.</p></div>") if demo else ""

    # by location
    group_of: dict[str, str] = {}
    location_html = ""
    if loc:
        per_location = Counter(r[loc] for r in rows)
        named = [k for k, _ in sorted(per_location.items(), key=lambda kv: (-kv[1], kv[0]))][:MAX_LOCATIONS]
        group_of = {k: (k if k in named else "other") for k in per_location}
        stats: dict[str, dict] = defaultdict(lambda: {"seqs": 0, "clustered": 0, "clusters": set()})
        for r, c in zip(rows, in_cluster):
            s = stats[group_of[r[loc]]]
            s["seqs"] += 1
            if c:
                s["clustered"] += 1
                s["clusters"].add(r["cluster_id"])
        n_other = sum(g == "other" for g in group_of.values())
        body = "".join(
            f"<tr><td class='gt'>{_e(f'Other sites ({n_other})' if name == 'other' else name)}</td><td class='num'>{s['seqs']}</td>"
            f"<td class='num'>{s['clustered']}</td><td class='num'>{_pct(s['clustered'], s['seqs'])}</td>"
            f"<td class='num'>{len(s['clusters'])}</td><td class='num'>{len(s['clusters'] & spanning)}</td></tr>"
            for name, s in sorted(stats.items(), key=lambda kv: (kv[0] == "other", -kv[1]["seqs"], kv[0])))
        body += (f"<tr class='total'><td class='gt'>All</td><td class='num'>{len(rows)}</td><td class='num'>{len(clustered)}</td>"
                 f"<td class='num'>{_pct(len(clustered), len(rows))}</td><td class='num'>{len(members)}</td><td class='num'>{len(spanning)}</td></tr>")
        location_html = (
            "<div><h2>By location</h2><table><tr><th>Location</th><th class='num'>Seqs</th><th class='num'>In cluster</th>"
            "<th class='num'>%</th><th class='num'>Clusters</th><th class='num'>Multi-site</th></tr>"
            f"{body}</table><p class='note'>A cluster counts at every location it includes (multi-site: it includes another location too), "
            f"so rows add up to more than the total.</p></div>")

    # collection months, largest clusters, cluster map
    if chart_months is not None and chart_months < 1 or active_months < 1:
        raise ValueError("--chart-months and --active-months must be at least 1")
    date_of = {r["sample_id"]: _parse_date(r.get(date_col)) for r in rows} if date_col else {}
    latest = max((d for d in date_of.values() if d), default=None)
    cluster_dates = {c: sorted(d for r in m if (d := date_of.get(r["sample_id"]))) for c, m in members.items()}
    active = {c for c, ds in cluster_dates.items() if ds and ds[-1] >= _months_before(latest, active_months)} if latest else set()
    if active_only and not latest:
        raise ValueError("--active-only needs collection dates: the run's table has no usable date column (see --date-column)")
    # the largest clusters of all time, unless only the active ones are wanted
    top = sorted((c for c in members if c in active or not active_only), key=lambda c: (-len(members[c]), c))
    key_of = dict(zip(top[:MAP_CLUSTERS], "abcdefghij"))  # the clusters on the map and the timelines, and their keys in the table
    shown = [r for c in key_of for r in members[c]]
    by = (lambda r: r[loc]) if loc else (lambda r: r["genotype"])
    order = [v for v, _ in Counter(by(r) for r in shown).most_common(len(CATEGORY_COLOURS))]
    palette = {v: CATEGORY_COLOURS[i] for i, v in enumerate(order)}
    value_of = {r["sample_id"]: by(r) if by(r) in palette else "other" for r in shown}
    palette["other"] = GREY
    new_ids = growth["new_ids"] if growth else set()

    months_html = ""
    range_text = ""  # what the data cover, for the header
    all_months = Counter((d.year, d.month) for d in date_of.values() if d)
    if all_months:
        first, last = _window_start(min(all_months), max(all_months), chart_from, chart_months), max(all_months)
        span = [(y, m) for y in range(first[0], last[0] + 1) for m in range(1, 13) if first <= (y, m) <= last]
        labels = [datetime.date(y, m, 1).strftime("%b") + (f" ’{str(y)[2:]}" if m == 1 or i == 0 else "") for i, (y, m) in enumerate(span)]
        hidden = sum(n for k, n in all_months.items() if k < first)
        range_text = f"Samples collected {datetime.date(*min(all_months), 1):%b %Y} to {datetime.date(*last, 1):%b %Y}"
        if hidden:  # only the two time charts are windowed; everything else counts every sequence
            range_text += f" · time charts from {datetime.date(*first, 1):%b %Y}"
        earlier = f"; {hidden} earlier not shown" if hidden else ""
        months_html = (f"<div><h2>Sequences by collection month</h2>{_bars([all_months.get(k, 0) for k in span], labels, BLUE, every=4, height=125)}"
                       f"<p class='note'>{sum(all_months.values()) - hidden} sequences, {datetime.date(*first, 1):%b %Y} to {datetime.date(*last, 1):%b %Y}{earlier}.</p></div>")
        lines = [(key_of[c], [(d, palette[value_of[r["sample_id"]]], r["sample_id"] in new_ids)
                              for r in members[c] if (d := date_of.get(r["sample_id"]))]) for c in key_of]
        lines = [(key, dots) for key, dots in lines if dots]
        if lines:
            months_html += (f"<div><h2>Cluster timelines</h2>{_timelines(lines, span, labels, height=125)}"
                            f"<p class='note'>Dots are sequences; rows and colours match the map"
                            f"{'; circled = new' if growth else ''}{'; ◂ = continues earlier' if hidden else ''}.</p></div>")
        months_html = f"<div class='stack'>{months_html}</div>"
    top_html = ""
    if members:
        body = ""
        for c in top[:MAP_CLUSTERS]:
            ds = cluster_dates[c]
            first_cell = f"<td>{ds[0]:%b} ’{ds[0]:%y}</td>" if ds else "<td>–</td>"
            dot = f"<i class='dot {'on' if c in active else 'off'}'></i>"  # filled: a sample in the last few months
            last_cell = f"<td class='last'>{dot}{ds[-1]:%b} ’{ds[-1]:%y}</td>" if ds else "<td>–</td>"
            sites = f"<td class='num'>{len({r[loc] for r in members[c]})}</td>" if loc else ""
            moved = "<td class='num'>–</td>"  # the column is always there; it only has values when there is an earlier report
            if growth:
                d = growth["change"][c]
                moved = ("<td><b class='new'>NEW</b></td>" if c in growth["new"] else f"<td class='num up'>+{d}</td>" if d > 0 else moved)
            body += (f"<tr><td class='key'>{key_of[c]}</td><td>{_e(name_of[c])}</td><td class='num'>{len(members[c])}</td>{moved}{sites}"
                     f"{first_cell}{last_cell}</tr>")
        site_head = "<th class='num'>Sites</th>" if loc else ""
        change_head = "<th class='num'>Change</th>"
        top_html = (f"<div><h2>Largest clusters</h2><table><tr><th></th><th>Label</th><th class='num'>Seqs</th>{change_head}{site_head}"
                    f"<th>First</th><th>Last sample</th></tr>{body}</table><p class='note'>Letters match the map"
                    f"{f'; ● = a sample in the last {active_months} months' if latest else ''}{'; active clusters only' if active_only else ''}.</p></div>")
    map_html = ""
    if members and links:
        drawn = _cluster_map(members, links, {s: palette[v] for s, v in value_of.items()}, key_of, top=MAP_CLUSTERS,
                             labelled=MAP_CLUSTERS, label_size=10, added=new_ids, height=160)
        if drawn:
            svg, _caption = drawn
            used = set(value_of.values())
            legend = "".join(f"<span><i style='background:{palette[v]}'></i>{_e('Other sites' if v == 'other' else v)}</span>"
                             for v in [*order, "other"] if v in used)
            map_html = (f"<div><h2>Cluster map</h2>{svg}<div class='legend'>{legend}</div><p class='note'>{len(key_of)} largest of "
                        f"{f'{len(top)} active clusters ({len(members)} in all)' if active_only else f'{len(members)} clusters'}{'' if loc else '; colour is subtype'}"
                        f"{'; circled = new' if growth else ''}.</p></div>")

    def columns(*blocks: str, widths: str = "1fr") -> str:
        blocks = tuple(b for b in blocks if b)
        template = f"repeat({len(blocks)},1fr)" if len(blocks) != 2 else widths
        return f"<div class='cols' style='grid-template-columns:{template}'>{''.join(blocks)}</div>" if blocks else ""

    sequence_block = columns(months_html, demo_html, widths="0.85fr 1.15fr")
    cluster_block = columns(location_html, f"<div class='stack'>{map_html}{top_html}</div>" if map_html or top_html else "", widths="1fr 1fr")

    # what changed in the registry
    if growth:
        bits = [f"{len(growth['new_ids'])} new sequences", f"{len(growth['new'])} new cluster(s)",
                f"{len(growth['grown'])} existing cluster(s) grew (+{sum(growth['change'][c] for c in growth['grown'])} sequences)"]
        since_text = f"the report of {_long_date(growth['date'])}" if growth["date"] else "the earlier run"
        since_html = f"<b>Since {since_text}:</b> {', '.join(bits)}."
    elif events_file.exists():
        text = (f"{kinds['created']} new cluster(s) named (first run in this registry, so every cluster is new)" if first_run else
                f"{kinds['created']} new, {joined} grown, {kinds['merged']} merged, {kinds['split']} split cluster(s)")
        since_html = f"<b>Since last report:</b> {text}."
    else:
        since_html = ("<b>Since last report:</b> not tracked. Run with a cluster registry (<code>--registry</code>) to report "
                 "new, grown, merged and split clusters.")

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{_e(title)}</title><style>{_css()}</style></head><body><div class="page">
<header>{_logo_html(logo)}<div><h1>{_e(title)}</h1><p>Hepatitis C virus genomic surveillance</p>{f'<p class="range">{range_text}</p>' if range_text else ''}</div>
<div class="meta"><b>{_long_date(report_date)}</b>Data to {_long_date(cutoff)}</div></header>
<section><div class="band">Sequence information</div><div class="funnel">{funnel_html}</div>{sequence_block}</section>
<section><div class="band">Cluster information</div><div class="kpis">{tiles_html}</div>{cluster_block}</section>
<section class="callout">{since_html}</section>
<footer><b>Method.</b> {_method_text(params)}Pipeline hcv-cluster {_e(info.get("pipeline_version", __version__))}.<br>
<b>Interpretation.</b> A cluster shows genetic similarity and is a lead for epidemiological follow-up, not proof of transmission. Counts are of sequences, not people. No patient identifiers are shown.</footer>
</div></body></html>"""
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    return out


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="hcv-report",
        description="Write a one-page A4 progress report (HTML) from a finished run's output folder. "
                    "Open it in a browser and print to PDF.",
    )
    p.add_argument("run_dir", type=Path, help="Output folder of a finished run")
    p.add_argument("--out", type=Path, default=None, help=f"Where to write the report (default: RUN_DIR/{REPORT_FILENAME})")
    p.add_argument("--title", default=DEFAULT_TITLE)
    p.add_argument("--logo", type=Path, default=None, help="PNG, JPG or SVG for the top of the page")
    p.add_argument("--location-column", default=None, help="Default: a column named location, site or setting")
    p.add_argument("--date-column", default=None, help="Collection date column (default: sample_collection_date, collection_date or DOC)")
    p.add_argument("--age-column", default=None, help="Default: a column named age")
    p.add_argument("--sex-column", default=None, help="Default: a column named sex or gender")
    p.add_argument("--diagnosis-column", default=None, help="Default: a column named how_diagnosed or reason_for_test")
    p.add_argument("--injecting-column", default=None, help="Default: a column named injecting_status")
    p.add_argument("--indigenous-column", default=None, help="Default: a column named indigenous_status")
    p.add_argument("--since", type=Path, default=None, metavar="PREV_RUN_DIR",
                   help="Output folder of the run behind the last report (same --registry): flags new and grown clusters")
    window = p.add_mutually_exclusive_group()
    window.add_argument("--chart-from", metavar="YYYY-MM-DD", type=datetime.date.fromisoformat, default=None,
                        help="First month of the two time charts (default: the whole study, at most the last 36 months)")
    window.add_argument("--chart-months", metavar="N", type=int, default=None,
                        help="Show the last N months in the two time charts (default: the whole study, at most 36)")
    p.add_argument("--active-months", metavar="N", type=int, default=6,
                   help="A cluster is active if it has a sample in the last N months of the data (default: 6)")
    p.add_argument("--active-only", action="store_true",
                   help="Limit the map, timelines and largest-clusters table to active clusters (default: largest of all time)")
    p.add_argument("--min-cell", type=int, default=5, help="Demographic groups with fewer sequences are pooled as 'Other' (default: 5)")
    p.add_argument("--report-date", metavar="YYYY-MM-DD", type=datetime.date.fromisoformat, default=None,
                   help="Date on the page (default: today)")
    a = p.parse_args(argv)
    if a.logo is not None and not a.logo.is_file():
        p.error(f"--logo: no such file: {a.logo}")
    try:
        written = build_report(
            a.run_dir, a.out, title=a.title, logo=a.logo, location_column=a.location_column, date_column=a.date_column,
            age_column=a.age_column, sex_column=a.sex_column, diagnosis_column=a.diagnosis_column,
            injecting_column=a.injecting_column, indigenous_column=a.indigenous_column, min_cell=a.min_cell, since=a.since,
            chart_from=a.chart_from, chart_months=a.chart_months, active_months=a.active_months, active_only=a.active_only,
            report_date=a.report_date,
        )
    except (ValueError, FileNotFoundError) as exc:
        p.error(str(exc))
    print(f"Report written to {written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
