import csv
import datetime
import json
from collections import Counter
import subprocess
import sys
from pathlib import Path

import pytest

from hcv_cluster import report


@pytest.fixture(autouse=True)
def _no_logo_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(report, "DEFAULT_LOGO", None)  # keeps the test pages small; the logo has its own test


def _write(path: Path, fieldnames: list[str], rows: list[dict]) -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _run_dir(tmp_path: Path, *, metadata: bool = True, registry: bool = True, name: str = "run",
             without: tuple[str, ...] = ()) -> Path:
    """A finished run: ten sequences received, two fail, eight are analysed (``without`` leaves some out, as an earlier run would)."""
    people = [  # sample, cluster, location, collected, age, sex, injecting, indigenous
        ("SMP-A1", "1a_C0001", "Clinic X", "5-Mar-24", "22", "M", "Last month", "Aboriginal"),
        ("SMP-A2", "1a_C0001", "Clinic X", "9-Mar-24", "31", "M", "Last 6 months", "Neither"),
        ("SMP-A3", "1a_C0001", "Clinic Y", "2-Apr-24", "", "F", "Unknown", "Torres Strait Islander"),
        ("SMP-B1", "1a_C0002", "Clinic X", "7-May-24", "40", "M", "Never", "Neither"),
        ("SMP-B2", "1a_C0002", "Clinic X", "9-May-24", "58", "M", "Last month", "Aboriginal"),
        ("SMP-S1", "1a_C0003", "Clinic Z", "1-Jun-24", "35", "F", "Never", "Neither"),
        ("SMP-S2", "1a_C0004", "Clinic Z", "2-Jun-24", "36", "F", "Last month", "Aboriginal"),
        ("SMP-S3", "1a_C0005", "Clinic W", "3-Jun-24", "29", "M", "Last 2 years", "Neither"),
    ]
    people = [p for p in people if p[0] not in without]
    size = Counter(p[1] for p in people)
    fields = ["sample_id", "genotype", "cluster_id", "cluster_label", "cluster_size"]
    rows = [{"sample_id": s, "genotype": "1a", "cluster_id": c, "cluster_label": f"HCV-{c[-4:]}" if size[c] > 1 else "",
             "cluster_size": size[c], "Location": loc, "sample_collection_date": d, "Age": age, "Gender": sex,
             "injecting_status": inj, "indigenous_status": ind}
            for s, c, loc, d, age, sex, inj, ind in people]
    if metadata:
        fields += ["Location", "sample_collection_date", "Age", "Gender", "injecting_status", "indigenous_status"]
    else:
        rows = [{k: r[k] for k in fields} for r in rows]
    run = tmp_path / name
    run.mkdir()
    _write(run / "metadata.surviving_samples.csv", fields, rows)
    _write(run / "sample_audit.csv", ["sample_id", "overall_status", "failed_stage"],
           [{"sample_id": p[0], "overall_status": "pass", "failed_stage": ""} for p in people]
           + [{"sample_id": "SMP-F1", "overall_status": "fail", "failed_stage": "genotyping"},
              {"sample_id": "SMP-F2", "overall_status": "fail", "failed_stage": "region_qc"}])
    _write(run / "links.csv", ["source", "target", "distance"], [
        {"source": "SMP-A1", "target": "SMP-A2", "distance": "0.01"}, {"source": "SMP-A3", "target": "SMP-A2", "distance": "0.02"},
        {"source": "SMP-B1", "target": "SMP-B2", "distance": "0.01"}])
    (run / "run_info.json").write_text(json.dumps({"created": "2026-10-28T01:00:00", "parameters": {"region": "e1-e2"}}))
    if registry:
        _write(run / "label_changes.csv", ["event", "label"], [{"event": "created", "label": "HCV-C0001"}, {"event": "created", "label": "HCV-C0002"}])
    return run


def test_report_summarises_a_run_and_shows_no_identifiers(tmp_path: Path) -> None:
    run = _run_dir(tmp_path)
    page = report.build_report(run, min_cell=3, title="Test update").read_text()
    assert (run / "report.html").exists()  # the default place
    assert "<b>10</b><span>sequenced" in page and "<b>8</b><span>passed QC" in page and "<b>5</b><span>in a cluster" in page
    assert "Sequence information" in page and "Cluster information" in page and "region" not in page.split("<footer>")[0].lower()
    assert "<b>2</b><span>clusters (2+ sequences)" in page and "<b>1</b><span>across 2+ locations" in page  # cluster A: Clinic X and Y
    assert all(f"Clinic {c}" in page for c in "XYZW") and "Other sites" not in page  # every location has its own row
    assert "Who is in the data" in page and "median 35, IQR" in page and "Male" in page
    assert "Samples collected Mar 2024 to Jun 2024" in page and "time charts from" not in page  # the header says what the data cover
    assert "Injecting drug use" in page and "Last 6 months" in page and "Indigenous status" in page
    assert "Aboriginal and/or TSI" in page and "Torres Strait Islander" not in page  # any Aboriginal or TSI answer is pooled
    assert "Sequences by collection month" in page and "Cluster timelines" in page and "Cluster map" in page
    assert "<td class='key'>a</td><td>HCV-0001</td>" in page and "<td class='key'>b</td><td>HCV-0002</td>" in page  # the map's keys
    assert "<b>–</b><span>new clusters since last report" in page and "first run in this registry" in page  # a first run has nothing to compare with
    assert "SMP-" not in page  # counts only, never a sample identifier
    assert 'onclick="window.print()"' in page and "print-color-adjust: exact" in page and "<script" not in page  # a Save as PDF button, no other script


def test_report_leaves_out_what_the_run_has_no_data_for(tmp_path: Path) -> None:
    run = _run_dir(tmp_path, metadata=False, registry=False)
    page = report.build_report(run).read_text()
    for absent in ("By location", "Who is in the data", "collection month", "across 2+ locations"):
        assert absent not in page
    assert "colour is subtype" in page and "not tracked" in page
    with pytest.raises(ValueError, match="'Ward' is not a column"):
        report.build_report(run, location_column="Ward")
    with pytest.raises(FileNotFoundError, match="is it a finished run"):
        report.build_report(tmp_path)


def test_report_flags_growth_since_an_earlier_run(tmp_path: Path) -> None:
    now = _run_dir(tmp_path)
    earlier = _run_dir(tmp_path, name="earlier", without=("SMP-A3", "SMP-B1", "SMP-B2"))  # cluster A had two members; B did not exist
    page = report.build_report(now, since=earlier).read_text()
    assert "<b>8 <em>+3</em></b><span>passed QC" in page  # eight sequences, three of them new
    assert "<b>1</b><span>new clusters since last report" in page and "<b>1</b><span>existing clusters that grew" in page
    assert "<b class='new'>NEW</b>" in page and "<td class='num up'>+1</td>" in page and "3 new sequences" in page
    assert page.count('r="5.2" fill="none"') == 6  # the three new sequences are circled on the map and on the timelines
    assert "<em>" not in report.build_report(now).read_text()  # without --since nothing is flagged
    unlabelled = _run_dir(tmp_path, name="unlabelled", without=("SMP-A1", "SMP-A2", "SMP-B1", "SMP-B2"))
    with pytest.raises(ValueError, match="persistent cluster labels"):
        report.build_report(unlabelled, since=earlier)


def test_report_finds_the_collection_date_column_by_name(tmp_path: Path) -> None:
    run = _run_dir(tmp_path)
    table = run / "metadata.surviving_samples.csv"
    rows = list(csv.DictReader(table.open()))
    _write(table, [("DOC" if f == "sample_collection_date" else f) for f in rows[0]],
           [{("DOC" if k == "sample_collection_date" else k): v for k, v in r.items()} for r in rows])
    assert "Sequences by collection month" in report.build_report(run).read_text()  # DOC is recognised
    assert "Sequences by collection month" not in report.build_report(run, date_column="Age").read_text()  # a named column wins (no dates in it)


def test_time_charts_show_at_most_36_months_unless_told_otherwise() -> None:
    first, latest = (2020, 1), (2026, 9)
    assert report._window_start(first, latest, None, None) == (2023, 10)  # the automatic cap: 36 months ending in Sep 2026
    assert report._window_start((2025, 6), latest, None, None) == (2025, 6)  # a short study is shown whole
    assert report._window_start(first, latest, None, 12) == (2025, 10)
    assert report._window_start(first, latest, datetime.date(2022, 3, 9), None) == (2022, 3)  # a date sets the month it falls in
    assert report._window_start(first, latest, datetime.date(2019, 1, 1), None) == first  # never before the first sample


def test_report_windows_the_charts_and_flags_active_clusters(tmp_path: Path) -> None:
    run = _run_dir(tmp_path)  # samples are collected 5 Mar to 3 Jun 2024; cluster A last on 2 Apr and B on 9 May
    windowed = report.build_report(run, chart_from=datetime.date(2024, 5, 1)).read_text()
    assert "3 earlier not shown" in windowed and "◂ = continues earlier" in windowed
    assert "Samples collected Mar 2024 to Jun 2024 · time charts from May 2024" in windowed  # the header keeps the two ranges apart  # the Mar/Apr sequences are off the chart
    assert "earlier not shown" not in report.build_report(run).read_text()
    page = report.build_report(run, active_months=1).read_text()  # active: a sample since 3 May
    assert "<i class='dot off'></i>Apr ’24" in page and "<i class='dot on'></i>May ’24" in page  # A is inactive, B active
    only = report.build_report(run, active_months=1, active_only=True).read_text()
    assert "<td class='key'>a</td><td>HCV-0002</td>" in only and "HCV-0001" not in only and "active clusters only" in only
    with pytest.raises(ValueError, match="needs collection dates"):
        report.build_report(_run_dir(tmp_path, metadata=False, name="nodates"), active_only=True)


def test_report_remembers_itself_in_the_registry_and_shows_growth_next_time(tmp_path: Path) -> None:
    registry = tmp_path / "reg"
    earlier = _run_dir(tmp_path, name="earlier", without=("SMP-A3", "SMP-B1", "SMP-B2"))
    now = _run_dir(tmp_path)
    assert "<em>" not in report.build_report(earlier, registry=registry).read_text()  # nothing to compare with yet
    assert not (registry / "last_report.json").exists()  # and not saved unless asked
    report.build_report(earlier, registry=registry, save_snapshot=True, report_date=datetime.date(2026, 7, 14))
    page = report.build_report(now, registry=registry).read_text()  # no --since: the saved report is the comparison
    assert "<b>8 <em>+3</em></b><span>passed QC" in page and "Since the report of 14 July 2026" in page
    table = now / "metadata.surviving_samples.csv"  # a run without persistent labels must not replace the saved report
    rows = list(csv.DictReader(table.open()))
    _write(table, list(rows[0]), [{**r, "cluster_label": ""} for r in rows])
    report.build_report(now, registry=registry, save_snapshot=True, report_date=datetime.date(2026, 11, 3))
    assert report.read_snapshot(registry)["date"] == "2026-07-14"


def test_report_has_the_h2seq_logo_and_title_unless_told_otherwise(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    logo = tmp_path / "logo.png"
    logo.write_bytes(b"\x89PNG fake")
    monkeypatch.setattr(report, "DEFAULT_LOGO", logo)
    run = _run_dir(tmp_path)
    page = report.build_report(run).read_text()
    assert "<img" in page and "H2Seq HCV sequencing" in page
    assert "<img" not in report.build_report(run, no_logo=True).read_text()


def test_pdf_needs_a_chrome_and_says_so_when_there_is_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fake = tmp_path / "chrome"
    fake.write_text("")
    monkeypatch.setenv("HCV_CHROME", str(fake))
    assert report.find_chrome() == str(fake)
    monkeypatch.setattr(report, "find_chrome", lambda: None)
    assert report.write_pdf(tmp_path / "r.html", tmp_path / "r.pdf") is None


def test_report_is_a_standalone_command(tmp_path: Path) -> None:
    run = _run_dir(tmp_path)
    out = tmp_path / "elsewhere" / "update.html"
    assert report.main([str(run), "--out", str(out), "--report-date", "2026-11-02"]) == 0
    assert "2 November 2026" in out.read_text()
    with pytest.raises(SystemExit):  # a bad date is refused rather than written
        report.main([str(run), "--report-date", "soon"])


def test_cluster_map_layout_does_not_depend_on_the_python_process() -> None:
    code = (
        "from hcv_cluster import viz\n"
        "rows = [{'sample_id': s, 'cluster_id': 'c'} for s in 'abcdefgh']\n"
        "edges = [{'source': a, 'target': b, 'distance': '0.01'} for a, b in "
        "[('a','b'),('b','c'),('c','d'),('d','e'),('e','f'),('a','h'),('g','h')]]\n"
        "p = viz.compute_network_layout(rows, edges, mode='components')\n"
        "print(round(p['a'][0], 4), round(p['h'][1], 4))\n"
    )
    outputs = {
        subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                       env={"PYTHONHASHSEED": str(seed), "PATH": ""}).stdout
        for seed in (1, 2, 3, 4)
    }
    assert len(outputs) == 1


def test_sequencing_status_adds_only_samples_that_never_reached_the_input(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    run = _run_dir(tmp_path)  # ten samples in the input
    sheet = _write(tmp_path / "status.csv", ["sample_id", "status"], [
        {"sample_id": "SMP-A1", "status": "sequenced"},     # in the input: the FASTA decides, the list changes nothing
        {"sample_id": "smp-f1", "status": "sequenced"},     # the input sample SMP-F1 in another letter case: not an extra
        {"sample_id": "X1_1a", "status": "not_attempted"}, {"sample_id": "X1_3a_TR2", "status": "Not attempted"},  # tags dropped: one sample
        {"sample_id": "X2", "status": "not-attempted"}, {"sample_id": "SMP-A2_1a", "status": "sequenced"},  # the tag is dropped before matching the input
        {"sample_id": "X3", "status": "pending"}, {"sample_id": "X4", "status": "Fail"},  # an old word, read as sequenced: X4 belongs in the FASTA
        {"sample_id": "X5", "status": "sequenced"},         # says sequenced, but is not in the input: counted, and warned about
        {"sample_id": "", "status": "Fail"},
    ])
    page = report.build_report(run, min_cell=3, sequencing_status=sheet).read_text()
    assert "<b>15</b><span>samples collected" in page and "<b>10</b><span>sequenced" in page  # 10 + the five that are not in the input
    warnings = capsys.readouterr().err
    assert "X4" in warnings and "X5" in warnings and "letter case" in warnings
    plain = report.build_report(run, min_cell=3).read_text()
    assert "<b>10</b><span>samples collected" in plain  # no sheet: as before
    with pytest.raises(ValueError, match="'sample_id' and 'status'"):
        report.read_sequencing_status(_write(tmp_path / "bad.csv", ["id", "state"], [{"id": "a", "state": "b"}]))
