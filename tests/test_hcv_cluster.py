from pathlib import Path

import pytest

import hcv_cluster
import hcv_workflow


# ---------------------------------------------------------------------------
# resolve_threshold
# ---------------------------------------------------------------------------


def test_resolve_threshold_returns_user_value() -> None:
    assert hcv_cluster.resolve_threshold("core-e2", 0.030) == 0.030


def test_resolve_threshold_core_e2_nohvr1_default() -> None:
    assert hcv_cluster.resolve_threshold("core-e2-nohvr1", None) == 0.03


def test_resolve_threshold_core_e2_full() -> None:
    assert hcv_cluster.resolve_threshold("core-e2", None) == 0.045


def test_resolve_threshold_ns5b() -> None:
    assert hcv_cluster.resolve_threshold("ns5b", None) == 0.015


def test_resolve_threshold_preset_expands_before_lookup() -> None:
    assert hcv_cluster.resolve_threshold("structural", None) == hcv_cluster.resolve_threshold(
        "core-e2-nohvr1", None
    )


def test_resolve_threshold_unevidenced_region_returns_fallback() -> None:
    assert hcv_cluster.resolve_threshold("e1-e2", None) == hcv_cluster.FALLBACK_THRESHOLD


def test_resolve_threshold_case_insensitive() -> None:
    assert hcv_cluster.resolve_threshold("NS5B", None) == hcv_cluster.resolve_threshold("ns5b", None)


def test_region_threshold_is_evidence_based() -> None:
    assert hcv_cluster.region_threshold_is_evidence_based("core-e2-nohvr1") is True
    assert hcv_cluster.region_threshold_is_evidence_based("core-e2") is True
    assert hcv_cluster.region_threshold_is_evidence_based("ns5b") is True
    assert hcv_cluster.region_threshold_is_evidence_based("e1-e2") is False
    assert hcv_cluster.region_threshold_is_evidence_based("cds") is False


def test_resolve_threshold_snp_matches_lamoury_p_distance() -> None:
    # This pipeline's SNP distance is uncorrected p-distance -- the same metric
    # Lamoury et al. 2015 used, so their region values apply directly (not as an
    # approximation the way they're reused for TN93's core-e2/ns5b defaults).
    assert hcv_cluster.resolve_threshold("core-e2-nohvr1", None, distance="snp") == 0.03
    assert hcv_cluster.resolve_threshold("core-e2", None, distance="snp") == 0.045
    assert hcv_cluster.resolve_threshold("ns5b", None, distance="snp") == 0.015


def test_resolve_threshold_snp_user_override_still_wins() -> None:
    assert hcv_cluster.resolve_threshold("core-e2", 0.02, distance="snp") == 0.02


def test_resolve_threshold_snp_unevidenced_region_returns_snp_fallback() -> None:
    assert hcv_cluster.resolve_threshold(
        "e1-e2", None, distance="snp"
    ) == hcv_cluster.FALLBACK_THRESHOLD_SNP


def test_region_threshold_is_evidence_based_snp() -> None:
    assert hcv_cluster.region_threshold_is_evidence_based("core-e2-nohvr1", distance="snp") is True
    assert hcv_cluster.region_threshold_is_evidence_based("e1-e2", distance="snp") is False


# ---------------------------------------------------------------------------
# compute_snp_distances / compute_snp_distances_detailed
# ---------------------------------------------------------------------------

_ALIGNED_FASTA_TWO_SEQS = """\
>seq1
ACGT-N
>seq2
ACGG-N
"""

_ALIGNED_FASTA_IDENTICAL = """\
>seqA
ACGTACGT
>seqB
ACGTACGT
"""


def test_compute_snp_distances_detailed_counts(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)

    results = hcv_cluster.compute_snp_distances_detailed(fasta)

    assert len(results) == 1
    id_i, id_j, snp_count, comparable_sites, normalized = results[0]
    assert id_i == "seq1"
    assert id_j == "seq2"
    # positions: A/A (same), C/C (same), G/G (same), T/G (SNP), -/- (gap, excluded), N/N (N, excluded)
    assert snp_count == 1
    assert comparable_sites == 4
    assert abs(normalized - 1 / 4) < 1e-9


def test_compute_snp_distances_gaps_and_ns_excluded(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)

    results = hcv_cluster.compute_snp_distances_detailed(fasta)

    _, _, snp_count, comparable_sites, _ = results[0]
    # gap (-) and N positions must not count
    assert comparable_sites == 4


def test_compute_snp_distances_no_self_pairs(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)

    rows = hcv_cluster.compute_snp_distances(fasta)

    for row in rows:
        assert row.source != row.target


def test_compute_snp_distances_identical_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_IDENTICAL)

    rows = hcv_cluster.compute_snp_distances(fasta)

    assert len(rows) == 1
    assert rows[0].distance == 0.0


# ---------------------------------------------------------------------------
# write_snp_csv
# ---------------------------------------------------------------------------


def test_write_snp_csv_header(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)
    out_csv = tmp_path / "snp.csv"

    hcv_cluster.write_snp_csv(out_csv, fasta)

    lines = out_csv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source,target,distance,snp_count,comparable_sites"


def test_write_snp_csv_distance_values(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)
    out_csv = tmp_path / "snp.csv"

    rows = hcv_cluster.write_snp_csv(out_csv, fasta)

    assert len(rows) == 1
    assert abs(rows[0].distance - 1 / 4) < 1e-9


def test_write_snp_csv_parseable_by_parse_distance_table(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)
    out_csv = tmp_path / "snp.csv"

    written_rows = hcv_cluster.write_snp_csv(out_csv, fasta)
    parsed_rows = hcv_workflow.parse_distance_table(out_csv)

    assert len(parsed_rows) == len(written_rows)
    for parsed, written in zip(parsed_rows, written_rows):
        assert parsed.source == written.source
        assert parsed.target == written.target
        assert abs(parsed.distance - written.distance) < 1e-9


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def test_build_parser_advanced_suppressed() -> None:
    parser = hcv_cluster.build_parser(show_advanced=False)
    # Find the run subparser and check that --min-coverage help is SUPPRESS
    run_subparser = None
    for action in parser._subparsers._actions:  # type: ignore[attr-defined]
        if hasattr(action, "_name_parser_map"):
            run_subparser = action._name_parser_map.get("run")
    assert run_subparser is not None
    min_cov_action = next(
        (a for a in run_subparser._actions if "--min-coverage" in getattr(a, "option_strings", [])),
        None,
    )
    assert min_cov_action is not None
    import argparse as _argparse
    assert min_cov_action.help == _argparse.SUPPRESS


def test_build_parser_advanced_visible() -> None:
    parser = hcv_cluster.build_parser(show_advanced=True)
    run_subparser = None
    for action in parser._subparsers._actions:  # type: ignore[attr-defined]
        if hasattr(action, "_name_parser_map"):
            run_subparser = action._name_parser_map.get("run")
    assert run_subparser is not None
    min_cov_action = next(
        (a for a in run_subparser._actions if "--min-coverage" in getattr(a, "option_strings", [])),
        None,
    )
    assert min_cov_action is not None
    import argparse as _argparse
    assert min_cov_action.help != _argparse.SUPPRESS


def test_main_help_advanced_returns_zero() -> None:
    rc = hcv_cluster.main(["--help-advanced"])
    assert rc == 0


def test_main_empty_argv_returns_zero() -> None:
    rc = hcv_cluster.main([])
    assert rc == 0


# ---------------------------------------------------------------------------
# main subcommand routing
# ---------------------------------------------------------------------------


def test_main_unknown_subcommand_returns_nonzero() -> None:
    import io
    import contextlib

    with contextlib.suppress(SystemExit):
        rc = hcv_cluster.main(["unknown-subcommand"])
        assert rc != 0


# ---------------------------------------------------------------------------
# write_snp_links_csv / SNP counts
# ---------------------------------------------------------------------------


def _detailed(pairs):
    """(source, target, snp_count, comparable_sites, distance) tuples."""
    return [
        (source, target, snp_count, sites, snp_count / sites)
        for source, target, snp_count, sites in pairs
    ]


def test_write_snp_links_csv_carries_absolute_counts(tmp_path: Path) -> None:
    out_csv = tmp_path / "snp_links.csv"

    links = hcv_cluster.write_snp_links_csv(
        out_csv, _detailed([("a", "b", 30, 2000)]), threshold=0.03
    )

    lines = out_csv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source,target,distance,snp_count,comparable_sites"
    assert lines[1].startswith("a,b,0.015,30,2000")
    assert links[0]["snp_count"] == 30


def test_write_snp_links_csv_thresholds_on_p_distance_by_default(tmp_path: Path) -> None:
    # Same 30 SNPs, but over few enough sites that p-distance exceeds the threshold.
    links = hcv_cluster.write_snp_links_csv(
        tmp_path / "l.csv",
        _detailed([("a", "b", 30, 2000), ("c", "d", 30, 500)]),
        threshold=0.03,
    )

    assert [(link["source"], link["target"]) for link in links] == [("a", "b")]


def test_snp_count_threshold_links_on_raw_count(tmp_path: Path) -> None:
    """The pair p-distance rejects is linked when the user asks for '<=30 SNPs'."""
    links = hcv_cluster.write_snp_links_csv(
        tmp_path / "l.csv",
        _detailed([("a", "b", 30, 2000), ("c", "d", 30, 500), ("e", "f", 31, 5000)]),
        threshold=0.03,
        snp_count_threshold=30,
    )

    assert [(link["source"], link["target"]) for link in links] == [("a", "b"), ("c", "d")]


def test_snp_count_threshold_of_zero_links_only_identical_pairs(tmp_path: Path) -> None:
    links = hcv_cluster.write_snp_links_csv(
        tmp_path / "l.csv",
        _detailed([("a", "b", 0, 2000), ("c", "d", 1, 2000)]),
        threshold=0.03,
        snp_count_threshold=0,
    )

    assert [(link["source"], link["target"]) for link in links] == [("a", "b")]


def test_negative_snp_count_threshold_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        hcv_cluster.write_snp_links_csv(
            tmp_path / "l.csv", _detailed([("a", "b", 1, 100)]), threshold=0.03,
            snp_count_threshold=-1,
        )


def test_snp_links_are_parseable_for_clustering(tmp_path: Path) -> None:
    out_csv = tmp_path / "snp_links.csv"
    hcv_cluster.write_snp_links_csv(
        out_csv, _detailed([("a", "b", 3, 1000)]), threshold=0.03
    )

    parsed = hcv_workflow.parse_links(out_csv)

    assert [(link.source, link.target) for link in parsed] == [("a", "b")]


def test_snp_count_summary_reports_range_and_median() -> None:
    summary = hcv_cluster._snp_count_summary(
        [
            {"snp_count": 5, "comparable_sites": 2000},
            {"snp_count": 30, "comparable_sites": 2100},
            {"snp_count": 12, "comparable_sites": 1900},
        ]
    )

    assert "5-30 SNPs" in summary
    assert "median 12" in summary
    assert "1900-2100 comparable sites" in summary


def test_snp_count_summary_handles_no_links() -> None:
    assert hcv_cluster._snp_count_summary([]) == "no linked pairs"


def test_merge_link_tables_preserves_snp_columns(tmp_path: Path) -> None:
    first = tmp_path / "gt1.csv"
    second = tmp_path / "gt2.csv"
    hcv_cluster.write_snp_links_csv(first, _detailed([("a", "b", 3, 1000)]), threshold=0.03)
    hcv_cluster.write_snp_links_csv(second, _detailed([("c", "d", 7, 1200)]), threshold=0.03)
    merged = tmp_path / "merged.csv"

    hcv_workflow.merge_link_tables([("1a", first), ("1b", second)], merged)

    lines = merged.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source,target,distance,snp_count,comparable_sites"
    assert lines[1].endswith(",3,1000")
    assert lines[2].endswith(",7,1200")


def test_merge_link_tables_drops_snp_columns_when_not_in_every_input(tmp_path: Path) -> None:
    with_counts = tmp_path / "snp.csv"
    hcv_cluster.write_snp_links_csv(with_counts, _detailed([("a", "b", 3, 1000)]), threshold=0.03)
    without_counts = tmp_path / "tn93.csv"
    hcv_workflow.write_links(without_counts, [hcv_workflow.LinkRow("c", "d", 0.01)])
    merged = tmp_path / "merged.csv"

    hcv_workflow.merge_link_tables([("1a", with_counts), ("1b", without_counts)], merged)

    assert merged.read_text(encoding="utf-8").splitlines()[0] == "source,target,distance"
