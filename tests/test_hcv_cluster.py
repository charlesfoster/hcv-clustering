from pathlib import Path

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
