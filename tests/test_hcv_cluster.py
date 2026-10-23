from pathlib import Path

import pytest

from hcv_cluster import cli, workflow


# ---------------------------------------------------------------------------
# resolve_threshold
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("region", "distance", "expected", "evidence_based"),
    [
        pytest.param("core-e2-nohvr1", "tn93", 0.03, True, id="tn93-core-e2-nohvr1"),
        pytest.param("core-e2", "tn93", 0.045, True, id="tn93-core-e2"),
        pytest.param("ns5b", "tn93", 0.015, True, id="tn93-ns5b"),
        pytest.param("NS5B", "tn93", 0.015, None, id="tn93-case-insensitive"),
        # E1-E2 defaults are explicit but provisional (not evidence-based).
        pytest.param("e1-e2", "tn93", 0.045, False, id="tn93-e1-e2-provisional"),
        pytest.param("e1-e2-nohvr1", "tn93", 0.03, False, id="tn93-e1-e2-nohvr1-provisional"),
        # Presets and equivalent expressions normalise before lookup.
        pytest.param("structural", "tn93", 0.03, None, id="tn93-preset-expands"),
        pytest.param("envelope", "tn93", 0.045, None, id="tn93-envelope"),
        pytest.param("e1+e2", "tn93", 0.045, None, id="tn93-e1+e2"),
        pytest.param("core-e1+e2", "tn93", 0.045, None, id="tn93-core-e1+e2"),
        pytest.param("core+e1+e2", "tn93", 0.045, None, id="tn93-core+e1+e2"),
        # The pipeline's SNP distance is uncorrected p-distance -- the metric
        # Lamoury et al. 2015 used -- so their region values apply directly.
        pytest.param("core-e2-nohvr1", "snp", 0.03, True, id="snp-core-e2-nohvr1"),
        pytest.param("core-e2", "snp", 0.045, None, id="snp-core-e2"),
        pytest.param("ns5b", "snp", 0.015, None, id="snp-ns5b"),
        pytest.param("e1-e2", "snp", 0.045, False, id="snp-e1-e2-provisional"),
        pytest.param("cds", "snp", None, False, id="snp-fallback"),
    ],
)
def test_resolve_threshold_defaults_and_evidence_flag(
    region: str, distance: str, expected: float | None, evidence_based: bool | None
) -> None:
    if expected is None:  # unevidenced SNP region -> SNP-specific fallback
        expected = cli.FALLBACK_THRESHOLD_SNP
    assert cli.resolve_threshold(region, None, distance=distance) == expected
    if evidence_based is not None:
        assert (
            cli.region_threshold_is_evidence_based(region, distance=distance)
            is evidence_based
        )


def test_resolve_threshold_unevidenced_region_returns_fallback() -> None:
    assert cli.resolve_threshold("cds", None) == cli.FALLBACK_THRESHOLD == 0.03  # not tied to e1-e2
    assert cli.region_threshold_is_evidence_based("cds") is False
    assert cli.region_threshold_is_explicit("cds") is False
    assert cli.region_threshold_is_explicit("e1-e2") is True


@pytest.mark.parametrize(
    ("region", "distance", "user_value"),
    [
        pytest.param("core-e2", "tn93", 0.030, id="tn93"),
        pytest.param("core-e2", "snp", 0.02, id="snp"),
    ],
)
def test_resolve_threshold_user_value_wins(region: str, distance: str, user_value: float) -> None:
    assert cli.resolve_threshold(region, user_value, distance=distance) == user_value


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

    results = cli.compute_snp_distances_detailed(fasta)

    assert len(results) == 1
    id_i, id_j, snp_count, comparable_sites, normalized = results[0]
    assert id_i == "seq1"
    assert id_j == "seq2"
    # positions: A/A (same), C/C (same), G/G (same), T/G (SNP), -/- (gap, excluded), N/N (N, excluded)
    assert snp_count == 1
    assert comparable_sites == 4
    assert abs(normalized - 1 / 4) < 1e-9


def test_compute_snp_distances_identical_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_IDENTICAL)

    rows = cli.compute_snp_distances(fasta)

    assert len(rows) == 1
    assert rows[0].distance == 0.0


# ---------------------------------------------------------------------------
# write_snp_csv
# ---------------------------------------------------------------------------


def test_write_snp_csv_parseable_by_parse_distance_table(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(_ALIGNED_FASTA_TWO_SEQS)
    out_csv = tmp_path / "snp.csv"

    written_rows = cli.write_snp_csv(out_csv, fasta)
    parsed_rows = workflow.parse_distance_table(out_csv)

    assert len(parsed_rows) == len(written_rows)
    for parsed, written in zip(parsed_rows, written_rows):
        assert parsed.source == written.source
        assert parsed.target == written.target
        assert abs(parsed.distance - written.distance) < 1e-9


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def test_run_parser_uses_e1_e2_and_tolerant_alignment_defaults() -> None:
    args = cli.build_parser().parse_args(["run", "--input", "samples.fasta"])
    assert args.region == "e1-e2"
    assert args.min_coverage == 0.7
    assert args.minimap2_preset == "asm20"


@pytest.mark.parametrize("argv", [["--help-advanced"], []], ids=["help-advanced", "empty-argv"])
def test_main_without_subcommand_returns_zero(argv: list[str]) -> None:
    assert cli.main(argv) == 0


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

    links = cli.write_snp_links_csv(
        out_csv, _detailed([("a", "b", 30, 2000)]), threshold=0.03
    )

    lines = out_csv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source,target,distance,snp_count,comparable_sites"
    assert lines[1].startswith("a,b,0.015,30,2000")
    assert links[0]["snp_count"] == 30


def test_write_snp_links_csv_thresholds_on_p_distance_by_default(tmp_path: Path) -> None:
    # Same 30 SNPs, but over few enough sites that p-distance exceeds the threshold.
    links = cli.write_snp_links_csv(
        tmp_path / "l.csv",
        _detailed([("a", "b", 30, 2000), ("c", "d", 30, 500)]),
        threshold=0.03,
    )

    assert [(link["source"], link["target"]) for link in links] == [("a", "b")]


@pytest.mark.parametrize(
    ("pairs", "snp_count_threshold", "expected"),
    [
        pytest.param(
            [("a", "b", 30, 2000), ("c", "d", 30, 500), ("e", "f", 31, 5000)],
            30,
            [("a", "b"), ("c", "d")],
            id="raw-count-overrides-p-distance",
        ),
        pytest.param(
            [("a", "b", 0, 2000), ("c", "d", 1, 2000)],
            0,
            [("a", "b")],
            id="zero-links-only-identical-pairs",
        ),
    ],
)
def test_snp_count_threshold_links_on_raw_count(
    tmp_path: Path, pairs: list[tuple], snp_count_threshold: int, expected: list[tuple]
) -> None:
    """A pair p-distance rejects is linked when the user asks for '<=30 SNPs'."""
    links = cli.write_snp_links_csv(
        tmp_path / "l.csv",
        _detailed(pairs),
        threshold=0.03,
        snp_count_threshold=snp_count_threshold,
    )

    assert [(link["source"], link["target"]) for link in links] == expected


def test_negative_snp_count_threshold_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        cli.write_snp_links_csv(
            tmp_path / "l.csv", _detailed([("a", "b", 1, 100)]), threshold=0.03,
            snp_count_threshold=-1,
        )


def test_snp_links_are_parseable_for_clustering(tmp_path: Path) -> None:
    out_csv = tmp_path / "snp_links.csv"
    cli.write_snp_links_csv(
        out_csv, _detailed([("a", "b", 3, 1000)]), threshold=0.03
    )

    parsed = workflow.parse_links(out_csv)

    assert [(link.source, link.target) for link in parsed] == [("a", "b")]


def test_merge_link_tables_preserves_snp_columns(tmp_path: Path) -> None:
    first = tmp_path / "gt1.csv"
    second = tmp_path / "gt2.csv"
    cli.write_snp_links_csv(first, _detailed([("a", "b", 3, 1000)]), threshold=0.03)
    cli.write_snp_links_csv(second, _detailed([("c", "d", 7, 1200)]), threshold=0.03)
    merged = tmp_path / "merged.csv"

    workflow.merge_link_tables([("1a", first), ("1b", second)], merged)

    lines = merged.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "source,target,distance,snp_count,comparable_sites"
    assert lines[1].endswith(",3,1000")
    assert lines[2].endswith(",7,1200")


def test_merge_link_tables_drops_snp_columns_when_not_in_every_input(tmp_path: Path) -> None:
    with_counts = tmp_path / "snp.csv"
    cli.write_snp_links_csv(with_counts, _detailed([("a", "b", 3, 1000)]), threshold=0.03)
    without_counts = tmp_path / "tn93.csv"
    workflow.write_links(without_counts, [workflow.LinkRow("c", "d", 0.01)])
    merged = tmp_path / "merged.csv"

    workflow.merge_link_tables([("1a", with_counts), ("1b", without_counts)], merged)

    assert merged.read_text(encoding="utf-8").splitlines()[0] == "source,target,distance"
