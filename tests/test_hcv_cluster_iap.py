import csv
import datetime
from pathlib import Path

import pytest

from hcv_cluster import cli, iap, stages


def test_time_difference_and_rate_follow_the_documented_example() -> None:
    first, second = iap.parse_date("1/7/1995"), iap.parse_date("1/12/1997")
    assert f"{iap.years_between(first, second):.6f}" == "2.416667"
    assert iap.years_between(second, first) == iap.years_between(first, second)
    assert f"{24 / iap.years_between(first, second):.6f}" == "9.931034"
    assert iap.parse_date("1997-12-01") == second  # ISO dates are read the same way


@pytest.mark.parametrize(
    "text, expected",
    [("5-Mar-24", datetime.date(2024, 3, 5)), ("29-Feb-24", datetime.date(2024, 2, 29)),
     ("19-feb-1995", datetime.date(1995, 2, 19)), ("1-Sept-98", datetime.date(1998, 9, 1)),
     ("12 DECEMBER 07", datetime.date(2007, 12, 12))],
)
def test_day_month_name_dates_are_read(text: str, expected: datetime.date) -> None:
    assert iap.parse_date(text) == expected


@pytest.mark.parametrize("text", ["5-Foo-24", "5-Ju-24", "Mar-24", "5-Mar-2"])
def test_malformed_month_name_dates_are_rejected(text: str) -> None:
    with pytest.raises(ValueError, match="is not d/m/yyyy"):
        iap.parse_date(text)


@pytest.mark.parametrize(
    "differences, years, expected",
    # Reference values from an independent sum of negative-binomial probabilities (stdlib only).
    [(5, 0, 0.544339), (10, 0, 0.334735), (24, 29 / 12, 0.296720), (10, 29 / 12, 0.834518),
     (54, 10, 0.454245), (100, 10, 0.000064)],
)
def test_relatedness_index_matches_the_negative_binomial_model(differences: int, years: float, expected: float) -> None:
    assert iap.relatedness_index([differences], [years])[0] == pytest.approx(expected, abs=1e-6)


def test_relatedness_index_is_undefined_beyond_22_years() -> None:
    inside, outside = iap.relatedness_index([50, 50], [22.0, 22.3])
    assert 0 < inside < 1
    assert outside != outside  # NaN: the model's probability parameter reaches 1 at about 22.2 years


def test_collection_dates_reject_bad_input() -> None:
    rows = [
        {"sample_id": "a", "sample_collection_date": "3/4/2001"},
        {"sample_id": "b", "sample_collection_date": "(missing)"},
        {"sample_id": "c", "sample_collection_date": "April 2001"},
        {"sample_id": "d", "sample_collection_date": "31/2/2001"},
    ]
    with pytest.raises(ValueError, match=r"c: 'April 2001'\n- d: '31/2/2001'"):
        iap.collection_dates(rows)
    with pytest.raises(ValueError, match="'seen_on' is not a metadata column"):
        iap.collection_dates(rows[:2], "seen_on")
    with pytest.raises(ValueError, match="needs --metadata"):
        iap.collection_dates([])
    assert iap.collection_dates(rows[:2]) == {"a": datetime.date(2001, 4, 3)}


def test_iap_pairs_are_written_per_subtype_and_merged(tmp_path: Path) -> None:
    base = "ACGTACGTAC" * 3
    fasta_1a = tmp_path / "1a.fasta"
    fasta_1a.write_text(
        f">a\n{base}\n"
        f">b\nT{base[1:]}\n"            # 1 SNP from a
        f">c\nT{base[1:10]}---{base[13:]}\n"  # same SNP as b plus a 3-nt deletion
        f">d__v1\n{base[:-1]}G\n"       # variant of sample d: dated through its base id
        f">undated\n{base}\n"
    )
    fasta_1b = tmp_path / "1b.fasta"
    fasta_1b.write_text(f">e\n{base}\n>f\nT{base[1:]}\n>g\n{base}\n")
    (tmp_path / "variants.csv").write_text(
        "sequence_id,base_sample_id,variant,n_variants,input_record,completeness,genotype,region_coverage,status,detail\n"
        "d__v1,d,v1,2,4,30,1a,,kept,\n"
    )
    (tmp_path / "metadata.csv").write_text(
        "sample_id,seen_on\n"
        "a,1/7/1995\nb,1/12/1997\nc,1/12/1997\nd,1/7/1995\ne,1/7/1995\nf,1/8/1995\ng,1/7/2020\nundated,\n"
    )

    tables = {}
    for subtype, fasta in (("1b", fasta_1b), ("1a", fasta_1a)):
        tables[subtype] = tmp_path / f"{subtype}.iap_pairs.csv"
        assert stages.main([
            "iap-pairs", "--input", str(fasta), "--genotype", subtype,
            "--metadata", str(tmp_path / "metadata.csv"), "--date-column", "seen_on",
            "--variants", str(tmp_path / "variants.csv"), "--out", str(tables[subtype]),
        ]) == 0
    merged = tmp_path / "iap_pairs.csv"
    assert stages.main(["iap-merge", "--pairs", *(f"{k}={v}" for k, v in tables.items()), "--out", str(merged)]) == 0

    rows_1a = list(csv.DictReader(tables["1a"].open()))
    rows_1b = list(csv.DictReader(tables["1b"].open()))
    assert list(rows_1a[0]) == [
        "seqid_1", "seqid_2", "raw_dist", "time_diff_years", "snp_rate", "subtype", "relatedness_index", "note",
    ]
    assert {r["subtype"] for r in rows_1a} == {"1a"} and [(r["seqid_1"], r["seqid_2"]) for r in rows_1b] == [("e", "f"), ("e", "g"), ("f", "g")]
    assert list(csv.DictReader(merged.open())) == rows_1a + rows_1b  # subtype order, no cross-subtype pairs
    scored = [float(r.pop("relatedness_index")) for r in rows_1a]  # every pair here is within the model's range
    assert all(0 <= index <= 1 for index in scored)
    # 25 years apart is beyond the model: no index, and the note says why.
    assert [(r["seqid_2"], r["relatedness_index"], r["note"]) for r in rows_1b if r["seqid_2"] == "g"] == [
        ("g", "", "too_far_apart"), ("g", "", "too_far_apart"),
    ]
    rows = rows_1a
    by_pair = {(r["seqid_1"], r["seqid_2"]): r for r in rows}
    assert set(by_pair) == {("a", "b"), ("a", "c"), ("a", "d__v1"), ("b", "c"), ("b", "d__v1"), ("c", "d__v1")}
    assert by_pair[("a", "b")] == {
        "seqid_1": "a", "seqid_2": "b", "raw_dist": "1", "time_diff_years": "2.416667", "snp_rate": "0.413793",
        "subtype": "1a", "note": "",
    }
    assert by_pair[("a", "c")]["raw_dist"] == "1"  # the deletion adds nothing to the count
    assert by_pair[("b", "c")] == {
        "seqid_1": "b", "seqid_2": "c", "raw_dist": "0", "time_diff_years": "0.000000", "snp_rate": "",
        "subtype": "1a", "note": "same_day",
    }
    assert by_pair[("a", "d__v1")]["time_diff_years"] == "0.000000"  # collected the same day as a


@pytest.mark.parametrize("region", ["e1-e2", "E1+E2", "envelope"])
def test_e1_e2_is_the_only_calibrated_region(region: str) -> None:
    iap.check_region(region)


@pytest.mark.parametrize(
    "region, strategy",
    [("core-e2", "fixed"), ("e1-e2-nohvr1", "fixed"), ("e1", "fixed"), ("e1-e2", "max-usable")],
)
def test_other_regions_abort_before_any_work(region: str, strategy: str, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="calibrated for the full e1-e2 region"):
        iap.check_region(region, strategy)

    (tmp_path / "metadata.csv").write_text("sample_id,sample_collection_date\na,1/7/1995\n")
    assert stages.main(["iap-check", "--metadata", str(tmp_path / "metadata.csv"), "--region", region,
                        "--region-strategy", strategy]) == 1
    assert cli.main(["run", "-i", str(tmp_path / "in.fasta"), "-o", str(tmp_path / "out"), "--iap-model",
                     "--metadata", str(tmp_path / "metadata.csv"), "--region", region,
                     "--region-strategy", strategy]) == 1
    assert not (tmp_path / "out").exists()


def test_iap_clusters_are_built_from_the_index_alone_and_compared_with_the_main_ones(tmp_path: Path) -> None:
    (tmp_path / "clusters.csv").write_text(
        "sample_id,genotype,cluster_id,source_cluster_id,cluster_size,cluster_label\n"
        + "".join(f"{s},1a,1a_{c},{c},{n},HCV-{c}\n" for s, c, n in [
            ("a", "C0001", 3), ("b", "C0001", 3), ("c", "C0001", 3), ("d", "C0002", 2), ("e", "C0002", 2),
            ("f", "C0003", 2), ("g", "C0003", 2), ("h", "C0004", 1), ("i", "C0005", 1), ("j", "C0006", 1)])
    )
    (tmp_path / "links.csv").write_text("source,target,distance\na,b,0.01\nb,c,0.02\nd,e,0.01\ng,f,0.03\n")
    (tmp_path / "iap_pairs.csv").write_text(
        "seqid_1,seqid_2,raw_dist,time_diff_years,snp_rate,subtype,relatedness_index,note\n"
        "a,b,3,1.000000,3.000000,1a,0.5,\n"
        "c,b,30,1.000000,30.000000,1a,0.01,\n"  # a main link the model does not support (stored the other way round)
        "d,e,5,1.000000,5.000000,1a,0.3,\n"
        "h,i,4,1.000000,4.000000,1a,0.4,\n"  # a link only the model makes
        "h,j,4,1.000000,4.000000,1a,0.001,\n"
        # f and g have no collection date, so they are not in this table
    )
    assert stages.main(["iap-support", "--clusters", str(tmp_path / "clusters.csv"), "--links", str(tmp_path / "links.csv"),
                        "--pairs", str(tmp_path / "iap_pairs.csv"), "--alpha", "0.05", "--outdir", str(tmp_path),
                        "--tables-dir", str(tmp_path / "tables")]) == 0

    links = {(r["source"], r["target"]): r for r in csv.DictReader((tmp_path / "links.iap_information.csv").open())}
    assert {k: (r["relatedness_index"], r["supported"], r["note"]) for k, r in links.items()} == {
        ("a", "b"): ("0.5", "TRUE", ""), ("b", "c"): ("0.01", "FALSE", ""),
        ("d", "e"): ("0.3", "TRUE", ""), ("g", "f"): ("", "FALSE", "no_collection_date"),
    }
    assert links[("b", "c")]["raw_dist"] == "30" and links[("a", "b")]["distance"] == "0.01"

    edges = {(r["source"], r["target"]): r for r in csv.DictReader((tmp_path / "links.iap.csv").open())}
    assert {k: (r["distance"], r["in_main"]) for k, r in edges.items()} == {
        ("a", "b"): ("0.01", "TRUE"), ("d", "e"): ("0.01", "TRUE"), ("h", "i"): ("", "FALSE"),
    }

    clusters = {r["sample_id"]: r for r in csv.DictReader((tmp_path / "clusters.iap.csv").open())}
    assert {s: (r["cluster_id"], r["cluster_size"]) for s, r in clusters.items()} == {
        "a": ("1a_IAP0001", "2"), "b": ("1a_IAP0001", "2"), "c": ("1a_IAP0002", "1"), "d": ("1a_IAP0003", "2"),
        "e": ("1a_IAP0003", "2"), "f": ("1a_IAP0004", "1"), "g": ("1a_IAP0005", "1"), "h": ("1a_IAP0006", "2"),
        "i": ("1a_IAP0006", "2"), "j": ("1a_IAP0007", "1"),
    }
    assert "base_sample_id" in clusters["a"]  # annotated like the main table

    comparison = {(r["set"], r["cluster_id"]): r for r in csv.DictReader((tmp_path / "clusters.comparison.csv").open())}
    assert {k: (r["relationship"], r["concordant"], r["counterparts"], r["unclustered_in_other"]) for k, r in comparison.items()} == {
        ("main", "1a_C0001"): ("split", "FALSE", "1a_IAP0001(2)", "1"),  # c is left out by the model
        ("main", "1a_C0002"): ("identical", "TRUE", "1a_IAP0003(2)", "0"),
        ("main", "1a_C0003"): ("not_clustered_in_other", "FALSE", "", "2"),  # f and g have no dates
        ("iap", "1a_IAP0001"): ("contained", "FALSE", "1a_C0001(2)", "0"),
        ("iap", "1a_IAP0003"): ("identical", "TRUE", "1a_C0002(2)", "0"),
        ("iap", "1a_IAP0006"): ("not_clustered_in_other", "FALSE", "", "2"),  # only the model links h and i
    }

    # Per-genotype tables, in the form the network plots read.
    assert [r["cluster_id"] for r in csv.DictReader((tmp_path / "tables" / "1a.iap.clusters.csv").open())][:2] == ["IAP0001", "IAP0001"]
    assert len(list(csv.DictReader((tmp_path / "tables" / "1a.iap.links.csv").open()))) == 3

    assert stages.main(["iap-support", "--clusters", str(tmp_path / "clusters.csv"), "--links", str(tmp_path / "links.csv"),
                        "--pairs", str(tmp_path / "iap_pairs.csv"), "--alpha", "1.5", "--outdir", str(tmp_path)]) == 1


def test_iap_labels_have_their_own_prefix_and_numbering_in_the_same_registry(tmp_path: Path) -> None:
    (tmp_path / "clusters.csv").write_text(
        "sample_id,genotype,cluster_id,source_cluster_id,cluster_size\n"
        + "".join(f"{s},1a,1a_{c},{c},2\n" for s, c in [("a", "C0001"), ("b", "C0001"), ("c", "C0002"), ("d", "C0002")])
    )
    (tmp_path / "links.csv").write_text("source,target,distance\na,b,0.01\nc,d,0.01\n")
    (tmp_path / "iap_pairs.csv").write_text(
        "seqid_1,seqid_2,raw_dist,time_diff_years,snp_rate,subtype,relatedness_index,note\n"
        "a,b,3,1.000000,3.000000,1a,0.5,\nc,d,3,1.000000,3.000000,1a,0.5,\n"
    )
    registry = tmp_path / "registry"
    common = ["--registry", str(registry), "--update-registry", "--label-date", "2026-10-23"]
    assert stages.main(["label", "--clusters", str(tmp_path / "clusters.csv"), "--out-clusters", str(tmp_path / "main.csv"),
                        "--out-labels", str(tmp_path / "l.csv"), "--out-changes", str(tmp_path / "c.csv"), *common]) == 0
    assert stages.main(["iap-support", "--clusters", str(tmp_path / "main.csv"), "--links", str(tmp_path / "links.csv"),
                        "--pairs", str(tmp_path / "iap_pairs.csv"), "--outdir", str(tmp_path), *common]) == 0

    main = {r["cluster_label"] for r in csv.DictReader((tmp_path / "main.csv").open())}
    iap_labels = {r["cluster_label"] for r in csv.DictReader((tmp_path / "clusters.iap.csv").open())}
    assert main == {"HCV-2610-001", "HCV-2610-002"}
    assert iap_labels == {"HCV_IAP-2610-001", "HCV_IAP-2610-002"}  # numbered from 001, not after the main labels
