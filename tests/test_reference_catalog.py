"""Validate DEFAULT_REFERENCE_CATALOG against ICTV's confirmed subtype table.

Most GenBank records for the older HCV reference genomes (D00944 = HC-J6,
D17763 = NZL1, D84262 = Th580, ...) carry no subtype string of their own, so the
catalog's subtype assignments cannot be checked against the records themselves.
ICTV Flaviviridae/Hepacivirus Table 1 is the authority that makes them auditable,
and reference_data/ictv_hcv_subtypes.tsv mirrors it.
"""

from pathlib import Path

import hcv_cluster_prep

ICTV_TABLE = Path(__file__).resolve().parent.parent / "reference_data" / "ictv_hcv_subtypes.tsv"

# Catalog entries that deliberately use a genome other than ICTV's first-listed
# exemplar. Each is a different isolate of the *same* ICTV-confirmed subtype, kept
# because it carries the complete-CDS annotation the prep stage needs to resolve
# core-E2 boundaries. Every one of these self-declares its subtype in its own
# GenBank record, so the assignment is independently verifiable.
#
# 1a and 1b predate the catalog expansion and anchor all existing alignments and
# cached coordinates; they are not swapped for ICTV's exemplars.
ICTV_EXEMPLAR_DEVIATIONS = {
    "1a": "EF407457.1",  # isolate 1003; /note="subtype: 1a"
    "1b": "EU781828.1",  # isolate TN28; DEFINITION "subtype 1b"
    "2j": "HM777359.1",  # isolate C1292; DEFINITION "genotype 2j"
    "4a": "DQ418789.1",  # isolate L835; DEFINITION "subtype 4a"
    "4f": "EU392174.1",  # isolate PS4; DEFINITION "genotype 4f"
    "4k": "EU392171.1",  # isolate PB65185; DEFINITION "genotype 4k"
    "6i": "DQ835762.1",  # isolate C-0159; DEFINITION "subtype 6i"
    "6j": "DQ835761.1",  # isolate C-0667; DEFINITION "subtype 6j"
    "6m": "DQ835766.1",  # isolate C-0192; DEFINITION "subtype 6m"
}

# ICTV subtypes with no configured reference. 6xa-6xj are recent designations the
# bundled genotyping panel does not call; a reference can be supplied at runtime
# via --reference-map.
KNOWN_UNCOVERED = {"6xa", "6xb", "6xc", "6xd", "6xe", "6xf", "6xg", "6xh", "6xi", "6xj"}


def load_ictv_table() -> dict[str, list[str]]:
    table: dict[str, list[str]] = {}
    for line in ICTV_TABLE.read_text().splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("subtype\t"):
            continue
        subtype, accessions = line.split("\t")
        table[subtype.strip()] = [a.strip() for a in accessions.split(",")]
    return table


def test_ictv_table_has_all_94_confirmed_subtypes() -> None:
    assert len(load_ictv_table()) == 94


def test_every_catalog_subtype_is_ictv_confirmed() -> None:
    ictv = load_ictv_table()
    unknown = sorted(set(hcv_cluster_prep.DEFAULT_REFERENCE_CATALOG) - set(ictv))
    assert unknown == [], f"catalog names subtypes ICTV does not confirm: {unknown}"


def test_catalog_accessions_match_ictv_or_a_pinned_deviation() -> None:
    """Every accession is either ICTV's exemplar or a documented substitute."""
    ictv = load_ictv_table()
    unexpected = []
    for subtype, spec in sorted(hcv_cluster_prep.DEFAULT_REFERENCE_CATALOG.items()):
        versionless = spec.accession.split(".")[0]
        if versionless in ictv[subtype]:
            continue
        if ICTV_EXEMPLAR_DEVIATIONS.get(subtype) == spec.accession:
            continue
        unexpected.append(f"{subtype}={spec.accession} (ICTV: {','.join(ictv[subtype])})")
    assert unexpected == [], (
        "catalog accessions diverge from ICTV without a pinned deviation: "
        + "; ".join(unexpected)
    )


def test_pinned_deviations_are_all_still_in_the_catalog() -> None:
    """A deviation left behind after an accession change would silently stop checking."""
    catalog = hcv_cluster_prep.DEFAULT_REFERENCE_CATALOG
    stale = [
        subtype
        for subtype, accession in ICTV_EXEMPLAR_DEVIATIONS.items()
        if subtype not in catalog or catalog[subtype].accession != accession
    ]
    assert stale == [], f"stale pinned deviations: {stale}"


def test_uncovered_subtypes_are_exactly_the_documented_ones() -> None:
    ictv = load_ictv_table()
    uncovered = set(ictv) - set(hcv_cluster_prep.DEFAULT_REFERENCE_CATALOG)
    assert uncovered == KNOWN_UNCOVERED


def test_reference_spec_genotype_matches_its_catalog_key() -> None:
    mismatched = [
        (key, spec.genotype)
        for key, spec in hcv_cluster_prep.DEFAULT_REFERENCE_CATALOG.items()
        if key != spec.genotype
    ]
    assert mismatched == []


def test_accessions_are_unique_across_subtypes() -> None:
    seen: dict[str, str] = {}
    duplicates = []
    for subtype, spec in sorted(hcv_cluster_prep.DEFAULT_REFERENCE_CATALOG.items()):
        if spec.accession in seen:
            duplicates.append(f"{spec.accession} used by {seen[spec.accession]} and {subtype}")
        seen[spec.accession] = subtype
    assert duplicates == []
