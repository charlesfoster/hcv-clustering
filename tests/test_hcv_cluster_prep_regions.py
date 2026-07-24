import pytest

import hcv_cluster_prep


def _synthetic_regions() -> dict[str, hcv_cluster_prep.RegionSegment]:
    # Loosely H77-shaped coordinates: core, e1, e2 contiguous, generous e2 length.
    return {
        "core": hcv_cluster_prep.RegionSegment("core", 342, 914, "test", ""),
        "e1": hcv_cluster_prep.RegionSegment("e1", 915, 1490, "test", ""),
        "e2": hcv_cluster_prep.RegionSegment("e2", 1491, 2768, "test", ""),
    }


def test_mask_hvr1_from_e2_trims_n_terminal_span() -> None:
    e2 = hcv_cluster_prep.RegionSegment("e2", 1491, 2768, "test", "detail")
    trimmed = hcv_cluster_prep.mask_hvr1_from_e2(e2)
    assert trimmed.start == 1491 + hcv_cluster_prep.HVR1_LENGTH_NT
    assert trimmed.end == e2.end
    assert trimmed.length == e2.length - hcv_cluster_prep.HVR1_LENGTH_NT


def test_mask_hvr1_from_e2_rejects_e2_shorter_than_mask() -> None:
    tiny_e2 = hcv_cluster_prep.RegionSegment("e2", 1491, 1491 + 10, "test", "")
    with pytest.raises(hcv_cluster_prep.HcvPrepError):
        hcv_cluster_prep.mask_hvr1_from_e2(tiny_e2)


def test_resolve_core_e2_nohvr1_returns_two_discontiguous_segments() -> None:
    regions = _synthetic_regions()
    segments = hcv_cluster_prep.resolve_core_e2_nohvr1(regions)
    assert len(segments) == 2
    core_e1_segment, e2_segment = segments
    assert core_e1_segment.start == regions["core"].start
    assert core_e1_segment.end == regions["e1"].end
    assert e2_segment.start == regions["e2"].start + hcv_cluster_prep.HVR1_LENGTH_NT
    assert e2_segment.end == regions["e2"].end


def test_core_e2_nohvr1_selection_excludes_hvr1_length_from_total() -> None:
    regions = _synthetic_regions()
    selection = hcv_cluster_prep.resolve_region_selection("core-e2-nohvr1", regions, strategy="fixed")
    full_span = regions["e2"].end - regions["core"].start + 1
    assert selection.length == full_span - hcv_cluster_prep.HVR1_LENGTH_NT
    # Overall start/end still span the full genomic range, including the masked gap.
    assert selection.start == regions["core"].start
    assert selection.end == regions["e2"].end


def test_core_e2_full_selection_is_unaffected_and_contiguous() -> None:
    regions = _synthetic_regions()
    selection = hcv_cluster_prep.resolve_region_selection("core-e2", regions, strategy="fixed")
    assert len(selection.segments) == 1
    assert selection.length == regions["e2"].end - regions["core"].start + 1


def test_structural_preset_matches_core_e2_nohvr1() -> None:
    regions = _synthetic_regions()
    masked = hcv_cluster_prep.resolve_region_selection("core-e2-nohvr1", regions, strategy="fixed")
    via_preset = hcv_cluster_prep.resolve_region_selection("structural", regions, strategy="fixed")
    assert via_preset.length == masked.length
    assert via_preset.segments == masked.segments


def test_extract_region_selection_concatenates_around_masked_hvr1() -> None:
    # 10 reference bases for core (1-10), 5 for e1 (11-15), 10 for e2 (16-25);
    # mask the first 4 bases of e2 (positions 16-19) to simulate HVR1 removal.
    regions = {
        "core": hcv_cluster_prep.RegionSegment("core", 1, 10, "test", ""),
        "e1": hcv_cluster_prep.RegionSegment("e1", 11, 15, "test", ""),
        "e2": hcv_cluster_prep.RegionSegment("e2", 16, 25, "test", ""),
    }
    reference_aligned = "A" * 25
    query_aligned = "C" * 25
    selection = hcv_cluster_prep.RegionSelection(
        expression="core-e2-nohvr1",
        strategy="fixed",
        segments=(
            hcv_cluster_prep.RegionSegment("core-e1", 1, 15, "test", ""),
            hcv_cluster_prep.RegionSegment("e2-nohvr1", 20, 25, "test", ""),
        ),
    )
    extracted = hcv_cluster_prep.extract_region_selection(
        reference_aligned, {"sample": query_aligned}, selection
    )
    assert len(extracted["sample"]) == 15 + 6
