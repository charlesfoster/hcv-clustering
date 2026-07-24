import hcv_cluster_prep


def test_validate_region_syntax_accepts_single_region() -> None:
    assert hcv_cluster_prep.validate_region_syntax("ns5b") is None


def test_validate_region_syntax_accepts_range() -> None:
    assert hcv_cluster_prep.validate_region_syntax("e1-e2") is None


def test_validate_region_syntax_accepts_preset() -> None:
    assert hcv_cluster_prep.validate_region_syntax("core-e2") is None
    assert hcv_cluster_prep.validate_region_syntax("core-e2-nohvr1") is None
    assert hcv_cluster_prep.validate_region_syntax("structural") is None
    assert hcv_cluster_prep.validate_region_syntax("cds") is None


def test_validate_region_syntax_accepts_union() -> None:
    assert hcv_cluster_prep.validate_region_syntax("core-e2+ns3") is None
    assert hcv_cluster_prep.validate_region_syntax("core-e2-nohvr1+ns3") is None


def test_validate_region_syntax_rejects_unknown_region() -> None:
    error = hcv_cluster_prep.validate_region_syntax("not-a-region")
    assert error is not None
    assert "not-a-region" in error


def test_validate_region_syntax_rejects_reversed_range() -> None:
    error = hcv_cluster_prep.validate_region_syntax("ns5b-core")
    assert error is not None
    assert "reversed" in error


def test_validate_region_syntax_rejects_empty() -> None:
    error = hcv_cluster_prep.validate_region_syntax("")
    assert error is not None
