import pytest

from hcv_cluster import prep


@pytest.mark.parametrize(
    "region",
    [
        pytest.param("ns5b", id="single-region"),
        pytest.param("e1-e2", id="range"),
        pytest.param("core-e2", id="preset-core-e2"),
        pytest.param("core-e2-nohvr1", id="preset-core-e2-nohvr1"),
        pytest.param("e1-e2-nohvr1", id="preset-e1-e2-nohvr1"),
        pytest.param("structural", id="preset-structural"),
        pytest.param("cds", id="preset-cds"),
        pytest.param("core-e2+ns3", id="union"),
        pytest.param("core-e2-nohvr1+ns3", id="union-nohvr1"),
    ],
)
def test_validate_region_syntax_accepts_valid_expressions(region: str) -> None:
    assert prep.validate_region_syntax(region) is None


@pytest.mark.parametrize(
    ("region", "message_fragment"),
    [
        pytest.param("not-a-region", "not-a-region", id="unknown-region"),
        pytest.param("ns5b-core", "reversed", id="reversed-range"),
        pytest.param("", "", id="empty"),
    ],
)
def test_validate_region_syntax_rejects_invalid_expressions(
    region: str, message_fragment: str
) -> None:
    error = prep.validate_region_syntax(region)
    assert error is not None
    assert message_fragment in error
