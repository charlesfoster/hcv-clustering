"""One release number, kept in step everywhere it appears (see docs/RELEASING.md)."""
import re
from pathlib import Path

import hcv_cluster

ROOT = Path(__file__).resolve().parents[1]


def _find(path: str, pattern: str) -> list[str]:
    return re.findall(pattern, (ROOT / path).read_text(), flags=re.M)


def test_version_is_the_same_everywhere() -> None:
    version = hcv_cluster.__version__
    assert re.fullmatch(r"\d+\.\d+\.\d+", version)
    assert _find("nextflow.config", r'^\s*version\s*=\s*"([^"]+)"') == [version]
    assert _find("conf/hcv_cluster.config", r'^\s*pipeline_version\s*=\s*"([^"]+)"') == [version]
    assert _find("pixi.toml", r'^version = "([^"]+)"') == [version]
    assert "dynamic = [\"version\"]" in (ROOT / "pyproject.toml").read_text()  # read from __init__
    assert _find("CHANGELOG.md", r"^## \[(\d+\.\d+\.\d+)\]")[0] == version  # newest entry first
