"""Per-user defaults for hcv-cluster, read from a small TOML file.

A lab sets these once and every run, from the command line or the GUI, picks them up::

    registry = "/shared/hcv/cluster_registry"   # where persistent cluster labels live
    generate_report = true                      # the GUI's "generate the report" box starts ticked
    threads = 12                                # MAFFT threads

``hcv-init`` writes this file for you by asking the questions.

The file is ``~/.hcv-cluster/config.toml``, or the path in ``HCV_CLUSTER_CONFIG``.
A missing file means no defaults. Nothing here is required: the command line and the GUI
always say what they use.
"""

from __future__ import annotations

import json
import os
import tomllib
from dataclasses import dataclass
from pathlib import Path

CONFIG_ENV = "HCV_CLUSTER_CONFIG"


@dataclass(frozen=True)
class Settings:
    registry: Path | None = None
    generate_report: bool = False
    threads: int | None = None


def home_dir() -> Path:
    """The folder hcv-cluster keeps its settings, and by default its registry, in."""
    return Path.home() / ".hcv-cluster"


def config_path() -> Path:
    return Path(os.environ.get(CONFIG_ENV) or home_dir() / "config.toml")


def default_registry() -> Path:
    """The registry folder ``hcv-init`` offers when none is set."""
    return home_dir() / "label_registry"


def _read(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        return tomllib.loads(path.read_text(encoding="utf-8"))
    except (tomllib.TOMLDecodeError, OSError) as exc:
        raise ValueError(f"cannot read the settings file {path}: {exc}") from exc


def load(path: Path | None = None) -> Settings:
    path = path or config_path()
    raw = _read(path)
    registry, report, threads = raw.get("registry"), raw.get("generate_report", False), raw.get("threads")
    if registry is not None and not (isinstance(registry, str) and registry.strip()):
        raise ValueError(f"{path}: registry must be a path in quotes")
    if not isinstance(report, bool):
        raise ValueError(f"{path}: generate_report must be true or false")
    if threads is not None and not (isinstance(threads, int) and not isinstance(threads, bool) and threads >= 1):
        raise ValueError(f"{path}: threads must be a whole number of at least 1")
    return Settings(Path(registry).expanduser() if registry else None, report, threads)


def save(*, registry: Path | None = None, generate_report: bool | None = None, threads: int | None = None,
         path: Path | None = None) -> Path:
    """Remember the given settings; any others already in the file are kept."""
    path = path or config_path()
    values = _read(path)
    for key, value in (("registry", None if registry is None else str(registry)),
                       ("generate_report", generate_report), ("threads", threads)):
        if value is not None:
            values[key] = value
    lines = []
    for key, value in values.items():
        lines.append(f"{key} = {'true' if value is True else 'false' if value is False else json.dumps(value)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
