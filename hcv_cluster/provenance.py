#!/usr/bin/env python3
"""Run provenance: what produced a set of results.

``run_info.json`` records the pipeline version and commit, the versions of the external
tools, the reference data in use, and the parameters, so a result can be traced back to
exactly the code and tools that made it.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path
from typing import Any, Mapping

from . import __version__, genotyping, prep

TOOL_VERSION_COMMANDS = {
    "tn93": ["tn93", "--version"],
    "minimap2": ["minimap2", "--version"],
    "mafft": ["mafft", "--version"],  # MAFFT prints its version on stderr
}


def tool_version(command: list[str]) -> str | None:
    """First line of a tool's version output, or None if it cannot be run."""
    try:
        done = subprocess.run(command, capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    text = (done.stdout or done.stderr).strip()
    return text.splitlines()[0] if text else None


def git_commit() -> str | None:
    """Commit of the running code: HCV_CLUSTER_COMMIT (set in the container), else git."""
    if os.environ.get("HCV_CLUSTER_COMMIT") not in (None, "", "unknown"):
        return os.environ["HCV_CLUSTER_COMMIT"]
    try:
        done = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return done.stdout.strip() or None if done.returncode == 0 else None


def git_dirty() -> bool | None:
    """Whether the checkout has uncommitted changes (None if unknown, e.g. in a container)."""
    try:
        done = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "status", "--porcelain"],
            capture_output=True, text=True, timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return bool(done.stdout.strip()) if done.returncode == 0 else None


def _sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16] if path.exists() else None


def collect(params: Mapping[str, Any] | None = None, extra: Mapping[str, Any] | None = None) -> dict[str, Any]:
    info: dict[str, Any] = {
        "pipeline": "hcv-cluster",
        "pipeline_version": __version__,
        "git_commit": git_commit(),
        "git_uncommitted_changes": None if os.environ.get("HCV_CLUSTER_COMMIT") else git_dirty(),
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "tools": {name: tool_version(command) for name, command in TOOL_VERSION_COMMANDS.items()},
        "reference_data": {
            "genotype_panel_sha256": _sha256(genotyping.DEFAULT_PANEL),
            "subtypes_with_reference": len(prep.DEFAULT_REFERENCE_CATALOG),
        },
        "parameters": dict(params or {}),
    }
    if extra:
        info.update(extra)
    return info


def write(path: Path, params: Mapping[str, Any] | None = None, extra: Mapping[str, Any] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(collect(params, extra), indent=2, sort_keys=True, default=str) + "\n")
    return path
