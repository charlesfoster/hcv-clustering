"""``hcv-init``: ask a few questions and write the settings file (see hcv_cluster.settings).

Run it once on each computer::

    pixi run hcv-init

Every question shows its default in [brackets]: press Enter to accept it. At the registry question,
type ``browse`` to pick the folder in a window instead of typing it. Folders that do not exist
yet are created. Anything given as an option is not asked, and ``--yes`` accepts every default
without asking (for scripts).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from . import labels, settings


def _ask(question: str, default: str = "") -> str:
    answer = input(f"{question}" + (f" [{default}]" if default else "") + ": ").strip()
    return answer or default


def _yes(question: str, default: bool) -> bool:
    while True:
        answer = input(f"{question} [{'Y/n' if default else 'y/N'}]: ").strip().lower()
        if not answer:
            return default
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False


def _choose_folder(start: Path | None) -> Path | None:
    """Pick a folder in a window (Tk, or the macOS dialog); None if no window can be shown or nothing was chosen."""
    start = start if start and start.exists() else Path.home()
    try:
        import tkinter
        from tkinter import filedialog

        root = tkinter.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        chosen = filedialog.askdirectory(parent=root, initialdir=str(start), mustexist=False, title="Choose the registry folder")
        root.destroy()
        return Path(chosen) if chosen else None
    except Exception:  # no Tk, or no display (a remote terminal)
        pass
    if sys.platform == "darwin":
        script = f'POSIX path of (choose folder with prompt "Choose the registry folder" default location POSIX file "{start}")'
        done = subprocess.run(["osascript", "-e", script], capture_output=True, text=True, check=False)
        if done.returncode == 0 and done.stdout.strip():
            return Path(done.stdout.strip())
    return None


def _as_path(raw: str) -> Path:
    """A typed path: ``~`` and environment variables like $HOME are expanded; an unset variable is an error."""
    expanded = os.path.expandvars(os.path.expanduser(raw))
    if "$" in expanded:
        raise ValueError(f"{raw!r} contains a variable that is not set; type the full path instead")
    return Path(expanded)


def _registry(current: Path | None, given: Path | None, ask: bool) -> Path | None:
    """Choose the registry folder, check it, and offer to create it (and any missing folders above it)."""
    chosen = given
    while True:
        if chosen is None:
            if not ask:
                return current
            raw = _ask("Registry folder ('none' for no registry)", str(current or settings.default_registry()))
            if raw.lower() in ("none", "no"):
                return None
            if raw.lower() in ("browse", "b", "?"):
                chosen = _choose_folder(current)
                if chosen is None:
                    print("  No folder was chosen (or no window could be shown here): type the path instead.")
                    continue
                print(f"  Chosen: {chosen}")
            else:
                try:
                    chosen = _as_path(raw)
                except ValueError as exc:
                    print(f"  {exc}")
                    continue
        try:
            new = labels.check_location(chosen, new_ok=True, strict=False)  # missing folders are fine: they are created below
        except labels.RegistryError as exc:
            print(f"  Cannot use that folder: {exc}")
            if not ask:
                raise SystemExit(2)
            chosen = None
            continue
        if not new:
            print(f"  Found a registry: {labels.describe(chosen)}.")
        else:
            missing = [parent for parent in (chosen, *chosen.parents) if not parent.exists()]
            note = f" (the folder{'s' if len(missing) > 1 else ''} will be created)" if missing else ""
            if not ask or _yes(f"  There is no registry in {chosen} yet{note}. Create an empty one now?", True):
                labels.save_registry(labels.Registry(), chosen)
                print(f"  Created an empty registry in {chosen}.")
            else:
                print("  Not created: the first run will need --new-registry (or the first-run box in the GUI).")
        return chosen


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="hcv-init", description="Write the hcv-cluster settings file by answering a few questions.")
    p.add_argument("--config", type=Path, default=None, help=f"Settings file to write (default: {settings.config_path()})")
    p.add_argument("--registry", type=Path, default=None, help="Registry folder (not asked if given)")
    p.add_argument("--generate-report", action=argparse.BooleanOptionalAction, default=None,
                   help="Start the GUI's 'generate the report' box ticked (not asked if given)")
    p.add_argument("--threads", type=int, default=None, help="MAFFT threads (not asked if given)")
    p.add_argument("--yes", action="store_true", help="Accept every default without asking")
    a = p.parse_args(argv)
    ask = not a.yes
    if ask:
        print("Press Enter at any question to accept the answer shown in [brackets].\n")

    path = a.config
    while path is None:
        try:
            path = _as_path(_ask("Settings file", str(settings.config_path())) if ask else str(settings.config_path()))
        except ValueError as exc:
            print(f"  {exc}")
    try:
        current = settings.load(path)
    except ValueError as exc:
        print(f"Note: {exc}\nStarting from empty settings.")
        current = settings.Settings()

    print("\nCluster label registry: remembers each cluster's name (e.g. HCV-2611-001) between runs, so use the same folder every time.\n"
          "  Type the folder's path, or type 'browse' to pick it in a window. Folders that do not exist are created.")
    registry = _registry(current.registry, a.registry, ask)
    report = a.generate_report
    if report is None:
        report = _yes("\nStart the GUI's 'Generate the one-page report after clustering' box ticked?", current.generate_report) if ask else current.generate_report
    threads = a.threads
    if threads is None:
        default = current.threads or min(os.cpu_count() or 1, 8)
        threads = int(_ask("\nMAFFT threads (more is faster, up to the cores in this computer)", str(default))) if ask else default
    if threads < 1:
        p.error("--threads must be at least 1")

    written = settings.save(registry=registry, generate_report=report, threads=threads, path=path)
    print(f"\nSaved {written}:\n" + "".join(f"  {line}\n" for line in written.read_text().splitlines()))
    if written != settings.config_path():
        print(f"hcv-cluster only reads this file if you set HCV_CLUSTER_CONFIG={written} in your environment.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
