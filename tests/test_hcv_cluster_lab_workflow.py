"""The lab workflow: settings, a registry that refuses a wrong location or wrong data, backups, and the report on demand."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from hcv_cluster import cli, config_init, labels, report, settings

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "assets" / "test" / "example.fasta"
DEFINITION = {"metric": "tn93", "region": "e1-e2", "threshold": 0.045}


def _members(prefix: str, n: int = 12) -> dict[tuple[str, str], str]:
    return {("1a", f"{prefix}{i}"): "1a_C0001" for i in range(n)}


def _registry(directory: Path) -> Path:
    labels.save_registry(labels.Registry(), directory)
    return directory


def test_settings_are_read_checked_and_saved(tmp_path: Path) -> None:
    path = Path(os.environ[settings.CONFIG_ENV])
    assert settings.load() == settings.Settings()  # no file: no defaults
    path.write_text('registry = "~/labels"\ngenerate_report = true\nthreads = 12\n')
    assert settings.load() == settings.Settings(Path("~/labels").expanduser(), True, 12)
    settings.save(registry=tmp_path / "reg")  # the other settings stay
    assert settings.load() == settings.Settings(tmp_path / "reg", True, 12)
    for bad in ("threads = 0", 'generate_report = "yes"', "registry = 5", "threads = ["):
        path.write_text(bad)
        with pytest.raises(ValueError, match="config.toml"):
            settings.load()


def test_a_registry_location_that_looks_wrong_is_refused(tmp_path: Path) -> None:
    good = _registry(tmp_path / "reg")
    assert labels.check_location(good, strict=True, new_ok=False) is False
    absent = tmp_path / "typo"
    assert labels.check_location(absent) is True  # lenient: a new registry may be created
    with pytest.raises(labels.RegistryError, match="no registry at"):
        labels.check_location(absent, strict=True, new_ok=False)  # strict: a typo, or a drive that is not mounted
    assert labels.check_location(absent, strict=True, new_ok=True) is True
    with pytest.raises(labels.RegistryError, match="does not exist"):
        labels.check_location(tmp_path / "no_parent" / "reg", strict=True, new_ok=True)
    stuff = tmp_path / "documents"
    stuff.mkdir()
    (stuff / "notes.txt").write_text("x")
    with pytest.raises(labels.RegistryError, match="does not look like"):
        labels.check_location(stuff)  # never write registry files among other files
    (good / labels.MEMBERSHIP_FILE).write_text("")
    with pytest.raises(labels.RegistryError, match="damaged"):
        labels.check_location(good)
    (tmp_path / "afile").write_text("x")
    with pytest.raises(labels.RegistryError, match="is a file"):
        labels.check_location(tmp_path / "afile")


def test_saving_labels_backs_up_the_registry_and_refuses_the_wrong_data(tmp_path: Path) -> None:
    reg = tmp_path / "reg"
    kwargs = dict(run_date="2026-11-03", outdir="out", prefix="HCV")
    labels.label_run(reg, "tn93", DEFINITION, _members("a"), run_id="r1", commit=True, **kwargs)
    assert not (reg / labels.BACKUP_DIR).exists()  # nothing to back up before the first save
    labels.label_run(reg, "tn93", DEFINITION, {**_members("a"), ("1a", "new"): "1a_C0001"}, run_id="r2", commit=True, **kwargs)
    assert (reg / labels.BACKUP_DIR / "r2" / labels.MEMBERSHIP_FILE).is_file()  # the state before run 2
    with pytest.raises(labels.RegistryError, match="only 0 of the 13 samples"):  # different data entirely
        labels.label_run(reg, "tn93", DEFINITION, _members("z"), run_id="r3", commit=True, **kwargs)
    assert "r3" not in {row["run_id"] for row in labels.load_registry(reg).runs}  # nothing was saved
    labels.label_run(reg, "tn93", DEFINITION, _members("z"), run_id="r3", commit=False, **kwargs)  # an unsaved run only warns
    labels.label_run(reg, "tn93", DEFINITION, _members("z"), run_id="r3", commit=True, allow_mismatch=True, **kwargs)


def test_default_outdir_is_dated_and_the_registry_offers_the_last_saved_run(tmp_path: Path) -> None:
    first = cli.default_outdir(tmp_path)
    assert first.name.startswith("results_20")
    first.mkdir()
    assert cli.default_outdir(tmp_path).name == f"{first.name}_2"  # never lands on an earlier run

    reg = tmp_path / "reg"
    done, test_run, gone = (tmp_path / name for name in ("done", "test", "gone"))
    for folder in (done, test_run):
        (folder / "by_genotype" / "1a").mkdir(parents=True)
        (folder / "by_genotype" / "1a" / "prep.aligned.fasta").write_text(">r\nA\n")
    kwargs = dict(run_date="2026-11-03", prefix="HCV")
    for run_id, folder, commit in (("r1", done, True), ("r2", test_run, False), ("r3", gone, True)):
        labels.label_run(reg, "tn93", DEFINITION, _members("a"), run_id=run_id, outdir=str(folder), commit=commit, **kwargs)
    assert labels.last_run_outdir(reg) == (done, "2026-11-03")  # skips the unsaved test run and the deleted folder


def test_run_uses_the_registry_in_the_settings_and_writes_the_report(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    if shutil.which("tn93") is None:
        pytest.skip("needs the clustering tools (pixi run -e nf test)")
    registry = tmp_path / "labels"
    settings.save(registry=registry)  # a registry that does not exist yet: the run must not create it by accident
    out = tmp_path / "out"
    argv = ["run", "-i", str(FIXTURE), "-o", str(out)]
    assert cli.main(argv) == 1 and "no registry at" in capsys.readouterr().err and not out.exists()
    assert cli.main([*argv, "--new-registry", "--generate-report"]) == 0  # first run: created, saved, reported
    assert (registry / labels.TRACKS_FILE).is_file() and (out / "report.html").is_file()
    assert (registry / labels.SNAPSHOT_FILE).is_file()  # the next report will show what changed since this one
    assert cli.main([*argv, "--no-registry"]) == 0 and "no cluster label registry" in capsys.readouterr().out


def test_gui_blocks_a_run_until_the_registry_location_is_right(tmp_path: Path) -> None:
    from streamlit.testing.v1 import AppTest

    settings.save(registry=tmp_path / "typo")
    app = AppTest.from_file(str(ROOT / "streamlit_app.py"), default_timeout=30).run()
    assert any("no registry at" in error.value for error in app.error)
    assert [b for b in app.button if b.label == "Run pipeline"][0].disabled


def test_init_writes_the_settings_from_options_without_asking(tmp_path: Path) -> None:
    reg = tmp_path / "labels"
    assert config_init.main(["--yes", "--registry", str(reg), "--generate-report", "--threads", "4"]) == 0
    assert settings.load() == settings.Settings(reg, True, 4)
    assert labels.check_location(reg, strict=True, new_ok=False) is False  # an empty registry was created, so the first run needs no flag
    with pytest.raises(SystemExit):  # a folder that is not a registry is refused, not written into
        (tmp_path / "docs").mkdir()
        (tmp_path / "docs" / "notes.txt").write_text("x")
        config_init.main(["--yes", "--registry", str(tmp_path / "docs")])


def test_init_asks_its_questions_and_creates_missing_folders(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reg = tmp_path / "new" / "folders" / "labels"  # none of the folders above it exist yet
    answers = iter([
        "",            # settings file: Enter keeps the default
        "$NOT_SET_X/labels",  # a variable that is not set: refused with a message, then it asks again
        str(reg), "",  # the registry folder, then Enter accepts the default (yes) to create it
        "y",           # start the report box ticked
        "",            # threads: Enter keeps the default
    ])
    monkeypatch.delenv("NOT_SET_X", raising=False)
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    assert config_init.main([]) == 0
    saved = settings.load()
    assert saved.registry == reg and saved.generate_report is True and saved.threads >= 1
    assert labels.check_location(reg, strict=True, new_ok=False) is False  # created, with the folders above it


def test_init_can_pick_the_registry_folder_in_a_window(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    picked = tmp_path / "picked"
    answers = iter(["", "browse", "", "", ""])  # Enter, browse (a window), Enter = create, Enter, Enter
    monkeypatch.setattr("builtins.input", lambda prompt="": next(answers))
    monkeypatch.setattr(config_init, "_choose_folder", lambda start: picked)
    assert config_init.main([]) == 0 and settings.load().registry == picked
    # no window available (a remote terminal): it says so and asks for a typed path instead
    answers = iter(["", "browse", str(tmp_path / "typed"), "", "", ""])
    monkeypatch.setattr(config_init, "_choose_folder", lambda start: None)
    assert config_init.main([]) == 0 and settings.load().registry == tmp_path / "typed"


def test_init_defaults_live_in_one_hidden_folder_in_the_home_directory(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(settings.CONFIG_ENV)  # as on a real computer: no override
    home = Path(os.environ["HOME"])
    assert settings.config_path() == home / ".hcv-cluster" / "config.toml"
    monkeypatch.setattr("builtins.input", lambda prompt="": "")  # Enter at every question
    assert config_init.main([]) == 0
    saved = settings.load()
    assert saved.registry == home / ".hcv-cluster" / "label_registry" and labels.check_location(saved.registry, new_ok=False, strict=True) is False
    assert (home / ".hcv-cluster" / "config.toml").is_file()
