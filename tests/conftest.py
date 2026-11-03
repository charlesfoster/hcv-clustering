import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


import pytest  # noqa: E402


@pytest.fixture(autouse=True)
def _no_user_settings(tmp_path_factory, monkeypatch):
    """Tests never read or write the developer's real ~/.hcv-cluster."""
    from hcv_cluster import settings

    home = tmp_path_factory.mktemp("home")
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv(settings.CONFIG_ENV, str(home / "config.toml"))


@pytest.fixture(autouse=True)
def _never_launch_a_browser(monkeypatch):
    """No test opens a real Chrome: with a temporary HOME it cannot find the login keychain and macOS shows a dialog."""
    from hcv_cluster import report

    def refuse(command):
        raise OSError("tests do not launch a browser")

    monkeypatch.setattr(report, "_launch", refuse)  # write_pdf then returns None, as when there is no Chrome
