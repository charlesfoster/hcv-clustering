from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_gui_renders_composite_and_small_multiple_networks(tmp_path: Path) -> None:
    results = tmp_path / "results"
    results.mkdir()
    (results / "clusters.csv").write_text(
        "sample_id,genotype,cluster_id,source_cluster_id,cluster_size\n"
        "A,1a,1a_C0001,C0001,2\n"
        "B,1a,1a_C0001,C0001,2\n"
        "C,3a,3a_C0001,C0001,1\n",
        encoding="utf-8",
    )
    (results / "links.csv").write_text(
        "source,target,distance\nA,B,0.01\n",
        encoding="utf-8",
    )
    (results / "metadata.csv").write_text(
        "sample_id,location,age_range,injecting_status,indigenous_status\n"
        "A,Prison A,0–30,yes,no\n"
        "B,Prison B,31–60,no,yes\n"
        "C,Prison A,61+,yes,no\n",
        encoding="utf-8",
    )

    app_path = Path(__file__).resolve().parents[1] / "hcv_cluster_gui.py"
    app = AppTest.from_file(str(app_path), default_timeout=15)
    app.session_state["last_run"] = {"outdir": str(results), "distance": "tn93"}
    app.session_state["cluster_view"] = {
        "outdir": str(results),
        "scope": "All genotypes",
        "metric": "tn93",
        "layout": "Packed clusters",
    }

    app.run()
    assert not app.exception
    assert len(app.get("plotly_chart")) == 1

    app.session_state["network_display_mode"] = "Small multiples"
    app.session_state["network_small_multiple_fields"] = [
        "location",
        "age_range",
        "injecting_status",
        "indigenous_status",
    ]
    app.run()
    assert not app.exception
    assert len(app.get("plotly_chart")) == 1
