"""Streamlit entry point: `pixi run hcv-cluster-gui` (or `streamlit run streamlit_app.py`).

Streamlit re-executes this script on every interaction, so the page is built by a
function call rather than at import time.
"""
from hcv_cluster import gui

gui.main()
