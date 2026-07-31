#!/usr/bin/env python3
"""Streamlit GUI for HCV genomic clustering. Thin wrapper around hcv_cluster.command_run."""

from __future__ import annotations

import contextlib
import csv
import io
import os
import signal
import threading
import time
from pathlib import Path

import streamlit as st

import assign_hcv_genotypes_from_fasta
import hcv_cluster
import hcv_cluster_prep
import hcv_cluster_viz

st.set_page_config(page_title="HCV Clustering", layout="wide")

REGION_PRESETS_UI = (
    ("core-e2-nohvr1", "Core–E2, HVR1 masked (default)"),
    ("core-e2", "Core–E2, full (HVR1 included)"),
    ("e1-e2", "Envelope (E1–E2)"),
    ("ns2-ns5b", "Nonstructural (NS2–NS5B)"),
    ("cds", "Whole coding sequence"),
    ("ns5b", "NS5B only"),
)


def _region_picker() -> tuple[str, str | None]:
    if "region_expr" not in st.session_state:
        st.session_state.region_expr = "core-e2-nohvr1"

    st.write("Region")
    preset_cols = st.columns(len(REGION_PRESETS_UI))
    for col, (expr, label) in zip(preset_cols, REGION_PRESETS_UI):
        if col.button(label, key=f"region_preset_{expr}", width="stretch"):
            st.session_state.region_expr = expr

    region = st.text_input(
        "Region expression",
        key="region_expr",
        help=(
            "Click a preset above, or type a custom expression: individual regions "
            "(core e1 e2 p7 ns2 ns3 ns4a ns4b ns5a ns5b), ranges (e1-e2, ns3-ns5b), "
            "or unions (core-e2-nohvr1+ns3). core-e2-nohvr1 (default) masks HVR1; "
            "core-e2 keeps it. Checked against known region names, not against any "
            "specific reference's annotations."
        ),
    )
    region_error = hcv_cluster_prep.validate_region_syntax(region) if region.strip() else "Region expression cannot be empty"
    if region_error:
        st.error(region_error)
    return region, region_error


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _metric_file_paths(genotype_dir: Path, distance: str, metric: str) -> tuple[Path, Path]:
    if distance == "both" and metric == "snp":
        return genotype_dir / "snp_links.csv", genotype_dir / "snp_clusters.csv"
    return genotype_dir / "links.csv", genotype_dir / "clusters.csv"


def _show_cluster_table(title: str, path: Path) -> None:
    if not path.exists():
        return
    rows = _read_csv_rows(path)
    st.subheader(title)
    st.dataframe(rows, width="stretch")
    st.download_button(
        f"Download {path.name}",
        data=path.read_bytes(),
        file_name=path.name,
        mime="text/csv",
        key=f"download_{path}",
    )


def _run_tab() -> None:
    st.header("Run clustering")

    uploaded = st.file_uploader("Input multi-FASTA of HCV consensus sequences", type=["fasta", "fa", "fna"])
    outdir = Path(st.text_input("Output directory", value="results_gui"))

    region, region_error = _region_picker()

    col1, col2, col3 = st.columns(3)
    with col1:
        distance = st.selectbox("Distance metric", ("tn93", "snp", "both"), index=0)
    with col2:
        default_threshold = hcv_cluster.resolve_threshold(region, None)
        override_threshold = st.checkbox("Override default threshold")
        if override_threshold:
            threshold = st.number_input(
                "Threshold", min_value=0.0, max_value=1.0, value=default_threshold, step=0.001, format="%.4f"
            )
        else:
            threshold = None
            st.caption(f"Default for '{region}': **{default_threshold:.4g}**")
    with col3:
        threads = st.number_input("MAFFT threads", min_value=1, value=1, step=1)

    verbose = st.checkbox("Verbose logging")
    dry_run = st.checkbox("Dry run (preview commands only, no execution)")

    with st.expander("Advanced options"):
        genotype_filter_raw = st.text_input(
            "Restrict to genotype(s) (comma-separated, blank = auto-detect all)",
            value="",
            help=(
                "Only process these detected genotypes; sequences assigned to any other "
                "genotype are skipped. Leave blank to process every genotype the "
                "auto-detection step finds. Comma-separated, e.g. '1a,3a'."
            ),
        )
        adv_col1, adv_col2 = st.columns(2)
        with adv_col1:
            min_coverage = st.number_input(
                "Min region coverage fraction",
                value=0.8,
                step=0.05,
                format="%.2f",
                help=(
                    "Minimum fraction of the selected region that must be covered by real "
                    "bases (not alignment gaps or Ns) for a sequence to pass QC and be "
                    "included in clustering. Checked after region extraction, so masked "
                    "spans like HVR1 don't count against it."
                ),
            )
            min_identity = st.number_input(
                "Min minimap2 identity",
                value=0.75,
                step=0.05,
                format="%.2f",
                help=(
                    "Minimum alignment identity to the best-matching reference during the "
                    "genotyping step. Sequences below this fail genotype assignment "
                    "(qc_fail_reason=low_identity) and never reach clustering."
                ),
            )
            min_query_coverage = st.number_input(
                "Min minimap2 query coverage",
                value=0.50,
                step=0.05,
                format="%.2f",
                help=(
                    "Minimum fraction of the input sequence that must align to the "
                    "best-matching reference during the genotyping step. Sequences below "
                    "this fail genotype assignment (qc_fail_reason=low_query_coverage)."
                ),
            )
            region_strategy = st.selectbox(
                "Region strategy",
                ("fixed", "max-usable"),
                index=0,
                help=(
                    "fixed (default): use the Region expression above exactly as given. "
                    "max-usable: instead of a fixed region, search for the widest "
                    "reproducible high-coverage window — useful when sequencing coverage "
                    "varies a lot between samples and a fixed region would fail too many "
                    "of them."
                ),
            )
            genotype_validation = st.selectbox(
                "Genotype validation",
                ("auto", "headers", "kmer", "none"),
                index=0,
                help=(
                    "A second sanity check that each sequence's assigned genotype is "
                    "correct, run before alignment. auto (default): use genotype tags in "
                    "FASTA headers if present (e.g. 'sample_id genotype=1a'), otherwise "
                    "fall back to comparing k-mers against reference sequences. headers: "
                    "require header tags and error if missing. kmer: validate by k-mer "
                    "comparison. none: skip this check entirely."
                ),
            )
            ambiguities = st.selectbox(
                "TN93 ambiguity handling",
                ("average", "resolve", "skip", "gapmm"),
                index=0,
                help=(
                    "How tn93 handles real IUPAC ambiguity codes (R, Y, W, S, K, M) "
                    "when computing distance. N is always masked to a gap before "
                    "clustering regardless of this setting, so low-depth-masked "
                    "positions can't distort distances (see "
                    "docs/threshold_rationale.md). average (default): proportional "
                    "treatment over the possible resolutions — e.g. R-A counts as "
                    "0.5 A-A + 0.5 G-A. resolve: pick the interpretation that "
                    "minimizes distance. skip: exclude ambiguous positions from the "
                    "comparison. gapmm: an alternative gap-aware handling mode — see "
                    "the tn93 tool's own documentation."
                ),
            )
        with adv_col2:
            refs_dir = Path(
                st.text_input(
                    "References cache dir",
                    value="refs",
                    help="Directory where downloaded GenBank/FASTA reference records are cached and reused across runs, so they aren't re-downloaded every time.",
                )
            )
            panel_fasta = Path(
                st.text_input(
                    "Genotyping panel FASTA",
                    value=str(assign_hcv_genotypes_from_fasta.DEFAULT_PANEL),
                    help="Reference panel used for the initial competitive-alignment genotyping step (minimap2 against every sequence in this file).",
                )
            )
            minimap2_exe = st.text_input(
                "minimap2 executable",
                value="minimap2",
                help="Path or command name for minimap2, if it isn't on your PATH under the default name.",
            )
            mafft_exe = st.text_input(
                "mafft executable",
                value="mafft",
                help="Path or command name for MAFFT, if it isn't on your PATH under the default name.",
            )
            tn93_exe = st.text_input(
                "tn93 executable",
                value="tn93",
                help="Path or command name for the tn93 tool, if it isn't on your PATH under the default name.",
            )
            ncbi_email = st.text_input(
                "NCBI email (for reference downloads)",
                value="",
                help="Email address sent to NCBI's E-utilities when downloading reference GenBank records — required by NCBI's usage policy, not used for anything else.",
            )
        force_download = st.checkbox(
            "Force re-download of cached reference records",
            help="Re-download reference GenBank/FASTA records even if already cached locally in the References cache dir.",
        )
        keep_temp = st.checkbox(
            "Keep temporary MAFFT files",
            help="Keep MAFFT's intermediate alignment files after the run instead of deleting them — useful for debugging alignment issues.",
        )
        keep_paf = st.checkbox(
            "Keep raw minimap2 PAF file",
            help="Keep the raw minimap2 PAF alignment file from the genotyping step instead of deleting it — useful for debugging genotype assignment.",
        )

    run_clicked = st.button(
        "Run pipeline", type="primary", disabled=uploaded is None or region_error is not None
    )

    if run_clicked:
        if uploaded is None:
            st.error("Upload an input FASTA first.")
        elif region_error is not None:
            st.error(region_error)
        else:
            outdir.mkdir(parents=True, exist_ok=True)
            input_path = outdir / "input.fasta"
            input_path.write_bytes(uploaded.getvalue())

            genotype_filter = [g.strip() for g in genotype_filter_raw.split(",") if g.strip()] or None

            argv = [
                "run",
                "--input", str(input_path),
                "--outdir", str(outdir),
                "--distance", distance,
                "--region", region,
                "--threads", str(threads),
                "--min-coverage", str(min_coverage),
                "--min-identity", str(min_identity),
                "--min-query-coverage", str(min_query_coverage),
                "--region-strategy", region_strategy,
                "--genotype-validation", genotype_validation,
                "--ambiguities", ambiguities,
                "--refs-dir", str(refs_dir),
                "--panel-fasta", str(panel_fasta),
                "--minimap2", minimap2_exe,
                "--mafft", mafft_exe,
                "--tn93", tn93_exe,
            ]
            if threshold is not None:
                argv.extend(["--threshold", str(threshold)])
            if verbose:
                argv.append("--verbose")
            if dry_run:
                argv.append("--dry-run")
            if genotype_filter:
                for genotype in genotype_filter:
                    argv.extend(["--genotype", genotype])
            if ncbi_email:
                argv.extend(["--ncbi-email", ncbi_email])
            if force_download:
                argv.append("--force-download")
            if keep_temp:
                argv.append("--keep-temp")
            if keep_paf:
                argv.append("--keep-paf")

            parser = hcv_cluster.build_parser()
            args = parser.parse_args(argv)

            log_buffer = io.StringIO()
            with st.spinner("Running pipeline..."):
                try:
                    with contextlib.redirect_stdout(log_buffer):
                        rc = args.func(args)
                except (RuntimeError, ValueError, hcv_cluster_prep.HcvPrepError) as exc:
                    rc = 1
                    log_buffer.write(f"ERROR: {exc}\n")

            with st.expander("Run log", expanded=rc != 0):
                st.code(log_buffer.getvalue() or "(no output)")

            if rc != 0:
                st.error(f"Pipeline failed (exit code {rc}). See run log above.")
            elif dry_run:
                st.info("Dry run complete — no files were written. See the commands above.")
            else:
                st.success(f"Done. Results written to: {outdir}")
                st.session_state.last_run = {"outdir": str(outdir), "distance": distance}

    if "last_run" in st.session_state:
        _render_results(Path(st.session_state.last_run["outdir"]), st.session_state.last_run["distance"])


def _render_results(outdir: Path, distance: str) -> None:
    st.subheader("Results")
    st.caption("clusters.csv + links.csv together are a MicrobeTrace-ready node/edge pair spanning all genotypes.")
    _show_cluster_table("Merged clusters (primary metric)", outdir / "clusters.csv")
    _show_cluster_table("Merged links (primary metric)", outdir / "links.csv")
    if distance == "both":
        _show_cluster_table("Merged clusters (SNP)", outdir / "clusters.snp.csv")
        _show_cluster_table("Merged links (SNP)", outdir / "links.snp.csv")

    by_genotype_dir = outdir / "by_genotype"
    genotypes = sorted(p.name for p in by_genotype_dir.iterdir() if p.is_dir()) if by_genotype_dir.exists() else []
    if not genotypes:
        return

    st.subheader("Cluster network")
    col_genotype, col_metric = st.columns(2)
    with col_genotype:
        selected_genotype = st.selectbox("Genotype", genotypes, key="view_genotype")
    with col_metric:
        if distance == "both":
            metric = st.selectbox("Distance metric", ("tn93", "snp"), key="view_metric")
        else:
            metric = distance
            st.write(f"Distance metric: **{metric.upper()}**")

    if st.button("View Clusters"):
        st.session_state.cluster_view = {
            "outdir": str(outdir),
            "genotype": selected_genotype,
            "metric": metric,
        }

    view = st.session_state.get("cluster_view")
    if not view or view["outdir"] != str(outdir):
        return

    genotype_dir = by_genotype_dir / view["genotype"]
    links_path, clusters_path = _metric_file_paths(genotype_dir, distance, view["metric"])
    if not (links_path.exists() and clusters_path.exists()):
        st.error(f"No {view['metric'].upper()} result files found for genotype {view['genotype']}.")
        return

    node_rows = _read_csv_rows(clusters_path)
    edge_rows = _read_csv_rows(links_path)
    if not node_rows:
        st.info("No sequences passed QC for this genotype.")
        return

    n_multi = len({row["cluster_id"] for row in node_rows if int(row["cluster_size"]) > 1})
    n_singleton = sum(1 for row in node_rows if int(row["cluster_size"]) == 1)

    hide_singletons = st.checkbox("Hide singletons", key="hide_singletons")
    plot_node_rows, plot_edge_rows = (
        hcv_cluster_viz.filter_singletons(node_rows, edge_rows) if hide_singletons else (node_rows, edge_rows)
    )
    if not plot_node_rows:
        st.info("No multi-member clusters for this genotype/metric.")
        return

    with st.spinner("Building network..."):
        fig = hcv_cluster_viz.build_network_figure(plot_node_rows, plot_edge_rows)
    st.plotly_chart(fig, width="stretch")

    try:
        png_bytes = fig.to_image(format="png", scale=2)
    except Exception as exc:
        st.caption(f"PNG export unavailable ({exc}); use the camera icon in the plot toolbar instead.")
    else:
        st.download_button(
            "Save image (PNG)",
            data=png_bytes,
            file_name=f"{view['genotype']}_{view['metric']}_network.png",
            mime="image/png",
        )

    st.caption(
        f"Genotype {view['genotype']} ({view['metric'].upper()}): "
        f"{len(node_rows)} sequence(s), {n_multi} cluster(s), {n_singleton} singleton(s)"
    )


def _refs_tab() -> None:
    st.header("Cache reference sequences")
    st.caption("Downloads and caches GenBank reference records used for alignment. Run once before the first clustering run.")

    genotype = st.text_input("Genotype to prepare, or 'all'", value="all")
    refs_dir = Path(st.text_input("References cache dir", value="refs", key="refs_tab_dir"))
    ncbi_email = st.text_input("NCBI email", value="", key="refs_tab_email")
    force_download = st.checkbox("Force re-download", key="refs_tab_force")
    verbose = st.checkbox("Verbose logging", key="refs_tab_verbose")

    if not st.button("Cache references"):
        return

    argv = ["refs", "--genotype", genotype, "--refs-dir", str(refs_dir)]
    if ncbi_email:
        argv.extend(["--ncbi-email", ncbi_email])
    if force_download:
        argv.append("--force-download")
    if verbose:
        argv.append("--verbose")

    parser = hcv_cluster.build_parser()
    args = parser.parse_args(argv)

    log_buffer = io.StringIO()
    with st.spinner("Caching references..."):
        try:
            with contextlib.redirect_stdout(log_buffer):
                rc = args.func(args)
        except (RuntimeError, ValueError, hcv_cluster_prep.HcvPrepError) as exc:
            rc = 1
            log_buffer.write(f"ERROR: {exc}\n")

    st.code(log_buffer.getvalue() or "(no output)")
    if rc == 0:
        st.success(f"References cached in: {refs_dir}")
    else:
        st.error(f"Failed (exit code {rc}).")


def _shutdown_server() -> None:
    def _terminate() -> None:
        time.sleep(0.75)
        os.kill(os.getpid(), signal.SIGTERM)

    threading.Thread(target=_terminate, daemon=True).start()


def _server_sidebar() -> None:
    with st.sidebar:
        st.header("Server")
        if not st.session_state.get("confirm_exit"):
            if st.button("Exit"):
                st.session_state.confirm_exit = True
                st.rerun()
        else:
            st.warning("Shut down the GUI server? This closes the app for all users.")
            col_yes, col_no = st.columns(2)
            if col_yes.button("Yes, shut down", type="primary"):
                st.info("Shutting down — you can close this browser tab.")
                _shutdown_server()
                st.stop()
            if col_no.button("Cancel"):
                st.session_state.confirm_exit = False
                st.rerun()


st.title("HCV Genomic Clustering")
_server_sidebar()

run_tab, refs_tab = st.tabs(["Run clustering", "Cache references"])
with run_tab:
    _run_tab()
with refs_tab:
    _refs_tab()
