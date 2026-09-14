#!/usr/bin/env python3
"""Streamlit GUI for HCV genomic clustering. Thin wrapper around hcv_cluster.command_run."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import math
import os
import signal
import tempfile
import threading
import time
from pathlib import Path
from typing import Any

import streamlit as st

import assign_hcv_genotypes_from_fasta
import hcv_cluster
import hcv_cluster_metadata
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


def _validate_metadata_upload(data: bytes) -> list[dict[str, Any]]:
    """Validate uploaded CSV bytes without writing into the selected result directory."""
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise hcv_cluster_metadata.MetadataValidationError(
            "Metadata CSV must be UTF-8 encoded"
        ) from exc
    try:
        reader = csv.DictReader(io.StringIO(text))
        original_fields = reader.fieldnames
        if original_fields is None:
            raise hcv_cluster_metadata.MetadataValidationError("Metadata CSV has no header row")
        cleaned_fields = [field.strip() for field in original_fields]
        source_rows = []
        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise hcv_cluster_metadata.MetadataValidationError(
                    f"Metadata row {row_number} has more values than the header"
                )
            source_rows.append(
                {
                    cleaned: row.get(original)
                    for original, cleaned in zip(original_fields, cleaned_fields, strict=True)
                }
            )
    except csv.Error as exc:
        raise hcv_cluster_metadata.MetadataValidationError(
            f"Could not parse metadata CSV: {exc}"
        ) from exc
    return hcv_cluster_metadata.normalize_metadata_rows(
        source_rows, fieldnames=cleaned_fields
    )


def _scope_metric_file_paths(
    outdir: Path, scope: str, distance: str, metric: str
) -> tuple[Path, Path]:
    """Resolve result files for either the merged graph or a genotype graph."""
    if scope == "All genotypes":
        if distance == "both" and metric == "snp":
            return outdir / "links.snp.csv", outdir / "clusters.snp.csv"
        return outdir / "links.csv", outdir / "clusters.csv"
    return _metric_file_paths(outdir / "by_genotype" / scope, distance, metric)


def _network_fingerprint(
    node_rows: list[dict[str, Any]], edge_rows: list[dict[str, Any]]
) -> str:
    """Identify graph topology for the session layout cache."""
    digest = hashlib.sha256()
    for row in sorted(node_rows, key=lambda item: str(item["sample_id"])):
        sample_id = str(row["sample_id"])
        genotype = str(row.get("genotype", ""))
        digest.update(f"n\0{sample_id}\0{genotype}\n".encode())
    for source, target in sorted(
        tuple(sorted((str(row["source"]), str(row["target"])))) for row in edge_rows
    ):
        digest.update(f"e\0{source}\0{target}\n".encode())
    return digest.hexdigest()


def _metadata_display_fields(node_rows: list[dict[str, Any]]) -> list[str]:
    excluded = {"sample_id", "cluster_id", "cluster_size", "source_cluster_id"}
    return [
        field
        for field in dict.fromkeys(field for row in node_rows for field in row)
        if field not in excluded
    ]


def _size_display_fields(
    fields: list[str], node_rows: list[dict[str, Any]]
) -> list[str]:
    """Return fields the renderer can size without an explicit category order."""
    result: list[str] = []
    for field in fields:
        values = [
            str(row.get(field, "")).strip()
            for row in node_rows
            if str(row.get(field, "")).strip()
            not in {"", hcv_cluster_metadata.MISSING_VALUE}
        ]
        if field.casefold().replace(" ", "_") == "age_range":
            result.append(field)
            continue
        try:
            if values and all(math.isfinite(float(value)) for value in values):
                result.append(field)
        except ValueError:
            pass
    return result


def _epidemiology_defaults(
    fields: list[str], node_rows: list[dict[str, Any]]
) -> dict[str, str | None]:
    """Choose conservative defaults only when the user requests the preset."""
    lowered = {field.casefold(): field for field in fields}
    ordered_size_fields = set(_size_display_fields(fields, node_rows))

    def first_named(*names: str) -> str | None:
        return next((lowered[name] for name in names if name in lowered), None)

    def low_cardinality(excluding: set[str]) -> str | None:
        for field in fields:
            if field in ordered_size_fields:
                continue
            values = {
                str(row.get(field) or hcv_cluster_metadata.MISSING_VALUE)
                for row in node_rows
            }
            values.discard(hcv_cluster_metadata.MISSING_VALUE)
            if field not in excluding and 1 < len(values) <= 8:
                return field
        return None

    color = first_named("location", "prison", "facility", "subtype", "genotype")
    shape = first_named("injecting_status", "injecting status", "indigenous_status", "indigenous status")
    if shape == color:
        shape = None
    shape = shape or low_cardinality({color} if color else set())
    size = first_named("age_range", "age range", "age")
    outline = first_named("indigenous_status", "indigenous status", "injecting_status", "injecting status")
    if outline in {color, shape}:
        outline = low_cardinality({value for value in (color, shape) if value})
    return {"color": color, "symbol": shape, "size": size, "outline": outline}


def _safe_plot_filename(value: str) -> str:
    cleaned = "".join(character if character.isalnum() or character in "-_" else "_" for character in value)
    return cleaned.strip("_") or "network"


def _figure_svg_bytes(figure: Any) -> bytes:
    """Use the shared vector-safety check and return data for Streamlit download."""
    with tempfile.TemporaryDirectory(prefix="hcv-cluster-svg-") as temp_dir:
        svg_path = Path(temp_dir) / "network.svg"
        hcv_cluster_viz.export_figure_svg(figure, svg_path)
        return svg_path.read_bytes()


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

    upload_col, metadata_col = st.columns(2)
    with upload_col:
        uploaded = st.file_uploader(
            "Input multi-FASTA of HCV consensus sequences", type=["fasta", "fa", "fna"]
        )
    with metadata_col:
        uploaded_metadata = st.file_uploader(
            "Optional sample metadata CSV",
            type=["csv"],
            help=(
                "Requires a sample_id column matching FASTA identifiers. Other columns "
                "are available for colour, shape, size, outline, and hover."
            ),
        )
    metadata_error: str | None = None
    if uploaded_metadata is not None:
        try:
            metadata_preview = _validate_metadata_upload(uploaded_metadata.getvalue())
        except hcv_cluster_metadata.MetadataValidationError as exc:
            metadata_error = str(exc)
            st.error(f"Metadata CSV is not valid: {metadata_error}")
        else:
            metadata_fields = hcv_cluster_metadata.metadata_fields(metadata_preview)
            st.success(
                f"Metadata ready: {len(metadata_preview)} sample(s), "
                f"{len(metadata_fields)} field(s)."
            )
    outdir = Path(st.text_input("Output directory", value="results_gui"))

    region, region_error = _region_picker()

    col1, col2, col3 = st.columns(3)
    with col1:
        distance = st.selectbox("Distance metric", ("tn93", "snp", "both"), index=0)
    with col2:
        default_threshold_tn93 = hcv_cluster.resolve_threshold(region, None, distance="tn93")
        default_threshold_snp = hcv_cluster.resolve_threshold(region, None, distance="snp")
        default_threshold = default_threshold_snp if distance == "snp" else default_threshold_tn93
        override_threshold = st.checkbox("Override default threshold")
        if override_threshold:
            threshold = st.number_input(
                "Threshold", min_value=0.0, max_value=1.0, value=default_threshold, step=0.001, format="%.4f"
            )
        else:
            threshold = None
            if distance == "both" and default_threshold_tn93 != default_threshold_snp:
                st.caption(
                    f"Default for '{region}': TN93 **{default_threshold_tn93:.4g}**, "
                    f"SNP **{default_threshold_snp:.4g}**"
                )
            else:
                metric_label = {"tn93": "TN93", "snp": "SNP", "both": "TN93 & SNP"}[distance]
                st.caption(f"Default for '{region}' ({metric_label}): **{default_threshold:.4g}**")
    with col3:
        threads = st.number_input("MAFFT threads", min_value=1, value=1, step=1)

    snp_count_threshold = None
    if distance in ("snp", "both"):
        use_snp_count = st.checkbox(
            "Link SNP pairs by SNP count instead of p-distance",
            help=(
                "Cluster on the raw number of differing bases (e.g. \"within 30 SNPs\") "
                "rather than the proportion. The SNP count is reported either way."
            ),
        )
        if use_snp_count:
            snp_count_threshold = int(
                st.number_input("Maximum SNPs between linked samples", min_value=0, value=30, step=1)
            )
            st.warning(
                "Pairs are compared over differing numbers of sites (N and gap positions "
                "are skipped), so the same SNP count means different divergence for "
                "different pairs, and no HCV clustering threshold in the literature is "
                "defined this way. Prefer the p-distance default for anything reportable."
            )

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
        reuse_alignments_from = st.text_input(
            "Reuse alignments from results directory",
            value="",
            help=(
                "Reuse unchanged samples from per-genotype prep.aligned.fasta files "
                "under this previous results directory, and send only new or changed "
                "samples to MAFFT. Leave blank to align every sample."
            ),
        )
        keep_paf = st.checkbox(
            "Keep raw minimap2 PAF file",
            help="Keep the raw minimap2 PAF alignment file from the genotyping step instead of deleting it — useful for debugging genotype assignment.",
        )

    run_clicked = st.button(
        "Run pipeline",
        type="primary",
        disabled=uploaded is None or region_error is not None or metadata_error is not None,
    )

    if run_clicked:
        if uploaded is None:
            st.error("Upload an input FASTA first.")
        elif region_error is not None:
            st.error(region_error)
        elif metadata_error is not None:
            st.error(f"Correct the metadata CSV before running: {metadata_error}")
        else:
            outdir.mkdir(parents=True, exist_ok=True)
            input_path = outdir / "input.fasta"
            input_path.write_bytes(uploaded.getvalue())
            metadata_path = None
            if uploaded_metadata is not None:
                metadata_path = outdir / "input_metadata.csv"
                metadata_path.write_bytes(uploaded_metadata.getvalue())

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
            if reuse_alignments_from.strip():
                argv.extend(["--reuse-alignments", reuse_alignments_from.strip()])
            if snp_count_threshold is not None:
                argv.extend(["--snp-count-threshold", str(snp_count_threshold)])
            if keep_paf:
                argv.append("--keep-paf")
            if metadata_path is not None:
                argv.extend(["--metadata", str(metadata_path)])

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
    has_merged_network = (outdir / "clusters.csv").exists() and (outdir / "links.csv").exists()
    # Preserve the original genotype-first GUI experience while making the merged
    # presentation available from the same selector.
    scopes = genotypes + (["All genotypes"] if has_merged_network else [])
    if not scopes:
        return

    st.subheader("Cluster network")
    st.caption(
        "The all-genotype view combines already-computed networks for presentation only. "
        "Clustering and genetic-distance comparisons remain genotype-stratified; no "
        "cross-genotype links are created."
    )
    col_scope, col_metric = st.columns(2)
    with col_scope:
        if st.session_state.get("view_scope") not in scopes:
            st.session_state.view_scope = scopes[0]
        selected_scope = st.selectbox("Network scope", scopes, key="view_scope")
    with col_metric:
        if distance == "both":
            metric = st.selectbox("Distance metric", ("tn93", "snp"), key="view_metric")
        else:
            metric = distance
            st.write(f"Distance metric: **{metric.upper()}**")

    selected_layout = "Spring layout"
    if selected_scope == "All genotypes":
        selected_layout = st.selectbox(
            "Combined layout",
            ("Packed clusters", "Grouped by genotype"),
            help=(
                "Both options lay out each connected cluster independently. Grouped mode "
                "places genotype networks in separate spatial regions."
            ),
            key="view_combined_layout",
        )

    if st.button("View / update network"):
        st.session_state.cluster_view = {
            "outdir": str(outdir),
            "scope": selected_scope,
            "metric": metric,
            "layout": selected_layout,
        }

    view = st.session_state.get("cluster_view")
    if not view or view["outdir"] != str(outdir):
        return

    # Old Streamlit sessions may still contain the pre-metadata genotype-only state.
    scope = view.get("scope", view.get("genotype"))
    if scope not in scopes:
        st.info("Choose an available network scope and click View / update network.")
        return
    layout_label = view.get("layout", "Spring layout")
    links_path, clusters_path = _scope_metric_file_paths(outdir, scope, distance, view["metric"])
    if not (links_path.exists() and clusters_path.exists()):
        st.error(f"No {view['metric'].upper()} result files found for {scope}.")
        return

    node_rows = _read_csv_rows(clusters_path)
    edge_rows = _read_csv_rows(links_path)
    if not node_rows:
        st.info("No sequences passed QC for this network scope.")
        return

    if scope != "All genotypes":
        for row in node_rows:
            row.setdefault("genotype", scope)

    metadata_rows: list[dict[str, Any]] = []
    metadata_path = outdir / "metadata.csv"
    if metadata_path.exists():
        try:
            metadata_rows = hcv_cluster_metadata.load_metadata_csv(metadata_path)
            node_rows, join_report = hcv_cluster_metadata.join_metadata(node_rows, metadata_rows)
        except hcv_cluster_metadata.MetadataValidationError as exc:
            st.error(
                f"Could not apply saved metadata ({metadata_path}): {exc}. "
                "Correct the source CSV and rerun clustering."
            )
            return
        if join_report.missing_metadata_sample_ids:
            preview = ", ".join(join_report.missing_metadata_sample_ids[:8])
            suffix = "…" if len(join_report.missing_metadata_sample_ids) > 8 else ""
            st.warning(
                f"{len(join_report.missing_metadata_sample_ids)} plotted sample(s) have no "
                f"metadata and are shown as {hcv_cluster_metadata.MISSING_VALUE}: {preview}{suffix}"
            )
        if scope == "All genotypes" and join_report.unknown_metadata_sample_ids:
            preview = ", ".join(join_report.unknown_metadata_sample_ids[:8])
            suffix = "…" if len(join_report.unknown_metadata_sample_ids) > 8 else ""
            st.info(
                f"{len(join_report.unknown_metadata_sample_ids)} metadata sample(s) are not "
                f"in this plotted result (for example, they may have failed QC): {preview}{suffix}"
            )

    n_multi = len({row["cluster_id"] for row in node_rows if int(row["cluster_size"]) > 1})
    n_singleton = sum(1 for row in node_rows if int(row["cluster_size"]) == 1)

    hide_singletons = st.checkbox("Hide singletons", key="hide_singletons")
    plot_node_rows, plot_edge_rows = (
        hcv_cluster_viz.filter_singletons(node_rows, edge_rows) if hide_singletons else (node_rows, edge_rows)
    )
    if not plot_node_rows:
        st.info("No multi-member clusters for this scope/metric.")
        return

    fields = _metadata_display_fields(plot_node_rows)
    st.markdown("#### Node appearance and hover")
    preset_col, reset_col, preset_help_col = st.columns([1, 1, 2])
    with preset_col:
        apply_epidemiology = st.button(
            "Apply epidemiology view",
            help="Choose sensible available metadata fields across several visual channels.",
        )
    with reset_col:
        reset_cluster = st.button(
            "Reset to cluster view",
            help="Restore cluster colours and remove all other node encodings.",
        )
    with preset_help_col:
        st.caption(
            "You can then change any dropdown independently. These controls only restyle "
            "the cached graph; they do not rerun clustering or move nodes."
        )

    size_fields = _size_display_fields(fields, plot_node_rows)
    channel_keys = {
        "color": "network_color_by",
        "symbol": "network_symbol_by",
        "size": "network_size_by",
        "outline": "network_outline_by",
    }
    if reset_cluster:
        for key in channel_keys.values():
            st.session_state[key] = None
    elif apply_epidemiology:
        defaults = _epidemiology_defaults(fields, plot_node_rows)
        for channel, key in channel_keys.items():
            st.session_state[key] = defaults[channel]
    for channel, key in channel_keys.items():
        options = size_fields if channel == "size" else fields
        if st.session_state.get(key) not in [None, *options]:
            st.session_state[key] = None

    appearance_columns = st.columns(4)
    with appearance_columns[0]:
        color_by = st.selectbox(
            "Node colour",
            [None, *fields],
            format_func=lambda value: "Cluster (default)" if value is None else value,
            key=channel_keys["color"],
        )
    with appearance_columns[1]:
        symbol_by = st.selectbox(
            "Node shape",
            [None, *fields],
            format_func=lambda value: "None" if value is None else value,
            key=channel_keys["symbol"],
        )
    with appearance_columns[2]:
        size_by = st.selectbox(
            "Node size",
            [None, *size_fields],
            format_func=lambda value: "None" if value is None else value,
            key=channel_keys["size"],
            help="Numeric fields and age_range are available for meaningful ordered sizes.",
        )
    with appearance_columns[3]:
        outline_by = st.selectbox(
            "Node outline",
            [None, *fields],
            format_func=lambda value: "None" if value is None else value,
            key=channel_keys["outline"],
        )

    hover_key = "network_hover_fields"
    if hover_key not in st.session_state:
        st.session_state[hover_key] = fields.copy()
    else:
        st.session_state[hover_key] = [
            field for field in st.session_state[hover_key] if field in fields
        ]
    hover_fields = st.multiselect(
        "Additional hover fields",
        fields,
        key=hover_key,
        help="Sample ID, cluster ID, and cluster size are always included.",
    )

    for warning in hcv_cluster_viz.get_encoding_warnings(
        plot_node_rows,
        color_by=color_by,
        symbol_by=symbol_by,
        size_by=size_by,
        outline_by=outline_by,
    ):
        st.warning(warning)

    # Build maps from the complete saved metadata/root-node universe, not only the
    # currently selected genotype. This keeps (for example) Prison A the same colour
    # while users switch repeatedly between combined and per-genotype views.
    metadata_field_names = set(hcv_cluster_metadata.metadata_fields(metadata_rows))
    all_genotype_rows = (
        _read_csv_rows(outdir / "clusters.csv")
        if (outdir / "clusters.csv").exists()
        else []
    )
    encoding_maps: dict[str, dict[str, str]] = {}
    for channel, field in (
        ("color", color_by),
        ("symbol", symbol_by),
        ("outline", outline_by),
    ):
        if field is None:
            continue
        values = [row.get(field) for row in plot_node_rows]
        if field in metadata_field_names:
            values.extend(row.get(field) for row in metadata_rows)
        if field == "genotype":
            values.extend(row.get("genotype") for row in all_genotype_rows)
        encoding_maps[f"{channel}:{field}"] = hcv_cluster_viz.build_category_mapping(
            values, channel=channel
        )

    if scope == "All genotypes":
        layout_mode = "components"
        group_by = "genotype" if layout_label == "Grouped by genotype" else None
    else:
        layout_mode = "spring"
        group_by = None
    cache_key = (
        str(outdir.resolve()),
        scope,
        view["metric"],
        layout_label,
        hide_singletons,
        _network_fingerprint(plot_node_rows, plot_edge_rows),
    )
    layout_cache = st.session_state.setdefault("network_layout_cache", {})
    if cache_key not in layout_cache:
        layout_cache[cache_key] = hcv_cluster_viz.compute_network_layout(
            plot_node_rows,
            plot_edge_rows,
            mode=layout_mode,
            group_by=group_by,
        )
        # Avoid retaining layouts for an unbounded number of result directories.
        if len(layout_cache) > 24:
            oldest_key = next(iter(layout_cache))
            if oldest_key != cache_key:
                layout_cache.pop(oldest_key)

    with st.spinner("Building network..."):
        try:
            fig = hcv_cluster_viz.build_network_figure(
                plot_node_rows,
                plot_edge_rows,
                positions=layout_cache[cache_key],
                group_by=group_by,
                color_by=color_by,
                symbol_by=symbol_by,
                size_by=size_by,
                outline_by=outline_by,
                hover_fields=hover_fields,
                encoding_maps=encoding_maps,
            )
        except ValueError as exc:
            st.error(f"Cannot apply the selected network appearance: {exc}")
            return
    st.plotly_chart(fig, width="stretch")

    scope_filename = "all_genotypes" if scope == "All genotypes" else scope
    filename_base = _safe_plot_filename(f"{scope_filename}_{view['metric']}_network")
    download_columns = st.columns(2)
    try:
        png_bytes = fig.to_image(format="png", scale=2)
    except Exception as exc:
        st.caption(f"PNG export unavailable ({exc}); use the camera icon in the plot toolbar instead.")
    else:
        with download_columns[0]:
            st.download_button(
                "Save image (PNG)",
                data=png_bytes,
                file_name=f"{filename_base}.png",
                mime="image/png",
            )
    try:
        svg_bytes = _figure_svg_bytes(fig)
    except Exception as exc:
        st.caption(f"SVG export unavailable ({exc}).")
    else:
        with download_columns[1]:
            st.download_button(
                "Save editable vector (SVG)",
                data=svg_bytes,
                file_name=f"{filename_base}.svg",
                mime="image/svg+xml",
            )

    scope_label = "All genotypes" if scope == "All genotypes" else f"Genotype {scope}"
    st.caption(
        f"{scope_label} ({view['metric'].upper()}): "
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
