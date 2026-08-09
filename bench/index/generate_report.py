#!/usr/bin/env python3
"""
z-fasta Index Benchmark Report Generator

Reads zebrac JSON joined with metadata JSONL from bench/index/results/,
produces Markdown report + PNG figures using pandas + matplotlib.

Usage:
    .venv/bin/python bench/index/generate_report.py [results_dir]

Defaults to bench/index/results/ (latest run via results/LATEST).
"""

import json
import argparse
from pathlib import Path
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import pandas as pd
import tabulate  # noqa: F401  -- needed by pd.to_markdown()

MARKDOWNLINT_DISABLE = "<!-- markdownlint-disable MD024 MD032 MD033 MD036 MD041 MD049 -->"


class ReportCounters:
    """Global table and figure numbers across the report."""

    def __init__(self) -> None:
        self.table = 1
        self.figure = 1

    def next_table(self) -> int:
        n = self.table
        self.table += 1
        return n

    def next_figure(self) -> int:
        n = self.figure
        self.figure += 1
        return n

# Styling.
COLORS = {
    "z-fasta-fai": "#F7A41D",  # Zig gold (peer-comparable)
    "z-fasta-fai-nodedup": "#FBC02D",
    "z-fasta-zfi": "#E65100",  # deep orange (first-class)
    "z-fasta-zfi-nodedup": "#EF6C00",
    "noodles": "#C45C26",  # rust bronze
    "rustbio-custom-index": "#8B3A2A",  # rust-bio benchmark lane
    "samtools": "#555555",  # C grey (GitHub linguist)
    "seqkit": "#009485",  # Go teal
    "fastahack": "#F34B7D",  # C++ pink (GitHub linguist)
    "pyfaidx": "#3572A5",  # Python blue (GitHub linguist)
}
DISPLAY_NAMES = {
    "z-fasta-fai": "z-fasta (.fai)",
    "z-fasta-fai-nodedup": "z-fasta (.fai) --no-dedup",
    "z-fasta-zfi": "z-fasta (.zfi)",
    "z-fasta-zfi-nodedup": "z-fasta (.zfi) --no-dedup",
    "samtools": "samtools",
    "seqkit": "seqkit",
    "fastahack": "fastahack",
    "pyfaidx": "pyfaidx",
    "noodles": "noodles",
    "rustbio-custom-index": "rust-bio",
}
EDGE_EXIT_COLUMNS = (
    "zfasta_exit",
    "zfasta_zfi_exit",
    "noodles_exit",
    "rustbio_exit",
    "samtools_exit",
    "seqkit_exit",
    "fastahack_exit",
    "pyfaidx_exit",
)
EDGE_EXIT_LABELS = {
    "zfasta_exit": "z-fasta (.fai)",
    "zfasta_zfi_exit": "z-fasta (.zfi)",
    "samtools_exit": "samtools",
    "seqkit_exit": "seqkit",
    "fastahack_exit": "fastahack",
    "pyfaidx_exit": "pyfaidx",
    "noodles_exit": "noodles",
    "rustbio_exit": "rust-bio",
}
EDGE_CONTRACT_CLASSES = {
    "fai_parity": "FAI parity",
    "zfi_messy": "ZFI messy support",
    "invalid_input": "Invalid input review",
}
EDGE_MESSY_CASES = {"uniform", "mixed_widths", "trailing_whitespace", "blank_lines", "mixed_crlf"}
DATASET_ORDER = ["Genome", "Transcriptome", "Proteome"]
# Real benchmark datasets (bench/shared/data; see bench/shared/download_data.sh).
REAL_DATASETS = [
    {
        "id": "Genome",
        "description": "Homo sapiens GRCh38 primary assembly (Ensembl release 113)",
    },
    {
        "id": "Transcriptome",
        "description": "Homo sapiens GENCODE v46 transcript sequences",
    },
    {
        "id": "Proteome",
        "description": "Homo sapiens UniProt reference proteome UP000005640",
    },
]

# Product section: format and dedup matrix on the record-rich real datasets.
PRODUCT_WORKLOADS = ("Transcriptome", "Proteome")
PRODUCT_COMPARISON_TOOLS = [
    "z-fasta-fai",
    "z-fasta-fai-nodedup",
    "z-fasta-zfi",
    "z-fasta-zfi-nodedup",
]

# Headline real-dataset chart: peer-comparable z-fasta (.fai) plus competitors.
HEADLINE_PERF_TOOLS = [
    "z-fasta-fai",
    "noodles",
    "rustbio-custom-index",
    "samtools",
    "seqkit",
    "fastahack",
    "pyfaidx",
]


def dataset_sort_key(name):
    try:
        return DATASET_ORDER.index(name), str(name)
    except ValueError:
        return len(DATASET_ORDER), str(name)


def display_tool(name):
    return DISPLAY_NAMES.get(name, name)


def edge_exit_columns(df: pd.DataFrame) -> list[str]:
    return [column for column in EDGE_EXIT_COLUMNS if column in df.columns]


def edge_exit_label(column: str) -> str:
    return EDGE_EXIT_LABELS.get(column, column)


def edge_contract_classes(df: pd.DataFrame) -> pd.Series:
    if "contract_class" in df.columns:
        return df["contract_class"]

    def classify(case: str) -> str:
        if case == "binary_data":
            return "invalid_input"
        if case in EDGE_MESSY_CASES and case != "uniform":
            return "zfi_messy"
        return "fai_parity"

    return df["test_case"].map(classify)


def pivot_with_display_names(pivot: pd.DataFrame) -> pd.DataFrame:
    out = pivot.copy()
    out.columns = [display_tool(str(c)) for c in out.columns]
    return out


def filter_tools(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    return df[df["tool"].isin(tools)].copy()


def _ordered_pivot(
    pivot: pd.DataFrame,
    tools: list[str],
    *,
    sort_index=None,
) -> pd.DataFrame:
    cols = [c for c in tools if c in pivot.columns]
    pivot = pivot[cols]
    if sort_index is not None:
        pivot = pivot.reindex(sorted(pivot.index, key=sort_index))
    else:
        pivot = pivot.sort_index()
    return pivot_with_display_names(pivot)


def md_tool_pivot_table(
    df: pd.DataFrame,
    tools: list[str],
    index_col: str,
    value_col: str,
    fmt,
    *,
    sort_index=dataset_sort_key,
    empty_msg: str = "_No data._",
    require_non_null: bool = False,
) -> str:
    work = filter_tools(df, tools)
    if require_non_null and value_col in work.columns:
        work = work[work[value_col].notna()].copy()
    if work.empty:
        return empty_msg
    work = work.copy()
    work["cell"] = work.apply(fmt, axis=1)
    pivot = work.pivot(index=index_col, columns="tool", values="cell")
    return _ordered_pivot(pivot, tools, sort_index=sort_index).to_markdown()


# Data loading.


def load_metadata(results_dir: Path, manifest: dict | None = None) -> pd.DataFrame | None:
    """Load metadata JSONL for the active run (manifest-linked when available)."""
    paths: list[Path] = []
    if manifest is not None:
        metadata_name = manifest.get("metadata")
        if not isinstance(metadata_name, str):
            return None
        candidate = results_dir / metadata_name
        if not candidate.is_file():
            return None
        paths = [candidate]
    else:
        latest = results_dir / "LATEST"
        if latest.exists():
            ts = latest.read_text().strip()
            candidate = results_dir / f"metadata_{ts}.jsonl"
            if candidate.exists():
                paths = [candidate]
    if manifest is None and not paths:
        found = sorted(results_dir.glob("metadata_*.jsonl"), reverse=True)
        if found:
            paths = [found[0]]

    rows = []
    for path in paths:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                rows.append(row)
    if not rows:
        return None
    return pd.DataFrame(rows)


def meta_rows_for_json(metadata_df: pd.DataFrame | None, path: Path) -> list[dict]:
    """Rows from metadata JSONL whose raw_json matches path."""
    if metadata_df is None or metadata_df.empty or "raw_json" not in metadata_df.columns:
        return []
    resolved = path.resolve()
    rows = []
    for _, row in metadata_df.iterrows():
        raw = Path(row["raw_json"])
        try:
            if raw.resolve() == resolved:
                rows.append(row.to_dict())
        except OSError:
            if raw == path:
                rows.append(row.to_dict())
    return rows


def load_zebrac_json(path: Path, metadata_df: pd.DataFrame | None = None) -> pd.DataFrame:
    with open(path) as f:
        data = json.load(f)
    if "zebrac_version" not in data:
        raise ValueError(f"expected zebrac JSON, got unsupported format: {path}")

    meta_rows = meta_rows_for_json(metadata_df, path)
    meta_by_command = {row.get("command"): row for row in meta_rows}
    rows = []
    for idx, result in enumerate(data.get("results", [])):
        command = result.get("command", "")
        meta = meta_by_command.get(command, {})
        if not meta and idx < len(meta_rows):
            meta = meta_rows[idx]
        wall = result.get("wall_time", {})
        peak_rss = result.get("peak_rss", {})
        mean_s = wall.get("mean", 0) / 1_000_000_000.0
        input_bytes = meta.get("input_bytes")
        throughput_mibs = None
        input_mib = None
        if input_bytes:
            input_mib = float(input_bytes) / (1024 * 1024)
            if mean_s > 0:
                throughput_mibs = input_mib / mean_s
        rows.append(
            {
                "tool": meta.get("tool", command),
                "mean": mean_s,
                "stddev": wall.get("std_dev", 0) / 1_000_000_000.0,
                "median": wall.get("median", 0) / 1_000_000_000.0,
                "min": wall.get("min", 0) / 1_000_000_000.0,
                "max": wall.get("max", 0) / 1_000_000_000.0,
                "peak_rss_mb": peak_rss.get("mean", 0) / (1024 * 1024),
                "peak_rss_stddev_mb": peak_rss.get("std_dev", 0) / (1024 * 1024),
                "peak_rss_min_mb": peak_rss.get("min", 0) / (1024 * 1024),
                "peak_rss_max_mb": peak_rss.get("max", 0) / (1024 * 1024),
                "minor_faults": result.get("minor_faults", {}).get("mean", 0),
                "major_faults": result.get("major_faults", {}).get("mean", 0),
                "input_mib": input_mib,
                "throughput_mibs": throughput_mibs,
            }
        )
    return pd.DataFrame(rows)


def discover_latest(results_dir: Path, prefix: str, manifest: dict | None = None) -> Path | None:
    """Find timestamped directory or file for the active run.

    When a manifest is present, never mix in files from another run.
    """
    section_key = {
        "perf": "real",
        "scale_size": "scale_size",
        "scale_seqs_budget": "scale_seqs_budget",
        "scale_seqs_fixed": "scale_seqs_fixed",
        "tests": "correctness",
    }.get(prefix)
    if manifest is not None and section_key:
        rel = manifest.get("sections", {}).get(section_key)
        if not rel:
            if section_key != "correctness" or manifest.get("schema_version") == "index-run.v3":
                return None
        else:
            candidate = results_dir / rel
            return candidate if candidate.exists() else None

    # Unmapped prefixes (e.g. tests_*.csv): pin to the run timestamp.
    if manifest is not None:
        ts = manifest.get("timestamp")
        if not ts:
            return None
        for candidate in (
            results_dir / f"{prefix}_{ts}",
            results_dir / f"{prefix}_{ts}.csv",
        ):
            if candidate.is_dir() or candidate.is_file():
                return candidate
        return None

    dirs = sorted(results_dir.glob(f"{prefix}_*"), reverse=True)
    dirs = [d for d in dirs if d.is_dir()]
    if dirs:
        return dirs[0]

    csvs = sorted(results_dir.glob(f"{prefix}_*.csv"), reverse=True)
    if csvs:
        return csvs[0]
    return None


def load_run_manifest(
    results_dir: Path, manifest_path: Path | None = None
) -> dict | None:
    """Load an explicit manifest, then LATEST, then the newest timestamped run."""
    if manifest_path is not None:
        return json.loads(manifest_path.read_text())
    latest = results_dir / "LATEST"
    if latest.is_file():
        ts = latest.read_text().strip()
        path = results_dir / f"run_{ts}.json"
        if path.is_file():
            return json.loads(path.read_text())
    paths = sorted(results_dir.glob("run_*.json"), reverse=True)
    if paths:
        return json.loads(paths[0].read_text())
    return None


def prune_stale_pngs(figures_dir: Path, allowed) -> list[str]:
    """Delete PNGs not in allowed; return removed filenames."""
    allowed_set = set(allowed)
    removed: list[str] = []
    if not figures_dir.is_dir():
        return removed
    for path in figures_dir.glob("*.png"):
        if path.name in allowed_set:
            continue
        path.unlink()
        removed.append(path.name)
    return removed


def enrich_manifest(manifest: dict | None) -> dict | None:
    """Apply non-evidentiary defaults without probing the live environment."""
    if not manifest:
        return manifest

    manifest.setdefault("runner", "zebrac")
    manifest.setdefault("mode", "warm")
    return manifest


def load_json_suite(
    results_dir: Path,
    prefix: str,
    metadata_df: pd.DataFrame | None,
    manifest: dict | None,
    *,
    param_col: str | None = None,
    dataset_from_stem: bool = False,
) -> pd.DataFrame | None:
    """Load a directory of zebrac JSON files into one DataFrame."""
    d = discover_latest(results_dir, prefix, manifest)
    if not d:
        return None

    frames = []
    for jf in sorted(d.glob("*.json")):
        frame = load_zebrac_json(jf, metadata_df)
        if dataset_from_stem:
            frame["dataset"] = jf.stem
        elif param_col:
            stem = jf.stem
            try:
                val = float(stem.replace("mb", ""))
            except ValueError as exc:
                raise ValueError(
                    f"cannot derive {param_col} from result filename: {jf.name}"
                ) from exc
            frame[param_col] = val
        frames.append(frame)

    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_perf(
    results_dir: Path,
    metadata_df: pd.DataFrame | None = None,
    manifest: dict | None = None,
) -> pd.DataFrame | None:
    """Load real-dataset perf results (dataset column from filename stem)."""
    return load_json_suite(
        results_dir,
        "perf",
        metadata_df,
        manifest,
        dataset_from_stem=True,
    )


def load_scaling(
    results_dir: Path,
    prefix: str,
    param_col: str,
    metadata_df: pd.DataFrame | None = None,
    manifest: dict | None = None,
) -> pd.DataFrame | None:
    """Load scale_size or scale_seqs directory of JSON files."""
    return load_json_suite(
        results_dir,
        prefix,
        metadata_df,
        manifest,
        param_col=param_col,
    )


def load_tests(results_dir: Path, manifest: dict | None = None) -> pd.DataFrame | None:
    """Load edge-case test CSV."""
    p = discover_latest(results_dir, "tests", manifest)
    if not p:
        return None
    return pd.read_csv(p)


# Figures.


def _save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def filter_headline_perf(df: pd.DataFrame) -> pd.DataFrame:
    return filter_tools(df, HEADLINE_PERF_TOOLS)


def product_workloads(manifest: dict | None) -> tuple[str, ...]:
    configured = (manifest or {}).get("product_workloads")
    if (
        isinstance(configured, list)
        and configured
        and all(isinstance(workload, str) for workload in configured)
    ):
        return tuple(configured)
    return PRODUCT_WORKLOADS


def filter_product_comparison(
    df: pd.DataFrame,
    manifest: dict | None = None,
) -> pd.DataFrame:
    work = filter_tools(df, PRODUCT_COMPARISON_TOOLS)
    return work[work["dataset"].isin(product_workloads(manifest))].copy()


def _scaling_param_label(param_col: str, value: float) -> str:
    v = int(value)
    if param_col == "size_mb":
        return f"{v} MiB"
    return f"{v:,}"


def build_ratio_comparisons(
    work: pd.DataFrame,
    value_col: str,
    *,
    baseline: str,
    peer_tools: list[str],
    group_col: str = "dataset",
    group_sort=None,
    label_group=None,
) -> pd.DataFrame:
    """Baseline vs each peer; ratio = peer value / baseline value.

    ``baseline`` is required: suites use different tool ids (index: ``z-fasta-fai``,
    get: ``z-fasta-default``, stats: ``z-fasta-full``).
    """
    if group_sort is None:
        group_sort = dataset_sort_key if group_col == "dataset" else None
    groups = sorted(work[group_col].unique(), key=group_sort) if group_sort else sorted(
        work[group_col].unique()
    )
    rows: list[dict] = []
    for group in groups:
        base = work[(work[group_col] == group) & (work["tool"] == baseline)]
        if base.empty:
            continue
        base_v = float(base[value_col].values[0])
        group_label = label_group(group) if label_group else group
        for tool in peer_tools:
            if tool == baseline:
                continue
            hit = work[(work[group_col] == group) & (work["tool"] == tool)]
            if hit.empty:
                continue
            peer_v = float(hit[value_col].values[0])
            ratio = peer_v / base_v if base_v > 0 else None
            rows.append(
                {
                    "dataset": group_label,
                    "tool": tool,
                    "competitor": display_tool(tool),
                    "zfasta_v": base_v,
                    "comp_v": peer_v,
                    "ratio": ratio,
                }
            )
    return pd.DataFrame(rows)


def build_headline_ratio_comparisons(work: pd.DataFrame, value_col: str) -> pd.DataFrame:
    return build_ratio_comparisons(
        work,
        value_col,
        baseline="z-fasta-fai",
        peer_tools=HEADLINE_PERF_TOOLS[1:],
    )


def build_scaling_ratio_comparisons(
    work: pd.DataFrame,
    param_col: str,
    value_col: str = "mean",
) -> pd.DataFrame:
    return build_ratio_comparisons(
        work,
        value_col,
        baseline="z-fasta-fai",
        peer_tools=HEADLINE_PERF_TOOLS[1:],
        group_col=param_col,
        label_group=lambda p: _scaling_param_label(param_col, p),
    )


def build_time_throughput_comparisons(
    work: pd.DataFrame,
    ratio_df: pd.DataFrame,
    *,
    baseline: str,
) -> pd.DataFrame:
    """Wall-time ratio rows enriched with throughput columns."""
    if ratio_df.empty:
        return ratio_df
    rows: list[dict] = []
    for row in ratio_df.itertuples(index=False):
        zf = work[(work["dataset"] == row.dataset) & (work["tool"] == baseline)]
        comp = work[(work["dataset"] == row.dataset) & (work["tool"] == row.tool)]
        zf_tp = (
            float(zf["throughput_mibs"].values[0])
            if not zf.empty and pd.notna(zf["throughput_mibs"].values[0])
            else None
        )
        comp_tp = (
            float(comp["throughput_mibs"].values[0])
            if not comp.empty and pd.notna(comp["throughput_mibs"].values[0])
            else None
        )
        tp_ratio = zf_tp / comp_tp if zf_tp and comp_tp and comp_tp > 0 else None
        rows.append(
            {
                "dataset": row.dataset,
                "tool": row.tool,
                "competitor": row.competitor,
                "zfasta_s": row.zfasta_v,
                "comp_s": row.comp_v,
                "speedup": row.ratio,
                "ratio": row.ratio,
                "zfasta_mibs": zf_tp,
                "comp_mibs": comp_tp,
                "throughput_x": tp_ratio,
            }
        )
    return pd.DataFrame(rows)


def build_headline_comparisons(work: pd.DataFrame) -> pd.DataFrame:
    return build_time_throughput_comparisons(
        work,
        build_headline_ratio_comparisons(work, "mean"),
        baseline="z-fasta-fai",
    )


def _format_speedup(ratio: float | None) -> str:
    if ratio is None or ratio <= 0:
        return "n/a"
    if ratio >= 100:
        return f"{ratio:.0f}x"
    if ratio >= 10:
        return f"{ratio:.1f}x"
    if ratio >= 2:
        return f"{ratio:.1f}x"
    if ratio >= 0.1:
        return f"{ratio:.2f}x"
    if ratio >= 0.01:
        return f"{ratio:.3f}x"
    return f"{ratio:.4f}x"


def _annotate_headline_comparisons(
    ax1,
    datasets: list[str],
    tools: list[str],
    width: float,
    bar_tops: dict[tuple[str, str], tuple[float, float]],
    comparisons: pd.DataFrame,
) -> None:
    """Speedup badges (z-fasta vs competitor) and per-dataset reference lines."""
    zf_color = COLORS.get("z-fasta-fai", "#F7A41D")
    badge_offsets: dict[tuple[str, str], int] = {}

    for ds in datasets:
        previous_y: float | None = None
        level = 0
        for tool in tools:
            key = (ds, tool)
            if key not in bar_tops:
                continue
            current_y = bar_tops[key][1]
            close = (
                previous_y is not None
                and max(previous_y, current_y) / min(previous_y, current_y) <= 1.2
            )
            level = 1 - level if close else 0
            badge_offsets[key] = 9 + level * 18
            previous_y = current_y

    for ds in datasets:
        zf_key = (ds, "z-fasta-fai")
        if zf_key not in bar_tops:
            continue
        zf_x, zf_y = bar_tops[zf_key]
        cluster_xs = [bar_tops[(ds, t)][0] for t in tools if (ds, t) in bar_tops]
        if not cluster_xs:
            continue
        span_lo = min(cluster_xs) - width * 0.65
        span_hi = max(cluster_xs) + width * 0.65
        ax1.hlines(
            zf_y,
            span_lo,
            span_hi,
            colors=zf_color,
            linestyles=(0, (4, 3)),
            linewidth=1.1,
            alpha=0.55,
            zorder=1,
        )
        ax1.annotate(
            "1x",
            (zf_x, zf_y),
            xytext=(0, badge_offsets[zf_key]),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=zf_color,
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor=zf_color,
                linewidth=1.1,
            ),
            zorder=5,
        )

    for row in comparisons.itertuples(index=False):
        key = (row.dataset, row.tool)
        if key not in bar_tops:
            continue
        xpos, mean_t = bar_tops[key]
        color = COLORS.get(row.tool, "#666666")
        ratio = getattr(row, "ratio", None)
        if ratio is None:
            ratio = getattr(row, "speedup", None)
        ax1.annotate(
            _format_speedup(ratio),
            (xpos, mean_t),
            xytext=(0, badge_offsets[key]),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            linespacing=0.95,
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor=color,
                linewidth=1.1,
                alpha=0.96,
            ),
            zorder=5,
        )


def _fig_throughput_dual_axis(
    work: pd.DataFrame,
    out: Path,
    *,
    tools: list[str],
    comparisons: pd.DataFrame,
    title: str,
    fig_note: str,
    bar_width: float,
) -> Path:
    """Grouped wall-time bars (log) with throughput lines on a twin axis."""
    datasets = sorted(work["dataset"].unique(), key=dataset_sort_key)
    tools = [t for t in tools if t in work["tool"].unique()]

    fig, ax1 = plt.subplots(figsize=(14, 7.2))
    ax2 = ax1.twinx()

    width = min(bar_width, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for i, tool in enumerate(tools):
        color = COLORS.get(tool, "#888888")
        xs: list[float] = []
        tputs: list[float] = []
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            mean_t = max(float(row["mean"].values[0]), 1e-6)
            std_t = max(float(row["stddev"].values[0]), 0.0)
            xpos = di + (i - len(tools) / 2 + 0.5) * width
            xs.append(xpos)
            ax1.bar(
                xpos,
                mean_t,
                width,
                color=color,
                alpha=0.88,
                yerr=std_t,
                capsize=2,
                error_kw={"elinewidth": 0.8, "capthick": 0.8, "ecolor": "#333333"},
                zorder=2,
            )
            bar_tops[(ds, tool)] = (xpos, mean_t)
            tput = row["throughput_mibs"].values[0]
            if pd.notna(tput) and float(tput) > 0:
                tputs.append(float(tput))

        if xs and tputs and len(xs) == len(tputs):
            ax2.plot(
                xs,
                tputs,
                color=color,
                marker="o",
                linestyle="-",
                linewidth=1.8,
                markersize=6,
                markeredgecolor="white",
                markeredgewidth=0.7,
                label=display_tool(tool),
                zorder=3,
            )

    if not comparisons.empty:
        _annotate_headline_comparisons(ax1, datasets, tools, width, bar_tops, comparisons)

    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets, fontsize=11)
    ax1.set_ylabel("Wall Time (s, log scale)", fontsize=10)
    ax2.set_ylabel("Throughput (MiB/s)", fontsize=10)
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.935,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    ax1.grid(axis="y", alpha=0.28, which="both")
    ax1.set_axisbelow(True)

    handles, _labels = ax2.get_legend_handles_labels()
    fig.legend(
        handles,
        _labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=len(tools),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.88)
    return _save(fig, out)


def fig_performance_throughput(df: pd.DataFrame, out: Path) -> Path:
    """Grouped wall-time bars (left, log) with throughput lines/markers (right)."""
    work = filter_headline_perf(df)
    return _fig_throughput_dual_axis(
        work,
        out,
        tools=HEADLINE_PERF_TOOLS,
        comparisons=build_headline_comparisons(work),
        title="Real-Dataset Indexing: Wall Time & Throughput",
        fig_note=(
            "Bar labels: z-fasta speedup vs each tool "
            "(competitor wall time / z-fasta wall time)"
        ),
        bar_width=0.095,
    )


def fig_scaling_headline_lines(
    df: pd.DataFrame,
    out: Path,
    *,
    param_col: str,
    xlabel: str,
    title: str,
    fig_note: str,
) -> Path:
    """Line chart: wall time vs scaling parameter (headline tools, no ratio badges)."""
    work = filter_headline_perf(df)
    tools = [t for t in HEADLINE_PERF_TOOLS if t in work["tool"].unique()]

    fig, ax = plt.subplots(figsize=(14, 7.2))

    for tool in tools:
        tdf = work[work["tool"] == tool].sort_values(param_col)
        if tdf.empty:
            continue
        xs = tdf[param_col].astype(float)
        ys = tdf["mean"].astype(float).clip(lower=1e-6)
        color = COLORS.get(tool, "#888888")
        linewidth = 2.4 if tool == "z-fasta-fai" else 1.8
        ax.plot(
            xs,
            ys,
            color=color,
            marker="o",
            linestyle="-",
            linewidth=linewidth,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=display_tool(tool),
            zorder=3 if tool == "z-fasta-fai" else 2,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Wall Time (s, log scale)", fontsize=10)
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.935,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    ax.grid(alpha=0.28, which="both")
    ax.set_axisbelow(True)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=min(len(handles), 8),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.88)
    return _save(fig, out)


def _fig_metric_bars(
    work: pd.DataFrame,
    out: Path,
    *,
    tools: list[str],
    comparisons: pd.DataFrame,
    value_col: str,
    std_col: str | None,
    ylabel: str,
    title: str,
    value_floor: float = 1e-3,
    fig_note: str | None = None,
    after_bars=None,
    bar_width: float = 0.095,
    legend_ncol: int | None = None,
) -> Path:
    """Grouped bars with z-fasta ratio badges."""
    datasets = sorted(work["dataset"].unique(), key=dataset_sort_key)
    tools = [t for t in tools if t in work["tool"].unique()]

    fig, ax = plt.subplots(figsize=(14, 7.2))
    width = min(bar_width, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for i, tool in enumerate(tools):
        color = COLORS.get(tool, "#888888")
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            val = max(float(row[value_col].values[0]), value_floor)
            std = 0.0
            if std_col and std_col in row.columns:
                raw_std = row[std_col].values[0]
                if pd.notna(raw_std):
                    std = max(float(raw_std), 0.0)
            xpos = di + (i - len(tools) / 2 + 0.5) * width
            bar_kw: dict = {
                "width": width,
                "color": color,
                "alpha": 0.88,
                "zorder": 2,
            }
            if std > 0:
                bar_kw.update(
                    {
                        "yerr": std,
                        "capsize": 2,
                        "error_kw": {
                            "elinewidth": 0.8,
                            "capthick": 0.8,
                            "ecolor": "#333333",
                        },
                    }
                )
            ax.bar(xpos, val, **bar_kw)
            bar_tops[(ds, tool)] = (xpos, val)

    if after_bars is not None:
        after_bars(ax, work, datasets, tools, width)

    if not comparisons.empty:
        _annotate_headline_comparisons(ax, datasets, tools, width, bar_tops, comparisons)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=10)
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    if fig_note:
        fig.text(
            0.5,
            0.935,
            fig_note,
            ha="center",
            va="top",
            fontsize=9,
            color="#444444",
            style="italic",
        )
    ax.grid(axis="y", alpha=0.28, which="both")
    ax.set_axisbelow(True)

    patches = [
        mpatches.Patch(color=COLORS.get(tool, "#888888"), label=display_tool(tool))
        for tool in tools
    ]
    extra_handles, extra_labels = ax.get_legend_handles_labels()
    handles = patches + list(extra_handles)
    labels = [display_tool(t) for t in tools] + list(extra_labels)
    ncol = legend_ncol if legend_ncol is not None else min(len(handles), 8)
    fig.legend(
        handles=handles,
        labels=labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=ncol,
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.88)
    return _save(fig, out)


def _fig_headline_metric_bars(
    df: pd.DataFrame,
    out: Path,
    *,
    value_col: str,
    std_col: str | None,
    ylabel: str,
    title: str,
    value_floor: float = 1e-3,
    fig_note: str | None = None,
    after_bars=None,
) -> Path:
    """Grouped bars for headline tools across real datasets."""
    work = filter_headline_perf(df)
    return _fig_metric_bars(
        work,
        out,
        tools=HEADLINE_PERF_TOOLS,
        comparisons=build_headline_ratio_comparisons(work, value_col),
        value_col=value_col,
        std_col=std_col,
        ylabel=ylabel,
        title=title,
        value_floor=value_floor,
        fig_note=fig_note,
        after_bars=after_bars,
    )


def _draw_input_file_reference(
    ax,
    work: pd.DataFrame,
    datasets: list[str],
    tools: list[str],
    width: float,
) -> None:
    """Dashed reference at FASTA file size (MiB) per dataset cluster."""
    for di, ds in enumerate(datasets):
        ref = work[(work["dataset"] == ds) & (work["tool"] == "z-fasta-fai")]
        if ref.empty or "input_mib" not in ref.columns:
            continue
        input_mib = ref["input_mib"].values[0]
        if pd.isna(input_mib) or float(input_mib) <= 0:
            continue
        input_mb = float(input_mib)
        cluster_xs = [di + (i - len(tools) / 2 + 0.5) * width for i in range(len(tools))]
        span_lo = min(cluster_xs) - width * 0.65
        span_hi = max(cluster_xs) + width * 0.65
        ax.hlines(
            input_mb,
            span_lo,
            span_hi,
            colors="#666666",
            linestyles=(0, (5, 4)),
            linewidth=1.0,
            alpha=0.65,
            zorder=1,
        )
    ax.plot([], [], color="#666666", linestyle="--", linewidth=1.0, label="input file size")


def fig_memory_rss(df: pd.DataFrame, out: Path) -> Path:
    """Peak RSS for headline tools across real datasets."""
    return _fig_headline_metric_bars(
        df,
        out,
        value_col="peak_rss_mb",
        std_col="peak_rss_stddev_mb",
        ylabel="Peak RSS (MiB, log scale)",
        title="Real-Dataset Indexing: Peak RSS",
        value_floor=1e-3,
        after_bars=_draw_input_file_reference,
        fig_note=(
            "Bar labels: competitor RSS / z-fasta RSS. "
            "Dashed grey = FASTA file size on disk (MiB)."
        ),
    )


def fig_page_faults(df: pd.DataFrame, out: Path) -> Path:
    """Minor page faults for headline tools across real datasets."""
    return _fig_headline_metric_bars(
        df,
        out,
        value_col="minor_faults",
        std_col=None,
        ylabel="Minor Page Faults (log scale)",
        title="Real-Dataset Indexing: Minor Page Faults",
        value_floor=1.0,
        fig_note="Bar labels: competitor faults / z-fasta faults.",
    )


def fig_product_comparison_tradeoff(
    df: pd.DataFrame,
    out: Path,
) -> Path | None:
    """Time vs peak RSS scatter, one facet per real dataset."""
    work = filter_tools(df, PRODUCT_COMPARISON_TOOLS)
    if work.empty:
        return None

    datasets = sorted(work["dataset"].unique(), key=dataset_sort_key)
    tools = [t for t in PRODUCT_COMPARISON_TOOLS if t in work["tool"].unique()]

    fig, axes = plt.subplots(1, len(datasets), figsize=(5 * len(datasets), 5.2))
    if len(datasets) == 1:
        axes = [axes]

    legend_handles: list = []
    legend_labels: list[str] = []

    def metric_ticks(values: list[float]) -> list[float]:
        low, high = min(values), max(values)
        if low == high:
            return [low]
        if high / low <= 1.5:
            return [low, (low * high) ** 0.5, high]
        return [low, high]

    for ax, ds in zip(axes, datasets):
        sub = work[work["dataset"] == ds]
        for tool in tools:
            row = sub[sub["tool"] == tool]
            if row.empty:
                continue
            color = COLORS.get(tool, "#888")
            label = display_tool(tool)
            ax.scatter(
                max(float(row["peak_rss_mb"].values[0]), 1e-3),
                max(float(row["mean"].values[0]), 1e-6),
                s=110,
                color=color,
                edgecolors="white",
                linewidths=0.6,
                zorder=3,
            )
            if label not in legend_labels:
                legend_handles.append(
                    mpatches.Patch(color=color, label=label)
                )
                legend_labels.append(label)

        ax.set_xscale("log")
        ax.set_yscale("log")
        rss_values = [max(float(value), 1e-3) for value in sub["peak_rss_mb"]]
        wall_values = [max(float(value), 1e-6) for value in sub["mean"]]
        ax.xaxis.set_major_locator(mticker.FixedLocator(metric_ticks(rss_values)))
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_locator(mticker.FixedLocator(metric_ticks(wall_values)))
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.4f"))
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_title(ds, fontsize=10, fontweight="bold")
        ax.grid(alpha=0.3, which="both")

    axes[0].set_ylabel("Wall Time (s, log scale)", fontsize=10)
    for ax in axes:
        ax.set_xlabel("Peak RSS (MiB, log scale)", fontsize=10)

    fig.suptitle(
        "Format and Dedup Comparison: Time vs Memory",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )

    fig.legend(
        handles=legend_handles,
        labels=legend_labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.12),
        ncol=len(legend_handles),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(top=0.88, bottom=0.22, wspace=0.28)
    return _save(fig, out)


def fig_edge_heatmap(df: pd.DataFrame, out: Path) -> Path:
    """Raw FAI/ZFI exit codes plus the composite contract result."""
    cases = df["test_case"].tolist()
    tool_cols = edge_exit_columns(df)
    tool_labels = [edge_exit_label(column) for column in tool_cols]

    n_tools = len(tool_cols)
    match_x = n_tools + 0.3

    fig, ax = plt.subplots(
        figsize=(max(8, n_tools * 1.15 + 2), max(7, len(cases) * 0.34))
    )

    for i, case in enumerate(cases):
        yi = len(cases) - i - 1
        for j, col in enumerate(tool_cols):
            val = df.iloc[i][col]
            missing = pd.isna(val) or int(val) == 127
            clr = "#9E9E9E" if missing else "#4CAF50" if int(val) == 0 else "#F44336"
            ax.add_patch(
                plt.Rectangle((j, yi), 1, 1, facecolor=clr, edgecolor="white", lw=2)
            )
            ax.text(
                j + 0.5,
                yi + 0.5,
                "n/a" if missing else str(int(val)),
                ha="center",
                va="center",
                fontsize=7,
                color="white",
                fontweight="bold",
            )

        m = df.iloc[i]["output_match"]
        if m == "MATCH":
            clr, label = "#4CAF50", "Y"
        elif m == "REVIEW":
            clr, label = "#7E57C2", "R"
        else:
            clr, label = "#FF9800", "N"
        ax.add_patch(
            plt.Rectangle(
                (match_x, yi), 1, 1, facecolor=clr, edgecolor="white", lw=2
            )
        )
        ax.text(
            match_x + 0.5,
            yi + 0.5,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            fontweight="bold",
        )

    ax.set_xlim(0, match_x + 1.2)
    ax.set_ylim(0, len(cases))
    ax.set_xticks([j + 0.5 for j in range(n_tools)] + [match_x + 0.5])
    ax.set_xticklabels(
        tool_labels + ["Contract"],
        fontweight="bold",
        rotation=90,
        ha="center",
        va="top",
    )
    ax.set_yticks([len(cases) - i - 0.5 for i in range(len(cases))])
    ax.set_yticklabels(cases, fontsize=8)
    ax.set_title(
        "Index Edge Cases: raw exit codes and contract result",
        fontweight="bold",
        pad=12,
    )

    patches = [
        mpatches.Patch(color="#4CAF50", label="Exit 0 (accepted)"),
        mpatches.Patch(color="#F44336", label="Non-zero exit"),
        mpatches.Patch(color="#9E9E9E", label="Not run"),
        mpatches.Patch(color="#7E57C2", label="Review"),
        mpatches.Patch(color="#FF9800", label="Contract mismatch"),
    ]
    ax.legend(handles=patches, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
    fig.subplots_adjust(bottom=0.28, right=0.82)
    return _save(fig, out)


# Markdown report.


def normalize_markdown(chunks: list[str]) -> str:
    """Join chunks with one blank line between blocks; avoid MD012 (multiple blanks)."""
    cleaned: list[str] = []
    for chunk in chunks:
        text = chunk.strip()
        if not text:
            continue
        while "\n\n\n" in text:
            text = text.replace("\n\n\n", "\n\n")
        cleaned.append(text)
    return "\n\n".join(cleaned) + "\n"


def short_tool_version(value: str | None, prefix: str) -> str:
    """Normalize captured tool output to a compact version label."""
    cleaned = (value or "unknown").split(" (", 1)[0]
    if cleaned.startswith(prefix):
        cleaned = cleaned[len(prefix) :]
    if cleaned != "unknown" and not cleaned.startswith("v"):
        cleaned = f"v{cleaned}"
    return cleaned


def readable_list(values: list[str] | tuple[str, ...]) -> str:
    """Join a short list for prose."""
    if len(values) < 2:
        return "".join(values)
    if len(values) == 2:
        return " and ".join(values)
    return f"{', '.join(values[:-1])}, and {values[-1]}"


def md_datasets() -> str:
    """Describe the three real datasets without exposing local storage paths."""
    lines = []
    for spec in REAL_DATASETS:
        lines.append(f"- **{spec['id']}:** {spec['description']}.")
    return "\n".join(lines)


def md_overview(manifest: dict | None) -> str:
    """Explain the benchmark's scope and the meaning of each z-fasta lane."""
    product_dataset_text = readable_list(list(product_workloads(manifest)))
    return "\n\n".join(
        [
            "This report presents the benchmark sections captured by the selected run: "
            "correctness, real-data performance and memory, format and dedup behavior, and "
            "synthetic scaling where available.",
            (
                "Cross-tool tables use `z-fasta index --emit-fai` because the peer tools "
                "write samtools-compatible `.fai` indexes. The z-fasta-only comparison "
                f"separately covers `.fai` and the default `.zfi`, with dedup on and "
                f"off, on {product_dataset_text}."
            ),
            (
                "The **rust-bio** label needs one clarification: rust-bio has no equivalent "
                "standalone FAI index command, so this benchmark lane uses the strict custom "
                "Rust FAI writer in the rust-bio wrapper. The label stays short in tables and "
                "figures, but the measured indexing code is custom."
            ),
            (
                "Correctness covers strict FAI compatibility and z-fasta's `.zfi` support "
                "for ragged wrapping, trailing spaces, blank lines, and mixed line endings."
            ),
        ]
    )


def md_run_provenance(manifest: dict | None) -> str:
    """Summarize the build and benchmark setup."""
    if not manifest:
        return "_No run information available._"

    build = manifest.get("build") or {}
    target = str(build.get("target", "unknown")).split(".", 1)[0]
    mode = manifest.get("mode", "warm")
    return "\n".join(
        (
            f"- **Build:** z-fasta "
            f"{short_tool_version(manifest.get('z_fasta'), 'z-fasta ')}; "
            f"Zig {build.get('zig_version', 'unknown')}; "
            f"{build.get('optimize', 'unknown')}; `{target}`.",
            f"- **Benchmark:** zebrac "
            f"{short_tool_version(manifest.get('zebrac'), 'zebrac ')}; "
            f"{mode} cache; {manifest.get('runs', '?')} measured samples after "
            f"{manifest.get('warmup', '?')} warmups; "
            f"{manifest.get('duration_ms', '?')} ms budget per command.",
        )
    )



def md_perf_headline_table(df: pd.DataFrame) -> str:
    """Wall-time table for headline tools (z-fasta default + competitors)."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "mean",
        lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
        empty_msg="_No performance data._",
    )


def md_throughput_headline_table(df: pd.DataFrame) -> str:
    """Throughput table for headline tools."""
    if "throughput_mibs" not in df.columns:
        return "_No throughput data._"
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "throughput_mibs",
        lambda row: f"{row['throughput_mibs']:.1f} MiB/s",
        empty_msg="_No throughput data._",
        require_non_null=True,
    )


def md_zfasta_vs_tools_table(df: pd.DataFrame) -> str:
    """z-fasta vs competitor comparisons (matches Figure 1 bar labels)."""
    work = filter_headline_perf(df)
    comparisons = build_headline_comparisons(work)
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Speedup",
        fmt_zf=lambda r: f"{r.zfasta_s:.4f}s",
        fmt_comp=lambda r: f"{r.comp_s:.4f}s",
        extra_columns=lambda r: {
            "z-fasta MiB/s": f"{r.zfasta_mibs:.1f}" if r.zfasta_mibs else "n/a",
            "Competitor MiB/s": f"{r.comp_mibs:.1f}" if r.comp_mibs else "n/a",
            "Throughput ratio": _format_speedup(r.throughput_x),
        },
    )


def md_zfasta_vs_ratio_table(
    comparisons: pd.DataFrame,
    *,
    group_label: str = "Dataset",
    vs_label: str = "z-fasta vs",
    zf_label: str,
    comp_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    extra_columns=None,
) -> str:
    """Generic baseline vs peer table (matches bar-label ratios on figures)."""
    if comparisons.empty:
        return "_No comparison data._"
    rows: list[dict] = []
    for row in comparisons.itertuples(index=False):
        entry = {
            group_label: row.dataset,
            vs_label: row.competitor,
            zf_label: fmt_zf(row),
            comp_label: fmt_comp(row),
            ratio_label: _format_speedup(row.ratio),
        }
        if extra_columns:
            entry.update(extra_columns(row))
        rows.append(entry)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_zfasta_vs_rss_table(df: pd.DataFrame) -> str:
    work = filter_headline_perf(df)
    comparisons = build_headline_ratio_comparisons(work, "peak_rss_mb")
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="RSS ratio",
        fmt_zf=lambda r: f"{r.zfasta_v:.2f} MiB",
        fmt_comp=lambda r: f"{r.comp_v:.2f} MiB",
    )


def md_zfasta_vs_faults_table(df: pd.DataFrame) -> str:
    work = filter_headline_perf(df)
    comparisons = build_headline_ratio_comparisons(work, "minor_faults")
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Faults ratio",
        fmt_zf=lambda r: f"{int(r.zfasta_v):,}",
        fmt_comp=lambda r: f"{int(r.comp_v):,}",
    )


def md_real_performance_result(df: pd.DataFrame) -> str:
    comparisons = build_headline_comparisons(filter_headline_perf(df))
    if comparisons.empty:
        return ""
    losses = comparisons[comparisons["ratio"] < 1.0]
    if losses.empty:
        slowest_win = comparisons.loc[comparisons["ratio"].idxmin()]
        datasets = sorted(df["dataset"].unique(), key=dataset_sort_key)
        return (
            f"**Measured result:** z-fasta is faster than every reported peer on "
            f"{readable_list(datasets)}. The narrowest lead is {slowest_win.ratio:.2f}x against "
            f"{slowest_win.competitor} on {slowest_win.dataset}."
        )
    loss = losses.loc[losses["ratio"].idxmin()]
    return (
        f"**Measured result:** {loss.competitor} is faster than z-fasta on {loss.dataset} "
        f"(z-fasta speedup {loss.ratio:.2f}x)."
    )


def md_performance_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Performance section body: captioned tables, figure, and reading notes."""
    throughput = md_throughput_headline_table(df)
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_cmp = nums.next_table()
    f_perf = nums.next_figure()
    blocks = [
        "Headline comparison on the three human reference datasets (see **Datasets**). "
        "**z-fasta (.fai)** is `index --emit-fai` with dedup. It constructs the "
        "peer-compatible output while discarding the emitted stream; peers write sidecar "
        "files. Output destinations and duplicate-name policies therefore are not identical. "
        "The product table isolates z-fasta's dedup cost. CLI default `index` writes `.zfi`; "
        "see **z-fasta Format and Dedup Comparison**.",
        md_real_performance_result(df),
        f"**Table {t_wall}:** Zebrac wall time per dataset (seconds, mean ± stddev, warm cache). "
        "Lower is better. Tool order: z-fasta, noodles, rust-bio, samtools, seqkit, "
        "fastahack, pyfaidx.",
        md_perf_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_tp}:</strong> Input throughput (MiB/s), calculated from "
        f"input size and wall time. Higher is better. Same tools and order "
        f"as Table {t_wall}.</summary>",
        "",
        throughput,
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "Speedup = competitor wall time / z-fasta wall time (higher = z-fasta faster). "
        f"Same ratios appear as bar labels on Figure {f_perf}. Throughput ratio uses MiB/s.</summary>",
        "",
        md_zfasta_vs_tools_table(df),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_perf}:** Dual-axis chart for Table {t_wall} and Table {t_tp}. Bars show wall time "
        "(left axis, log scale). Lines and markers show throughput (right axis). "
        "Legend colors: gold = z-fasta; bronze = noodles; red-brown = rust-bio; grey = "
        "samtools; teal = seqkit; pink = fastahack; blue = pyfaidx. Bar labels show "
        f"z-fasta speedup vs each tool (see Table {t_cmp}); gold dashed line = z-fasta wall time "
        "per dataset.",
        f"![Figure {f_perf}: real-dataset wall time and throughput](results/figures/performance.png)",
        (
            f"**Reading Figure {f_perf}**\n"
            "- **Bars (left):** zebrac mean wall time. Error bars are one standard deviation.\n"
            "- **Lines and markers (right):** input MiB/s for the same tool at each dataset.\n"
            "- **Legend order:** z-fasta, noodles, rust-bio, samtools, seqkit, "
            "fastahack, pyfaidx.\n"
            "- **Bar labels:** `1x` on z-fasta (baseline); competitor labels = z-fasta speedup "
            f"(competitor time / z-fasta time). Border color matches the tool. Dashed gold line "
            f"= z-fasta wall time for that dataset. Details in Table {t_cmp}.\n"
            "- **z-fasta only:** this section uses the peer-comparable `.fai` lane; format and "
            "dedup tradeoffs are in **z-fasta Format and Dedup Comparison**."
        ),
    ]
    return "\n\n".join(blocks)


def md_scaling_headline_table(
    df: pd.DataFrame, param_col: str, param_label: str
) -> str:
    """Wall-time pivot for headline tools on a scaling sweep."""
    work = filter_headline_perf(df)
    if work.empty:
        return "_No scaling data._"
    pivot = work.pivot_table(
        index=param_col, columns="tool", values="mean", aggfunc="first"
    )
    pivot = _ordered_pivot(pivot, HEADLINE_PERF_TOOLS)
    pivot.index = [_scaling_param_label(param_col, v) for v in pivot.index]
    pivot.index.name = param_label
    return pivot.to_markdown(floatfmt=".4f")


def md_scaling_zfasta_rss_table(
    df: pd.DataFrame, param_col: str, param_label: str
) -> str:
    """Peak RSS for the z-fasta FAI lane across one scaling sweep."""
    work = df[df["tool"] == "z-fasta-fai"].sort_values(param_col)
    if work.empty:
        return "_No z-fasta scaling memory data._"
    rows = []
    for row in work.itertuples(index=False):
        rows.append(
            {
                param_label: _scaling_param_label(param_col, getattr(row, param_col)),
                "Peak RSS (MiB)": f"{row.peak_rss_mb:.2f}",
                "RSS spread (MiB)": (
                    f"{row.peak_rss_max_mb - row.peak_rss_min_mb:.2f}"
                ),
            }
        )
    return pd.DataFrame(rows).to_markdown(index=False)


def md_zfasta_vs_scaling_table(
    df: pd.DataFrame, param_col: str, param_label: str
) -> str:
    work = filter_headline_perf(df)
    comparisons = build_scaling_ratio_comparisons(work, param_col, "mean")
    return md_zfasta_vs_ratio_table(
        comparisons,
        group_label=param_label,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Time ratio",
        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
    )


def md_scaling_result(comparisons: pd.DataFrame) -> str:
    """Summarize scaling losses from the display labels in the ratio table."""
    losses = comparisons[comparisons["ratio"] < 1.0]
    if losses.empty:
        return "**Measured result:** z-fasta is fastest at every measured point."
    grouped: dict[str, list[str]] = {}
    for row in losses.sort_values("ratio").head(3).itertuples(index=False):
        grouped.setdefault(row.dataset, []).append(row.competitor)
    names = "; ".join(
        f"{readable_list(peers)} at {point}"
        for point, peers in grouped.items()
    )
    suffix = "; additional losses remain in the table" if len(losses) > 3 else ""
    return f"**Measured result:** faster than z-fasta: {names}{suffix}."


def md_scaling_section(
    df: pd.DataFrame,
    nums: ReportCounters,
    *,
    intro: str,
    subsection: str,
    param_col: str,
    param_label: str,
    figure_path: str,
    figure_title: str,
    table_x_label: str,
    reading_x_label: str,
) -> str:
    """Shared template for file-size and sequence-count scaling sections."""
    t_wall = nums.next_table()
    t_cmp = nums.next_table()
    t_rss = nums.next_table()
    f_scale = nums.next_figure()
    comparisons = build_scaling_ratio_comparisons(
        filter_headline_perf(df), param_col, "mean"
    )
    result = md_scaling_result(comparisons)
    blocks = [
        intro,
        result,
        subsection,
        "<details>",
        f"<summary><strong>Table {t_wall}:</strong> Wall time per {table_x_label} "
        "(zebrac mean, seconds). Tool order: z-fasta, noodles, rust-bio, samtools, "
        "seqkit, fastahack, pyfaidx.</summary>",
        "",
        md_scaling_headline_table(df, param_col, param_label),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor at each "
        f"{table_x_label}. Time ratio = competitor wall time / z-fasta wall time.</summary>",
        "",
        md_zfasta_vs_scaling_table(df, param_col, param_label),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_rss}:</strong> z-fasta peak RSS across the same "
        f"{table_x_label} sweep.</summary>",
        "",
        md_scaling_zfasta_rss_table(df, param_col, param_label),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_scale}:** Table {t_wall} as lines (log scales on both axes). "
        "One line per tool; absolute wall times only.",
        f"![Figure {f_scale}: {figure_title}]({figure_path})",
        (
            f"**Reading Figure {f_scale}**\n"
            f"- Lines: mean wall time vs {reading_x_label}. Both axes use log scale.\n"
            f"- Lower on the chart = faster at that {table_x_label}.\n"
            "- Legend order: z-fasta, noodles, rust-bio, samtools, seqkit, fastahack, "
            "pyfaidx.\n"
            f"- Time ratios are in Table {t_cmp}; the chart does not show them."
        ),
    ]
    return "\n\n".join(blocks)


def md_product_effect_table(df: pd.DataFrame) -> str:
    """Dedup-on vs no-dedup wall time and RSS within each output format."""
    formats = (
        (".fai", "z-fasta-fai", "z-fasta-fai-nodedup"),
        (".zfi", "z-fasta-zfi", "z-fasta-zfi-nodedup"),
    )
    rows = []
    for dataset in sorted(df["dataset"].unique(), key=dataset_sort_key):
        for output, dedup_tool, no_dedup_tool in formats:
            dedup = df[(df["dataset"] == dataset) & (df["tool"] == dedup_tool)]
            no_dedup = df[(df["dataset"] == dataset) & (df["tool"] == no_dedup_tool)]
            if dedup.empty or no_dedup.empty:
                continue
            dedup_row = dedup.iloc[0]
            no_dedup_row = no_dedup.iloc[0]
            no_dedup_wall = float(no_dedup_row["mean"])
            no_dedup_rss = float(no_dedup_row["peak_rss_mb"])
            wall_ratio = (
                float(dedup_row["mean"]) / no_dedup_wall
                if no_dedup_wall > 0
                else None
            )
            rss_ratio = (
                float(dedup_row["peak_rss_mb"]) / no_dedup_rss
                if no_dedup_rss > 0
                else None
            )
            rows.append(
                {
                    "Dataset": dataset,
                    "Output": output,
                    "Dedup wall": f"{dedup_row['mean']:.4f}s ±{dedup_row['stddev']:.4f}",
                    "No-dedup wall": (
                        f"{no_dedup_row['mean']:.4f}s ±{no_dedup_row['stddev']:.4f}"
                    ),
                    "Wall ratio": f"{wall_ratio:.2f}x" if wall_ratio else "n/a",
                    "Dedup RSS": f"{dedup_row['peak_rss_mb']:.2f} MiB",
                    "No-dedup RSS": f"{no_dedup_row['peak_rss_mb']:.2f} MiB",
                    "RSS ratio": f"{rss_ratio:.2f}x" if rss_ratio else "n/a",
                }
            )
    if not rows:
        return "_No complete format and dedup pairs._"
    return pd.DataFrame(rows).to_markdown(index=False)


def md_product_comparison_section(
    df: pd.DataFrame,
    nums: ReportCounters,
) -> str:
    """Thin z-fasta-only view of output format and dedup cost."""
    t_effect = nums.next_table()
    f_trade = nums.next_figure()
    datasets = sorted(df["dataset"].unique(), key=dataset_sort_key)
    dataset_text = readable_list([f"**{dataset}**" for dataset in datasets])

    effects = []
    for dataset in datasets:
        for output, dedup_tool, no_dedup_tool in (
            (".fai", "z-fasta-fai", "z-fasta-fai-nodedup"),
            (".zfi", "z-fasta-zfi", "z-fasta-zfi-nodedup"),
        ):
            dedup = df[(df["dataset"] == dataset) & (df["tool"] == dedup_tool)]
            no_dedup = df[(df["dataset"] == dataset) & (df["tool"] == no_dedup_tool)]
            if dedup.empty or no_dedup.empty:
                continue
            effects.append(
                (
                    float(dedup.iloc[0]["peak_rss_mb"])
                    / float(no_dedup.iloc[0]["peak_rss_mb"]),
                    dataset,
                    output,
                )
            )
    if effects:
        rss_ratio, rss_dataset, rss_output = max(effects)
        result = (
            f"**Measured result:** the largest dedup memory cost is {rss_ratio:.2f}x on "
            f"{rss_dataset} (`{rss_output}`)."
        )
    else:
        result = ""

    blocks = [
        (
            f"This focused z-fasta comparison uses {dataset_text}. It isolates output format "
            "(`.fai` or `.zfi`) and duplicate-name tracking (dedup or `--no-dedup`). Peer "
            "indexers stay in the cross-tool performance sections, where they answer a "
            "different question."
        ),
        result,
        (
            f"**Table {t_effect}:** Measured wall time and peak RSS for each matched pair. "
            "Each ratio is dedup divided by no-dedup within the same dataset and output "
            "format. Values above `1x` are the measured dedup cost."
        ),
        md_product_effect_table(df),
        '<div style="margin: 1.5em 0"></div>',
        (
            f"**Figure {f_trade}:** The same measured configurations plotted by peak RSS "
            "and wall time. Both axes use log scale; lower-left is faster with less process "
            "memory."
        ),
        (
            f"![Figure {f_trade}: format and dedup time versus memory]"
            "(results/figures/product_comparison.png)"
        ),
        (
            f"**Reading Figure {f_trade}**\n"
            f"- Panels: {readable_list(datasets)}.\n"
            "- Gold points: `.fai`; orange points: `.zfi`.\n"
            "- Lighter points: `--no-dedup`; darker points: dedup.\n"
            f"- Exact measurements and within-format ratios are in Table {t_effect}."
        ),
    ]
    return "\n\n".join(blocks)


def md_memory_rss_headline_table(df: pd.DataFrame) -> str:
    """Peak RSS table for headline tools (same order as Figure 1)."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "peak_rss_mb",
        lambda row: f"{row['peak_rss_mb']:.2f} MiB",
        empty_msg="_No memory data._",
    )


def md_page_faults_headline_table(df: pd.DataFrame) -> str:
    """Minor page-fault table for headline tools."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "minor_faults",
        lambda row: f"{int(row['minor_faults']):,}",
        empty_msg="_No page-fault data._",
    )


def md_major_faults_headline_table(df: pd.DataFrame) -> str:
    """Major page-fault table for headline tools."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "major_faults",
        lambda row: f"{int(row['major_faults']):,}",
        empty_msg="_No page-fault data._",
    )


def md_memory_usage_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Memory usage section: headline peak RSS on real datasets."""
    t_rss = nums.next_table()
    t_cmp = nums.next_table()
    f_mem = nums.next_figure()
    comparisons = build_headline_ratio_comparisons(
        filter_headline_perf(df), "peak_rss_mb"
    )
    lower = comparisons[comparisons["ratio"] < 1.0]
    if lower.empty:
        result = "**Measured result:** no reported peer uses less peak RSS than z-fasta."
    else:
        best = lower.loc[lower["ratio"].idxmin()]
        result = (
            f"**Measured result:** {best.competitor} uses less peak RSS on {best.dataset} "
            f"({best.comp_v:.2f} MiB vs {best.zfasta_v:.2f} MiB). The headline z-fasta "
            "lane includes default duplicate-name tracking."
        )
    blocks = [
        "Same zebrac runs as **Performance: Real Datasets**. z-fasta (.fai) only; "
        "other product lanes are in **z-fasta Format and Dedup Comparison**.",
        result,
        (
            "### Peak RSS\n\n"
            "zebrac starts a new process for each sample and records peak RSS when it "
            "exits. That is the most RAM the process had in use at once (Linux "
            "`ru_maxrss` via `getrusage`). Table and figure show the mean across "
            "samples.\n\n"
            "Use this to compare tools on the same file and host. Grey dashed lines show file "
            "size as context; output and dedup state can still affect memory. A low bar does "
            "not always mean less work, so wall time and page faults provide the surrounding "
            "context."
        ),
        f"**Table {t_rss}:** Peak RSS (MiB, zebrac mean). Same tool order as Performance.",
        md_memory_rss_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "RSS ratio = competitor peak RSS / z-fasta peak RSS. Same ratios as bar labels on "
        f"Figure {f_mem}.</summary>",
        "",
        md_zfasta_vs_rss_table(df),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_mem}:** Table {t_rss} as grouped bars (log scale). Bar labels = "
        f"RSS ratio (see Table {t_cmp}). Gold dashed line = z-fasta RSS; grey dashed = "
        "FASTA file size on disk.",
        f"![Figure {f_mem}: real-dataset peak RSS](results/figures/memory.png)",
        (
            f"**Reading Figure {f_mem}**\n"
            "- Bars: mean peak RSS. `1x` on z-fasta; other labels = competitor / z-fasta.\n"
            "- Gold dashed line: z-fasta RSS for that dataset.\n"
            "- Grey dashed line: input FASTA size.\n"
            f"- Details in Table {t_cmp}."
        ),
    ]
    return "\n\n".join(blocks)


def md_page_faults_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Page-fault section: headline minor faults on real datasets."""
    major = md_major_faults_headline_table(df)
    t_minor = nums.next_table()
    t_cmp = nums.next_table()
    t_major = nums.next_table()
    f_pf = nums.next_figure()
    blocks = [
        "Same zebrac samples as the real-dataset performance and memory sections.",
        (
            "### Page faults\n\n"
            "A **minor** fault maps a page without reading disk. A **major** fault reads "
            "from disk. zebrac reports both per run, like wall time. "
            f"Table {t_major} lists major faults separately."
        ),
        f"**Table {t_minor}:** Minor page faults (zebrac mean). Same tool order as Performance.",
        md_page_faults_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "Faults ratio = competitor minor faults / z-fasta minor faults. Same ratios as bar "
        f"labels on Figure {f_pf}.</summary>",
        "",
        md_zfasta_vs_faults_table(df),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_major}:</strong> Major page faults.</summary>",
        "",
        major,
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_pf}:** Table {t_minor} as grouped bars (log scale). Bar labels = "
        f"Faults ratio (see Table {t_cmp}). Gold dashed line = z-fasta minor faults.",
        f"![Figure {f_pf}: real-dataset minor page faults](results/figures/page_faults.png)",
        (
            f"**Reading Figure {f_pf}**\n"
            "- Bars: mean minor faults. `1x` on z-fasta; other labels = competitor / z-fasta.\n"
            "- Gold dashed line: z-fasta count for that dataset.\n"
            "- Lower bar = fewer minor faults.\n"
            f"- Details in Table {t_cmp}."
        ),
    ]
    return "\n\n".join(blocks)


def load_messy_index(
    results_dir: Path, manifest: dict | None = None
) -> pd.DataFrame | None:
    """Load messy FASTA indexing results pinned to the run timestamp."""
    if not manifest or not manifest.get("timestamp"):
        return None
    configured = (manifest.get("sections") or {}).get("messy")
    if not configured:
        configured = (manifest.get("messy") or {}).get("directory")
    d = results_dir / (configured or f"messy_{manifest['timestamp']}")
    if not d.is_dir():
        return None
    meta_by_command: dict[str, str] = {}
    metadata_name = (manifest.get("messy") or {}).get("metadata") or "metadata.jsonl"
    meta_path = d / metadata_name
    if meta_path.exists():
        for line in meta_path.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            meta_by_command[row["command"]] = row["tool"]

    rows = []
    for f in sorted(d.glob("*.json")):
        data = json.loads(f.read_text())
        variant = f.stem
        for r in data.get("results", []):
            cmd = r.get("command", "")
            tool = meta_by_command.get(cmd)
            if not tool:
                tool = cmd.split()[0].split("/")[-1]
                if tool == "faidx":
                    tool = "pyfaidx"
            if tool == "rustbio-custom-index":
                tool = "rust-bio"
            if tool == "z-fasta":
                tool = (
                    "z-fasta (.fai)"
                    if "--emit-fai" in cmd
                    else "z-fasta (.zfi)"
                )
            failed = r.get("failed_sample_count", 0)
            samples = r.get("sample_count", 0)
            success = failed == 0 and samples > 0
            rows.append({"variant": variant, "tool": tool, "success": success})
    if not rows:
        return None
    return pd.DataFrame(rows)


def md_messy_table(df: pd.DataFrame) -> str:
    """Markdown compatibility matrix for messy FASTA variants."""
    messy_order = [
        "z-fasta (.zfi)",
        "z-fasta (.fai)",
        "noodles",
        "rust-bio",
        "samtools",
        "seqkit",
        "fastahack",
        "pyfaidx",
    ]
    tools = [t for t in messy_order if t in df["tool"].values]
    # append any unexpected tools at end
    for t in df["tool"].unique():
        if t not in tools and t != "seqtk":
            tools.append(t)
    variants = sorted(df["variant"].unique())
    rows = []
    for variant in variants:
        row: dict = {"Variant": variant}
        for tool in tools:
            sub = df[(df["variant"] == variant) & (df["tool"] == tool)]
            if sub.empty:
                row[tool] = "n/a"
            elif tool.startswith("z-fasta (.zfi)"):
                row[tool] = "ok" if sub.iloc[0]["success"] else "fail"
            else:
                row[tool] = "exit 0" if sub.iloc[0]["success"] else "nonzero"
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_edge_case_summary_table(df: pd.DataFrame) -> str:
    """Summarize comparability by the rule used for each case class."""
    classes = edge_contract_classes(df)
    meanings = {
        "fai_parity": "z-fasta (.fai) vs samtools acceptance and FAI bytes",
        "zfi_messy": "z-fasta (.zfi) side table; FAI peers expected incompatible",
        "invalid_input": "raw behavior shown; no parity claim for binary input",
    }
    rows = []
    for class_name, label in EDGE_CONTRACT_CLASSES.items():
        subset = df[classes == class_name]
        if subset.empty:
            continue
        matches = int((subset["output_match"] == "MATCH").sum())
        reviews = int((subset["output_match"] == "REVIEW").sum())
        result = (
            f"{reviews}/{len(subset)} review"
            if class_name == "invalid_input"
            else f"{matches}/{len(subset)} pass"
        )
        rows.append(
            {
                "Contract basis": label,
                "Cases": len(subset),
                "Result": result,
                "Comparison": meanings[class_name],
            }
        )
    return pd.DataFrame(rows).to_markdown(index=False)


def md_edge_case_section(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Edge-case correctness: summary table, heatmap figure, and reading notes."""
    t_summary = nums.next_table()
    f_edge = nums.next_figure()
    blocks = [
        (
            "Each heatmap row is one structural FASTA edge case. The z-fasta columns "
            "separate the "
            "samtools-compatible FAI attempt from the production `.zfi` index. The "
            "**Contract** column marks whether z-fasta met that test's rules."
        ),
        (
            "**Contract meaning:** for structural cases, `z-fasta (.fai)` must agree with "
            "samtools on FAI acceptance and bytes when both accept. If one accepts and the "
            "other rejects, the contract is `N`. For the four named messy fixtures, the "
            "`.zfi` side-table path is the contract because samtools cannot represent those "
            "layouts; the FAI rejection pair plus valid `.zfi` side table is `Y`."
        ),
        "**Scoring rules:**",
        (
            "- *Samtools parity* cases: when both z-fasta and samtools accept a file, their "
            "`.fai` output must match. When both reject it, that also counts as a match. "
            "A production `.zfi` success does not hide an FAI exit-code mismatch."
        ),
        (
            "- *z-fasta-only* messy cases: z-fasta must index files where samtools rejects "
            "non-uniform wrapping (side-table `.zfi`)."
        ),
        (
            f"**Table {t_summary}:** Comparability summary by contract class. Raw exit codes "
            "are diagnostic; `Result` applies the class-specific comparison rule."
        ),
        md_edge_case_summary_table(df),
        '<div style="margin: 1.5em 0"></div>',
        (
            f"**Figure {f_edge}:** Raw per-test exit codes and the composite contract result. "
            "Green = exit 0; red = non-zero exit; gray = tool not run; purple = review; "
            "orange = contract mismatch."
        ),
        (
            f"![Figure {f_edge}: edge-case exit codes](results/figures/edge_cases.png)"
        ),
        (
            f"**Reading Figure {f_edge}**\n"
            "- Rows: structural FASTA edge cases.\n"
            "- Tool columns show the numeric exit code; `n/a` means the optional tool was not run.\n"
            "- **Contract:** `Y` = match, `R` = review, `N` = mismatch.\n"
            "- Cell colors match the legend on the chart.\n"
            f"- Aggregate counts are in Table {t_summary}."
        ),
    ]
    return "\n\n".join(blocks)


def md_messy_section(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Messy FASTA compatibility: matrix table and reading notes."""
    t_matrix = nums.next_table()
    if "z-fasta (.zfi)" in set(df["tool"]):
        zf_lane = "`z-fasta (.zfi)` via default `index`"
    elif "z-fasta (.fai)" in set(df["tool"]):
        zf_lane = "historical `z-fasta (.fai)` data"
    else:
        zf_lane = "no z-fasta lane"
    blocks = [
        (
            "Proteome-derived messy FASTA fixtures. Each cell is zebrac with "
            "`--allow-failures` (repeated samples, not a single exit check)."
        ),
        (
            f"**z-fasta lane:** This table must use {zf_lane}. The `.fai` lane is not the "
            "messy compatibility claim; variable widths, trailing whitespace, and mixed "
            "line endings require the production `.zfi` side table."
        ),
        (
            "**What z-fasta handles:** irregular line wrapping, trailing whitespace on "
            "sequence lines, blank lines between records, and mixed CRLF/LF. Those are "
            "common export glitches, not arbitrary byte corruption. z-fasta records true "
            "sequence boundaries in a side-table `.zfi` instead of assuming every line in a "
            "record has the same width."
        ),
        (
            "**What we still reject:** empty files, headers without sequence, malformed "
            "headers, and other cases where indexing would be meaningless (see **Edge Case "
            "Correctness** above). We do not claim to repair truncated or binary-corrupted FASTA."
        ),
        (
            "`mixed_crlf` deliberately alternates LF and CRLF sequence-line endings. It is "
            "not a uniformly CRLF-wrapped control file."
        ),
        (
            f"**Table {t_matrix}:** Measured process result per variant and tool from the "
            "messy benchmark data. The z-fasta `.zfi` lane is `ok` only when every sample "
            "exits 0. Peer cells report only whether every measured sample exited 0; they do "
            "not claim sequence-retrieval compatibility."
        ),
        md_messy_table(df),
        (
            f"**Reading Table {t_matrix}**\n"
            "- Rows: messy variant names (`all_messy`, `mixed_crlf`, `mixed_widths`, "
            "`trailing_whitespace`).\n"
            "- Columns: production `z-fasta (.zfi)` followed by FAI peers.\n"
            "- `exit 0` and `nonzero` are measured process results, not output validation."
        ),
    ]
    return "\n\n".join(blocks)


def md_tools_tested(tools: dict[str, str] | None = None, z_fasta: str | None = None) -> str:
    """List each tool once, with z-fasta modes nested under the product version."""
    tools = tools or {}

    zf_version = short_tool_version(z_fasta, "z-fasta ")
    lines = [
        f"- **z-fasta ({zf_version})**\n"
        "  - `.fai`, dedup: peer-comparable cross-tool and scaling lane.\n"
        "  - `.fai`, no dedup: isolates duplicate-name tracking cost.\n"
        "  - `.zfi`, dedup: default product mode.\n"
        "  - `.zfi`, no dedup: default format without duplicate-name tracking."
    ]
    specs = (
        ("noodles", "noodles", "noodles-fasta ", "Rust noodles-fasta wrapper"),
        ("rustbio", "rust-bio", "rust-bio ", "Rust FASTA library"),
        ("samtools", "samtools", "samtools ", "reference FAI indexer"),
        ("seqkit", "seqkit", "seqkit ", "Go FASTA toolkit"),
        ("fastahack", "fastahack", "fastahack ", "C++ FAI indexer"),
        ("pyfaidx", "pyfaidx", "", "Python FAI indexer"),
    )
    for key, label, prefix, description in specs:
        if key not in tools:
            continue
        version = short_tool_version(tools[key], prefix)
        lines.append(f"- **{label} ({version})**: {description}.")
    return "\n".join(lines)


# Main.


def main():
    parser = argparse.ArgumentParser(
        description="Generate index benchmark report and figures from raw results."
    )
    parser.add_argument(
        "results_dir",
        nargs="?",
        type=Path,
        help="Results directory, defaults to bench/index/results.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Use this run manifest instead of results/LATEST.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Explicit report path, defaults to bench/index/REPORT.md.",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    results_dir = args.results_dir if args.results_dir else script_dir / "results"
    manifest = enrich_manifest(load_run_manifest(results_dir, args.manifest))
    if args.output is not None:
        report_path = args.output.resolve()
        final_figures_dir = report_path.parent / "results" / "figures"
    else:
        report_path = script_dir / "REPORT.md"
        final_figures_dir = results_dir / "figures"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    final_figures_dir.mkdir(parents=True, exist_ok=True)
    figure_stage = tempfile.TemporaryDirectory(
        prefix=".index-report-figures-",
        dir=report_path.parent,
    )
    figures_dir = Path(figure_stage.name)

    report_lines: list[str] = []

    def section(title, level=2):
        report_lines.append(f"{'#' * level} {title}")

    report_lines.append("# z-fasta Index Benchmark Report")
    report_lines.append("_Auto-generated by `generate_report.py` from zebrac results._")

    section("Overview")
    report_lines.append(md_overview(manifest))

    if manifest:
        section("Run Provenance")
        report_lines.append(md_run_provenance(manifest))

    section("Datasets")
    report_lines.append(md_datasets())

    section("Tools Tested")
    report_lines.append(
        md_tools_tested(
            (manifest or {}).get("tools"),
            (manifest or {}).get("z_fasta"),
        )
    )

    generated_figs: list[str] = []
    has_results = False
    metadata_df = load_metadata(results_dir, manifest)
    nums = ReportCounters()

    test_df = load_tests(results_dir, manifest)
    if test_df is not None and len(test_df):
        has_results = True
        section("Edge Case Correctness")
        report_lines.append(md_edge_case_section(test_df, nums))
        p = fig_edge_heatmap(test_df, figures_dir / "edge_cases.png")
        generated_figs.append(str(p))

    messy_df = load_messy_index(results_dir, manifest)
    if messy_df is not None and len(messy_df):
        has_results = True
        section("Messy FASTA Compatibility")
        report_lines.append(md_messy_section(messy_df, nums))

    perf_df = load_perf(results_dir, metadata_df, manifest)
    if perf_df is not None and len(perf_df):
        has_results = True
        section("Performance: Real Datasets")
        report_lines.append(md_performance_real_datasets(perf_df, nums))

        p = fig_performance_throughput(perf_df, figures_dir / "performance.png")
        generated_figs.append(str(p))

        section("Memory Usage: Real Datasets")
        report_lines.append(md_memory_usage_real_datasets(perf_df, nums))
        p = fig_memory_rss(perf_df, figures_dir / "memory.png")
        generated_figs.append(str(p))
    product_df = (
        filter_product_comparison(perf_df, manifest)
        if perf_df is not None and len(perf_df)
        else None
    )
    if product_df is not None and not product_df.empty:
        section("z-fasta Format and Dedup Comparison")
        report_lines.append(md_product_comparison_section(product_df, nums))
        p = fig_product_comparison_tradeoff(
            product_df,
            figures_dir / "product_comparison.png",
        )
        if p:
            generated_figs.append(str(p))

    if perf_df is not None and len(perf_df):
        section("Page Faults: Real Datasets")
        report_lines.append(md_page_faults_real_datasets(perf_df, nums))
        p = fig_page_faults(perf_df, figures_dir / "page_faults.png")
        generated_figs.append(str(p))

    size_df = load_scaling(results_dir, "scale_size", "size_mb", metadata_df, manifest)
    if size_df is not None and len(size_df):
        has_results = True
        section("Scaling: File Size")
        report_lines.append(
            md_scaling_section(
                size_df,
                nums,
                intro=(
                    "Synthetic FASTA files from 1 MiB to 1000 MiB with 100 sequences per file. "
                    "Uses the zebrac warm-cache setup in **Run Provenance** and the "
                    "peer-comparable z-fasta `.fai` lane."
                ),
                subsection=(
                    "### Wall time vs file size\n\n"
                    "Mean wall time (seconds, zebrac) at each file size. **Lower is better.** "
                    "Record count stays fixed at 100 while total bytes on disk grow."
                ),
                param_col="size_mb",
                param_label="File Size",
                figure_path="results/figures/scaling_size.png",
                figure_title="scaling by file size",
                table_x_label="file size",
                reading_x_label="file size (MiB)",
            )
        )
        p = fig_scaling_headline_lines(
            size_df,
            figures_dir / "scaling_size.png",
            param_col="size_mb",
            xlabel="File Size (MiB, log scale)",
            title="Scaling: Wall Time vs File Size",
            fig_note="One line per tool. Lower wall time is better.",
        )
        generated_figs.append(str(p))

    seq_budget_df = load_scaling(
        results_dir, "scale_seqs_budget", "seq_count", metadata_df, manifest
    )
    if seq_budget_df is not None and len(seq_budget_df):
        has_results = True
        section("Scaling: Sequence Count (Constant File Size)")
        report_lines.append(
            md_scaling_section(
                seq_budget_df,
                nums,
                intro=(
                    "Synthetic FASTA holding about **50 MiB** total per file; sequence length "
                    "shrinks as record count rises. Uses the zebrac warm-cache setup in "
                    "**Run Provenance** and the peer-comparable z-fasta `.fai` lane."
                ),
                subsection=(
                    "### Wall time vs sequence count (constant file size)\n\n"
                    "Mean wall time (seconds, zebrac) at each record count. **Lower is better.** "
                    "Isolates per-record and header overhead at near-constant file size."
                ),
                param_col="seq_count",
                param_label="Sequences",
                figure_path="results/figures/scaling_seqs_budget.png",
                figure_title="constant-size sequence scaling",
                table_x_label="record count",
                reading_x_label="sequence count",
            )
        )
        p = fig_scaling_headline_lines(
            seq_budget_df,
            figures_dir / "scaling_seqs_budget.png",
            param_col="seq_count",
            xlabel="Sequences (log scale)",
            title="Scaling: Wall Time vs Sequence Count (Constant ~50 MiB)",
            fig_note="One line per tool. Lower wall time is better.",
        )
        generated_figs.append(str(p))

    seq_fixed_df = load_scaling(
        results_dir, "scale_seqs_fixed", "seq_count", metadata_df, manifest
    )
    if seq_fixed_df is not None and len(seq_fixed_df):
        has_results = True
        section("Scaling: Sequence Count (Fixed Length)")
        report_lines.append(
            md_scaling_section(
                seq_fixed_df,
                nums,
                intro=(
                    "Synthetic FASTA with **1024 bp** per sequence; total file size grows with "
                    "record count (100k to 1M). Uses the zebrac warm-cache setup in "
                    "**Run Provenance** and the peer-comparable z-fasta `.fai` lane."
                ),
                subsection=(
                    "### Wall time vs sequence count (fixed length)\n\n"
                    "Mean wall time (seconds, zebrac) at each record count. **Lower is better.** "
                    "Bytes and index entries grow together; the 1M point is about 1 GiB on disk."
                ),
                param_col="seq_count",
                param_label="Sequences",
                figure_path="results/figures/scaling_seqs_fixed.png",
                figure_title="fixed-length sequence scaling",
                table_x_label="record count",
                reading_x_label="sequence count",
            )
        )
        p = fig_scaling_headline_lines(
            seq_fixed_df,
            figures_dir / "scaling_seqs_fixed.png",
            param_col="seq_count",
            xlabel="Sequences (log scale)",
            title="Scaling: Wall Time vs Sequence Count (1024 bp / record)",
            fig_note="One line per tool. Lower wall time is better.",
        )
        generated_figs.append(str(p))

    # Write report.
    if not has_results:
        raise SystemExit(
            "No benchmark result data found. Run at least one benchmark or correctness "
            "section first."
        )

    published_figs = []
    for staged_path in map(Path, generated_figs):
        destination = final_figures_dir / staged_path.name
        staged_path.replace(destination)
        published_figs.append(str(destination))

    report_temp = report_path.with_name(f".{report_path.name}.tmp")
    report_temp.write_text(f"{MARKDOWNLINT_DISABLE}\n\n{normalize_markdown(report_lines)}")
    report_temp.replace(report_path)
    prune_stale_pngs(final_figures_dir, [Path(path).name for path in published_figs])
    figure_stage.cleanup()

    print(f"Report: {report_path}")
    print(f"Figures: {len(published_figs)}")
    for f in published_figs:
        print(f"  {f}")


if __name__ == "__main__":
    main()
