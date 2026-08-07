#!/usr/bin/env python3
"""Generate stats benchmark REPORT.md and figures from zebrac JSON.

Mirrors bench/get and bench/index report grammar (overview, field matrix,
provenance, wall/RSS/faults, and scaling). Color map is STATS_COLORS below.
Prose follows plan/WRITING.md (ASCII; no emojis).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
PROJECT_ROOT = SCRIPT_DIR.parent.parent
FIGURES_DIR = RESULTS_DIR / "figures"

# Facet / table order matches index REPORT (Genome | Transcriptome | Proteome).
DATASET_ORDER = ["Genome", "Transcriptome", "Proteome"]

BASELINE = "z-fasta-full"

FULL_TOOLS = [
    "z-fasta-full",
    "noodles",
    "rustbio",
    "seqkit",
    "seqtk",
]

SCALING_TOOLS = [
    "z-fasta-full",
    "z-fasta-full-fai",
    "noodles",
    "rustbio",
    "seqkit",
    "seqtk",
]

# Cleaner legend/table labels for scaling.
SCALING_DISPLAY = {
    "z-fasta-full": "z-fasta",
    "z-fasta-full-fai": "z-fasta (fai)",
    "seqkit": "seqkit",
    "noodles": "noodles",
    "rustbio": "rust-bio",
    "seqtk": "seqtk (comp)",
}

REFERENCE_TOOLS = frozenset({"seqtk"})

STATS_COLORS = {
    "z-fasta-full": "#F7A41D",
    "z-fasta-full-fai": "#FFB74D",
    "seqkit": "#1565C0",
    "noodles": "#C45C26",
    "rustbio": "#8B3A2A",
    "seqtk": "#6A1B9A",
}

STATS_DISPLAY = {
    "z-fasta-full": "z-fasta full",
    "z-fasta-full-fai": "z-fasta full (fai)",
    "seqkit": "seqkit",
    "noodles": "noodles",
    "rustbio": "rust-bio",
    "seqtk": "seqtk (comp)",
}

# In peer-only full section, shorten the baseline label.
FULL_SECTION_DISPLAY = {
    **STATS_DISPLAY,
    "z-fasta-full": "z-fasta",
}

# Gap between wall | RSS | faults metric facets (GET RC / messy style).
FACET_WSPACE = 0.28

REPORT_FIGURES = frozenset(
    {
        "perf_full.png",
        "scaling_size.png",
        "scaling_size_slope.png",
        "scaling_seqs_fixed.png",
        "scaling_seqs_fixed_slope.png",
    }
)

SIZE_MB_ORDER = [1, 5, 10, 25, 50, 100, 250, 500, 1000]
SEQ_COUNT_ORDER = [100000, 250000, 500000, 1000000]


def load_index_report():
    path = SCRIPT_DIR.parent / "index" / "generate_report.py"
    spec = importlib.util.spec_from_file_location("index_generate_report", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    mod.COLORS.update(STATS_COLORS)
    mod.DISPLAY_NAMES.update(STATS_DISPLAY)
    # Shared ratio tables use _format_speedup; keep stats REPORT WRITING.md-safe.
    mod._format_speedup = _format_ratio
    return mod


def resolve_verify_status(manifest: dict, results_dir: Path) -> tuple[bool, str | None]:
    """Return (skipped, pass_count). Prefer verify_<ts>.log when present.

    Resume with `--skip-tests` (alias `--skip-verify`) can rewrite the manifest and clear `verify_pass`
    even when `verify_<timestamp>.log` already recorded a green correctness run.
    """
    ts = str(manifest.get("timestamp") or "")
    log = results_dir / f"verify_{ts}.log" if ts else None
    if log is not None and log.is_file():
        text = log.read_text(encoding="utf-8", errors="replace")
        if "ALL PASSED" in text:
            for line in text.splitlines():
                if line.startswith("Results:"):
                    parts = line.split()
                    if len(parts) >= 2:
                        return False, parts[1]
            return False, str(manifest.get("verify_pass") or "unknown")
    skipped = bool(manifest.get("verify_skipped", False))
    raw = manifest.get("verify_pass")
    if raw is None:
        return skipped, None
    return skipped, str(raw)


def tools_in_run(df: pd.DataFrame | None, order: list[str]) -> list[str]:
    if df is None or df.empty:
        return []
    present = set(df["tool"].unique())
    return [t for t in order if t in present]


def peer_tools(tools: list[str], baseline: str = BASELINE) -> list[str]:
    """All non-baseline tools (including hatched seqtk reference) for ratio badges."""
    return [t for t in tools if t != baseline]


def dataset_sort_key(name):
    try:
        return DATASET_ORDER.index(name), str(name)
    except ValueError:
        return len(DATASET_ORDER), str(name)


def load_zebrac_json(path: Path, metadata_df: pd.DataFrame | None, ir) -> pd.DataFrame:
    data = json.loads(path.read_text())
    meta_rows = ir.meta_rows_for_json(metadata_df, path)
    meta_by_command = {row.get("command"): row for row in meta_rows}

    rows = []
    for idx, result in enumerate(data.get("results", [])):
        command = result.get("command", "")
        meta = meta_by_command.get(command, {})
        if not meta and idx < len(meta_rows):
            meta = meta_rows[idx]
        wall = result.get("wall_time", {})
        peak = result.get("peak_rss", {})
        minor = result.get("minor_faults", {})
        mean_s = float(wall.get("mean", 0) or 0) / 1_000_000_000.0
        std_s = float(wall.get("std_dev", 0) or 0) / 1_000_000_000.0
        rss = float(peak.get("mean", 0) or 0) / (1024.0 * 1024.0)
        rss_std = float(peak.get("std_dev", 0) or 0) / (1024.0 * 1024.0)
        faults = float(minor.get("mean", 0) or 0)
        faults_std = float(minor.get("std_dev", 0) or 0)
        workload = meta.get("workload", path.stem.split("__", 1)[0])
        section = meta.get("section", "")
        tool = meta.get("tool")
        if not tool:
            # Fallback from filename: Genome__z-fasta-full.json
            stem = path.stem
            if "__" in stem:
                tool = stem.split("__", 1)[1]
            else:
                tool = "unknown"
        input_bytes = meta.get("input_bytes")
        input_mib = None
        throughput_mibs = None
        if input_bytes is not None:
            input_mib = float(input_bytes) / (1024.0 * 1024.0)
            if mean_s > 0:
                throughput_mibs = input_mib / mean_s
        # REAL sections: workload is the dataset name.
        dataset = workload if workload in DATASET_ORDER else None
        rows.append(
            {
                "tool": tool,
                "section": section,
                "workload": workload,
                "dataset": dataset,
                "mean": mean_s,
                "stddev": std_s,
                "peak_rss_mb": rss,
                "peak_rss_stddev_mb": rss_std,
                "minor_faults": faults,
                "minor_faults_stddev": faults_std,
                "input_bytes": input_bytes,
                "input_mib": input_mib,
                "throughput_mibs": throughput_mibs,
                "command": command,
                "raw_json": str(path),
            }
        )
    return pd.DataFrame(rows)


def load_section(results_dir: Path, manifest: dict | None, key: str, ir) -> pd.DataFrame | None:
    return ir.load_section_frames(
        results_dir,
        manifest,
        key,
        load_json=lambda p, m: load_zebrac_json(p, m, ir),
        load_metadata=ir.load_metadata,
    )


def enrich_size(df: pd.DataFrame | None) -> pd.DataFrame | None:
    if df is None or df.empty:
        return df
    out = df.copy()
    out["size_mb"] = out["workload"].astype(str).str.replace("mb", "", regex=False).astype(float)
    return out


def enrich_seqs(df: pd.DataFrame | None) -> pd.DataFrame | None:
    if df is None or df.empty:
        return df
    out = df.copy()
    out["seq_count"] = out["workload"].astype(int)
    return out


def _std_col(value_col: str) -> str | None:
    if value_col == "mean":
        return "stddev"
    if value_col == "peak_rss_mb":
        return "peak_rss_stddev_mb"
    if value_col == "minor_faults":
        return "minor_faults_stddev"
    return None


def _format_ratio(ratio: float | None) -> str:
    """ASCII ratio label (WRITING.md: no multiplication sign)."""
    if ratio is None:
        return "n/a"
    if abs(ratio - 1.0) < 1e-9:
        return "1x"
    if ratio >= 100:
        return f"{ratio:.0f}x"
    if ratio >= 10:
        return f"{ratio:.1f}x"
    if ratio >= 1:
        return f"{ratio:.2f}x"
    return f"{ratio:.3f}x"


def _bar_patches(tools: list[str], ir, display_map: dict[str, str] | None = None) -> list:
    patches = []
    names = display_map or {}
    for tool in tools:
        color = ir.COLORS.get(tool, "#888888")
        label = names.get(tool) or ir.display_tool(tool)
        kw: dict = {
            "facecolor": color,
            "edgecolor": color if tool in REFERENCE_TOOLS else "none",
            "label": label,
            "alpha": 0.75 if tool in REFERENCE_TOOLS else 0.88,
        }
        if tool in REFERENCE_TOOLS:
            kw["hatch"] = "///"
        patches.append(mpatches.Patch(**kw))
    return patches


def _annotate_ratios(
    ax,
    dataset: str,
    tools: list[str],
    bar_tops: dict[tuple[str, str], tuple[float, float]],
    comparisons: pd.DataFrame,
    baseline: str,
    width: float,
    *,
    rotation: int = 0,
    fontsize: int = 9,
    pad: float = 0.28,
    y_offset: int = 9,
) -> None:
    base_color = STATS_COLORS.get(baseline, "#F7A41D")
    base_key = (dataset, baseline)
    if base_key in bar_tops:
        bx, by = bar_tops[base_key]
        xs = [bar_tops[(dataset, t)][0] for t in tools if (dataset, t) in bar_tops]
        if xs:
            ax.hlines(
                by,
                min(xs) - width * 0.65,
                max(xs) + width * 0.65,
                colors=base_color,
                linestyles=(0, (4, 3)),
                linewidth=1.0,
                alpha=0.45,
                zorder=1,
            )
        ax.annotate(
            "1x",
            (bx, by),
            xytext=(0, y_offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=rotation,
            fontsize=fontsize,
            fontweight="bold",
            color=base_color,
            bbox=dict(
                boxstyle=f"round,pad={pad}",
                facecolor="white",
                edgecolor=base_color,
                linewidth=0.9,
            ),
            zorder=5,
        )
    if comparisons.empty:
        return
    for row in comparisons.itertuples(index=False):
        key = (dataset, row.tool)
        if key not in bar_tops:
            continue
        xpos, mean_t = bar_tops[key]
        color = STATS_COLORS.get(row.tool, "#666666")
        ax.annotate(
            _format_ratio(row.ratio),
            (xpos, mean_t),
            xytext=(0, y_offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            rotation=rotation,
            fontsize=fontsize,
            fontweight="bold",
            color="#111111",
            bbox=dict(
                boxstyle=f"round,pad={pad}",
                facecolor="white",
                edgecolor=color,
                linewidth=0.9,
                alpha=0.96,
            ),
            zorder=5,
        )


def _draw_metric_facet(
    ax,
    work: pd.DataFrame,
    tools: list[str],
    datasets: list[str],
    *,
    value_col: str,
    ylabel: str,
    value_floor: float,
    log_y: bool,
    baseline: str,
    annotate: bool,
    ir,
) -> None:
    """One metric facet: datasets on x-axis, colored bars per tool (GET RC layout)."""
    std_col = _std_col(value_col)
    peers = peer_tools(tools, baseline)
    width = min(0.095, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for ti, tool in enumerate(tools):
        color = ir.COLORS.get(tool, "#888888")
        is_ref = tool in REFERENCE_TOOLS
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            val = max(float(row[value_col].iloc[0]), value_floor)
            std = 0.0
            if std_col and std_col in row.columns and pd.notna(row[std_col].iloc[0]):
                std = max(float(row[std_col].iloc[0]), 0.0)
            xpos = di + (ti - len(tools) / 2 + 0.5) * width
            bar_kw: dict = {
                "width": width,
                "color": color,
                "alpha": 0.75 if is_ref else 0.88,
                "zorder": 2,
            }
            if is_ref:
                bar_kw["hatch"] = "///"
                bar_kw["edgecolor"] = color
                bar_kw["linewidth"] = 0.6
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

    if annotate and peers:
        for ds in datasets:
            ds_work = work[work["dataset"] == ds]
            present = [t for t in tools if (ds, t) in bar_tops]
            if not present:
                continue
            comparisons = ir.build_ratio_comparisons(
                ds_work,
                value_col,
                baseline=baseline,
                peer_tools=[t for t in peers if t in present],
            )
            _annotate_ratios(
                ax,
                ds,
                present,
                bar_tops,
                comparisons,
                baseline,
                width,
                rotation=90,
                fontsize=7,
                pad=0.18,
                y_offset=6,
            )

    if log_y:
        ax.set_yscale("log")
        lo, hi = ax.get_ylim()
        if hi > lo > 0:
            ax.set_ylim(lo, hi * 2.0)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(axis="y", alpha=0.28, which="both")
    ax.set_axisbelow(True)


def fig_metric_facets(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    ir,
    *,
    title: str,
    fig_note: str,
    baseline: str = BASELINE,
    annotate: bool = True,
    display_map: dict[str, str] | None = None,
) -> Path:
    """Datasets on x-axis; facets = wall time, peak RSS, page faults (GET RC style)."""
    filtered = ir.filter_tools(work, tools)
    legend_tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in DATASET_ORDER if d in filtered["dataset"].unique()]
    facets = [
        ("mean", "Wall Time (s)", True, 1e-6),
        ("peak_rss_mb", "Peak RSS (MB)", False, 0.1),
        ("minor_faults", "Minor Page Faults", True, 1.0),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6.8), squeeze=False)
    for ax, (value_col, ylabel, log_y, floor) in zip(axes[0], facets):
        _draw_metric_facet(
            ax,
            filtered,
            legend_tools,
            datasets,
            value_col=value_col,
            ylabel=ylabel,
            value_floor=floor,
            log_y=log_y,
            baseline=baseline,
            annotate=annotate,
            ir=ir,
        )

    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.93,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    patches = _bar_patches(legend_tools, ir, display_map)
    fig.legend(
        handles=patches,
        labels=[p.get_label() for p in patches],
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.06),
        ncol=min(len(legend_tools), 5),
        frameon=False,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.18, top=0.86, wspace=FACET_WSPACE)
    return ir._save(fig, out)


def _scaling_style(tool: str, *, baseline: str = BASELINE) -> dict:
    is_ref = tool in REFERENCE_TOOLS
    is_fai = tool == "z-fasta-full-fai"
    return {
        "color": STATS_COLORS.get(tool, "#888888"),
        "marker": "D" if is_fai else "o",
        "linestyle": "--" if is_ref else "-",
        "linewidth": 2.4 if tool == baseline else (1.9 if is_fai else 1.8),
        "markersize": 5.5 if is_fai else 6,
        "alpha": 0.92 if is_fai else (0.75 if is_ref else 0.95),
        "zorder": 3 if tool == baseline else (2.6 if is_fai else 2),
    }


def _scaling_legend_patches(tools: list[str], ir, display_map: dict[str, str] | None = None) -> list:
    """Line handles matching scaling plot markers/styles (not bar patches)."""
    names = display_map or {}
    handles = []
    for tool in tools:
        st = _scaling_style(tool)
        label = names.get(tool) or ir.display_tool(tool)
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=st["color"],
                marker=st["marker"],
                linestyle=st["linestyle"],
                linewidth=st["linewidth"],
                markersize=st["markersize"],
                markeredgecolor="white",
                markeredgewidth=0.6,
                alpha=st["alpha"],
                label=label,
            )
        )
    return handles


SCALING_METRICS = (
    ("mean", "Wall time (s)", True, 1e-6),
    ("peak_rss_mb", "Peak RSS (MB)", True, 0.1),  # log so peers stay visible vs mmap full
    ("minor_faults", "Minor page faults", True, 1.0),
)


def _plot_scaling_abs_lines(ax, df: pd.DataFrame, tools: list[str], param_col: str, value_col: str, floor: float, log_y: bool, baseline: str) -> None:
    for tool in tools:
        tdf = df[df["tool"] == tool].sort_values(param_col)
        if tdf.empty:
            continue
        st = _scaling_style(tool, baseline=baseline)
        ax.plot(
            tdf[param_col].astype(float),
            tdf[value_col].astype(float).clip(lower=floor),
            marker=st["marker"],
            linestyle=st["linestyle"],
            color=st["color"],
            linewidth=st["linewidth"],
            markersize=st["markersize"],
            markeredgecolor="white",
            markeredgewidth=0.6,
            alpha=st["alpha"],
            zorder=st["zorder"],
            label=tool,
        )
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.grid(alpha=0.28, which="both")
    ax.set_axisbelow(True)


def _finish_scaling_figure(fig, axes_row, tools, ir, title: str, note: str, xlabel: str) -> None:
    for ax in axes_row:
        ax.set_xlabel(xlabel, fontsize=10)
        ax.tick_params(axis="both", labelsize=8.5)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.18, top=0.84, wspace=FACET_WSPACE)
    # Title + note left-aligned to first panel y-axis.
    pos = axes_row[0].get_position()
    fig.text(pos.x0, 0.965, title, ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(pos.x0, 0.922, note, ha="left", va="top", fontsize=8.5, color="#444444", style="italic")
    saved = dict(ir.DISPLAY_NAMES)
    ir.DISPLAY_NAMES.update(SCALING_DISPLAY)
    patches = _scaling_legend_patches(tools, ir, SCALING_DISPLAY)
    ir.DISPLAY_NAMES.clear()
    ir.DISPLAY_NAMES.update(saved)
    fig.legend(
        handles=patches,
        labels=[p.get_label() for p in patches],
        fontsize=8.5,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.08),
        ncol=min(len(tools), 6),
        frameon=False,
        columnspacing=1.1,
        handletextpad=0.4,
    )


def fig_scaling_abs_facets(
    df: pd.DataFrame,
    out: Path,
    *,
    param_col: str,
    xlabel: str,
    title: str,
    tools: list[str],
    ir,
    baseline: str = BASELINE,
) -> Path:
    """1x3 absolute facets: wall | RSS | faults (log-log lines)."""
    work = ir.filter_tools(df, tools)
    tools = [t for t in tools if t in work["tool"].unique()]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.4), squeeze=False)
    for ax, (col, ylab, log_y, floor) in zip(axes[0], SCALING_METRICS):
        _plot_scaling_abs_lines(ax, work, tools, param_col, col, floor, log_y, baseline)
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_title(ylab, fontsize=11, fontweight="bold", pad=8)
    _finish_scaling_figure(
        fig,
        axes[0],
        tools,
        ir,
        title,
        "Log-log absolute values. Diamond = z-fasta (fai). Dashed = seqtk (comp) reference.",
        xlabel,
    )
    return ir._save(fig, out)


def fig_scaling_slope_ratios(
    df: pd.DataFrame,
    out: Path,
    *,
    param_col: str,
    title: str,
    tools: list[str],
    ir,
    baseline: str = BASELINE,
) -> Path:
    """Min->max slope chart of x vs full (does the gap grow with scale?)."""
    work = ir.filter_tools(df, tools)
    tools = [t for t in tools if t in work["tool"].unique() and t != baseline]
    xs = sorted(work[param_col].unique())
    if len(xs) < 2:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "Need ≥2 scale points", ha="center")
        return ir._save(fig, out)
    x0, x1 = xs[0], xs[-1]
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.2), squeeze=False)
    for ax, (col, ylab, _log_y, floor) in zip(axes[0], SCALING_METRICS):
        b0 = work[(work["tool"] == baseline) & (work[param_col] == x0)]
        b1 = work[(work["tool"] == baseline) & (work[param_col] == x1)]
        if b0.empty or b1.empty:
            continue
        base0 = float(b0[col].iloc[0])
        base1 = float(b1[col].iloc[0])
        for tool in tools:
            r0 = work[(work["tool"] == tool) & (work[param_col] == x0)]
            r1 = work[(work["tool"] == tool) & (work[param_col] == x1)]
            if r0.empty or r1.empty:
                continue
            v0 = max(float(r0[col].iloc[0]), floor) / max(base0, floor)
            v1 = max(float(r1[col].iloc[0]), floor) / max(base1, floor)
            st = _scaling_style(tool, baseline=baseline)
            ax.plot(
                [0, 1],
                [v0, v1],
                color=st["color"],
                linewidth=2.0,
                marker=st["marker"],
                markersize=7,
                markeredgecolor="white",
                markeredgewidth=0.7,
                alpha=0.92,
            )
            ax.text(
                1.04,
                v1,
                SCALING_DISPLAY.get(tool, ir.display_tool(tool)),
                fontsize=7.5,
                va="center",
                color=st["color"],
            )
        ax.axhline(1.0, color=STATS_COLORS.get(baseline, "#F7A41D"), linestyle=(0, (4, 3)), linewidth=1.2, alpha=0.7)
        ax.set_yscale("log")
        ax.set_xlim(-0.08, 1.55)
        if param_col == "size_mb":
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"{int(x0)} MB", f"{int(x1)} MB"])
        else:
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"{int(x0):,}", f"{int(x1):,}"])
        ax.set_ylabel("x vs z-fasta", fontsize=10)
        ax.set_title(ylab, fontsize=11, fontweight="bold", pad=8)
        ax.tick_params(axis="both", labelsize=8.5)
        ax.grid(axis="y", alpha=0.28, which="both")
        ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.12, top=0.84, wspace=0.32)
    pos = axes[0, 0].get_position()
    fig.text(pos.x0, 0.965, title, ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(
        pos.x0,
        0.922,
        "Left = smallest fixture; right = largest. Rising = peer worsens vs z-fasta with scale. Gold dashed = 1x.",
        ha="left",
        va="top",
        fontsize=8.5,
        color="#444444",
        style="italic",
    )
    return ir._save(fig, out)



def fmt_wall(row) -> str:
    mean = float(row["mean"])
    std = float(row["stddev"]) if pd.notna(row.get("stddev")) else 0.0
    if mean < 0.01:
        return f"{mean * 1000:.2f}±{std * 1000:.2f} ms"
    return f"{mean:.3f}±{std:.3f} s"


def fmt_rss(row) -> str:
    mean = float(row["peak_rss_mb"])
    std = float(row["peak_rss_stddev_mb"]) if pd.notna(row.get("peak_rss_stddev_mb")) else 0.0
    return f"{mean:.1f}±{std:.1f} MB"


def fmt_faults(row) -> str:
    mean = float(row["minor_faults"])
    std = float(row["minor_faults_stddev"]) if pd.notna(row.get("minor_faults_stddev")) else 0.0
    return f"{mean:.0f}±{std:.0f}"


def _join(parts: list[str]) -> str:
    return "\n".join(parts)


def md_ratio_details(
    work: pd.DataFrame,
    tools: list[str],
    value_col: str,
    *,
    baseline: str,
    table_num: int,
    zf_label: str,
    comp_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    ir,
    summary: str,
) -> str:
    peers = peer_tools(tools, baseline)
    comparisons = ir.build_ratio_comparisons(
        work,
        value_col,
        baseline=baseline,
        peer_tools=peers,
    )
    body = ir.md_zfasta_vs_ratio_table(
        comparisons,
        zf_label=zf_label,
        comp_label=comp_label,
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
    )
    return _join(
        [
            f"<details><summary><strong>Table {table_num}:</strong> {summary}</summary>",
            "",
            body,
            "",
            "</details>",
        ]
    )


def md_metric_tables(
    *,
    heading: str,
    intro: str,
    work: pd.DataFrame,
    tools: list[str],
    value_col: str,
    fmt,
    ir,
    nums,
    baseline: str,
    ratio_summary: str,
    zf_label: str,
    comp_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    display_map: dict[str, str] | None = None,
) -> str:
    """Heading + wide table + collapsible ratio details (no figure)."""
    tools = tools_in_run(work, tools)
    t_main = nums.next_table()
    saved = dict(ir.DISPLAY_NAMES)
    if display_map:
        ir.DISPLAY_NAMES.update(display_map)
    pivot = ir.md_tool_pivot_table(
        work,
        tools,
        "dataset",
        value_col,
        fmt,
        sort_index=dataset_sort_key,
        empty_msg="_No data._",
    )
    if isinstance(pivot, str):
        pivot = pivot.replace("| nan ", "|  ").replace("| nan |", "|  |")
    ir.DISPLAY_NAMES.clear()
    ir.DISPLAY_NAMES.update(saved)

    t_ratio = nums.next_table()
    ratio_block = md_ratio_details(
        work,
        tools,
        value_col,
        baseline=baseline,
        table_num=t_ratio,
        zf_label=zf_label,
        comp_label=comp_label,
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
        ir=ir,
        summary=ratio_summary,
    )
    return _join(
        [
            f"## {heading}",
            "",
            intro,
            "",
            f"**Table {t_main}:** Mean ± stddev by dataset and tool.",
            "",
            pivot,
            "",
            ratio_block,
            "",
            '<div style="margin: 1.5em 0"></div>',
        ]
    )


def md_combined_metric_figure(
    *,
    nums,
    fig_rel: str,
    fig_caption: str,
    reading: list[str],
) -> str:
    """Single wall|RSS|faults figure after metric tables (GET RC grammar)."""
    fnum = nums.next_figure()
    reading_md = "\n".join(f"- {b}" for b in reading)
    return _join(
        [
            f"**Figure {fnum}:** {fig_caption}.",
            "",
            f"![Figure {fnum}]({fig_rel})",
            "",
            f"**Reading Figure {fnum}**",
            "",
            reading_md,
            "",
            '<div style="margin: 1.5em 0"></div>',
        ]
    )


def md_overview(manifest: dict, results_dir: Path) -> str:
    verify_skipped, verify = resolve_verify_status(manifest, results_dir)
    if verify_skipped:
        verify_line = "Verify was skipped for this run (`--skip-tests` / `--skip-verify`)."
    elif verify is not None:
        verify_line = (
            f"Correctness reported **{verify}** passing checks before perf "
            "(see `bench/stats/run.sh`)."
        )
    else:
        verify_line = "Verify status was not recorded in the manifest."
    return _join(
        [
            "This report times `z-fasta stats` against clean-FASTA peers on the shared REAL "
            "datasets (Genome, Transcriptome, Proteome) and on synthetic size / sequence-count sweeps.",
            "",
            "**What is timed**",
            "",
            "- **Full stats (peers):** `z-fasta stats` with `.zfi` (lengths + composition scan) vs "
            "noodles / rust-bio wrappers, seqkit `stats -a`, and seqtk `comp` (nucleotide only; "
            "hatched reference).",
            "- **Scaling:** complete stats and composition peers: wall / RSS / page faults vs "
            "file size and vs sequence count (absolute facets + min->max x slopes; x tables in details).",
            "",
            "**Index policy:** sidecars are preloaded once. Timed commands only run `stats` / peer "
            "tools; index *build* is outside this wall (see `bench/index/REPORT.md`).",
            "",
            verify_line,
            "",
            "Correctness values (N50, GC, ...) are gated by verify against BioPython "
            "(`bench/stats/oracle.py`), not re-tabulated here.",
        ]
    )


def md_field_matrix() -> str:
    """Capability overview for complete stats."""
    return _join(
        [
            "`z-fasta stats` reads length metadata from the selected sidecar, then scans "
            "sequence bytes for composition. Peers re-parse the FASTA every time.",
            "",
            "### Assembly metrics (lengths / index)",
            "",
            "| Metric | z-fasta | noodles | rust-bio | seqkit `-a` | seqtk `comp` |",
            "| --- | --- | --- | --- | --- | --- |",
            "| Sequences / total bases | yes | yes | yes | yes | totals |",
            "| Min / max / mean / median | yes | yes | yes | yes | - |",
            "| N50 / L50 | yes | yes | yes | yes | - |",
            "| N90 / L90 / AU | yes | yes | yes | - | - |",
            "| Shortest / longest names | yes | yes | yes | - | - |",
            "",
            "### Composition (sequence scan)",
            "",
            "| Metric | z-fasta | noodles | rust-bio | seqkit `-a` | seqtk `comp` |",
            "| --- | --- | --- | --- | --- | --- |",
            "| GC / GC skew / soft-mask | yes | yes | yes | GC | partial |",
            "| A / C / G / T / N / Other | yes | yes | yes | - | yes |",
            "| Top amino acids (protein) | yes | yes | yes | - | - |",
            "",
            "**How to read the benches:** Full-stats peers are fair for composition work. "
            "seqkit also reports Q1/Q3/gaps/Q20 (FASTA/Q QC); those are out of scope for z-fasta. "
            "Wrappers are clean-FASTA peers only (no messy / side-table path).",
            "",
            "Oracle: [`bench/stats/oracle.py`](oracle.py). Verify: "
            "[`bench/stats/run.sh`](run.sh) correctness (`run_tests`).",
        ]
    )


def md_run_provenance(manifest: dict, ir) -> str:
    tools = manifest.get("tools") or {}
    tool_lines = []
    for name in ("seqkit", "seqtk", "noodles", "rustbio", "samtools"):
        if name in tools and tools[name]:
            tool_lines.append(f"- **{name}:** {tools[name]}")
    sections = manifest.get("sections") or {}
    sec_lines = [
        f"- `{key}` -> `{sections[key]}`"
        for key in ("perf_full", "scale_size", "scale_seqs_fixed")
        if key in sections
    ]
    parts = [
        f"- **Timestamp:** `{manifest.get('timestamp', '')}`",
        f"- **Runner:** {manifest.get('runner', 'zebrac')} ({manifest.get('mode', 'warm')})",
        f"- **zebrac:** {manifest.get('zebrac', '')}",
        f"- **z-fasta:** {manifest.get('z_fasta', '')}",
        f"- **Samples:** runs={manifest.get('runs')}, warmup={manifest.get('warmup')}, "
        f"duration_ms={manifest.get('duration_ms')}",
        f"- **Metadata:** `{manifest.get('metadata', '')}`",
        f"- **Index preload:** {manifest.get('index_preload', True)} "
        "(build time not included in stats wall)",
        "",
        "**Sections**",
        "",
    ]
    parts += sec_lines or ["- _(none)_"]
    parts += ["", "**Peer versions**", ""]
    parts += tool_lines or ["- _(none recorded)_"]
    parts += ["", ir.md_data_used(PROJECT_ROOT)]
    return _join(parts)


def _scaling_pivot(df: pd.DataFrame, tools: list[str], param_col: str, param_order: list, xlabel: str, fmt, ir) -> str:
    work = ir.filter_tools(df, tools).copy()
    work["cell"] = work.apply(fmt, axis=1)
    pivot = work.pivot(index=param_col, columns="tool", values="cell")
    pivot = pivot.reindex([p for p in param_order if p in pivot.index])
    cols = [c for c in tools if c in pivot.columns]
    pivot = pivot[cols]
    pivot = pivot.rename(columns={c: SCALING_DISPLAY.get(c, ir.display_tool(c)) for c in pivot.columns})
    if param_col == "size_mb":
        pivot.index = [f"{int(v)} MB" for v in pivot.index]
    else:
        pivot.index = [f"{int(v):,}" for v in pivot.index]
    pivot.index.name = xlabel
    return pivot.to_markdown()


def md_scaling_section(
    df: pd.DataFrame,
    *,
    param_col: str,
    param_order: list,
    xlabel: str,
    section_title: str,
    fig_stem: str,
    fig_title: str,
    sort_key,
    ir,
    nums,
    figures_dir: Path,
) -> str:
    """Absolute + slope figures; wall table open; RSS/faults + x tables in details."""
    tools = tools_in_run(df, SCALING_TOOLS)
    saved_names = dict(ir.DISPLAY_NAMES)
    ir.DISPLAY_NAMES.update(SCALING_DISPLAY)
    abs_name = f"{fig_stem}.png"
    slope_name = f"{fig_stem}_slope.png"

    fig_scaling_abs_facets(
        df,
        figures_dir / abs_name,
        param_col=param_col,
        xlabel=xlabel,
        title=fig_title,
        tools=tools,
        ir=ir,
    )
    fig_scaling_slope_ratios(
        df,
        figures_dir / slope_name,
        param_col=param_col,
        title=f"{fig_title}: x gap (smallest to largest)",
        tools=tools,
        ir=ir,
    )

    t_wall = nums.next_table()
    wall_md = _scaling_pivot(df, tools, param_col, param_order, xlabel, fmt_wall, ir)

    peers = peer_tools(tools, BASELINE)
    label_group = (lambda v: f"{int(v)} MB" if param_col == "size_mb" else f"{int(v):,}")
    t_wall_ratio = nums.next_table()
    wall_ratio = ir.md_zfasta_vs_ratio_table(
        ir.build_ratio_comparisons(
            df,
            "mean",
            baseline=BASELINE,
            peer_tools=peers,
            group_col=param_col,
            group_sort=sort_key,
            label_group=label_group,
        ),
        group_label=xlabel,
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="Time x",
        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
    )

    t_rss = nums.next_table()
    rss_md = _scaling_pivot(df, tools, param_col, param_order, xlabel, fmt_rss, ir)
    t_rss_ratio = nums.next_table()
    rss_ratio = ir.md_zfasta_vs_ratio_table(
        ir.build_ratio_comparisons(
            df,
            "peak_rss_mb",
            baseline=BASELINE,
            peer_tools=peers,
            group_col=param_col,
            group_sort=sort_key,
            label_group=label_group,
        ),
        group_label=xlabel,
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="RSS x",
        fmt_zf=lambda r: f"{r.zfasta_v:.1f} MB",
        fmt_comp=lambda r: f"{r.comp_v:.1f} MB",
    )

    t_pf = nums.next_table()
    pf_md = _scaling_pivot(df, tools, param_col, param_order, xlabel, fmt_faults, ir)
    t_pf_ratio = nums.next_table()
    pf_ratio = ir.md_zfasta_vs_ratio_table(
        ir.build_ratio_comparisons(
            df,
            "minor_faults",
            baseline=BASELINE,
            peer_tools=peers,
            group_col=param_col,
            group_sort=sort_key,
            label_group=label_group,
        ),
        group_label=xlabel,
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="Faults x",
        fmt_zf=lambda r: f"{r.zfasta_v:.0f}",
        fmt_comp=lambda r: f"{r.comp_v:.0f}",
    )

    f_abs = nums.next_figure()
    f_slope = nums.next_figure()

    ir.DISPLAY_NAMES.clear()
    ir.DISPLAY_NAMES.update(saved_names)

    return _join(
        [
            f"## {section_title}",
            "",
            "Synthetic FASTAs under `bench/shared/cache/scaling/` (generated on demand). "
            "Indexes preloaded; timed work is complete `stats` / peers only. "
            "Two figures per sweep: absolute wall / RSS / faults, then whether the "
            "x vs z-fasta gap grows from the smallest fixture to the largest. "
            "Per-point x values stay in the collapsible tables.",
            "",
            f"**Table {t_wall}:** Mean ± stddev wall time by {xlabel} and tool.",
            "",
            wall_md,
            "",
            f"<details><summary><strong>Table {t_wall_ratio}:</strong> "
            "Time x = peer / z-fasta at each point.</summary>",
            "",
            wall_ratio,
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_rss}:</strong> Peak RSS (MB) by {xlabel}.</summary>",
            "",
            rss_md,
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_rss_ratio}:</strong> RSS x = peer / z-fasta.</summary>",
            "",
            rss_ratio,
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_pf}:</strong> Minor page faults by {xlabel}.</summary>",
            "",
            pf_md,
            "",
            "</details>",
            "",
            f"<details><summary><strong>Table {t_pf_ratio}:</strong> Faults x = peer / z-fasta.</summary>",
            "",
            pf_ratio,
            "",
            "</details>",
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
            f"**Figure {f_abs}:** Absolute wall / RSS / page faults vs {xlabel} (log-log).",
            "",
            f"![Figure {f_abs}](results/figures/{abs_name})",
            "",
            f"**Reading Figure {f_abs}**",
            "",
            "- **Facets:** wall time | peak RSS | minor page faults (log-log).",
            f"- **X:** {xlabel}. Gold = z-fasta (`.zfi`); diamond = z-fasta (`.fai`); dashed = seqtk (comp).",
            "- RSS is log-scaled so peer lines stay readable next to mmap-full.",
            "- Per-point x vs z-fasta: open the Time / RSS / Faults x tables above.",
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
            f"**Figure {f_slope}:** x vs z-fasta at the smallest vs largest fixture.",
            "",
            f"![Figure {f_slope}](results/figures/{slope_name})",
            "",
            f"**Reading Figure {f_slope}**",
            "",
            "- Left = smallest fixture; right = largest. Gold dashed = 1x.",
            "- Rising = peer gets relatively worse with scale; falling = relatively better.",
            "",
            '<div style="margin: 1.5em 0"></div>',
        ]
    )



def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", nargs="?", default=str(RESULTS_DIR))
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    ir = load_index_report()
    # Keep index helper dataset sort aligned with stats facet order.
    ir.DATASET_ORDER = list(DATASET_ORDER)

    manifest = ir.load_run_manifest(results_dir)
    if not manifest:
        raise SystemExit(f"No run_*.json under {results_dir}")
    incomplete = ir.is_incomplete_run(
        manifest,
        required_sections=("perf_full", "scale_size", "scale_seqs_fixed"),
        skip_flags=("skip_full", "skip_scale_size", "skip_scale_seqs_fixed"),
    )
    if incomplete and not args.allow_incomplete:
        raise SystemExit(
            "Refusing to overwrite REPORT.md from incomplete run "
            "(missing sections or skip_* flags). Pass --allow-incomplete for a draft."
        )

    full_df = load_section(results_dir, manifest, "perf_full", ir)
    size_df = enrich_size(load_section(results_dir, manifest, "scale_size", ir))
    seq_df = enrich_seqs(load_section(results_dir, manifest, "scale_seqs_fixed", ir))

    if full_df is None and size_df is None and seq_df is None:
        raise SystemExit("No section data found")

    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    ir.prune_stale_pngs(figures_dir, REPORT_FIGURES)

    nums = ir.ReportCounters()
    report_lines: list[str] = [
        "<!-- markdownlint-disable MD024 MD032 MD033 MD036 MD041 MD049 -->",
        "# z-fasta Stats Benchmark Report",
        "_Auto-generated by `bench/stats/generate_report.py` from zebrac results._",
    ]
    if incomplete:
        report_lines.append(
            "> **DRAFT:** incomplete manifest (skipped sections and/or low sample counts). "
            "Numbers may be partial."
        )

    report_lines.extend(
        [
            "## Overview",
            md_overview(manifest, results_dir),
            "## What each tool reports",
            md_field_matrix(),
            "## Run Provenance",
            md_run_provenance(manifest, ir),
        ]
    )

    # --- Full peers ---
    if full_df is not None and not full_df.empty:
        full_tools = tools_in_run(full_df, FULL_TOOLS)
        sample_n = manifest.get("runs", "?")
        fig_metric_facets(
            full_df,
            figures_dir / "perf_full.png",
            full_tools,
            ir,
            title="Full stats: wall time, peak RSS, page faults",
            fig_note=(
                f"X = Genome / Transcriptome / Proteome. Error bars = zebrac stddev "
                f"(n={sample_n}). Hatched = seqtk (comp) reference (omitted on Proteome). "
                "Rotated labels = peer / z-fasta (all facets)."
            ),
            baseline=BASELINE,
            annotate=True,
            display_map=FULL_SECTION_DISPLAY,
        )

        report_lines.append("## Performance: Full stats")
        report_lines.append(
            "Peers re-parse the FASTA. z-fasta full loads `.zfi`, emits length metrics from the "
            "index, then mmap-scans sequence bytes for composition. Baseline for ratios: "
            "**z-fasta full**. **Figure facets** are wall time, peak RSS, and minor page "
            "faults; **x-axis** is dataset; **bar color** is tool (GET RC layout)."
        )

        t_wall = nums.table
        report_lines.append(
            md_metric_tables(
                heading="Wall time",
                intro=(
                    "Zebrac mean wall time (seconds) for one `stats` / peer process with "
                    "stdout discarded."
                ),
                work=full_df,
                tools=full_tools,
                value_col="mean",
                fmt=fmt_wall,
                ir=ir,
                nums=nums,
                baseline=BASELINE,
                ratio_summary="Time x = peer wall / z-fasta full. Same ratios as bar labels.",
                zf_label="z-fasta",
                comp_label="Peer",
                ratio_label="Time x",
                fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                fmt_comp=lambda r: f"{r.comp_v:.4f}s",
                display_map=FULL_SECTION_DISPLAY,
            )
        )
        t_rss = nums.table
        report_lines.append(
            md_metric_tables(
                heading="Memory Usage",
                intro="Peak RSS from zebrac samples for the same full-stats commands.",
                work=full_df,
                tools=full_tools,
                value_col="peak_rss_mb",
                fmt=fmt_rss,
                ir=ir,
                nums=nums,
                baseline=BASELINE,
                ratio_summary="RSS x = peer peak RSS / z-fasta full.",
                zf_label="z-fasta",
                comp_label="Peer",
                ratio_label="RSS x",
                fmt_zf=lambda r: f"{r.zfasta_v:.1f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.1f} MB",
                display_map=FULL_SECTION_DISPLAY,
            )
        )
        t_pf = nums.table
        report_lines.append(
            md_metric_tables(
                heading="Page Faults",
                intro="Minor page faults for the same full-stats commands.",
                work=full_df,
                tools=full_tools,
                value_col="minor_faults",
                fmt=fmt_faults,
                ir=ir,
                nums=nums,
                baseline=BASELINE,
                ratio_summary="Faults x = peer minor faults / z-fasta full.",
                zf_label="z-fasta",
                comp_label="Peer",
                ratio_label="Faults x",
                fmt_zf=lambda r: f"{r.zfasta_v:.0f}",
                fmt_comp=lambda r: f"{r.comp_v:.0f}",
                display_map=FULL_SECTION_DISPLAY,
            )
        )
        report_lines.append(
            md_combined_metric_figure(
                nums=nums,
                fig_rel="results/figures/perf_full.png",
                fig_caption=(
                    f"Tables {t_wall}, {t_rss}, and {t_pf} as grouped bars "
                    "(datasets on x-axis; tool color from legend)"
                ),
                reading=[
                    "**Facets:** wall time (log y) | peak RSS (linear) | minor page faults (log y).",
                    "**X-axis:** Genome, Transcriptome, Proteome.",
                    "**Bar colors:** one lane per tool (legend below figure).",
                    "**Bar labels (rotated):** `1x` on z-fasta; peers = peer / z-fasta. "
                    "Same ratio rules on wall, RSS, and page faults.",
                    "**Hatched bars:** seqtk `comp` (composition-only reference; omitted on Proteome).",
                    "Gold = z-fasta full (`.zfi`). Error bars = 1σ across measured samples.",
                ],
            )
        )

    # --- Scaling (composition peers only) ---
    if size_df is not None and not size_df.empty:
        report_lines.append(
            md_scaling_section(
                size_df,
                param_col="size_mb",
                param_order=SIZE_MB_ORDER,
                xlabel="File size (MB)",
                section_title="Scaling: file size",
                fig_stem="scaling_size",
                fig_title="Scaling by file size",
                sort_key=lambda v: float(v),
                ir=ir,
                nums=nums,
                figures_dir=figures_dir,
            )
        )
    if seq_df is not None and not seq_df.empty:
        report_lines.append(
            md_scaling_section(
                seq_df,
                param_col="seq_count",
                param_order=SEQ_COUNT_ORDER,
                xlabel="Sequence count",
                section_title="Scaling: sequence count (fixed 1024 bp)",
                fig_stem="scaling_seqs_fixed",
                fig_title="Scaling by sequence count (1024 bp/seq)",
                sort_key=lambda v: int(v),
                ir=ir,
                nums=nums,
                figures_dir=figures_dir,
            )
        )

    report_path = SCRIPT_DIR / "REPORT.md"
    report_path.write_text(ir.normalize_markdown(report_lines))
    print(f"Wrote {report_path}")
    print(f"Figures under {figures_dir}")


if __name__ == "__main__":
    main()
