#!/usr/bin/env python3
"""Generate stats benchmark REPORT.md and figures from zebrac JSON.

Mirrors bench/get and bench/index report grammar (Overview, field matrix,
Provenance, wall/RSS/faults, mode 2x2, scaling). Colors and labels follow
plan/stats-bench.md. Prose follows plan/WRITING.md (ASCII; no emojis).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
PROJECT_ROOT = SCRIPT_DIR.parent.parent
FIGURES_DIR = RESULTS_DIR / "figures"

# Facet / table order matches plan/stats-bench.md (Genome | Proteome | Transcriptome).
DATASET_ORDER = ["Genome", "Proteome", "Transcriptome"]

BASELINE = "z-fasta-full"

FULL_TOOLS = [
    "z-fasta-full",
    "noodles",
    "rustbio",
    "seqkit",
    "seqtk",
]

MODE_TOOLS = [
    "z-fasta-full",
    "z-fasta-full-fai",
    "z-fasta-indexed-zfi",
    "z-fasta-indexed-fai",
]

SCALING_TOOLS = [
    "z-fasta-full",
    "z-fasta-full-fai",
    "z-fasta-indexed-zfi",
    "z-fasta-indexed-fai",
    "noodles",
    "rustbio",
    "seqkit",
    "seqtk",
]

REFERENCE_TOOLS = frozenset({"seqtk"})

STATS_COLORS = {
    "z-fasta-full": "#F7A41D",
    "z-fasta-full-fai": "#FFB74D",
    "z-fasta-indexed-zfi": "#C67A00",
    "z-fasta-indexed-fai": "#E65100",
    "seqkit": "#1565C0",
    "noodles": "#C45C26",
    "rustbio": "#8B3A2A",
    "seqtk": "#6A1B9A",
}

STATS_DISPLAY = {
    "z-fasta-full": "z-fasta full",
    "z-fasta-full-fai": "z-fasta full (fai)",
    "z-fasta-indexed-zfi": "z-fasta indexed (zfi)",
    "z-fasta-indexed-fai": "z-fasta indexed (fai)",
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

GROUPED_BAR_WSPACE = 0.10

REPORT_FIGURES = frozenset(
    {
        "perf_full_wall.png",
        "perf_full_rss.png",
        "perf_full_faults.png",
        "perf_mode_wall.png",
        "perf_mode_rss.png",
        "perf_mode_faults.png",
        "scaling_size.png",
        "scaling_seqs_fixed.png",
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
    return mod


def load_manifest(results_dir: Path) -> dict | None:
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


def is_incomplete(manifest: dict | None) -> bool:
    if not manifest:
        return True
    sections = manifest.get("sections") or {}
    for key in ("perf_full", "perf_mode", "scale_size", "scale_seqs_fixed"):
        if key not in sections:
            return True
    for flag in ("skip_full", "skip_mode", "skip_scale_size", "skip_scale_seqs_fixed"):
        if manifest.get(flag):
            return True
    try:
        if int(manifest.get("warmup", 0)) < 1:
            return True
        if int(manifest.get("runs", 0)) < 3:
            return True
    except (TypeError, ValueError):
        return True
    return False


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


def load_zebrac_json(path: Path, metadata_df: pd.DataFrame | None) -> pd.DataFrame:
    data = json.loads(path.read_text())
    meta_rows: list[dict] = []
    if metadata_df is not None and "raw_json" in metadata_df.columns:
        resolved = path.resolve()
        for _, row in metadata_df.iterrows():
            raw = Path(row["raw_json"])
            try:
                if raw.resolve() == resolved:
                    meta_rows.append(row.to_dict())
            except OSError:
                if raw == path:
                    meta_rows.append(row.to_dict())
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
    rel = (manifest or {}).get("sections", {}).get(key)
    if not rel:
        return None
    section_dir = results_dir / rel
    if not section_dir.is_dir():
        return None
    metadata = ir.load_metadata(results_dir, manifest)
    frames = []
    for jf in sorted(section_dir.glob("*.json")):
        frame = load_zebrac_json(jf, metadata)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


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
    if ratio is None:
        return "n/a"
    if abs(ratio - 1.0) < 1e-9:
        return "1×"
    if ratio >= 100:
        return f"{ratio:.0f}×"
    if ratio >= 10:
        return f"{ratio:.1f}×"
    if ratio >= 1:
        return f"{ratio:.2f}×"
    return f"{ratio:.3f}×"


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
                linewidth=1.1,
                alpha=0.55,
                zorder=1,
            )
        ax.annotate(
            "1×",
            (bx, by),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=base_color,
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor=base_color,
                linewidth=1.1,
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
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
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


def fig_dataset_grouped_bars(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    value_col: str,
    ylabel: str,
    title: str,
    fig_note: str,
    ir,
    *,
    log_y: bool = True,
    value_floor: float = 1e-6,
    baseline: str = BASELINE,
    annotate: bool = True,
    display_map: dict[str, str] | None = None,
) -> Path:
    """1x3 Genome|Proteome|Transcriptome facets; one tool cluster per panel."""
    filtered = ir.filter_tools(work, tools)
    tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in DATASET_ORDER if d in filtered["dataset"].unique()]
    std_col = _std_col(value_col)
    peers = peer_tools(tools, baseline)

    fig, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(max(10.5, 3.4 * len(datasets)), 7.2),
        sharey=False,
        gridspec_kw={"wspace": GROUPED_BAR_WSPACE},
    )
    if len(datasets) == 1:
        axes = [axes]

    ylab_set = False
    for ax, ds in zip(axes, datasets):
        ds_work = filtered[filtered["dataset"] == ds]
        present = [t for t in tools if t in ds_work["tool"].unique()]
        if not present:
            ax.set_visible(False)
            continue
        if not ylab_set:
            ax.set_ylabel(ylabel, fontsize=10, labelpad=2)
            ylab_set = True
        n = len(present)
        width = min(0.18, 0.78 / max(1, n))
        bar_tops: dict[tuple[str, str], tuple[float, float]] = {}
        for ti, tool in enumerate(present):
            row = ds_work[ds_work["tool"] == tool]
            if row.empty:
                continue
            val = max(float(row[value_col].iloc[0]), value_floor)
            std = 0.0
            if std_col and std_col in row.columns and pd.notna(row[std_col].iloc[0]):
                std = max(float(row[std_col].iloc[0]), 0.0)
            xpos = (ti - n / 2 + 0.5) * width
            color = ir.COLORS.get(tool, "#888888")
            bar_kw: dict = {
                "width": width * 0.92,
                "color": color,
                "alpha": 0.75 if tool in REFERENCE_TOOLS else 0.88,
                "zorder": 2,
            }
            if tool in REFERENCE_TOOLS:
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
            comparisons = ir.build_ratio_comparisons(
                ds_work,
                value_col,
                baseline=baseline,
                peer_tools=peers,
            )
            # build_ratio_comparisons labels dataset from group; ensure match
            _annotate_ratios(ax, ds, present, bar_tops, comparisons, baseline, width)

        if log_y:
            ax.set_yscale("log")
        ax.set_xticks([])
        ax.set_xlabel("")
        ax.set_title(ds, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.28, which="both")
        ax.set_axisbelow(True)
        ax.margins(x=0.12)

    patches = _bar_patches(tools, ir, display_map)
    labels = [p.get_label() for p in patches]
    fig.subplots_adjust(left=0.06, right=0.995, bottom=0.16, top=0.84)
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.97)
    fig.text(0.5, 0.915, fig_note, ha="center", va="top", fontsize=9, color="#444444", style="italic")
    fig.legend(
        handles=patches,
        labels=labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.12),
        bbox_transform=fig.transFigure,
        ncol=min(len(labels), 4),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def fig_scaling_lines(
    df: pd.DataFrame,
    out: Path,
    *,
    param_col: str,
    xlabel: str,
    title: str,
    fig_note: str,
    tools: list[str],
    ir,
    baseline: str = BASELINE,
) -> Path:
    work = ir.filter_tools(df, tools)
    tools = [t for t in tools if t in work["tool"].unique()]
    fig, ax = plt.subplots(figsize=(14, 7.2))
    for tool in tools:
        tdf = work[work["tool"] == tool].sort_values(param_col)
        if tdf.empty:
            continue
        xs = tdf[param_col].astype(float)
        ys = tdf["mean"].astype(float).clip(lower=1e-6)
        color = ir.COLORS.get(tool, "#888888")
        linewidth = 2.4 if tool == baseline else 1.8
        linestyle = "--" if tool in REFERENCE_TOOLS else "-"
        ax.plot(
            xs,
            ys,
            color=color,
            marker="o",
            linestyle=linestyle,
            linewidth=linewidth,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=ir.display_tool(tool),
            zorder=3 if tool == baseline else 2,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Wall Time (s, log scale)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=20)
    fig.text(0.5, 0.945, fig_note, ha="center", va="top", fontsize=9, color="#444444", style="italic")
    ax.grid(alpha=0.28, which="both")
    ax.set_axisbelow(True)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=min(len(handles), 4),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.22, top=0.88)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


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


def md_metric_block(
    *,
    heading: str,
    intro: str,
    work: pd.DataFrame,
    tools: list[str],
    value_col: str,
    fmt,
    ir,
    nums,
    fig_rel: str,
    fig_caption: str,
    reading: list[str],
    baseline: str,
    ratio_summary: str,
    zf_label: str,
    comp_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    display_map: dict[str, str] | None = None,
) -> str:
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
    # Missing tool/dataset cells come through as the literal "nan" from pandas.
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
    fnum = nums.next_figure()
    reading_md = "\n".join(f"- {b}" for b in reading)
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
            "",
            f"**Figure {fnum}:** {fig_caption}",
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


def md_overview(manifest: dict) -> str:
    verify = manifest.get("verify_pass")
    verify_skipped = manifest.get("verify_skipped", False)
    if verify_skipped:
        verify_line = "Verify was skipped for this run (`--skip-verify`)."
    elif verify is not None:
        verify_line = (
            f"L2 verify reported **{verify}** passing checks before perf "
            "(see `bench/stats/verify.sh`)."
        )
    else:
        verify_line = "Verify status was not recorded in the manifest."
    return _join(
        [
            "This report times `z-fasta stats` against clean-FASTA peers on the shared REAL "
            "datasets (Genome, Proteome, Transcriptome) and on synthetic size / sequence-count sweeps.",
            "",
            "**What is timed**",
            "",
            "- **Full stats (peers):** `z-fasta stats` with `.zfi` (lengths + composition scan) vs "
            "noodles / rust-bio wrappers, seqkit `stats -a`, and seqtk `comp` (nucleotide only; "
            "hatched reference).",
            "- **z-fasta modes (2x2):** full and `--index-only`, each with `.zfi` and with `.fai` "
            "(`.zfi` stashed). Same stats surface for both index formats; messy side tables need `.zfi`.",
            "- **Scaling:** wall time vs file size and vs fixed-length sequence count for the four "
            "z-fasta lanes plus peers.",
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
    return _join(
        [
            "Static capability matrix. Indexed = `z-fasta stats --index-only` (lengths only). "
            "Full = lengths + composition scan.",
            "",
            "| Field / capability | z-fasta full | z-fasta indexed | seqkit `-a` | noodles wrap | rustbio wrap | seqtk `comp` | BioPython oracle |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |",
            "| Sequences / sum_len | yes | yes | yes | yes | yes | yes (rows) | yes |",
            "| min / max / mean / median | yes | yes | yes | yes | yes | no | yes |",
            "| N50 / L50 / N90 / L90 / AU | yes | yes | N50 only | yes | yes | no | yes |",
            "| GC / GC skew / soft-mask | yes (nuc) | no | GC | yes | yes | partial | yes |",
            "| A/C/G/T/N/Other | yes (nuc) | no | no | yes | yes | yes | yes |",
            "| Top AA / lowercase | yes (aa) | no | no | yes | yes | no | yes |",
            "| Q1 / Q3 / gaps / Q20 | no | no | yes | no | no | no | out of scope |",
            "",
            "Peer field gaps: seqkit exposes Q1/Q3/gaps/Q20 (FASTA/Q QC); z-fasta exposes "
            "N90/L90/AU/GC skew in one place. Wrappers are clean-FASTA comparison peers only.",
            "",
            "Oracle: [`bench/stats/oracle.py`](oracle.py). Verify entry: "
            "[`bench/stats/verify.sh`](verify.sh). Layout twins are generated under "
            "`bench/stats/data/verify/layout_twins/`.",
        ]
    )


def md_run_provenance(manifest: dict, ir) -> str:
    tools = manifest.get("tools") or {}
    tool_lines = []
    for name in ("seqkit", "seqtk", "noodles", "rustbio", "samtools"):
        if name in tools and tools[name]:
            tool_lines.append(f"- **{name}:** {tools[name]}")
    sections = manifest.get("sections") or {}
    sec_lines = [f"- `{k}` -> `{v}`" for k, v in sections.items()]
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


def md_composition_tax(mode_df: pd.DataFrame, nums) -> str:
    """(full - indexed_zfi) / full wall fraction per dataset."""
    rows = []
    for ds in DATASET_ORDER:
        full = mode_df[(mode_df["dataset"] == ds) & (mode_df["tool"] == "z-fasta-full")]
        idx = mode_df[(mode_df["dataset"] == ds) & (mode_df["tool"] == "z-fasta-indexed-zfi")]
        if full.empty or idx.empty:
            continue
        fv = float(full["mean"].iloc[0])
        iv = float(idx["mean"].iloc[0])
        frac = (fv - iv) / fv if fv > 0 else None
        rows.append(
            {
                "Dataset": ds,
                "full (s)": f"{fv:.4f}",
                "indexed zfi (s)": f"{iv:.4f}",
                "composition fraction": f"{frac:.2%}" if frac is not None else "n/a",
            }
        )
    if not rows:
        return ""
    t = nums.next_table()
    return _join(
        [
            f"**Table {t}:** Composition tax. Fraction of full wall explained by the mmap scan: "
            "`(full - indexed(zfi)) / full`.",
            "",
            pd.DataFrame(rows).to_markdown(index=False),
        ]
    )


def prune_stale_figures(figures_dir: Path) -> None:
    if not figures_dir.is_dir():
        return
    for path in figures_dir.glob("*.png"):
        if path.name not in REPORT_FIGURES:
            path.unlink()


def md_scaling_section(
    df: pd.DataFrame,
    *,
    param_col: str,
    param_order: list,
    xlabel: str,
    section_title: str,
    fig_name: str,
    fig_title: str,
    sort_key,
    ir,
    nums,
    figures_dir: Path,
) -> str:
    tools = tools_in_run(df, SCALING_TOOLS)
    png = figures_dir / fig_name
    fig_scaling_lines(
        df,
        png,
        param_col=param_col,
        xlabel=xlabel,
        title=fig_title,
        fig_note="Absolute wall times (log-log). Dashed = seqtk reference. Thick gold = z-fasta full.",
        tools=tools,
        ir=ir,
    )
    t_main = nums.next_table()
    work = ir.filter_tools(df, tools).copy()
    work["cell"] = work.apply(fmt_wall, axis=1)
    pivot = work.pivot(index=param_col, columns="tool", values="cell")
    pivot = pivot.reindex([p for p in param_order if p in pivot.index])
    cols = [c for c in tools if c in pivot.columns]
    pivot = pivot[cols]
    pivot = pivot.rename(columns={c: ir.display_tool(c) for c in pivot.columns})
    if param_col == "size_mb":
        pivot.index = [f"{int(v)} MB" for v in pivot.index]
    else:
        pivot.index = [f"{int(v):,}" for v in pivot.index]
    pivot.index.name = xlabel

    t_ratio = nums.next_table()
    peers = peer_tools(tools, BASELINE)
    comparisons = ir.build_ratio_comparisons(
        df,
        "mean",
        baseline=BASELINE,
        peer_tools=peers,
        group_col=param_col,
        group_sort=sort_key,
        label_group=(
            lambda v: f"{int(v)} MB" if param_col == "size_mb" else f"{int(v):,}"
        ),
    )
    fnum = nums.next_figure()
    return _join(
        [
            f"## {section_title}",
            "",
            "Synthetic FASTAs under `bench/stats/data/` (generated on demand). "
            "Indexes preloaded; timed work is `stats` / peers only.",
            "",
            f"**Table {t_main}:** Mean ± stddev wall time by {xlabel.lower()} and tool.",
            "",
            pivot.to_markdown(),
            "",
            f"<details><summary><strong>Table {t_ratio}:</strong> "
            "Time × = other / z-fasta full at each point.</summary>",
            "",
            ir.md_zfasta_vs_ratio_table(
                comparisons,
                group_label=xlabel,
                zf_label="z-fasta full",
                comp_label="Other",
                ratio_label="Time ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                fmt_comp=lambda r: f"{r.comp_v:.4f}s",
            ),
            "",
            "</details>",
            "",
            '<div style="margin: 1.5em 0"></div>',
            "",
            f"**Figure {fnum}:** Wall time vs {xlabel.lower()} (log-log lines).",
            "",
            f"![Figure {fnum}](results/figures/{fig_name})",
            "",
            f"**Reading Figure {fnum}**",
            "",
            f"- X = {xlabel.lower()}; Y = wall seconds (both log).",
            "- Thick gold = z-fasta full; lighter golds = full (fai) / indexed lanes.",
            "- Dashed = seqtk (comp) reference.",
            "- Indexed lines stay near the floor; full tracks composition cost as size or N grows.",
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

    manifest = load_manifest(results_dir)
    if not manifest:
        raise SystemExit(f"No run_*.json under {results_dir}")
    incomplete = is_incomplete(manifest)
    if incomplete and not args.allow_incomplete:
        raise SystemExit(
            "Refusing to overwrite REPORT.md from incomplete run "
            "(missing sections or skip_* flags). Pass --allow-incomplete for a draft."
        )

    full_df = load_section(results_dir, manifest, "perf_full", ir)
    mode_df = load_section(results_dir, manifest, "perf_mode", ir)
    size_df = enrich_size(load_section(results_dir, manifest, "scale_size", ir))
    seq_df = enrich_seqs(load_section(results_dir, manifest, "scale_seqs_fixed", ir))

    if full_df is None and mode_df is None and size_df is None and seq_df is None:
        raise SystemExit("No section data found")

    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    prune_stale_figures(figures_dir)

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
            md_overview(manifest),
            "## What each tool reports",
            md_field_matrix(),
            "## Run Provenance",
            md_run_provenance(manifest, ir),
        ]
    )

    # --- Full peers ---
    if full_df is not None and not full_df.empty:
        full_tools = tools_in_run(full_df, FULL_TOOLS)
        fig_note = (
            "Bar labels = peer / z-fasta. Hatched = seqtk (comp) reference. Error bars = 1σ."
        )
        fig_dataset_grouped_bars(
            full_df,
            figures_dir / "perf_full_wall.png",
            full_tools,
            "mean",
            "Wall time (s, log)",
            "Full stats: wall time",
            fig_note,
            ir,
            log_y=True,
            display_map=FULL_SECTION_DISPLAY,
        )
        fig_dataset_grouped_bars(
            full_df,
            figures_dir / "perf_full_rss.png",
            full_tools,
            "peak_rss_mb",
            "Peak RSS (MB)",
            "Full stats: peak RSS",
            "Bar labels = peer / z-fasta. Error bars = 1σ.",
            ir,
            log_y=False,
            value_floor=1e-3,
            display_map=FULL_SECTION_DISPLAY,
        )
        fig_dataset_grouped_bars(
            full_df,
            figures_dir / "perf_full_faults.png",
            full_tools,
            "minor_faults",
            "Minor page faults (log)",
            "Full stats: minor page faults",
            "Bar labels = peer / z-fasta. Error bars = 1σ.",
            ir,
            log_y=True,
            value_floor=1.0,
            display_map=FULL_SECTION_DISPLAY,
        )

        report_lines.append("## Performance: Full stats")
        report_lines.append(
            "Peers re-parse the FASTA. z-fasta full loads `.zfi`, emits length metrics from the "
            "index, then mmap-scans sequence bytes for composition. Baseline for ratios: "
            "**z-fasta full**."
        )
        report_lines.append(
            md_metric_block(
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
                fig_rel="results/figures/perf_full_wall.png",
                fig_caption=(
                    "Table above as grouped bars (log y). Three panels: Genome, Proteome, "
                    "Transcriptome."
                ),
                reading=[
                    "Panels are Genome / Proteome / Transcriptome (`sharey=False`).",
                    "Bars = zebrac mean wall; error bars = 1σ across measured samples.",
                    "Gold = z-fasta full (`.zfi`). Hatched purple = seqtk `comp` "
                    "(composition-only reference; omitted on Proteome).",
                    "`1×` on z-fasta; other labels = peer wall / z-fasta wall. "
                    "Details in the ratio table.",
                ],
                baseline=BASELINE,
                ratio_summary="Time × = peer wall / z-fasta full. Same ratios as bar labels.",
                zf_label="z-fasta",
                comp_label="Peer",
                ratio_label="Time ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                fmt_comp=lambda r: f"{r.comp_v:.4f}s",
                display_map=FULL_SECTION_DISPLAY,
            )
        )
        report_lines.append(
            md_metric_block(
                heading="Memory Usage",
                intro="Peak RSS from zebrac samples for the same full-stats commands.",
                work=full_df,
                tools=full_tools,
                value_col="peak_rss_mb",
                fmt=fmt_rss,
                ir=ir,
                nums=nums,
                fig_rel="results/figures/perf_full_rss.png",
                fig_caption="Peak RSS (linear y). Same tools and datasets as wall time.",
                reading=[
                    "Same panel layout as wall time.",
                    "Bar labels = peer RSS / z-fasta RSS.",
                    "Hatched = seqtk reference when present.",
                ],
                baseline=BASELINE,
                ratio_summary="RSS × = peer peak RSS / z-fasta full.",
                zf_label="z-fasta",
                comp_label="Peer",
                ratio_label="RSS ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.1f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.1f} MB",
                display_map=FULL_SECTION_DISPLAY,
            )
        )
        report_lines.append(
            md_metric_block(
                heading="Page Faults",
                intro="Minor page faults for the same full-stats commands.",
                work=full_df,
                tools=full_tools,
                value_col="minor_faults",
                fmt=fmt_faults,
                ir=ir,
                nums=nums,
                fig_rel="results/figures/perf_full_faults.png",
                fig_caption="Minor page faults (log y).",
                reading=[
                    "Same panel layout as wall time.",
                    "Bar labels = peer faults / z-fasta faults.",
                ],
                baseline=BASELINE,
                ratio_summary="Faults × = peer minor faults / z-fasta full.",
                zf_label="z-fasta",
                comp_label="Peer",
                ratio_label="Faults ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.0f}",
                fmt_comp=lambda r: f"{r.comp_v:.0f}",
                display_map=FULL_SECTION_DISPLAY,
            )
        )

    # --- Mode 2x2 ---
    if mode_df is not None and not mode_df.empty:
        mode_tools = tools_in_run(mode_df, MODE_TOOLS)
        report_lines.append("## z-fasta Mode Comparison")
        report_lines.append(
            "Four lanes: **full** and **indexed** (`--index-only`), each with `.zfi` and with "
            "`.fai` (`.zfi` stashed outside the timed window). Indexed skips the composition "
            "mmap scan. `.fai` is first-class for users who skip `.zfi` (messy side tables still "
            "need `.zfi`)."
        )
        tax = md_composition_tax(mode_df, nums)
        if tax:
            report_lines.append(tax)

        fig_dataset_grouped_bars(
            mode_df,
            figures_dir / "perf_mode_wall.png",
            mode_tools,
            "mean",
            "Wall time (s, log)",
            "Mode comparison: wall time",
            "Bar labels = other / z-fasta full. Error bars = 1σ.",
            ir,
            log_y=True,
            annotate=True,
        )
        fig_dataset_grouped_bars(
            mode_df,
            figures_dir / "perf_mode_rss.png",
            mode_tools,
            "peak_rss_mb",
            "Peak RSS (MB)",
            "Mode comparison: peak RSS",
            "Bar labels = other / z-fasta full. Error bars = 1σ.",
            ir,
            log_y=False,
            value_floor=1e-3,
        )
        fig_dataset_grouped_bars(
            mode_df,
            figures_dir / "perf_mode_faults.png",
            mode_tools,
            "minor_faults",
            "Minor page faults (log)",
            "Mode comparison: minor page faults",
            "Bar labels = other / z-fasta full. Error bars = 1σ.",
            ir,
            log_y=True,
            value_floor=1.0,
        )

        report_lines.append(
            md_metric_block(
                heading="Mode wall time",
                intro="Wall time for the four z-fasta lanes on REAL datasets.",
                work=mode_df,
                tools=mode_tools,
                value_col="mean",
                fmt=fmt_wall,
                ir=ir,
                nums=nums,
                fig_rel="results/figures/perf_mode_wall.png",
                fig_caption="Mode wall times (log y). Baseline gold = full (`.zfi`).",
                reading=[
                    "Four solid gold shades: full, full (fai), indexed (zfi), indexed (fai).",
                    "`1×` on full; other labels = lane / full.",
                    "Indexed lanes should be near-instant (lengths only). FAI tax shows most "
                    "on many-record files.",
                ],
                baseline=BASELINE,
                ratio_summary="Time × = other lane / z-fasta full.",
                zf_label="z-fasta full",
                comp_label="Other",
                ratio_label="Time ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                fmt_comp=lambda r: f"{r.comp_v:.4f}s",
            )
        )
        report_lines.append(
            md_metric_block(
                heading="Mode memory",
                intro="Peak RSS for the four z-fasta lanes.",
                work=mode_df,
                tools=mode_tools,
                value_col="peak_rss_mb",
                fmt=fmt_rss,
                ir=ir,
                nums=nums,
                fig_rel="results/figures/perf_mode_rss.png",
                fig_caption="Mode peak RSS (linear y).",
                reading=[
                    "Same four lanes as mode wall.",
                    "Bar labels = other RSS / full RSS.",
                ],
                baseline=BASELINE,
                ratio_summary="RSS × = other lane / z-fasta full.",
                zf_label="z-fasta full",
                comp_label="Other",
                ratio_label="RSS ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.1f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.1f} MB",
            )
        )
        report_lines.append(
            md_metric_block(
                heading="Mode page faults",
                intro="Minor page faults for the four z-fasta lanes.",
                work=mode_df,
                tools=mode_tools,
                value_col="minor_faults",
                fmt=fmt_faults,
                ir=ir,
                nums=nums,
                fig_rel="results/figures/perf_mode_faults.png",
                fig_caption="Mode minor page faults (log y).",
                reading=[
                    "Same four lanes as mode wall.",
                    "Bar labels = other faults / full faults.",
                ],
                baseline=BASELINE,
                ratio_summary="Faults × = other lane / z-fasta full.",
                zf_label="z-fasta full",
                comp_label="Other",
                ratio_label="Faults ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.0f}",
                fmt_comp=lambda r: f"{r.comp_v:.0f}",
            )
        )

    if size_df is not None and not size_df.empty:
        report_lines.append(
            md_scaling_section(
                size_df,
                param_col="size_mb",
                param_order=SIZE_MB_ORDER,
                xlabel="File size",
                section_title="Scaling: file size",
                fig_name="scaling_size.png",
                fig_title="Stats scaling by file size",
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
                fig_name="scaling_seqs_fixed.png",
                fig_title="Stats scaling by sequence count (1024 bp/seq)",
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
