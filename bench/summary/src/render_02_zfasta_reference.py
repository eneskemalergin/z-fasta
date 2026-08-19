#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np

from representation_common import (
    BASELINE_TOOL, DATASETS, DATASET_DETAILS, DATASET_IDENTITIES, MARKERS, ROOT,
    TASKS, TOOLS, Theme,
    Y_LIMITS, clean_axis, configure, extract, plot_color, provenance,
    row_geometry, save, tool_legend, tool_versions,
)


def render(frame, theme: Theme, meta, versions):
    fig = plt.figure(figsize=(16, 6.75))
    grid = fig.add_gridspec(
        1, 3, width_ratios=(0.92, 4.7, 4.7), left=0.025, right=0.975,
        top=0.735, bottom=0.21, wspace=0.012,
    )
    labels = fig.add_subplot(grid[0, 0])
    wall = fig.add_subplot(grid[0, 1])
    rss = fig.add_subplot(grid[0, 2], sharey=wall)
    y_map, boundaries = row_geometry()

    labels.set_xlim(0, 1)
    labels.set_ylim(*Y_LIMITS)
    labels.axis("off")
    for task in TASKS:
        center = np.mean([y_map[(task, dataset)] for dataset in DATASETS])
        labels.text(
            0.08, center, task, rotation=90, ha="center", va="center",
            fontsize=11.5, fontweight="bold", color=theme.section,
        )
        labels.plot(
            [0.19, 0.19], [center - 1.15, center + 1.15],
            color=theme.grid, linewidth=0.8,
        )
        for dataset in DATASETS:
            row_y = y_map[(task, dataset)]
            labels.text(0.97, row_y, dataset, ha="right", va="center", fontsize=9.2)

    ranges = {"wall": (0.8, 128), "rss": (0.06, 512)}
    ticks = {
        "wall": [1, 2, 4, 8, 16, 32, 64, 128],
        "rss": [0.1, 0.3, 1, 3, 10, 30, 100, 300],
    }
    for ax, metric, title in ((wall, "wall", "WALL TIME"), (rss, "rss", "PEAK RSS")):
        for task in TASKS:
            for dataset in DATASETS:
                center = y_map[(task, dataset)]
                subset = frame[(frame["task"] == task) & (frame["dataset"] == dataset)].copy()
                baseline_tool = BASELINE_TOOL[task]
                baseline = float(subset[subset["tool"].astype(str) == baseline_tool].iloc[0][metric])
                peers = subset[subset["tool"].astype(str) != baseline_tool]
                ordered = [tool for tool in TOOLS if tool in set(peers["tool"].astype(str))]
                offsets = np.linspace(0.28, -0.28, max(len(ordered), 1))
                for offset, tool in zip(offsets, ordered, strict=True):
                    row = peers[peers["tool"].astype(str) == tool].iloc[0]
                    ratio = float(row[metric]) / baseline
                    plotted = min(max(ratio, ranges[metric][0]), ranges[metric][1])
                    secondary = row["comparison"] != "complete"
                    color = plot_color(task, tool, theme)
                    y = center + offset
                    ax.plot([1, plotted], [y, y], color=color, linewidth=1.15, alpha=0.5, linestyle="--" if secondary else "-", zorder=2)
                    ax.scatter(plotted, y, s=52 if tool.startswith("z-fasta") else 34, marker=MARKERS[tool], facecolor=theme.background if secondary else color, edgecolor=color, linewidth=1.15, zorder=4)
                    if ratio > ranges[metric][1]:
                        ax.annotate(f"{ratio:.0f}x, off scale", (plotted, y), xytext=(-4, 0), textcoords="offset points", ha="right", va="center", fontsize=6.5, fontweight="bold", color=color)
                baseline_color = plot_color(task, baseline_tool, theme)
                ax.scatter(1, center, s=66, marker=MARKERS[baseline_tool], facecolor=baseline_color, edgecolor=theme.ink, linewidth=1, zorder=5)
        ax.axvline(1, color="#F7A41D", linewidth=1.5, linestyle=(0, (3, 3)), zorder=1)
        ax.set_xscale("log")
        ax.set_xlim(*ranges[metric])
        ax.xaxis.set_major_locator(FixedLocator(ticks[metric]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}x"))
        ax.set_ylim(*Y_LIMITS)
        ax.set_yticks([])
        ax.set_title("")
        ax.set_xlabel("peer / z-fasta", fontsize=9)
        clean_axis(ax, theme)
        for row_y in y_map.values():
            ax.axhline(row_y, color=theme.grid, linewidth=0.45, alpha=0.42, zorder=0)
    wall_box = wall.get_position()
    rss_box = rss.get_position()
    wall_one_x = fig.transFigure.inverted().transform(wall.transData.transform((1, 0)))[0]
    fig.text(
        0.145, 0.956, "Speed and memory across real FASTA shapes",
        ha="left", va="top", fontsize=22, fontweight="bold", color=theme.ink,
    )
    fig.text(
        0.145, 0.902,
        "FASTA index construction | 1 kbp positional extraction | complete statistics output",
        ha="left", va="top", fontsize=10.2, color=theme.muted,
    )
    fig.text(
        0.985, 0.952, f"z-fasta v{meta.zfasta_version}",
        ha="right", va="top", fontsize=10.5, fontweight="bold", color=theme.ink,
    )
    fig.text(
        0.985, 0.925, "github.com/eneskemalergin/z-fasta",
        ha="right", va="top", fontsize=7.2, fontweight="bold", color=theme.section,
    )
    fig.text(
        0.985, 0.894, "REAL FASTA BENCHMARK SUMMARY",
        ha="right", va="top", fontsize=7.4, fontweight="bold", color=theme.muted,
    )
    fig.add_artist(plt.Line2D(
        [0.145, 0.985], [0.852, 0.852], transform=fig.transFigure,
        color=theme.grid, linewidth=0.8,
    ))
    fig.text(wall_one_x, 0.805, "WALL TIME", ha="left", va="center", fontsize=13, fontweight="bold", color=theme.section)
    fig.text(wall_one_x, 0.768, "Right of 1x means z-fasta is faster", ha="left", va="center", fontsize=8.0, color=theme.muted)
    fig.text(rss_box.x1, 0.805, "PEAK RSS", ha="right", va="center", fontsize=13, fontweight="bold", color=theme.section)
    fig.text(rss_box.x1, 0.768, "Right of 1x means z-fasta uses less memory", ha="right", va="center", fontsize=8.0, color=theme.muted)
    tool_legend(fig, theme, x=0.565, y=0.803, loc="center", ncol=5, fontsize=6.7, versions=versions)
    footer_top = 0.135
    for x in (0.415, 0.705):
        fig.add_artist(plt.Line2D(
            [x, x], [0.027, footer_top], transform=fig.transFigure,
            color=theme.grid, linewidth=0.8,
        ))
    fig.text(0.025, 0.116, "HUMAN DATASETS", ha="left", va="center", fontsize=7.0, fontweight="bold", color=theme.section)
    for y, dataset in zip((0.087, 0.060, 0.033), DATASETS, strict=True):
        details = DATASET_DETAILS[dataset].replace("  ·  ", " | ")
        fig.text(
            0.025, y, f"{dataset}  {DATASET_IDENTITIES[dataset]} | {details}",
            ha="left", va="center", fontsize=6.6, color=theme.ink,
        )
    fig.text(0.560, 0.116, "RATIO REFERENCE", ha="center", va="center", fontsize=7.0, fontweight="bold", color=theme.section)
    fig.text(0.560, 0.087, "INDEX  z-fasta .fai = 1x", ha="center", va="center", fontsize=6.8, fontweight="bold", color=theme.ink)
    fig.text(0.560, 0.060, "GET + STATS  z-fasta .zfi = 1x", ha="center", va="center", fontsize=6.8, fontweight="bold", color=theme.ink)
    fig.text(0.560, 0.033, "Open + dashed marks = reference or partial work", ha="center", va="center", fontsize=6.6, color=theme.muted)
    fig.text(0.975, 0.116, "MEASUREMENT", ha="right", va="center", fontsize=7.0, fontweight="bold", color=theme.section)
    fig.text(0.975, 0.087, f"zebrac {meta.zebrac_version} | {meta.samples} samples | {meta.warmups} warmups", ha="right", va="center", fontsize=6.8, fontweight="bold", color=theme.ink)
    fig.text(0.975, 0.060, f"{meta.duration_ms} ms budget | Zig {meta.zig_version}", ha="right", va="center", fontsize=6.8, color=theme.ink)
    fig.text(0.975, 0.033, f"Warm-cache benchmark | {meta.benchmark_date}", ha="right", va="center", fontsize=6.6, color=theme.muted)
    return fig


if __name__ == "__main__":
    data = extract()
    run_meta = provenance()
    peer_versions = tool_versions(run_meta.zfasta_version)
    icon = ROOT / "assets/icon.svg"
    for theme_name in ("light", "dark"):
        active_theme = configure(theme_name)
        save(
            render(data, active_theme, run_meta, peer_versions),
            f"02-zfasta-reference-{theme_name}",
            embedded_svg=icon,
            embedded_box=(51.84, 16.5, 88, 67.1),
            render_png=False,
        )
