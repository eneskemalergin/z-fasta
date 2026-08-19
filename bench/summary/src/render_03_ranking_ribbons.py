#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

from representation_common import (
    DATASETS, DATASET_DETAILS, DATASET_IDENTITIES, ROOT, TASKS, Theme,
    configure, extract, format_multiple, format_rss, format_time, plot_color,
    provenance, save, tool_legend, tool_versions,
)


def ascii_multiple(value: float) -> str:
    return format_multiple(value).replace("×", "x")


def render(frame, theme: Theme, meta, versions):
    fig, axes = plt.subplots(3, 3, figsize=(16, 8.6))
    fig.subplots_adjust(
        left=0.075, right=0.975, top=0.720, bottom=0.180,
        wspace=0.16, hspace=0.12,
    )

    for row_i, task in enumerate(TASKS):
        for col_i, dataset in enumerate(DATASETS):
            ax = axes[row_i, col_i]
            subset = frame[(frame["task"] == task) & (frame["dataset"] == dataset)].copy()
            complete = subset[subset["comparison"] == "complete"].copy()
            references = subset[subset["comparison"] != "complete"].copy()
            complete["time_rank"] = complete["wall"].rank(method="min")
            complete["rss_rank"] = complete["rss"].rank(method="min")
            complete["time_multiple"] = complete["wall"] / complete["wall"].min()
            complete["rss_multiple"] = complete["rss"] / complete["rss"].min()
            references = references.sort_values("wall").reset_index(drop=True)
            references["time_rank"] = len(complete) + 0.8 + references.index * 0.8
            references = references.sort_values("rss").reset_index(drop=True)
            references["rss_rank"] = len(complete) + 0.8 + references.index * 0.8
            subset = pd.concat([complete, references], ignore_index=True)
            n = max(float(subset["time_rank"].max()), float(subset["rss_rank"].max()))

            for record in subset.itertuples():
                tool = str(record.tool)
                is_zfasta = tool.startswith("z-fasta")
                secondary = record.comparison != "complete"
                color = plot_color(task, tool, theme)
                ax.plot(
                    [0, 1], [record.time_rank, record.rss_rank],
                    color=color, linewidth=3.2 if is_zfasta else 1.45,
                    alpha=1 if is_zfasta else 0.78,
                    linestyle="--" if secondary else "-",
                )
                ax.scatter(
                    [0, 1], [record.time_rank, record.rss_rank],
                    facecolor=theme.background if secondary else color,
                    edgecolor=color, s=36 if is_zfasta else 23,
                    linewidth=1.1, zorder=3,
                )
                if secondary:
                    time_value = f"{format_time(record.wall)} | ref"
                    rss_value = f"{format_rss(record.rss)} | ref"
                else:
                    time_value = f"{format_time(record.wall)} | {ascii_multiple(record.time_multiple)}"
                    rss_value = f"{format_rss(record.rss)} | {ascii_multiple(record.rss_multiple)}"
                weight = "bold" if is_zfasta else "normal"
                text_color = theme.ink if is_zfasta else theme.muted
                ax.text(
                    -0.035, record.time_rank, f"{record.display}  {time_value}",
                    ha="right", va="center", fontsize=6.0,
                    color=text_color, fontweight=weight,
                )
                ax.text(
                    1.055, record.rss_rank, f"{rss_value}  {record.display}",
                    ha="left", va="center", fontsize=6.0,
                    color=text_color, fontweight=weight,
                )

            ax.set_xlim(-0.68, 1.68)
            ax.set_ylim(n + 0.65, 0.35)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["TIME | x BEST", "RSS | x BEST"] if row_i == 2 else ["", ""])
            ax.set_yticks([])
            ax.spines[:].set_visible(False)
            ax.tick_params(length=0, labelsize=7.2, colors=theme.muted)

    for row_i, task in enumerate(TASKS):
        left_box = axes[row_i, 0].get_position()
        center = (left_box.y0 + left_box.y1) / 2
        fig.text(
            0.025, center, task, rotation=90, ha="center", va="center",
            fontsize=11.5, fontweight="bold", color=theme.section,
        )
        fig.add_artist(plt.Line2D(
            [0.042, 0.042], [left_box.y0 + 0.015, left_box.y1 - 0.015],
            transform=fig.transFigure, color=theme.grid, linewidth=0.8,
        ))

    for col_i, dataset in enumerate(DATASETS):
        box = axes[0, col_i].get_position()
        center = (box.x0 + box.x1) / 2
        fig.text(
            center, 0.744, dataset.upper(), ha="center", va="center",
            fontsize=9.5, fontweight="bold", color=theme.ink,
        )

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
    tool_legend(
        fig, theme, x=0.565, y=0.800, loc="center",
        ncol=5, fontsize=6.7, versions=versions,
    )

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
    fig.text(0.560, 0.116, "RANKING RIBBONS", ha="center", va="center", fontsize=7.0, fontweight="bold", color=theme.section)
    fig.text(0.560, 0.087, "Rank 1 = best complete lane for that metric", ha="center", va="center", fontsize=6.8, fontweight="bold", color=theme.ink)
    fig.text(0.560, 0.060, "Labels = absolute value | x best complete", ha="center", va="center", fontsize=6.8, fontweight="bold", color=theme.ink)
    fig.text(0.560, 0.033, "Open + dashed lanes = reference or partial work, not ranked", ha="center", va="center", fontsize=6.6, color=theme.muted)
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
            f"03-ranking-ribbons-{theme_name}",
            embedded_svg=icon,
            embedded_box=(51.84, 16.5, 88, 67.1),
            render_png=False,
        )
