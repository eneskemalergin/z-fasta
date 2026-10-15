#!/usr/bin/env python3

from representation_common import (
    DATASETS, TASKS, Theme, draw_footer, draw_header, format_multiple,
    format_rss, format_time, main, plot_color, tool_legend,
)

import matplotlib.pyplot as plt
import pandas as pd



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
            for metric in ("time", "rss"):
                ranks = complete[f"{metric}_rank"]
                complete[f"{metric}_position"] = ranks + 0.8 * complete.groupby(ranks).cumcount()
            complete["time_multiple"] = complete["wall"] / complete["wall"].min()
            complete["rss_multiple"] = complete["rss"] / complete["rss"].min()
            references = references.sort_values("wall").reset_index(drop=True)
            references["time_rank"] = len(complete) + 0.8 + references.index * 0.8
            references = references.sort_values("rss").reset_index(drop=True)
            references["rss_rank"] = len(complete) + 0.8 + references.index * 0.8
            references["time_position"] = references["time_rank"]
            references["rss_position"] = references["rss_rank"]
            subset = pd.concat([complete, references], ignore_index=True)
            n = max(float(subset["time_position"].max()), float(subset["rss_position"].max()))

            for record in subset.itertuples():
                tool = str(record.tool)
                is_zfasta = tool.startswith("z-fasta")
                secondary = record.comparison != "complete"
                color = plot_color(task, tool, theme)
                ax.plot(
                    [0, 1], [record.time_position, record.rss_position],
                    color=color, linewidth=3.2 if is_zfasta else 1.45,
                    alpha=1 if is_zfasta else 0.78,
                    linestyle="--" if secondary else "-",
                )
                ax.scatter(
                    [0, 1], [record.time_position, record.rss_position],
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
                    -0.035, record.time_position, f"{record.display}  {time_value}",
                    ha="right", va="center", fontsize=6.0,
                    color=text_color, fontweight=weight,
                )
                ax.text(
                    1.055, record.rss_position, f"{rss_value}  {record.display}",
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

    draw_header(fig, theme, meta)
    tool_legend(
        fig, theme, x=0.565, y=0.800, loc="center",
        ncol=5, fontsize=6.7, versions=versions,
    )

    draw_footer(fig, theme, meta, ("RANKING RIBBONS", "Rank 1 = best complete lane for that metric", "Labels = absolute value | x best complete", "Open + dashed lanes = reference or partial work, not ranked"))
    return fig


if __name__ == "__main__":
    main(render, "03-ranking-ribbons")
