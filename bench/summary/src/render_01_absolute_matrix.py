#!/usr/bin/env python3

from representation_common import (
    DATASETS, TASKS, Theme, Y_LIMITS, clean_axis, draw_footer, draw_header,
    draw_points, draw_row_labels, main, row_geometry, tool_legend,
)

import matplotlib.pyplot as plt
import numpy as np



def render(frame, theme: Theme, meta, versions):
    fig = plt.figure(figsize=(16, 6.75))
    grid = fig.add_gridspec(
        1, 3, width_ratios=(0.92, 4.7, 4.7), left=0.025, right=0.975,
        top=0.735, bottom=0.21, wspace=0.012,
    )
    labels = fig.add_subplot(grid[0, 0])
    wall = fig.add_subplot(grid[0, 1])
    rss = fig.add_subplot(grid[0, 2], sharey=wall)
    y_map, _ = row_geometry()

    draw_row_labels(labels, y_map, theme)
    for ax, metric, unit in (
        (wall, "wall", "seconds"),
        (rss, "rss", "MB"),
    ):
        for task in TASKS:
            for dataset in DATASETS:
                subset = frame[(frame["task"] == task) & (frame["dataset"] == dataset)]
                draw_points(ax, subset, metric, y_map[(task, dataset)], theme)
        ax.set_xscale("log")
        ax.set_ylim(*Y_LIMITS)
        ax.set_yticks([])
        ax.set_title("")
        ax.set_xlabel(unit, fontsize=9)
        clean_axis(ax, theme)
        for row_y in y_map.values():
            ax.axhline(row_y, color=theme.grid, linewidth=0.45, alpha=0.42, zorder=0)

    wall_box = wall.get_position()
    rss_box = rss.get_position()
    draw_header(fig, theme, meta)
    fig.text(wall_box.x0, 0.805, "WALL TIME", ha="left", va="center", fontsize=13, fontweight="bold", color=theme.section)
    fig.text(wall_box.x0, 0.768, "Lower values mean less elapsed time", ha="left", va="center", fontsize=8.0, color=theme.muted)
    fig.text(rss_box.x1, 0.805, "PEAK RSS", ha="right", va="center", fontsize=13, fontweight="bold", color=theme.section)
    fig.text(rss_box.x1, 0.768, "Lower values mean less process memory", ha="right", va="center", fontsize=8.0, color=theme.muted)
    tool_legend(
        fig, theme, x=0.565, y=0.803, loc="center",
        ncol=5, fontsize=6.7, versions=versions,
    )

    draw_footer(fig, theme, meta, ("ABSOLUTE VALUES", "Points = zebrac mean", "Horizontal whiskers = 1 standard deviation", "Open marks = reference or partial work"))
    return fig


if __name__ == "__main__":
    main(render, "01-absolute-matrix")
