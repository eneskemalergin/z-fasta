#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "bench/index/results/perf_20260818_070220"
GET = ROOT / "bench/get/results/perf_pos_20260818_070223"
STATS = ROOT / "bench/stats/results/perf_full_20260818_073417"

DATASETS = ["Genome", "Transcriptome", "Proteome"]
TASKS = ["INDEX", "GET", "STATS"]
TOOLS = [
    "z-fasta .zfi", "z-fasta .fai", "noodles", "rust-bio", "samtools",
    "seqkit", "seqtk", "fastahack", "pyfaidx",
]

MARKERS = {
    "z-fasta .zfi": "D", "z-fasta .fai": "s", "noodles": "o",
    "rust-bio": "^", "samtools": "v", "seqkit": "P", "seqtk": "*",
    "fastahack": "X", "pyfaidx": "h",
}

TASK_COLORS_LIGHT = {
    "INDEX": {
        "z-fasta .fai": "#F7A41D", "noodles": "#C45C26",
        "rust-bio": "#8B3A2A", "samtools": "#555555",
        "seqkit": "#009485", "fastahack": "#F34B7D", "pyfaidx": "#3572A5",
    },
    "GET": {
        "z-fasta .zfi": "#F7A41D", "z-fasta .fai": "#FFC45C",
        "noodles": "#C45C26", "rust-bio": "#8B3A2A",
        "samtools": "#555555", "seqtk": "#6A1B9A",
    },
    "STATS": {
        "z-fasta .zfi": "#F7A41D", "z-fasta .fai": "#FFB74D",
        "noodles": "#C45C26", "rust-bio": "#8B3A2A",
        "seqkit": "#1565C0", "seqtk": "#6A1B9A",
    },
}

TASK_COLORS_DARK = {
    "INDEX": {
        "z-fasta .fai": "#F7A41D", "noodles": "#E07A3F",
        "rust-bio": "#D46A54", "samtools": "#B8C0C5",
        "seqkit": "#25BDB0", "fastahack": "#FF6B9C", "pyfaidx": "#6CA7D9",
    },
    "GET": {
        "z-fasta .zfi": "#F7A41D", "z-fasta .fai": "#FFD078",
        "noodles": "#E07A3F", "rust-bio": "#D46A54",
        "samtools": "#B8C0C5", "seqtk": "#C58AF9",
    },
    "STATS": {
        "z-fasta .zfi": "#F7A41D", "z-fasta .fai": "#FFD078",
        "noodles": "#E07A3F", "rust-bio": "#D46A54",
        "seqkit": "#58A6FF", "seqtk": "#C58AF9",
    },
}


@dataclass(frozen=True)
class Theme:
    name: str
    background: str
    ink: str
    muted: str
    grid: str
    section: str
    colors: dict[str, dict[str, str]]


THEMES = {
    "light": Theme(
        name="light", background="#FFFFFF", ink="#1F2328", muted="#59636E",
        grid="#D0D7DE", section="#B66C11", colors=TASK_COLORS_LIGHT,
    ),
    "dark": Theme(
        name="dark", background="#0D1117", ink="#F0F3F6", muted="#9DA7B0",
        grid="#30363D", section="#F7A41D", colors=TASK_COLORS_DARK,
    ),
}

DATASET_DETAILS = {
    "Genome": "3.15 GB  ·  194 sequences",
    "Transcriptome": "481 MB  ·  254K sequences",
    "Proteome": "13.7 MB  ·  20.7K sequences",
}

DATASET_IDENTITIES = {
    "Genome": "GRCh38 primary (Ensembl 113)",
    "Transcriptome": "GENCODE v46",
    "Proteome": "UniProt UP000005640",
}

TASK_DETAILS = {
    "INDEX": "FAI build",
    "GET": "1 kbp positional",
    "STATS": "complete stats",
}

BASELINE_TOOL = {
    "INDEX": "z-fasta .fai",
    "GET": "z-fasta .zfi",
    "STATS": "z-fasta .zfi",
}

Y_LIMITS = (-1.72, 8.72)


@dataclass(frozen=True)
class Provenance:
    zfasta_version: str
    zig_version: str
    zebrac_version: str
    samples: int
    warmups: int
    duration_ms: int
    index_run: str
    get_run: str
    stats_run: str
    benchmark_date: str


def provenance() -> Provenance:
    zon = (ROOT / "build.zig.zon").read_text()
    version_match = re.search(r'\.version\s*=\s*"([^"]+)"', zon)
    zig_match = re.search(r'\.minimum_zig_version\s*=\s*"([^"]+)"', zon)
    if version_match is None or zig_match is None:
        raise ValueError("build.zig.zon is missing version metadata")

    manifests = [
        load_json(INDEX / "Genome.json"),
        load_json(GET / "Genome_1kbp_mid.json"),
        load_json(STATS / "Genome__z-fasta-zfi.json"),
    ]
    configs = [manifest["config"] for manifest in manifests]
    config_fields = ("min_samples", "max_samples", "warmup", "duration_ms")
    if any(
        config[field] != configs[0][field]
        for config in configs[1:]
        for field in config_fields
    ):
        raise ValueError("selected benchmark runs do not share one sampling configuration")
    zebrac_versions = {manifest["zebrac_version"] for manifest in manifests}
    if len(zebrac_versions) != 1:
        raise ValueError("selected benchmark runs do not share one zebrac version")

    config = configs[0]
    run_ids = {
        "index": INDEX.name.removeprefix("perf_"),
        "get": GET.name.removeprefix("perf_pos_"),
        "stats": STATS.name.removeprefix("perf_full_"),
    }
    run_dates = {
        datetime.strptime(run_id[:8], "%Y%m%d").strftime("%Y-%m-%d")
        for run_id in run_ids.values()
    }
    if len(run_dates) != 1:
        raise ValueError("selected benchmark runs do not share one benchmark date")
    return Provenance(
        zfasta_version=version_match.group(1),
        zig_version=zig_match.group(1),
        zebrac_version=zebrac_versions.pop(),
        samples=config["max_samples"],
        warmups=config["warmup"],
        duration_ms=config["duration_ms"],
        index_run=run_ids["index"],
        get_run=run_ids["get"],
        stats_run=run_ids["stats"],
        benchmark_date=run_dates.pop(),
    )


def tool_versions(zfasta_version: str) -> dict[str, str]:
    source = (ROOT / "tools/versions.sh").read_text()
    assignments = dict(re.findall(r'^([A-Z0-9_]+)="([^"]+)"$', source, re.MULTILINE))
    required = {
        "NOODLES_VERSION", "RUSTBIO_VERSION", "SAMTOOLS_VERSION",
        "SEQKIT_VERSION", "SEQTK_DISPLAY_VERSION", "FASTAHACK_VERSION",
        "PYFAIDX_VERSION",
    }
    missing = required - assignments.keys()
    if missing:
        raise ValueError(f"tools/versions.sh is missing: {', '.join(sorted(missing))}")
    return {
        "z-fasta .zfi": zfasta_version,
        "z-fasta .fai": zfasta_version,
        "noodles": assignments["NOODLES_VERSION"],
        "rust-bio": assignments["RUSTBIO_VERSION"],
        "samtools": assignments["SAMTOOLS_VERSION"],
        "seqkit": assignments["SEQKIT_VERSION"],
        "seqtk": assignments["SEQTK_DISPLAY_VERSION"],
        "fastahack": assignments["FASTAHACK_VERSION"],
        "pyfaidx": assignments["PYFAIDX_VERSION"],
    }


def configure(name: str) -> Theme:
    theme = THEMES[name]
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "svg.fonttype": "none",
        "figure.facecolor": theme.background,
        "axes.facecolor": theme.background,
        "savefig.facecolor": theme.background,
        "text.color": theme.ink,
        "axes.labelcolor": theme.muted,
        "xtick.color": theme.muted,
        "ytick.color": theme.ink,
        "axes.edgecolor": theme.grid,
    })
    return theme


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def classify_command(command: str) -> str | None:
    if "/zig-out/bin/z-fasta" in command:
        return "z-fasta"
    for token, label in (
        ("/tools/bin/noodles", "noodles"),
        ("/tools/bin/rustbio", "rust-bio"),
        ("/tools/bin/samtools", "samtools"),
        ("/tools/bin/seqkit", "seqkit"),
        ("/tools/bin/seqtk", "seqtk"),
        ("/tools/bin/fastahack", "fastahack"),
        ("/tools/bin/faidx", "pyfaidx"),
    ):
        if token in command:
            return label
    return None


def observation(
    task: str,
    dataset: str,
    tool: str,
    result: dict,
    comparison: str = "complete",
    display: str | None = None,
) -> dict:
    return {
        "task": task,
        "dataset": dataset,
        "tool": tool,
        "display": display or tool,
        "comparison": comparison,
        "wall": result["wall_time"]["mean"] / 1e9,
        "wall_std": result["wall_time"]["std_dev"] / 1e9,
        "rss": result["peak_rss"]["mean"] / 1e6,
        "rss_std": result["peak_rss"]["std_dev"] / 1e6,
    }


def extract() -> pd.DataFrame:
    rows: list[dict] = []
    for dataset in DATASETS:
        for result in load_json(INDEX / f"{dataset}.json")["results"]:
            command = result["command"]
            tool = classify_command(command)
            if tool == "z-fasta":
                if "index\\ --emit-fai" not in command or "--no-dedup" in command:
                    continue
                tool = "z-fasta .fai"
            if tool is not None:
                display = "SeqKit faidx" if tool == "seqkit" else tool
                rows.append(observation("INDEX", dataset, tool, result, display=display))

    for dataset in DATASETS:
        for result in load_json(GET / f"{dataset}_1kbp_mid.json")["results"]:
            command = result["command"]
            tool = classify_command(command)
            if tool == "z-fasta":
                tool = "z-fasta .zfi" if "/zfi/" in command else "z-fasta .fai"
            if tool in {"z-fasta .zfi", "z-fasta .fai", "noodles", "rust-bio", "samtools", "seqtk"}:
                comparison = "reference" if tool == "seqtk" else "complete"
                display = "seqtk scan ref" if tool == "seqtk" else tool
                rows.append(observation("GET", dataset, tool, result, comparison, display))

    for path in sorted(STATS.glob("*.json")):
        dataset, raw_tool = path.stem.split("__", maxsplit=1)
        tool = {"z-fasta-zfi": "z-fasta .zfi", "z-fasta-fai": "z-fasta .fai", "rustbio": "rust-bio"}.get(raw_tool, raw_tool)
        if tool not in {"z-fasta .zfi", "z-fasta .fai", "noodles", "rust-bio", "seqkit", "seqtk"}:
            continue
        comparison = "partial" if tool in {"seqkit", "seqtk"} else "complete"
        display = f"{tool} partial" if comparison == "partial" else tool
        rows.append(observation("STATS", dataset, tool, load_json(path)["results"][0], comparison, display))

    frame = pd.DataFrame(rows)
    frame["task"] = pd.Categorical(frame["task"], TASKS, ordered=True)
    frame["dataset"] = pd.Categorical(frame["dataset"], DATASETS, ordered=True)
    frame["tool"] = pd.Categorical(frame["tool"], TOOLS, ordered=True)
    return frame.sort_values(["task", "dataset", "tool"]).reset_index(drop=True)


def plot_color(task: str, tool: str, theme: Theme) -> str:
    return theme.colors[task][tool]


def row_geometry() -> tuple[dict[tuple[str, str], float], list[float]]:
    positions: dict[tuple[str, str], float] = {}
    boundaries: list[float] = []
    y = 8.0
    for task_i, task in enumerate(TASKS):
        for dataset in DATASETS:
            positions[(task, dataset)] = y
            y -= 1.0
        if task_i < len(TASKS) - 1:
            next_top = y - 0.55
            boundaries.append((positions[(task, DATASETS[-1])] + next_top) / 2)
            y = next_top
    return positions, boundaries


def clean_axis(ax: plt.Axes, theme: Theme, grid: bool = True) -> None:
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(theme.grid)
    ax.tick_params(axis="both", length=0, labelsize=8)
    if grid:
        ax.grid(axis="x", color=theme.grid, linewidth=0.7, linestyle=(0, (2, 4)))
        ax.set_axisbelow(True)


def draw_points(ax: plt.Axes, subset: pd.DataFrame, metric: str, y: float, theme: Theme, spread: float = 0.28) -> None:
    task = str(subset.iloc[0]["task"])
    ordered = [tool for tool in TOOLS if tool in set(subset["tool"].astype(str))]
    offsets = np.linspace(spread, -spread, max(len(ordered), 1))
    for offset, tool in zip(offsets, ordered, strict=True):
        row = subset[subset["tool"].astype(str) == tool].iloc[0]
        value = float(row[metric])
        error = float(row[f"{metric}_std"])
        is_zfasta = tool.startswith("z-fasta")
        secondary = row["comparison"] != "complete"
        color = plot_color(task, tool, theme)
        ax.errorbar(
            value, y + offset, xerr=error if error > 0 else None,
            marker=MARKERS[tool], markersize=7.5 if is_zfasta else 5.5,
            markerfacecolor=theme.background if secondary else color,
            markeredgecolor=theme.ink if is_zfasta else color, markeredgewidth=1.0,
            color=color, ecolor=mcolors.to_rgba(color, 0.55), elinewidth=1,
            capsize=1.5, linewidth=0, zorder=4 if is_zfasta else 3,
        )


def tool_legend(
    fig: plt.Figure,
    theme: Theme,
    y: float = 0.02,
    x: float = 0.5,
    loc: str = "lower center",
    ncol: int = 9,
    fontsize: float = 8,
    versions: dict[str, str] | None = None,
) -> None:
    specs = [
        ("z-fasta .zfi", "z-fasta .zfi", "GET", False),
        ("z-fasta .fai", "z-fasta .fai", "INDEX", False),
        ("noodles", "noodles", "INDEX", False),
        ("rust-bio", "rust-bio", "INDEX", False),
        ("samtools", "samtools", "INDEX", False),
        ("seqkit", "SeqKit faidx", "INDEX", False),
        ("seqkit", "SeqKit partial", "STATS", True),
        ("seqtk", "seqtk reference", "GET", True),
        ("fastahack", "fastahack", "INDEX", False),
        ("pyfaidx", "pyfaidx", "INDEX", False),
    ]
    handles = []
    for tool, label, task, secondary in specs:
        color = plot_color(task, tool, theme)
        display = f"{label} {versions[tool]}" if versions is not None else label
        handles.append(Line2D(
            [0], [0], marker=MARKERS[tool], color=color if secondary else "none",
            markerfacecolor=theme.background if secondary else color, markeredgecolor=color,
            linestyle="--" if secondary else "none", markeredgewidth=1,
            markersize=7, label=display,
        ))
    fig.legend(
        handles=handles, loc=loc, bbox_to_anchor=(x, y), ncol=ncol,
        frameon=False, fontsize=fontsize, handletextpad=0.35, columnspacing=1.0,
    )


def format_time(value: float) -> str:
    return f"{value * 1000:.3g} ms" if value < 0.1 else f"{value:.3g} s"


def format_rss(value: float) -> str:
    return f"{value:,.0f} MB" if value >= 100 else f"{value:.3g} MB"


def format_multiple(value: float) -> str:
    return f"{value:.1f}×" if value < 10 else f"{value:.0f}×"


def save(
    fig: plt.Figure,
    stem: str,
    embedded_svg: Path | None = None,
    embedded_box: tuple[float, float, float, float] | None = None,
    render_png: bool = True,
) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    svg_path = OUT / f"{stem}.svg"
    png_path = OUT / f"{stem}.png"
    fig.savefig(svg_path, metadata={"Creator": stem})
    if embedded_svg is None and render_png:
        fig.savefig(png_path, dpi=150)
    elif embedded_svg is not None:
        if embedded_box is None:
            raise ValueError("embedded_box is required with embedded_svg")
        x, y, width, height = embedded_box
        import cairosvg
        icon_source = embedded_svg.read_text()
        icon_source = re.sub(
            r'font-family="ui-monospace[^"]*"',
            'font-family="Liberation Mono"',
            icon_source,
        )
        outlined = cairosvg.svg2svg(bytestring=icon_source.encode("utf-8"))
        source = ET.fromstring(outlined)
        for element in source.iter():
            if "id" in element.attrib:
                element.attrib["id"] = f"zfasta-logo-{element.attrib['id']}"
            for attribute, value in tuple(element.attrib.items()):
                if attribute.endswith("href") and value.startswith("#"):
                    element.attrib[attribute] = f"#zfasta-logo-{value[1:]}"
        vector = "".join(ET.tostring(child, encoding="unicode") for child in source)
        image = (
            f'<svg id="zfasta-logo" x="{x:g}" y="{y:g}" '
            f'width="{width:g}" height="{height:g}" viewBox="0 0 160 122" '
            f'preserveAspectRatio="xMidYMid meet">{vector}</svg>\n'
        )
        document = svg_path.read_text()
        svg_path.write_text(document.replace("</svg>", image + "</svg>"))
        if render_png:
            cairosvg.svg2png(
                url=str(svg_path), write_to=str(png_path),
                output_width=int(fig.get_figwidth() * 150),
                output_height=int(fig.get_figheight() * 150),
            )
    plt.close(fig)
