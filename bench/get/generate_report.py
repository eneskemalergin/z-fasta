#!/usr/bin/env python3
"""
z-fasta GET Benchmark Report Generator

Reads raw hyperfine JSON + CSV data from bench/get/results/,
produces Markdown report + PNG figures using pandas + matplotlib.

Usage:
    .venv/bin/python bench/get/generate_report.py [results_dir]

Defaults to bench/get/results/ (latest timestamped files).
"""

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import tabulate  # noqa: F401

# ── Styling ────────────────────────────────────────────────────────
COLORS = {
    "z-fasta": "#2E7D32",
    "samtools": "#1565C0",
    "bedtools": "#6A1B9A",
    "seqkit": "#E65100",
    "fastahack": "#7B1FA2",
    "seqtk": "#00838F",
    "pyfaidx": "#F7A41D",
}
TOOL_ORDER = ["z-fasta", "samtools", "bedtools", "seqkit", "fastahack", "seqtk", "pyfaidx"]
RC_COLORS = {
    "z-fasta-forward": "#1B5E20",
    "z-fasta-rc": "#2E7D32",
    "z-fasta-sibling": "#66BB6A",
    "samtools": "#1565C0",
    "bedtools": "#8E24AA",
    "seqtk": "#00838F",
}
TOOL_DISPLAY = {
    "z-fasta": "z-fasta",
    "samtools": "samtools",
    "bedtools": "bedtools",
    "seqkit": "seqkit",
    "fastahack": "fastahack",
    "seqtk": "seqtk",
    "pyfaidx": "pyfaidx",
}

matplotlib.rcParams.update({
    "figure.facecolor": "#FBFAF7",
    "axes.facecolor": "#FBFAF7",
    "axes.edgecolor": "#D7D1C7",
    "axes.labelcolor": "#2B2B2B",
    "axes.titleweight": "bold",
    "xtick.color": "#3A3A3A",
    "ytick.color": "#3A3A3A",
    "grid.color": "#D9D3CA",
    "grid.alpha": 0.45,
    "legend.frameon": False,
    "font.size": 10,
})


def tool_sort_key(name):
    # Strip prefixes to find the base tool name
    base = name
    for t in TOOL_ORDER:
        if t in name:
            base = t
            break
    try:
        return TOOL_ORDER.index(base)
    except ValueError:
        return 99


def tool_color(name):
    for t in TOOL_ORDER:
        if t in name:
            return COLORS.get(t, "#888")
    return "#888"


def rc_tool_family(name: str) -> str:
    if name in {"forward", "bed-forward"} or name.endswith("-forward"):
        return "z-fasta-forward"
    if "samtools" in name:
        return "samtools"
    if "bedtools" in name:
        return "bedtools"
    if "seqtk" in name:
        return "seqtk"
    if name in {"zf-rc", "rc", "bed-zf-rc", "bed-rc", "bed-honor-strand-zf-rc", "bed-honor-strand-rc"} or name.endswith("-zf-rc"):
        return "z-fasta-rc"
    if name in {"reverse-only", "complement-only", "rc-annotate", "bed-reverse-only", "bed-complement-only"}:
        return "z-fasta-sibling"
    return "z-fasta-sibling"


def rc_color(name: str) -> str:
    return RC_COLORS.get(rc_tool_family(name), "#777777")


def _style_axes(ax):
    ax.grid(axis="y", alpha=0.45)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CFC7BC")
    ax.spines["bottom"].set_color("#CFC7BC")


def _strip_tool_prefix(name: str) -> str:
    """Extract base tool name from hyperfine label like '10mb_full_z-fasta'."""
    for base in TOOL_ORDER:
        if name.endswith(base):
            return base
        if base in name:
            return base
    return name


def _human_bench(name: str) -> str:
    """Convert raw benchmark stem to human-readable string.

    Examples:
        'size_10mb_100bp'    -> '10 MB / 100 bp'
        'size_100mb_10kbp'   -> '100 MB / 10 kbp'
        '10mb_full'          -> '10 MB (full)'
        '100bp'              -> '100 bp'
        'Genome_1kb'         -> 'Genome: 1 kbp region'
        'Genome_full'        -> 'Genome: full seq'
    """
    # single-region latency: size_Xmb_Ybp or size_Xmb_Ykbp
    m = re.match(r"size_(\d+)mb_(\d+)(k?)bp$", name)
    if m:
        mb, n, k = m.group(1), m.group(2), m.group(3)
        region = f"{n} kbp" if k else f"{n} bp"
        return f"{mb} MB / {region}"
    # full-sequence: Xmb_full
    m = re.match(r"(\d+)mb_full$", name)
    if m:
        return f"{m.group(1)} MB (full)"
    # region size only: Xbp or Xkbp
    m = re.match(r"(\d+)(k?)bp$", name)
    if m:
        n, k = m.group(1), m.group(2)
        return f"{n} kbp" if k else f"{n} bp"
    # real dataset: Name_region
    m = re.match(r"([A-Za-z]+)_(full|[\d]+(?:kb|mb))$", name)
    if m:
        ds, sz = m.group(1), m.group(2)
        label = {"full": "full seq", "1kb": "1 kbp region", "1mb": "1 Mbp region"}.get(sz, sz)
        return f"{ds}: {label}"
    m = re.match(r"(\d+)regions_(default|stranded)$", name)
    if m:
        count, mode = m.group(1), m.group(2)
        return f"{count} regions ({'default' if mode == 'default' else 'stranded'})"
    return name


def _parse_region_amount(label: str) -> int:
    match = re.match(r"(\d+)\s+(bp|kbp|Mbp)", label)
    if not match:
        return 0
    value = int(match.group(1))
    unit = match.group(2)
    if unit == "kbp":
        return value * 1_000
    if unit == "Mbp":
        return value * 1_000_000
    return value


def workload_sort_key(label: str):
    m = re.match(r"(\d+) MB / (\d+\s+(?:bp|kbp))$", label)
    if m:
        return (1, int(m.group(1)), _parse_region_amount(m.group(2)))

    m = re.match(r"(\d+) MB \(full\)$", label)
    if m:
        return (2, int(m.group(1)))

    m = re.match(r"(\d+)\s+(bp|kbp|Mbp)$", label)
    if m:
        return (3, _parse_region_amount(label))

    m = re.match(r"(\d+)\s+region[s]?$", label)
    if m:
        return (4, int(m.group(1)))

    m = re.match(r"(\d+) regions \((default|stranded)\)$", label)
    if m:
        return (5, int(m.group(1)), 0 if m.group(2) == "default" else 1)

    m = re.match(r"(Genome|Transcriptome|Proteome):\s+(.*)$", label)
    if m:
        dataset_order = {"Genome": 0, "Transcriptome": 1, "Proteome": 2}
        detail = m.group(2)
        detail_order = {"1 kbp region": 0, "1 Mbp region": 1, "full seq": 2}
        return (6, dataset_order.get(m.group(1), 99), detail_order.get(detail, 99), detail)

    return (99, label)


# ══════════════════════════════════════════════════════════════════════
#  Data Loading
# ══════════════════════════════════════════════════════════════════════


def load_hyperfine_json(path: Path) -> pd.DataFrame:
    with open(path) as f:
        data = json.load(f)
    rows = []
    for r in data.get("results", []):
        rows.append({
            "tool": r["command"],
            "mean": r["mean"],
            "stddev": r["stddev"],
            "median": r["median"],
            "min": r["min"],
            "max": r["max"],
        })
    return pd.DataFrame(rows)


def _extract_timestamp(name: str, prefix: str) -> str | None:
    match = re.match(rf"^{re.escape(prefix)}_(\d{{8}}_\d{{6}})(?:\.csv)?$", name)
    if match:
        return match.group(1)
    return None


def load_latest_run_manifest(results_dir: Path) -> dict | None:
    manifests = sorted(results_dir.glob("run_*.json"), reverse=True)
    if not manifests:
        return None
    with open(manifests[0]) as f:
        return json.load(f)


def infer_latest_bundle_timestamp(results_dir: Path) -> str | None:
    prefixes = ["single", "fullseq", "scale_region", "real", "bed", "multi", "memory", "rc_review"]
    coverage: dict[str, set[str]] = {}
    for prefix in prefixes:
        for path in results_dir.glob(f"{prefix}_*"):
            timestamp = _extract_timestamp(path.name, prefix)
            if timestamp is None:
                continue
            coverage.setdefault(timestamp, set()).add(prefix)
    if not coverage:
        return None
    return max(coverage.items(), key=lambda item: (len(item[1]), item[0]))[0]


def discover_dataset(results_dir: Path, prefix: str, bundle_timestamp: str | None = None) -> Path | None:
    if bundle_timestamp is not None:
        exact_dir = results_dir / f"{prefix}_{bundle_timestamp}"
        if exact_dir.exists():
            return exact_dir
        exact_csv = results_dir / f"{prefix}_{bundle_timestamp}.csv"
        if exact_csv.exists():
            return exact_csv
    return discover_latest(results_dir, prefix)


def discover_latest(results_dir: Path, prefix: str) -> Path | None:
    dirs = sorted(results_dir.glob(f"{prefix}_*"), reverse=True)
    dirs = [d for d in dirs if d.is_dir()]
    if dirs:
        return dirs[0]
    csvs = sorted(results_dir.glob(f"{prefix}_*.csv"), reverse=True)
    if csvs:
        return csvs[0]
    return None


def load_dir_jsons(results_dir: Path, prefix: str, bundle_timestamp: str | None = None) -> pd.DataFrame | None:
    d = discover_dataset(results_dir, prefix, bundle_timestamp)
    if not d or not d.is_dir():
        return None
    frames = []
    for jf in sorted(d.glob("*.json")):
        df = load_hyperfine_json(jf)
        df["benchmark"] = jf.stem
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_memory(results_dir: Path, bundle_timestamp: str | None = None) -> pd.DataFrame | None:
    p = discover_dataset(results_dir, "memory", bundle_timestamp)
    if not p:
        return None
    df = pd.read_csv(p)
    # CSV columns: tool,time_s,mem_kb,major_faults,minor_faults,region_type
    if "mem_kb" in df.columns:
        df["mem_mb"] = df["mem_kb"] / 1024.0
    return df


def load_multi(results_dir: Path, bundle_timestamp: str | None = None) -> pd.DataFrame | None:
    """Load multi-region benchmark results (bench [6]).

    Files come from results/multi_<timestamp>/<N>regions.json where N is
    the region count (1, 10, 50, 100).
    """
    d = discover_dataset(results_dir, "multi", bundle_timestamp)
    if not d or not d.is_dir():
        return None
    frames = []
    for jf in sorted(d.glob("*.json")):
        df = load_hyperfine_json(jf)
        # stem is e.g. "10regions" or "100regions"
        stem = jf.stem  # e.g. "10regions"
        count_str = stem.replace("regions", "")
        try:
            count = int(count_str)
        except ValueError:
            count = 0
        df["region_count"] = count
        df["benchmark"] = stem
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_bed(results_dir: Path, bundle_timestamp: str | None = None) -> pd.DataFrame | None:
    d = discover_dataset(results_dir, "bed", bundle_timestamp)
    if not d or not d.is_dir():
        return None
    frames = []
    for jf in sorted(d.glob("*.json")):
        df = load_hyperfine_json(jf)
        df["benchmark"] = jf.stem
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_rc_review(results_dir: Path, bundle_timestamp: str | None = None) -> pd.DataFrame | None:
    d = discover_dataset(results_dir, "rc_review", bundle_timestamp)
    if not d or not d.is_dir():
        return None
    json_dir = d / "json"
    if not json_dir.is_dir():
        return None
    frames = []
    for jf in sorted(json_dir.glob("*.json")):
        df = load_hyperfine_json(jf)
        df["slice"] = jf.stem
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_rc_rss(results_dir: Path, bundle_timestamp: str | None = None) -> pd.DataFrame | None:
    d = discover_dataset(results_dir, "rc_review", bundle_timestamp)
    if not d or not d.is_dir():
        return None
    rss_path = d / "rss.tsv"
    if not rss_path.is_file():
        return None
    return pd.read_csv(rss_path, sep="\t")


# ══════════════════════════════════════════════════════════════════════
#  Figures
# ══════════════════════════════════════════════════════════════════════

def _save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _bar_labels(ax, bars, vals_ms):
    """Annotate bars with ms values."""
    for bar, v in zip(bars, vals_ms):
        if v > 0:
            label = f"{v:.2f}" if v < 10 else f"{v:.1f}"
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.03,
                    label, ha="center", va="bottom", fontsize=7, rotation=45)


def _prep(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["tool_base"] = df["tool"].apply(_strip_tool_prefix)
    return df


def fig_latency_panels(df: pd.DataFrame, out: Path) -> Path:
    """Two-panel plot: x-axis = file size, bars = tools, subplots = region size.

    Handles benchmark names like 'size_10mb_100bp', 'size_50mb_10kbp'.
    """
    df = _prep(df)
    df["file_mb"] = df["benchmark"].str.extract(r"size_(\d+)mb").astype(float)
    df["region_tag"] = df["benchmark"].str.extract(r"_(\d+k?bp)$").iloc[:, 0]

    region_tags = sorted(
        df["region_tag"].dropna().unique(),
        key=lambda x: int(re.sub(r"k?bp", "", x)) * (1000 if "k" in x else 1),
    )
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    n = len(region_tags)

    fig, axes = plt.subplots(1, max(n, 1), figsize=(5 * max(n, 1) + 1, 5), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, rtag in zip(axes, region_tags):
        rdf = df[df["region_tag"] == rtag]
        file_sizes = sorted(rdf["file_mb"].dropna().unique())
        x = range(len(file_sizes))
        width = 0.8 / max(len(tools), 1)

        for i, t in enumerate(tools):
            tdf = rdf[rdf["tool_base"] == t]
            vals_ms = []
            for mb in file_sizes:
                row = tdf[tdf["file_mb"] == mb]
                vals_ms.append(row["mean"].values[0] * 1000 if len(row) else 0)
            offset = (i - len(tools) / 2 + 0.5) * width
            bars = ax.bar([xi + offset for xi in x], vals_ms, width,
                          label=t, color=COLORS.get(t, "#888"), alpha=0.9)
            _bar_labels(ax, bars, vals_ms)

        rtag_label = rtag.replace("kbp", " kbp").replace("bp", " bp")
        ax.set_xticks(list(x))
        ax.set_xticklabels([f"{int(mb)} MB" for mb in file_sizes])
        ax.set_xlabel("File Size")
        ax.set_ylabel("Time (ms)")
        ax.set_title(f"Region size: {rtag_label}")
        ax.legend(fontsize=8)
        _style_axes(ax)

    fig.suptitle("Single-Region Extraction Latency", fontweight="bold", fontsize=13)
    plt.tight_layout()
    return _save(fig, out)


def fig_fullseq(df: pd.DataFrame, out: Path) -> Path:
    """Full-sequence extraction: x-axis = file size (MB), bars by tool."""
    df = _prep(df)
    df["file_mb"] = df["benchmark"].str.extract(r"(\d+)mb").astype(float)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    file_sizes = sorted(df["file_mb"].dropna().unique())

    fig, ax = plt.subplots(figsize=(max(7, len(file_sizes) * 2), 5))
    x = range(len(file_sizes))
    width = 0.8 / max(len(tools), 1)

    for i, t in enumerate(tools):
        tdf = df[df["tool_base"] == t]
        vals_ms = []
        for mb in file_sizes:
            row = tdf[tdf["file_mb"] == mb]
            vals_ms.append(row["mean"].values[0] * 1000 if len(row) else 0)
        offset = (i - len(tools) / 2 + 0.5) * width
        bars = ax.bar([xi + offset for xi in x], vals_ms, width,
                      label=t, color=COLORS.get(t, "#888"), alpha=0.9)
        _bar_labels(ax, bars, vals_ms)

    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{int(mb)} MB" for mb in file_sizes])
    ax.set_xlabel("File Size")
    ax.set_ylabel("Time (ms)")
    ax.set_title("Full-Sequence Extraction by File Size")
    ax.legend(fontsize=8)
    _style_axes(ax)
    return _save(fig, out)


def fig_scaling_region(df: pd.DataFrame, out: Path) -> Path:
    """Region-size scaling: x-axis = region bp (log), lines by tool."""
    df = _prep(df)
    df["region_bp"] = df["benchmark"].str.extract(r"^(\d+)").astype(float)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]

    fig, ax = plt.subplots(figsize=(9, 5))
    for t in tools:
        tdf = df[df["tool_base"] == t].sort_values("region_bp")
        ax.plot(tdf["region_bp"], tdf["mean"] * 1000, "o-", label=t,
                color=COLORS.get(t, "#888"), linewidth=2, markersize=6)

    ax.set_xscale("log")
    ax.set_xlabel("Region Size")
    ax.set_ylabel("Time (ms)")
    ax.set_title("GET Latency vs Region Size")
    ax.legend(fontsize=9)
    _style_axes(ax)
    # Human-readable x ticks
    unique_bps = sorted(df["region_bp"].dropna().unique())
    ax.set_xticks(unique_bps)
    tick_labels = []
    for v in unique_bps:
        if v >= 1_000_000:
            tick_labels.append(f"{int(v/1_000_000)} Mbp")
        elif v >= 1_000:
            tick_labels.append(f"{int(v/1_000)} kbp")
        else:
            tick_labels.append(f"{int(v)} bp")
    ax.set_xticklabels(tick_labels, fontsize=8)
    return _save(fig, out)


def fig_speedup_vs(df: pd.DataFrame, baseline: str, title: str, out: Path) -> Path | None:
    """Horizontal bar chart: speedup of z-fasta vs a named baseline tool."""
    df = _prep(df)
    df["bench_human"] = df["benchmark"].apply(_human_bench)
    benchmarks = df["bench_human"].unique()

    rows = []
    for b in benchmarks:
        bdf = df[df["bench_human"] == b]
        ref = bdf[bdf["tool_base"] == baseline]
        zf = bdf[bdf["tool_base"] == "z-fasta"]
        if len(ref) and len(zf) and zf["mean"].values[0] > 0:
            speedup = ref["mean"].values[0] / zf["mean"].values[0]
            rows.append({"bench": b, "speedup": speedup})

    if not rows:
        return None
    sdf = pd.DataFrame(rows).sort_values("speedup")
    colors = [COLORS["z-fasta"] if v >= 1 else "#C62828" for v in sdf["speedup"]]

    fig, ax = plt.subplots(figsize=(7, max(3, len(rows) * 0.7)))
    bars = ax.barh(sdf["bench"], sdf["speedup"], color=colors, height=0.5, alpha=0.9)
    for bar, val in zip(bars, sdf["speedup"]):
        sign = "faster" if val >= 1 else "slower"
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                f"{val:.1f}x {sign}", va="center", fontweight="bold", fontsize=9)
    ax.axvline(x=1, color="gray", linestyle="--", alpha=0.6, linewidth=1.5)
    ax.set_xlabel(f"Speedup (z-fasta vs {baseline}; >1 = z-fasta wins)")
    ax.set_title(title)
    ax.set_xlim(0, max(sdf["speedup"].max() * 1.3, 2))
    ax.grid(axis="x", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return _save(fig, out)


def fig_real_datasets(df: pd.DataFrame, out: Path) -> Path:
    """Real dataset benchmarks grouped by dataset+region, bars by tool."""
    df = _prep(df)
    df["bench_human"] = df["benchmark"].apply(_human_bench)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    benchmarks = sorted(df["bench_human"].unique(), key=workload_sort_key)

    fig, ax = plt.subplots(figsize=(max(10, len(benchmarks) * 1.8), 5))
    x = range(len(benchmarks))
    width = 0.8 / max(len(tools), 1)

    for i, t in enumerate(tools):
        tdf = df[df["tool_base"] == t]
        vals_ms = []
        for b in benchmarks:
            row = tdf[tdf["bench_human"] == b]
            vals_ms.append(row["mean"].values[0] * 1000 if len(row) else 0)
        offset = (i - len(tools) / 2 + 0.5) * width
        bars = ax.bar([xi + offset for xi in x], vals_ms, width,
                      label=t, color=COLORS.get(t, "#888"), alpha=0.9)
        _bar_labels(ax, bars, vals_ms)

    ax.set_xticks(list(x))
    ax.set_xticklabels(benchmarks, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Time (ms)")
    ax.set_title("Real Dataset Extraction Performance")
    ax.legend(fontsize=8)
    _style_axes(ax)
    return _save(fig, out)


def fig_bed_batch(df: pd.DataFrame, out: Path) -> Path:
    df = _prep(df)
    df["bench_human"] = df["benchmark"].apply(_human_bench)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    benchmarks = sorted(df["bench_human"].unique(), key=workload_sort_key)

    fig, ax = plt.subplots(figsize=(max(10, len(benchmarks) * 1.8), 5))
    x = range(len(benchmarks))
    width = 0.8 / max(len(tools), 1)

    for i, t in enumerate(tools):
        tdf = df[df["tool_base"] == t]
        vals_ms = []
        for b in benchmarks:
            row = tdf[tdf["bench_human"] == b]
            vals_ms.append(row["mean"].values[0] * 1000 if len(row) else 0)
        offset = (i - len(tools) / 2 + 0.5) * width
        bars = ax.bar([xi + offset for xi in x], vals_ms, width,
                      label=t, color=COLORS.get(t, "#888"), alpha=0.9)
        _bar_labels(ax, bars, vals_ms)

    ax.set_xticks(list(x))
    ax.set_xticklabels(benchmarks, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Time (ms)")
    ax.set_title("BED Batch Extraction Performance")
    ax.legend(fontsize=8)
    _style_axes(ax)
    return _save(fig, out)


def fig_rc_slice(df: pd.DataFrame, slice_name: str, title: str, out: Path) -> Path | None:
    sdf = df[df["slice"] == slice_name].copy()
    if sdf.empty:
        return None

    order_map = {
        "orientation_small": ["forward", "zf-rc", "samtools-rc", "reverse-only", "complement-only", "rc-annotate"],
        "orientation_medium": ["forward", "zf-rc", "samtools-rc", "reverse-only", "complement-only", "rc-annotate"],
        "orientation_large": ["forward", "zf-rc", "samtools-rc", "reverse-only", "complement-only", "rc-annotate"],
        "full_sequence": ["forward", "zf-rc", "samtools-rc", "rc-annotate"],
        "bed_batch": ["bed-forward", "bed-zf-rc", "bed-honor-strand-zf-rc", "bed-bedtools-seqtk-rc", "bed-honor-strand-bedtools-seqtk-rc", "bed-reverse-only", "bed-complement-only"],
        "multi_region": ["multi10-forward", "multi10-zf-rc", "multi10-samtools-rc", "multi50-forward", "multi50-zf-rc", "multi50-samtools-rc", "multi100-forward", "multi100-zf-rc", "multi100-samtools-rc"],
    }
    ordered_tools = [name for name in order_map.get(slice_name, []) if name in set(sdf["tool"])]
    if not ordered_tools:
        ordered_tools = list(sdf["tool"])

    sdf["tool"] = pd.Categorical(sdf["tool"], categories=ordered_tools, ordered=True)
    sdf = sdf.sort_values("tool")

    labels = [rc_tool_label(name) for name in sdf["tool"].astype(str)]
    vals_ms = (sdf["mean"] * 1000).tolist()
    colors = [rc_color(name) for name in sdf["tool"].astype(str)]

    fig, ax = plt.subplots(figsize=(max(9, len(labels) * 1.15), 4.8))
    bars = ax.bar(range(len(labels)), vals_ms, color=colors, alpha=0.95, edgecolor="#EDE7DF", linewidth=0.8)
    _bar_labels(ax, bars, vals_ms)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=8)
    ax.set_ylabel("Time (ms)")
    ax.set_title(title)
    _style_axes(ax)
    return _save(fig, out)


def fig_rc_orientation_overview(df: pd.DataFrame, out: Path) -> Path | None:
    slices = [
        ("orientation_small", "100 bp"),
        ("orientation_medium", "10 kbp"),
        ("orientation_large", "250 kbp"),
        ("full_sequence", "full seq"),
    ]
    mode_order = ["forward", "zf-rc", "samtools-rc", "reverse-only", "complement-only", "rc-annotate"]
    mode_labels = {
        "forward": "forward",
        "zf-rc": "z-fasta --rc",
        "samtools-rc": "samtools -i",
        "reverse-only": "reverse-only",
        "complement-only": "complement-only",
        "rc-annotate": "--rc + annotate",
    }

    plotted = False
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharey=False)
    axes = axes.flatten()
    for ax, (slice_name, label) in zip(axes, slices):
        sdf = df[df["slice"] == slice_name].copy()
        present_modes = [mode for mode in mode_order if mode in set(sdf["tool"])]
        if not present_modes:
            ax.set_visible(False)
            continue
        plotted = True
        sdf["tool"] = pd.Categorical(sdf["tool"], categories=present_modes, ordered=True)
        sdf = sdf.sort_values("tool")
        vals_ms = (sdf["mean"] * 1000).tolist()
        bars = ax.bar(
            range(len(present_modes)),
            vals_ms,
            color=[rc_color(name) for name in present_modes],
            edgecolor="#EDE7DF",
            linewidth=0.8,
        )
        _bar_labels(ax, bars, vals_ms)
        ax.set_xticks(range(len(present_modes)))
        ax.set_xticklabels([mode_labels[name] for name in present_modes], rotation=18, ha="right", fontsize=8)
        ax.set_title(label)
        ax.set_ylabel("Time (ms)")
        _style_axes(ax)

    if not plotted:
        plt.close(fig)
        return None

    fig.suptitle("RC Positional Comparison", fontsize=14, fontweight="bold")
    plt.tight_layout()
    return _save(fig, out)


def fig_rc_speedup(df: pd.DataFrame, out: Path) -> Path | None:
    slices = [
        ("orientation_small", "100 bp"),
        ("orientation_medium", "10 kbp"),
        ("orientation_large", "250 kbp"),
        ("full_sequence", "full seq"),
    ]
    rows = []
    for slice_name, label in slices:
        sdf = df[df["slice"] == slice_name]
        zf = sdf[sdf["tool"] == "zf-rc"]
        st = sdf[sdf["tool"] == "samtools-rc"]
        if zf.empty or st.empty:
            continue
        rows.append({"slice": label, "speedup": st["mean"].values[0] / zf["mean"].values[0]})
    if not rows:
        return None
    speed_df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    bars = ax.barh(speed_df["slice"], speed_df["speedup"], color=RC_COLORS["z-fasta-rc"], alpha=0.95)
    for bar, value in zip(bars, speed_df["speedup"]):
        ax.text(value + 0.05, bar.get_y() + bar.get_height() / 2, f"{value:.1f}x", va="center", fontsize=9, fontweight="bold")
    ax.axvline(1, color="#7C7468", linestyle="--", linewidth=1.2)
    ax.set_xlabel("samtools -i time / z-fasta --rc time")
    ax.set_title("RC Speedup vs samtools")
    _style_axes(ax)
    return _save(fig, out)


def fig_rc_batch_compare(df: pd.DataFrame, out: Path) -> Path | None:
    multi_df = df[df["slice"] == "multi_region"].copy()
    bed_df = df[df["slice"] == "bed_batch"].copy()
    if multi_df.empty and bed_df.empty:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

    if not multi_df.empty:
        ax = axes[0]
        groups = []
        for prefix in ["multi10", "multi50", "multi100"]:
            present = [f"{prefix}-forward", f"{prefix}-zf-rc", f"{prefix}-samtools-rc"]
            if any(name in set(multi_df["tool"]) for name in present):
                groups.append(prefix.replace("multi", "") + " regions")
        x_positions = range(len(groups))
        width = 0.22
        series = [
            ("forward", "z-fasta forward"),
            ("zf-rc", "z-fasta --rc"),
            ("samtools-rc", "samtools -i"),
        ]
        for offset_idx, (suffix, label) in enumerate(series):
            vals = []
            for group in groups:
                prefix = f"multi{group.split()[0]}"
                row = multi_df[multi_df["tool"] == f"{prefix}-{suffix}"]
                vals.append(row["mean"].values[0] * 1000 if len(row) else 0)
            bars = ax.bar([x + (offset_idx - 1) * width for x in x_positions], vals, width=width,
                          label=label, color=rc_color(f"multi10-{suffix}"), edgecolor="#EDE7DF", linewidth=0.8)
            _bar_labels(ax, bars, vals)
        ax.set_xticks(list(x_positions))
        ax.set_xticklabels(groups)
        ax.set_ylabel("Time (ms)")
        ax.set_title("Multi-region RC")
        _style_axes(ax)
        ax.legend(fontsize=8)
    else:
        axes[0].set_visible(False)

    if not bed_df.empty:
        ax = axes[1]
        group_defs = [
            ("BED default", ["bed-forward", "bed-zf-rc", "bed-bedtools-seqtk-rc"]),
            ("BED honor-strand", ["bed-forward", "bed-honor-strand-zf-rc", "bed-honor-strand-bedtools-seqtk-rc"]),
        ]
        active_groups = [(label, tools) for label, tools in group_defs if any(tool in set(bed_df["tool"]) for tool in tools)]
        x_positions = range(len(active_groups))
        width = 0.22
        series_labels = [
            (0, "z-fasta forward"),
            (1, "z-fasta RC path"),
            (2, "bedtools + seqtk"),
        ]
        for offset_idx, legend_label in series_labels:
            vals = []
            colors = []
            for _, tools in active_groups:
                tool_name = tools[offset_idx]
                row = bed_df[bed_df["tool"] == tool_name]
                vals.append(row["mean"].values[0] * 1000 if len(row) else 0)
                colors.append(rc_color(tool_name))
            bars = ax.bar([x + (offset_idx - 1) * width for x in x_positions], vals, width=width,
                          label=legend_label, color=colors, edgecolor="#EDE7DF", linewidth=0.8)
            _bar_labels(ax, bars, vals)
        ax.set_xticks(list(x_positions))
        ax.set_xticklabels([label for label, _ in active_groups])
        ax.set_ylabel("Time (ms)")
        ax.set_title("BED RC vs external composition")
        _style_axes(ax)
        ax.legend(fontsize=8)
    else:
        axes[1].set_visible(False)

    fig.suptitle("RC Batch Comparison", fontsize=14, fontweight="bold")
    plt.tight_layout()
    return _save(fig, out)


# ══════════════════════════════════════════════════════════════════════
#  Markdown
# ══════════════════════════════════════════════════════════════════════

def _strip_tool_prefix(name: str) -> str:
    """Extract base tool name from hyperfine label like '10mb_full_z-fasta'."""
    for base in TOOL_ORDER:
        if name.endswith(base):
            return base
        if base in name:
            return base
    return name


def format_timing(mean: float | None, stddev: float | None) -> str:
    if mean is None:
        return "-"
    if stddev is None:
        return f"{mean * 1000:.2f} ms"
    return f"{mean * 1000:.2f} ms +/- {stddev * 1000:.2f}"


def format_ratio(baseline: float | None, zf: float | None) -> str:
    if baseline is None or zf is None or zf <= 0:
        return "-"
    ratio = baseline / zf
    if ratio >= 1:
        return f"**{ratio:.1f}x**"
    return f"{ratio:.2f}x"


def tool_display_name(tool: str) -> str:
    return TOOL_DISPLAY.get(tool, tool)


def ordered_present_tools(df: pd.DataFrame, preferred: list[str] | None = None) -> list[str]:
    present = set(df["tool_base"])
    order = preferred if preferred is not None else TOOL_ORDER
    return [tool for tool in order if tool in present]


def md_benchmark_summary_table(
    df: pd.DataFrame,
    group_col: str,
    *,
    human_names: bool = True,
    group_name: str = "Workload",
    baseline_tool: str = "samtools",
    preferred_tools: list[str] | None = None,
    baseline_ratio_exempt: set[str] | None = None,
) -> str:
    df = df.copy()
    df["tool_base"] = df["tool"].apply(_strip_tool_prefix)
    if human_names:
        df[group_col] = df[group_col].apply(_human_bench)

    groups = sorted(df[group_col].unique(), key=workload_sort_key)
    tools = ordered_present_tools(df, preferred_tools)
    rows = []
    for group in groups:
        gdf = df[df[group_col] == group]
        row = {group_name: group}
        zf_mean = None
        for tool in tools:
            match = gdf[gdf["tool_base"] == tool]
            if len(match):
                mean = match["mean"].values[0]
                stddev = match["stddev"].values[0]
                row[tool_display_name(tool)] = format_timing(mean, stddev)
                if tool == "z-fasta":
                    zf_mean = mean
            else:
                row[tool_display_name(tool)] = "-"

        baseline_match = gdf[gdf["tool_base"] == baseline_tool]
        baseline_mean = baseline_match["mean"].values[0] if len(baseline_match) else None
        exempt_baseline_ratio = (
            baseline_ratio_exempt is not None
            and (
                group in baseline_ratio_exempt
                or str(group).endswith("_stranded")
                or str(group).endswith("(stranded)")
            )
        )
        if exempt_baseline_ratio:
            row[f"{baseline_tool} / z-fasta"] = "n/a (strandless)"
        else:
            row[f"{baseline_tool} / z-fasta"] = format_ratio(baseline_mean, zf_mean)

        competitors = []
        for tool in tools:
            if tool == "z-fasta":
                continue
            if exempt_baseline_ratio and tool == baseline_tool:
                continue
            match = gdf[gdf["tool_base"] == tool]
            if len(match):
                competitors.append((tool, match["mean"].values[0]))
        if competitors and zf_mean is not None:
            best_tool, best_mean = min(competitors, key=lambda item: item[1])
            row["best alt / z-fasta"] = f"{format_ratio(best_mean, zf_mean)} {tool_display_name(best_tool)}"
        else:
            row["best alt / z-fasta"] = "-"
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_memory_table(df: pd.DataFrame) -> str:
    # CSV columns: tool, time_s, mem_kb, major_faults, minor_faults, region_type
    df = df.copy()
    col_map = {"tool": "Tool", "region_type": "Region", "time_s": "Time (s)",
               "mem_kb": "MaxRSS (KB)", "mem_mb": "MaxRSS (MB)",
               "major_faults": "Maj. Faults", "minor_faults": "Min. Faults"}
    want = ["tool", "region_type", "time_s", "mem_kb", "major_faults"]
    cols = [c for c in want if c in df.columns]
    display = df[cols].rename(columns=col_map)
    if "MaxRSS (KB)" in display.columns:
        display["MaxRSS (KB)"] = display["MaxRSS (KB)"].apply(lambda v: f"{int(v):,}")
    return display.to_markdown(index=False, floatfmt=".4f")


def md_memory_summary_table(df: pd.DataFrame) -> str:
    rows = []
    for region in ["small", "large", "full"]:
        rdf = df[df["region_type"] == region]
        zf = rdf[rdf["tool"] == "z-fasta"]
        sam = rdf[rdf["tool"] == "samtools"]
        row = {"Region": region}
        zf_mem = zf["mem_kb"].values[0] if len(zf) else None
        sam_mem = sam["mem_kb"].values[0] if len(sam) else None
        row["z-fasta RSS"] = f"{int(zf_mem):,} KB" if zf_mem is not None else "-"
        row["samtools RSS"] = f"{int(sam_mem):,} KB" if sam_mem is not None else "-"
        row["samtools / z-fasta RSS"] = format_ratio(sam_mem, zf_mem)
        competitors = rdf[rdf["tool"] != "z-fasta"]
        if len(competitors) and zf_mem is not None:
            best = competitors.loc[competitors["mem_kb"].idxmin()]
            row["lowest alt / z-fasta RSS"] = f"{format_ratio(best['mem_kb'], zf_mem)} {tool_display_name(best['tool'])}"
        else:
            row["lowest alt / z-fasta RSS"] = "-"
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def rc_slice_label(name: str) -> str:
    return {
        "orientation_small": "100 bp region",
        "orientation_medium": "10 kbp region",
        "orientation_large": "large region",
        "full_sequence": "full sequence",
        "multi_region": "multi-region",
        "bed_batch": "BED batches",
    }.get(name, name.replace("_", " "))


def rc_tool_label(name: str) -> str:
    mapping = {
        "forward": "z-fasta forward",
        "zf-rc": "z-fasta --rc",
        "rc": "z-fasta --rc",
        "samtools-rc": "samtools -i",
        "reverse-only": "z-fasta --reverse-only",
        "complement-only": "z-fasta --complement-only",
        "rc-annotate": "z-fasta --rc --annotate-rc",
        "multi10-forward": "10-region forward",
        "multi10-zf-rc": "10-region z-fasta --rc",
        "multi10-rc": "10-region z-fasta --rc",
        "multi10-samtools-rc": "10-region samtools -i",
        "multi50-forward": "50-region forward",
        "multi50-zf-rc": "50-region z-fasta --rc",
        "multi50-rc": "50-region z-fasta --rc",
        "multi50-samtools-rc": "50-region samtools -i",
        "multi100-forward": "100-region forward",
        "multi100-zf-rc": "100-region z-fasta --rc",
        "multi100-rc": "100-region z-fasta --rc",
        "multi100-samtools-rc": "100-region samtools -i",
        "bed-forward": "BED forward",
        "bed-zf-rc": "BED z-fasta --rc",
        "bed-rc": "BED z-fasta --rc",
        "bed-reverse-only": "BED --reverse-only",
        "bed-complement-only": "BED --complement-only",
        "bed-bedtools-seqtk-rc": "BED bedtools + seqtk rc",
        "bed-honor-strand-zf-rc": "BED honor-strand + z-fasta --rc",
        "bed-honor-strand-rc": "BED honor-strand + z-fasta --rc",
        "bed-honor-strand-bedtools-seqtk-rc": "BED honor-strand + bedtools + seqtk rc",
        "zf_rc": "z-fasta --rc",
        "multi_b_zf_rc": "multi-region z-fasta --rc",
        "multi_b_rc": "multi-region z-fasta --rc",
        "bed_honor_strand_zf_rc": "BED honor-strand + z-fasta --rc",
        "bed_honor_strand_rc": "BED honor-strand + z-fasta --rc",
        "samtools_rc": "samtools -i",
        "multi_b_samtools_rc": "multi-region samtools -i",
        "bed_honor_strand_bedtools_seqtk_rc": "BED honor-strand + bedtools + seqtk rc",
        "rc_annotate": "z-fasta --rc --annotate-rc",
        "multi_b_forward": "multi-region forward",
        "bed_forward": "BED forward",
    }
    return mapping.get(name, name.replace("_", " "))


def md_rc_slice_table(df: pd.DataFrame, slice_name: str) -> str:
    sdf = df[df["slice"] == slice_name].copy()
    if sdf.empty:
        return "_No RC benchmark data available._"
    sdf["Mode"] = sdf["tool"].apply(rc_tool_label)
    sdf["Mean (ms)"] = sdf["mean"].apply(lambda v: f"{v * 1000:.2f}")
    sdf["Stddev (ms)"] = sdf["stddev"].apply(lambda v: f"{v * 1000:.2f}")
    sdf = sdf[["Mode", "Mean (ms)", "Stddev (ms)"]]
    return sdf.to_markdown(index=False)


def md_rc_summary_table(df: pd.DataFrame) -> str:
    rows = []
    for slice_name in ["orientation_small", "orientation_medium", "orientation_large", "full_sequence"]:
        sdf = df[df["slice"] == slice_name]
        if sdf.empty:
            continue
        row = {"Workload": rc_slice_label(slice_name)}
        forward_mean = None
        rc_mean = None
        samtools_mean = None
        for tool_name, col_name in [
            ("forward", "z-fasta forward"),
            ("zf-rc", "z-fasta --rc"),
            ("samtools-rc", "samtools -i"),
            ("reverse-only", "reverse-only"),
            ("complement-only", "complement-only"),
            ("rc-annotate", "--rc + annotate"),
        ]:
            match = sdf[sdf["tool"] == tool_name]
            row[col_name] = format_timing(match["mean"].values[0], match["stddev"].values[0]) if len(match) else "-"
            if tool_name == "forward" and len(match):
                forward_mean = match["mean"].values[0]
            if tool_name == "zf-rc" and len(match):
                rc_mean = match["mean"].values[0]
            if tool_name == "samtools-rc" and len(match):
                samtools_mean = match["mean"].values[0]
        row["RC / forward"] = format_ratio(rc_mean, forward_mean)
        row["samtools / z-fasta --rc"] = format_ratio(samtools_mean, rc_mean)
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_rc_batch_table(df: pd.DataFrame) -> str:
    rows = []
    multi_df = df[df["slice"] == "multi_region"]
    for count in [10, 50, 100]:
        forward = multi_df[multi_df["tool"] == f"multi{count}-forward"]
        rc = multi_df[multi_df["tool"] == f"multi{count}-zf-rc"]
        sam = multi_df[multi_df["tool"] == f"multi{count}-samtools-rc"]
        if len(forward) or len(rc) or len(sam):
            forward_mean = forward["mean"].values[0] if len(forward) else None
            rc_mean = rc["mean"].values[0] if len(rc) else None
            sam_mean = sam["mean"].values[0] if len(sam) else None
            rows.append({
                "Workload": f"multi-region ({count})",
                "z-fasta forward": format_timing(forward_mean, forward["stddev"].values[0] if len(forward) else None),
                "z-fasta --rc": format_timing(rc_mean, rc["stddev"].values[0] if len(rc) else None),
                "baseline": format_timing(sam_mean, sam["stddev"].values[0] if len(sam) else None),
                "RC / forward": format_ratio(rc_mean, forward_mean),
                "baseline / z-fasta --rc": format_ratio(sam_mean, rc_mean),
            })

    bed_df = df[df["slice"] == "bed_batch"]
    for label, forward_name, rc_name, baseline_name in [
        ("BED default", "bed-forward", "bed-zf-rc", "bed-bedtools-seqtk-rc"),
        ("BED honor-strand", "bed-forward", "bed-honor-strand-zf-rc", "bed-honor-strand-bedtools-seqtk-rc"),
    ]:
        forward = bed_df[bed_df["tool"] == forward_name]
        rc = bed_df[bed_df["tool"] == rc_name]
        baseline = bed_df[bed_df["tool"] == baseline_name]
        if len(forward) or len(rc) or len(baseline):
            forward_mean = forward["mean"].values[0] if len(forward) else None
            rc_mean = rc["mean"].values[0] if len(rc) else None
            baseline_mean = baseline["mean"].values[0] if len(baseline) else None
            rows.append({
                "Workload": label,
                "z-fasta forward": format_timing(forward_mean, forward["stddev"].values[0] if len(forward) else None),
                "z-fasta --rc": format_timing(rc_mean, rc["stddev"].values[0] if len(rc) else None),
                "baseline": format_timing(baseline_mean, baseline["stddev"].values[0] if len(baseline) else None),
                "RC / forward": format_ratio(rc_mean, forward_mean),
                "baseline / z-fasta --rc": format_ratio(baseline_mean, rc_mean),
            })
    return pd.DataFrame(rows).to_markdown(index=False)


def md_rc_rss_table(df: pd.DataFrame) -> str:
    rdf = df.copy()
    rdf["Mode"] = rdf["label"].apply(rc_tool_label)
    rdf["Elapsed"] = rdf["elapsed"]
    rdf["MaxRSS (KB)"] = rdf["maxrss_kb"].apply(lambda v: f"{int(v):,}")
    zf_rc_row = rdf[rdf["label"].isin(["zf_rc", "rc", "z-fasta --rc"]) ]
    zf_rc_mem = zf_rc_row["maxrss_kb"].values[0] if len(zf_rc_row) else None
    rdf["vs z-fasta --rc"] = rdf["maxrss_kb"].apply(lambda value: format_ratio(value, zf_rc_mem))
    return rdf[["Mode", "Elapsed", "MaxRSS (KB)", "vs z-fasta --rc"]].to_markdown(index=False)


def rc_speedup_text(df: pd.DataFrame, slice_name: str, zf_name: str, baseline_name: str) -> str | None:
    sdf = df[df["slice"] == slice_name]
    zf = sdf[sdf["tool"] == zf_name]
    baseline = sdf[sdf["tool"] == baseline_name]
    if zf.empty or baseline.empty:
        return None
    ratio = baseline["mean"].values[0] / zf["mean"].values[0]
    return f"**{ratio:.1f}x**"


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    script_dir = Path(__file__).resolve().parent
    results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else script_dir / "results"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    run_manifest = load_latest_run_manifest(results_dir)
    bundle_timestamp = run_manifest["timestamp"] if run_manifest else None

    report = []
    report.append("<!-- markdownlint-disable MD012 MD024 -->\n")
    report.append("# z-fasta GET Benchmark Report\n")
    report.append("_Auto-generated by `generate_report.py`_\n")
    if bundle_timestamp:
        report.append("\n## Run Context\n")
        report.append(f"> Report anchored to run bundle `{bundle_timestamp}` so all sections come from one benchmark family instead of mixing timestamps.\n")
        if run_manifest is not None:
            produced = ", ".join(run_manifest.get("produced", []))
            report.append(f"> Runs: `{run_manifest.get('runs')}` | Warmup: `{run_manifest.get('warmup')}` | Mode: `{run_manifest.get('bench_mode')}` | Produced: `{produced}`\n")
    else:
        inferred = infer_latest_bundle_timestamp(results_dir)
        report.append("\n## Run Context\n")
        if inferred:
            report.append(f"> No run manifest is present yet. This report is using the latest available slice per section, so sections may come from different timestamps. The most complete recent timestamp family is `{inferred}`.\n")
        else:
            report.append("> No run manifest is present yet. This report is using the latest available slice per section, so sections may come from different timestamps.\n")

    # ── 1. Single-region latency ───────────────────────────────────
    single_df = load_dir_jsons(results_dir, "single", bundle_timestamp)
    if single_df is not None and len(single_df):
        report.append("\n## Single-Region Extraction Latency\n")
        report.append("> All measurements on warm cache. "
                      "Each benchmark extracts a region from the **middle** of the target file.\n")
        report.append(md_benchmark_summary_table(single_df, "benchmark"))
        report.append("")
        fig_latency_panels(single_df, figures_dir / "single_latency.png")
        report.append("\n![Single-Region Latency](results/figures/single_latency.png)\n")

        report.append("> **Interpretation:** z-fasta small-region extraction is startup-dominated: "
                  "the measured path is O(1) index lookup + `pread`, with ~1-2 ms "
                  "end-to-end CLI latency on this host. seqkit carries a ~14 ms Go runtime "
                  "startup cost that dominates at small regions. fastahack and samtools are "
                  "in the same low-millisecond range for indexed single-region access.\n")

        p = fig_speedup_vs(single_df, "samtools", "z-fasta GET Speedup vs samtools",
                           figures_dir / "single_speedup.png")
        if p:
            report.append("\n![Single-Region Speedup](results/figures/single_speedup.png)\n")

    # ── 2. Full-sequence extraction ────────────────────────────────
    fullseq_df = load_dir_jsons(results_dir, "fullseq", bundle_timestamp)
    if fullseq_df is not None and len(fullseq_df):
        report.append("\n## Full-Sequence Extraction\n")
        report.append(md_benchmark_summary_table(fullseq_df, "benchmark"))
        report.append("")
        fig_fullseq(fullseq_df, figures_dir / "fullseq.png")
        report.append("\n![Full-Sequence Extraction](results/figures/fullseq.png)\n")
        report.append("> **Note:** fastahack is faster than z-fasta for large full-sequence "
                      "extraction (≥50 MB). fastahack uses a simpler unbuffered write path "
                      "while z-fasta's buffered I/O adds overhead proportional to sequence "
                      "length. This is an optimization target for a future release.\n")

    # ── 3. Region-size scaling ─────────────────────────────────────
    scale_df = load_dir_jsons(results_dir, "scale_region", bundle_timestamp)
    if scale_df is not None and len(scale_df):
        report.append("\n## Region-Size Scaling\n")
        report.append(md_benchmark_summary_table(scale_df, "benchmark"))
        report.append("")
        fig_scaling_region(scale_df.copy(), figures_dir / "scaling_region.png")
        report.append("\n![Region-Size Scaling](results/figures/scaling_region.png)\n")
        report.append("> **Interpretation:** Latency is O(1) up to ~10 kbp. Process startup "
                      "and index lookup dominate below that threshold. Above ~100 kbp it grows linearly "
                      "with region size as the I/O read cost dominates. fastahack's output "
                      "advantage grows with region size due to its simpler write path.\n")

    # ── 4. Multi-region extraction ────────────────────────────────
    multi_df = load_multi(results_dir, bundle_timestamp)
    if multi_df is not None and len(multi_df):
        report.append("\n## Multi-Region Extraction\n")
        report.append("> Time to extract N regions in a single CLI call, loading the index once and streaming all results to stdout. "
                      "z-fasta uses O(1) offset lookup per region while samtools and seqtk pay higher per-call overhead.\n")
        multi_df = multi_df.copy()
        multi_df["benchmark"] = multi_df["region_count"].apply(lambda value: f"{value} region" if value == 1 else f"{value} regions")
        report.append(md_benchmark_summary_table(multi_df, "benchmark", human_names=False))
        report.append("")
        report.append("> **Note:** seqtk does not use an index for region extraction. It performs a full-file scan for every call regardless of region count, so its time is roughly constant here and reflects file I/O cost, not region count.\n")

    # ── 5. BED batch extraction ──────────────────────────────────
    bed_df = load_bed(results_dir, bundle_timestamp)
    if bed_df is not None and len(bed_df):
        report.append("\n## BED Batch Extraction\n")
        report.append(
            "> Synthetic BED batches of 100, 1K, 10K, and 100K regions on a single indexed FASTA. "
            "`z-fasta` runs with `--bed` and an explicit chunk size that forces multi-batch processing. "
            "`samtools faidx -r` ignores strand, so the stranded rows are not apples-to-apples against samtools: "
            "they show the extra work `z-fasta` and `bedtools getfasta -s` do for strand handling while samtools remains a forward-only fetch baseline.\n"
        )
        stranded_bed_rows = {
            "100 regions (stranded)",
            "1000 regions (stranded)",
            "10000 regions (stranded)",
            "100000 regions (stranded)",
        }
        report.append(md_benchmark_summary_table(
            bed_df,
            "benchmark",
            preferred_tools=["z-fasta", "samtools", "bedtools"],
            baseline_ratio_exempt=stranded_bed_rows,
        ))
        report.append("")
        report.append(
            "> **Interpretation:** Treat the stranded `samtools / z-fasta` ratio as an orientation-overhead reference, not as a fair feature-parity comparison. "
            "For strand-aware equivalence, compare `z-fasta --honor-strand --rc` against `bedtools getfasta -s | seqtk seq -r` in the RC section.\n"
        )
        fig_bed_batch(bed_df, figures_dir / "bed_batch.png")
        report.append("\n![BED Batch Extraction](results/figures/bed_batch.png)\n")

    # ── 6. Reverse/complement orientation ───────────────────────
    rc_df = load_rc_review(results_dir, bundle_timestamp)
    rc_rss_df = load_rc_rss(results_dir, bundle_timestamp)
    if rc_df is not None and len(rc_df):
        report.append("\n## Reverse / Complement Orientation\n")
        report.append(
            "> Correctness lives in `verify_rc.sh`; this section tracks the runtime and RSS cost of shipped orientation modes and compares them against external RC-capable baselines where those paths exist.\n"
        )

        report.append("\n### Positional RC Summary\n")
        report.append(md_rc_summary_table(rc_df))
        report.append("")
        pos_fig = fig_rc_orientation_overview(rc_df, figures_dir / "rc_orientation_overview.png")
        if pos_fig:
            report.append("\n![RC Positional Comparison](results/figures/rc_orientation_overview.png)\n")

        speedup_fig = fig_rc_speedup(rc_df, figures_dir / "rc_speedup.png")
        if speedup_fig:
            report.append("\n![RC Speedup vs samtools](results/figures/rc_speedup.png)\n")

        report.append("\n### Batch RC Summary\n")
        report.append(md_rc_batch_table(rc_df))
        report.append("")
        batch_fig = fig_rc_batch_compare(rc_df, figures_dir / "rc_batch_compare.png")
        if batch_fig:
            report.append("\n![RC Batch Comparison](results/figures/rc_batch_compare.png)\n")

        large_vs_samtools = rc_speedup_text(rc_df, "orientation_large", "zf-rc", "samtools-rc")
        bed_vs_baseline = rc_speedup_text(rc_df, "bed_batch", "bed-honor-strand-zf-rc", "bed-honor-strand-bedtools-seqtk-rc")
        if large_vs_samtools or bed_vs_baseline:
            report.append("\n### RC Takeaways\n")
            if large_vs_samtools:
                report.append(f"- Large-region `--rc` stays in the same low-millisecond class as forward extraction and remains {large_vs_samtools} faster than the local `samtools faidx -i` baseline on the synthetic review slice.")
            if bed_vs_baseline:
                report.append(f"- BED `--honor-strand --rc` remains a single-pass path and is {bed_vs_baseline} faster than the local `bedtools getfasta -s | seqtk seq -r` composition on the review slice.")
            report.append("- The inline RC review slice stays quick by default so these comparisons can be rerun in the normal edit loop without a separate benchmark entrypoint.")

        if rc_rss_df is not None and len(rc_rss_df):
            report.append("\n### RC RSS Snapshot\n")
            report.append(md_rc_rss_table(rc_rss_df))
            report.append("")

    # ── 7. Real datasets ──────────────────────────────────────────
    real_df = load_dir_jsons(results_dir, "real", bundle_timestamp)
    if real_df is not None and len(real_df):
        report.append("\n## Real Dataset Extraction\n")
        report.append(md_benchmark_summary_table(real_df, "benchmark"))
        report.append("")
        fig_real_datasets(real_df, figures_dir / "real_get.png")
        report.append("\n![Real Dataset Extraction](results/figures/real_get.png)\n")
    else:
        report.append("\n## Real Dataset Extraction\n")
        report.append("> _No real-dataset results found for the selected run bundle. Run `bench/shared/download_data.sh` and rerun benchmarks without `--skip-real` for a coherent real-data section._\n")

    # ── 8. Memory ────────────────────────────────────────────────
    mem_df = load_memory(results_dir, bundle_timestamp)
    if mem_df is not None and len(mem_df):
        report.append("\n## Memory Usage\n")
        report.append("> RSS measured with `/usr/bin/time -v`. For mmap-based tools (z-fasta, samtools, fastahack) RSS reflects file pages mapped by the OS, not heap allocation. `Major Faults` should be ~0 on warm cache.\n")
        report.append(md_memory_summary_table(mem_df))
        report.append("")
        report.append("\n### Detailed Memory Table\n")
        report.append(md_memory_table(mem_df))
        report.append("")

    # ── Write report ──────────────────────────────────────────────
    report_text = "\n".join(report) + "\n"
    report_text = re.sub(r"\n{3,}", "\n\n", report_text)
    report_text = report_text.rstrip() + "\n"
    report_path = script_dir / "REPORT.md"
    report_path.write_text(report_text)
    print(f"Report written to {report_path}")
    print(f"Figures in {figures_dir}")


if __name__ == "__main__":
    main()
