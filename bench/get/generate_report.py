#!/usr/bin/env python3
"""
z-fasta GET Benchmark Report Generator

Reads raw hyperfine JSON + CSV data from bench/get/results/,
produces Markdown report + PNG figures using pandas + matplotlib.

Usage:
    python3 bench/get/generate_report.py [results_dir]

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
    "seqkit": "#E65100",
    "fastahack": "#7B1FA2",
    "seqtk": "#00838F",
    "pyfaidx": "#F7A41D",
}
TOOL_ORDER = ["z-fasta", "samtools", "seqkit", "fastahack", "seqtk", "pyfaidx"]


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
        'size_10mb_100bp'    → '10 MB / 100 bp'
        'size_100mb_10kbp'   → '100 MB / 10 kbp'
        '10mb_full'          → '10 MB (full)'
        '100bp'              → '100 bp'
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
    return name


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


def discover_latest(results_dir: Path, prefix: str) -> Path | None:
    dirs = sorted(results_dir.glob(f"{prefix}_*"), reverse=True)
    dirs = [d for d in dirs if d.is_dir()]
    if dirs:
        return dirs[0]
    csvs = sorted(results_dir.glob(f"{prefix}_*.csv"), reverse=True)
    if csvs:
        return csvs[0]
    return None


def load_dir_jsons(results_dir: Path, prefix: str) -> pd.DataFrame | None:
    d = discover_latest(results_dir, prefix)
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


def load_memory(results_dir: Path) -> pd.DataFrame | None:
    p = discover_latest(results_dir, "memory")
    if not p:
        return None
    df = pd.read_csv(p)
    # CSV columns: tool,time_s,mem_kb,major_faults,minor_faults,region_type
    if "mem_kb" in df.columns:
        df["mem_mb"] = df["mem_kb"] / 1024.0
    return df


def load_multi(results_dir: Path) -> pd.DataFrame | None:
    """Load multi-region benchmark results (bench [6]).

    Files come from results/multi_<timestamp>/<N>regions.json where N is
    the region count (1, 10, 50, 100).
    """
    d = discover_latest(results_dir, "multi")
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
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

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
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
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
    ax.grid(alpha=0.3)
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
                f"{val:.1f}× {sign}", va="center", fontweight="bold", fontsize=9)
    ax.axvline(x=1, color="gray", linestyle="--", alpha=0.6, linewidth=1.5)
    ax.set_xlabel(f"Speedup (z-fasta vs {baseline}; >1 = z-fasta wins)")
    ax.set_title(title)
    ax.set_xlim(0, max(sdf["speedup"].max() * 1.3, 2))
    ax.grid(axis="x", alpha=0.3)
    return _save(fig, out)


def fig_real_datasets(df: pd.DataFrame, out: Path) -> Path:
    """Real dataset benchmarks grouped by dataset+region, bars by tool."""
    df = _prep(df)
    df["bench_human"] = df["benchmark"].apply(_human_bench)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    benchmarks = list(df["bench_human"].unique())

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
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
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


def md_table(df: pd.DataFrame, group_col: str, human_names: bool = True) -> str:
    def fmt(row):
        return f"{row['mean']:.4f}s ±{row['stddev']:.4f}"
    df = df.copy()
    df["tool_base"] = df["tool"].apply(_strip_tool_prefix)
    df["cell"] = df.apply(fmt, axis=1)
    if human_names:
        df[group_col] = df[group_col].apply(_human_bench)
    pivot = df.pivot_table(index=group_col, columns="tool_base", values="cell", aggfunc="first")
    cols = sorted(pivot.columns, key=tool_sort_key)
    return pivot[cols].to_markdown()


def md_speedup_table(df: pd.DataFrame, group_col: str,
                     include_fastahack: bool = False) -> str:
    df = df.copy()
    df["tool_base"] = df["tool"].apply(_strip_tool_prefix)
    df["bench_human"] = df[group_col].apply(_human_bench)
    groups = sorted(df["bench_human"].unique())
    rows = []
    for g in groups:
        gdf = df[df["bench_human"] == g]
        sam = gdf[gdf["tool_base"] == "samtools"]
        zf = gdf[gdf["tool_base"] == "z-fasta"]
        if not len(sam) or not len(zf):
            continue
        sam_t = sam["mean"].values[0]
        zf_t = zf["mean"].values[0]
        row = {"Benchmark": g}
        if zf_t > 0:
            row["z-fasta vs samtools"] = f"**{sam_t / zf_t:.1f}×**"
        else:
            row["z-fasta vs samtools"] = "N/A"
        if include_fastahack:
            fh = gdf[gdf["tool_base"] == "fastahack"]
            if len(fh) and zf_t > 0:
                fh_t = fh["mean"].values[0]
                ratio = fh_t / zf_t
                if ratio >= 1:
                    row["z-fasta vs fastahack"] = f"**{ratio:.1f}×** (z-fasta wins)"
                else:
                    row["z-fasta vs fastahack"] = f"{ratio:.2f}× (fastahack faster)"
        rows.append(row)
    if not rows:
        return "_No comparison data available._"
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


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    script_dir = Path(__file__).resolve().parent
    results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else script_dir / "results"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    report = []
    report.append("# z-fasta GET Benchmark Report\n")
    report.append("_Auto-generated by `generate_report.py`_\n")

    # ── 1. Single-region latency ───────────────────────────────────
    single_df = load_dir_jsons(results_dir, "single")
    if single_df is not None and len(single_df):
        report.append("\n## Single-Region Extraction Latency\n")
        report.append("> All measurements on warm cache. "
                      "Each benchmark extracts a region from the **middle** of the target file.\n")
        report.append(md_table(single_df, "benchmark"))
        report.append("")
        fig_latency_panels(single_df, figures_dir / "single_latency.png")
        report.append("\n![Single-Region Latency](results/figures/single_latency.png)\n")

        report.append("> **Interpretation:** z-fasta achieves sub-millisecond latency via O(1) "
                      "index lookup + `pread`. samtools is ~2–2.6× slower due to full `.fai` "
                      "parse overhead on every call. seqkit carries a ~14 ms Go runtime startup "
                      "cost that dominates at small regions. fastahack matches samtools.\n")

        report.append("\n### Speedup vs samtools\n")
        report.append(md_speedup_table(single_df, "benchmark"))
        report.append("")
        p = fig_speedup_vs(single_df, "samtools", "z-fasta GET Speedup vs samtools",
                           figures_dir / "single_speedup.png")
        if p:
            report.append("\n![Single-Region Speedup](results/figures/single_speedup.png)\n")

    # ── 2. Full-sequence extraction ────────────────────────────────
    fullseq_df = load_dir_jsons(results_dir, "fullseq")
    if fullseq_df is not None and len(fullseq_df):
        report.append("\n## Full-Sequence Extraction\n")
        report.append(md_table(fullseq_df, "benchmark"))
        report.append("")
        fig_fullseq(fullseq_df, figures_dir / "fullseq.png")
        report.append("\n![Full-Sequence Extraction](results/figures/fullseq.png)\n")

        report.append("\n### Speedup vs samtools (and fastahack)\n")
        report.append(md_speedup_table(fullseq_df, "benchmark", include_fastahack=True))
        report.append("")
        report.append("> **Note:** fastahack is faster than z-fasta for large full-sequence "
                      "extraction (≥50 MB). fastahack uses a simpler unbuffered write path "
                      "while z-fasta's buffered I/O adds overhead proportional to sequence "
                      "length. This is an optimization target for a future release.\n")

    # ── 3. Region-size scaling ─────────────────────────────────────
    scale_df = load_dir_jsons(results_dir, "scale_region")
    if scale_df is not None and len(scale_df):
        report.append("\n## Region-Size Scaling\n")

        scale_df = _prep(scale_df)
        scale_df["region_bp"] = scale_df["benchmark"].str.extract(r"^(\d+)").astype(float)

        pivot = scale_df.pivot_table(
            index="region_bp", columns="tool_base", values="mean", aggfunc="first"
        )
        cols = sorted(pivot.columns, key=tool_sort_key)
        pivot = pivot[cols]
        tick_labels = []
        for v in pivot.index:
            if v >= 1_000_000:
                tick_labels.append(f"{int(v/1_000_000)} Mbp")
            elif v >= 1_000:
                tick_labels.append(f"{int(v/1_000)} kbp")
            else:
                tick_labels.append(f"{int(v)} bp")
        pivot.index = tick_labels
        pivot.index.name = "Region Size"
        report.append(pivot.to_markdown(floatfmt=".4f"))
        report.append("")
        fig_scaling_region(scale_df.copy(), figures_dir / "scaling_region.png")
        report.append("\n![Region-Size Scaling](results/figures/scaling_region.png)\n")
        report.append("> **Interpretation:** Latency is O(1) up to ~10 kbp. Process startup "
                      "and index lookup dominate below that threshold. Above ~100 kbp it grows linearly "
                      "with region size as the I/O read cost dominates. fastahack's output "
                      "advantage grows with region size due to its simpler write path.\n")

    # ── 4. Real datasets ──────────────────────────────────────────
    real_df = load_dir_jsons(results_dir, "real")
    if real_df is not None and len(real_df):
        report.append("\n## Real Dataset Extraction\n")
        report.append(md_table(real_df, "benchmark"))
        report.append("")
        fig_real_datasets(real_df, figures_dir / "real_get.png")
        report.append("\n![Real Dataset Extraction](results/figures/real_get.png)\n")
    else:
        report.append("\n## Real Dataset Extraction\n")
        report.append("> _No real-dataset results found. "
                      "Run `bench/shared/download_data.sh` then re-run benchmarks "
                      "without `--skip-real`._\n")

    # ── 5. Memory ────────────────────────────────────────────────
    mem_df = load_memory(results_dir)
    if mem_df is not None and len(mem_df):
        report.append("\n## Memory Usage\n")
        report.append("> RSS measured with `/usr/bin/time -v`. For mmap-based tools "
                      "(z-fasta, samtools, fastahack) RSS reflects file pages mapped by the OS, "
                      "not heap allocation. `Major Faults` should be ~0 on warm cache.\n")
        report.append(md_memory_table(mem_df))
        report.append("")

    # ── 6. Multi-region extraction ────────────────────────────────
    multi_df = load_multi(results_dir)
    if multi_df is not None and len(multi_df):
        report.append("\n## Multi-Region Extraction (v0.2.4)\n")
        report.append("> Time to extract N regions in a single CLI call, loading the index once "
                      "and streaming all results to stdout. z-fasta uses O(1) offset lookup per "
                      "region. samtools re-parses its FAI header per call but otherwise uses the "
                      "same mmap strategy. seqtk scan time dominates because it has no index.\n")

        # Build a tidy table: rows = region count, columns = tool, values = mean time
        multi_df2 = multi_df.copy()
        multi_df2["tool_base"] = multi_df2["tool"].apply(_strip_tool_prefix)
        multi_df2 = multi_df2.sort_values("region_count")

        pivot = multi_df2.pivot_table(
            index="region_count", columns="tool_base", values="mean", aggfunc="first"
        )
        cols = sorted(pivot.columns, key=tool_sort_key)
        pivot = pivot[cols]

        # Format as ms
        def fmt_ms(v):
            if pd.isna(v):
                return "N/A"
            return f"{v * 1000:.1f} ms"

        fmt_pivot = pivot.map(fmt_ms)
        fmt_pivot.index.name = "Regions"
        report.append(fmt_pivot.to_markdown())
        report.append("")

        # Speedup vs samtools table
        rows_spd = []
        for rc in sorted(multi_df2["region_count"].unique()):
            gdf = multi_df2[multi_df2["region_count"] == rc]
            sam = gdf[gdf["tool_base"] == "samtools"]
            zf = gdf[gdf["tool_base"] == "z-fasta"]
            if len(sam) and len(zf):
                sam_t = sam["mean"].values[0]
                zf_t = zf["mean"].values[0]
                speedup = sam_t / zf_t if zf_t > 0 else float("nan")
                rows_spd.append({
                    "Regions": rc,
                    "z-fasta (ms)": f"{zf_t * 1000:.1f}",
                    "samtools (ms)": f"{sam_t * 1000:.1f}",
                    "Speedup": f"**{speedup:.1f}x**" if not pd.isna(speedup) else "N/A",
                })
        if rows_spd:
            report.append("\n### Speedup vs samtools\n")
            report.append(pd.DataFrame(rows_spd).to_markdown(index=False))
            report.append("")
            report.append("> **Note:** seqtk does not use an index for region extraction. It "
                          "performs a full-file scan for every call regardless of region count, "
                          "so its time is roughly constant here and reflects file I/O cost, not "
                          "region count. It is listed for reference only.\n")

    # ── Write report ──────────────────────────────────────────────
    report_text = "\n".join(report) + "\n"
    report_path = script_dir / "REPORT.md"
    report_path.write_text(report_text)
    print(f"Report written to {report_path}")
    print(f"Figures in {figures_dir}")


if __name__ == "__main__":
    main()
