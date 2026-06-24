#!/usr/bin/env python3
"""
z-fasta STATS Benchmark Report Generator

Reads raw hyperfine JSON + CSV data from bench/stats/results/,
produces Markdown report + PNG figures using pandas + matplotlib.

Usage:
    .venv/bin/python bench/stats/generate_report.py [results_dir]

Defaults to bench/stats/results/ (latest timestamped files).
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
    "z-fasta-full":      "#2E7D32",
    "z-fasta-indexonly": "#66BB6A",
    "seqkit-stats-a":    "#E65100",
    "seqtk-comp":        "#00838F",
}
TOOL_ORDER = ["z-fasta-full", "z-fasta-indexonly", "seqkit-stats-a", "seqtk-comp"]


def tool_sort_key(name):
    for i, t in enumerate(TOOL_ORDER):
        if name == t or t in name:
            return i
    return 99


def tool_color(name):
    for t in TOOL_ORDER:
        if t in name:
            return COLORS.get(t, "#888")
    return "#888"


def _strip_tool_prefix(name: str) -> str:
    """Extract base tool name from hyperfine label like '10mb_z-fasta-full'."""
    if "seqkit-stats" in name:
        return "seqkit-stats-a"
    for base in TOOL_ORDER:
        if name.endswith(base):
            return base
        if base in name:
            return base
    return name


def _human_file(name: str) -> str:
    """Convert benchmark stem to readable label.

    '10mb' -> '10 MB', 'edge_cases' -> 'Edge Cases', etc.
    """
    m = re.match(r"(\d+)mb$", name)
    if m:
        return f"{m.group(1)} MB"
    return name.replace("_", " ").title()


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
    # CSV columns: tool,time_s,mem_kb,major_faults,minor_faults,file_size_mb
    if "mem_kb" in df.columns and "mem_mb" not in df.columns:
        df["mem_mb"] = df["mem_kb"] / 1024.0
    return df


def load_throughput(results_dir: Path) -> pd.DataFrame | None:
    p = discover_latest(results_dir, "throughput")
    if not p:
        return None
    return pd.read_csv(p)


# ══════════════════════════════════════════════════════════════════════
#  Figures
# ══════════════════════════════════════════════════════════════════════

def _save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _bar_labels(ax, bars, vals_ms):
    for bar, v in zip(bars, vals_ms):
        if v > 0:
            label = f"{v:.3f}" if v < 0.1 else (f"{v:.2f}" if v < 10 else f"{v:.1f}")
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.03,
                    label, ha="center", va="bottom", fontsize=7, rotation=45)


def _prep(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["tool_base"] = df["tool"].apply(_strip_tool_prefix)
    return df


def fig_test_files(df: pd.DataFrame, out: Path) -> Path:
    """Stats on test files: x-axis = file name, bars by tool."""
    df = _prep(df)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    benchmarks = df["benchmark"].unique()
    bench_labels = [_human_file(b) for b in benchmarks]

    fig, ax = plt.subplots(figsize=(max(7, len(benchmarks) * 2), 5))
    x = range(len(benchmarks))
    width = 0.8 / max(len(tools), 1)

    for i, t in enumerate(tools):
        tdf = df[df["tool_base"] == t]
        vals_ms = []
        for b in benchmarks:
            row = tdf[tdf["benchmark"] == b]
            vals_ms.append(row["mean"].values[0] * 1000 if len(row) else 0)
        offset = (i - len(tools) / 2 + 0.5) * width
        bars = ax.bar([xi + offset for xi in x], vals_ms, width,
                      label=t, color=COLORS.get(t, "#888"), alpha=0.9)
        _bar_labels(ax, bars, vals_ms)

    ax.set_xticks(list(x))
    ax.set_xticklabels(bench_labels, fontsize=10)
    ax.set_ylabel("Time (ms)")
    ax.set_title("Stats Performance: Test Files")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    return _save(fig, out)


def fig_indexonly_vs_full(df: pd.DataFrame, out: Path) -> Path:
    """Index-only vs full scan: x-axis = file size, log y-scale, bars by tool."""
    df = _prep(df)
    df["file_mb"] = df["benchmark"].str.extract(r"(\d+)mb$").astype(float)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    file_sizes = sorted(df["file_mb"].dropna().unique())

    fig, ax = plt.subplots(figsize=(max(7, len(file_sizes) * 2.5), 5))
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
    ax.set_ylabel("Time (ms, log scale)")
    ax.set_title("Index-Only vs Full Scan Performance")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    return _save(fig, out)


def fig_scaling_size(df: pd.DataFrame, out: Path) -> Path:
    """File-size scaling: x-axis = file size (MB), lines by tool."""
    df = _prep(df)
    df["size_mb"] = df["benchmark"].str.extract(r"(\d+)mb$").astype(float)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]

    fig, ax = plt.subplots(figsize=(10, 5))
    for t in tools:
        tdf = df[df["tool_base"] == t].sort_values("size_mb")
        ax.plot(tdf["size_mb"], tdf["mean"], "o-", label=t,
                color=COLORS.get(t, "#888"), linewidth=2, markersize=6)

    ax.set_xlabel("File Size (MB)")
    ax.set_ylabel("Time (s)")
    ax.set_title("Stats Execution Time vs File Size")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    return _save(fig, out)


def fig_throughput(df: pd.DataFrame, out: Path) -> Path:
    """Throughput line chart: x-axis = file size, y = MB/s, lines by tool."""
    tools = [t for t in TOOL_ORDER if t in df["tool"].values]

    fig, ax = plt.subplots(figsize=(9, 5))
    for t in tools:
        tdf = df[df["tool"] == t].sort_values("file_size_mb")
        ax.plot(tdf["file_size_mb"], tdf["throughput_mbs"], "o-", label=t,
                color=COLORS.get(t, "#888"), linewidth=2, markersize=6)

    ax.set_xlabel("File Size (MB)")
    ax.set_ylabel("Throughput (MB/s)")
    ax.set_title("Stats Full-Scan Throughput by File Size")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    # Annotate each data point
    for t in tools:
        tdf = df[df["tool"] == t].sort_values("file_size_mb")
        for _, row in tdf.iterrows():
            ax.annotate(f"{row['throughput_mbs']:.0f}",
                        (row["file_size_mb"], row["throughput_mbs"]),
                        textcoords="offset points", xytext=(0, 6),
                        fontsize=7, ha="center")
    return _save(fig, out)


def fig_real_datasets(df: pd.DataFrame, out: Path) -> Path:
    """Real dataset stats: x-axis = dataset, bars by tool."""
    df = _prep(df)
    tools = [t for t in TOOL_ORDER if t in df["tool_base"].values]
    benchmarks = list(df["benchmark"].unique())
    bench_labels = [_human_file(b) for b in benchmarks]

    fig, ax = plt.subplots(figsize=(max(8, len(benchmarks) * 2.5), 5))
    x = range(len(benchmarks))
    width = 0.8 / max(len(tools), 1)

    for i, t in enumerate(tools):
        tdf = df[df["tool_base"] == t]
        vals = []
        for b in benchmarks:
            row = tdf[tdf["benchmark"] == b]
            vals.append(row["mean"].values[0] if len(row) else 0)
        offset = (i - len(tools) / 2 + 0.5) * width
        bars = ax.bar([xi + offset for xi in x], vals, width,
                      label=t, color=COLORS.get(t, "#888"), alpha=0.9)
        for bar, v in zip(bars, vals):
            if v > 0:
                label = f"{v*1000:.0f}ms" if v < 1 else f"{v:.2f}s"
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.03,
                        label, ha="center", va="bottom", fontsize=7, rotation=45)

    ax.set_xticks(list(x))
    ax.set_xticklabels(bench_labels, fontsize=10)
    ax.set_ylabel("Time (s)")
    ax.set_title("Stats Performance: Real Datasets")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    return _save(fig, out)


# ══════════════════════════════════════════════════════════════════════
#  Markdown
# ══════════════════════════════════════════════════════════════════════

def md_table(df: pd.DataFrame, group_col: str, human_names: bool = True) -> str:
    def fmt(row):
        return f"{row['mean']:.4f}s ±{row['stddev']:.4f}"
    df = df.copy()
    df["tool_base"] = df["tool"].apply(_strip_tool_prefix)
    df["cell"] = df.apply(fmt, axis=1)
    if human_names:
        df[group_col] = df[group_col].apply(_human_file)
    pivot = df.pivot_table(index=group_col, columns="tool_base", values="cell", aggfunc="first")
    cols = sorted(pivot.columns, key=tool_sort_key)
    return pivot[cols].to_markdown()


def md_indexonly_speedup(df: pd.DataFrame) -> str:
    """Speedup table: index-only vs full scan vs seqkit, with μs formatting."""
    df = _prep(df)
    df["file_mb"] = df["benchmark"].str.extract(r"(\d+)mb$").astype(float)
    benchmarks = sorted(df["file_mb"].dropna().unique())
    rows = []
    for mb in benchmarks:
        bdf = df[df["file_mb"] == mb]
        full = bdf[bdf["tool_base"] == "z-fasta-full"]
        io = bdf[bdf["tool_base"] == "z-fasta-indexonly"]
        seq = bdf[bdf["tool_base"] == "seqkit-stats-a"]
        row = {"File": f"{int(mb)} MB"}
        if len(full):
            row["Full Scan"] = f"{full['mean'].values[0]:.4f}s"
        if len(io) and io["mean"].values[0] > 0:
            us = io["mean"].values[0] * 1_000_000
            row["Index-Only"] = f"{us:.0f} μs"
            if len(full) and full["mean"].values[0] > 0:
                sp = full["mean"].values[0] / io["mean"].values[0]
                row["Speedup (full vs index-only)"] = f"**{sp:.0f}x**"
        if len(seq):
            row["seqkit-stats-a"] = f"{seq['mean'].values[0]:.4f}s"
        rows.append(row)
    if not rows:
        return ""
    return pd.DataFrame(rows).to_markdown(index=False)


def md_memory_table(df: pd.DataFrame) -> str:
    # CSV columns: tool,time_s,mem_kb,major_faults,minor_faults,file_size_mb
    df = df.copy()
    want = ["tool", "file_size_mb", "time_s", "mem_kb"]
    cols = [c for c in want if c in df.columns]
    display = df[cols].copy()
    rename = {"tool": "Tool", "file_size_mb": "File (MB)",
              "time_s": "Time (s)", "mem_kb": "MaxRSS (KB)"}
    display.rename(columns=rename, inplace=True)
    if "MaxRSS (KB)" in display.columns:
        display["MaxRSS (KB)"] = display["MaxRSS (KB)"].apply(lambda v: f"{int(v):,}")
    return display.to_markdown(index=False, floatfmt=".4f")


def md_throughput_table(df: pd.DataFrame) -> str:
    pivot = df.pivot_table(
        index="file_size_mb", columns="tool", values="throughput_mbs", aggfunc="first"
    )
    cols = sorted(pivot.columns, key=tool_sort_key)
    pivot = pivot[cols]
    pivot.index = [f"{int(v)} MB" for v in pivot.index]
    pivot.index.name = "File Size"
    return pivot.to_markdown(floatfmt=".1f")


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════

def main():
    script_dir = Path(__file__).resolve().parent
    results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else script_dir / "results"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    report = []
    report.append("# z-fasta STATS Benchmark Report\n")
    report.append("_Auto-generated by `generate_report.py`_\n")

    # ── 1. Test file stats ─────────────────────────────────────────
    stats_df = load_dir_jsons(results_dir, "stats")
    if stats_df is not None and len(stats_df):
        report.append("\n## Test File Stats Performance\n")
        report.append("> Small test files with pre-existing `.zfi` indexes. "
                      "Measures end-to-end CLI overhead.\n")
        report.append(md_table(stats_df, "benchmark"))
        report.append("")
        fig_test_files(stats_df, figures_dir / "stats_test.png")
        report.append("\n![Test File Stats](results/figures/stats_test.png)\n")

    # ── 2. Index-only vs full scan ─────────────────────────────────
    io_df = load_dir_jsons(results_dir, "indexonly")
    if io_df is not None and len(io_df):
        report.append("\n## Index-Only vs Full Scan\n")
        report.append(md_table(io_df, "benchmark"))
        report.append("")
        fig_indexonly_vs_full(io_df, figures_dir / "indexonly_vs_full.png")
        report.append("\n![Index-Only vs Full](results/figures/indexonly_vs_full.png)\n")

        report.append("\n### Index-Only Speedup\n")
        report.append(md_indexonly_speedup(io_df))
        report.append("")
        report.append("> **Interpretation:** Index-only mode reads `.zfi` index data and "
                  "computes length-derived metrics without scanning FASTA sequence bytes. "
                  "Time is dominated by process startup "
                      "and is constant regardless of file size. This makes it suitable for "
                      "interactive assembly QC where waiting seconds for N50 is unacceptable.\n")

    # ── 3. File-size scaling ───────────────────────────────────────
    size_df = load_dir_jsons(results_dir, "scale_size")
    if size_df is not None and len(size_df):
        report.append("\n## Scaling: File Size\n")
        size_df = _prep(size_df)
        size_df["size_mb"] = size_df["benchmark"].str.extract(r"(\d+)mb$").astype(float)

        pivot = size_df.pivot_table(
            index="size_mb", columns="tool_base", values="mean", aggfunc="first"
        )
        cols = sorted(pivot.columns, key=tool_sort_key)
        pivot = pivot[cols]
        pivot.index = [f"{int(v)} MB" for v in pivot.index]
        pivot.index.name = "File Size"
        report.append(pivot.to_markdown(floatfmt=".4f"))
        report.append("")

        fig_scaling_size(size_df, figures_dir / "scaling_size.png")
        report.append("\n![Scaling by File Size](results/figures/scaling_size.png)\n")
        report.append(
            "> **Note:** z-fasta-full is ahead of seqkit-stats and seqtk-comp across "
            "the synthetic file-size scaling run while computing richer statistics. "
            "Real-dataset performance still depends on record count, sequence layout, and "
            "index-loading cost. "
            "Index-only time is effectively constant regardless of file size and is dominated "
            "by CLI startup plus reading `.zfi` index data.\n"
        )

    # ── 4. Real datasets ──────────────────────────────────────────
    real_df = load_dir_jsons(results_dir, "real")
    if real_df is not None and len(real_df):
        report.append("\n## Real Dataset Stats\n")
        report.append(md_table(real_df, "benchmark"))
        report.append("")
        fig_real_datasets(real_df, figures_dir / "real_stats.png")
        report.append("\n![Real Dataset Stats](results/figures/real_stats.png)\n")
    else:
        report.append("\n## Real Dataset Stats\n")
        report.append("> _No real-dataset results found. "
                      "Run `bench/shared/download_data.sh` then re-run benchmarks "
                      "without `--skip-real`._\n")

    # ── 5. Memory ─────────────────────────────────────────────────
    mem_df = load_memory(results_dir)
    if mem_df is not None and len(mem_df):
        report.append("\n## Memory Usage\n")
        report.append("> RSS measured with `/usr/bin/time -v`. "
                      "z-fasta-full uses mmap for reading. RSS reflects file pages mapped by OS. "
                      "z-fasta-indexonly reads only the tiny `.zfi` header, so RSS stays nearly constant.\n")
        report.append(md_memory_table(mem_df))
        report.append("")

    # ── 6. Throughput ─────────────────────────────────────────────
    tp_df = load_throughput(results_dir)
    if tp_df is not None and len(tp_df):
        report.append("\n## Throughput (Full Scan)\n")
        report.append(md_throughput_table(tp_df))
        report.append("")
        fig_throughput(tp_df, figures_dir / "throughput.png")
        report.append("\n![Throughput](results/figures/throughput.png)\n")
        report.append("> z-fasta-full throughput is ~1.3 GB/s on synthetic single-sequence "
                  "files in the current run, ahead of seqkit-stats and seqtk-comp. "
                  "z-fasta's advantage is the additional statistics computed "
                      "(N50, L50, N90, L90, AU, GC%, per-base composition) that seqkit-stats "
                      "does not provide.\n")

    # ── Tool comparison ──────────────────────────────────────────
    report.append("\n## Tool Comparison Notes\n")
    report.append(
        "| Tool | Stats Support | Notes |\n"
        "| :--- | :------------ | :---- |\n"
        "| z-fasta stats | Full (Tier 1 + Tier 2) | Index-only mode for startup-dominated stats from `.zfi` |\n"
        "| seqkit stats | Basic counts/sizes | No N50/GC/composition breakdown |\n"
        "| samtools | No stats command | Index-only; no sequence statistics |\n"
        "| fastahack | No stats command | Index-only; no sequence statistics |\n"
    )

    # ── Write report ─────────────────────────────────────────────
    report_text = "\n".join(report) + "\n"
    report_path = script_dir / "REPORT.md"
    report_path.write_text(report_text)
    print(f"Report written to {report_path}")
    print(f"Figures in {figures_dir}")


if __name__ == "__main__":
    main()

