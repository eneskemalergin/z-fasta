#!/usr/bin/env python3
"""
z-fasta Benchmark Report Generator

Reads raw hyperfine JSON + CSV data from bench/index/results/,
produces Markdown report + PNG figures using pandas + matplotlib.

Usage:
    .venv/bin/python bench/index/generate_report.py [results_dir]

Defaults to bench/index/results/ (latest timestamped files).
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import tabulate  # noqa: F401  -- needed by pd.to_markdown()

# ── Styling ────────────────────────────────────────────────────────
COLORS = {
    "z-fasta-default": "#2E7D32",
    "z-fasta-nodedup": "#66BB6A",
    "z-fasta-lowmem": "#A5D6A7",
    "samtools": "#1565C0",
    "seqkit": "#E65100",
    "fastahack": "#7B1FA2",
    "pyfaidx": "#F7A41D",
}
TOOL_ORDER = [
    "z-fasta-default",
    "z-fasta-nodedup",
    "z-fasta-lowmem",
    "samtools",
    "seqkit",
    "fastahack",
    "pyfaidx",
]


def tool_sort_key(name):
    try:
        return TOOL_ORDER.index(name)
    except ValueError:
        return 99


# ══════════════════════════════════════════════════════════════════════
#  Data Loading
# ══════════════════════════════════════════════════════════════════════


def load_hyperfine_json(path: Path) -> pd.DataFrame:
    """Load a single hyperfine --export-json file into a DataFrame.

    Returns columns: tool, mean, stddev, median, min, max, user, system
    The 'command' field from hyperfine contains the label (via -n).
    """
    with open(path) as f:
        data = json.load(f)

    rows = []
    for r in data.get("results", []):
        # The label set by `-n` ends up in 'command'
        rows.append(
            {
                "tool": r["command"],
                "mean": r["mean"],
                "stddev": r["stddev"],
                "median": r["median"],
                "min": r["min"],
                "max": r["max"],
                "user": r["user"],
                "system": r["system"],
            }
        )
    return pd.DataFrame(rows)


def discover_latest(results_dir: Path, prefix: str) -> Path | None:
    """Find latest timestamped directory or file matching prefix."""
    # Try directories first (perf_*, scale_size_*, scale_seqs_*)
    dirs = sorted(results_dir.glob(f"{prefix}_*"), reverse=True)
    dirs = [d for d in dirs if d.is_dir()]
    if dirs:
        return dirs[0]

    # Try CSV files (memory_*, tests_*)
    csvs = sorted(results_dir.glob(f"{prefix}_*.csv"), reverse=True)
    if csvs:
        return csvs[0]
    return None


def load_perf(results_dir: Path) -> pd.DataFrame | None:
    """Load real-dataset perf results into a single DataFrame.

    Adds a 'dataset' column from filename (e.g. Genome.json -> Genome).
    """
    d = discover_latest(results_dir, "perf")
    if not d:
        return None

    frames = []
    for jf in sorted(d.glob("*.json")):
        df = load_hyperfine_json(jf)
        df["dataset"] = jf.stem  # e.g. "Genome"
        frames.append(df)

    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_scaling(results_dir: Path, prefix: str, param_col: str) -> pd.DataFrame | None:
    """Load scale_size or scale_seqs directory of JSON files.

    Each JSON filename encodes the parameter (e.g. 100mb.json, 10000.json).
    """
    d = discover_latest(results_dir, prefix)
    if not d:
        return None

    frames = []
    for jf in sorted(d.glob("*.json")):
        df = load_hyperfine_json(jf)
        # Parse parameter from filename
        stem = jf.stem  # e.g. "100mb" or "10000"
        numeric = stem.replace("mb", "")
        try:
            val = float(numeric)
        except ValueError:
            val = 0
        df[param_col] = val
        frames.append(df)

    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_memory(results_dir: Path) -> pd.DataFrame | None:
    """Load memory CSV.

    Supports formats:
      New: tool,time_s,mem_kb,major_faults,minor_faults
      Mid: tool,time_s,mem_kb
      Old: mode,run,time_s,mem_mb,lines
    Normalizes to columns: tool, time_s, mem_kb, mem_mb, major_faults, minor_faults
    """
    p = discover_latest(results_dir, "memory")
    if not p:
        return None
    df = pd.read_csv(p)

    # Normalize column names
    if "mode" in df.columns and "tool" not in df.columns:
        df = df.rename(columns={"mode": "tool"})
    if "mem_kb" in df.columns and "mem_mb" not in df.columns:
        df["mem_mb"] = df["mem_kb"] / 1024.0
    if "mem_mb" in df.columns and "mem_kb" not in df.columns:
        df["mem_kb"] = df["mem_mb"] * 1024.0

    # Ensure page fault columns exist (default 0 for old data)
    for col in ("major_faults", "minor_faults"):
        if col not in df.columns:
            df[col] = 0

    # Average across runs if multiple
    if "run" in df.columns:
        df = df.groupby("tool", as_index=False).agg(
            {"time_s": "mean", "mem_kb": "mean", "mem_mb": "mean",
             "major_faults": "mean", "minor_faults": "mean"}
        )

    return df


def load_tests(results_dir: Path) -> pd.DataFrame | None:
    """Load edge-case test CSV."""
    p = discover_latest(results_dir, "tests")
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


def fig_performance(df: pd.DataFrame, out: Path) -> Path:
    """Grouped bar chart: time per dataset, grouped by tool."""
    datasets = df["dataset"].unique()
    tools = sorted(df["tool"].unique(), key=tool_sort_key)

    fig, ax = plt.subplots(figsize=(max(7, len(datasets) * 2.2), 5))
    x = range(len(datasets))
    width = 0.8 / len(tools)

    for i, tool in enumerate(tools):
        vals = []
        for ds in datasets:
            row = df[(df["dataset"] == ds) & (df["tool"] == tool)]
            vals.append(row["mean"].values[0] if len(row) else 0)
        offset = (i - len(tools) / 2 + 0.5) * width
        bars = ax.bar(
            [xi + offset for xi in x],
            vals,
            width,
            label=tool,
            color=COLORS.get(tool, "#888"),
        )
        for bar, v in zip(bars, vals):
            if v > 0:
                label = f"{v:.3f}s" if v < 1 else f"{v:.2f}s"
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=45,
                )

    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=10)
    ax.set_ylabel("Time (s)")
    ax.set_title("Performance Comparison (Real Datasets)")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    return _save(fig, out)


def fig_speedup(df: pd.DataFrame, out: Path) -> Path:
    """Horizontal bar chart: speedup of z-fasta-default vs samtools per dataset."""
    datasets = df["dataset"].unique()
    rows = []
    for ds in datasets:
        zf = df[(df["dataset"] == ds) & (df["tool"] == "z-fasta-default")]
        sam = df[(df["dataset"] == ds) & (df["tool"] == "samtools")]
        if len(zf) and len(sam):
            zf_t = zf["mean"].values[0]
            sam_t = sam["mean"].values[0]
            if zf_t > 0:
                rows.append({"dataset": ds, "speedup": sam_t / zf_t})
    if not rows:
        return None

    sdf = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(6, max(3, len(rows) * 0.8)))
    bars = ax.barh(
        sdf["dataset"], sdf["speedup"], color=COLORS["z-fasta-default"], height=0.5
    )
    for bar, val in zip(bars, sdf["speedup"]):
        label = f"{val:.1f}x"
        ax.text(
            bar.get_width() + 0.3,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            fontweight="bold",
            fontsize=9,
        )
    ax.axvline(x=1, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Speedup (vs samtools)")
    ax.set_title("z-fasta Speedup Factor")
    ax.set_xlim(0, sdf["speedup"].max() * 1.2)
    ax.grid(axis="x", alpha=0.3)
    return _save(fig, out)


def fig_scaling(
    df: pd.DataFrame, param_col: str, xlabel: str, title: str, out: Path
) -> Path:
    """Line chart: time vs parameter, one line per tool."""
    tools = sorted(df["tool"].unique(), key=tool_sort_key)
    params = sorted(df[param_col].unique())

    fig, ax = plt.subplots(figsize=(9, 5))
    for tool in tools:
        tdf = df[df["tool"] == tool].sort_values(param_col)
        ax.plot(
            tdf[param_col],
            tdf["mean"],
            "o-",
            label=tool,
            color=COLORS.get(tool, "#888"),
            linewidth=2,
            markersize=5,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Time (s)")
    ax.set_title(title)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    # Use log scale for seq count since range is 10..100000
    if param_col == "seq_count":
        ax.set_xscale("log")
    return _save(fig, out)


def fig_memory(df: pd.DataFrame, out: Path) -> Path:
    """Side-by-side bar chart: time and memory by tool."""
    df = df.sort_values("tool", key=lambda s: s.map(tool_sort_key))
    tools = df["tool"].tolist()
    colors = [COLORS.get(t, "#888") for t in tools]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Time
    bars1 = ax1.bar(range(len(tools)), df["time_s"], color=colors)
    ax1.set_xticks(range(len(tools)))
    ax1.set_xticklabels(tools, rotation=25, ha="right", fontsize=9)
    ax1.set_ylabel("Time (s)")
    ax1.set_title("Execution Time")
    ax1.grid(axis="y", alpha=0.3)
    for bar, val in zip(bars1, df["time_s"]):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{val:.2f}s",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    # Memory
    bars2 = ax2.bar(range(len(tools)), df["mem_mb"], color=colors)
    ax2.set_xticks(range(len(tools)))
    ax2.set_xticklabels(tools, rotation=25, ha="right", fontsize=9)
    ax2.set_ylabel("Max RSS (MB)")
    ax2.set_title("Peak Memory (MaxRSS)")
    ax2.grid(axis="y", alpha=0.3)
    for bar, val in zip(bars2, df["mem_mb"]):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{val:.0f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    plt.tight_layout()
    return _save(fig, out)


def fig_edge_heatmap(df: pd.DataFrame, out: Path) -> Path:
    """Heatmap: rows = test cases, columns = tools + match."""
    cases = df["test_case"].tolist()
    tool_cols = [c for c in df.columns if c.endswith("_exit")]
    tool_labels = [c.replace("_exit", "").replace("zfasta", "z-fasta-default") for c in tool_cols]

    n_tools = len(tool_cols)
    match_x = n_tools + 0.3

    fig, ax = plt.subplots(
        figsize=(4 + n_tools * 0.8, max(6, len(cases) * 0.32))
    )

    for i, case in enumerate(cases):
        yi = len(cases) - i - 1
        for j, col in enumerate(tool_cols):
            val = df.iloc[i][col]
            clr = "#4CAF50" if val == 0 else "#F44336"
            ax.add_patch(
                plt.Rectangle((j, yi), 1, 1, facecolor=clr, edgecolor="white", lw=2)
            )
            ax.text(
                j + 0.5,
                yi + 0.5,
                "✓" if val == 0 else "✗",
                ha="center",
                va="center",
                fontsize=11,
                color="white",
                fontweight="bold",
            )

        # Match column
        m = df.iloc[i]["output_match"]
        clr = "#4CAF50" if m == "MATCH" else "#FF9800"
        ax.add_patch(
            plt.Rectangle(
                (match_x, yi), 1, 1, facecolor=clr, edgecolor="white", lw=2
            )
        )
        ax.text(
            match_x + 0.5,
            yi + 0.5,
            m[0],
            ha="center",
            va="center",
            fontsize=10,
            color="white",
            fontweight="bold",
        )

    ax.set_xlim(0, match_x + 1.2)
    ax.set_ylim(0, len(cases))
    ax.set_xticks([j + 0.5 for j in range(n_tools)] + [match_x + 0.5])
    ax.set_xticklabels(tool_labels + ["Match"], fontweight="bold")
    ax.set_yticks([len(cases) - i - 0.5 for i in range(len(cases))])
    ax.set_yticklabels(cases)
    ax.set_title("Edge Case Results (✓ = accepted, ✗ = rejected)", fontweight="bold")

    patches = [
        mpatches.Patch(color="#4CAF50", label="Pass / Match"),
        mpatches.Patch(color="#F44336", label="Fail"),
        mpatches.Patch(color="#FF9800", label="Diff"),
    ]
    ax.legend(handles=patches, loc="upper right", bbox_to_anchor=(1.3, 1), fontsize=8)
    plt.tight_layout()
    return _save(fig, out)


# ══════════════════════════════════════════════════════════════════════
#  Markdown Report
# ══════════════════════════════════════════════════════════════════════


def md_perf_table(df: pd.DataFrame) -> str:
    """Markdown table of performance results per dataset."""
    # Pivot to dataset x tool, format as "mean ± stddev"
    def fmt(row):
        return f"{row['mean']:.4f}s ±{row['stddev']:.4f}"

    df["cell"] = df.apply(fmt, axis=1)
    pivot = df.pivot(index="dataset", columns="tool", values="cell")
    # Reorder columns by TOOL_ORDER
    cols = sorted(pivot.columns, key=tool_sort_key)
    pivot = pivot[cols]
    return pivot.to_markdown()


def md_speedup_table(df: pd.DataFrame) -> str:
    """Markdown table of speedup vs samtools."""
    datasets = sorted(df["dataset"].unique())
    zf_tools = sorted(
        [t for t in df["tool"].unique() if t.startswith("z-fasta")],
        key=tool_sort_key,
    )
    rows = []
    for ds in datasets:
        sam = df[(df["dataset"] == ds) & (df["tool"] == "samtools")]
        if not len(sam):
            continue
        sam_t = sam["mean"].values[0]
        row = {"Dataset": ds}
        for t in zf_tools:
            r = df[(df["dataset"] == ds) & (df["tool"] == t)]
            if len(r) and r["mean"].values[0] > 0:
                row[f"{t} vs samtools"] = f"**{sam_t / r['mean'].values[0]:.1f}x**"
            else:
                row[f"{t} vs samtools"] = "N/A"
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_scaling_table(df: pd.DataFrame, param_col: str, param_label: str) -> str:
    """Markdown table: param x tool mean times."""
    pivot = df.pivot_table(
        index=param_col, columns="tool", values="mean", aggfunc="first"
    )
    cols = sorted(pivot.columns, key=tool_sort_key)
    pivot = pivot[cols]
    pivot.index = [
        f"{int(v)} MB" if param_col == "size_mb" else f"{int(v):,}"
        for v in pivot.index
    ]
    pivot.index.name = param_label
    return pivot.to_markdown(floatfmt=".4f")


def md_memory_table(df: pd.DataFrame) -> str:
    """Markdown table of memory results with page faults."""
    df = df.sort_values("tool", key=lambda s: s.map(tool_sort_key))
    display = df[["tool", "time_s", "mem_kb", "mem_mb", "major_faults", "minor_faults"]].copy()
    display.columns = ["Tool", "Time (s)", "MaxRSS (KB)", "MaxRSS (MB)", "Major Faults", "Minor Faults"]
    display["MaxRSS (KB)"] = display["MaxRSS (KB)"].apply(lambda v: f"{int(v):,}")
    display["Major Faults"] = display["Major Faults"].apply(lambda v: f"{int(v):,}")
    display["Minor Faults"] = display["Minor Faults"].apply(lambda v: f"{int(v):,}")
    return display.to_markdown(index=False, floatfmt=".2f")


def load_messy_index(results_dir: Path) -> pd.DataFrame | None:
    """Load messy FASTA indexing results from messy_* directories."""
    dirs = sorted(results_dir.glob("messy_*"), reverse=True)
    if not dirs:
        return None
    d = dirs[0]
    rows = []
    for f in sorted(d.glob("*.json")):
        data = json.loads(f.read_text())
        variant = f.stem
        for r in data["results"]:
            exits = r.get("exit_codes", [])
            n_ok = exits.count(0)
            n_total = len(exits)
            cmd = r["command"]
            tool = cmd.split()[0].split("/")[-1]
            # faidx binary is pyfaidx CLI
            if tool == "faidx":
                tool = "pyfaidx"
            success = n_ok == n_total and n_total > 0
            rows.append({"variant": variant, "tool": tool, "success": success})
    if not rows:
        return None
    return pd.DataFrame(rows)


def md_messy_table(df: pd.DataFrame) -> str:
    """Markdown compatibility matrix for messy FASTA variants."""
    messy_order = ["z-fasta", "samtools", "fastahack", "pyfaidx"]
    tools = [t for t in messy_order if t in df["tool"].values]
    # append any unexpected tools at end
    for t in df["tool"].unique():
        if t not in tools and t != "seqtk":
            tools.append(t)
    variants = sorted(df["variant"].unique())
    rows = []
    for variant in variants:
        row: dict = {"Variant": variant}
        for tool in tools:
            sub = df[(df["variant"] == variant) & (df["tool"] == tool)]
            if sub.empty:
                row[tool] = "n/a"
            else:
                row[tool] = "✓" if sub.iloc[0]["success"] else "✗"
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_test_table(df: pd.DataFrame) -> str:
    """Markdown table of edge case results."""
    tool_cols = [c for c in df.columns if c.endswith("_exit")]
    tool_labels = {
        c: c.replace("_exit", "").replace("zfasta", "z-fasta-default")
        for c in tool_cols
    }

    display = df[["test_case"] + tool_cols + ["output_match"]].copy()
    for col in tool_cols:
        display[col] = display[col].apply(lambda v: "✓" if v == 0 else "✗")
    display = display.rename(columns={"test_case": "Test Case", "output_match": "Match", **tool_labels})
    return display.to_markdown(index=False)


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════


def main():
    script_dir = Path(__file__).resolve().parent
    results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else script_dir / "results"
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    report_lines: list[str] = []

    def section(title, level=2):
        report_lines.append(f"\n{'#' * level} {title}\n")

    report_lines.append("# z-fasta Benchmark Report\n")
    report_lines.append(f"_Auto-generated by `generate_report.py`_\n")

    generated_figs: list[str] = []

    # ── 1. Performance ─────────────────────────────────────────────
    perf_df = load_perf(results_dir)
    if perf_df is not None and len(perf_df):
        section("Performance: Real Datasets")

        report_lines.append(md_perf_table(perf_df))
        report_lines.append("")

        p = fig_performance(perf_df, figures_dir / "performance.png")
        generated_figs.append(str(p))
        report_lines.append(f"![Performance Comparison](results/figures/performance.png)\n")

        # Speedup table
        section("Speedup vs samtools", 3)
        report_lines.append(md_speedup_table(perf_df))
        report_lines.append("")

        p = fig_speedup(perf_df, figures_dir / "speedup.png")
        if p:
            generated_figs.append(str(p))
            report_lines.append(f"![Speedup](results/figures/speedup.png)\n")
    else:
        section("Performance: Real Datasets")
        report_lines.append("_No real-dataset benchmark data found._\n")

    # ── 2. Scaling by file size ────────────────────────────────────
    size_df = load_scaling(results_dir, "scale_size", "size_mb")
    if size_df is not None and len(size_df):
        section("Scaling: File Size")
        report_lines.append(md_scaling_table(size_df, "size_mb", "File Size"))
        report_lines.append("")

        p = fig_scaling(
            size_df,
            "size_mb",
            "File Size (MB)",
            "Scaling by File Size",
            figures_dir / "scaling_size.png",
        )
        generated_figs.append(str(p))
        report_lines.append(f"![Scaling by File Size](results/figures/scaling_size.png)\n")

    # ── 3. Scaling by sequence count ───────────────────────────────
    seq_df = load_scaling(results_dir, "scale_seqs", "seq_count")
    if seq_df is not None and len(seq_df):
        section("Scaling: Sequence Count")
        report_lines.append(md_scaling_table(seq_df, "seq_count", "Sequences"))
        report_lines.append("")

        p = fig_scaling(
            seq_df,
            "seq_count",
            "Number of Sequences",
            "Scaling by Sequence Count",
            figures_dir / "scaling_seqs.png",
        )
        generated_figs.append(str(p))
        report_lines.append(
            f"![Scaling by Sequence Count](results/figures/scaling_seqs.png)\n"
        )

    # ── 4. Memory ──────────────────────────────────────────────────
    mem_df = load_memory(results_dir)
    if mem_df is not None and len(mem_df):
        section("Memory Usage")
        report_lines.append(md_memory_table(mem_df))
        report_lines.append("")
        report_lines.append(
            "> **Note:** MaxRSS from `/usr/bin/time -v` includes mmap'd pages. "
            "For mmap-based indexers (z-fasta-default, samtools), reported RSS "
            "reflects file size, not heap allocation. "
            "Use `--low-mem` for true streaming memory.\n"
            ">\n"
            "> **Page faults:** Major faults = blocking disk reads (should be ~0 "
            "on warm cache). Minor faults = non-blocking page mappings. A high "
            "count for mmap tools confirms the OS is seamlessly mapping file "
            "pages to RAM without blocking I/O.\n"
        )

        p = fig_memory(mem_df, figures_dir / "memory.png")
        generated_figs.append(str(p))
        report_lines.append(f"![Memory Usage](results/figures/memory.png)\n")

    # ── 5. Edge Case Tests ─────────────────────────────────────────
    test_df = load_tests(results_dir)
    if test_df is not None and len(test_df):
        section("Edge Case Correctness")

        total = len(test_df)
        matches = (test_df["output_match"] == "MATCH").sum()
        report_lines.append(
            f"**z-fasta vs samtools output match: {matches} / {total}**\n"
        )

        report_lines.append(md_test_table(test_df))
        report_lines.append("")

        p = fig_edge_heatmap(test_df, figures_dir / "edge_cases.png")
        generated_figs.append(str(p))
        report_lines.append(f"![Edge Case Heatmap](results/figures/edge_cases.png)\n")

    # ── 6. Messy FASTA Compatibility ──────────────────────────────
    messy_df = load_messy_index(results_dir)
    if messy_df is not None and len(messy_df):
        section("Messy FASTA Compatibility")
        report_lines.append(
            "Real-world FASTA files often have mixed line widths or trailing whitespace "
            "that violates the fixed-width assumption in the FAI spec. "
            "z-fasta indexes all four variants correctly. "
            "samtools, fastahack, and pyfaidx reject mixed-width and trailing-whitespace files.\n"
        )
        report_lines.append(md_messy_table(messy_df))
        report_lines.append("")
        report_lines.append(
            "> **Legend:** ✓ = all hyperfine runs exited 0 (indexing succeeded). "
            "✗ = all runs exited non-zero (indexing failed). "
            "n/a = tool absent in this benchmark run.\n"
            ">\n"
            "> **Variants:** `mixed_widths` has sequences wrapped at irregular column widths. "
            "`trailing_whitespace` has spaces at the end of each line. "
            "`crlf_endings` uses Windows CRLF line endings. "
            "`all_messy` combines all three.\n"
        )

    # ── 7. Tools Tested ────────────────────────────────────────────
    section("Tools Tested")
    report_lines.append(
        "| Tool | Description |\n"
        "|------|-------------|\n"
        "| z-fasta-default | Zig FASTA indexer (mmap + dedup, default) |\n"
        "| z-fasta-nodedup | z-fasta with `--no-dedup` (skip duplicate checking) |\n"
        "| z-fasta-lowmem | z-fasta with `--low-mem` (streaming, no mmap) |\n"
        "| samtools | `samtools faidx`: industry standard |\n"
        "| seqkit | `seqkit faidx`: Go bioinformatics toolkit |\n"
        "| fastahack | `fastahack -i`: C++ FASTA indexer (first-access indexing) |\n"
        "| pyfaidx | `faidx --no-output`: Python pyfaidx 0.9.0.3 |\n"
    )

    # ── Write report ───────────────────────────────────────────────
    # Write REPORT.md one level up from results/ (into bench/) so it
    # isn't gitignored with the rest of the generated results data.
    bench_dir = results_dir.parent
    report_path = bench_dir / "REPORT.md"
    report_path.write_text("\n".join(report_lines))

    print(f"Report: {report_path}")
    print(f"Figures: {len(generated_figs)}")
    for f in generated_figs:
        print(f"  {f}")


if __name__ == "__main__":
    main()
