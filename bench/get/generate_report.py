#!/usr/bin/env python3
"""Generate GET benchmark REPORT.md and figures from zebrac JSON."""

from __future__ import annotations

import argparse
import importlib.util
import json
import statistics
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results"
PROJECT_ROOT = SCRIPT_DIR.parent.parent

DATASET_ORDER = ["Genome", "Transcriptome", "Proteome"]

# Chart/table order: z-fasta first, then Rust peers, then samtools/bedtools, then reference lanes.
TOOL_ORDER = [
    "z-fasta-default",
    "z-fasta-fai",
    "z-fasta-chunk-all",
    "z-fasta-chunk-1",
    "z-fasta-rc",
    "z-fasta-complement-only",
    "z-fasta-reverse-only",
    "noodles",
    "rustbio-custom-get",
    "samtools",
    "bedtools",
    "seqtk-reference",
    "fastahack",
    "fastahack-reference",
    "noodles-rc",
    "rustbio-custom-get-rc",
]

# Positional charts: indexed tools + reference seqtk + fastahack on full_seq.
POS_CHART_TOOLS = [
    "z-fasta-default",
    "noodles",
    "rustbio-custom-get",
    "samtools",
    "fastahack",
    "seqtk-reference",
]

POS_TABLE_TOOLS = POS_CHART_TOOLS

# Back-compat alias
POS_HEADLINE_TOOLS = POS_TABLE_TOOLS

REGION_LABELS = {
    "100bp": "100 bp",
    "1kbp_mid": "1 kbp",
    "10kbp": "10 kbp",
    "full_seq": "full seq",
}

REFERENCE_TOOLS = frozenset({"seqtk-reference", "fastahack-reference"})

# RC chart only: plain z-fasta is a no-transform baseline (hatch like seqtk)
# so it does not blend with the three transform modes.
RC_BASELINE_TOOLS = frozenset({"z-fasta-default"})

POS_EXCLUDE = frozenset(
    {
        "z-fasta-fai",
        "z-fasta-chunk-all",
        "z-fasta-chunk-1",
        "z-fasta-rc",
        "noodles-rc",
        "rustbio-custom-get-rc",
    }
)

MODE_POS_TOOLS = ["z-fasta-default", "z-fasta-fai"]

# BED batch: 1 kbp per row (see run.sh / verify bed fixtures).
BED_BASES_PER_ROW = 1000
BED_ROW_ORDER = [10, 100, 1000, 10000]

# Headline BED chart/table peers (chunk modes are mode-comparison only).
BED_HEADLINE_TOOLS = [
    "z-fasta-default",
    "noodles",
    "rustbio-custom-get",
    "samtools",
    "bedtools",
]
BED_MODE_TOOLS = ["z-fasta-default", "z-fasta-chunk-all", "z-fasta-chunk-1"]
BED_DEFAULT_CHUNK = 4096

# Multi-region headline peers (reference loops omitted from default timed runs).
MULTI_HEADLINE_TOOLS = [
    "z-fasta-default",
    "noodles",
    "rustbio-custom-get",
    "samtools",
]
MULTI_CHART_TOOLS = MULTI_HEADLINE_TOOLS
MULTI_TABLE_TOOLS = MULTI_CHART_TOOLS
MULTI_N_ORDER = [1, 10, 100, 1000]
MULTI_REPORT_FIGURES = frozenset(
    {
        "perf_multi_wall.png",
        "perf_multi_rss.png",
        "perf_multi_faults.png",
    }
)
BED_REPORT_FIGURES = frozenset(
    {
        "perf_bed_wall.png",
        "perf_bed_rss.png",
        "perf_bed_faults.png",
    }
)

# RC overhead: nucleotide REAL datasets only (Proteome excluded; see md_rc_section).
RC_DATASET_ORDER = ["Genome", "Transcriptome"]
RC_EXCLUDED_DATASETS = frozenset({"Proteome"})

# Plain + RC-capable peers on a fixed 1 kbp region per RC dataset.
RC_HEADLINE_TOOLS = [
    "z-fasta-default",
    "z-fasta-rc",
    "z-fasta-complement-only",
    "z-fasta-reverse-only",
    "noodles-rc",
    "rustbio-custom-get-rc",
    "seqtk-reference",
]

# Consistent ratio wording in tables, figure notes, and reading sections.
LABEL_PEER_RATIO = "peer / z-fasta"
LABEL_PEER_WALL_RATIO = "peer wall time ÷ z-fasta wall time"
LABEL_PEER_RSS_RATIO = "peer peak RSS ÷ z-fasta peak RSS"
LABEL_PEER_FAULTS_RATIO = "peer minor faults ÷ z-fasta minor faults"

# Horizontal gap between Genome / Transcriptome / Proteome facet panels.
GROUPED_BAR_WSPACE = 0.10
# Gap between RC / messy multi-facet columns.
FACET_WSPACE = 0.28
# BED grouped-bar bottom margin and legend gap (two-line x labels; never rotated).
BED_FIG_BOTTOM = 0.28
BED_FIG_LEGEND_GAP = 0.055

# Figures always written when the corresponding section has data.
CORE_REPORT_FIGURES = frozenset(
    {
        "perf_pos_wall.png",
        "perf_pos_rss.png",
        "perf_pos_faults.png",
        "mode_pos.png",
        "messy_ratio_headline.png",
        "messy_ratio_summary.png",
        "perf_rc.png",
    }
)

# Back-compat alias (tests / docs).
GET_REPORT_FIGURES = CORE_REPORT_FIGURES

REGION_ORDER = ["100bp", "1kbp_mid", "10kbp", "full_seq"]

MESSY_PERF_JSON = SCRIPT_DIR / "messy_perf.json"


def format_md_table_cell(value: object) -> object:
    """Render pipe-containing strings without breaking markdown table columns (MD056)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return value
    if not isinstance(value, str):
        return value
    if "|" not in value:
        return value
    return value.replace("|", "&#124;")


def df_to_markdown(df: pd.DataFrame, **kwargs) -> str:
    """Render a DataFrame as markdown with pipe-safe cell values."""
    if df.empty:
        return "_No data._"
    safe = df.copy()
    for col in safe.columns:
        safe[col] = safe[col].map(format_md_table_cell)
    if kwargs.get("index", True) is not False:
        if isinstance(safe.index, pd.MultiIndex):
            safe.index = pd.MultiIndex.from_tuples(
                [tuple(format_md_table_cell(part) for part in row) for row in safe.index],
                names=[
                    format_md_table_cell(name) if isinstance(name, str) else name
                    for name in safe.index.names
                ],
            )
        else:
            safe.index = safe.index.map(format_md_table_cell)
    return safe.to_markdown(**kwargs)


def load_messy_perf_config() -> dict:
    data = json.loads(MESSY_PERF_JSON.read_text())
    spans = data["spans"]
    span_order = list(spans.keys())
    headline = data.get("headline_span", span_order[0] if span_order else "")
    return {
        "variants": list(data["variants"]),
        "span_order": span_order,
        "headline_span": headline,
        "spans": spans,
    }


def messy_span_bases(region_full: str) -> int:
    coords = region_full.rsplit(":", 1)[-1]
    start, end = coords.split("-", 1)
    return int(end) - int(start) + 1


MESSY_LAYOUT_TOOLS = ["z-fasta-uniform", "z-fasta-messy"]

GET_COLORS = {
    "z-fasta-default": "#F7A41D",
    "z-fasta-uniform": "#9E9E9E",
    "z-fasta-messy": "#F7A41D",
    "z-fasta-fai": "#E65100",
    "z-fasta-chunk-all": "#FFF59D",
    "z-fasta-chunk-1": "#FBC02D",
    # RC transform modes: darker / mid / pale shades of brand gold (#F7A41D).
    "z-fasta-rc": "#C67A00",
    "z-fasta-complement-only": "#F5B041",
    "z-fasta-reverse-only": "#FFE082",
    "noodles": "#C45C26",
    "noodles-rc": "#C45C26",
    "rustbio-custom-get": "#8B3A2A",
    "rustbio-custom-get-rc": "#8B3A2A",
    "samtools": "#555555",
    "bedtools": "#2E7D32",
    "seqtk-reference": "#6A1B9A",
    "fastahack": "#F34B7D",
    "fastahack-reference": "#F34B7D",
}

GET_DISPLAY = {
    "z-fasta-default": "z-fasta",
    "z-fasta-uniform": "z-fasta (uniform)",
    "z-fasta-messy": "z-fasta (messy)",
    "z-fasta-fai": "z-fasta (.fai)",
    "z-fasta-chunk-all": "z-fasta --chunk-size -1",
    "z-fasta-chunk-1": "z-fasta --chunk-size 1",
    "z-fasta-rc": "z-fasta --rc",
    "z-fasta-complement-only": "z-fasta --complement-only",
    "z-fasta-reverse-only": "z-fasta --reverse-only",
    "noodles": "noodles",
    "noodles-rc": "noodles --rc",
    "rustbio-custom-get": "rust-bio",
    "rustbio-custom-get-rc": "rust-bio --rc",
    "samtools": "samtools",
    "bedtools": "bedtools",
    "seqtk-reference": "seqtk (ref)",
    "fastahack": "fastahack",
    "fastahack-reference": "fastahack (ref)",
}

GET_PEERS = [
    {
        "key": "samtools",
        "label": "samtools",
        "language": "C",
        "command": "samtools faidx region",
        "version_from": "`samtools --version`",
        "pin": None,
    },
    {
        "key": "bedtools",
        "label": "bedtools",
        "language": "C++",
        "command": "bedtools getfasta",
        "version_from": "`bedtools --version`",
        "pin": None,
    },
    {
        "key": "noodles",
        "label": "noodles",
        "language": "Rust",
        "command": "tools/noodles_wrapper get",
        "version_from": "noodles-fasta crate in `tools/noodles_wrapper/Cargo.toml`",
        "pin": "0.61",
    },
    {
        "key": "rustbio",
        "label": "rust-bio",
        "language": "Rust",
        "command": "tools/rustbio_wrapper get",
        "version_from": "bio crate in `tools/rustbio_wrapper/Cargo.toml`",
        "pin": "2.2",
    },
    {
        "key": "fastahack",
        "label": "fastahack",
        "language": "C++",
        "command": "fastahack region",
        "version_from": "directory pin `tools/fastahack-1.0.0/`",
        "pin": "1.0.0",
    },
    {
        "key": "seqtk",
        "label": "seqtk",
        "language": "C",
        "command": "seqtk subseq (reference loop)",
        "version_from": "`seqtk 2>&1`",
        "pin": None,
    },
]


def load_index_report():
    path = SCRIPT_DIR.parent / "index" / "generate_report.py"
    spec = importlib.util.spec_from_file_location("index_generate_report", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    mod.COLORS.update(GET_COLORS)
    mod.DISPLAY_NAMES.update(GET_DISPLAY)
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
        return False
    sections = manifest.get("sections") or {}
    for key in ("perf_pos", "perf_multi", "perf_bed", "perf_rc"):
        if key not in sections:
            return True
    for flag in ("skip_pos", "skip_multi", "skip_bed", "skip_rc"):
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


def load_zebrac_json(path: Path, metadata_df: pd.DataFrame | None) -> pd.DataFrame:
    data = json.loads(path.read_text())
    meta_rows = []
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
        mean_s = wall.get("mean", 0) / 1_000_000_000.0
        std_s = wall.get("std_dev", 0) / 1_000_000_000.0
        out_bases = meta.get("output_bases")
        throughput_mbp_s = None
        if out_bases and mean_s > 0:
            throughput_mbp_s = float(out_bases) / mean_s / 1_000_000.0
        workload = meta.get("workload", path.stem)
        dataset = meta.get("dataset")
        region = meta.get("region")
        if dataset is None and isinstance(workload, str) and "/" in workload:
            dataset, region = workload.split("/", 1)
        elif dataset is None and "_" in path.stem:
            parts = path.stem.split("_", 1)
            if len(parts) == 2:
                dataset, region = parts
        rows.append(
            {
                "tool": meta.get("tool", ""),
                "section": meta.get("section", ""),
                "workload": workload,
                "dataset": dataset,
                "region": region,
                "mean": mean_s,
                "stddev": std_s,
                "peak_rss_mb": peak.get("mean", 0) / (1024 * 1024),
                "peak_rss_stddev_mb": peak.get("std_dev", 0) / (1024 * 1024),
                "minor_faults": result.get("minor_faults", {}).get("mean", 0),
                "minor_faults_stddev": result.get("minor_faults", {}).get("std_dev", 0),
                "major_faults": result.get("major_faults", {}).get("mean", 0),
                "output_bases": out_bases,
                "throughput_mibs": throughput_mbp_s,
                "input_bytes": meta.get("input_bytes"),
            }
        )
    return pd.DataFrame(rows)


def get_dataset_sort_key(name: str) -> tuple:
    try:
        return DATASET_ORDER.index(name), str(name)
    except ValueError:
        return len(DATASET_ORDER), str(name)


def tools_in_run(df: pd.DataFrame, *, exclude: frozenset[str] = frozenset()) -> list[str]:
    present = set(df["tool"].astype(str))
    return [t for t in TOOL_ORDER if t in present and t not in exclude]


def speedup_peer_tools(tools: list[str]) -> list[str]:
    """Peers for speedup bar labels (exclude z-fasta baseline and reference lanes)."""
    return [t for t in tools if t != "z-fasta-default" and t not in REFERENCE_TOOLS]


def manifest_sample_count(manifest: dict | None) -> str:
    if not manifest:
        return "?"
    runs = manifest.get("runs")
    return str(runs) if runs is not None else "?"


def region_sort_key(value: str) -> tuple:
    try:
        return (REGION_ORDER.index(value), value)
    except ValueError:
        return (len(REGION_ORDER), value)


def parse_multi_n(text: str) -> int | None:
    if not isinstance(text, str):
        return None
    part = text.split("/", 1)[-1]
    if part.startswith("N="):
        try:
            return int(part[2:])
        except ValueError:
            return None
    return None


def multi_n_label(n: int) -> str:
    return f"N={n:,}"


def multi_ds_n_label(dataset: str, n: int) -> str:
    return f"{dataset} / {multi_n_label(n)}"


def multi_ds_n_sort_key(label: str) -> tuple:
    if " / " not in str(label):
        return (len(DATASET_ORDER), 0, str(label))
    ds, rest = str(label).split(" / ", 1)
    n = int(rest.replace("N=", "").replace(",", ""))
    return get_dataset_sort_key(ds), n


def n_values_for_dataset(work: pd.DataFrame, dataset: str) -> list[int]:
    allowed = set(MULTI_N_ORDER)
    vals = work.loc[work["dataset"] == dataset, "n"].dropna().unique()
    return [n for n in MULTI_N_ORDER if n in allowed and n in {int(v) for v in vals}]


def multi_chart_tools(df: pd.DataFrame) -> list[str]:
    present = set(df["tool"].astype(str))
    tools = [t for t in MULTI_HEADLINE_TOOLS if t in present]
    for t in REFERENCE_TOOLS:
        if t in present and t not in tools:
            tools.append(t)
    return tools


def multi_table_tools(df: pd.DataFrame) -> list[str]:
    return multi_chart_tools(df)


def md_multi_intro_blurb(work: pd.DataFrame, table_tools: list[str]) -> str:
    datasets = real_datasets_in_work(work)
    n_list = ", ".join(multi_n_label(n) for n in MULTI_N_ORDER)
    ref_note = (
        " Reference loops (seqtk, fastahack) are omitted from default timed runs."
        if not any(t in REFERENCE_TOOLS for t in work["tool"].astype(str))
        else " Reference loops run only when `GET_MULTI_REFERENCE=1`."
    )
    tool_tail = (
        ", seqtk (ref), fastahack (ref)."
        if any(t in REFERENCE_TOOLS for t in table_tools)
        else "."
    )
    genome_ns = n_values_for_dataset(work, "Genome") if "Genome" in datasets else []
    genome_note = ""
    if genome_ns:
        missing = [n for n in MULTI_N_ORDER if n not in genome_ns]
        if missing:
            genome_note = (
                f" Genome omits {', '.join(multi_n_label(n) for n in missing)} "
                "(fewer eligible 1 kbp regions in the fixture)."
            )
    return (
        "One `get` invocation fetches **N regions × 1 kbp each** (e.g. N=100 → 100 kbp ≈ "
        f"0.1 Mb output). Timed on {', '.join(datasets)} with **{n_list}** (log-spaced). "
        "Three figure panels (Genome, Transcriptome, Proteome); panel width follows each "
        "dataset's N count. At **N≥16**, z-fasta uses `lookup_full_map` (N=100 and N=1,000 "
        f"here; N=1 and N=10 use the lighter index path (see `getter.zig`)).{genome_note}"
        f"{ref_note} Chart/table tool order: z-fasta, noodles, rust-bio, samtools{tool_tail}"
    )


def md_multi_figure_reading_base(
    chart_tools: list[str],
    *,
    y_note: str,
    ratio_label: str,
    t_cmp: int,
    f_num: int,
) -> str:
    ref_lines = ""
    if any(t in REFERENCE_TOOLS for t in chart_tools):
        ref_lines = (
            "- **Legend order:** z-fasta, noodles, rust-bio, samtools, seqtk (ref), "
            "fastahack (ref).\n"
            "- **Hatched bars:** reference loops (one subprocess per region).\n"
        )
    else:
        ref_lines = "- **Legend order:** z-fasta, noodles, rust-bio, samtools.\n"
    return (
        f"**Reading Figure {f_num}**\n"
        "- Three panels: Genome, Transcriptome, Proteome (N on x-axis).\n"
        f"- {y_note}\n"
        f"{ref_lines}"
        f"- **Bar labels:** `1×` on z-fasta; other labels = {ratio_label}. "
        f"Details in Table {t_cmp}."
    )


def parse_bed_rows(workload: str) -> int | None:
    if not isinstance(workload, str):
        return None
    part = workload.split("/", 1)[-1]
    if part.startswith("rows="):
        try:
            return int(part[5:].split("/")[0])
        except ValueError:
            return None
    return None


def real_datasets_in_work(work: pd.DataFrame) -> list[str]:
    if "dataset" not in work.columns or work["dataset"].isna().all():
        return ["Transcriptome"]
    present = set(work["dataset"].dropna().astype(str))
    ordered = [d for d in DATASET_ORDER if d in present]
    return ordered or sorted(present)


def bed_output_mb(row_count: int) -> float:
    return row_count * BED_BASES_PER_ROW / 1_000_000.0


def fmt_bed_output_mb(row_count: int) -> str:
    mb = bed_output_mb(row_count)
    if mb >= 10 and abs(mb - round(mb)) < 0.05:
        return f"{int(round(mb))} Mb"
    if mb >= 1:
        return f"{mb:.1f} Mb"
    return f"{mb:.2f} Mb"


def bed_row_sort_key(value) -> tuple:
    if isinstance(value, (int, float)) and not pd.isna(value):
        return (0, int(value))
    text = str(value)
    n = parse_bed_rows(f"rows={text.split()[0].replace(',', '')}")
    if n is not None:
        return (0, n)
    return (1, text)


def bed_row_table_label(row_count: int) -> str:
    return f"{row_count:,} rows ({fmt_bed_output_mb(row_count)})"


def bed_row_axis_label(row_count: int) -> str:
    return f"{row_count:,} rows\n({fmt_bed_output_mb(row_count)})"


def sort_comparisons_by_tool_order(
    comparisons: pd.DataFrame,
    tools: list[str],
    *,
    group_sort_key=None,
) -> pd.DataFrame:
    """Sort comparison rows by group then chart tool order (not alphabetical tool id)."""
    if comparisons.empty:
        return comparisons
    peers = speedup_peer_tools(tools)
    tool_rank = {t: i for i, t in enumerate(peers)}
    out = comparisons.copy()
    if group_sort_key is not None:
        out["_group_sort"] = out["dataset"].map(group_sort_key)
        out["_tool_sort"] = out["tool"].map(tool_rank)
        out = out.sort_values(["_group_sort", "_tool_sort"]).drop(
            columns=["_group_sort", "_tool_sort"]
        )
    else:
        out["_tool_sort"] = out["tool"].map(tool_rank)
        out = out.sort_values(["dataset", "_tool_sort"]).drop(columns="_tool_sort")
    return out


def bed_work_sorted(work: pd.DataFrame) -> pd.DataFrame:
    out = work.copy()
    out["bed_rows"] = out["rows"].astype(int)
    out = out.sort_values("bed_rows")
    return out


def bed_tools_present(work: pd.DataFrame, tool_list: list[str]) -> list[str]:
    present = set(work["tool"].astype(str))
    return [t for t in tool_list if t in present]


def bed_display_tool(tool: str, ir, *, mode_table: bool = False) -> str:
    if mode_table and tool == "z-fasta-default":
        return f"z-fasta (--chunk-size {BED_DEFAULT_CHUNK})"
    return ir.display_tool(tool)


def bed_ds_row_label(dataset: str, row_count: int) -> str:
    return f"{dataset} / {bed_row_table_label(row_count)}"


def bed_ds_row_sort_key(label: str) -> tuple:
    ds, rest = str(label).split(" / ", 1)
    row_part = rest.split(" rows", 1)[0].replace(",", "")
    return get_dataset_sort_key(ds), int(row_part)


def bed_rows_for_dataset(work: pd.DataFrame, dataset: str) -> list[int]:
    vals = work.loc[work["dataset"] == dataset, "bed_rows"].dropna().unique()
    present = {int(v) for v in vals}
    return [n for n in BED_ROW_ORDER if n in present]


def bed_chart_tools(df: pd.DataFrame) -> list[str]:
    present = set(df["tool"].astype(str))
    return [t for t in BED_HEADLINE_TOOLS if t in present]


def md_bed_intro_blurb(work: pd.DataFrame, *, t_mode: int | None = None) -> str:
    datasets = real_datasets_in_work(work)
    rows_present = [
        n for n in BED_ROW_ORDER if n in {int(v) for v in work["bed_rows"].unique()}
    ]
    row_list = ", ".join(bed_row_table_label(n) for n in rows_present)
    mode_blurb = ""
    if t_mode is not None:
        mode_blurb = (
            f" Table {t_mode} compares z-fasta `--chunk-size` settings when those lanes ran; "
            f"headline chart uses **default `--chunk-size {BED_DEFAULT_CHUNK}` only**."
        )
    return (
        f"`z-fasta get --bed` on {', '.join(datasets)}. Each BED row fetches **1 kbp** "
        f"(sweep: {row_list}). Three figure panels (Genome, Transcriptome, Proteome); "
        f"grouped bars per row count. Headline z-fasta is plain `get --bed` with CLI default "
        f"**`--chunk-size {BED_DEFAULT_CHUNK}`** (internal batching; not `--chunk-size -1` "
        f"read-whole-file or `--chunk-size 1`).{mode_blurb} Chart/table tool order: z-fasta, "
        "noodles, rust-bio, samtools, bedtools."
    )


def md_bed_figure_reading_base(
    *,
    y_note: str,
    ratio_label: str,
    t_cmp: int,
    f_num: int,
    include_chunk_note: bool = False,
) -> str:
    chunk_line = ""
    if include_chunk_note:
        chunk_line = (
            f"- **Headline z-fasta:** default `--chunk-size {BED_DEFAULT_CHUNK}` only.\n"
        )
    return (
        f"**Reading Figure {f_num}**\n"
        "- Three panels: Genome, Transcriptome, Proteome (BED row count on x-axis).\n"
        f"- {y_note}\n"
        "- **Legend order:** z-fasta, noodles, rust-bio, samtools, bedtools.\n"
        f"{chunk_line}"
        f"- **Bar labels:** `1×` on z-fasta; other labels = {ratio_label}. "
        f"Details in Table {t_cmp}."
    )


def rc_tools_present(work: pd.DataFrame, tool_list: list[str]) -> list[str]:
    present = set(work["tool"].astype(str))
    return [t for t in tool_list if t in present]


def enrich_rc(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    legacy = out["workload"].astype(str) == "1kbp"
    out.loc[legacy, "dataset"] = "Transcriptome"
    out.loc[legacy, "region"] = "1kbp_mid"
    missing_ds = out["dataset"].isna() if "dataset" in out.columns else pd.Series(True, index=out.index)
    if missing_ds.any():
        out.loc[missing_ds, "dataset"] = (
            out.loc[missing_ds, "workload"].astype(str).str.split("/", n=1).str[0]
        )
    if "region" not in out.columns:
        out["region"] = None
    missing_region = out["region"].isna() | (out["region"].astype(str) == "")
    if missing_region.any():
        split = out.loc[missing_region, "workload"].astype(str).str.split("/", n=1, expand=True)
        if split.shape[1] == 2:
            out.loc[missing_region, "region"] = split[1]
    out["region"] = out["region"].fillna("1kbp_mid")
    return out


def rc_chart_work(df: pd.DataFrame) -> pd.DataFrame:
    work = enrich_rc(df)
    if work.empty or "dataset" not in work.columns:
        return work
    return work[~work["dataset"].astype(str).isin(RC_EXCLUDED_DATASETS)]


def rc_tools_for_dataset(work: pd.DataFrame, tool_list: list[str]) -> list[str]:
    present = set(work["tool"].astype(str))
    return [t for t in tool_list if t in present]


def rc_peer_tools(tools: list[str]) -> list[str]:
    """RC section peer labels: compare RC lanes to z-fasta --rc (not plain)."""
    return [
        t
        for t in tools
        if t
        not in {
            "z-fasta-default",
            "z-fasta-rc",
            "z-fasta-complement-only",
            "z-fasta-reverse-only",
        }
        and t not in REFERENCE_TOOLS
    ]


def rc_metric_ratio(
    work: pd.DataFrame,
    dataset: str,
    tool: str,
    baseline: str,
    value_col: str,
) -> float | None:
    sub = work[work["dataset"] == dataset]
    base = sub.loc[sub["tool"] == baseline, value_col]
    peer = sub.loc[sub["tool"] == tool, value_col]
    if base.empty or peer.empty or float(base.iloc[0]) <= 0:
        return None
    return float(peer.iloc[0]) / float(base.iloc[0])


def rc_zfasta_overhead_ratio(work: pd.DataFrame, dataset: str | None = None) -> float | None:
    sub = work
    if dataset is not None:
        sub = work[work["dataset"] == dataset]
    return rc_metric_ratio(sub, dataset or str(sub["dataset"].iloc[0]), "z-fasta-rc", "z-fasta-default", "mean")


def rc_annotation_pairs(tools: list[str]) -> list[tuple[str, str]]:
    """(tool, baseline) pairs for RC facet bar labels."""
    pairs: list[tuple[str, str]] = []
    for tool in tools:
        if tool == "z-fasta-default":
            continue
        if tool in {"z-fasta-rc", "z-fasta-complement-only", "z-fasta-reverse-only"}:
            pairs.append((tool, "z-fasta-default"))
        else:
            pairs.append((tool, "z-fasta-rc"))
    return pairs


def _rc_annotate_facet(
    ax,
    work: pd.DataFrame,
    tools: list[str],
    datasets: list[str],
    *,
    value_col: str,
    bar_tops: dict[tuple[str, str], tuple[float, float]],
    width: float,
    ir,
) -> None:
    """Rotated ratio badges on every comparable RC bar (all metric facets)."""
    plain_color = ir.COLORS.get("z-fasta-default", "#F7A41D")
    rc_color = ir.COLORS.get("z-fasta-rc", "#C67A00")

    for ds in datasets:
        ds_work = work[work["dataset"] == ds]
        plain_key = (ds, "z-fasta-default")
        if plain_key in bar_tops:
            plain_x, plain_y = bar_tops[plain_key]
            cluster_xs = [bar_tops[(ds, t)][0] for t in tools if (ds, t) in bar_tops]
            if cluster_xs:
                ax.hlines(
                    plain_y,
                    min(cluster_xs) - width * 0.65,
                    max(cluster_xs) + width * 0.65,
                    colors=plain_color,
                    linestyles=(0, (4, 3)),
                    linewidth=1.0,
                    alpha=0.45,
                    zorder=1,
                )
            ax.annotate(
                "1×",
                (plain_x, plain_y),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=7,
                fontweight="bold",
                color=plain_color,
                bbox=dict(
                    boxstyle="round,pad=0.18",
                    facecolor="white",
                    edgecolor=plain_color,
                    linewidth=0.9,
                ),
                zorder=5,
            )

        for tool, baseline in rc_annotation_pairs(tools):
            key = (ds, tool)
            if key not in bar_tops:
                continue
            ratio = rc_metric_ratio(ds_work, ds, tool, baseline, value_col)
            if ratio is None:
                continue
            xpos, y = bar_tops[key]
            edge = (
                rc_color
                if baseline == "z-fasta-default" and tool == "z-fasta-rc"
                else ir.COLORS.get(tool, "#666666")
            )
            ax.annotate(
                ir._format_speedup(ratio),
                (xpos, y),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=7,
                fontweight="bold",
                color="#111111",
                bbox=dict(
                    boxstyle="round,pad=0.18",
                    facecolor="white",
                    edgecolor=edge,
                    linewidth=0.9,
                    alpha=0.96,
                ),
                zorder=5,
            )


def parse_messy_workload(workload: str) -> tuple[str | None, str | None, str | None]:
    text = str(workload)
    if text.count("/") >= 2:
        variant, span_id, layout = text.rsplit("/", 2)
        if layout in {"uniform", "messy"}:
            return variant, span_id, layout
    if "/" in text:
        base, layout = text.rsplit("/", 1)
        if layout not in {"uniform", "messy"}:
            return None, None, None
    else:
        base, layout = text, "messy"
    if ":" not in base:
        return None, None, None
    variant, region = base.split(":", 1)
    return variant, region, layout


def enrich_messy(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    parsed = out["workload"].map(parse_messy_workload)
    out["variant"] = parsed.map(lambda p: p[0] if p else None)
    out["span_id"] = parsed.map(lambda p: p[1] if p else None)
    out["layout"] = parsed.map(lambda p: p[2] if p else None)
    cfg = load_messy_perf_config()
    out["region"] = out["span_id"]
    out["region_full"] = out["span_id"].map(
        lambda s: cfg["spans"].get(str(s), {}).get("region") if pd.notna(s) else None
    )
    return out.dropna(subset=["variant", "span_id", "layout"])


def filter_messy_perf_work(work: pd.DataFrame) -> pd.DataFrame:
    """Keep only variant/span pairs defined in messy_perf.json."""
    cfg = load_messy_perf_config()
    allowed_variants = set(cfg["variants"])
    allowed_spans = set(cfg["spans"])
    return work[
        work["variant"].astype(str).isin(allowed_variants)
        & work["span_id"].astype(str).isin(allowed_spans)
    ].copy()


def messy_variants_in_work(work: pd.DataFrame) -> list[str]:
    cfg = load_messy_perf_config()
    present = set(work["variant"].astype(str))
    return [v for v in cfg["variants"] if v in present]


def messy_span_ids_in_work(work: pd.DataFrame) -> list[str]:
    cfg = load_messy_perf_config()
    present = set(work["span_id"].astype(str))
    return [s for s in cfg["span_order"] if s in present]


def messy_layout_value(
    work: pd.DataFrame,
    variant: str,
    span_id: str,
    layout: str,
    metric: str,
) -> float | None:
    row = work[
        (work["variant"] == variant)
        & (work["span_id"] == span_id)
        & (work["layout"] == layout)
    ]
    if row.empty or metric not in row.columns:
        return None
    val = row[metric].iloc[0]
    if pd.isna(val):
        return None
    return float(val)


def messy_overhead_ratio(
    work: pd.DataFrame,
    variant: str,
    span_id: str,
    metric: str,
) -> float | None:
    uniform = messy_layout_value(work, variant, span_id, "uniform", metric)
    messy = messy_layout_value(work, variant, span_id, "messy", metric)
    if uniform is None or messy is None or uniform <= 0:
        return None
    return messy / uniform


def fmt_us(seconds: float) -> str:
    return f"{seconds * 1_000_000:.1f} µs"


def fmt_us_pm(mean_s: float, std_s: float) -> str:
    return f"{fmt_us(mean_s)} ±{std_s * 1_000_000:.1f}"


def messy_region_catalog_table() -> pd.DataFrame:
    cfg = load_messy_perf_config()
    rows = []
    headline = cfg["headline_span"]
    for variant in cfg["variants"]:
        for span_id in cfg["span_order"]:
            span = cfg["spans"][span_id]
            region_full = span["region"]
            rows.append(
                {
                    "variant": variant,
                    "span": span_id,
                    "region": region_full,
                    "bases": messy_span_bases(region_full),
                    "purpose": span["purpose"],
                    "headline": "yes" if span_id == headline else "",
                }
            )
    return pd.DataFrame(rows)


def _messy_ratio_ylim(ratios: list[float], *, min_span: float = 0.18) -> tuple[float, float]:
    """Ratio y-limits that always show parity (1.0) and avoid micro-zoom."""
    if not ratios:
        return 0.90, 1.10
    lo_data = min(ratios)
    hi_data = max(ratios)
    lo = min(lo_data, 1.0)
    hi = max(hi_data, 1.0)
    pad = max(0.04, (hi_data - lo_data) * 0.25)
    ymin = lo - pad
    ymax = hi + pad
    if ymax - ymin < min_span:
        ymax = max(hi_data + 0.04, 1.0 + min_span * 0.55)
        ymin = min(1.0 - min_span * 0.45, lo_data - 0.04)
    return max(0.5, ymin), min(2.5, ymax)


def messy_variant_span_ratios(
    work: pd.DataFrame,
    variant: str,
    metric: str,
    span_order: list[str],
) -> list[float]:
    ratios: list[float] = []
    for span_id in span_order:
        ratio = messy_overhead_ratio(work, variant, span_id, metric)
        if ratio is not None:
            ratios.append(ratio)
    return ratios


def fig_messy_ratio_headline(
    work: pd.DataFrame,
    out: Path,
    variants: list[str],
    ir,
) -> Path:
    """Headline span: messy/uniform ratio panels with parity at 1.0×."""
    cfg = load_messy_perf_config()
    headline = cfg["headline_span"]
    metrics = [
        ("mean", "Wall time ×"),
        ("peak_rss_mb", "Peak RSS ×"),
        ("minor_faults", "Minor faults ×"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5), squeeze=False)
    color = ir.COLORS.get("z-fasta-messy", "#F7A41D")

    for col, (metric, ylabel) in enumerate(metrics):
        ax = axes[0, col]
        ratios: list[float] = []
        labels: list[str] = []
        for variant in variants:
            ratio = messy_overhead_ratio(work, variant, headline, metric)
            if ratio is None:
                continue
            ratios.append(ratio)
            labels.append(variant.replace("_", "\n"))
        xpos = list(range(len(ratios)))
        bars = ax.bar(xpos, ratios, color=color, alpha=0.92, width=0.62, zorder=2)
        ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.1, zorder=1)
        for bar, ratio in zip(bars, ratios):
            ax.annotate(
                f"{ratio:.3f}×",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
        ymin, ymax = _messy_ratio_ylim(ratios)
        ax.set_ylim(ymin, ymax)
        ax.grid(axis="y", alpha=0.28)
        ax.set_axisbelow(True)

    fig.suptitle(
        f"Messy FASTA GET: headline span `{headline}` (messy / uniform)",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.93,
        "Proteome-derived fixtures (~20k seqs). Dashed line = parity (1.0×).",
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    fig.subplots_adjust(bottom=0.14, top=0.84, wspace=FACET_WSPACE)
    return ir._save(fig, out)


def fig_messy_ratio_summary(
    work: pd.DataFrame,
    out: Path,
    variants: list[str],
    ir,
) -> Path:
    """Layout-variant summary: median ratio over spans, min-max whiskers, headline marker."""
    cfg = load_messy_perf_config()
    headline = cfg["headline_span"]
    span_order = messy_span_ids_in_work(work)
    metrics = [
        ("mean", "Wall time ×"),
        ("peak_rss_mb", "Peak RSS ×"),
        ("minor_faults", "Minor faults ×"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.8), squeeze=False)
    color = ir.COLORS.get("z-fasta-messy", "#F7A41D")
    headline_color = "#1565C0"

    all_ratios: list[float] = []
    for variant in variants:
        for metric, _ in metrics:
            all_ratios.extend(messy_variant_span_ratios(work, variant, metric, span_order))
    shared_ylim = _messy_ratio_ylim(all_ratios)

    for col, (metric, ylabel) in enumerate(metrics):
        ax = axes[0, col]
        labels = [variant.replace("_", "\n") for variant in variants]
        xpos = list(range(len(variants)))
        medians: list[float] = []
        yerr_lo: list[float] = []
        yerr_hi: list[float] = []
        headline_vals: list[float] = []

        for variant in variants:
            ratios = messy_variant_span_ratios(work, variant, metric, span_order)
            if not ratios:
                medians.append(float("nan"))
                yerr_lo.append(0.0)
                yerr_hi.append(0.0)
                headline_vals.append(float("nan"))
                continue
            med = statistics.median(ratios)
            medians.append(med)
            yerr_lo.append(med - min(ratios))
            yerr_hi.append(max(ratios) - med)
            headline_ratio = messy_overhead_ratio(work, variant, headline, metric)
            headline_vals.append(headline_ratio if headline_ratio is not None else med)

        bars = ax.bar(xpos, medians, color=color, alpha=0.92, width=0.62, zorder=2)
        ax.errorbar(
            xpos,
            medians,
            yerr=[yerr_lo, yerr_hi],
            fmt="none",
            ecolor="#333333",
            elinewidth=1.0,
            capsize=4,
            capthick=1.0,
            zorder=3,
        )
        ax.plot(
            xpos,
            headline_vals,
            linestyle="none",
            marker="D",
            markersize=7,
            color=headline_color,
            markeredgecolor="white",
            markeredgewidth=0.8,
            zorder=4,
        )
        ax.axhline(1.0, color="#666666", linestyle="--", linewidth=1.1, zorder=1)
        for bar, med in zip(bars, medians):
            if pd.isna(med):
                continue
            ax.annotate(
                f"{med:.3f}×",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
        ax.set_ylim(shared_ylim)
        ax.grid(axis="y", alpha=0.28)
        ax.set_axisbelow(True)

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=color, alpha=0.92),
        plt.Line2D([0], [0], color="#333333", marker="|", linestyle="none", markersize=8),
        plt.Line2D(
            [0],
            [0],
            linestyle="none",
            marker="D",
            markersize=7,
            color=headline_color,
            markeredgecolor="white",
        ),
    ]
    fig.legend(
        handles,
        [
            "Median over spans",
            "Min-max across spans",
            f"Headline span `{headline}`",
        ],
        loc="lower center",
        ncol=3,
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
    )

    fig.suptitle(
        "Messy FASTA GET: layout overhead (median over spans)",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.93,
        "Overhead is driven by layout defect type, not extraction length (1k-40k spans agree). "
        "Dashed line = parity (1.0×).",
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    fig.subplots_adjust(bottom=0.20, top=0.84, wspace=FACET_WSPACE)
    return ir._save(fig, out)


def messy_detail_table(work: pd.DataFrame) -> pd.DataFrame:
    cfg = load_messy_perf_config()
    rows = []
    for variant in messy_variants_in_work(work):
        for span_id in cfg["span_order"]:
            if span_id not in set(work["span_id"].astype(str)):
                continue
            u_row = work[
                (work["variant"] == variant)
                & (work["span_id"] == span_id)
                & (work["layout"] == "uniform")
            ]
            m_row = work[
                (work["variant"] == variant)
                & (work["span_id"] == span_id)
                & (work["layout"] == "messy")
            ]
            if u_row.empty or m_row.empty:
                continue
            region_full = cfg["spans"][span_id]["region"]
            rows.append(
                {
                    "variant": variant,
                    "span": span_id,
                    "region": region_full,
                    "bases": messy_span_bases(region_full),
                    "uniform wall": fmt_us_pm(
                        float(u_row["mean"].iloc[0]), float(u_row["stddev"].iloc[0])
                    ),
                    "messy wall": fmt_us_pm(
                        float(m_row["mean"].iloc[0]), float(m_row["stddev"].iloc[0])
                    ),
                    "wall ×": (
                        f"{messy_overhead_ratio(work, variant, span_id, 'mean'):.3f}×"
                        if messy_overhead_ratio(work, variant, span_id, "mean") is not None
                        else "n/a"
                    ),
                    "RSS ×": (
                        f"{messy_overhead_ratio(work, variant, span_id, 'peak_rss_mb'):.3f}×"
                        if messy_overhead_ratio(work, variant, span_id, "peak_rss_mb") is not None
                        else "n/a"
                    ),
                    "faults ×": (
                        f"{messy_overhead_ratio(work, variant, span_id, 'minor_faults'):.3f}×"
                        if messy_overhead_ratio(work, variant, span_id, "minor_faults") is not None
                        else "n/a"
                    ),
                }
            )
    return pd.DataFrame(rows)


def enrich_multi(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "dataset" in out.columns:
        legacy = out["dataset"].isna() & out["workload"].astype(str).str.match(r"^N=\d+$", na=False)
        out.loc[legacy, "dataset"] = "Transcriptome"
    out["n"] = out.apply(
        lambda r: parse_multi_n(
            r["region"]
            if pd.notna(r.get("region")) and str(r.get("region", "")).startswith("N=")
            else r["workload"]
        ),
        axis=1,
    )
    missing_ds = out["dataset"].isna() if "dataset" in out.columns else pd.Series(True, index=out.index)
    if missing_ds.any():
        out.loc[missing_ds, "dataset"] = (
            out.loc[missing_ds, "workload"].astype(str).str.split("/", n=1).str[0]
        )
    out = out.dropna(subset=["n"])
    return out[out["n"].isin(MULTI_N_ORDER)]


def enrich_bed(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["rows"] = out.apply(
        lambda r: parse_bed_rows(
            r["region"]
            if pd.notna(r.get("region")) and str(r.get("region", "")).startswith("rows=")
            else r["workload"]
        ),
        axis=1,
    )
    out = out.dropna(subset=["rows"])
    out["rows"] = out["rows"].astype(int)
    out["bed_rows"] = out["rows"]
    if "dataset" in out.columns:
        legacy = out["dataset"].isna() & out["workload"].astype(str).str.match(r"^rows=", na=False)
        out.loc[legacy, "dataset"] = "Transcriptome"
    missing_ds = out["dataset"].isna() if "dataset" in out.columns else pd.Series(True, index=out.index)
    if missing_ds.any():
        out.loc[missing_ds, "dataset"] = (
            out.loc[missing_ds, "workload"].astype(str).str.split("/", n=1).str[0]
        )
    return out


def tool_legend_blurb(tools: list[str], ir) -> str:
    parts = []
    for tool in tools:
        color_word = {
            "z-fasta-default": "gold",
            "z-fasta-fai": "deep orange",
            "z-fasta-chunk-all": "pale gold",
            "noodles": "bronze",
            "rustbio-custom-get": "red-brown",
            "samtools": "grey",
            "bedtools": "green",
            "pyfaidx": "blue",
            "seqtk": "purple",
            "fastahack": "pink",
            "fastahack-reference": "pink",
        }.get(tool, "grey")
        parts.append(f"{color_word} = {ir.display_tool(tool)}")
    return "; ".join(parts)


def build_throughput_comparisons(
    work: pd.DataFrame, tools: list[str], ir, *, group_col: str = "dataset"
) -> pd.DataFrame:
    peers = speedup_peer_tools(tools)
    ratio = ir.build_ratio_comparisons(
        work, "mean", peer_tools=peers, group_col=group_col
    )
    return ir.build_time_throughput_comparisons(work, ratio)


def fig_throughput_dual_axis(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    comparisons: pd.DataFrame,
    title: str,
    fig_note: str,
    ir,
    *,
    throughput_ylabel: str = "Output Throughput (Mbp/s)",
    dataset_order: list[str] | None = None,
) -> Path:
    if dataset_order:
        datasets = [d for d in dataset_order if d in work["dataset"].unique()]
    else:
        datasets = sorted(work["dataset"].unique(), key=get_dataset_sort_key)
    tools = [t for t in tools if t in work["tool"].unique()]

    fig, ax1 = plt.subplots(figsize=(14, 7.2))
    ax2 = ax1.twinx()
    width = min(0.095, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for i, tool in enumerate(tools):
        color = ir.COLORS.get(tool, "#888888")
        xs: list[float] = []
        tputs: list[float] = []
        linestyle = "--" if tool in REFERENCE_TOOLS else "-"
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            mean_t = max(float(row["mean"].values[0]), 1e-6)
            std_t = max(float(row["stddev"].values[0]), 0.0)
            xpos = di + (i - len(tools) / 2 + 0.5) * width
            xs.append(xpos)
            ax1.bar(
                xpos,
                mean_t,
                width,
                color=color,
                alpha=0.88,
                yerr=std_t,
                capsize=2,
                error_kw={"elinewidth": 0.8, "capthick": 0.8, "ecolor": "#333333"},
                zorder=2,
            )
            bar_tops[(ds, tool)] = (xpos, mean_t)
            tput = row["throughput_mibs"].values[0]
            if pd.notna(tput) and float(tput) > 0:
                tputs.append(float(tput))
        if xs and tputs and len(xs) == len(tputs):
            ax2.plot(
                xs,
                tputs,
                color=color,
                marker="o",
                linestyle=linestyle,
                linewidth=1.8,
                markersize=6,
                markeredgecolor="white",
                markeredgewidth=0.7,
                label=ir.display_tool(tool),
                zorder=3,
            )

    if not comparisons.empty:
        ir._annotate_headline_comparisons(ax1, datasets, tools, width, bar_tops, comparisons)

    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets, fontsize=11)
    ax1.set_ylabel("Wall Time (s, log scale)", fontsize=10)
    ax2.set_ylabel(throughput_ylabel, fontsize=10)
    ax1.set_title(title, fontsize=12, fontweight="bold", pad=20)
    fig.text(
        0.3125,
        0.945,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    ax1.grid(axis="y", alpha=0.28, which="both")
    ax1.set_axisbelow(True)
    handles, _labels = ax2.get_legend_handles_labels()
    fig.legend(
        handles,
        _labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=len(tools),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.92)
    return ir._save(fig, out)


def bed_filter_work(df: pd.DataFrame) -> pd.DataFrame:
    work = enrich_bed(df)
    if work.empty:
        return work
    return work[~work["workload"].astype(str).str.contains("/stdin", na=False)]


def md_bed_zfasta_vs_metric_table(
    work: pd.DataFrame,
    tools: list[str],
    ir,
    *,
    value_col: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
) -> str:
    work = bed_work_sorted(work)
    headline_tools = bed_tools_present(work, tools)
    tagged = work[work["tool"].isin(headline_tools)].copy()
    tagged["ds_row"] = tagged.apply(
        lambda r: bed_ds_row_label(r["dataset"], int(r["bed_rows"])),
        axis=1,
    )
    comparisons = ir.build_ratio_comparisons(
        tagged,
        value_col,
        peer_tools=speedup_peer_tools(headline_tools),
        group_col="ds_row",
        group_sort=bed_ds_row_sort_key,
    )
    if comparisons.empty:
        return "_No z-fasta comparison data._"
    comparisons = sort_comparisons_by_tool_order(
        comparisons,
        headline_tools,
        group_sort_key=bed_ds_row_sort_key,
    )
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        group_label="Dataset / BED rows",
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
    )


def md_bed_detail_pivot(
    work: pd.DataFrame,
    tools: list[str],
    ir,
    *,
    formatter,
    mode_table: bool = False,
) -> str:
    rows = []
    filtered = ir.filter_tools(work, tools)
    tools = bed_tools_present(filtered, tools)
    for ds in real_datasets_in_work(filtered):
        ds_work = filtered[filtered["dataset"] == ds]
        for n in bed_rows_for_dataset(filtered, ds):
            for _, row in ds_work[ds_work["bed_rows"] == n].iterrows():
                formatted = formatter(row)
                if formatted is None:
                    continue
                rows.append(
                    {
                        "dataset_row": bed_ds_row_label(ds, n),
                        "tool": bed_display_tool(row["tool"], ir, mode_table=mode_table),
                        "value": formatted,
                    }
                )
    detail = pd.DataFrame(rows)
    if detail.empty:
        return "_No BED data._"
    detail = detail.pivot(index="dataset_row", columns="tool", values="value")
    detail = detail.reindex(index=sorted(detail.index, key=bed_ds_row_sort_key))
    detail.index.name = "Dataset / BED rows"
    if mode_table:
        tool_labels = [
            bed_display_tool(t, ir, mode_table=True)
            for t in tools
            if bed_display_tool(t, ir, mode_table=True) in detail.columns
        ]
    else:
        tool_labels = [
            ir.display_tool(t) for t in tools if ir.display_tool(t) in detail.columns
        ]
    return detail.reindex(columns=tool_labels).pipe(df_to_markdown)


def md_bed_zfasta_vs_table(work: pd.DataFrame, tools: list[str], ir) -> str:
    work = bed_work_sorted(work)
    headline_tools = bed_tools_present(work, tools)
    tagged = work[work["tool"].isin(headline_tools)].copy()
    tagged["ds_row"] = tagged.apply(
        lambda r: bed_ds_row_label(r["dataset"], int(r["bed_rows"])),
        axis=1,
    )
    comparisons = ir.build_ratio_comparisons(
        tagged,
        "mean",
        peer_tools=speedup_peer_tools(headline_tools),
        group_col="ds_row",
        group_sort=bed_ds_row_sort_key,
    )
    if comparisons.empty:
        return "_No z-fasta comparison data._"
    comparisons = sort_comparisons_by_tool_order(
        comparisons,
        headline_tools,
        group_sort_key=bed_ds_row_sort_key,
    )
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        group_label="Dataset / BED rows",
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="Speedup",
        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
    )


def fig_metric_bars(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    value_col: str,
    std_col: str | None,
    ylabel: str,
    title: str,
    ir,
    *,
    fig_note: str | None = None,
    after_bars=None,
) -> Path:
    filtered = ir.filter_tools(work, tools)
    comparisons = ir.build_ratio_comparisons(
        filtered, value_col, peer_tools=speedup_peer_tools(tools), group_col="dataset"
    )
    return ir._fig_metric_bars(
        filtered,
        out,
        tools=tools,
        comparisons=comparisons,
        value_col=value_col,
        std_col=std_col,
        ylabel=ylabel,
        title=title,
        fig_note=fig_note,
        after_bars=after_bars,
    )


def pos_chart_tools(df: pd.DataFrame) -> list[str]:
    present = set(df["tool"].astype(str))
    return [t for t in POS_CHART_TOOLS if t in present]


def pos_table_tools(df: pd.DataFrame) -> list[str]:
    present = set(df["tool"].astype(str))
    return [t for t in POS_TABLE_TOOLS if t in present]


def _pos_std_col(value_col: str) -> str | None:
    if value_col == "mean":
        return "stddev"
    if value_col == "peak_rss_mb":
        return "peak_rss_stddev_mb"
    if value_col == "minor_faults":
        return "minor_faults_stddev"
    return None


def regions_for_dataset(work: pd.DataFrame, dataset: str) -> list[str]:
    regs = sorted(
        work.loc[work["dataset"] == dataset, "region"].dropna().unique(),
        key=region_sort_key,
    )
    return [r for r in regs if r in REGION_ORDER]


def _grouped_bar_patches(tools: list[str], ir) -> list:
    patches = []
    for tool in tools:
        color = ir.COLORS.get(tool, "#888888")
        patch_kw: dict = {
            "facecolor": color,
            "edgecolor": color if tool in REFERENCE_TOOLS else "none",
            "label": ir.display_tool(tool),
            "alpha": 0.75 if tool in REFERENCE_TOOLS else 0.88,
        }
        if tool in REFERENCE_TOOLS:
            patch_kw["hatch"] = "///"
        patches.append(mpatches.Patch(**patch_kw))
    return patches


def _finish_grouped_bar_figure(
    fig,
    patches: list,
    tool_labels: list[str],
    *,
    rotate_xt: bool,
) -> None:
    """Shared legend/margins for positional and multi 3-panel grouped bars."""
    if rotate_xt:
        bottom = 0.20
        legend_gap = 0.045
    else:
        bottom = 0.15
        legend_gap = 0.045
    fig.subplots_adjust(left=0.035, right=0.995, bottom=bottom, top=0.86)
    fig.legend(
        handles=patches,
        labels=tool_labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, bottom - legend_gap),
        bbox_transform=fig.transFigure,
        ncol=min(len(tool_labels), 6),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )


def fig_pos_grouped_bars(
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
    annotate_comparisons: bool = False,
) -> Path:
    """Three side-by-side panels (Genome / Transcriptome / Proteome); grouped bars per region."""
    filtered = ir.filter_tools(work, tools)
    tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in DATASET_ORDER if d in filtered["dataset"].unique()]
    std_col = _pos_std_col(value_col)
    peer_tools = speedup_peer_tools(tools)
    region_counts = [len(regions_for_dataset(filtered, d)) for d in datasets]
    width_ratios = [max(1, n) for n in region_counts]
    inches_per_region = 1.25
    fig_w = inches_per_region * sum(width_ratios)

    fig, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(fig_w, 7.2),
        sharey=False,
        gridspec_kw={"width_ratios": width_ratios, "wspace": GROUPED_BAR_WSPACE},
    )
    if len(datasets) == 1:
        axes = [axes]

    ylab_set = False
    for ax, ds in zip(axes, datasets):
        ds_work = filtered[filtered["dataset"] == ds]
        regions = regions_for_dataset(filtered, ds)
        if not regions:
            ax.set_visible(False)
            continue
        if not ylab_set:
            ax.set_ylabel(ylabel, fontsize=10, labelpad=2)
            ylab_set = True
        present_tools = [t for t in tools if t in ds_work["tool"].unique()]
        width = min(0.11, 0.82 / max(1, len(present_tools)))
        x = list(range(len(regions)))
        bar_tops_by_reg: dict[tuple[str, str, str], tuple[float, float]] = {}

        for ti, tool in enumerate(present_tools):
            color = ir.COLORS.get(tool, "#888888")
            for ri, reg in enumerate(regions):
                row = ds_work[(ds_work["tool"] == tool) & (ds_work["region"] == reg)]
                if row.empty:
                    continue
                val = max(float(row[value_col].iloc[0]), value_floor)
                std = 0.0
                if std_col and std_col in row.columns:
                    raw = row[std_col].iloc[0]
                    if pd.notna(raw):
                        std = max(float(raw), 0.0)
                xpos = ri + (ti - len(present_tools) / 2 + 0.5) * width
                bar_kw: dict = {
                    "width": width,
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
                bar_tops_by_reg[(ds, reg, tool)] = (xpos, val)

        if annotate_comparisons and peer_tools:
            for reg in regions:
                reg_work = ds_work[ds_work["region"] == reg]
                reg_tops = {
                    (ds, t): bar_tops_by_reg[(ds, reg, t)]
                    for t in present_tools
                    if (ds, reg, t) in bar_tops_by_reg
                }
                if (ds, "z-fasta-default") not in reg_tops:
                    continue
                comparisons = ir.build_ratio_comparisons(
                    reg_work,
                    value_col,
                    peer_tools=peer_tools,
                )
                if not comparisons.empty:
                    ir._annotate_headline_comparisons(
                        ax, [ds], present_tools, width, reg_tops, comparisons
                    )

        if log_y:
            ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels([REGION_LABELS.get(r, r) for r in regions], fontsize=10)
        ax.set_title(ds, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.28, which="both")
        ax.set_axisbelow(True)

    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.97)
    fig.text(
        0.5,
        0.915,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    patches = _grouped_bar_patches(tools, ir)
    _finish_grouped_bar_figure(
        fig,
        patches,
        [ir.display_tool(t) for t in tools],
        rotate_xt=False,
    )
    return ir._save(fig, out)


def fig_multi_grouped_bars(
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
    annotate_comparisons: bool = False,
) -> Path:
    """Three side-by-side panels (Genome / Transcriptome / Proteome); grouped bars per N."""
    filtered = ir.filter_tools(work, tools)
    tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in DATASET_ORDER if d in filtered["dataset"].unique()]
    std_col = _pos_std_col(value_col)
    peer_tools = speedup_peer_tools(tools)
    n_counts = [len(n_values_for_dataset(filtered, d)) for d in datasets]
    width_ratios = [max(1, n) for n in n_counts]
    # ~1.25" per N cluster (positional uses 5.2" for 4 regions); panel width ∝ its N count.
    inches_per_n = 1.25
    fig_w = inches_per_n * sum(width_ratios)
    max_n = max(width_ratios)
    rotate_xt = max_n > 5

    fig, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(fig_w, 7.2),
        sharey=False,
        gridspec_kw={"width_ratios": width_ratios, "wspace": GROUPED_BAR_WSPACE},
    )
    if len(datasets) == 1:
        axes = [axes]

    ylab_set = False
    for ax, ds in zip(axes, datasets):
        ds_work = filtered[filtered["dataset"] == ds]
        n_vals = n_values_for_dataset(filtered, ds)
        if not n_vals:
            ax.set_visible(False)
            continue
        if not ylab_set:
            ax.set_ylabel(ylabel, fontsize=10, labelpad=2)
            ylab_set = True
        present_tools = [t for t in tools if t in ds_work["tool"].unique()]
        width = min(0.11, 0.82 / max(1, len(present_tools)))
        x = list(range(len(n_vals)))
        bar_tops_by_n: dict[tuple[str, int, str], tuple[float, float]] = {}

        for ti, tool in enumerate(present_tools):
            color = ir.COLORS.get(tool, "#888888")
            for ni, n in enumerate(n_vals):
                row = ds_work[(ds_work["tool"] == tool) & (ds_work["n"] == n)]
                if row.empty:
                    continue
                val = max(float(row[value_col].iloc[0]), value_floor)
                std = 0.0
                if std_col and std_col in row.columns:
                    raw = row[std_col].iloc[0]
                    if pd.notna(raw):
                        std = max(float(raw), 0.0)
                xpos = ni + (ti - len(present_tools) / 2 + 0.5) * width
                bar_kw: dict = {
                    "width": width,
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
                bar_tops_by_n[(ds, n, tool)] = (xpos, val)

        if annotate_comparisons and peer_tools:
            for n in n_vals:
                n_work = ds_work[ds_work["n"] == n]
                n_tops = {
                    (ds, t): bar_tops_by_n[(ds, n, t)]
                    for t in present_tools
                    if (ds, n, t) in bar_tops_by_n
                }
                if (ds, "z-fasta-default") not in n_tops:
                    continue
                comparisons = ir.build_ratio_comparisons(
                    n_work,
                    value_col,
                    peer_tools=peer_tools,
                )
                if not comparisons.empty:
                    ir._annotate_headline_comparisons(
                        ax, [ds], present_tools, width, n_tops, comparisons
                    )

        if log_y:
            ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [multi_n_label(n) for n in n_vals],
            fontsize=10 if not rotate_xt else 9,
            rotation=45 if rotate_xt else 0,
            ha="right" if rotate_xt else "center",
        )
        ax.set_title(ds, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.28, which="both")
        ax.set_axisbelow(True)

    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.97)
    fig.text(
        0.5,
        0.915,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    patches = _grouped_bar_patches(tools, ir)
    _finish_grouped_bar_figure(
        fig,
        patches,
        [ir.display_tool(t) for t in tools],
        rotate_xt=rotate_xt,
    )
    return ir._save(fig, out)


def fig_bed_grouped_bars(
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
    annotate_comparisons: bool = False,
) -> Path:
    """Three side-by-side panels (Genome / Transcriptome / Proteome); grouped bars per BED row count."""
    filtered = ir.filter_tools(work, tools)
    tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in DATASET_ORDER if d in filtered["dataset"].unique()]
    std_col = _pos_std_col(value_col)
    peer_tools = speedup_peer_tools(tools)
    row_counts_per_ds = [len(bed_rows_for_dataset(filtered, d)) for d in datasets]
    width_ratios = [max(1, n) for n in row_counts_per_ds]
    inches_per_row = 1.25
    fig_w = inches_per_row * sum(width_ratios)

    fig, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(fig_w, 7.2),
        sharey=False,
        gridspec_kw={"width_ratios": width_ratios, "wspace": GROUPED_BAR_WSPACE},
    )
    if len(datasets) == 1:
        axes = [axes]

    ylab_set = False
    for ax, ds in zip(axes, datasets):
        ds_work = filtered[filtered["dataset"] == ds]
        row_counts = bed_rows_for_dataset(filtered, ds)
        if not row_counts:
            ax.set_visible(False)
            continue
        if not ylab_set:
            ax.set_ylabel(ylabel, fontsize=10, labelpad=2)
            ylab_set = True
        present_tools = [t for t in tools if t in ds_work["tool"].unique()]
        width = min(0.11, 0.82 / max(1, len(present_tools)))
        x = list(range(len(row_counts)))
        bar_tops_by_row: dict[tuple[str, int, str], tuple[float, float]] = {}

        for ti, tool in enumerate(present_tools):
            color = ir.COLORS.get(tool, "#888888")
            for ri, n in enumerate(row_counts):
                row = ds_work[(ds_work["tool"] == tool) & (ds_work["bed_rows"] == n)]
                if row.empty:
                    continue
                val = max(float(row[value_col].iloc[0]), value_floor)
                std = 0.0
                if std_col and std_col in row.columns:
                    raw = row[std_col].iloc[0]
                    if pd.notna(raw):
                        std = max(float(raw), 0.0)
                xpos = ri + (ti - len(present_tools) / 2 + 0.5) * width
                bar_kw: dict = {
                    "width": width,
                    "color": color,
                    "alpha": 0.88,
                    "zorder": 2,
                }
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
                bar_tops_by_row[(ds, n, tool)] = (xpos, val)

        if annotate_comparisons and peer_tools:
            for n in row_counts:
                n_work = ds_work[ds_work["bed_rows"] == n]
                n_tops = {
                    (ds, t): bar_tops_by_row[(ds, n, t)]
                    for t in present_tools
                    if (ds, n, t) in bar_tops_by_row
                }
                if (ds, "z-fasta-default") not in n_tops:
                    continue
                comparisons = ir.build_ratio_comparisons(
                    n_work,
                    value_col,
                    peer_tools=peer_tools,
                )
                if not comparisons.empty:
                    ir._annotate_headline_comparisons(
                        ax, [ds], present_tools, width, n_tops, comparisons
                    )

        if log_y:
            ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [bed_row_axis_label(n) for n in row_counts],
            fontsize=10,
            ha="center",
        )
        ax.set_title(ds, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.28, which="both")
        ax.set_axisbelow(True)

    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.97)
    fig.text(
        0.5,
        0.915,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    patches = _grouped_bar_patches(tools, ir)
    fig.subplots_adjust(left=0.05, right=0.995, bottom=BED_FIG_BOTTOM, top=0.86)
    fig.legend(
        handles=patches,
        labels=[ir.display_tool(t) for t in tools],
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, BED_FIG_BOTTOM - BED_FIG_LEGEND_GAP),
        bbox_transform=fig.transFigure,
        ncol=min(len(tools), 6),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    return ir._save(fig, out)


def _rc_draw_facet(
    ax,
    work: pd.DataFrame,
    tools: list[str],
    datasets: list[str],
    *,
    value_col: str,
    std_col: str | None,
    ylabel: str,
    value_floor: float,
    log_y: bool,
    ir,
) -> None:
    """One RC facet: datasets on x-axis, colored bars per tool (mode_pos layout)."""
    width = min(0.095, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for ti, tool in enumerate(tools):
        color = ir.COLORS.get(tool, "#888888")
        is_styled = tool in REFERENCE_TOOLS or tool in RC_BASELINE_TOOLS
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            val = max(float(row[value_col].iloc[0]), value_floor)
            std = 0.0
            if std_col and std_col in row.columns and pd.notna(row[std_col].iloc[0]):
                std = max(float(row[std_col].iloc[0]), 0.0)
            xpos = di + (ti - len(tools) / 2 + 0.5) * width
            bar_kw: dict = {
                "width": width,
                "color": color,
                "alpha": 0.75 if is_styled else 0.88,
                "zorder": 2,
            }
            if is_styled:
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

    _rc_annotate_facet(
        ax,
        work,
        tools,
        datasets,
        value_col=value_col,
        bar_tops=bar_tops,
        width=width,
        ir=ir,
    )

    if log_y:
        ax.set_yscale("log")
        lo, hi = ax.get_ylim()
        if hi > lo > 0:
            ax.set_ylim(lo, hi * 2.0)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(axis="y", alpha=0.28, which="both")
    ax.set_axisbelow(True)


def fig_rc_overhead(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    ir,
    *,
    sample_n: int | str,
) -> Path:
    """RC overhead: datasets on x-axis; facets = wall time, peak RSS, page faults."""
    filtered = ir.filter_tools(rc_chart_work(work), tools)
    legend_tools = [t for t in tools if t in filtered["tool"].unique()]
    datasets = [d for d in RC_DATASET_ORDER if d in filtered["dataset"].unique()]
    facets = [
        ("mean", "stddev", "Wall Time (s)", True, 1e-6),
        ("peak_rss_mb", "peak_rss_stddev_mb", "Peak RSS (MB)", False, 0.1),
        ("minor_faults", "minor_faults_stddev", "Minor Page Faults", True, 1.0),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6.8), squeeze=False)
    for ax, (value_col, std_col, ylabel, log_y, floor) in zip(axes[0], facets):
        _rc_draw_facet(
            ax,
            filtered,
            legend_tools,
            datasets,
            value_col=value_col,
            std_col=std_col,
            ylabel=ylabel,
            value_floor=floor,
            log_y=log_y,
            ir=ir,
        )

    fig.suptitle(
        "RC GET: wall time, peak RSS, page faults (1 kbp mid slice)",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.93,
        (
            f"Grouped bars per dataset (Genome, Transcriptome). "
            f"Error bars = zebrac stddev (n={sample_n}). "
            "Hatched = plain z-fasta (no transform) and seqtk (ref). "
            "Rotated bar labels: z-fasta `--rc` / `--complement-only` / `--reverse-only` "
            "vs plain; RC peers and seqtk (ref) vs z-fasta `--rc` (all facets)."
        ),
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    patches = []
    for tool in legend_tools:
        color = ir.COLORS.get(tool, "#888888")
        is_styled = tool in REFERENCE_TOOLS or tool in RC_BASELINE_TOOLS
        patch_kw: dict = {
            "facecolor": color,
            "edgecolor": color if is_styled else "none",
            "label": ir.display_tool(tool),
            "alpha": 0.75 if is_styled else 0.88,
        }
        if is_styled:
            patch_kw["hatch"] = "///"
        patches.append(mpatches.Patch(**patch_kw))
    fig.legend(
        handles=patches,
        labels=[ir.display_tool(t) for t in legend_tools],
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.06),
        ncol=min(len(legend_tools), 4),
        frameon=False,
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.18, top=0.86, wspace=FACET_WSPACE)
    return ir._save(fig, out)


def fig_scaling_metric_lines(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    param_col: str,
    value_col: str,
    xlabel: str,
    ylabel: str,
    title: str,
    fig_note: str,
    ir,
    *,
    log_x: bool = True,
    log_y: bool = True,
    value_floor: float = 1e-6,
) -> Path:
    filtered = ir.filter_tools(work, tools)
    tools = [t for t in tools if t in filtered["tool"].unique()]
    fig, ax = plt.subplots(figsize=(14, 7.2))
    for tool in tools:
        tdf = filtered[filtered["tool"] == tool].sort_values(param_col)
        if tdf.empty:
            continue
        xs = tdf[param_col].astype(float)
        ys = tdf[value_col].astype(float).clip(lower=value_floor)
        std_col = _pos_std_col(value_col)
        yerr = None
        if std_col and std_col in tdf.columns:
            yerr = tdf[std_col].astype(float).clip(lower=0.0).values
        color = ir.COLORS.get(tool, "#888888")
        linewidth = 2.4 if tool == "z-fasta-default" else 1.8
        linestyle = "--" if tool in REFERENCE_TOOLS else "-"
        ax.errorbar(
            xs,
            ys,
            yerr=yerr,
            color=color,
            marker="o",
            linestyle=linestyle,
            linewidth=linewidth,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=ir.display_tool(tool),
            capsize=3,
            elinewidth=0.9,
            zorder=3 if tool == "z-fasta-default" else 2,
        )
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=20)
    fig.text(
        0.3125,
        0.945,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    ax.grid(alpha=0.28, which="both")
    ax.set_axisbelow(True)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=min(len(handles), 6),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.92)
    return ir._save(fig, out)


def fig_scaling_lines(
    work: pd.DataFrame,
    out: Path,
    tools: list[str],
    param_col: str,
    xlabel: str,
    title: str,
    fig_note: str,
    ir,
    *,
    label_fn=None,
) -> Path:
    filtered = ir.filter_tools(work, tools)
    tools = [t for t in tools if t in filtered["tool"].unique()]
    fig, ax = plt.subplots(figsize=(14, 7.2))
    for tool in tools:
        tdf = filtered[filtered["tool"] == tool].sort_values(param_col)
        if tdf.empty:
            continue
        xs = tdf[param_col].astype(float)
        ys = tdf["mean"].astype(float).clip(lower=1e-6)
        yerr = None
        if "stddev" in tdf.columns:
            yerr = tdf["stddev"].astype(float).clip(lower=0.0).values
        color = ir.COLORS.get(tool, "#888888")
        linewidth = 2.4 if tool == "z-fasta-default" else 1.8
        linestyle = "--" if tool in REFERENCE_TOOLS else "-"
        ax.errorbar(
            xs,
            ys,
            yerr=yerr,
            color=color,
            marker="o",
            linestyle=linestyle,
            linewidth=linewidth,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=ir.display_tool(tool),
            capsize=3,
            elinewidth=0.9,
            zorder=3 if tool == "z-fasta-default" else 2,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Wall Time (s, log scale)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=20)
    fig.text(
        0.3125,
        0.945,
        fig_note,
        ha="center",
        va="top",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    ax.grid(alpha=0.28, which="both")
    ax.set_axisbelow(True)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=len(tools),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.92)
    return ir._save(fig, out)


def md_wall_table(df: pd.DataFrame, tools: list[str], index_col: str, ir) -> str:
    return ir.md_tool_pivot_table(
        df,
        tools,
        index_col,
        "mean",
        lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
        sort_index=get_dataset_sort_key if index_col == "dataset" else None,
        empty_msg="_No performance data._",
    )


def fmt_mbp_s(value: float) -> str:
    if value >= 10:
        return f"{value:.1f} Mbp/s"
    if value >= 0.01:
        return f"{value:.2f} Mbp/s"
    return f"{value:.3f} Mbp/s"


def md_throughput_table(df: pd.DataFrame, tools: list[str], index_col: str, ir) -> str:
    if "throughput_mibs" not in df.columns:
        return "_No throughput data._"
    return ir.md_tool_pivot_table(
        df,
        tools,
        index_col,
        "throughput_mibs",
        lambda row: fmt_mbp_s(float(row["throughput_mibs"])),
        sort_index=get_dataset_sort_key if index_col == "dataset" else None,
        empty_msg="_No throughput data._",
        require_non_null=True,
    )


def md_zfasta_vs_table(work: pd.DataFrame, tools: list[str], ir, *, group_col: str = "dataset") -> str:
    comparisons = build_throughput_comparisons(work, tools, ir, group_col=group_col)
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="Speedup",
        fmt_zf=lambda r: f"{r.zfasta_s:.4f}s",
        fmt_comp=lambda r: f"{r.comp_s:.4f}s",
        extra_columns=lambda r: {
            "z-fasta Mbp/s": fmt_mbp_s(r.zfasta_mibs) if r.zfasta_mibs else "n/a",
            "Peer Mbp/s": fmt_mbp_s(r.comp_mibs) if r.comp_mibs else "n/a",
            "Throughput ×": ir._format_speedup(r.throughput_x),
        },
    )


def md_zfasta_vs_metric_table(work: pd.DataFrame, tools: list[str], value_col: str, ir, *, unit: str, ratio_label: str) -> str:
    filtered = ir.filter_tools(work, tools)
    comparisons = ir.build_ratio_comparisons(
        filtered, value_col, peer_tools=speedup_peer_tools(tools), group_col="dataset"
    )
    if value_col == "peak_rss_mb":
        fmt_zf = lambda r: f"{r.zfasta_v:.2f} MB"
        fmt_comp = lambda r: f"{r.comp_v:.2f} MB"
    else:
        fmt_zf = lambda r: f"{int(r.zfasta_v):,}"
        fmt_comp = lambda r: f"{int(r.comp_v):,}"
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
    )


def md_perf_dual_axis_section(
    work: pd.DataFrame,
    nums,
    tools: list[str],
    ir,
    *,
    intro: str,
    wall_caption: str,
    throughput_summary: str,
    comparison_summary: str,
    figure_path: str,
    figure_caption: str,
    reading_notes: str,
    figure_title: str,
    fig_note: str,
    out_figure: Path,
    dataset_order: list[str] | None = None,
) -> str:
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_cmp = nums.next_table()
    f_perf = nums.next_figure()
    fig_throughput_dual_axis(
        work,
        out_figure,
        tools,
        build_throughput_comparisons(work, tools, ir),
        figure_title,
        fig_note,
        ir,
        dataset_order=dataset_order,
    )
    blocks = [
        intro,
        f"**Table {t_wall}:** {wall_caption}",
        md_wall_table(work, tools, "dataset", ir),
        "<details>",
        f"<summary><strong>Table {t_tp}:</strong> {throughput_summary}</summary>",
        "",
        md_throughput_table(work, tools, "dataset", ir),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> {comparison_summary}</summary>",
        "",
        md_zfasta_vs_table(work, tools, ir),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_perf}:** {figure_caption}",
        f"![Figure {f_perf}: {figure_title}]({figure_path})",
        reading_notes.format(f_perf=f_perf, t_cmp=t_cmp, t_wall=t_wall, t_tp=t_tp),
    ]
    return "\n\n".join(blocks)


def md_memory_section(work: pd.DataFrame, nums, tools: list[str], ir, *, context: str, figure_path: str) -> str:
    t_rss = nums.next_table()
    t_cmp = nums.next_table()
    f_mem = nums.next_figure()
    return "\n\n".join(
        [
            context,
            (
                "### Peak RSS\n\n"
                "zebrac starts a new process for each sample and records peak RSS when it "
                "exits. That is the most RAM the process had in use at once (Linux "
                "`ru_maxrss` via `getrusage`). Table and figure show the mean across samples."
            ),
            f"**Table {t_rss}:** Peak RSS (MB, zebrac mean). Same tool order as Performance.",
            ir.md_tool_pivot_table(
                work,
                tools,
                "dataset",
                "peak_rss_mb",
                lambda row: f"{row['peak_rss_mb']:.2f} MB",
                sort_index=get_dataset_sort_key,
            ),
            "<details>",
            f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
            f"RSS × = peer peak RSS ÷ z-fasta peak RSS. Same ratios as bar labels on "
            f"Figure {f_mem}.</summary>",
            "",
            md_zfasta_vs_metric_table(work, tools, "peak_rss_mb", ir, unit="MB", ratio_label="RSS ×"),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            f"**Figure {f_mem}:** Table {t_rss} as grouped bars (log scale). Bar labels = "
            f"RSS × (see Table {t_cmp}). Gold dashed line = z-fasta RSS per dataset.",
            f"![Figure {f_mem}: peak RSS]({figure_path})",
            (
                f"**Reading Figure {f_mem}**\n"
                "- Bars: mean peak RSS. `1×` on z-fasta; other labels = peer / z-fasta.\n"
                "- Gold dashed line: z-fasta RSS for that dataset.\n"
                f"- Details in Table {t_cmp}."
            ),
        ]
    )


def md_page_faults_section(work: pd.DataFrame, nums, tools: list[str], ir, *, context: str, figure_path: str) -> str:
    t_minor = nums.next_table()
    t_cmp = nums.next_table()
    t_major = nums.next_table()
    f_pf = nums.next_figure()
    blocks = [
        context,
        (
            "### Page faults\n\n"
            "A **minor** fault maps a page without reading disk. A **major** fault reads "
            "from disk. zebrac reports both per run, like wall time."
        ),
        f"**Table {t_minor}:** Minor page faults (zebrac mean). Same tool order as Performance.",
        ir.md_tool_pivot_table(
            work,
            tools,
            "dataset",
            "minor_faults",
            lambda row: f"{int(row['minor_faults']):,}",
        ),
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
        f"Faults × = peer minor faults / z-fasta minor faults. Same ratios as bar "
        f"labels on Figure {f_pf}.</summary>",
        "",
        md_zfasta_vs_metric_table(work, tools, "minor_faults", ir, unit="", ratio_label="Faults ×"),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_major}:</strong> Major page faults.</summary>",
        "",
        ir.md_tool_pivot_table(
            work,
            tools,
            "dataset",
            "major_faults",
            lambda row: f"{int(row['major_faults']):,}",
        ),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_pf}:** Table {t_minor} as grouped bars (log scale). Bar labels = "
        f"Faults × (see Table {t_cmp}).",
        f"![Figure {f_pf}: page faults]({figure_path})",
        (
            f"**Reading Figure {f_pf}**\n"
            "- Bars: mean minor page faults. `1×` on z-fasta; other labels = peer / z-fasta.\n"
            f"- Details in Table {t_cmp}."
        ),
    ]
    return "\n\n".join(blocks)


def md_overview(manifest: dict | None) -> str:
    lines = [
        "This report compares `z-fasta get` to other FASTA region fetch tools on preloaded "
        "`.zfi` and `.fai` sidecars. Wall time, peak RSS, and page faults come from "
        "the same zebrac samples for every tool in each section.",
        "Other timed tools are **peers**: they perform the same GET work (named regions, "
        "multi-region batches, or BED batch fetch) on the same fixtures. Charts label "
        f"peer ratios as **{LABEL_PEER_RATIO}**. **Reference lanes** (seqtk positional, "
        "multi-region subprocess loops) use different invocation patterns and appear "
        "hatched when present.",
        "**Performance:** positional extraction (wall time, RSS, page faults vs region "
        "size), multi-region batching (N=1, 10, 100, 1,000; one get with N×1 kbp regions), "
        "and BED batch (`get --bed`; 10–10,000 rows × 1 kbp; three-panel grouped bars per "
        "metric) on each REAL dataset present "
        "in this run, RC overhead on nucleotide datasets, and z-fasta messy vs uniform "
        "GET pairs on verify fixtures.",
        "**Correctness:** `bench/get/verify.sh` (not timed here) checks byte-identical output "
        "against golden fixtures.",
        "**Tool order in charts:** z-fasta, noodles, rust-bio, samtools, fastahack (on "
        "`full seq`), then reference lanes (seqtk positional; multi seqtk/fastahack loops "
        "omitted from default timed runs). **bedtools** "
        "is timed in **Performance: BED batch** (`bedtools getfasta`); it is not a "
        "positional or multi-region peer because those workloads use region/name syntax, "
        "not BED files.",
        "**Multi-region reference loops** (seqtk, fastahack: one subprocess per region) "
        "are **omitted from default timed runs** because wall time grows linearly with N "
        "and largely duplicates positional seqtk reference data. Set `GET_MULTI_REFERENCE=1` "
        "in `run.sh` for a one-off reference sweep.",
        "**pyfaidx** is omitted from timed runs: Python per-call overhead makes it much "
        "slower on these fixtures, which would dominate charts and lengthen reruns. "
        "Enable it in `run.sh` if needed.",
        "**fastahack** positional timing uses `fastahack -r` on **`full seq`** rows "
        "(whole-entry fetch on bounded sequences; see fixture notes in `run.sh`).",
        "**seqtk (ref)** positional timing pipes one region via stdin on each call with "
        "no persistent index. It appears in charts with hatched bars. On large indexed "
        "FASTA files (e.g. Genome), seqtk rescans the whole file each call, so log-scale "
        "wall-time charts compress indexed tools near the bottom; compare indexed peers "
        "within each region, not bar height vs seqtk.",
    ]
    if manifest and not manifest.get("verify_skipped") and manifest.get("verify_pass"):
        lines.append(
            f"**Verify:** `bench/get/verify.sh` passed **{manifest['verify_pass']}** checks before this perf run."
        )
    if manifest and manifest.get("skip_multi"):
        lines.append("_This run skipped multi-region scaling (`skip_multi=true`)._")
    if manifest and manifest.get("skip_bed"):
        lines.append("_This run skipped BED batch (`skip_bed=true`)._")
    if manifest and manifest.get("skip_rc"):
        lines.append("_This run skipped RC overhead (`skip_rc=true`)._")
    if manifest and manifest.get("skip_messy"):
        lines.append("_This run skipped messy FASTA (`skip_messy=true`)._")
    elif manifest and manifest.get("messy_runs"):
        lines.append(
            f"_Messy GET used higher zebrac sampling ({manifest['messy_runs']} runs, "
            f"{manifest.get('messy_warmup', '?')} warmup, "
            f"{manifest.get('messy_duration_ms', '?')} ms) because timings are sub-millisecond._"
        )
    return "\n\n".join(lines)


def md_run_provenance(manifest: dict | None, project_root: Path, ir) -> str:
    if not manifest:
        return "_No run manifest found._"

    ts = manifest.get("timestamp", "unknown")
    runner = manifest.get("runner", "zebrac")
    mode = manifest.get("mode", "warm")
    runs = manifest.get("runs", "?")
    warmup = manifest.get("warmup", "?")
    duration = manifest.get("duration_ms", "?")
    zebrac = manifest.get("zebrac", "unknown")
    zfasta = manifest.get("z_fasta", "unknown")
    tools = manifest.get("tools") or {}
    sections = manifest.get("sections") or {}
    metadata = manifest.get("metadata") or f"metadata_{ts}.jsonl"

    lines = [
        f"Run **`{ts}`** used **{runner}** in **{mode}** mode: {runs} measured samples, "
        f"{warmup} warmup passes, {duration} ms minimum per sample.",
        f"- **Subject:** {zfasta} (Zig; default `.zfi` get, `.fai` lane, chunk modes in tables)",
        f"- **Runner:** {zebrac}",
    ]

    artifacts: list[str] = []
    for key, label in [
        ("perf_pos", "positional extraction"),
        ("perf_multi", "multi-region scaling"),
        ("perf_bed", "BED batch"),
        ("perf_rc", "RC overhead"),
        ("messy", "messy FASTA"),
    ]:
        if sections.get(key):
            artifacts.append(f"`results/{sections[key]}/` {label}")
    if artifacts:
        lines.append(
            f"- **Artifacts:** {'; '.join(artifacts)}; metadata in `results/{metadata}`"
        )
    else:
        lines.append(f"- **Artifacts:** metadata in `results/{metadata}`")
    if manifest.get("messy_runs") and not manifest.get("skip_messy"):
        lines.append(
            f"- **Messy GET zebrac:** {manifest['messy_runs']} samples, "
            f"{manifest.get('messy_warmup', '?')} warmup, "
            f"{manifest.get('messy_duration_ms', '?')} ms (headline sections use "
            f"{runs}/{warmup}/{duration} ms)."
        )

    lines.append("")
    lines.append(ir.md_data_used(project_root))

    measured = []
    for spec in GET_PEERS:
        ver = tools.get(spec["key"])
        if not ver:
            continue
        ver = str(ver).splitlines()[0]
        pin_text = f"; pin {spec['pin']}" if spec.get("pin") else ""
        measured.append(
            f"- **{spec['label']}** ({spec['language']}): {ver} via `{spec['command']}` "
            f"(from {spec['version_from']}{pin_text})"
        )
    if measured:
        lines.append(
            "**Peers in this run** (same GET task; versions captured at benchmark time; "
            "vendored pins in `bench/shared/tools.sh` / `tools/`):\n" + "\n".join(measured)
        )

    lines.append(
        "Indexes (`.zfi` / `.fai`) are built once in preload before timed GET commands. "
        "Timed commands do not delete sidecars."
    )
    return "\n\n".join(lines)


def md_pos_zfasta_vs_table(work: pd.DataFrame, tools: list[str], ir) -> str:
    """z-fasta vs each peer for every dataset/region row (includes reference lanes)."""
    peers = [t for t in tools if t != "z-fasta-default"]
    tagged = work.copy()
    tagged["ds_region"] = tagged.apply(
        lambda r: f"{r['dataset']} / {REGION_LABELS.get(r['region'], r['region'])}",
        axis=1,
    )
    comparisons = ir.build_ratio_comparisons(
        tagged,
        "mean",
        peer_tools=peers,
        group_col="ds_region",
        group_sort=None,
    )
    if comparisons.empty:
        return "_No z-fasta comparison data._"
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        group_label="Dataset / region",
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="Speedup",
        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
    )


def md_pos_zfasta_vs_metric_table(
    work: pd.DataFrame,
    tools: list[str],
    value_col: str,
    ir,
    *,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
) -> str:
    peers = [t for t in tools if t != "z-fasta-default"]
    tagged = work.copy()
    tagged["ds_region"] = tagged.apply(
        lambda r: f"{r['dataset']} / {REGION_LABELS.get(r['region'], r['region'])}",
        axis=1,
    )
    comparisons = ir.build_ratio_comparisons(
        tagged,
        value_col,
        peer_tools=peers,
        group_col="ds_region",
        group_sort=None,
    )
    if comparisons.empty:
        return "_No z-fasta comparison data._"
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        group_label="Dataset / region",
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
    )


def md_pos_detail_pivot(
    work: pd.DataFrame,
    tools: list[str],
    ir,
    *,
    value_col: str,
    formatter,
) -> str:
    rows = []
    for _, row in ir.filter_tools(work, tools).iterrows():
        formatted = formatter(row)
        if formatted is None:
            continue
        rows.append(
            {
                "dataset_region": (
                    f"{row['dataset']} / {REGION_LABELS.get(row['region'], row['region'])}"
                ),
                "tool": ir.display_tool(row["tool"]),
                "value": formatted,
            }
        )
    detail = pd.DataFrame(rows)
    if detail.empty:
        return "_No data._"
    detail = detail.pivot(index="dataset_region", columns="tool", values="value")
    tool_labels = [ir.display_tool(t) for t in tools if ir.display_tool(t) in detail.columns]
    return detail.reindex(columns=tool_labels).pipe(df_to_markdown)


def md_pos_wall_pivot_table(work: pd.DataFrame, tools: list[str], ir) -> str:
    return ir.md_tool_pivot_table(
        work,
        tools,
        "dataset",
        "mean",
        lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
        sort_index=get_dataset_sort_key,
        empty_msg="_No performance data._",
    )


def md_pos_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
    manifest: dict | None = None,
) -> str:
    chart_tools = pos_chart_tools(df)
    table_tools = pos_table_tools(df)
    if not chart_tools:
        return "_No positional chart data._"

    chart_work = df[df["tool"].isin(chart_tools)].copy()
    table_work = df[df["tool"].isin(table_tools)].copy()
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_cmp = nums.next_table()
    f_wall = nums.next_figure()
    sample_n = manifest_sample_count(manifest)

    fig_pos_grouped_bars(
        chart_work,
        figures_dir / "perf_pos_wall.png",
        chart_tools,
        "mean",
        "Wall Time (s)",
        "Positional GET: Wall Time vs Region Size",
        f"Error bars = zebrac stddev (n={sample_n}). Hatched bars = reference lanes (seqtk). "
        f"Bar labels = {LABEL_PEER_RATIO} (indexed peers).",
        ir,
        log_y=True,
        annotate_comparisons=True,
    )

    wall_md = md_pos_detail_pivot(
        table_work,
        table_tools,
        ir,
        value_col="mean",
        formatter=lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
    )
    throughput_md = md_pos_detail_pivot(
        table_work,
        table_tools,
        ir,
        value_col="throughput_mibs",
        formatter=lambda row: (
            fmt_mbp_s(float(row["throughput_mibs"]))
            if pd.notna(row.get("throughput_mibs"))
            else None
        ),
    )

    return "\n\n".join(
        [
            "Positional GET on preloaded indexes. Three panels (Genome, Transcriptome, "
            "Proteome) with grouped bars per region size. **`full seq`** is a "
            "whole-entry fetch on a bounded sequence (Genome ≤1 kbp; Transcriptome/Proteome "
            "medium entries, not titin). Tool order: z-fasta, noodles, rust-bio, samtools, "
            "fastahack (when present), seqtk (ref).",
            (
                f"**Table {t_wall}:** Wall time (mean ± stddev, seconds). Lower is better. "
                "Same tool order as the figure."
            ),
            wall_md,
            "<details>",
            (
                f"<summary><strong>Table {t_tp}:</strong> Output throughput (Mbp/s) from "
                f"metadata <code>output_bases</code> and wall time. Higher is better. Same "
                f"tools and order as Table {t_wall}.</summary>"
            ),
            "",
            throughput_md,
            "</details>",
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                "Speedup = peer wall time ÷ z-fasta wall time (higher = z-fasta faster). "
                f"Includes reference lanes. Same ratios as bar labels on "
                f"Figure {f_wall}.</summary>"
            ),
            "",
            md_pos_zfasta_vs_table(table_work, table_tools, ir),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_wall}:** Table {t_wall} as grouped bars (log y). Hatched bars = "
                f"seqtk (ref). Bar labels = {LABEL_PEER_RATIO} (indexed peers)."
            ),
            "![positional wall time](results/figures/perf_pos_wall.png)",
            (
                f"**Reading Figure {f_wall}**\n"
                "- Three panels: Genome, Transcriptome, Proteome.\n"
                "- **Bars:** zebrac mean wall time. Error bars are one standard deviation.\n"
                "- **Legend order:** z-fasta, noodles, rust-bio, samtools, fastahack, seqtk (ref).\n"
                "- **Hatched bars:** seqtk (ref) has no index and rescans the FASTA each call.\n"
                "- **On large FASTA files:** seqtk (ref) can dominate log-scale y-axis; "
                "compare indexed tools within each region, not bar height vs seqtk.\n"
                "- **Bar labels:** `1×` on z-fasta; other labels = "
                f"{LABEL_PEER_RATIO} (indexed peers only). Details in Table {t_cmp}."
            ),
        ]
    )


def md_pos_memory_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
) -> str:
    chart_tools = pos_chart_tools(df)
    if not chart_tools:
        return "_No positional memory data._"
    chart_work = df[df["tool"].isin(chart_tools)].copy()
    t_rss = nums.next_table()
    t_cmp = nums.next_table()
    f_rss = nums.next_figure()

    fig_pos_grouped_bars(
        chart_work,
        figures_dir / "perf_pos_rss.png",
        chart_tools,
        "peak_rss_mb",
        "Peak RSS (MB)",
        "Positional GET: Peak RSS vs Region Size",
        "Error bars = zebrac stddev when non-zero. Hatched bars = reference lanes. "
        f"Bar labels = RSS × ({LABEL_PEER_RATIO}).",
        ir,
        log_y=False,
        value_floor=0.1,
        annotate_comparisons=True,
    )

    rss_md = md_pos_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        value_col="peak_rss_mb",
        formatter=lambda row: f"{row['peak_rss_mb']:.2f} MB",
    )

    return "\n\n".join(
        [
            "Same zebrac samples as **Performance: Positional extraction**.",
            (
                "### Peak RSS\n\n"
                "zebrac starts a new process for each sample and records peak RSS when it "
                "exits (`ru_maxrss`). Table and figure show the mean across samples."
            ),
            f"**Table {t_rss}:** Peak RSS (MB, zebrac mean). Same tool order as Performance.",
            rss_md,
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                f"RSS × = peer peak RSS ÷ z-fasta peak RSS. Same ratios as bar labels "
                f"on Figure {f_rss}.</summary>"
            ),
            "",
            md_pos_zfasta_vs_metric_table(
                chart_work,
                chart_tools,
                "peak_rss_mb",
                ir,
                ratio_label="RSS ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.2f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.2f} MB",
            ),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_rss}:** Table {t_rss} as grouped bars. Bar labels = RSS × "
                f"({LABEL_PEER_RATIO}; see Table {t_cmp})."
            ),
            "![positional RSS](results/figures/perf_pos_rss.png)",
            (
                f"**Reading Figure {f_rss}**\n"
                "- Same three-panel layout as positional wall time; linear y-axis (MB).\n"
                "- `1×` on z-fasta; other labels = "
                f"{LABEL_PEER_RATIO}.\n"
                f"- Details in Table {t_cmp}."
            ),
        ]
    )


def md_pos_page_faults_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
) -> str:
    chart_tools = pos_chart_tools(df)
    if not chart_tools:
        return "_No positional page-fault data._"
    chart_work = df[df["tool"].isin(chart_tools)].copy()
    t_minor = nums.next_table()
    t_cmp = nums.next_table()
    t_major = nums.next_table()
    f_pf = nums.next_figure()

    fig_pos_grouped_bars(
        chart_work,
        figures_dir / "perf_pos_faults.png",
        chart_tools,
        "minor_faults",
        "Minor Page Faults",
        "Positional GET: Minor Page Faults vs Region Size",
        "Error bars = zebrac stddev. Hatched bars = reference lanes. "
        f"Bar labels = Faults × ({LABEL_PEER_RATIO}).",
        ir,
        log_y=True,
        value_floor=1.0,
        annotate_comparisons=True,
    )

    minor_md = md_pos_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        value_col="minor_faults",
        formatter=lambda row: f"{int(row['minor_faults']):,}",
    )
    major_md = md_pos_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        value_col="major_faults",
        formatter=lambda row: f"{int(row['major_faults']):,}",
    )

    return "\n\n".join(
        [
            "Same zebrac samples again.",
            (
                "### Page faults\n\n"
                "A **minor** fault maps a page without reading disk. A **major** fault reads "
                "from disk. zebrac reports both per run, like wall time."
            ),
            f"**Table {t_minor}:** Minor page faults (zebrac mean). Same tool order as Performance.",
            minor_md,
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                f"Faults × = peer minor faults / z-fasta minor faults. Same ratios as "
                f"bar labels on Figure {f_pf}.</summary>"
            ),
            "",
            md_pos_zfasta_vs_metric_table(
                chart_work,
                chart_tools,
                "minor_faults",
                ir,
                ratio_label="Faults ×",
                fmt_zf=lambda r: f"{int(r.zfasta_v):,}",
                fmt_comp=lambda r: f"{int(r.comp_v):,}",
            ),
            "</details>",
            "<details>",
            f"<summary><strong>Table {t_major}:</strong> Major page faults.</summary>",
            "",
            major_md,
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_pf}:** Table {t_minor} as grouped bars (log y). "
                f"Bar labels = Faults × ({LABEL_PEER_RATIO}; see Table {t_cmp})."
            ),
            "![positional page faults](results/figures/perf_pos_faults.png)",
            (
                f"**Reading Figure {f_pf}**\n"
                "- Same three-panel layout as positional wall time.\n"
                "- `1×` on z-fasta; other labels = "
                f"{LABEL_PEER_RATIO}.\n"
                f"- Details in Table {t_cmp}."
            ),
        ]
    )


def md_mode_comparison(df: pd.DataFrame, nums, ir, figures_dir: Path) -> str:
    sub = ir.filter_tools(df[df["region"] == "1kbp_mid"], MODE_POS_TOOLS)
    if sub.empty:
        return "_No z-fasta mode comparison data._"
    chart_ds = DATASET_ORDER
    sub = sub[sub["dataset"].isin(chart_ds)]
    t_wall = nums.next_table()
    t_cmp = nums.next_table()
    f_mode = nums.next_figure()
    fig_metric_bars(
        sub,
        figures_dir / "mode_pos.png",
        MODE_POS_TOOLS,
        "mean",
        "stddev",
        "Wall Time (s, log scale)",
        "z-fasta GET: default vs .fai lane (1 kbp_mid)",
        ir,
        fig_note="Bar labels: .fai lane wall time / default wall time.",
    )
    return "\n\n".join(
        [
            "z-fasta default (`.zfi`) vs `.fai` lane on **1kbp_mid**. "
            "The `.fai` lane stashes `.zfi` once per region outside the zebrac sample "
            "(see `run.sh`); native `.fai` tools do not pay that cost. "
            "Transcriptome **~7.6×** is mostly **index load**, not the 1 kbp fetch: the "
            "`.fai` path parses **~254k** tab-separated lines and mmaps a **~33 MiB** "
            "`.fai` on every zebrac sample, while `.zfi` mmaps a **~9.7 MiB** binary "
            "index with no line parsing. Genome (**194** seqs) and Proteome (**~21k**) pay "
            "far less parse overhead. Single-region GET uses `records_only` (no hash map); "
            "lookup is O(n) on both paths, so startup parse cost dominates on large catalogs.",
            f"**Table {t_wall}:** Wall time per dataset (Genome, Transcriptome, Proteome).",
            md_wall_table(sub, MODE_POS_TOOLS, "dataset", ir),
            "<details>",
            f"<summary><strong>Table {t_cmp}:</strong> default vs `.fai` lane.</summary>",
            "",
            ir.md_zfasta_vs_ratio_table(
                ir.build_ratio_comparisons(sub, "mean", peer_tools=MODE_POS_TOOLS[1:]),
                zf_label="z-fasta",
                comp_label="z-fasta (.fai)",
                ratio_label="Time ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                fmt_comp=lambda r: f"{r.comp_v:.4f}s",
            ),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            f"**Figure {f_mode}:** Wall time for Table {t_wall}.",
            f"![Figure {f_mode}: z-fasta mode comparison](results/figures/mode_pos.png)",
            (
                f"**Reading Figure {f_mode}**\n"
                "- Grouped bars: wall time per dataset.\n"
                f"- Bar labels match Table {t_cmp} (`.fai` lane / default)."
            ),
        ]
    )


def md_multi_detail_pivot(
    work: pd.DataFrame,
    tools: list[str],
    ir,
    *,
    value_col: str,
    formatter,
) -> str:
    rows = []
    filtered = ir.filter_tools(work, tools)
    for ds in real_datasets_in_work(filtered):
        ds_work = filtered[filtered["dataset"] == ds]
        for n in n_values_for_dataset(ds_work, ds):
            for _, row in ds_work[ds_work["n"] == n].iterrows():
                formatted = formatter(row)
                if formatted is None:
                    continue
                rows.append(
                    {
                        "dataset_n": multi_ds_n_label(ds, n),
                        "tool": ir.display_tool(row["tool"]),
                        "value": formatted,
                    }
                )
    detail = pd.DataFrame(rows)
    if detail.empty:
        return "_No data._"
    detail = detail.pivot(index="dataset_n", columns="tool", values="value")
    detail = detail.reindex(
        index=sorted(detail.index, key=multi_ds_n_sort_key),
    )
    detail.index.name = "Dataset / N"
    tool_labels = [ir.display_tool(t) for t in tools if ir.display_tool(t) in detail.columns]
    return detail.reindex(columns=tool_labels).pipe(df_to_markdown)


def md_multi_zfasta_vs_table(work: pd.DataFrame, tools: list[str], ir) -> str:
    peers = [t for t in tools if t != "z-fasta-default"]
    tagged = work.copy()
    tagged["ds_n"] = tagged.apply(
        lambda r: multi_ds_n_label(r["dataset"], int(r["n"])),
        axis=1,
    )
    comparisons = ir.build_ratio_comparisons(
        tagged,
        "mean",
        peer_tools=peers,
        group_col="ds_n",
        group_sort=multi_ds_n_sort_key,
    )
    if comparisons.empty:
        return "_No z-fasta comparison data._"
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        group_label="Dataset / N",
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label="Speedup",
        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
    )


def md_multi_zfasta_vs_metric_table(
    work: pd.DataFrame,
    tools: list[str],
    value_col: str,
    ir,
    *,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
) -> str:
    peers = [t for t in tools if t != "z-fasta-default"]
    tagged = work.copy()
    tagged["ds_n"] = tagged.apply(
        lambda r: multi_ds_n_label(r["dataset"], int(r["n"])),
        axis=1,
    )
    comparisons = ir.build_ratio_comparisons(
        tagged,
        value_col,
        peer_tools=peers,
        group_col="ds_n",
        group_sort=multi_ds_n_sort_key,
    )
    if comparisons.empty:
        return "_No z-fasta comparison data._"
    return ir.md_zfasta_vs_ratio_table(
        comparisons,
        group_label="Dataset / N",
        zf_label="z-fasta",
        comp_label="Peer",
        ratio_label=ratio_label,
        fmt_zf=fmt_zf,
        fmt_comp=fmt_comp,
    )


def md_multi_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
    manifest: dict | None = None,
) -> str:
    work = enrich_multi(df)
    chart_tools = multi_chart_tools(work)
    table_tools = multi_table_tools(work)
    if not chart_tools:
        return "_No multi-region chart data._"

    chart_work = work[work["tool"].isin(chart_tools)].copy()
    table_work = work[work["tool"].isin(table_tools)].copy()
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_cmp = nums.next_table()
    f_wall = nums.next_figure()
    sample_n = manifest_sample_count(manifest)
    ref_fig_note = (
        " Hatched bars = reference loops (seqtk, fastahack)."
        if any(t in REFERENCE_TOOLS for t in chart_tools)
        else ""
    )
    cmp_ref_blurb = (
        " Includes reference loops when run with `GET_MULTI_REFERENCE=1`."
        if any(t in REFERENCE_TOOLS for t in table_tools)
        else ""
    )

    fig_multi_grouped_bars(
        chart_work,
        figures_dir / "perf_multi_wall.png",
        chart_tools,
        "mean",
        "Wall Time (s)",
        "Multi-region GET: Wall Time vs N",
        (
            f"1 kbp per region. Error bars = zebrac stddev (n={sample_n}).{ref_fig_note} "
            f"Bar labels = {LABEL_PEER_RATIO} (indexed peers)."
        ),
        ir,
        log_y=True,
        annotate_comparisons=True,
    )

    wall_md = md_multi_detail_pivot(
        table_work,
        table_tools,
        ir,
        value_col="mean",
        formatter=lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
    )
    throughput_md = md_multi_detail_pivot(
        table_work,
        table_tools,
        ir,
        value_col="throughput_mibs",
        formatter=lambda row: (
            fmt_mbp_s(float(row["throughput_mibs"]))
            if pd.notna(row.get("throughput_mibs"))
            else None
        ),
    )

    return "\n\n".join(
        [
            md_multi_intro_blurb(work, table_tools),
            (
                f"**Table {t_wall}:** Wall time (mean ± stddev, seconds). Rows sorted by "
                "dataset then N. Lower is better. Same tool order as the figure."
            ),
            wall_md,
            "<details>",
            (
                f"<summary><strong>Table {t_tp}:</strong> Output throughput (Mbp/s) from "
                f"metadata <code>output_bases</code> ÷ wall time. Higher is better. Same "
                f"rows and tool order as Table {t_wall}.</summary>"
            ),
            "",
            throughput_md,
            "</details>",
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                "Speedup = peer wall time ÷ z-fasta wall time (higher = z-fasta faster)."
                f"{cmp_ref_blurb} Same ratios as bar labels on "
                f"Figure {f_wall}.</summary>"
            ),
            "",
            md_multi_zfasta_vs_table(table_work, table_tools, ir),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_wall}:** Table {t_wall} as grouped bars (log y).{ref_fig_note} "
                f"Bar labels = {LABEL_PEER_RATIO} (indexed peers)."
            ),
            "![multi-region wall time](results/figures/perf_multi_wall.png)",
            md_multi_figure_reading_base(
                chart_tools,
                y_note=(
                    "**Bars:** zebrac mean wall time. Error bars = one standard deviation."
                ),
                ratio_label=f"{LABEL_PEER_RATIO} (indexed peers only)",
                t_cmp=t_cmp,
                f_num=f_wall,
            ),
        ]
    )


def md_multi_memory_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
) -> str:
    work = enrich_multi(df)
    chart_tools = multi_chart_tools(work)
    if not chart_tools:
        return "_No multi-region memory data._"
    chart_work = work[work["tool"].isin(chart_tools)].copy()
    t_rss = nums.next_table()
    t_cmp = nums.next_table()
    f_rss = nums.next_figure()
    ref_fig_note = (
        " Hatched bars = reference loops."
        if any(t in REFERENCE_TOOLS for t in chart_tools)
        else ""
    )

    fig_multi_grouped_bars(
        chart_work,
        figures_dir / "perf_multi_rss.png",
        chart_tools,
        "peak_rss_mb",
        "Peak RSS (MB)",
        "Multi-region GET: Peak RSS vs N",
        (
            "Error bars = zebrac stddev when non-zero."
            f"{ref_fig_note} Bar labels = RSS × (peer peak RSS ÷ z-fasta peak RSS)."
        ),
        ir,
        log_y=False,
        value_floor=0.1,
        annotate_comparisons=True,
    )

    rss_md = md_multi_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        value_col="peak_rss_mb",
        formatter=lambda row: f"{row['peak_rss_mb']:.2f} MB",
    )

    return "\n\n".join(
        [
            "Same zebrac samples as **Performance: Multi-region scaling**.",
            (
                "### Peak RSS\n\n"
                "zebrac starts a new process for each sample and records peak RSS when it "
                "exits (`ru_maxrss`). Table and figure show the mean across samples."
            ),
            f"**Table {t_rss}:** Peak RSS (MB, zebrac mean). Same tool order as Performance.",
            rss_md,
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                f"RSS × = peer peak RSS ÷ z-fasta peak RSS. Same ratios as bar labels "
                f"on Figure {f_rss}.</summary>"
            ),
            "",
            md_multi_zfasta_vs_metric_table(
                chart_work,
                chart_tools,
                "peak_rss_mb",
                ir,
                ratio_label="RSS ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.2f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.2f} MB",
            ),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_rss}:** Table {t_rss} as grouped bars (linear y, MB). "
                f"Bar labels = RSS × (see Table {t_cmp})."
            ),
            "![multi-region RSS](results/figures/perf_multi_rss.png)",
            md_multi_figure_reading_base(
                chart_tools,
                y_note=(
                    "**Bars:** zebrac mean peak RSS. Error bars when stddev is non-zero. "
                    "Linear y-axis (MB)."
                ),
                ratio_label="RSS × (peer peak RSS ÷ z-fasta peak RSS)",
                t_cmp=t_cmp,
                f_num=f_rss,
            ),
        ]
    )


def md_multi_page_faults_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
) -> str:
    work = enrich_multi(df)
    chart_tools = multi_chart_tools(work)
    if not chart_tools:
        return "_No multi-region page-fault data._"
    chart_work = work[work["tool"].isin(chart_tools)].copy()
    t_minor = nums.next_table()
    t_cmp = nums.next_table()
    t_major = nums.next_table()
    f_pf = nums.next_figure()
    ref_fig_note = (
        " Hatched bars = reference loops."
        if any(t in REFERENCE_TOOLS for t in chart_tools)
        else ""
    )

    fig_multi_grouped_bars(
        chart_work,
        figures_dir / "perf_multi_faults.png",
        chart_tools,
        "minor_faults",
        "Minor Page Faults",
        "Multi-region GET: Minor Page Faults vs N",
        (
            "Error bars = zebrac stddev."
            f"{ref_fig_note} Bar labels = Faults × (peer minor faults ÷ z-fasta minor faults)."
        ),
        ir,
        log_y=True,
        value_floor=1.0,
        annotate_comparisons=True,
    )

    minor_md = md_multi_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        value_col="minor_faults",
        formatter=lambda row: f"{int(row['minor_faults']):,}",
    )
    major_md = md_multi_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        value_col="major_faults",
        formatter=lambda row: f"{int(row['major_faults']):,}",
    )

    major_all_zero = (
        chart_work["major_faults"].fillna(0).eq(0).all()
        if "major_faults" in chart_work.columns and not chart_work.empty
        else True
    )
    major_summary = (
        f"<summary><strong>Table {t_major}:</strong> Major page faults (all zero for "
        "this run).</summary>"
        if major_all_zero
        else f"<summary><strong>Table {t_major}:</strong> Major page faults.</summary>"
    )

    return "\n\n".join(
        [
            "Same zebrac samples as **Performance: Multi-region scaling**.",
            (
                "### Page faults\n\n"
                "A **minor** fault maps a page without reading disk. A **major** fault reads "
                "from disk. zebrac reports both per run, like wall time."
            ),
            f"**Table {t_minor}:** Minor page faults (zebrac mean). Same tool order as Performance.",
            minor_md,
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                f"Faults × = peer minor faults ÷ z-fasta minor faults. Same ratios as "
                f"bar labels on Figure {f_pf}.</summary>"
            ),
            "",
            md_multi_zfasta_vs_metric_table(
                chart_work,
                chart_tools,
                "minor_faults",
                ir,
                ratio_label="Faults ×",
                fmt_zf=lambda r: f"{int(r.zfasta_v):,}",
                fmt_comp=lambda r: f"{int(r.comp_v):,}",
            ),
            "</details>",
            "<details>",
            major_summary,
            "",
            major_md,
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_pf}:** Table {t_minor} as grouped bars (log y). "
                f"Bar labels = Faults × (see Table {t_cmp})."
            ),
            "![multi-region page faults](results/figures/perf_multi_faults.png)",
            md_multi_figure_reading_base(
                chart_tools,
                y_note="**Bars:** zebrac mean minor page faults. Log y-axis.",
                ratio_label="Faults × (peer minor faults ÷ z-fasta minor faults)",
                t_cmp=t_cmp,
                f_num=f_pf,
            ),
        ]
    )


def md_bed_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
    manifest: dict | None = None,
) -> str:
    work = bed_filter_work(df)
    chart_tools = bed_chart_tools(work)
    if not chart_tools:
        return "_No BED batch chart data._"

    chart_work = work[work["tool"].isin(chart_tools)].copy()
    mode_tools = bed_tools_present(work, BED_MODE_TOOLS)
    mode_work = ir.filter_tools(work, mode_tools) if mode_tools else pd.DataFrame()
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_cmp = nums.next_table()
    t_mode = nums.next_table() if not mode_work.empty else None
    f_wall = nums.next_figure()
    sample_n = manifest_sample_count(manifest)

    fig_bed_grouped_bars(
        chart_work,
        figures_dir / "perf_bed_wall.png",
        chart_tools,
        "mean",
        "Wall Time (s)",
        "BED batch GET: Wall Time vs Row Count",
        (
            f"1 kbp per row. Error bars = zebrac stddev (n={sample_n}). "
            f"Bar labels = {LABEL_PEER_RATIO} (indexed peers + bedtools)."
        ),
        ir,
        log_y=True,
        annotate_comparisons=True,
    )

    wall_md = md_bed_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        formatter=lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
    )
    throughput_md = md_bed_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        formatter=lambda row: (
            fmt_mbp_s(float(row["throughput_mibs"]))
            if pd.notna(row.get("throughput_mibs"))
            else None
        ),
    )

    blocks = [
        md_bed_intro_blurb(work, t_mode=t_mode),
        (
            f"**Table {t_wall}:** Wall time (mean ± stddev, seconds). Rows sorted by "
            "dataset then BED row count. Lower is better. Same tool order as the figure."
        ),
        wall_md,
        "<details>",
        (
            f"<summary><strong>Table {t_tp}:</strong> Output throughput (Mbp/s) from "
            f"metadata <code>output_bases</code> ÷ wall time. Higher is better. Same "
            f"rows and tool order as Table {t_wall}.</summary>"
        ),
        "",
        throughput_md,
        "</details>",
        "<details>",
        (
            f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
            "Speedup = peer wall time ÷ z-fasta wall time (higher = z-fasta faster). "
            f"Same ratios as bar labels on Figure {f_wall}.</summary>"
        ),
        "",
        md_bed_zfasta_vs_table(chart_work, chart_tools, ir),
        "</details>",
    ]
    if t_mode is not None:
        blocks.extend(
            [
                "<details>",
                (
                    f"<summary><strong>Table {t_mode}:</strong> z-fasta chunk-size comparison "
                    f"(default {BED_DEFAULT_CHUNK} vs -1 vs 1; not in Figure {f_wall}).</summary>"
                ),
                "",
                md_bed_detail_pivot(
                    mode_work,
                    mode_tools,
                    ir,
                    mode_table=True,
                    formatter=lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
                ),
                "</details>",
            ]
        )
    blocks.extend(
        [
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_wall}:** Table {t_wall} as grouped bars (log y). "
                f"Bar labels = {LABEL_PEER_RATIO} (indexed peers + bedtools)."
            ),
            "![bed batch wall time](results/figures/perf_bed_wall.png)",
            md_bed_figure_reading_base(
                y_note=(
                    "**Bars:** zebrac mean wall time. Error bars = one standard deviation."
                ),
                ratio_label=f"{LABEL_PEER_RATIO} (indexed peers + bedtools)",
                t_cmp=t_cmp,
                f_num=f_wall,
                include_chunk_note=True,
            ),
        ]
    )
    return "\n\n".join(blocks)


def md_bed_memory_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
) -> str:
    work = bed_filter_work(df)
    chart_tools = bed_chart_tools(work)
    if not chart_tools:
        return "_No BED batch memory data._"
    chart_work = work[work["tool"].isin(chart_tools)].copy()
    t_rss = nums.next_table()
    t_cmp = nums.next_table()
    f_rss = nums.next_figure()

    fig_bed_grouped_bars(
        chart_work,
        figures_dir / "perf_bed_rss.png",
        chart_tools,
        "peak_rss_mb",
        "Peak RSS (MB)",
        "BED batch GET: Peak RSS vs Row Count",
        (
            "Error bars = zebrac stddev when non-zero. "
            "Bar labels = RSS × (peer peak RSS ÷ z-fasta peak RSS)."
        ),
        ir,
        log_y=False,
        value_floor=0.1,
        annotate_comparisons=True,
    )

    rss_md = md_bed_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        formatter=lambda row: f"{row['peak_rss_mb']:.2f} MB",
    )

    return "\n\n".join(
        [
            "Same zebrac samples as **Performance: BED batch**.",
            (
                "### Peak RSS\n\n"
                "zebrac starts a new process for each sample and records peak RSS when it "
                "exits (`ru_maxrss`). Table and figure show the mean across samples."
            ),
            f"**Table {t_rss}:** Peak RSS (MB, zebrac mean). Same tool order as Performance.",
            rss_md,
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                f"RSS × = peer peak RSS ÷ z-fasta peak RSS. Same ratios as bar labels "
                f"on Figure {f_rss}.</summary>"
            ),
            "",
            md_bed_zfasta_vs_metric_table(
                chart_work,
                chart_tools,
                ir,
                value_col="peak_rss_mb",
                ratio_label="RSS ×",
                fmt_zf=lambda r: f"{r.zfasta_v:.2f} MB",
                fmt_comp=lambda r: f"{r.comp_v:.2f} MB",
            ),
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_rss}:** Table {t_rss} as grouped bars (linear y, MB). "
                f"Bar labels = RSS × (see Table {t_cmp})."
            ),
            "![bed batch RSS](results/figures/perf_bed_rss.png)",
            md_bed_figure_reading_base(
                y_note=(
                    "**Bars:** zebrac mean peak RSS. Error bars when stddev is non-zero. "
                    "Linear y-axis (MB)."
                ),
                ratio_label="RSS × (peer peak RSS ÷ z-fasta peak RSS)",
                t_cmp=t_cmp,
                f_num=f_rss,
            ),
        ]
    )


def md_bed_page_faults_section(df: pd.DataFrame, nums, ir, figures_dir: Path) -> str:
    work = bed_filter_work(df)
    chart_tools = bed_chart_tools(work)
    if not chart_tools:
        return "_No BED batch page-fault data._"
    chart_work = work[work["tool"].isin(chart_tools)].copy()
    t_minor = nums.next_table()
    t_cmp = nums.next_table()
    t_major = nums.next_table()
    f_pf = nums.next_figure()

    fig_bed_grouped_bars(
        chart_work,
        figures_dir / "perf_bed_faults.png",
        chart_tools,
        "minor_faults",
        "Minor Page Faults",
        "BED batch GET: Minor Page Faults vs Row Count",
        (
            "Error bars = zebrac stddev. "
            "Bar labels = Faults × (peer minor faults ÷ z-fasta minor faults)."
        ),
        ir,
        log_y=True,
        value_floor=1.0,
        annotate_comparisons=True,
    )

    minor_md = md_bed_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        formatter=lambda row: f"{int(row['minor_faults']):,}",
    )
    major_md = md_bed_detail_pivot(
        chart_work,
        chart_tools,
        ir,
        formatter=lambda row: f"{int(row['major_faults']):,}",
    )
    major_all_zero = (
        chart_work["major_faults"].fillna(0).eq(0).all()
        if "major_faults" in chart_work.columns and not chart_work.empty
        else True
    )
    major_summary = (
        f"<summary><strong>Table {t_major}:</strong> Major page faults (all zero for "
        "this run).</summary>"
        if major_all_zero
        else f"<summary><strong>Table {t_major}:</strong> Major page faults.</summary>"
    )

    return "\n\n".join(
        [
            "Same zebrac samples as **Performance: BED batch**.",
            (
                "### Page faults\n\n"
                "A **minor** fault maps a page without reading disk. A **major** fault reads "
                "from disk. zebrac reports both per run, like wall time."
            ),
            f"**Table {t_minor}:** Minor page faults (zebrac mean). Same tool order as Performance.",
            minor_md,
            "<details>",
            (
                f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each peer. "
                f"Faults × = peer minor faults ÷ z-fasta minor faults. Same ratios as "
                f"bar labels on Figure {f_pf}.</summary>"
            ),
            "",
            md_bed_zfasta_vs_metric_table(
                chart_work,
                chart_tools,
                ir,
                value_col="minor_faults",
                ratio_label="Faults ×",
                fmt_zf=lambda r: f"{int(r.zfasta_v):,}",
                fmt_comp=lambda r: f"{int(r.comp_v):,}",
            ),
            "</details>",
            "<details>",
            major_summary,
            "",
            major_md,
            "</details>",
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_pf}:** Table {t_minor} as grouped bars (log y). "
                f"Bar labels = Faults × (see Table {t_cmp})."
            ),
            "![bed batch page faults](results/figures/perf_bed_faults.png)",
            md_bed_figure_reading_base(
                y_note="**Bars:** zebrac mean minor page faults. Log y-axis.",
                ratio_label="Faults × (peer minor faults ÷ z-fasta minor faults)",
                t_cmp=t_cmp,
                f_num=f_pf,
            ),
        ]
    )


def md_rc_detail_pivot(
    work: pd.DataFrame,
    tools: list[str],
    ir,
    *,
    value_col: str,
    formatter,
) -> str:
    rows = []
    for _, row in ir.filter_tools(work, tools).iterrows():
        formatted = formatter(row)
        if formatted is None:
            continue
        rows.append(
            {
                "dataset": row["dataset"],
                "tool": ir.display_tool(row["tool"]),
                "value": formatted,
            }
        )
    detail = pd.DataFrame(rows)
    if detail.empty:
        return "_No data._"
    detail = detail.pivot(index="dataset", columns="tool", values="value")
    tool_labels = [ir.display_tool(t) for t in tools if ir.display_tool(t) in detail.columns]
    ordered_ds = [d for d in RC_DATASET_ORDER if d in detail.index]
    return df_to_markdown(detail.reindex(index=ordered_ds, columns=tool_labels))


def md_rc_self_overhead_table(work: pd.DataFrame) -> str:
    rows = []
    for ds in RC_DATASET_ORDER:
        sub = work[work["dataset"] == ds]
        if sub.empty:
            continue
        overhead = rc_zfasta_overhead_ratio(sub, ds)
        if overhead is None:
            continue
        rows.append({"dataset": ds, "z-fasta --rc / plain": f"{overhead:.3f}×"})
    if not rows:
        return "_No z-fasta RC self-overhead data._"
    return df_to_markdown(pd.DataFrame(rows), index=False)


def md_rc_section(
    df: pd.DataFrame,
    nums,
    ir,
    figures_dir: Path,
    manifest: dict | None = None,
) -> str:
    work = rc_chart_work(df)
    if work.empty:
        return "_No RC data._"
    tools = rc_tools_for_dataset(work, RC_HEADLINE_TOOLS)
    if not tools:
        return "_No RC data._"
    datasets = [d for d in RC_DATASET_ORDER if d in set(work["dataset"].astype(str))]
    t_wall = nums.next_table()
    t_rss = nums.next_table()
    t_pf = nums.next_table()
    t_tp = nums.next_table()
    t_self = nums.next_table()
    t_cmp = nums.next_table()
    f_rc = nums.next_figure()
    sample_n = manifest_sample_count(manifest)
    rc_peers = rc_peer_tools(tools)

    fig_rc_overhead(work, figures_dir / "perf_rc.png", tools, ir, sample_n=sample_n)

    wall_md = md_rc_detail_pivot(
        work,
        tools,
        ir,
        value_col="mean",
        formatter=lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
    )
    rss_md = md_rc_detail_pivot(
        work,
        tools,
        ir,
        value_col="peak_rss_mb",
        formatter=lambda row: f"{row['peak_rss_mb']:.2f} MB",
    )
    pf_md = md_rc_detail_pivot(
        work,
        tools,
        ir,
        value_col="minor_faults",
        formatter=lambda row: f"{int(row['minor_faults']):,}",
    )
    tp_md = md_rc_detail_pivot(
        work,
        tools,
        ir,
        value_col="throughput_mibs",
        formatter=lambda row: (
            fmt_mbp_s(float(row["throughput_mibs"]))
            if pd.notna(row.get("throughput_mibs"))
            else None
        ),
    )

    return "\n\n".join(
        [
            "Orientation transforms on a fixed **1 kbp** mid slice (`1kbp_mid`) on "
            f"**{', '.join(datasets)}** only. **Proteome is excluded:** protein records have "
            "no reverse-complement semantics, so `--rc` and `--complement-only` do not apply; "
            "`--reverse-only` only permutes amino-acid characters and adds negligible overhead "
            "vs plain GET, so it is not a meaningful benchmark axis. **X-axis = dataset** "
            "(Genome, Transcriptome); **bar color = tool/lane** (same layout as z-fasta mode "
            "comparison). Figure facets are wall time, peak RSS, and minor page faults. Primary "
            "comparison: z-fasta plain vs `--rc`. RC-capable peers use native `--rc`; seqtk (ref) "
            "uses `subseq` + `seq -r` (no index).",
            f"**Table {t_wall}:** Wall time (mean ± stddev, seconds).",
            wall_md,
            f"**Table {t_rss}:** Peak RSS (MB).",
            rss_md,
            f"**Table {t_pf}:** Minor page faults.",
            pf_md,
            "<details>",
            (
                f"<summary><strong>Table {t_tp}:</strong> Output throughput (Mbp/s).</summary>"
            ),
            "",
            tp_md,
            "</details>",
            f"**Table {t_self}:** z-fasta RC self-overhead (`--rc` / plain).",
            md_rc_self_overhead_table(work),
            *(
                [
                    "<details>",
                    (
                        f"<summary><strong>Table {t_cmp}:</strong> RC peers vs z-fasta --rc. "
                        f"Ratio = {LABEL_PEER_WALL_RATIO} (baseline = z-fasta --rc).</summary>"
                    ),
                    "",
                    ir.md_zfasta_vs_ratio_table(
                        ir.build_ratio_comparisons(
                            rc_chart_work(work),
                            "mean",
                            baseline="z-fasta-rc",
                            peer_tools=rc_peers,
                        ),
                        zf_label="z-fasta --rc",
                        comp_label="Peer",
                        ratio_label="Time ×",
                        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
                        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
                    ),
                    "</details>",
                ]
                if rc_peers
                else []
            ),
            '<div style="margin: 1.5em 0"></div>',
            (
                f"**Figure {f_rc}:** Tables {t_wall}, {t_rss}, and {t_pf} as grouped bars "
                "(datasets on x-axis; tool color from legend)."
            ),
            "![RC overhead](results/figures/perf_rc.png)",
            (
                f"**Reading Figure {f_rc}**\n"
                "- **Facets:** wall time (log y) | peak RSS (linear) | minor page faults (log y).\n"
                "- **X-axis:** Genome, Transcriptome (Proteome omitted; see section intro).\n"
                "- **Bar colors:** one lane per tool (legend below figure).\n"
                "- **Bar labels (rotated):** plain = `1×` baseline; z-fasta `--rc`, "
                "`--complement-only`, `--reverse-only` vs plain; noodles/rust-bio `--rc` "
                "and seqtk (ref) vs z-fasta `--rc`. Same ratio rules on wall time, peak RSS, "
                "and page faults.\n"
                "- **Hatched bars:** plain z-fasta (no transform) and seqtk (ref)."
            ),
        ]
    )


def md_messy_section(df: pd.DataFrame, nums, ir, figures_dir: Path, manifest: dict | None) -> str:
    work = filter_messy_perf_work(enrich_messy(df))
    if work.empty:
        return "_No messy FASTA data._"
    cfg = load_messy_perf_config()
    variants = messy_variants_in_work(work)
    headline = cfg["headline_span"]
    t_catalog = nums.next_table()
    t_headline = nums.next_table()
    t_detail = nums.next_table()
    f_headline = nums.next_figure()
    f_summary = nums.next_figure()

    fig_messy_ratio_headline(work, figures_dir / "messy_ratio_headline.png", variants, ir)
    fig_messy_ratio_summary(work, figures_dir / "messy_ratio_summary.png", variants, ir)

    catalog = messy_region_catalog_table()
    headline_rows = []
    for variant in variants:
        for metric, label in (
            ("mean", "wall time ×"),
            ("peak_rss_mb", "peak RSS ×"),
            ("minor_faults", "minor faults ×"),
        ):
            ratio = messy_overhead_ratio(work, variant, headline, metric)
            u_val = messy_layout_value(work, variant, headline, "uniform", metric)
            m_val = messy_layout_value(work, variant, headline, "messy", metric)
            if ratio is None or u_val is None or m_val is None:
                continue
            if metric == "mean":
                uniform_fmt = fmt_us(u_val)
                messy_fmt = fmt_us(m_val)
            elif metric == "peak_rss_mb":
                uniform_fmt = f"{u_val:.2f} MB"
                messy_fmt = f"{m_val:.2f} MB"
            else:
                uniform_fmt = f"{int(u_val):,}"
                messy_fmt = f"{int(m_val):,}"
            headline_rows.append(
                {
                    "variant": variant,
                    "span": headline,
                    "metric": label,
                    "uniform": uniform_fmt,
                    "messy": messy_fmt,
                    "messy / uniform": f"{ratio:.3f}×",
                }
            )
    headline_df = pd.DataFrame(headline_rows)
    detail_df = messy_detail_table(work)

    zebrac_note = ""
    if manifest and manifest.get("messy_runs"):
        zebrac_note = (
            f" Zebrac used **{manifest['messy_runs']}** samples, "
            f"**{manifest.get('messy_warmup', '?')}** warmup, "
            f"**{manifest.get('messy_duration_ms', '?')} ms** minimum per lane."
        )

    n_spans = len(cfg["span_order"])
    n_variants = len(cfg["variants"])
    return "\n\n".join(
        [
            "Paired `get` on proteome-derived messy fixtures in `bench/shared/messy_perf/` "
            f"({n_variants} layout variants, ~20k sequences each). Each workload extracts the same "
            "TITIN sub-regions from a **uniform** (`validate --fix`) and **messy** (side-table "
            f"`.zfi`) copy of each file. Spans are defined in `bench/get/messy_perf.json`. Index "
            f"build stays outside zebrac.{zebrac_note} Peers are omitted.",
            (
                "**Why ratio charts:** Figures plot **messy / uniform**; dashed line = parity (1.0×). "
                f"Figure {f_headline} quotes the headline span `{headline}`. "
                f"Figure {f_summary} summarizes all spans (median bar, min-max whiskers, "
                f"`{headline}` diamond) to show overhead tracks layout type, not extraction length. "
                "Collapsible Table below has per-span microseconds and ratios."
            ),
            f"**Table {t_catalog}:** Perf spans × variants ({n_variants * n_spans} timed pairs).",
            df_to_markdown(catalog, index=False),
            (
                f"**Reading Table {t_catalog}**\n"
                "- **span** = short id from `messy_perf.json`.\n"
                "- **region** = full `z-fasta get` region argument (`&#124;` = `|` in UniProt ids).\n"
                f"- **headline** = `{headline}` (same slice on every variant file)."
            ),
            f"**Table {t_headline}:** Headline span `{headline}` messy / uniform.",
            df_to_markdown(headline_df, index=False) if not headline_df.empty else "_No headline data._",
            '<div style="margin: 1.5em 0"></div>',
            f"**Figure {f_headline}:** Headline span `{headline}` ratios (Table {t_headline}).",
            "![messy ratio headline](results/figures/messy_ratio_headline.png)",
            (
                f"**Reading Figure {f_headline}**\n"
                "- Three panels: wall time, peak RSS, minor page faults (messy / uniform).\n"
                f"- One bar per layout variant at the headline slice `{headline}`."
            ),
            f"**Figure {f_summary}:** Layout overhead stable across spans (see Table {t_detail}).",
            "![messy ratio summary](results/figures/messy_ratio_summary.png)",
            (
                f"**Reading Figure {f_summary}**\n"
                "- Bar = median messy / uniform ratio over all perf spans.\n"
                "- Whiskers = min-max across spans (tight whiskers = span-insensitive overhead).\n"
                f"- Diamond = headline span `{headline}` (links back to Figure {f_headline})."
            ),
            "<details>",
            (
                f"<summary><strong>Table {t_detail}:</strong> All spans: wall (µs), RSS ×, "
                f"faults × (per-span drill-down).</summary>"
            ),
            "",
            df_to_markdown(detail_df, index=False) if not detail_df.empty else "_No detail data._",
            "</details>",
        ]
    )


def expected_report_figures(
    multi_df: pd.DataFrame | None,
    bed_df: pd.DataFrame | None,
) -> frozenset[str]:
    """PNG filenames referenced by the current REPORT for the loaded run."""
    names = set(CORE_REPORT_FIGURES)
    if multi_df is not None and not multi_df.empty:
        names |= MULTI_REPORT_FIGURES
    if bed_df is not None and not bed_df.empty:
        names |= BED_REPORT_FIGURES
    return frozenset(names)


def prune_stale_figures(figures_dir: Path, allowed: frozenset[str]) -> list[str]:
    """Remove PNGs not referenced by the current report layout."""
    removed: list[str] = []
    if not figures_dir.is_dir():
        return removed
    for path in figures_dir.glob("*.png"):
        if path.name in allowed:
            continue
        path.unlink()
        removed.append(path.name)
    return removed


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate GET benchmark report.")
    parser.add_argument("results_dir", nargs="?", type=Path, default=RESULTS_DIR)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()

    ir = load_index_report()
    results_dir = args.results_dir
    manifest = load_manifest(results_dir)
    if is_incomplete(manifest) and not args.allow_incomplete:
        ts = (manifest or {}).get("timestamp", "unknown")
        raise SystemExit(
            f"Refusing to overwrite REPORT.md from incomplete run {ts}. "
            "Pass --allow-incomplete for drafts."
        )

    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    nums = ir.ReportCounters()
    report_lines: list[str] = [
        "<!-- markdownlint-disable MD024 MD032 MD033 MD036 MD041 MD049 -->",
        "# z-fasta GET Benchmark Report",
        "_Auto-generated by `bench/get/generate_report.py` from zebrac results._",
    ]
    if is_incomplete(manifest):
        report_lines.append(
            "> **DRAFT: incomplete run.** Missing headline sections or below minimum "
            "zebrac sample settings. Do not publish until a full run completes without "
            "missing sections or `skip_*` flags on headline workloads."
        )
    report_lines.extend(
        [
        "## Overview",
        md_overview(manifest),
        "## Run Provenance",
        md_run_provenance(manifest, PROJECT_ROOT, ir),
        ]
    )

    pos_df = load_section(results_dir, manifest, "perf_pos", ir)
    multi_df = load_section(results_dir, manifest, "perf_multi", ir)
    bed_df = load_section(results_dir, manifest, "perf_bed", ir)
    rc_df = load_section(results_dir, manifest, "perf_rc", ir)
    messy_df = load_section(results_dir, manifest, "messy", ir)

    if pos_df is not None and not pos_df.empty:
        report_lines.append("## Performance: Positional extraction")
        report_lines.append(md_pos_section(pos_df, nums, ir, figures_dir, manifest))
        report_lines.append("## Memory Usage: Positional extraction")
        report_lines.append(md_pos_memory_section(pos_df, nums, ir, figures_dir))
        report_lines.append("## Page Faults: Positional extraction")
        report_lines.append(md_pos_page_faults_section(pos_df, nums, ir, figures_dir))
    else:
        report_lines.append("## Performance: Positional extraction")
        report_lines.append("_No positional zebrac data found._")

    if multi_df is not None and not multi_df.empty:
        report_lines.append("## Performance: Multi-region scaling")
        report_lines.append(md_multi_section(multi_df, nums, ir, figures_dir, manifest))
        report_lines.append("## Memory Usage: Multi-region scaling")
        report_lines.append(md_multi_memory_section(multi_df, nums, ir, figures_dir))
        report_lines.append("## Page Faults: Multi-region scaling")
        report_lines.append(md_multi_page_faults_section(multi_df, nums, ir, figures_dir))
    elif manifest and manifest.get("skip_multi"):
        report_lines.append("## Performance: Multi-region scaling")
        report_lines.append("_This run skipped multi-region scaling (`skip_multi=true`)._")
    elif manifest and not manifest.get("skip_multi"):
        report_lines.append("## Performance: Multi-region scaling")
        report_lines.append("_No multi-region zebrac data found._")

    if bed_df is not None and not bed_df.empty:
        report_lines.append("## Performance: BED batch")
        report_lines.append(md_bed_section(bed_df, nums, ir, figures_dir, manifest))
        report_lines.append("## Memory Usage: BED batch")
        report_lines.append(md_bed_memory_section(bed_df, nums, ir, figures_dir))
        report_lines.append("## Page Faults: BED batch")
        report_lines.append(md_bed_page_faults_section(bed_df, nums, ir, figures_dir))
    elif manifest and manifest.get("skip_bed"):
        report_lines.append("## Performance: BED batch")
        report_lines.append("_This run skipped BED batch (`skip_bed=true`)._")
    elif manifest and not manifest.get("skip_bed"):
        report_lines.append("## Performance: BED batch")
        report_lines.append("_No BED batch zebrac data found._")

    if rc_df is not None and not rc_df.empty:
        report_lines.append("## RC overhead")
        report_lines.append(md_rc_section(rc_df, nums, ir, figures_dir, manifest))
    elif manifest and manifest.get("skip_rc"):
        report_lines.append("## RC overhead")
        report_lines.append("_This run skipped RC overhead (`skip_rc=true`)._")
    elif manifest and not manifest.get("skip_rc"):
        report_lines.append("## RC overhead")
        report_lines.append("_No RC zebrac data found._")

    if pos_df is not None and not pos_df.empty:
        report_lines.append("## z-fasta Mode Comparison")
        report_lines.append(md_mode_comparison(pos_df, nums, ir, figures_dir))

    if messy_df is not None and not messy_df.empty:
        report_lines.append("## Messy FASTA GET (z-fasta only)")
        report_lines.append(md_messy_section(messy_df, nums, ir, figures_dir, manifest))
    elif manifest and manifest.get("skip_messy"):
        report_lines.append("## Messy FASTA GET (z-fasta only)")
        report_lines.append("_This run skipped messy FASTA (`skip_messy=true`)._")

    report_path = SCRIPT_DIR / "REPORT.md"
    allowed_figures = expected_report_figures(multi_df, bed_df)
    prune_stale_figures(figures_dir, allowed_figures)
    report_path.write_text(ir.normalize_markdown(report_lines))
    print(report_path)


if __name__ == "__main__":
    main()
