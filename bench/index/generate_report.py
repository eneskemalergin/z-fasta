#!/usr/bin/env python3
"""
z-fasta Index Benchmark Report Generator

Reads zebrac JSON joined with metadata JSONL from bench/index/results/,
produces Markdown report + PNG figures using pandas + matplotlib.

Usage:
    .venv/bin/python bench/index/generate_report.py [results_dir]

Defaults to bench/index/results/ (latest run via results/LATEST).
"""

import json
import argparse
import shutil
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import tabulate  # noqa: F401  -- needed by pd.to_markdown()

MARKDOWNLINT_DISABLE = "<!-- markdownlint-disable MD024 MD032 MD033 MD036 MD041 MD049 -->"


class ReportCounters:
    """Global table and figure numbers across the report."""

    def __init__(self) -> None:
        self.table = 1
        self.figure = 1

    def next_table(self) -> int:
        n = self.table
        self.table += 1
        return n

    def next_figure(self) -> int:
        n = self.figure
        self.figure += 1
        return n

# ── Styling ────────────────────────────────────────────────────────
COLORS = {
    "z-fasta-fai": "#F7A41D",  # Zig gold (peer-comparable)
    "z-fasta-fai-nodedup": "#FBC02D",
    "z-fasta-zfi": "#E65100",  # deep orange (first-class)
    "z-fasta-zfi-nodedup": "#EF6C00",
    "noodles": "#C45C26",  # rust bronze
    "rustbio-custom-index": "#8B3A2A",  # rust red-brown
    "samtools": "#555555",  # C grey (GitHub linguist)
    "seqkit": "#009485",  # Go teal
    "fastahack": "#F34B7D",  # C++ pink (GitHub linguist)
    "pyfaidx": "#3572A5",  # Python blue (GitHub linguist)
}
DISPLAY_NAMES = {
    "z-fasta-fai": "z-fasta (.fai)",
    "z-fasta-fai-nodedup": "z-fasta (.fai) --no-dedup",
    "z-fasta-zfi": "z-fasta (.zfi)",
    "z-fasta-zfi-nodedup": "z-fasta (.zfi) --no-dedup",
    "samtools": "samtools",
    "seqkit": "seqkit",
    "fastahack": "fastahack",
    "pyfaidx": "pyfaidx",
    "noodles": "noodles",
    "rustbio-custom-index": "rust-bio",
}
EDGE_EXIT_COLUMNS = (
    "zfasta_exit",
    "zfasta_zfi_exit",
    "samtools_exit",
    "seqkit_exit",
    "fastahack_exit",
    "noodles_exit",
    "rustbio_exit",
)
EDGE_EXIT_LABELS = {
    "zfasta_exit": "z-fasta (.fai)",
    "zfasta_zfi_exit": "z-fasta (.zfi)",
    "samtools_exit": "samtools",
    "seqkit_exit": "seqkit",
    "fastahack_exit": "fastahack",
    "noodles_exit": "noodles",
    "rustbio_exit": "rust-bio",
}
EDGE_CONTRACT_CLASSES = {
    "fai_parity": "FAI parity",
    "zfi_messy": "ZFI messy support",
    "invalid_input": "Invalid input review",
}
EDGE_MESSY_CASES = {"uniform", "mixed_widths", "trailing_whitespace", "blank_lines", "mixed_crlf"}
TOOL_ORDER = [
    "z-fasta-fai",
    "z-fasta-fai-nodedup",
    "z-fasta-zfi",
    "z-fasta-zfi-nodedup",
    "samtools",
    "seqkit",
    "fastahack",
    "pyfaidx",
    "noodles",
    "rustbio-custom-index",
]
DATASET_ORDER = ["Genome", "Transcriptome", "Proteome"]
SECTION_WORKLOADS = {
    "real": ("Genome", "Transcriptome", "Proteome"),
    "scale_size": (
        "1mb",
        "5mb",
        "10mb",
        "25mb",
        "50mb",
        "100mb",
        "250mb",
        "500mb",
        "1000mb",
    ),
    "scale_seqs_budget": ("1000", "10000", "100000", "250000"),
    "scale_seqs_fixed": ("100000", "250000", "500000", "1000000"),
}

# Pinned versions for vendored or wrapper tools (keep in sync with tools.sh /
# tools/*/Cargo.toml / tools/fastahack-*). install_tools.sh verifies these.
VERSION_PINS = {
    "pyfaidx": "0.9.0.3",
    "seqkit": "2.13.0",
    "fastahack": "1.0.0",
    "noodles": "0.61",
    "rustbio": "2.2",
}

# Competitor metadata for provenance text (command + version source).
COMPETITOR_TOOLS = [
    {
        "key": "samtools",
        "label": "samtools",
        "language": "C",
        "command": "samtools faidx",
        "version_from": "`samtools --version`",
        "pin": None,
    },
    {
        "key": "seqkit",
        "label": "seqkit",
        "language": "Go",
        "command": "seqkit faidx",
        "version_from": "`seqkit version`",
        "pin": VERSION_PINS["seqkit"],
    },
    {
        "key": "fastahack",
        "label": "fastahack",
        "language": "C++",
        "command": "fastahack -i",
        "version_from": "directory pin `tools/fastahack-1.0.0/`",
        "pin": VERSION_PINS["fastahack"],
    },
    {
        "key": "pyfaidx",
        "label": "pyfaidx",
        "language": "Python",
        "command": "faidx --no-output",
        "version_from": "`faidx --version` / `pyfaidx.__version__`",
        "pin": VERSION_PINS["pyfaidx"],
    },
    {
        "key": "noodles",
        "label": "noodles",
        "language": "Rust",
        "command": "tools/noodles_wrapper index",
        "version_from": "noodles-fasta crate in `tools/noodles_wrapper/Cargo.toml`",
        "pin": VERSION_PINS["noodles"],
    },
    {
        "key": "rustbio",
        "label": "rust-bio",
        "language": "Rust",
        "command": "tools/rustbio_wrapper index",
        "version_from": "bio crate in `tools/rustbio_wrapper/Cargo.toml`",
        "pin": VERSION_PINS["rustbio"],
    },
]

# Real benchmark datasets (bench/shared/data; see bench/shared/download_data.sh).
REAL_DATASETS = [
    {
        "id": "Genome",
        "file": "REAL_Genome.fa",
        "description": "Homo sapiens GRCh38 primary assembly (Ensembl release 113)",
        "entry_noun": "sequences",
        "len_unit": "bp",
    },
    {
        "id": "Transcriptome",
        "file": "REAL_Transcriptome.fa",
        "description": "Homo sapiens GENCODE v46 transcript sequences",
        "entry_noun": "transcripts",
        "len_unit": "bp",
    },
    {
        "id": "Proteome",
        "file": "REAL_Proteome.fasta",
        "description": "Homo sapiens UniProt reference proteome UP000005640",
        "entry_noun": "proteins",
        "len_unit": "residues",
    },
]

# Product section: format and dedup matrix plus headline competitor and industry reference.
PRODUCT_COMPARISON_TOOLS = [
    "z-fasta-fai",
    "z-fasta-fai-nodedup",
    "z-fasta-zfi",
    "z-fasta-zfi-nodedup",
    "noodles",
    "samtools",
]

# Headline real-dataset chart: peer-comparable z-fasta (.fai) plus competitors.
HEADLINE_PERF_TOOLS = [
    "z-fasta-fai",
    "noodles",
    "rustbio-custom-index",
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


def dataset_sort_key(name):
    try:
        return DATASET_ORDER.index(name), str(name)
    except ValueError:
        return len(DATASET_ORDER), str(name)


def display_tool(name):
    return DISPLAY_NAMES.get(name, name)


def edge_exit_columns(df: pd.DataFrame) -> list[str]:
    return [column for column in EDGE_EXIT_COLUMNS if column in df.columns]


def edge_exit_label(column: str) -> str:
    return EDGE_EXIT_LABELS.get(column, column)


def edge_contract_classes(df: pd.DataFrame) -> pd.Series:
    if "contract_class" in df.columns:
        return df["contract_class"]

    def classify(case: str) -> str:
        if case == "binary_data":
            return "invalid_input"
        if case in EDGE_MESSY_CASES and case != "uniform":
            return "zfi_messy"
        return "fai_parity"

    return df["test_case"].map(classify)


def pivot_with_display_names(pivot: pd.DataFrame) -> pd.DataFrame:
    out = pivot.copy()
    out.columns = [display_tool(str(c)) for c in out.columns]
    return out


def filter_tools(df: pd.DataFrame, tools: list[str]) -> pd.DataFrame:
    return df[df["tool"].isin(tools)].copy()


def _ordered_pivot(
    pivot: pd.DataFrame,
    tools: list[str],
    *,
    sort_index=None,
) -> pd.DataFrame:
    cols = [c for c in tools if c in pivot.columns]
    pivot = pivot[cols]
    if sort_index is not None:
        pivot = pivot.reindex(sorted(pivot.index, key=sort_index))
    else:
        pivot = pivot.sort_index()
    return pivot_with_display_names(pivot)


def md_tool_pivot_table(
    df: pd.DataFrame,
    tools: list[str],
    index_col: str,
    value_col: str,
    fmt,
    *,
    sort_index=dataset_sort_key,
    empty_msg: str = "_No data._",
    require_non_null: bool = False,
) -> str:
    work = filter_tools(df, tools)
    if require_non_null and value_col in work.columns:
        work = work[work[value_col].notna()].copy()
    if work.empty:
        return empty_msg
    work = work.copy()
    work["cell"] = work.apply(fmt, axis=1)
    pivot = work.pivot(index=index_col, columns="tool", values="cell")
    return _ordered_pivot(pivot, tools, sort_index=sort_index).to_markdown()


# ══════════════════════════════════════════════════════════════════════
#  Data Loading
# ══════════════════════════════════════════════════════════════════════


def load_metadata(results_dir: Path, manifest: dict | None = None) -> pd.DataFrame | None:
    """Load metadata JSONL for the active run (manifest-linked when available)."""
    paths: list[Path] = []
    if manifest and manifest.get("metadata"):
        candidate = results_dir / manifest["metadata"]
        if candidate.exists():
            paths = [candidate]
    if not paths:
        latest = results_dir / "LATEST"
        if latest.exists():
            ts = latest.read_text().strip()
            candidate = results_dir / f"metadata_{ts}.jsonl"
            if candidate.exists():
                paths = [candidate]
    if not paths:
        found = sorted(results_dir.glob("metadata_*.jsonl"), reverse=True)
        if found:
            paths = [found[0]]

    rows = []
    for path in paths:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                row["_metadata_file"] = str(path)
                rows.append(row)
    if not rows:
        return None
    return pd.DataFrame(rows)


def meta_rows_for_json(metadata_df: pd.DataFrame | None, path: Path) -> list[dict]:
    """Rows from metadata JSONL whose raw_json matches path."""
    if metadata_df is None or metadata_df.empty or "raw_json" not in metadata_df.columns:
        return []
    resolved = path.resolve()
    rows = []
    for _, row in metadata_df.iterrows():
        raw = Path(row["raw_json"])
        try:
            if raw.resolve() == resolved:
                rows.append(row.to_dict())
        except OSError:
            if raw == path:
                rows.append(row.to_dict())
    return rows


def load_zebrac_json(path: Path, metadata_df: pd.DataFrame | None = None) -> pd.DataFrame:
    with open(path) as f:
        data = json.load(f)

    meta_rows = meta_rows_for_json(metadata_df, path)
    meta_by_command = {row.get("command"): row for row in meta_rows}
    rows = []
    for idx, result in enumerate(data.get("results", [])):
        command = result.get("command", "")
        meta = meta_by_command.get(command, {})
        if not meta and idx < len(meta_rows):
            meta = meta_rows[idx]
        wall = result.get("wall_time", {})
        peak_rss = result.get("peak_rss", {})
        mean_s = wall.get("mean", 0) / 1_000_000_000.0
        input_bytes = meta.get("input_bytes")
        throughput_mibs = None
        input_mib = None
        if input_bytes:
            input_mib = float(input_bytes) / (1024 * 1024)
            if mean_s > 0:
                throughput_mibs = input_mib / mean_s
        rows.append(
            {
                "tool": meta.get("tool", command),
                "command": command,
                "metadata_schema": meta.get("schema_version"),
                "metadata_suite": meta.get("suite"),
                "metadata_section": meta.get("section"),
                "metadata_workload": meta.get("workload"),
                "metadata_command": meta.get("command"),
                "mean": mean_s,
                "stddev": wall.get("std_dev", 0) / 1_000_000_000.0,
                "median": wall.get("median", 0) / 1_000_000_000.0,
                "min": wall.get("min", 0) / 1_000_000_000.0,
                "max": wall.get("max", 0) / 1_000_000_000.0,
                "runner": "zebrac",
                "sample_count": result.get("sample_count"),
                "failed_sample_count": result.get("failed_sample_count", 0),
                "peak_rss_mb": peak_rss.get("mean", 0) / (1024 * 1024),
                "peak_rss_stddev_mb": peak_rss.get("std_dev", 0) / (1024 * 1024),
                "minor_faults": result.get("minor_faults", {}).get("mean", 0),
                "major_faults": result.get("major_faults", {}).get("mean", 0),
                "cpu_cycles": result.get("cpu_cycles", {}).get("mean", 0),
                "instructions": result.get("instructions", {}).get("mean", 0),
                "cache_references": result.get("cache_references", {}).get("mean", 0),
                "cache_misses": result.get("cache_misses", {}).get("mean", 0),
                "branch_misses": result.get("branch_misses", {}).get("mean", 0),
                "input_bytes": input_bytes,
                "input_mib": input_mib,
                "throughput_mibs": throughput_mibs,
            }
        )
    return pd.DataFrame(rows)


def load_result_json(path: Path, metadata_df: pd.DataFrame | None = None) -> pd.DataFrame:
    with open(path) as f:
        data = json.load(f)
    if "zebrac_version" not in data:
        raise ValueError(f"expected zebrac JSON, got unsupported format: {path}")
    return load_zebrac_json(path, metadata_df)


def discover_latest(results_dir: Path, prefix: str, manifest: dict | None = None) -> Path | None:
    """Find timestamped directory or file for the active run.

    When a manifest is present, never fall back to a newer foreign run's artifacts.
    """
    section_key = {
        "perf": "real",
        "scale_size": "scale_size",
        "scale_seqs_budget": "scale_seqs_budget",
        "scale_seqs_fixed": "scale_seqs_fixed",
    }.get(prefix)
    if manifest is not None and section_key:
        rel = manifest.get("sections", {}).get(section_key)
        if not rel:
            return None
        candidate = results_dir / rel
        return candidate if candidate.exists() else None

    # Unmapped prefixes (e.g. tests_*.csv): pin to the run timestamp.
    if manifest is not None:
        ts = manifest.get("timestamp")
        if not ts:
            return None
        for candidate in (
            results_dir / f"{prefix}_{ts}",
            results_dir / f"{prefix}_{ts}.csv",
        ):
            if candidate.is_dir() or candidate.is_file():
                return candidate
        return None

    dirs = sorted(results_dir.glob(f"{prefix}_*"), reverse=True)
    dirs = [d for d in dirs if d.is_dir()]
    if dirs:
        return dirs[0]

    csvs = sorted(results_dir.glob(f"{prefix}_*.csv"), reverse=True)
    if csvs:
        return csvs[0]
    return None


def load_run_manifest(results_dir: Path) -> dict | None:
    """Load run_<ts>.json via LATEST, else newest run_*.json."""
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


def is_incomplete_run(
    manifest: dict | None,
    *,
    required_sections: tuple[str, ...] = (),
    skip_flags: tuple[str, ...] = ("skip_real", "skip_scaling", "skip_size"),
    min_warmup: int = 5,
    min_runs: int = 5,
) -> bool:
    """True when required sections are missing, skip flags are set, or samples are thin."""
    if not manifest:
        return True
    sections = manifest.get("sections") or {}
    for key in required_sections:
        if key not in sections:
            return True
    for flag in skip_flags:
        if manifest.get(flag):
            return True
    try:
        if int(manifest.get("warmup", 0)) < min_warmup:
            return True
        if int(manifest.get("runs", 0)) < min_runs:
            return True
    except (TypeError, ValueError):
        return True
    return False


def expected_index_cells(manifest: dict) -> set[tuple[str, str, str]]:
    """Exact section, workload, and tool cells required by this manifest."""
    tools = manifest.get("tools") or {}
    peers = ["samtools"]
    for key, tool_id in (
        ("seqkit", "seqkit"),
        ("fastahack", "fastahack"),
        ("pyfaidx", "pyfaidx"),
        ("noodles", "noodles"),
        ("rustbio", "rustbio-custom-index"),
    ):
        if tools.get(key):
            peers.append(tool_id)

    cells: set[tuple[str, str, str]] = set()
    for section, workloads in SECTION_WORKLOADS.items():
        if section not in (manifest.get("sections") or {}):
            continue
        tool_ids = (
            [
                "z-fasta-fai",
                "z-fasta-fai-nodedup",
                "z-fasta-zfi",
                "z-fasta-zfi-nodedup",
            ]
            if section == "real"
            else ["z-fasta-fai"]
        )
        tool_ids.extend(peers)
        cells.update((section, workload, tool) for workload in workloads for tool in tool_ids)
    return cells


def publication_issues(manifest: dict | None, results_dir: Path) -> list[str]:
    """Return reasons the active run cannot replace the tracked index report."""
    if not manifest:
        return ["missing run manifest"]

    issues: list[str] = []
    if manifest.get("schema_version") != "index-run.v2":
        issues.append("manifest is not index-run.v2")

    source = manifest.get("source") or {}
    if not source.get("git_commit") or not isinstance(source.get("dirty"), bool):
        issues.append("missing source identity")
    build = manifest.get("build") or {}
    for field in ("zig_version", "target", "optimize", "binary_path"):
        if not build.get(field):
            issues.append(f"missing build identity field: {field}")
    if manifest.get("cache_condition") != "warm":
        issues.append("cache condition is not warm")

    correctness = manifest.get("correctness") or {}
    correctness_artifact = correctness.get("artifact")
    correctness_checks = correctness.get("checks")
    correctness_reviews = correctness.get("reviews")
    if (
        correctness.get("status") != "pass"
        or correctness_checks != 24
        or not isinstance(correctness_reviews, int)
    ):
        issues.append("index correctness did not pass 24 checks")
    elif not correctness_artifact or not (results_dir / correctness_artifact).is_file():
        issues.append("index correctness artifact is missing")
    else:
        try:
            correctness_rows = pd.read_csv(results_dir / correctness_artifact)
            outcomes = correctness_rows.get("output_match")
        except (OSError, pd.errors.ParserError):
            outcomes = None
        if (
            outcomes is None
            or len(outcomes) != correctness_checks + correctness_reviews
            or int((outcomes == "MATCH").sum()) != correctness_checks
            or int((outcomes == "REVIEW").sum()) != correctness_reviews
            or not bool(outcomes.isin(("MATCH", "REVIEW")).all())
        ):
            issues.append("index correctness artifact does not match its recorded result")

    expected = expected_index_cells(manifest)
    recorded_cells = manifest.get("expected_cells") or []
    try:
        recorded = {
            (cell["section"], str(cell["workload"]), cell["tool"])
            for cell in recorded_cells
        }
    except (KeyError, TypeError):
        recorded = set()
    if len(recorded_cells) != len(recorded) or recorded != expected:
        issues.append("manifest expected-cell set does not match the product matrix")

    metadata_name = manifest.get("metadata")
    metadata_path = results_dir / metadata_name if isinstance(metadata_name, str) else None
    if metadata_path is None or not metadata_path.is_file():
        metadata = None
        issues.append("manifest-linked benchmark metadata is missing")
    else:
        try:
            metadata = load_metadata(results_dir, manifest)
        except (OSError, ValueError, json.JSONDecodeError):
            metadata = None
            issues.append("benchmark metadata is unreadable")
    actual: set[tuple[str, str, str]] = set()
    try:
        minimum_samples = int(manifest.get("runs", 0) or 0)
    except (TypeError, ValueError):
        minimum_samples = 0
        issues.append("measured sample requirement is invalid")
    for section, workloads in SECTION_WORKLOADS.items():
        rel = (manifest.get("sections") or {}).get(section)
        if not rel:
            continue
        section_dir = results_dir / rel
        expected_files = {f"{workload}.json" for workload in workloads}
        actual_files = (
            {path.name for path in section_dir.glob("*.json")}
            if section_dir.is_dir()
            else set()
        )
        if actual_files != expected_files:
            issues.append(f"{section} result files do not match the required workloads")
            continue
        for workload in workloads:
            result_path = section_dir / f"{workload}.json"
            try:
                frame = load_result_json(result_path, metadata)
            except (OSError, ValueError, json.JSONDecodeError):
                issues.append(f"unreadable result file: {section}/{workload}")
                continue
            if frame.empty or "tool" not in frame.columns:
                issues.append(f"empty result file: {section}/{workload}")
                continue
            metadata_rows = meta_rows_for_json(metadata, result_path)
            metadata_matches = (
                len(metadata_rows) == len(frame)
                and bool((frame["metadata_schema"] == "bench.meta.v1").all())
                and bool((frame["metadata_suite"] == "index").all())
                and bool((frame["metadata_section"] == section).all())
                and bool((frame["metadata_workload"].astype(str) == workload).all())
                and bool((frame["metadata_command"] == frame["command"]).all())
            )
            if not metadata_matches:
                issues.append(f"metadata does not match result cell: {section}/{workload}")
            if frame["tool"].duplicated().any():
                issues.append(f"duplicate measured cells: {section}/{workload}")
            for row in frame.itertuples(index=False):
                actual.add((section, workload, row.tool))
                sample_count = row.sample_count
                failed_count = row.failed_sample_count
                try:
                    samples_valid = (
                        not pd.isna(sample_count)
                        and int(sample_count) >= minimum_samples
                    )
                except (TypeError, ValueError, OverflowError):
                    samples_valid = False
                try:
                    failures_valid = not pd.isna(failed_count) and int(failed_count) == 0
                except (TypeError, ValueError, OverflowError):
                    failures_valid = False
                if not samples_valid:
                    issues.append(f"insufficient successful samples: {section}/{workload}/{row.tool}")
                if not failures_valid:
                    issues.append(f"failed samples present: {section}/{workload}/{row.tool}")
    if actual != expected:
        issues.append("measured cells do not match the manifest product matrix")
    return issues


def load_section_frames(
    results_dir: Path,
    manifest: dict | None,
    key: str,
    *,
    load_json,
    load_metadata,
) -> pd.DataFrame | None:
    """Concat zebrac JSON frames for one manifest section key."""
    rel = (manifest or {}).get("sections", {}).get(key)
    if not rel:
        return None
    section_dir = results_dir / rel
    if not section_dir.is_dir():
        return None
    metadata = load_metadata(results_dir, manifest)
    frames = []
    for jf in sorted(section_dir.glob("*.json")):
        frame = load_json(jf, metadata)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def prune_stale_pngs(figures_dir: Path, allowed) -> list[str]:
    """Delete PNGs not in allowed; return removed filenames."""
    allowed_set = set(allowed)
    removed: list[str] = []
    if not figures_dir.is_dir():
        return removed
    for path in figures_dir.glob("*.png"):
        if path.name in allowed_set:
            continue
        path.unlink()
        removed.append(path.name)
    return removed


def enrich_manifest(manifest: dict | None, results_dir: Path) -> dict | None:
    """Fill missing provenance fields and resolve tool versions from live binaries."""
    if not manifest:
        return manifest

    manifest.setdefault("runner", "zebrac")
    manifest.setdefault("mode", "warm")
    if manifest.get("schema_version") == "index-run.v2":
        return manifest

    project_root = results_dir.parent.parent.parent
    versions = resolve_tool_versions(manifest, project_root)

    if versions.get("z-fasta"):
        manifest["z_fasta"] = versions["z-fasta"]
    if versions.get("zebrac"):
        manifest["zebrac"] = versions["zebrac"]

    tools = dict(manifest.get("tools") or {})
    for key in ("samtools", "seqkit", "fastahack", "pyfaidx", "noodles", "rustbio"):
        if versions.get(key):
            tools[key] = versions[key]
    manifest["tools"] = tools
    return manifest


def command_first_line(cmd: list[str]) -> str | None:
    try:
        out = subprocess.run(cmd, check=True, capture_output=True)
        text = (out.stdout or out.stderr or b"").decode(errors="replace").strip()
        if not text:
            return None
        return text.splitlines()[0]
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def probe_live_versions(project_root: Path) -> dict[str, str]:
    """Query installed benchmark binaries for version strings."""
    tools_dir = project_root / "tools"
    versions: dict[str, str] = {}

    zfasta = project_root / "zig-out" / "bin" / "z-fasta"
    if zfasta.is_file():
        line = command_first_line([str(zfasta), "--version"])
        if line:
            versions["z-fasta"] = line

    zebrac = tools_dir / "zebrac"
    if zebrac.is_file():
        line = command_first_line([str(zebrac), "--version"])
        if line:
            versions["zebrac"] = line

    line = command_first_line(["samtools", "--version"])
    if line:
        versions["samtools"] = line

    seqkit = tools_dir / "seqkit"
    if seqkit.is_file():
        line = command_first_line([str(seqkit), "version"])
        if line:
            versions["seqkit"] = line

    for faidx in (project_root / ".venv" / "bin" / "faidx",):
        candidates = [faidx]
        system_faidx = shutil.which("faidx")
        if system_faidx:
            candidates.append(Path(system_faidx))
        for candidate in candidates:
            if candidate.is_file():
                line = command_first_line([str(candidate), "--version"])
                if line:
                    versions["pyfaidx"] = line
                    break
        if "pyfaidx" in versions:
            break

    fastahack = tools_dir / "fastahack-1.0.0" / "fastahack"
    if fastahack.is_file():
        versions["fastahack"] = (
            f"fastahack {VERSION_PINS['fastahack']} (directory pin)"
        )

    noodles = tools_dir / "noodles_wrapper" / "target" / "release" / "noodles_wrapper"
    if noodles.is_file():
        versions["noodles"] = f"noodles-fasta {VERSION_PINS['noodles']} (wrapper)"

    rustbio = tools_dir / "rustbio_wrapper" / "target" / "release" / "rustbio_wrapper"
    if rustbio.is_file():
        versions["rustbio"] = (
            f"rust-bio {VERSION_PINS['rustbio']} (custom indexer wrapper)"
        )

    return versions


def resolve_tool_versions(manifest: dict | None, project_root: Path) -> dict[str, str]:
    """Merge manifest-captured versions with live probes."""
    merged = probe_live_versions(project_root)
    if not manifest:
        return merged

    for key, val in (manifest.get("tools") or {}).items():
        if val:
            merged[key] = val
    if manifest.get("z_fasta"):
        merged["z-fasta"] = manifest["z_fasta"]
    if manifest.get("zebrac"):
        merged["zebrac"] = manifest["zebrac"]
    return merged


def load_json_suite(
    results_dir: Path,
    prefix: str,
    metadata_df: pd.DataFrame | None,
    manifest: dict | None,
    *,
    param_col: str | None = None,
    dataset_from_stem: bool = False,
) -> pd.DataFrame | None:
    """Load a directory of zebrac JSON files into one DataFrame."""
    d = discover_latest(results_dir, prefix, manifest)
    if not d:
        return None

    frames = []
    for jf in sorted(d.glob("*.json")):
        frame = load_result_json(jf, metadata_df)
        if dataset_from_stem:
            frame["dataset"] = jf.stem
        elif param_col:
            stem = jf.stem
            try:
                val = float(stem.replace("mb", ""))
            except ValueError:
                val = 0.0
            frame[param_col] = val
        frames.append(frame)

    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def load_perf(
    results_dir: Path,
    metadata_df: pd.DataFrame | None = None,
    manifest: dict | None = None,
) -> pd.DataFrame | None:
    """Load real-dataset perf results (dataset column from filename stem)."""
    return load_json_suite(
        results_dir,
        "perf",
        metadata_df,
        manifest,
        dataset_from_stem=True,
    )


def load_scaling(
    results_dir: Path,
    prefix: str,
    param_col: str,
    metadata_df: pd.DataFrame | None = None,
    manifest: dict | None = None,
) -> pd.DataFrame | None:
    """Load scale_size or scale_seqs directory of JSON files."""
    return load_json_suite(
        results_dir,
        prefix,
        metadata_df,
        manifest,
        param_col=param_col,
    )


def load_tests(results_dir: Path, manifest: dict | None = None) -> pd.DataFrame | None:
    """Load edge-case test CSV."""
    p = discover_latest(results_dir, "tests", manifest)
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


def filter_headline_perf(df: pd.DataFrame) -> pd.DataFrame:
    return filter_tools(df, HEADLINE_PERF_TOOLS)


def filter_product_comparison(df: pd.DataFrame) -> pd.DataFrame:
    return filter_tools(df, PRODUCT_COMPARISON_TOOLS)


def _scaling_param_label(param_col: str, value: float) -> str:
    v = int(value)
    if param_col == "size_mb":
        return f"{v} MB"
    return f"{v:,}"


def build_ratio_comparisons(
    work: pd.DataFrame,
    value_col: str,
    *,
    baseline: str,
    peer_tools: list[str],
    group_col: str = "dataset",
    group_sort=None,
    label_group=None,
) -> pd.DataFrame:
    """Baseline vs each peer; ratio = peer value / baseline value.

    ``baseline`` is required: suites use different tool ids (index: ``z-fasta-fai``,
    get: ``z-fasta-default``, stats: ``z-fasta-full``).
    """
    if group_sort is None:
        group_sort = dataset_sort_key if group_col == "dataset" else None
    groups = sorted(work[group_col].unique(), key=group_sort) if group_sort else sorted(
        work[group_col].unique()
    )
    rows: list[dict] = []
    for group in groups:
        base = work[(work[group_col] == group) & (work["tool"] == baseline)]
        if base.empty:
            continue
        base_v = float(base[value_col].values[0])
        group_label = label_group(group) if label_group else group
        for tool in peer_tools:
            if tool == baseline:
                continue
            hit = work[(work[group_col] == group) & (work["tool"] == tool)]
            if hit.empty:
                continue
            peer_v = float(hit[value_col].values[0])
            ratio = peer_v / base_v if base_v > 0 else None
            rows.append(
                {
                    "dataset": group_label,
                    "tool": tool,
                    "competitor": display_tool(tool),
                    "zfasta_v": base_v,
                    "comp_v": peer_v,
                    "ratio": ratio,
                }
            )
    return pd.DataFrame(rows)


def build_headline_ratio_comparisons(work: pd.DataFrame, value_col: str) -> pd.DataFrame:
    return build_ratio_comparisons(
        work,
        value_col,
        baseline="z-fasta-fai",
        peer_tools=HEADLINE_PERF_TOOLS[1:],
    )


def build_product_ratio_comparisons(
    work: pd.DataFrame,
    value_col: str,
    *,
    baseline: str = "z-fasta-fai",
) -> pd.DataFrame:
    return build_ratio_comparisons(
        work,
        value_col,
        baseline=baseline,
        peer_tools=[t for t in PRODUCT_COMPARISON_TOOLS if t != baseline],
    )


def build_scaling_ratio_comparisons(
    work: pd.DataFrame,
    param_col: str,
    value_col: str = "mean",
) -> pd.DataFrame:
    return build_ratio_comparisons(
        work,
        value_col,
        baseline="z-fasta-fai",
        peer_tools=HEADLINE_PERF_TOOLS[1:],
        group_col=param_col,
        label_group=lambda p: _scaling_param_label(param_col, p),
    )


def build_time_throughput_comparisons(
    work: pd.DataFrame,
    ratio_df: pd.DataFrame,
    *,
    baseline: str,
) -> pd.DataFrame:
    """Wall-time ratio rows enriched with throughput columns."""
    if ratio_df.empty:
        return ratio_df
    rows: list[dict] = []
    for row in ratio_df.itertuples(index=False):
        zf = work[(work["dataset"] == row.dataset) & (work["tool"] == baseline)]
        comp = work[(work["dataset"] == row.dataset) & (work["tool"] == row.tool)]
        zf_tp = (
            float(zf["throughput_mibs"].values[0])
            if not zf.empty and pd.notna(zf["throughput_mibs"].values[0])
            else None
        )
        comp_tp = (
            float(comp["throughput_mibs"].values[0])
            if not comp.empty and pd.notna(comp["throughput_mibs"].values[0])
            else None
        )
        tp_ratio = zf_tp / comp_tp if zf_tp and comp_tp and comp_tp > 0 else None
        rows.append(
            {
                "dataset": row.dataset,
                "tool": row.tool,
                "competitor": row.competitor,
                "zfasta_s": row.zfasta_v,
                "comp_s": row.comp_v,
                "speedup": row.ratio,
                "ratio": row.ratio,
                "zfasta_mibs": zf_tp,
                "comp_mibs": comp_tp,
                "throughput_x": tp_ratio,
            }
        )
    return pd.DataFrame(rows)


def build_headline_comparisons(work: pd.DataFrame) -> pd.DataFrame:
    return build_time_throughput_comparisons(
        work,
        build_headline_ratio_comparisons(work, "mean"),
        baseline="z-fasta-fai",
    )


def build_product_comparisons(work: pd.DataFrame) -> pd.DataFrame:
    return build_time_throughput_comparisons(
        work,
        build_product_ratio_comparisons(work, "mean"),
        baseline="z-fasta-fai",
    )


def _format_speedup(ratio: float | None) -> str:
    if ratio is None or ratio <= 0:
        return "n/a"
    if ratio >= 100:
        return f"{ratio:.0f}×"
    if ratio >= 10:
        return f"{ratio:.1f}×"
    if ratio >= 2:
        return f"{ratio:.1f}×"
    if ratio >= 0.1:
        return f"{ratio:.2f}×"
    if ratio >= 0.01:
        return f"{ratio:.3f}×"
    return f"{ratio:.4f}×"


def _annotate_headline_comparisons(
    ax1,
    datasets: list[str],
    tools: list[str],
    width: float,
    bar_tops: dict[tuple[str, str], tuple[float, float]],
    comparisons: pd.DataFrame,
) -> None:
    """Speedup badges (z-fasta vs competitor) and per-dataset reference lines."""
    zf_color = COLORS.get("z-fasta-fai", "#F7A41D")

    for ds in datasets:
        zf_key = (ds, "z-fasta-fai")
        if zf_key not in bar_tops:
            continue
        zf_x, zf_y = bar_tops[zf_key]
        cluster_xs = [bar_tops[(ds, t)][0] for t in tools if (ds, t) in bar_tops]
        if not cluster_xs:
            continue
        span_lo = min(cluster_xs) - width * 0.65
        span_hi = max(cluster_xs) + width * 0.65
        ax1.hlines(
            zf_y,
            span_lo,
            span_hi,
            colors=zf_color,
            linestyles=(0, (4, 3)),
            linewidth=1.1,
            alpha=0.55,
            zorder=1,
        )
        ax1.annotate(
            "1×",
            (zf_x, zf_y),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=zf_color,
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor=zf_color,
                linewidth=1.1,
            ),
            zorder=5,
        )

    for row in comparisons.itertuples(index=False):
        key = (row.dataset, row.tool)
        if key not in bar_tops:
            continue
        xpos, mean_t = bar_tops[key]
        color = COLORS.get(row.tool, "#666666")
        ratio = getattr(row, "ratio", None)
        if ratio is None:
            ratio = getattr(row, "speedup", None)
        ax1.annotate(
            _format_speedup(ratio),
            (xpos, mean_t),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            linespacing=0.95,
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor=color,
                linewidth=1.1,
                alpha=0.96,
            ),
            zorder=5,
        )


def _fig_throughput_dual_axis(
    work: pd.DataFrame,
    out: Path,
    *,
    tools: list[str],
    comparisons: pd.DataFrame,
    title: str,
    fig_note: str,
    bar_width: float,
) -> Path:
    """Grouped wall-time bars (log) with throughput lines on a twin axis."""
    datasets = sorted(work["dataset"].unique(), key=dataset_sort_key)
    tools = [t for t in tools if t in work["tool"].unique()]

    fig, ax1 = plt.subplots(figsize=(14, 7.2))
    ax2 = ax1.twinx()

    width = min(bar_width, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for i, tool in enumerate(tools):
        color = COLORS.get(tool, "#888888")
        xs: list[float] = []
        tputs: list[float] = []
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
                linestyle="-",
                linewidth=1.8,
                markersize=6,
                markeredgecolor="white",
                markeredgewidth=0.7,
                label=display_tool(tool),
                zorder=3,
            )

    if not comparisons.empty:
        _annotate_headline_comparisons(ax1, datasets, tools, width, bar_tops, comparisons)

    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(datasets, fontsize=11)
    ax1.set_ylabel("Wall Time (s, log scale)", fontsize=10)
    ax2.set_ylabel("Throughput (MiB/s)", fontsize=10)
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
    return _save(fig, out)


def fig_performance_throughput(df: pd.DataFrame, out: Path) -> Path:
    """Grouped wall-time bars (left, log) with throughput lines/markers (right)."""
    work = filter_headline_perf(df)
    return _fig_throughput_dual_axis(
        work,
        out,
        tools=HEADLINE_PERF_TOOLS,
        comparisons=build_headline_comparisons(work),
        title="Real-Dataset Indexing: Wall Time & Throughput",
        fig_note=(
            "Bar labels: z-fasta speedup vs each tool "
            "(competitor wall time ÷ z-fasta wall time)"
        ),
        bar_width=0.095,
    )


def fig_scaling_headline_lines(
    df: pd.DataFrame,
    out: Path,
    *,
    param_col: str,
    xlabel: str,
    title: str,
    fig_note: str,
) -> Path:
    """Line chart: wall time vs scaling parameter (headline tools, no ratio badges)."""
    work = filter_headline_perf(df)
    tools = [t for t in HEADLINE_PERF_TOOLS if t in work["tool"].unique()]

    fig, ax = plt.subplots(figsize=(14, 7.2))

    for tool in tools:
        tdf = work[work["tool"] == tool].sort_values(param_col)
        if tdf.empty:
            continue
        xs = tdf[param_col].astype(float)
        ys = tdf["mean"].astype(float).clip(lower=1e-6)
        color = COLORS.get(tool, "#888888")
        linewidth = 2.4 if tool == "z-fasta-fai" else 1.8
        ax.plot(
            xs,
            ys,
            color=color,
            marker="o",
            linestyle="-",
            linewidth=linewidth,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=display_tool(tool),
            zorder=3 if tool == "z-fasta-fai" else 2,
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
        ncol=min(len(handles), 8),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.92)
    return _save(fig, out)


def _fig_metric_bars(
    work: pd.DataFrame,
    out: Path,
    *,
    tools: list[str],
    comparisons: pd.DataFrame,
    value_col: str,
    std_col: str | None,
    ylabel: str,
    title: str,
    value_floor: float = 1e-3,
    fig_note: str | None = None,
    after_bars=None,
    bar_width: float = 0.095,
    legend_ncol: int | None = None,
) -> Path:
    """Grouped bars with z-fasta ratio badges."""
    datasets = sorted(work["dataset"].unique(), key=dataset_sort_key)
    tools = [t for t in tools if t in work["tool"].unique()]

    fig, ax = plt.subplots(figsize=(14, 7.2))
    width = min(bar_width, 0.80 / max(1, len(tools)))
    x = list(range(len(datasets)))
    bar_tops: dict[tuple[str, str], tuple[float, float]] = {}

    for i, tool in enumerate(tools):
        color = COLORS.get(tool, "#888888")
        for di, ds in enumerate(datasets):
            row = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if row.empty:
                continue
            val = max(float(row[value_col].values[0]), value_floor)
            std = 0.0
            if std_col and std_col in row.columns:
                raw_std = row[std_col].values[0]
                if pd.notna(raw_std):
                    std = max(float(raw_std), 0.0)
            xpos = di + (i - len(tools) / 2 + 0.5) * width
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
            bar_tops[(ds, tool)] = (xpos, val)

    if after_bars is not None:
        after_bars(ax, work, datasets, tools, width)

    if not comparisons.empty:
        _annotate_headline_comparisons(ax, datasets, tools, width, bar_tops, comparisons)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=20)
    if fig_note:
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
    ax.grid(axis="y", alpha=0.28, which="both")
    ax.set_axisbelow(True)

    patches = [
        mpatches.Patch(color=COLORS.get(tool, "#888888"), label=display_tool(tool))
        for tool in tools
    ]
    extra_handles, extra_labels = ax.get_legend_handles_labels()
    handles = patches + list(extra_handles)
    labels = [display_tool(t) for t in tools] + list(extra_labels)
    ncol = legend_ncol if legend_ncol is not None else min(len(handles), 8)
    fig.legend(
        handles=handles,
        labels=labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.15),
        ncol=ncol,
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(bottom=0.20, top=0.92)
    return _save(fig, out)


def _fig_headline_metric_bars(
    df: pd.DataFrame,
    out: Path,
    *,
    value_col: str,
    std_col: str | None,
    ylabel: str,
    title: str,
    value_floor: float = 1e-3,
    fig_note: str | None = None,
    after_bars=None,
) -> Path:
    """Grouped bars for headline tools across real datasets."""
    work = filter_headline_perf(df)
    return _fig_metric_bars(
        work,
        out,
        tools=HEADLINE_PERF_TOOLS,
        comparisons=build_headline_ratio_comparisons(work, value_col),
        value_col=value_col,
        std_col=std_col,
        ylabel=ylabel,
        title=title,
        value_floor=value_floor,
        fig_note=fig_note,
        after_bars=after_bars,
    )


def _draw_input_file_reference(
    ax,
    work: pd.DataFrame,
    datasets: list[str],
    tools: list[str],
    width: float,
) -> None:
    """Dashed reference at FASTA file size (MiB) per dataset cluster."""
    for di, ds in enumerate(datasets):
        ref = work[(work["dataset"] == ds) & (work["tool"] == "z-fasta-fai")]
        if ref.empty or "input_mib" not in ref.columns:
            continue
        input_mib = ref["input_mib"].values[0]
        if pd.isna(input_mib) or float(input_mib) <= 0:
            continue
        input_mb = float(input_mib)
        cluster_xs = [di + (i - len(tools) / 2 + 0.5) * width for i in range(len(tools))]
        span_lo = min(cluster_xs) - width * 0.65
        span_hi = max(cluster_xs) + width * 0.65
        ax.hlines(
            input_mb,
            span_lo,
            span_hi,
            colors="#666666",
            linestyles=(0, (5, 4)),
            linewidth=1.0,
            alpha=0.65,
            zorder=1,
        )
    ax.plot([], [], color="#666666", linestyle="--", linewidth=1.0, label="input file size")


def fig_memory_rss(df: pd.DataFrame, out: Path) -> Path:
    """Peak RSS for headline tools across real datasets."""
    return _fig_headline_metric_bars(
        df,
        out,
        value_col="peak_rss_mb",
        std_col="peak_rss_stddev_mb",
        ylabel="Peak RSS (MB, log scale)",
        title="Real-Dataset Indexing: Peak RSS",
        value_floor=1e-3,
        after_bars=_draw_input_file_reference,
        fig_note=(
            "Bar labels: competitor RSS / z-fasta RSS. "
            "Dashed grey = FASTA file size on disk (MiB)."
        ),
    )


def fig_page_faults(df: pd.DataFrame, out: Path) -> Path:
    """Minor page faults for headline tools across real datasets."""
    return _fig_headline_metric_bars(
        df,
        out,
        value_col="minor_faults",
        std_col=None,
        ylabel="Minor Page Faults (log scale)",
        title="Real-Dataset Indexing: Minor Page Faults",
        value_floor=1.0,
        fig_note="Bar labels: competitor faults / z-fasta faults.",
    )


def fig_product_comparison_throughput(df: pd.DataFrame, out: Path) -> Path | None:
    """Grouped wall-time bars with throughput lines for product lanes."""
    work = filter_product_comparison(df)
    if work.empty:
        return None
    return _fig_throughput_dual_axis(
        work,
        out,
        tools=PRODUCT_COMPARISON_TOOLS,
        comparisons=build_product_comparisons(work),
        title="Format and Dedup Comparison: Wall Time & Throughput",
        fig_note="Bar labels: wall time ratio vs z-fasta (.fai) (each value ÷ baseline)",
        bar_width=0.11,
    )


def fig_product_comparison_memory(df: pd.DataFrame, out: Path) -> Path | None:
    """Peak RSS for format and dedup product lanes."""
    work = filter_product_comparison(df)
    if work.empty:
        return None
    return _fig_metric_bars(
        work,
        out,
        tools=PRODUCT_COMPARISON_TOOLS,
        comparisons=build_product_ratio_comparisons(work, "peak_rss_mb"),
        value_col="peak_rss_mb",
        std_col="peak_rss_stddev_mb",
        ylabel="Peak RSS (MB, log scale)",
        title="Format and Dedup Comparison: Peak RSS",
        value_floor=1e-3,
        fig_note=(
            "Bar labels: peak RSS ratio vs z-fasta (.fai). "
            "Dashed grey = FASTA file size on disk."
        ),
        after_bars=_draw_input_file_reference,
        bar_width=0.11,
        legend_ncol=6,
    )


def _dataset_facet_title(ds: str, stats_by_id: dict[str, dict]) -> str:
    """Facet title with on-disk size and record count when available."""
    entry = stats_by_id.get(ds)
    if not entry or not entry.get("present"):
        return ds
    size = format_file_size(entry["file_bytes"])
    records = entry.get("records")
    if records is not None:
        return f"{ds}\n{size} · {records:,} {entry['entry_noun']}"
    return f"{ds}\n{size}"


def fig_product_comparison_tradeoff(
    df: pd.DataFrame,
    out: Path,
    *,
    stats_by_id: dict[str, dict],
) -> Path | None:
    """Time vs peak RSS scatter, one facet per real dataset."""
    work = filter_product_comparison(df)
    if work.empty:
        return None

    datasets = sorted(work["dataset"].unique(), key=dataset_sort_key)
    tools = [t for t in PRODUCT_COMPARISON_TOOLS if t in work["tool"].unique()]

    fig, axes = plt.subplots(1, len(datasets), figsize=(5 * len(datasets), 5.2))
    if len(datasets) == 1:
        axes = [axes]

    legend_handles: list = []
    legend_labels: list[str] = []

    for ax, ds in zip(axes, datasets):
        sub = work[work["dataset"] == ds]
        for tool in tools:
            row = sub[sub["tool"] == tool]
            if row.empty:
                continue
            color = COLORS.get(tool, "#888")
            label = display_tool(tool)
            scatter = ax.scatter(
                max(float(row["peak_rss_mb"].values[0]), 1e-3),
                max(float(row["mean"].values[0]), 1e-6),
                s=110,
                color=color,
                edgecolors="white",
                linewidths=0.6,
                zorder=3,
            )
            if label not in legend_labels:
                legend_handles.append(
                    mpatches.Patch(color=color, label=label)
                )
                legend_labels.append(label)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(_dataset_facet_title(ds, stats_by_id), fontsize=10, fontweight="bold")
        ax.grid(alpha=0.3, which="both")

    axes[0].set_ylabel("Wall Time (s, log scale)", fontsize=10)
    for ax in axes:
        ax.set_xlabel("Peak RSS (MB, log scale)", fontsize=10)

    fig.suptitle(
        "Format and Dedup Comparison: Time vs Memory Tradeoff",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )

    middle = axes[1] if len(axes) > 1 else axes[0]
    box = middle.get_position()
    fig.legend(
        handles=legend_handles,
        labels=legend_labels,
        fontsize=9,
        loc="upper center",
        bbox_to_anchor=(box.x0 + box.width / 2, box.y0 - 0.14),
        ncol=len(legend_handles),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    fig.subplots_adjust(top=0.88, bottom=0.22, wspace=0.28)
    return _save(fig, out)


def enrich_counter_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add IPC, cache miss rate, and related derived metrics."""
    w = df.copy()
    cycles = pd.to_numeric(w.get("cpu_cycles"), errors="coerce")
    instr = pd.to_numeric(w.get("instructions"), errors="coerce")
    refs = pd.to_numeric(w.get("cache_references"), errors="coerce")
    misses = pd.to_numeric(w.get("cache_misses"), errors="coerce")
    branches = pd.to_numeric(w.get("branch_misses"), errors="coerce")
    if "input_mib" not in w.columns or w["input_mib"].isna().all():
        bytes_col = pd.to_numeric(w.get("input_bytes"), errors="coerce")
        w["input_mib"] = bytes_col / (1024 * 1024)
    mib = pd.to_numeric(w["input_mib"], errors="coerce")
    w["ipc"] = instr / cycles.where(cycles > 0)
    w["cache_miss_rate"] = misses / refs.where(refs > 0)
    w["branch_miss_rate"] = branches / instr.where(instr > 0)
    w["instructions_per_mib"] = instr / mib.where(mib > 0)
    return w


def fig_zebrac_instructions_per_mib(df: pd.DataFrame, out: Path) -> Path | None:
    """Instructions per MiB of input; total CPU work per byte indexed."""
    required = {"instructions", "input_bytes"}
    if not required.issubset(df.columns) and "input_mib" not in df.columns:
        return None
    work = enrich_counter_metrics(filter_headline_perf(df))
    if work["instructions_per_mib"].isna().all():
        return None
    return _fig_metric_bars(
        work,
        out,
        tools=HEADLINE_PERF_TOOLS,
        comparisons=build_headline_ratio_comparisons(work, "instructions_per_mib"),
        value_col="instructions_per_mib",
        std_col=None,
        ylabel="Instructions per MiB (log scale)",
        title="Real-Dataset Indexing: CPU Work per MiB",
        value_floor=1e3,
        fig_note=(
            "Bar labels: competitor instructions/MiB ÷ z-fasta. "
            "Higher = more CPU work than z-fasta. Lower bars are better."
        ),
    )


def fig_edge_heatmap(df: pd.DataFrame, out: Path) -> Path:
    """Raw FAI/ZFI exit codes plus the composite contract result."""
    cases = df["test_case"].tolist()
    tool_cols = edge_exit_columns(df)
    tool_labels = [edge_exit_label(column) for column in tool_cols]

    n_tools = len(tool_cols)
    match_x = n_tools + 0.3

    fig, ax = plt.subplots(
        figsize=(max(8, n_tools * 1.15 + 2), max(7, len(cases) * 0.34))
    )

    for i, case in enumerate(cases):
        yi = len(cases) - i - 1
        for j, col in enumerate(tool_cols):
            val = df.iloc[i][col]
            missing = pd.isna(val) or int(val) == 127
            clr = "#9E9E9E" if missing else "#4CAF50" if int(val) == 0 else "#F44336"
            ax.add_patch(
                plt.Rectangle((j, yi), 1, 1, facecolor=clr, edgecolor="white", lw=2)
            )
            ax.text(
                j + 0.5,
                yi + 0.5,
                "n/a" if missing else str(int(val)),
                ha="center",
                va="center",
                fontsize=7,
                color="white",
                fontweight="bold",
            )

        m = df.iloc[i]["output_match"]
        if m == "MATCH":
            clr, label = "#4CAF50", "Y"
        elif m == "REVIEW":
            clr, label = "#7E57C2", "R"
        else:
            clr, label = "#FF9800", "N"
        ax.add_patch(
            plt.Rectangle(
                (match_x, yi), 1, 1, facecolor=clr, edgecolor="white", lw=2
            )
        )
        ax.text(
            match_x + 0.5,
            yi + 0.5,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            fontweight="bold",
        )

    ax.set_xlim(0, match_x + 1.2)
    ax.set_ylim(0, len(cases))
    ax.set_xticks([j + 0.5 for j in range(n_tools)] + [match_x + 0.5])
    ax.set_xticklabels(
        tool_labels + ["Contract"],
        fontweight="bold",
        rotation=90,
        ha="center",
        va="top",
    )
    ax.set_yticks([len(cases) - i - 0.5 for i in range(len(cases))])
    ax.set_yticklabels(cases, fontsize=8)
    ax.set_title(
        "Index Edge Cases: raw exit codes and contract result",
        fontweight="bold",
        pad=12,
    )

    patches = [
        mpatches.Patch(color="#4CAF50", label="Exit 0 (accepted)"),
        mpatches.Patch(color="#F44336", label="Non-zero exit"),
        mpatches.Patch(color="#9E9E9E", label="Not run"),
        mpatches.Patch(color="#7E57C2", label="Review"),
        mpatches.Patch(color="#FF9800", label="Contract mismatch"),
    ]
    ax.legend(handles=patches, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
    fig.subplots_adjust(bottom=0.28, right=0.82)
    return _save(fig, out)


# ══════════════════════════════════════════════════════════════════════
#  Markdown Report
# ══════════════════════════════════════════════════════════════════════


def normalize_markdown(chunks: list[str]) -> str:
    """Join chunks with one blank line between blocks; avoid MD012 (multiple blanks)."""
    cleaned: list[str] = []
    for chunk in chunks:
        text = chunk.strip()
        if not text:
            continue
        while "\n\n\n" in text:
            text = text.replace("\n\n\n", "\n\n")
        cleaned.append(text)
    return "\n\n".join(cleaned) + "\n"


def format_file_size(num_bytes: int) -> str:
    gib = num_bytes / (1024**3)
    if gib >= 1.0:
        return f"{gib:.2f} GiB"
    mib = num_bytes / (1024**2)
    if mib >= 1.0:
        return f"{mib:.1f} MiB"
    return f"{num_bytes:,} B"


def format_total_length(num_bases: int, unit: str) -> str:
    if unit == "bp":
        if num_bases >= 1_000_000_000:
            return f"{num_bases / 1e9:.2f} Gbp"
        if num_bases >= 1_000_000:
            return f"{num_bases / 1e6:.0f} Mbp"
        return f"{num_bases:,} bp"
    if num_bases >= 1_000_000:
        return f"{num_bases / 1e6:.1f} M {unit}"
    return f"{num_bases:,} {unit}"


def parse_fai_stats(fai_path: Path) -> tuple[int, int]:
    records = 0
    total_len = 0
    with open(fai_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].isdigit():
                records += 1
                total_len += int(parts[1])
    return records, total_len


def parse_zfi_stats(zfi_path: Path) -> tuple[int, int]:
    """Return (record_count, file_bytes) from a ZFI header."""
    data = zfi_path.read_bytes()
    if len(data) < 16 or data[:3] != b"ZFI":
        return 0, len(data)
    count = int.from_bytes(data[4:8], "little")
    return count, len(data)


def load_index_format_stats(data_dir: Path) -> list[dict]:
    """On-disk `.zfi` vs `.fai` sizes for real benchmark fixtures."""
    rows: list[dict] = []
    for spec in REAL_DATASETS:
        fasta = data_dir / spec["file"]
        if not fasta.is_file():
            continue
        zfi_path = data_dir / f"{spec['file']}.zfi"
        fai_path = data_dir / f"{spec['file']}.fai"
        row: dict = {
            "dataset": spec["id"],
            "fasta_bytes": fasta.stat().st_size,
        }
        if zfi_path.is_file():
            zfi_records, zfi_bytes = parse_zfi_stats(zfi_path)
            row["zfi_bytes"] = zfi_bytes
            row["zfi_records"] = zfi_records
            row["zfi_bytes_per_rec"] = zfi_bytes / zfi_records if zfi_records else None
        if fai_path.is_file():
            fai_records, _ = parse_fai_stats(fai_path)
            fai_bytes = fai_path.stat().st_size
            row["fai_bytes"] = fai_bytes
            row["fai_records"] = fai_records
            row["fai_bytes_per_line"] = fai_bytes / fai_records if fai_records else None
        if row.get("zfi_bytes") and row.get("fai_bytes"):
            row["fai_over_zfi"] = row["fai_bytes"] / row["zfi_bytes"]
        rows.append(row)
    return rows


def load_real_dataset_stats(data_dir: Path) -> list[dict]:
    """Read on-disk sizes and FAI entry counts for real benchmark fixtures."""
    stats: list[dict] = []
    for spec in REAL_DATASETS:
        fasta = data_dir / spec["file"]
        entry: dict = dict(spec)
        entry["present"] = fasta.is_file()
        if not entry["present"]:
            stats.append(entry)
            continue
        entry["file_bytes"] = fasta.stat().st_size
        fai = data_dir / f"{spec['file']}.fai"
        if fai.is_file():
            records, total_len = parse_fai_stats(fai)
            entry["records"] = records
            entry["total_len"] = total_len
        stats.append(entry)
    return stats


def md_data_used(project_root: Path) -> str:
    data_dir = project_root / "bench" / "shared" / "data"
    lines = [
        "**Data used** (human reference files in `bench/shared/data/`; "
        "fetch with `bench/shared/download_data.sh`):",
    ]
    any_present = False
    for entry in load_real_dataset_stats(data_dir):
        if not entry.get("present"):
            lines.append(
                f"- **{entry['id']}** (`{entry['file']}`): {entry['description']} "
                "(file not found on this machine)"
            )
            continue
        any_present = True
        records = entry.get("records")
        total_len = entry.get("total_len")
        size = format_file_size(entry["file_bytes"])
        if records is not None and total_len is not None:
            length = format_total_length(total_len, entry["len_unit"])
            lines.append(
                f"- **{entry['id']}** (`{entry['file']}`): {entry['description']}. "
                f"{records:,} {entry['entry_noun']}, {length} total, {size} on disk."
            )
        else:
            lines.append(
                f"- **{entry['id']}** (`{entry['file']}`): {entry['description']}. "
                f"{size} on disk (index `.fai` not found; run indexing to count entries)."
            )
    if not any_present:
        lines.append(
            "_No real dataset files present. Run `bash bench/shared/download_data.sh`._"
        )
    return "\n".join(lines)


def md_overview(manifest: dict | None) -> str:
    """Opening narrative: what the index suite measures and what this report contains."""
    lines = [
        "This report compares `z-fasta index` to FASTA indexers that build or use "
        "FAI-style `.fai` files. Wall time, throughput, peak RSS, and hardware counters "
        "all come from the same zebrac samples for every tool, so the tables and charts "
        "are directly comparable.",
        (
            "**Why `.fai` in cross-tool tables:** peer tools (samtools, noodles, seqkit, …) "
            "emit text `.fai`. Headline **Performance**, **Memory**, **Page Faults**, and "
            "**Scaling** therefore time **`z-fasta (.fai)`** (`index --emit-fai`) so the "
            "comparison is format-fair. That lane is not the CLI default."
        ),
        (
            "**Why `.zfi` is first-class:** CLI default `index` writes binary `.zfi` (embedded "
            "names, optional side tables for messy wrap). That is the production path get and "
            "stats are built around. **z-fasta Format and Dedup Comparison** times the four "
            "product lanes: `.fai` vs `.zfi`, with and without `--no-dedup`, with noodles "
            "and samtools as references. All z-fasta lanes use the same bounded reader."
        ),
        "**Performance and scaling:** three human reference downloads (Genome, "
        "Transcriptome, Proteome; see **Data used** under Run Provenance), synthetic "
        "file-size scaling from 1 MB through 1 GB, two sequence-count sweeps (bounded "
        "~50 MiB total vs fixed 1 kbp per record through 1M sequences), and speedup "
        "vs peers on the `.fai` lane. Zebrac counter columns are for profiling, not "
        "headline marketing numbers.",
        "**Correctness:** edge-case fixtures (structurally invalid FASTA) and messy "
        "variants (ragged wrap, trailing spaces, CRLF exports). z-fasta indexes some "
        "messy inputs with a side-table `.zfi`; samtools, noodles, and rust-bio follow "
        "strict FAI line-width rules.",
        "**On-disk `.zfi`:** after Format and Dedup Comparison, index file sizes and "
        "load/integrity notes. Build timings for every product lane live in that comparison, "
        "not a second wall table.",
        "**Charts** highlight peer-comparable `z-fasta (.fai)` plus samtools, seqkit, "
        "noodles, and rust-bio. **Product tables** list all four z-fasta lanes when present.",
    ]
    if manifest and manifest.get("skip_real"):
        lines.append("_This run skipped real datasets (`skip_real=true`)._")
    if manifest and manifest.get("skip_scaling"):
        lines.append("_This run skipped sequence-count scaling (`skip_scaling=true`)._")
    if manifest and manifest.get("skip_size"):
        lines.append("_This run skipped file-size scaling (`skip_size=true`)._")
    return "\n\n".join(lines)


def md_run_provenance(manifest: dict | None, project_root: Path) -> str:
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
    source = manifest.get("source") or {}
    build = manifest.get("build") or {}
    correctness = manifest.get("correctness") or {}

    lines = [
        f"Run **`{ts}`** used **{runner}** in **{mode}** mode: {runs} measured samples, "
        f"{warmup} warmup passes, {duration} ms minimum per sample.",
    ]

    bullets = [
        (
            f"- **Subject:** {zfasta} (Zig; peer tables use `z-fasta (.fai)`; Format "
            "and Dedup Comparison covers `.fai` / `.zfi` × dedup / `--no-dedup`)"
        ),
        f"- **Runner:** {zebrac}",
        (
            f"- **Source:** commit `{source.get('git_commit', 'unknown')}`; "
            f"dirty={source.get('dirty', 'unknown')}"
        ),
        (
            f"- **Build:** Zig {build.get('zig_version', 'unknown')}; "
            f"target `{build.get('target', 'unknown')}`; "
            f"optimize `{build.get('optimize', 'unknown')}`; "
            f"binary `{build.get('binary_path', 'unknown')}`"
        ),
        (
            f"- **Correctness:** {correctness.get('status', 'unknown')}; "
            f"{correctness.get('checks', '?')} checks; artifact "
            f"`results/{correctness.get('artifact', 'unknown')}`"
        ),
    ]

    artifacts: list[str] = []
    if sections.get("real"):
        artifacts.append(f"`results/{sections['real']}/` real datasets")
    if sections.get("scale_size"):
        artifacts.append(f"`results/{sections['scale_size']}/` file-size scaling")
    if sections.get("scale_seqs_budget"):
        artifacts.append(
            f"`results/{sections['scale_seqs_budget']}/` sequence-count scaling (bounded bytes)"
        )
    if sections.get("scale_seqs_fixed"):
        artifacts.append(
            f"`results/{sections['scale_seqs_fixed']}/` sequence-count scaling (fixed length)"
        )
    if artifacts:
        bullets.append(
            f"- **Artifacts:** {'; '.join(artifacts)}; metadata in `results/{metadata}`"
        )

    blocks = [lines[0], "\n".join(bullets), md_data_used(project_root)]

    measured: list[str] = []
    for spec in COMPETITOR_TOOLS:
        ver = tools.get(spec["key"])
        if not ver:
            continue
        pin = spec["pin"]
        pin_text = f"; pin {pin}" if pin else ""
        measured.append(
            f"- **{spec['label']}** ({spec['language']}): {ver} via `{spec['command']}` "
            f"(from {spec['version_from']}{pin_text})"
        )

    if measured:
        blocks.append(
            "**Competitors in this run** (versions captured at benchmark time; "
            "vendored pins in `bench/shared/tools.sh` / `tools/`):\n"
            + "\n".join(measured)
        )
    else:
        blocks.append(
            "_Competitor binaries were not resolved on this machine. "
            "Re-run `bash bench/index/run.sh` to capture versions in the manifest._"
        )

    blocks.append(
        "Index-file cleanup runs inside each measured command because zebrac has no "
        "separate prepare hook. Every tool pays the same small reset overhead."
    )
    return "\n\n".join(blocks)


def _headline_pivot(pivot: pd.DataFrame) -> pd.DataFrame:
    return _ordered_pivot(pivot, HEADLINE_PERF_TOOLS, sort_index=dataset_sort_key)


def md_perf_headline_table(df: pd.DataFrame) -> str:
    """Wall-time table for headline tools (z-fasta default + competitors)."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "mean",
        lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
        empty_msg="_No performance data._",
    )


def md_throughput_headline_table(df: pd.DataFrame) -> str:
    """Throughput table for headline tools."""
    if "throughput_mibs" not in df.columns:
        return "_No throughput data._"
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "throughput_mibs",
        lambda row: f"{row['throughput_mibs']:.1f} MiB/s",
        empty_msg="_No throughput data._",
        require_non_null=True,
    )


def md_zfasta_vs_tools_table(df: pd.DataFrame) -> str:
    """z-fasta vs competitor comparisons (matches Figure 1 bar labels)."""
    work = filter_headline_perf(df)
    comparisons = build_headline_comparisons(work)
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Speedup",
        fmt_zf=lambda r: f"{r.zfasta_s:.4f}s",
        fmt_comp=lambda r: f"{r.comp_s:.4f}s",
        extra_columns=lambda r: {
            "z-fasta MiB/s": f"{r.zfasta_mibs:.1f}" if r.zfasta_mibs else "n/a",
            "Competitor MiB/s": f"{r.comp_mibs:.1f}" if r.comp_mibs else "n/a",
            "Throughput ×": _format_speedup(r.throughput_x),
        },
    )


def md_zfasta_vs_ratio_table(
    comparisons: pd.DataFrame,
    *,
    group_label: str = "Dataset",
    vs_label: str = "z-fasta vs",
    zf_label: str,
    comp_label: str,
    ratio_label: str,
    fmt_zf,
    fmt_comp,
    extra_columns=None,
) -> str:
    """Generic baseline vs peer table (matches bar-label ratios on figures)."""
    if comparisons.empty:
        return "_No comparison data._"
    rows: list[dict] = []
    for row in comparisons.itertuples(index=False):
        entry = {
            group_label: row.dataset,
            vs_label: row.competitor,
            zf_label: fmt_zf(row),
            comp_label: fmt_comp(row),
            ratio_label: _format_speedup(row.ratio),
        }
        if extra_columns:
            entry.update(extra_columns(row))
        rows.append(entry)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_zfasta_vs_rss_table(df: pd.DataFrame) -> str:
    work = filter_headline_perf(df)
    comparisons = build_headline_ratio_comparisons(work, "peak_rss_mb")
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="RSS ×",
        fmt_zf=lambda r: f"{r.zfasta_v:.2f} MB",
        fmt_comp=lambda r: f"{r.comp_v:.2f} MB",
    )


def md_zfasta_vs_faults_table(df: pd.DataFrame) -> str:
    work = filter_headline_perf(df)
    comparisons = build_headline_ratio_comparisons(work, "minor_faults")
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Faults ×",
        fmt_zf=lambda r: f"{int(r.zfasta_v):,}",
        fmt_comp=lambda r: f"{int(r.comp_v):,}",
    )


def md_performance_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Performance section body: captioned tables, figure, and reading notes."""
    throughput = md_throughput_headline_table(df)
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_cmp = nums.next_table()
    f_perf = nums.next_figure()
    blocks = [
        "Headline comparison on the three human reference datasets (see **Data used**). "
        "**z-fasta (.fai)** is `index --emit-fai` (bounded reader + dedup, FAI to stdout) so "
        "peers are format-fair. CLI default `index` writes `.zfi`; see **z-fasta Format and "
        "Dedup Comparison** and "
        "**z-fasta Production Index (.zfi)**.",
        f"**Table {t_wall}:** Zebrac wall time per dataset (seconds, mean ± stddev, warm cache). "
        "Lower is better. Tool order: z-fasta, noodles, rust-bio, samtools, seqkit, "
        "fastahack, pyfaidx.",
        md_perf_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_tp}:</strong> Input throughput (MiB/s) from metadata "
        f"<code>input_bytes</code> and wall time. Higher is better. Same tools and order "
        f"as Table {t_wall}.</summary>",
        "",
        throughput,
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "Speedup = competitor wall time ÷ z-fasta wall time (higher = z-fasta faster). "
        f"Same ratios appear as bar labels on Figure {f_perf}. Throughput × uses MiB/s.</summary>",
        "",
        md_zfasta_vs_tools_table(df),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_perf}:** Dual-axis chart for Table {t_wall} and Table {t_tp}. Bars show wall time "
        "(left axis, log scale). Lines and markers show throughput (right axis). "
        "Legend colors: gold = z-fasta; bronze = noodles; red-brown = rust-bio; grey = "
        "samtools; teal = seqkit; pink = fastahack; blue = pyfaidx. Bar labels show "
        f"z-fasta speedup vs each tool (see Table {t_cmp}); gold dashed line = z-fasta wall time "
        "per dataset.",
        f"![Figure {f_perf}: real-dataset wall time and throughput](results/figures/performance.png)",
        (
            f"**Reading Figure {f_perf}**\n"
            "- **Bars (left):** zebrac mean wall time. Error bars are one standard deviation.\n"
            "- **Lines and markers (right):** input MiB/s for the same tool at each dataset.\n"
            "- **Legend order:** z-fasta, noodles, rust-bio, samtools, seqkit, fastahack, pyfaidx.\n"
            "- **Bar labels:** `1×` on z-fasta (baseline); competitor labels = z-fasta speedup "
            f"(competitor time ÷ z-fasta time). Border color matches the tool. Dashed gold line "
            f"= z-fasta wall time for that dataset. Details in Table {t_cmp}.\n"
            "- **z-fasta only:** this section uses the peer-comparable `.fai` lane; format × "
            "dedup tradeoffs are in **z-fasta Format and Dedup Comparison** (before Edge Case Correctness)."
        ),
    ]
    return "\n\n".join(blocks)


def md_scaling_headline_table(
    df: pd.DataFrame, param_col: str, param_label: str
) -> str:
    """Wall-time pivot for headline tools on a scaling sweep."""
    work = filter_headline_perf(df)
    if work.empty:
        return "_No scaling data._"
    pivot = work.pivot_table(
        index=param_col, columns="tool", values="mean", aggfunc="first"
    )
    pivot = _ordered_pivot(pivot, HEADLINE_PERF_TOOLS)
    pivot.index = [_scaling_param_label(param_col, v) for v in pivot.index]
    pivot.index.name = param_label
    return pivot.to_markdown(floatfmt=".4f")


def md_zfasta_vs_scaling_table(
    df: pd.DataFrame, param_col: str, param_label: str
) -> str:
    work = filter_headline_perf(df)
    comparisons = build_scaling_ratio_comparisons(work, param_col, "mean")
    return md_zfasta_vs_ratio_table(
        comparisons,
        group_label=param_label,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Time ×",
        fmt_zf=lambda r: f"{r.zfasta_v:.4f}s",
        fmt_comp=lambda r: f"{r.comp_v:.4f}s",
    )


def md_scaling_section(
    df: pd.DataFrame,
    nums: ReportCounters,
    *,
    intro: str,
    subsection: str,
    param_col: str,
    param_label: str,
    figure_path: str,
    figure_title: str,
    table_x_label: str,
    reading_x_label: str,
) -> str:
    """Shared template for file-size and sequence-count scaling sections."""
    t_wall = nums.next_table()
    t_cmp = nums.next_table()
    f_scale = nums.next_figure()
    blocks = [
        intro,
        subsection,
        "<details>",
        f"<summary><strong>Table {t_wall}:</strong> Wall time per {table_x_label} "
        "(zebrac mean, seconds). Same tool order as Performance.</summary>",
        "",
        md_scaling_headline_table(df, param_col, param_label),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor at each "
        f"{table_x_label}. Time × = competitor wall time / z-fasta wall time.</summary>",
        "",
        md_zfasta_vs_scaling_table(df, param_col, param_label),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_scale}:** Table {t_wall} as lines (log scales on both axes). "
        "One line per tool; absolute wall times only.",
        f"![Figure {f_scale}: {figure_title}]({figure_path})",
        (
            f"**Reading Figure {f_scale}**\n"
            f"- Lines: mean wall time vs {reading_x_label}. Both axes use log scale.\n"
            f"- Lower on the chart = faster at that {table_x_label}.\n"
            "- Legend colors match **Performance: Real Datasets**.\n"
            f"- Time × ratios are in Table {t_cmp}; the chart does not show them."
        ),
    ]
    return "\n\n".join(blocks)


def md_product_wall_time_table(df: pd.DataFrame) -> str:
    return md_tool_pivot_table(
        df,
        PRODUCT_COMPARISON_TOOLS,
        "dataset",
        "mean",
        lambda row: f"{row['mean']:.4f}s ±{row['stddev']:.4f}",
        empty_msg="_No product comparison data._",
    )


def md_product_throughput_table(df: pd.DataFrame) -> str:
    if "throughput_mibs" not in df.columns:
        return "_No throughput data._"
    return md_tool_pivot_table(
        df,
        PRODUCT_COMPARISON_TOOLS,
        "dataset",
        "throughput_mibs",
        lambda row: f"{row['throughput_mibs']:.1f} MiB/s",
        empty_msg="_No throughput data._",
        require_non_null=True,
    )


def md_product_rss_table(df: pd.DataFrame) -> str:
    return md_tool_pivot_table(
        df,
        PRODUCT_COMPARISON_TOOLS,
        "dataset",
        "peak_rss_mb",
        lambda row: f"{row['peak_rss_mb']:.2f} MB",
        empty_msg="_No product comparison data._",
    )


def md_product_vs_default_table(
    df: pd.DataFrame,
    *,
    value_col: str,
    baseline_label: str,
    ratio_label: str,
    fmt_baseline,
    fmt_other,
) -> str:
    work = filter_product_comparison(df)
    comparisons = build_product_ratio_comparisons(work, value_col)
    return md_zfasta_vs_ratio_table(
        comparisons,
        vs_label="vs z-fasta (.fai)",
        zf_label=baseline_label,
        comp_label="Other",
        ratio_label=ratio_label,
        fmt_zf=fmt_baseline,
        fmt_comp=fmt_other,
    )


def md_product_time_vs_default_table(df: pd.DataFrame) -> str:
    work = filter_product_comparison(df)
    comparisons = build_product_comparisons(work)
    return md_zfasta_vs_ratio_table(
        comparisons,
        vs_label="vs z-fasta (.fai)",
        zf_label="z-fasta (.fai)",
        comp_label="Other",
        ratio_label="Time ×",
        fmt_zf=lambda r: f"{r.zfasta_s:.4f}s",
        fmt_comp=lambda r: f"{r.comp_s:.4f}s",
        extra_columns=lambda r: {
            "(.fai) MiB/s": (
                f"{r.zfasta_mibs:.1f}" if r.zfasta_mibs is not None else "n/a"
            ),
            "Other MiB/s": (
                f"{r.comp_mibs:.1f}" if r.comp_mibs is not None else "n/a"
            ),
            "Throughput ×": (
                f"{r.throughput_x:.2f}×"
                if r.throughput_x is not None and r.throughput_x > 0
                else "n/a"
            ),
        },
    )


def md_product_comparison_section(
    df: pd.DataFrame,
    nums: ReportCounters,
) -> str:
    """Format and dedup section: four z-fasta lanes plus noodles and samtools."""
    t_wall = nums.next_table()
    t_tp = nums.next_table()
    t_time_cmp = nums.next_table()
    t_rss = nums.next_table()
    t_rss_cmp = nums.next_table()
    f_time = nums.next_figure()
    f_mem = nums.next_figure()
    f_trade = nums.next_figure()

    blocks = [
        (
            "Compares every z-fasta index product lane on the same human reference files "
            "as **Performance: Real Datasets**, with **noodles** and **samtools** as "
            "references. Axes: **output** (`.fai` via `--emit-fai` vs `.zfi` via default "
            "`index`) and **dedup** (on vs `--no-dedup`). All four z-fasta commands use "
            "the same bounded reader; labels use `z-fasta (.fai)` / `z-fasta (.zfi)` plus "
            "`--no-dedup` (never `--emit-fai` in the short name). "
            "Time × / RSS × baseline is **`z-fasta (.fai)`** (peer-comparable lane)."
        ),
        (
            "### Wall time and throughput\n\n"
            "Mean zebrac wall time per dataset. **Lower is better.** Tool order follows "
            "Format and Dedup Comparison columns (four z-fasta lanes, then noodles, samtools). "
            "Bar labels and Time × use each value divided by `z-fasta (.fai)` "
            "(values below `1×` are faster than that baseline)."
        ),
        (
            f"**Table {t_wall}:** Wall time (seconds, mean ± stddev). Same tool order as "
            "the product comparison figures."
        ),
        md_product_wall_time_table(df),
        "<details>",
        (
            f"<summary><strong>Table {t_tp}:</strong> Input throughput (MiB/s) from metadata "
            f"<code>input_bytes</code> and wall time. Higher is better. Same tools and order "
            f"as Table {t_wall}.</summary>"
        ),
        "",
        md_product_throughput_table(df),
        "</details>",
        "<details>",
        (
            f"<summary><strong>Table {t_time_cmp}:</strong> Each configuration vs "
            "`z-fasta (.fai)`. Time × = other wall time ÷ baseline wall time. Throughput × = "
            f"baseline MiB/s ÷ other MiB/s. Same ratios as bar labels on Figure {f_time}.</summary>"
        ),
        "",
        md_product_time_vs_default_table(df),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        (
            f"**Figure {f_time}:** Dual-axis chart for Table {t_wall} and Table {t_tp}. Bars "
            "show wall time (left axis, log scale). Lines show throughput (right axis). "
            f"Bar labels show Time × vs `z-fasta (.fai)` (see Table {t_time_cmp}); gold dashed "
            "line = baseline wall time per dataset."
        ),
        (
            f"![Figure {f_time}: format and dedup comparison wall time and throughput]"
            "(results/figures/product_comparison_performance.png)"
        ),
        (
            f"**Reading Figure {f_time}**\n"
            "- **Bars (left):** zebrac mean wall time. Error bars are one standard deviation.\n"
            "- **Lines and markers (right):** input MiB/s for the same tool at each dataset.\n"
            "- **Legend:** four z-fasta lanes (`(.fai)` / `(.zfi)` × dedup / `--no-dedup`), "
            "then noodles, samtools.\n"
            "- **Bar labels:** `1×` on `z-fasta (.fai)`; other labels = Time × (other ÷ baseline).\n"
            f"- Details in Table {t_time_cmp}."
        ),
        (
            "### Peak RSS\n\n"
            "zebrac peak RSS when the indexer process exits (`ru_maxrss`). **Lower is better** "
            "for the same file. Grey dashed lines on the figure mark FASTA size on disk. "
            "Differences between z-fasta lanes reflect output and dedup state; input reading "
            "is identical."
        ),
        (
            f"**Table {t_rss}:** Peak RSS (MB, zebrac mean). Same tool order as Table {t_wall}."
        ),
        md_product_rss_table(df),
        "<details>",
        (
            f"<summary><strong>Table {t_rss_cmp}:</strong> Each configuration vs "
            "`z-fasta (.fai)`. RSS × = other peak RSS ÷ baseline peak RSS. Same ratios as "
            f"bar labels on Figure {f_mem}.</summary>"
        ),
        "",
        md_product_vs_default_table(
            df,
            value_col="peak_rss_mb",
            baseline_label="z-fasta (.fai)",
            ratio_label="RSS ×",
            fmt_baseline=lambda r: f"{r.zfasta_v:.2f} MB",
            fmt_other=lambda r: f"{r.comp_v:.2f} MB",
        ),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        (
            f"**Figure {f_mem}:** Table {t_rss} as grouped bars (log scale). Bar labels = RSS × "
            f"(see Table {t_rss_cmp}). Gold dashed line = `z-fasta (.fai)` RSS; grey dashed = "
            "FASTA file size on disk."
        ),
        (
            f"![Figure {f_mem}: format and dedup comparison peak RSS]"
            "(results/figures/product_comparison_memory.png)"
        ),
        (
            f"**Reading Figure {f_mem}**\n"
            "- **Legend:** same four z-fasta lanes plus noodles and samtools as Figure "
            f"{f_time}.\n"
            "- Bars: mean peak RSS. `1×` on `z-fasta (.fai)`; other labels = other ÷ baseline.\n"
            "- Gold dashed line: baseline RSS for that dataset.\n"
            "- Grey dashed line: file size from `input_bytes`.\n"
            f"- Details in Table {t_rss_cmp}."
        ),
        (
            "### Time vs memory tradeoff\n\n"
            "Each panel is one real dataset. Facet titles show on-disk size and record count "
            "from **Data used**. One point per tool; both axes use log scale. **Lower-left** "
            "is faster with less RAM; **upper-right** is slower with more RAM."
        ),
        (
            f"**Figure {f_trade}:** Scatter of wall time (y) vs peak RSS (x) for Format and "
            "Dedup Comparison tools. Shared legend sits below the Transcriptome panel."
        ),
        (
            f"![Figure {f_trade}: format and dedup comparison time vs memory]"
            "(results/figures/product_comparison_tradeoff.png)"
        ),
        (
            f"**Reading Figure {f_trade}**\n"
            "- One facet per dataset (Genome, Transcriptome, Proteome).\n"
            "- Facet title: dataset name, file size, and entry count when available.\n"
            f"- Colors match Figure {f_time}: gold family = `.fai` lanes; orange family = "
            "`.zfi` lanes; bronze = noodles; grey = samtools.\n"
            "- Legend is shared under the middle facet."
        ),
    ]
    return "\n\n".join(blocks)


def md_zfi_size_prose(index_stats: list[dict]) -> str:
    """Prose summary of on-disk `.zfi` vs `.fai` from benchmark fixtures."""
    if not index_stats:
        return (
            "No index files found under `bench/shared/data/`. Run "
            "`bench/shared/download_data.sh` to generate sidecars."
        )

    lines = [
        "Each uniform entry is a fixed 40-byte `IndexRecord` (`src/index_format.zig`). "
        "FAI lines grow with sequence name length plus five tab-separated numeric fields. "
        "Non-uniform (messy) records add side-table bytes on top of the 40-byte header.",
    ]
    by_id = {row["dataset"]: row for row in index_stats}

    tx = by_id.get("Transcriptome")
    if tx and tx.get("fai_over_zfi") and tx.get("zfi_records"):
        lines.append(
            f"On **Transcriptome** ({tx['zfi_records']:,} entries), `.zfi` is "
            f"{format_file_size(tx['zfi_bytes'])} vs `.fai` "
            f"{format_file_size(tx['fai_bytes'])} "
            f"(FAI / ZFI = {tx['fai_over_zfi']:.2f}x; about "
            f"{tx['zfi_bytes_per_rec']:.0f} B per record vs "
            f"{tx['fai_bytes_per_line']:.0f} B per FAI line). Long GENCODE names favor "
            "the binary layout."
        )

    genome = by_id.get("Genome")
    if genome and genome.get("fai_over_zfi") and genome.get("zfi_records"):
        if genome["fai_over_zfi"] < 1.0:
            lines.append(
                f"On **Genome** ({genome['zfi_records']:,} sequences), short chromosome names "
                f"make `.fai` smaller (FAI / ZFI = {genome['fai_over_zfi']:.2f}x; "
                f"~{genome['fai_bytes_per_line']:.0f} B per FAI line vs "
                f"~{genome['zfi_bytes_per_rec']:.0f} B per record). Fixed 40-byte records "
                "carry overhead when names are only a few characters."
            )
        else:
            lines.append(
                f"On **Genome** ({genome['zfi_records']:,} sequences), "
                f"FAI / ZFI = {genome['fai_over_zfi']:.2f}x."
            )

    proteome = by_id.get("Proteome")
    if proteome and proteome.get("fai_over_zfi") and proteome.get("zfi_records"):
        if abs(proteome["fai_over_zfi"] - 1.0) < 0.05:
            lines.append(
                f"On **Proteome** ({proteome['zfi_records']:,} entries), sizes are about "
                f"equal (FAI / ZFI = {proteome['fai_over_zfi']:.2f}x) because UniProt IDs "
                "are already compact in text form."
            )

    if any(
        row.get("fai_bytes") and not row.get("zfi_bytes") for row in index_stats
    ):
        lines.append(
            "`.zfi` sidecars were missing on disk when this report was built. "
            "Run `bash bench/index/run.sh` (real-dataset section) or "
            "`./zig-out/bin/z-fasta index` on each REAL_* file before regenerating."
        )

    return " ".join(lines)


def md_zfi_index_size_table(stats: list[dict]) -> str:
    """On-disk `.zfi` vs `.fai` size comparison from benchmark fixtures."""
    if not stats:
        return "_No index files found under `bench/shared/data/`._"

    headers = [
        "Dataset",
        "Records",
        "`.zfi`",
        "`.fai`",
        "FAI / ZFI",
        "ZFI B/rec",
        "FAI B/line",
    ]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in sorted(stats, key=lambda r: dataset_sort_key(r["dataset"])):
        records = row.get("zfi_records") or row.get("fai_records")
        rec_s = f"{records:,}" if records else "n/a"
        zfi_s = format_file_size(row["zfi_bytes"]) if row.get("zfi_bytes") else "n/a"
        fai_s = format_file_size(row["fai_bytes"]) if row.get("fai_bytes") else "n/a"
        ratio = row.get("fai_over_zfi")
        ratio_s = f"{ratio:.2f}x" if ratio else "n/a"
        zfi_bpr = row.get("zfi_bytes_per_rec")
        fai_bpl = row.get("fai_bytes_per_line")
        zfi_bpr_s = f"{zfi_bpr:.0f}" if zfi_bpr else "n/a"
        fai_bpl_s = f"{fai_bpl:.0f}" if fai_bpl else "n/a"
        lines.append(
            f"| {row['dataset']} | {rec_s} | {zfi_s} | {fai_s} | {ratio_s} | "
            f"{zfi_bpr_s} | {fai_bpl_s} |"
        )
    return "\n".join(lines)


def md_zfi_production_section(
    perf_df: pd.DataFrame,
    data_dir: Path,
    nums: ReportCounters,
) -> str:
    """On-disk `.zfi` size, load path, and integrity."""
    _ = perf_df
    index_stats = load_index_format_stats(data_dir)
    t_size = nums.next_table()

    blocks: list[str] = [
        (
            "Cross-tool **Performance** tables use **`z-fasta (.fai)`** so competitors that "
            "write text `.fai` stay format-fair. CLI default `index` writes binary `.zfi`. "
            "Build wall time and RSS for every format × dedup lane are in **z-fasta Format "
            "and Dedup Comparison**; this section is about the on-disk artifact."
        ),
        (
            "### On-disk size\n\n"
            + md_zfi_size_prose(index_stats)
            + f"\n\n**Table {t_size}:** Index file size on the three real datasets "
            "(`bench/shared/data/`). FAI / ZFI above 1.0 means the text `.fai` is larger."
        ),
        md_zfi_index_size_table(index_stats),
        (
            "### Load path\n\n"
            "`loadIndex` opens `.zfi` first, then falls back to `.fai`. For `.zfi`, the "
            "record array is a zero-copy view of mmap'd index bytes, and `lookup_full_map` "
            "builds a pointer hash over the embedded name blob. For `.fai`, the loader mmaps "
            "the text index, parses lines into a record array, and points "
            "`name_offset` / `name_len` into the mmap'd `.fai` (`tryLoadFai` in "
            "`index_format.zig`). `.zfi` entries store `name_offset` / `name_len` into the "
            "mmap'd FASTA as well. `stats --index-only` and small `get` batches can load "
            "either format in `.records_only` mode and skip the name hash; `.fai` stats scans "
            "use `.stats_scan` (on-demand name reads)."
        ),
        (
            "### Messy FASTA\n\n"
            "FAI assumes one `line_bases` / `line_bytes` pair per record. z-fasta indexes "
            "irregular wrap, trailing whitespace, blank lines, and mixed CRLF/LF using "
            "per-record side tables (see **Messy FASTA Compatibility**). Emit `.fai` only "
            "with `index --emit-fai` when every record is FAI-representable; otherwise use "
            "default `.zfi`. Load still falls back to an existing `.fai` for samtools-"
            "compatible text indexes."
        ),
        (
            "### Integrity\n\n"
            "`.zfi` stores `source_size` and, on current writers, a `ZFID` embedded FASTA "
            "mtime. Load rejects size or embedded-mtime mismatch. Legacy `.zfi` without "
            "`ZFID` keeps the weaker index-file mtime check. `.fai` identity remains "
            "mtime-only (no embedded source fields)."
        ),
    ]
    return "\n\n".join(blocks)


def md_memory_rss_headline_table(df: pd.DataFrame) -> str:
    """Peak RSS table for headline tools (same order as Figure 1)."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "peak_rss_mb",
        lambda row: f"{row['peak_rss_mb']:.2f} MB",
        empty_msg="_No memory data._",
    )


def md_page_faults_headline_table(df: pd.DataFrame) -> str:
    """Minor page-fault table for headline tools."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "minor_faults",
        lambda row: f"{int(row['minor_faults']):,}",
        empty_msg="_No page-fault data._",
    )


def md_major_faults_headline_table(df: pd.DataFrame) -> str:
    """Major page-fault table for headline tools."""
    return md_tool_pivot_table(
        df,
        HEADLINE_PERF_TOOLS,
        "dataset",
        "major_faults",
        lambda row: f"{int(row['major_faults']):,}",
        empty_msg="_No page-fault data._",
    )


def md_memory_usage_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Memory usage section: headline peak RSS on real datasets."""
    t_rss = nums.next_table()
    t_cmp = nums.next_table()
    f_mem = nums.next_figure()
    blocks = [
        "Same zebrac runs as **Performance: Real Datasets**. z-fasta (.fai) only; "
        "other product lanes are in **z-fasta Format and Dedup Comparison**.",
        (
            "### Peak RSS\n\n"
            "zebrac starts a new process for each sample and records peak RSS when it "
            "exits. That is the most RAM the process had in use at once (Linux "
            "`ru_maxrss` via `getrusage`). Table and figure show the mean across "
            "samples.\n\n"
            "Use this to compare tools on the same file and host. z-fasta reads through a "
            "bounded buffer, so its input path does not make RSS track FASTA size. Grey "
            "dashed lines show file size as context; output and dedup state can still affect "
            "memory. A low bar does not always mean less work, so page faults and CPU counters "
            "remain separate evidence."
        ),
        f"**Table {t_rss}:** Peak RSS (MB, zebrac mean). Same tool order as Performance.",
        md_memory_rss_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "RSS × = competitor peak RSS / z-fasta peak RSS. Same ratios as bar labels on "
        f"Figure {f_mem}.</summary>",
        "",
        md_zfasta_vs_rss_table(df),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_mem}:** Table {t_rss} as grouped bars (log scale). Bar labels = "
        f"RSS × (see Table {t_cmp}). Gold dashed line = z-fasta RSS; grey dashed = "
        "FASTA file size on disk.",
        f"![Figure {f_mem}: real-dataset peak RSS](results/figures/memory.png)",
        (
            f"**Reading Figure {f_mem}**\n"
            "- Bars: mean peak RSS. `1×` on z-fasta; other labels = competitor / z-fasta.\n"
            "- Gold dashed line: z-fasta RSS for that dataset.\n"
            "- Grey dashed line: file size from `input_bytes`.\n"
            f"- Details in Table {t_cmp}."
        ),
    ]
    return "\n\n".join(blocks)


def md_page_faults_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Page-fault section: headline minor faults on real datasets."""
    major = md_major_faults_headline_table(df)
    t_minor = nums.next_table()
    t_cmp = nums.next_table()
    t_major = nums.next_table()
    f_pf = nums.next_figure()
    blocks = [
        "Same zebrac samples again. Full counter dump is in **Zebrac Counters** below.",
        (
            "### Page faults\n\n"
            "A **minor** fault maps a page without reading disk. A **major** fault reads "
            "from disk. zebrac reports both per run, like wall time. "
            f"Table {t_major} lists major faults separately."
        ),
        f"**Table {t_minor}:** Minor page faults (zebrac mean). Same tool order as Performance.",
        md_page_faults_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "Faults × = competitor minor faults / z-fasta minor faults. Same ratios as bar "
        f"labels on Figure {f_pf}.</summary>",
        "",
        md_zfasta_vs_faults_table(df),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_major}:</strong> Major page faults.</summary>",
        "",
        major,
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_pf}:** Table {t_minor} as grouped bars (log scale). Bar labels = "
        f"Faults × (see Table {t_cmp}). Gold dashed line = z-fasta minor faults.",
        f"![Figure {f_pf}: real-dataset minor page faults](results/figures/page_faults.png)",
        (
            f"**Reading Figure {f_pf}**\n"
            "- Bars: mean minor faults. `1×` on z-fasta; other labels = competitor / z-fasta.\n"
            "- Gold dashed line: z-fasta count for that dataset.\n"
            "- Lower bar = fewer minor faults.\n"
            f"- Details in Table {t_cmp}."
        ),
    ]
    return "\n\n".join(blocks)


def _fmt_ipc(v: float | None) -> str:
    if v is None or pd.isna(v):
        return "n/a"
    return f"{float(v):.3f}"


def _fmt_rate_pct(v: float | None) -> str:
    if v is None or pd.isna(v):
        return "n/a"
    return f"{float(v) * 100:.3f}%"


def _fmt_instr_per_mib(v: float | None) -> str:
    if v is None or pd.isna(v):
        return "n/a"
    v = float(v)
    if v >= 1e9:
        return f"{v / 1e9:.2f}B"
    if v >= 1e6:
        return f"{v / 1e6:.2f}M"
    if v >= 1e3:
        return f"{v / 1e3:.1f}K"
    return f"{v:.0f}"


def md_zebrac_instr_per_mib_headline_table(df: pd.DataFrame) -> str:
    """Instructions per MiB pivot for headline tools."""
    work = enrich_counter_metrics(filter_headline_perf(df))
    if work.empty or "instructions_per_mib" not in work.columns:
        return "_No counter data._"
    work["cell"] = work["instructions_per_mib"].apply(_fmt_instr_per_mib)
    pivot = work.pivot(index="dataset", columns="tool", values="cell")
    return _headline_pivot(pivot).to_markdown()


def md_zebrac_secondary_derived_table(df: pd.DataFrame) -> str:
    """IPC, cache miss rate, and branch miss rate for headline tools."""
    work = enrich_counter_metrics(filter_headline_perf(df))
    if work.empty or "ipc" not in work.columns:
        return "_No counter data._"

    rows: list[dict] = []
    for ds in sorted(work["dataset"].unique(), key=dataset_sort_key):
        for tool in HEADLINE_PERF_TOOLS:
            hit = work[(work["dataset"] == ds) & (work["tool"] == tool)]
            if hit.empty:
                continue
            r = hit.iloc[0]
            rows.append(
                {
                    "Dataset": ds,
                    "Tool": display_tool(tool),
                    "IPC": _fmt_ipc(r["ipc"]),
                    "Cache miss %": _fmt_rate_pct(r["cache_miss_rate"]),
                    "Branch miss / instr %": _fmt_rate_pct(r["branch_miss_rate"]),
                }
            )
    if not rows:
        return "_No counter data._"
    return pd.DataFrame(rows).to_markdown(index=False)


def md_zfasta_vs_instr_per_mib_table(df: pd.DataFrame) -> str:
    work = enrich_counter_metrics(filter_headline_perf(df))
    comparisons = build_headline_ratio_comparisons(work, "instructions_per_mib")
    return md_zfasta_vs_ratio_table(
        comparisons,
        zf_label="z-fasta",
        comp_label="Competitor",
        ratio_label="Work ×",
        fmt_zf=lambda r: _fmt_instr_per_mib(r.zfasta_v),
        fmt_comp=lambda r: _fmt_instr_per_mib(r.comp_v),
    )


def md_zebrac_counter_table(df: pd.DataFrame) -> str:
    cols = [
        "dataset",
        "tool",
        "peak_rss_mb",
        "minor_faults",
        "major_faults",
        "cpu_cycles",
        "instructions",
        "cache_references",
        "cache_misses",
        "branch_misses",
    ]
    have_cols = [c for c in cols if c in df.columns]
    if "peak_rss_mb" not in have_cols:
        return ""
    display = df[have_cols].copy()
    display["_dataset_order"] = display["dataset"].map(dataset_sort_key)
    display["_tool_order"] = display["tool"].map(tool_sort_key)
    display = display.sort_values(["_dataset_order", "_tool_order"]).drop(
        columns=["_dataset_order", "_tool_order"]
    )
    display = display.rename(columns={
        "dataset": "Dataset",
        "tool": "Tool",
        "peak_rss_mb": "Peak RSS (MB)",
        "minor_faults": "Minor Faults",
        "major_faults": "Major Faults",
        "cpu_cycles": "CPU Cycles",
        "instructions": "Instructions",
        "cache_references": "Cache References",
        "cache_misses": "Cache Misses",
        "branch_misses": "Branch Misses",
    })
    if "Tool" in display.columns:
        display["Tool"] = display["Tool"].map(display_tool)
    int_cols = (
        "Minor Faults",
        "Major Faults",
        "CPU Cycles",
        "Instructions",
        "Cache References",
        "Cache Misses",
        "Branch Misses",
    )
    for col in int_cols:
        if col in display.columns:
            display[col] = display[col].apply(lambda v: f"{int(v):,}" if pd.notna(v) else "")
    return display.to_markdown(index=False, floatfmt=".2f")


def md_zebrac_counters_real_datasets(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Hardware counters: CPU work per MiB figure, secondary metrics in tables."""
    t_work = nums.next_table()
    t_cmp = nums.next_table()
    t_other = nums.next_table()
    t_raw = nums.next_table()
    f_work = nums.next_figure()
    blocks = [
        "Same zebrac samples as **Performance: Real Datasets**. z-fasta (.fai) only; "
        "other product lanes are in **z-fasta Format and Dedup Comparison**.",
        (
            "### CPU work per MiB\n\n"
            "**Instructions per MiB** is retired CPU instructions divided by input size "
            "in MiB (zebrac counters and metadata `input_bytes`). **Lower is better**: "
            "less work to index the same file.\n\n"
            "IPC and cache metrics are in Table {other} (collapsed). They help profiling "
            "but measure per-cycle or cache behavior, not total work per byte."
        ).format(other=t_other),
        (
            f"**Table {t_work}:** Instructions per MiB (zebrac mean). Same tool order as "
            "Performance. Lower = less CPU work."
        ),
        md_zebrac_instr_per_mib_headline_table(df),
        "<details>",
        f"<summary><strong>Table {t_cmp}:</strong> z-fasta vs each competitor. "
        "Work × = competitor instructions/MiB ÷ z-fasta instructions/MiB. Same ratios "
        f"as bar labels on Figure {f_work}.</summary>",
        "",
        md_zfasta_vs_instr_per_mib_table(df),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_other}:</strong> Other derived metrics: IPC "
        "(instructions / CPU cycles), cache miss rate (cache misses / cache references), "
        "branch miss rate (branch misses / instructions). Profiling detail only.</summary>",
        "",
        md_zebrac_secondary_derived_table(df),
        "</details>",
        "<details>",
        f"<summary><strong>Table {t_raw}:</strong> Full zebrac counter dump (all tools, "
        "including z-fasta product lanes).</summary>",
        "",
        md_zebrac_counter_table(df),
        "</details>",
        '<div style="margin: 1.5em 0"></div>',
        f"**Figure {f_work}:** Table {t_work} as grouped bars (log scale). Bar labels = "
        f"Work × (see Table {t_cmp}). Gold dashed line = z-fasta instructions/MiB.",
        f"![Figure {f_work}: real-dataset CPU work per MiB](results/figures/zebrac_instructions_per_mib.png)",
        (
            f"**Reading Figure {f_work}**\n"
            "- Bars: mean instructions per MiB (log scale). Lower bar = less CPU work.\n"
            "- `1×` on z-fasta; other labels = competitor work / z-fasta work "
            f"(Work ×; see Table {t_cmp}).\n"
            "- Gold dashed line: z-fasta value for that dataset.\n"
            f"- Details in Table {t_cmp}."
        ),
    ]
    return "\n\n".join(blocks)


def load_messy_index(
    results_dir: Path, manifest: dict | None = None
) -> pd.DataFrame | None:
    """Load messy FASTA indexing results pinned to the run timestamp."""
    if not manifest or not manifest.get("timestamp"):
        return None
    d = results_dir / f"messy_{manifest['timestamp']}"
    if not d.is_dir():
        return None
    meta_by_command: dict[str, str] = {}
    meta_path = d / "metadata.jsonl"
    if meta_path.exists():
        for line in meta_path.read_text().splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            meta_by_command[row["command"]] = row["tool"]

    rows = []
    for f in sorted(d.glob("*.json")):
        data = json.loads(f.read_text())
        variant = f.stem
        for r in data.get("results", []):
            cmd = r.get("command", "")
            tool = meta_by_command.get(cmd)
            if not tool:
                tool = cmd.split()[0].split("/")[-1]
                if tool == "faidx":
                    tool = "pyfaidx"
            if tool == "rustbio-custom-index":
                tool = "rust-bio"
            if tool == "z-fasta":
                tool = (
                    "z-fasta (.fai)"
                    if "--emit-fai" in cmd
                    else "z-fasta (.zfi)"
                )
            failed = r.get("failed_sample_count", 0)
            samples = r.get("sample_count", 0)
            success = failed == 0 and samples > 0
            rows.append({"variant": variant, "tool": tool, "success": success})
    if not rows:
        return None
    return pd.DataFrame(rows)


def md_messy_table(df: pd.DataFrame) -> str:
    """Markdown compatibility matrix for messy FASTA variants."""
    messy_order = [
        "z-fasta (.zfi)",
        "z-fasta (.fai)",
        "samtools",
        "noodles",
        "rust-bio",
        "fastahack",
        "pyfaidx",
    ]
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
            elif tool.startswith("z-fasta (.zfi)"):
                row[tool] = "ok" if sub.iloc[0]["success"] else "fail"
            else:
                row[tool] = "not compatible"
        rows.append(row)
    return pd.DataFrame(rows).to_markdown(index=False)


def md_edge_case_summary_table(df: pd.DataFrame) -> str:
    """Summarize comparability by the rule used for each case class."""
    classes = edge_contract_classes(df)
    meanings = {
        "fai_parity": "z-fasta (.fai) vs samtools acceptance and FAI bytes",
        "zfi_messy": "z-fasta (.zfi) side table; FAI peers expected incompatible",
        "invalid_input": "raw behavior shown; no parity claim for binary input",
    }
    rows = []
    for class_name, label in EDGE_CONTRACT_CLASSES.items():
        subset = df[classes == class_name]
        if subset.empty:
            continue
        matches = int((subset["output_match"] == "MATCH").sum())
        rows.append(
            {
                "Contract basis": label,
                "Cases": len(subset),
                "Matches": f"{matches}/{len(subset)}",
                "Comparison": meanings[class_name],
            }
        )
    return pd.DataFrame(rows).to_markdown(index=False)


def md_edge_case_section(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Edge-case correctness: summary table, heatmap figure, and reading notes."""
    t_summary = nums.next_table()
    f_edge = nums.next_figure()
    blocks = [
        (
            "Structural edge-case fixtures from `bench/index/run.sh` (edge_cases/). Each "
            "row in the heatmap is one test file. The z-fasta columns separate the "
            "samtools-compatible FAI attempt from the production `.zfi` index. The "
            "**Contract** column marks whether z-fasta met that test's rules."
        ),
        (
            "**Contract meaning:** for structural cases, `z-fasta (.fai)` must agree with "
            "samtools on FAI acceptance and bytes when both accept. If one accepts and the "
            "other rejects, the contract is `N`. For the four named messy fixtures, the "
            "`.zfi` side-table path is the contract because samtools cannot represent those "
            "layouts; the FAI rejection pair plus valid `.zfi` side table is `Y`."
        ),
        "**Scoring rules:**",
        (
            "- *Samtools parity* cases: when both z-fasta and samtools accept a file, their "
            "`.fai` output must match. When both reject it, that also counts as a match. "
            "A production `.zfi` success does not hide an FAI exit-code mismatch."
        ),
        (
            "- *z-fasta-only* messy cases: z-fasta must index files where samtools rejects "
            "non-uniform wrapping (side-table `.zfi`)."
        ),
        (
            f"**Table {t_summary}:** Comparability summary by contract class. Raw exit codes "
            "are diagnostic; `Matches` applies the class-specific comparison rule."
        ),
        md_edge_case_summary_table(df),
        '<div style="margin: 1.5em 0"></div>',
        (
            f"**Figure {f_edge}:** Raw per-test exit codes and the composite contract result. "
            "Green = exit 0; red = non-zero exit; gray = tool not run; orange = contract "
            "mismatch."
        ),
        (
            f"![Figure {f_edge}: edge-case exit codes](results/figures/edge_cases.png)"
        ),
        (
            f"**Reading Figure {f_edge}**\n"
            "- Rows: test case names from `bench/index/run.sh` edge_cases/.\n"
            "- Tool columns show the numeric exit code; `n/a` means the optional tool was not run.\n"
            "- **Contract:** `Y` = MATCH, `N` = mismatch vs test rules.\n"
            "- Green / red / gray / orange cells match the legend on the chart.\n"
            f"- Aggregate counts are in Table {t_summary}."
        ),
    ]
    return "\n\n".join(blocks)


def md_messy_section(df: pd.DataFrame, nums: ReportCounters) -> str:
    """Messy FASTA compatibility: matrix table and reading notes."""
    t_matrix = nums.next_table()
    if "z-fasta (.zfi)" in set(df["tool"]):
        zf_lane = "`z-fasta (.zfi)` via default `index`"
    elif "z-fasta (.fai)" in set(df["tool"]):
        zf_lane = "historical `z-fasta (.fai)` data"
    else:
        zf_lane = "no z-fasta lane"
    blocks = [
        (
            "Proteome-derived fixtures in `bench/shared/cache/messy_perf/`. Each cell is zebrac "
            "with `--allow-failures` (repeated samples, not a single exit check)."
        ),
        (
            f"**z-fasta lane:** This table must use {zf_lane}. The `.fai` lane is not the "
            "messy compatibility claim; variable widths, trailing whitespace, and mixed "
            "line endings require the production `.zfi` side table."
        ),
        (
            "**What z-fasta handles:** irregular line wrapping, trailing whitespace on "
            "sequence lines, blank lines between records, and mixed CRLF/LF. Those are "
            "common export glitches, not arbitrary byte corruption. z-fasta records true "
            "sequence boundaries in a side-table `.zfi` instead of assuming every line in a "
            "record has the same width."
        ),
        (
            "**What we still reject:** empty files, headers without sequence, malformed "
            "headers, and other cases where indexing would be meaningless (see **Edge Case "
            "Correctness** above). We do not claim to repair truncated or binary-corrupted FASTA."
        ),
        (
            "`mixed_crlf` deliberately alternates LF and CRLF sequence-line endings. It is "
            "not a uniformly CRLF-wrapped control file."
        ),
        (
            f"**Table {t_matrix}:** Indexing success per variant and tool from the messy "
            "benchmark CSV. The z-fasta `.zfi` lane is `ok` only when every sample exits 0. "
            "FAI peers are `not compatible` for these variable-width, whitespace, and mixed-"
            "ending layouts; an exit 0 alone does not prove that their index reproduces the "
            "sequence correctly."
        ),
        md_messy_table(df),
        (
            f"**Reading Table {t_matrix}**\n"
            "- Rows: messy variant names (`all_messy`, `mixed_crlf`, `mixed_widths`, "
            "`trailing_whitespace`).\n"
            "- Columns: production `z-fasta (.zfi)` followed by FAI peers.\n"
            "- `not compatible` is a layout contract result, not a raw process exit code.\n"
            "- Refresh with `bash bench/index/run.sh` (messy zebrac section)."
        ),
    ]
    return "\n\n".join(blocks)


def md_tools_tested(tools: dict[str, str] | None = None, z_fasta: str | None = None) -> str:
    """Bullet list of tools in the index suite."""
    tools = tools or {}

    def version_note(key: str) -> str:
        if tools.get(key):
            return f" ({tools[key]})"
        pin = VERSION_PINS.get(key)
        if pin:
            return f" (pin {pin})"
        return ""

    zf = f" ({z_fasta})" if z_fasta else ""
    return (
        f"- **z-fasta (.fai){zf}:** `index --emit-fai` (bounded reader, dedup). Peer-comparable "
        "lane for cross-tool tables and scaling.\n"
        "- **z-fasta (.fai) --no-dedup:** `index --emit-fai --no-dedup`.\n"
        "- **z-fasta (.zfi):** `index` (bounded reader, dedup). CLI default; writes `.zfi` "
        "on disk.\n"
        "- **z-fasta (.zfi) --no-dedup:** `index --no-dedup`.\n"
        f"- **samtools{version_note('samtools')}:** `samtools faidx`, industry reference.\n"
        f"- **seqkit{version_note('seqkit')}:** `seqkit faidx`, Go toolkit.\n"
        f"- **fastahack{version_note('fastahack')}:** `fastahack -i`, indexes on first access.\n"
        f"- **pyfaidx{version_note('pyfaidx')}:** `faidx --no-output`, Python wrapper.\n"
        f"- **noodles{version_note('noodles')}:** `tools/noodles_wrapper index`, "
        "Tier 2 noodles-fasta wrapper.\n"
        f"- **rust-bio{version_note('rustbio')}:** `tools/rustbio_wrapper index`. "
        "Custom CLI because rust-bio has no standalone index command. Strict FAI rules "
        "aligned with noodles and samtools.\n"
    )


# ══════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════


def main():
    parser = argparse.ArgumentParser(
        description="Generate index benchmark report and figures from raw results."
    )
    parser.add_argument(
        "results_dir",
        nargs="?",
        type=Path,
        help="Results directory, defaults to bench/index/results.",
    )
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Allow smoke or partial runs to overwrite the tracked report.",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    results_dir = args.results_dir if args.results_dir else script_dir / "results"
    project_root = results_dir.parent.parent.parent
    manifest = enrich_manifest(load_run_manifest(results_dir), results_dir)
    publication_failures = publication_issues(manifest, results_dir)
    incomplete = bool(publication_failures) or is_incomplete_run(
        manifest,
        required_sections=(
            "real",
            "scale_size",
            "scale_seqs_budget",
            "scale_seqs_fixed",
        ),
        skip_flags=("skip_real", "skip_scaling", "skip_size"),
    )
    if incomplete and not args.allow_incomplete:
        timestamp = manifest.get("timestamp", "unknown") if manifest else "unknown"
        runs = manifest.get("runs", "unknown") if manifest else "unknown"
        skip_real = manifest.get("skip_real", "unknown") if manifest else "unknown"
        skip_scaling = manifest.get("skip_scaling", "unknown") if manifest else "unknown"
        skip_size = manifest.get("skip_size", "unknown") if manifest else "unknown"
        shown_issues = publication_failures[:8]
        remaining_issue_count = len(publication_failures) - len(shown_issues)
        issue_text = "; ".join(shown_issues) or "section or sample settings"
        if remaining_issue_count:
            issue_text += f"; {remaining_issue_count} more"
        raise SystemExit(
            "Refusing to overwrite tracked benchmark report from incomplete run "
            f"{timestamp} (runs={runs}, skip_real={skip_real}, "
            f"skip_scaling={skip_scaling}, skip_size={skip_size}; "
            f"issues={issue_text}). "
            "Re-run the full suite "
            "or pass --allow-incomplete for local smoke analysis."
        )
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    report_lines: list[str] = []

    def section(title, level=2):
        report_lines.append(f"{'#' * level} {title}")

    report_lines.append("# z-fasta Index Benchmark Report")
    report_lines.append("_Auto-generated by `generate_report.py` from zebrac results._")
    if incomplete:
        report_lines.append(
            "> **DRAFT: incomplete run.** Required sections, source/build identity, exact "
            "cells, sample validity, or correctness evidence are incomplete. Do not publish "
            "until a full run completes."
        )

    section("Overview")
    report_lines.append(md_overview(manifest))

    section("Run Provenance")
    report_lines.append(md_run_provenance(manifest, project_root))

    generated_figs: list[str] = []
    metadata_df = load_metadata(results_dir, manifest)

    perf_df = load_perf(results_dir, metadata_df, manifest)
    nums = ReportCounters()
    if perf_df is not None and len(perf_df):
        section("Performance: Real Datasets")
        report_lines.append(md_performance_real_datasets(perf_df, nums))

        p = fig_performance_throughput(perf_df, figures_dir / "performance.png")
        generated_figs.append(str(p))

        section("Memory Usage: Real Datasets")
        report_lines.append(md_memory_usage_real_datasets(perf_df, nums))
        p = fig_memory_rss(perf_df, figures_dir / "memory.png")
        generated_figs.append(str(p))

        section("Page Faults: Real Datasets")
        report_lines.append(md_page_faults_real_datasets(perf_df, nums))
        p = fig_page_faults(perf_df, figures_dir / "page_faults.png")
        generated_figs.append(str(p))

        section("Zebrac Counters: Real Datasets")
        report_lines.append(md_zebrac_counters_real_datasets(perf_df, nums))
        p = fig_zebrac_instructions_per_mib(
            perf_df, figures_dir / "zebrac_instructions_per_mib.png"
        )
        if p:
            generated_figs.append(str(p))
    else:
        section("Performance: Real Datasets")
        report_lines.append("_No real-dataset zebrac data found._")

    size_df = load_scaling(results_dir, "scale_size", "size_mb", metadata_df, manifest)
    if size_df is not None and len(size_df):
        section("Scaling: File Size")
        report_lines.append(
            md_scaling_section(
                size_df,
                nums,
                intro=(
                    "Synthetic FASTA files from 1 MB to 1 GB (one sequence per file). "
                    "Same zebrac warm-cache setup as **Performance: Real Datasets**. "
                    "z-fasta (.fai) only; other product lanes are in **z-fasta Format and "
                    "Dedup Comparison**."
                ),
                subsection=(
                    "### Wall time vs file size\n\n"
                    "Mean wall time (seconds, zebrac) at each file size. **Lower is better.** "
                    "Sequence count stays near constant while total bytes on disk grow."
                ),
                param_col="size_mb",
                param_label="File Size",
                figure_path="results/figures/scaling_size.png",
                figure_title="scaling by file size",
                table_x_label="file size",
                reading_x_label="file size (MB)",
            )
        )
        p = fig_scaling_headline_lines(
            size_df,
            figures_dir / "scaling_size.png",
            param_col="size_mb",
            xlabel="File Size (MB, log scale)",
            title="Scaling: Wall Time vs File Size",
            fig_note="One line per tool. Lower wall time is better.",
        )
        generated_figs.append(str(p))

    seq_budget_df = load_scaling(
        results_dir, "scale_seqs_budget", "seq_count", metadata_df, manifest
    )
    if seq_budget_df is not None and len(seq_budget_df):
        section("Scaling: Sequence Count (Bounded Bytes)")
        report_lines.append(
            md_scaling_section(
                seq_budget_df,
                nums,
                intro=(
                    "Synthetic FASTA holding about **50 MiB** total per file; sequence length "
                    "shrinks as record count rises. Same zebrac warm-cache setup as "
                    "**Performance: Real Datasets**. z-fasta (.fai) only; other product lanes "
                    "are in **z-fasta Format and Dedup Comparison**."
                ),
                subsection=(
                    "### Wall time vs sequence count (bounded bytes)\n\n"
                    "Mean wall time (seconds, zebrac) at each record count. **Lower is better.** "
                    "Isolates per-record and header overhead at near-constant file size."
                ),
                param_col="seq_count",
                param_label="Sequences",
                figure_path="results/figures/scaling_seqs_budget.png",
                figure_title="bounded-bytes sequence scaling",
                table_x_label="record count",
                reading_x_label="sequence count",
            )
        )
        p = fig_scaling_headline_lines(
            seq_budget_df,
            figures_dir / "scaling_seqs_budget.png",
            param_col="seq_count",
            xlabel="Sequences (log scale)",
            title="Scaling: Wall Time vs Sequence Count (Bounded ~50 MiB)",
            fig_note="One line per tool. Lower wall time is better.",
        )
        generated_figs.append(str(p))

    seq_fixed_df = load_scaling(
        results_dir, "scale_seqs_fixed", "seq_count", metadata_df, manifest
    )
    if seq_fixed_df is not None and len(seq_fixed_df):
        section("Scaling: Sequence Count (Fixed Length)")
        report_lines.append(
            md_scaling_section(
                seq_fixed_df,
                nums,
                intro=(
                    "Synthetic FASTA with **1024 bp** per sequence; total file size grows with "
                    "record count (100k to 1M). Same zebrac warm-cache setup as "
                    "**Performance: Real Datasets**. z-fasta (.fai) only; other product lanes "
                    "are in **z-fasta Format and Dedup Comparison**."
                ),
                subsection=(
                    "### Wall time vs sequence count (fixed length)\n\n"
                    "Mean wall time (seconds, zebrac) at each record count. **Lower is better.** "
                    "Bytes and index entries grow together; the 1M point is about 1 GiB on disk."
                ),
                param_col="seq_count",
                param_label="Sequences",
                figure_path="results/figures/scaling_seqs_fixed.png",
                figure_title="fixed-length sequence scaling",
                table_x_label="record count",
                reading_x_label="sequence count",
            )
        )
        p = fig_scaling_headline_lines(
            seq_fixed_df,
            figures_dir / "scaling_seqs_fixed.png",
            param_col="seq_count",
            xlabel="Sequences (log scale)",
            title="Scaling: Wall Time vs Sequence Count (1024 bp / record)",
            fig_note="One line per tool. Lower wall time is better.",
        )
        generated_figs.append(str(p))

    test_df = load_tests(results_dir, manifest)
    if perf_df is not None and len(perf_df) and not filter_product_comparison(perf_df).empty:
        section("z-fasta Format and Dedup Comparison")
        report_lines.append(md_product_comparison_section(perf_df, nums))
        stats_by_id = {
            entry["id"]: entry
            for entry in load_real_dataset_stats(project_root / "bench" / "shared" / "data")
        }
        p = fig_product_comparison_throughput(
            perf_df, figures_dir / "product_comparison_performance.png"
        )
        if p:
            generated_figs.append(str(p))
        p = fig_product_comparison_memory(
            perf_df, figures_dir / "product_comparison_memory.png"
        )
        if p:
            generated_figs.append(str(p))
        p = fig_product_comparison_tradeoff(
            perf_df,
            figures_dir / "product_comparison_tradeoff.png",
            stats_by_id=stats_by_id,
        )
        if p:
            generated_figs.append(str(p))

    if perf_df is not None and len(perf_df):
        index_stats = load_index_format_stats(project_root / "bench" / "shared" / "data")
        if index_stats:
            section("z-fasta Production Index (.zfi)")
            report_lines.append(
                md_zfi_production_section(
                    perf_df,
                    project_root / "bench" / "shared" / "data",
                    nums,
                )
            )

    if test_df is not None and len(test_df):
        section("Edge Case Correctness")
        report_lines.append(md_edge_case_section(test_df, nums))
        p = fig_edge_heatmap(test_df, figures_dir / "edge_cases.png")
        generated_figs.append(str(p))

    messy_df = load_messy_index(results_dir, manifest)
    if messy_df is not None and len(messy_df):
        section("Messy FASTA Compatibility")
        report_lines.append(md_messy_section(messy_df, nums))

    section("Tools Tested")
    report_lines.append(
        md_tools_tested(
            (manifest or {}).get("tools"),
            (manifest or {}).get("z_fasta"),
        )
    )

    # ── Write report ───────────────────────────────────────────────
    if not generated_figs and not args.allow_incomplete:
        raise SystemExit(
            "Refusing to overwrite tracked benchmark report because no benchmark "
            "result data was found. Run the benchmark suite first or pass "
            "--allow-incomplete for local report drafts."
        )

    # Write REPORT.md one level up from results/ (into bench/) so it
    # isn't gitignored with the rest of the generated results data.
    bench_dir = results_dir.parent
    report_path = bench_dir / "REPORT.md"
    report_path.write_text(f"{MARKDOWNLINT_DISABLE}\n\n{normalize_markdown(report_lines)}")

    print(f"Report: {report_path}")
    print(f"Figures: {len(generated_figs)}")
    for f in generated_figs:
        print(f"  {f}")


if __name__ == "__main__":
    main()
