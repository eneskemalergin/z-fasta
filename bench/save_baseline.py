#!/usr/bin/env python3
"""
Capture a machine-readable performance baseline from bench results.

Writes a timestamped directory under bench/baselines/ (gitignored) containing:
  manifest.json   git rev, host, CPU, optimize mode, suite timestamps
  index.json      condensed index metrics (real datasets + key scaling points)
  get.json        condensed get metrics
  stats.json      condensed stats metrics
  compare.json    flat key->seconds map for diffing against a future baseline

Usage:
    .venv/bin/python bench/save_baseline.py
    .venv/bin/python bench/save_baseline.py --label v0.2.9-pre-fix
    .venv/bin/python bench/compare_baseline.py bench/baselines/<old> bench/baselines/<new>
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

BENCH_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = BENCH_ROOT.parent
BASELINES_DIR = BENCH_ROOT / "baselines"

INDEX_BUNDLE_PREFIXES = ("perf", "scale_size", "scale_seqs", "memory")
GET_BUNDLE_PREFIXES = (
    "single",
    "fullseq",
    "scale_region",
    "real",
    "bed",
    "multi",
    "memory",
    "rc_review",
)
STATS_BUNDLE_PREFIXES = (
    "real",
    "indexonly",
    "stats",
    "memory",
    "scale_size",
    "throughput",
)


def _extract_timestamp(name: str, prefix: str) -> str | None:
    match = re.match(rf"^{re.escape(prefix)}_(\d{{8}}_\d{{6}})(?:\.csv)?$", name)
    return match.group(1) if match else None


def infer_bundle_timestamp(results_dir: Path, prefixes: tuple[str, ...]) -> str | None:
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


def discover_dataset(
    results_dir: Path,
    prefix: str,
    bundle_timestamp: str | None = None,
) -> Path | None:
    if bundle_timestamp is not None:
        exact_dir = results_dir / f"{prefix}_{bundle_timestamp}"
        if exact_dir.is_dir():
            return exact_dir
        exact_csv = results_dir / f"{prefix}_{bundle_timestamp}.csv"
        if exact_csv.is_file():
            return exact_csv
        return None
    dirs = sorted(
        (d for d in results_dir.glob(f"{prefix}_*") if d.is_dir()),
        reverse=True,
    )
    if dirs:
        return dirs[0]
    csvs = sorted(results_dir.glob(f"{prefix}_*.csv"), reverse=True)
    return csvs[0] if csvs else None


def load_hyperfine(path: Path) -> list[dict]:
    with open(path) as f:
        data = json.load(f)
    rows = []
    for r in data.get("results", []):
        rows.append(
            {
                "tool": r["command"],
                "mean_s": r["mean"],
                "stddev_s": r.get("stddev"),
                "median_s": r.get("median"),
                "user_s": r.get("user"),
                "system_s": r.get("system"),
            }
        )
    return rows


def load_perf_dir(perf_dir: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for jf in sorted(perf_dir.glob("*.json")):
        out[jf.stem] = load_hyperfine(jf)
    return out


def load_memory_csv(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "tool": row["tool"],
                    "time_s": float(row["time_s"]),
                    "mem_kb": int(row["mem_kb"]),
                    "major_faults": int(row["major_faults"]),
                    "minor_faults": int(row["minor_faults"]),
                }
            )
    return rows


def git_info() -> dict:
    def run(*args: str) -> str:
        try:
            return (
                subprocess.check_output(["git", *args], cwd=PROJECT_ROOT, text=True)
                .strip()
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return ""

    rev = run("rev-parse", "HEAD")
    short = run("rev-parse", "--short", "HEAD")
    branch = run("rev-parse", "--abbrev-ref", "HEAD")
    describe = run("describe", "--tags", "--always", "--dirty")
    dirty = bool(run("status", "--porcelain"))
    return {
        "commit": rev,
        "commit_short": short,
        "branch": branch,
        "describe": describe,
        "dirty": dirty,
    }


def machine_info() -> dict:
    cpu = platform.processor() or "unknown"
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass

    mem_gb = None
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    mem_gb = round(kb / 1024 / 1024, 1)
                    break
    except OSError:
        pass

    return {
        "hostname": platform.node(),
        "os": platform.platform(),
        "kernel": platform.release(),
        "cpu": cpu,
        "logical_cpus": os_cpu_count(),
        "ram_gib": mem_gb,
    }


def os_cpu_count() -> int | None:
    try:
        with open("/proc/cpuinfo") as f:
            return sum(1 for line in f if line.startswith("processor"))
    except OSError:
        return None


def zfasta_binary_info() -> dict:
    binary = PROJECT_ROOT / "zig-out" / "bin" / "z-fasta"
    info: dict = {"path": str(binary)}
    if binary.is_file():
        info["size_bytes"] = binary.stat().st_size
        try:
            out = subprocess.check_output([str(binary), "--version"], text=True).strip()
            info["version"] = out
        except (subprocess.CalledProcessError, OSError):
            pass
    return info


def pick_zfasta(rows: list[dict], mode: str = "z-fasta-default") -> dict | None:
    for r in rows:
        if r["tool"] == mode or r["tool"].endswith(mode):
            return r
    for r in rows:
        if "z-fasta" in r["tool"]:
            return r
    return None


def flatten_compare(suite_name: str, suite_data: dict) -> dict[str, float]:
    flat: dict[str, float] = {}
    for section in (
        "real_datasets",
        "scale_size_mb",
        "scale_seqs",
        "single_region",
        "full_sequence",
        "multi_region",
        "bed_batch",
        "test_files",
        "indexonly_mb",
    ):
        bucket = suite_data.get(section)
        if not isinstance(bucket, dict):
            continue
        for param, rows in bucket.items():
            if not isinstance(rows, list):
                continue
            for row in rows:
                key = f"{suite_name}.{section}.{param}.{row['tool']}"
                flat[key] = row["mean_s"]
    return flat


def anchor_bundle_timestamp(results_dir: Path, anchor_prefix: str) -> str | None:
    anchor = discover_dataset(results_dir, anchor_prefix, None)
    if anchor is None:
        return None
    return _extract_timestamp(anchor.name, anchor_prefix)


def bundle_warnings(
    suite_name: str,
    results_dir: Path,
    prefixes: tuple[str, ...],
    bundle_timestamp: str | None,
    timestamps: dict[str, str],
) -> list[str]:
    warnings: list[str] = []
    if bundle_timestamp is None:
        return warnings
    for prefix in prefixes:
        if prefix in timestamps:
            continue
        if discover_dataset(results_dir, prefix, bundle_timestamp) is None:
            warnings.append(
                f"{suite_name}: no {prefix} results for bundle {bundle_timestamp}"
            )
    return warnings


def collect_index() -> dict:
    results = BENCH_ROOT / "index" / "results"
    bundle_ts = anchor_bundle_timestamp(results, "perf") or infer_bundle_timestamp(
        results, INDEX_BUNDLE_PREFIXES
    )

    data: dict = {"timestamps": {}, "bundle_timestamp": bundle_ts}
    perf = discover_dataset(results, "perf", bundle_ts)
    scale_size = discover_dataset(results, "scale_size", bundle_ts)
    scale_seqs = discover_dataset(results, "scale_seqs", bundle_ts)
    memory = discover_dataset(results, "memory", bundle_ts)

    if perf:
        data["timestamps"]["perf"] = perf.name
        data["real_datasets"] = load_perf_dir(perf)
        zf = pick_zfasta(data["real_datasets"].get("Genome", []))
        if zf:
            data["headline_genome_index_s"] = zf["mean_s"]
    if scale_size and scale_size.is_dir():
        data["timestamps"]["scale_size"] = scale_size.name
        all_scale = load_perf_dir(scale_size)
        keep = ("1mb", "10mb", "100mb", "1000mb")
        data["scale_size_mb"] = {k: v for k, v in all_scale.items() if k in keep}
    if scale_seqs and scale_seqs.is_dir():
        data["timestamps"]["scale_seqs"] = scale_seqs.name
        all_seqs = load_perf_dir(scale_seqs)
        keep = ("10", "1000", "100000")
        data["scale_seqs"] = {k: v for k, v in all_seqs.items() if k in keep}
    if memory and memory.suffix == ".csv":
        data["timestamps"]["memory"] = memory.name
        data["memory"] = load_memory_csv(memory)

    data["bundle_warnings"] = bundle_warnings(
        "index", results, INDEX_BUNDLE_PREFIXES, bundle_ts, data["timestamps"]
    )
    return data


def load_json_dir(
    results: Path,
    prefix: str,
    inner_key: str,
    data: dict,
    bundle_ts: str | None,
) -> None:
    d = discover_dataset(results, prefix, bundle_ts)
    if not d or not d.is_dir():
        return
    data["timestamps"][prefix] = d.name
    section: dict[str, list[dict]] = {}
    for jf in sorted(d.glob("*.json")):
        section[jf.stem] = load_hyperfine(jf)
    data[inner_key] = section


def collect_get() -> dict:
    results = BENCH_ROOT / "get" / "results"
    data: dict = {"timestamps": {}}

    bundle_ts = anchor_bundle_timestamp(results, "real") or infer_bundle_timestamp(
        results, GET_BUNDLE_PREFIXES
    )
    manifest_files = sorted(results.glob("run_*.json"), reverse=True)
    if manifest_files:
        with open(manifest_files[0]) as f:
            manifest = json.load(f)
        data["run_manifest"] = manifest
        data["timestamps"]["run"] = manifest_files[0].name
        manifest_ts = manifest.get("timestamp")
        if manifest_ts:
            bundle_ts = manifest_ts

    data["bundle_timestamp"] = bundle_ts

    load_json_dir(results, "single", "single_region", data, bundle_ts)
    load_json_dir(results, "real", "real_datasets", data, bundle_ts)
    load_json_dir(results, "multi", "multi_region", data, bundle_ts)
    load_json_dir(results, "fullseq", "full_sequence", data, bundle_ts)
    load_json_dir(results, "bed", "bed_batch", data, bundle_ts)

    mem = discover_dataset(results, "memory", bundle_ts)
    if mem and mem.suffix == ".csv":
        data["timestamps"]["memory"] = mem.name
        data["memory"] = load_memory_csv(mem)

    data["bundle_warnings"] = bundle_warnings(
        "get", results, GET_BUNDLE_PREFIXES, bundle_ts, data["timestamps"]
    )
    return data


def collect_stats() -> dict:
    results = BENCH_ROOT / "stats" / "results"
    bundle_ts = anchor_bundle_timestamp(results, "real") or infer_bundle_timestamp(
        results, STATS_BUNDLE_PREFIXES
    )

    data: dict = {"timestamps": {}, "bundle_timestamp": bundle_ts}

    real = discover_dataset(results, "real", bundle_ts)
    if real and real.is_dir():
        data["timestamps"]["real"] = real.name
        data["real_datasets"] = load_perf_dir(real)
        zf = pick_zfasta(data["real_datasets"].get("Genome", []), "z-fasta-full")
        if zf:
            data["headline_genome_stats_full_s"] = zf["mean_s"]

    indexonly = discover_dataset(results, "indexonly", bundle_ts)
    if indexonly and indexonly.is_dir():
        data["timestamps"]["indexonly"] = indexonly.name
        all_io = load_perf_dir(indexonly)
        data["indexonly_mb"] = {k: v for k, v in all_io.items() if k in ("10mb", "100mb")}

    stats_tests = discover_dataset(results, "stats", bundle_ts)
    if stats_tests and stats_tests.is_dir():
        data["timestamps"]["test_files"] = stats_tests.name
        data["test_files"] = load_perf_dir(stats_tests)

    mem = discover_dataset(results, "memory", bundle_ts)
    if mem and mem.suffix == ".csv":
        data["timestamps"]["memory"] = mem.name
        data["memory"] = load_memory_csv(mem)

    data["bundle_warnings"] = bundle_warnings(
        "stats", results, STATS_BUNDLE_PREFIXES, bundle_ts, data["timestamps"]
    )
    return data


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Save a performance baseline snapshot.")
    parser.add_argument(
        "--label",
        help="Optional suffix for the baseline directory name",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Override output directory (default: bench/baselines/<timestamp>_<rev>)",
    )
    args = parser.parse_args()

    git = git_info()
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    rev = git.get("commit_short") or "unknown"
    if git.get("dirty"):
        rev += "-dirty"

    if args.out_dir:
        out_dir = args.out_dir
    else:
        name = f"{ts}_{rev}"
        if args.label:
            safe = re.sub(r"[^A-Za-z0-9._-]+", "-", args.label)
            name += f"_{safe}"
        out_dir = BASELINES_DIR / name

    index_data = collect_index()
    get_data = collect_get()
    stats_data = collect_stats()

    if not index_data.get("real_datasets"):
        print(
            "ERROR: no index perf results found. Run bench/index/run_benchmarks.sh first.",
            file=sys.stderr,
        )
        return 1

    all_warnings = (
        index_data.get("bundle_warnings", [])
        + get_data.get("bundle_warnings", [])
        + stats_data.get("bundle_warnings", [])
    )

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git": git,
        "machine": machine_info(),
        "z_fasta": zfasta_binary_info(),
        "suites": {
            "index": {
                "bundle_timestamp": index_data.get("bundle_timestamp"),
                "artifacts": index_data.get("timestamps", {}),
            },
            "get": {
                "bundle_timestamp": get_data.get("bundle_timestamp"),
                "artifacts": get_data.get("timestamps", {}),
            },
            "stats": {
                "bundle_timestamp": stats_data.get("bundle_timestamp"),
                "artifacts": stats_data.get("timestamps", {}),
            },
        },
        "headlines": {
            "genome_index_s": index_data.get("headline_genome_index_s"),
            "genome_stats_full_s": stats_data.get("headline_genome_stats_full_s"),
        },
        "bundle_warnings": all_warnings,
    }

    compare: dict[str, float] = {}
    compare.update(flatten_compare("index", index_data))
    compare.update(flatten_compare("get", get_data))
    compare.update(flatten_compare("stats", stats_data))

    write_json(out_dir / "manifest.json", manifest)
    write_json(out_dir / "index.json", index_data)
    write_json(out_dir / "get.json", get_data)
    write_json(out_dir / "stats.json", stats_data)
    write_json(out_dir / "compare.json", compare)

    latest_ptr = BASELINES_DIR / "LATEST"
    latest_ptr.write_text(str(out_dir.name) + "\n")

    print(f"Baseline saved: {out_dir}")
    if manifest["headlines"].get("genome_index_s") is not None:
        print(f"  genome index: {manifest['headlines']['genome_index_s']:.4f} s")
    if all_warnings:
        print("Warnings (mixed or partial bench runs):", file=sys.stderr)
        for warning in all_warnings:
            print(f"  - {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
