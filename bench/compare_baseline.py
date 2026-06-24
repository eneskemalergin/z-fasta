#!/usr/bin/env python3
"""
Compare two bench/baselines snapshots (compare.json flat maps).

Usage:
    .venv/bin/python bench/compare_baseline.py bench/baselines/<old> bench/baselines/<new>
    .venv/bin/python bench/compare_baseline.py --threshold 10
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load_compare(path: Path) -> dict[str, float]:
    with open(path) as f:
        return json.load(f)


def load_manifest(path: Path) -> dict | None:
    manifest_path = path / "manifest.json"
    if not manifest_path.is_file():
        return None
    with open(manifest_path) as f:
        return json.load(f)


def print_machine_note(a_dir: Path, b_dir: Path) -> None:
    ma = load_manifest(a_dir)
    mb = load_manifest(b_dir)
    if not ma or not mb:
        return
    host_a = ma.get("machine", {}).get("hostname")
    host_b = mb.get("machine", {}).get("hostname")
    cpu_a = ma.get("machine", {}).get("cpu")
    cpu_b = mb.get("machine", {}).get("cpu")
    if host_a and host_b and (host_a != host_b or cpu_a != cpu_b):
        print(
            "Note: baselines were captured on different hosts or CPUs; "
            "timing deltas are indicative only.\n"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two baseline snapshots.")
    parser.add_argument("baseline_a", type=Path, help="Older baseline directory")
    parser.add_argument("baseline_b", type=Path, help="Newer baseline directory")
    parser.add_argument(
        "--threshold",
        type=float,
        default=5.0,
        help="Regression threshold in percent (default: 5)",
    )
    parser.add_argument(
        "--z-fasta-only",
        action="store_true",
        help="Only show z-fasta tool rows",
    )
    args = parser.parse_args()

    a = load_compare(args.baseline_a / "compare.json")
    b = load_compare(args.baseline_b / "compare.json")

    print(f"Comparing:\n  A: {args.baseline_a}\n  B: {args.baseline_b}\n")
    print_machine_note(args.baseline_a, args.baseline_b)

    keys = sorted(set(a) | set(b))
    regressions: list[tuple[str, float, float, float]] = []
    improvements: list[tuple[str, float, float, float]] = []
    missing_in_a: list[str] = []
    missing_in_b: list[str] = []

    print(f"{'key':<80} {'A (s)':>10} {'B (s)':>10} {'delta%':>10}")
    print("-" * 112)

    for key in keys:
        if args.z_fasta_only and "z-fasta" not in key:
            continue
        va = a.get(key)
        vb = b.get(key)
        if va is None:
            missing_in_a.append(key)
            continue
        if vb is None:
            missing_in_b.append(key)
            continue
        if va <= 0:
            continue
        delta_pct = ((vb - va) / va) * 100.0
        print(f"{key:<80} {va:10.4f} {vb:10.4f} {delta_pct:+10.1f}")
        if delta_pct > args.threshold:
            regressions.append((key, va, vb, delta_pct))
        elif delta_pct < -args.threshold:
            improvements.append((key, va, vb, delta_pct))

    print()
    if missing_in_a:
        print(f"Keys only in B ({len(missing_in_a)}):")
        for key in missing_in_a[:10]:
            print(f"  {key}")
        if len(missing_in_a) > 10:
            print(f"  ... and {len(missing_in_a) - 10} more")
        print()
    if missing_in_b:
        print(f"Keys only in A ({len(missing_in_b)}):")
        for key in missing_in_b[:10]:
            print(f"  {key}")
        if len(missing_in_b) > 10:
            print(f"  ... and {len(missing_in_b) - 10} more")
        print()

    if regressions:
        print(f"REGRESSIONS (>{args.threshold:.0f}% slower):")
        for key, va, vb, d in sorted(regressions, key=lambda x: -x[3]):
            print(f"  {key}: {va:.4f}s → {vb:.4f}s ({d:+.1f}%)")
    else:
        print(f"No regressions above {args.threshold:.0f}% threshold.")

    if improvements:
        print(f"\nIMPROVEMENTS (>{args.threshold:.0f}% faster):")
        for key, va, vb, d in sorted(improvements, key=lambda x: x[3])[:10]:
            print(f"  {key}: {va:.4f}s → {vb:.4f}s ({d:+.1f}%)")

    return 1 if regressions else 0


if __name__ == "__main__":
    raise SystemExit(main())
