#!/usr/bin/env python3
"""
verify_stats.py — Verify z-fasta stats output against BioPython.

Uses BioPython (>=1.80) to independently compute all stats and compares
against z-fasta stats output. Runs on non-REAL test files only.

Usage:
    python bench/stats/verify_stats.py
    # Or with explicit venv:
    .venv/bin/python bench/stats/verify_stats.py
"""

import math
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Config — paths are relative to project root (script resolves its own location)
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCH_ROOT = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(BENCH_ROOT)

ZFASTA = os.path.join(PROJECT_ROOT, "zig-out", "bin", "z-fasta")
TEST_DIR = os.path.join(PROJECT_ROOT, "tests", "data")
TEST_FILES = [
    os.path.join(TEST_DIR, "simple.fasta"),
    os.path.join(TEST_DIR, "proteome.fasta"),
    os.path.join(TEST_DIR, "single.fasta"),
    os.path.join(TEST_DIR, "edge_cases.fasta"),
    os.path.join(TEST_DIR, "mixed_widths.fasta"),
]

PASS_COUNT = 0
FAIL_COUNT = 0


def passed(label: str):
    global PASS_COUNT
    PASS_COUNT += 1
    print(f"  PASS: {label}")


def failed(label: str, expected, got):
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(f"  FAIL: {label} — expected {expected}, got {got}")


# ---------------------------------------------------------------------------
# Independent stats computation via BioPython
# ---------------------------------------------------------------------------

def load_indexed_names(fasta_path: str) -> list[str]:
    """Load sequence names from .fai to know which sequences the index contains.

    Both samtools and z-fasta skip empty sequences and deduplicate names.
    BioPython sees all raw records. We filter to match what the index has.
    """
    fai_path = fasta_path + ".fai"
    if not os.path.isfile(fai_path):
        return []  # Fall back to all records
    names = []
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                names.append(parts[0])
    return names


def compute_stats(fasta_path: str) -> dict:
    """Compute all stats independently using BioPython."""
    all_records = list(SeqIO.parse(fasta_path, "fasta"))
    indexed_names = load_indexed_names(fasta_path)

    if indexed_names:
        # Filter to only sequences present in the index (by name, first occurrence)
        seen = set()
        records = []
        for r in all_records:
            if r.id in indexed_names and r.id not in seen and len(r.seq) > 0:
                records.append(r)
                seen.add(r.id)
    else:
        records = [r for r in all_records if len(r.seq) > 0]

    num_seqs = len(records)
    lengths = [len(r.seq) for r in records]
    total_bases = sum(lengths)

    # Sort descending for N50/L50/N90/L90/AU
    sorted_desc = sorted(lengths, reverse=True)

    # Mean (integer division, matching Zig)
    mean = total_bases // num_seqs if num_seqs > 0 else 0

    # Median (matching Zig: integer average for even count)
    if num_seqs == 0:
        median = 0
    elif num_seqs % 2 == 1:
        median = sorted_desc[num_seqs // 2]
    else:
        median = (sorted_desc[num_seqs // 2 - 1] + sorted_desc[num_seqs // 2]) // 2

    # N50/L50/N90/L90/AU — matching Zig's ceiling thresholds
    threshold_50 = (total_bases + 1) // 2  # ceiling
    threshold_90 = (total_bases * 9 + 9) // 10  # ceiling

    bases_seen = 0
    au_sum = 0
    n50 = l50 = n90 = l90 = 0
    found_n50 = found_n90 = False

    for i, length in enumerate(sorted_desc):
        bases_seen += length
        au_sum += length * length

        if not found_n50 and bases_seen >= threshold_50:
            n50 = length
            l50 = i + 1
            found_n50 = True
        if not found_n90 and bases_seen >= threshold_90:
            n90 = length
            l90 = i + 1
            found_n90 = True

    au = au_sum // total_bases if total_bases > 0 else 0

    # Shortest / longest with names
    shortest_idx = lengths.index(min(lengths))
    longest_idx = lengths.index(max(lengths))
    shortest_name = records[shortest_idx].id
    longest_name = records[longest_idx].id

    # Composition: byte-level counting matching z-fasta
    counts = Counter()
    lowercase_count = 0
    comp_total = 0
    for r in records:
        seq_str = str(r.seq)
        for ch in seq_str:
            counts[ch] += 1
            comp_total += 1
            if ch.islower():
                lowercase_count += 1

    # Detect type (matching z-fasta: >90% ACGTNU = nucleotide)
    nuc_chars = sum(counts.get(c, 0) for c in "ACGTNUacgtnu")
    is_nucleotide = (nuc_chars * 10 > comp_total * 9) if comp_total > 0 else True

    result = {
        "num_seqs": num_seqs,
        "total_bases": total_bases,
        "shortest_len": min(lengths),
        "shortest_name": shortest_name,
        "longest_len": max(lengths),
        "longest_name": longest_name,
        "mean": mean,
        "median": median,
        "n50": n50,
        "l50": l50,
        "n90": n90,
        "l90": l90,
        "au": au,
        "is_nucleotide": is_nucleotide,
        "counts": counts,
        "comp_total": comp_total,
        "lowercase_count": lowercase_count,
    }

    if is_nucleotide:
        a = counts.get("A", 0) + counts.get("a", 0)
        c = counts.get("C", 0) + counts.get("c", 0)
        g = counts.get("G", 0) + counts.get("g", 0)
        t = counts.get("T", 0) + counts.get("t", 0)
        n = counts.get("N", 0) + counts.get("n", 0)
        acgt = a + c + g + t
        other = comp_total - a - c - g - t - n

        gc_pct = (g + c) / acgt * 100.0 if acgt > 0 else 0.0
        gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        soft_pct = lowercase_count / comp_total * 100.0 if comp_total > 0 else 0.0

        result.update({
            "a_count": a, "c_count": c, "g_count": g, "t_count": t, "n_count": n,
            "other_count": other,
            "a_pct": a / comp_total * 100.0 if comp_total > 0 else 0.0,
            "c_pct": c / comp_total * 100.0 if comp_total > 0 else 0.0,
            "g_pct": g / comp_total * 100.0 if comp_total > 0 else 0.0,
            "t_pct": t / comp_total * 100.0 if comp_total > 0 else 0.0,
            "n_pct": n / comp_total * 100.0 if comp_total > 0 else 0.0,
            "gc_pct": gc_pct,
            "gc_skew": gc_skew,
            "soft_pct": soft_pct,
        })
    else:
        soft_pct = lowercase_count / comp_total * 100.0 if comp_total > 0 else 0.0
        result["soft_pct"] = soft_pct

    return result


# ---------------------------------------------------------------------------
# Parse z-fasta stats output
# ---------------------------------------------------------------------------

def strip_commas(s: str) -> int:
    """Remove commas and convert to int."""
    return int(s.replace(",", ""))


def parse_zfasta_output(output: str) -> dict:
    """Parse the z-fasta stats text output into a dict."""
    result = {}
    lines = output.strip().split("\n")

    for line in lines:
        line = line.strip()

        m = re.match(r"Sequences:\s+(.+)", line)
        if m:
            result["num_seqs"] = strip_commas(m.group(1))

        m = re.match(r"Total bases:\s+(.+)", line)
        if m:
            result["total_bases"] = strip_commas(m.group(1))

        m = re.match(r"Shortest:\s+([\d,]+)\s+\((.+)\)", line)
        if m:
            result["shortest_len"] = strip_commas(m.group(1))
            result["shortest_name"] = m.group(2)

        m = re.match(r"Longest:\s+([\d,]+)\s+\((.+)\)", line)
        if m:
            result["longest_len"] = strip_commas(m.group(1))
            result["longest_name"] = m.group(2)

        m = re.match(r"Mean:\s+(.+)", line)
        if m:
            result["mean"] = strip_commas(m.group(1))

        m = re.match(r"Median:\s+(.+)", line)
        if m:
            result["median"] = strip_commas(m.group(1))

        m = re.match(r"N50:\s+(.+)", line)
        if m:
            result["n50"] = strip_commas(m.group(1))

        m = re.match(r"L50:\s+(.+)", line)
        if m:
            result["l50"] = strip_commas(m.group(1))

        m = re.match(r"N90:\s+(.+)", line)
        if m:
            result["n90"] = strip_commas(m.group(1))

        m = re.match(r"L90:\s+(.+)", line)
        if m:
            result["l90"] = strip_commas(m.group(1))

        m = re.match(r"AU:\s+(.+)", line)
        if m:
            result["au"] = strip_commas(m.group(1))

        m = re.match(r"Type:\s+(.+)", line)
        if m:
            result["type"] = m.group(1).strip()

        # Composition lines
        m = re.match(r"A:\s+([\d.]+)%", line)
        if m:
            result["a_pct"] = float(m.group(1))

        m = re.match(r"C:\s+([\d.]+)%", line)
        if m:
            result["c_pct"] = float(m.group(1))

        m = re.match(r"G:\s+([\d.]+)%", line)
        if m:
            result["g_pct"] = float(m.group(1))

        m = re.match(r"T:\s+([\d.]+)%", line)
        if m:
            result["t_pct"] = float(m.group(1))

        m = re.match(r"N:\s+([\d.]+)%", line)
        if m:
            result["n_pct"] = float(m.group(1))

        m = re.match(r"GC:\s+([\d.]+)%", line)
        if m:
            result["gc_pct"] = float(m.group(1))

        m = re.match(r"GC skew:\s+([+-]?[\d.]+)", line)
        if m:
            result["gc_skew"] = float(m.group(1))

        m = re.match(r"Soft-masked:\s+([\d.]+)%", line)
        if m:
            result["soft_pct"] = float(m.group(1))

        m = re.match(r"Lowercase:\s+([\d.]+)%", line)
        if m:
            result["soft_pct"] = float(m.group(1))

        m = re.match(r"N content:\s+(.+)", line)
        if m:
            result["n_content"] = strip_commas(m.group(1))

    return result


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare_float(label: str, expected: float, got: float, tolerance: float = 0.015):
    """Compare floats within tolerance (default allows for rounding at 2dp)."""
    if abs(expected - got) <= tolerance:
        passed(label)
    else:
        failed(label, f"{expected:.4f}", f"{got:.4f}")


def compare_int(label: str, expected: int, got: int):
    if expected == got:
        passed(label)
    else:
        failed(label, expected, got)


def compare_str(label: str, expected: str, got: str):
    if expected == got:
        passed(label)
    else:
        failed(label, expected, got)


def verify_file(fasta_path: str):
    """Run z-fasta stats and compare against BioPython computation."""
    filename = os.path.basename(fasta_path)
    print(f"\n=== {filename} ===")

    # Compute expected values with BioPython
    expected = compute_stats(fasta_path)

    # Run z-fasta stats
    proc = subprocess.run(
        [ZFASTA, "stats", fasta_path],
        capture_output=True, text=True, timeout=10,
    )
    if proc.returncode != 0:
        failed("z-fasta stats execution", "exit 0", f"exit {proc.returncode}: {proc.stderr}")
        return

    got = parse_zfasta_output(proc.stdout)

    # Compare Tier 1 stats (exact integer matches)
    compare_int("num_seqs", expected["num_seqs"], got.get("num_seqs", -1))
    compare_int("total_bases", expected["total_bases"], got.get("total_bases", -1))
    compare_int("shortest_len", expected["shortest_len"], got.get("shortest_len", -1))
    compare_str("shortest_name", expected["shortest_name"], got.get("shortest_name", "?"))
    compare_int("longest_len", expected["longest_len"], got.get("longest_len", -1))
    compare_str("longest_name", expected["longest_name"], got.get("longest_name", "?"))
    compare_int("mean", expected["mean"], got.get("mean", -1))
    compare_int("median", expected["median"], got.get("median", -1))
    compare_int("n50", expected["n50"], got.get("n50", -1))
    compare_int("l50", expected["l50"], got.get("l50", -1))
    compare_int("n90", expected["n90"], got.get("n90", -1))
    compare_int("l90", expected["l90"], got.get("l90", -1))
    compare_int("au", expected["au"], got.get("au", -1))

    # Compare type
    expected_type = "Nucleotide" if expected["is_nucleotide"] else "Protein"
    compare_str("type", expected_type, got.get("type", "?"))

    # Compare Tier 2 composition
    if expected["is_nucleotide"]:
        compare_float("a_pct", expected["a_pct"], got.get("a_pct", -1))
        compare_float("c_pct", expected["c_pct"], got.get("c_pct", -1))
        compare_float("g_pct", expected["g_pct"], got.get("g_pct", -1))
        compare_float("t_pct", expected["t_pct"], got.get("t_pct", -1))
        compare_float("n_pct", expected["n_pct"], got.get("n_pct", -1))
        compare_float("gc_pct", expected["gc_pct"], got.get("gc_pct", -1))
        compare_float("gc_skew", expected["gc_skew"], got.get("gc_skew", -1), tolerance=0.002)
        compare_float("soft_pct", expected["soft_pct"], got.get("soft_pct", -1))
        # N content: compare count
        compare_int("n_content", expected["n_count"], got.get("n_content", -1))
    else:
        compare_float("soft_pct (protein)", expected["soft_pct"], got.get("soft_pct", -1))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("z-fasta stats verification against BioPython")
    print("=" * 50)

    # Check z-fasta binary exists
    if not os.path.isfile(ZFASTA):
        print(f"ERROR: {ZFASTA} not found. Run 'zig build' first.")
        sys.exit(1)

    for fasta_path in TEST_FILES:
        if not os.path.isfile(fasta_path):
            print(f"\nSKIP: {fasta_path} not found")
            continue
        verify_file(fasta_path)

    print(f"\n{'=' * 50}")
    print(f"Results: {PASS_COUNT} passed, {FAIL_COUNT} failed")
    print(f"{'=' * 50}")

    if FAIL_COUNT > 0:
        print("SOME TESTS FAILED")
        sys.exit(1)
    else:
        print("ALL PASSED")


if __name__ == "__main__":
    main()
