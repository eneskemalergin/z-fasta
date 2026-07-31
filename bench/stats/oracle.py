#!/usr/bin/env python3
"""Stats oracles for bench/stats/run.sh correctness (run_tests) only. Extend in place; no sibling scripts."""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from collections.abc import Callable
from pathlib import Path

from Bio import SeqIO

IUPAC_NUC = "ACGTURYSWKMBDHVNacgturyswkmbdhvn"

AA_CODES = "ARNDCEQGHILKMFPSTWYV"
AA_NAMES = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartate",
    "C": "Cysteine",
    "E": "Glutamate",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
}


def format_size(num_bytes: int) -> str:
    gb = num_bytes / (1024**3)
    mb = num_bytes / (1024**2)
    kb = num_bytes / 1024
    if gb >= 1.0:
        return f"{gb:.1f} GB"
    if mb >= 1.0:
        return f"{mb:.1f} MB"
    if kb >= 1.0:
        return f"{kb:.1f} KB"
    return f"{num_bytes} B"


def load_indexed_names(fasta: str) -> list[str]:
    fai = f"{fasta}.fai"
    if not Path(fai).is_file():
        return []
    names: list[str] = []
    for line in Path(fai).read_text().splitlines():
        parts = line.split("\t")
        if parts:
            names.append(parts[0])
    return names


def indexed_records(fasta: str, *, dedup: bool = True) -> list:
    all_records = list(SeqIO.parse(fasta, "fasta"))
    if not dedup:
        return [r for r in all_records if len(r.seq) > 0]
    indexed = load_indexed_names(fasta)
    if not indexed:
        return [r for r in all_records if len(r.seq) > 0]
    want = set(indexed)
    seen: set[str] = set()
    out = []
    for rec in all_records:
        if rec.id in want and rec.id not in seen and len(rec.seq) > 0:
            out.append(rec)
            seen.add(rec.id)
    return out


def detect_type(counts: Counter, total: int) -> str:
    if total == 0:
        return "nucleotide"
    nuc = sum(counts.get(c, 0) for c in IUPAC_NUC)
    return "nucleotide" if nuc * 10 > total * 9 else "protein"


def top_amino_acids(counts: Counter, total: int) -> list[dict]:
    entries: list[dict] = []
    for i, code in enumerate(AA_CODES):
        cnt = counts.get(code, 0) + counts.get(code.lower(), 0)
        entries.append({"code": code, "count": cnt, "ord": i})
    # Match stats.zig: sort by count desc; ties keep amino_acid_names order.
    entries.sort(key=lambda e: (-e["count"], e["ord"]))
    out: list[dict] = []
    for entry in entries[:3]:
        out.append(
            {
                "code": entry["code"],
                "pct": entry["count"] / total * 100.0 if total else 0.0,
                "name": AA_NAMES[entry["code"]],
            }
        )
    return out


def name_duplicate_extras(names: list[str]) -> int:
    """Source-level extras: sum(k - 1) for each name occurring k times."""
    return sum(count - 1 for count in Counter(names).values() if count > 1)


def compute_expected(fasta: str, *, dedup: bool = True) -> dict:
    all_records = list(SeqIO.parse(fasta, "fasta"))
    records = indexed_records(fasta, dedup=dedup)
    lengths = [len(r.seq) for r in records]
    num_seqs = len(lengths)
    total_bases = sum(lengths)
    sorted_desc = sorted(lengths, reverse=True)

    mean = total_bases // num_seqs if num_seqs else 0
    if num_seqs == 0:
        median = 0
        median_half_up = 0
    elif num_seqs % 2 == 1:
        median = sorted_desc[num_seqs // 2]
        median_half_up = median
    else:
        central_sum = sorted_desc[num_seqs // 2 - 1] + sorted_desc[num_seqs // 2]
        median = central_sum // 2  # z-fasta integer floor
        median_half_up = (central_sum + 1) // 2  # seqkit-style round half up

    threshold_50 = (total_bases + 1) // 2
    threshold_90 = (total_bases * 9 + 9) // 10
    bases_seen = au_sum = 0
    n50 = l50 = n90 = l90 = 0
    found_n50 = found_n90 = False
    for i, length in enumerate(sorted_desc):
        bases_seen += length
        au_sum += length * length
        if not found_n50 and bases_seen >= threshold_50:
            n50, l50, found_n50 = length, i + 1, True
        if not found_n90 and bases_seen >= threshold_90:
            n90, l90, found_n90 = length, i + 1, True
    au = au_sum // total_bases if total_bases else 0

    shortest_idx = lengths.index(min(lengths))
    longest_idx = lengths.index(max(lengths))

    counts: Counter = Counter()
    lowercase = 0
    comp_total = 0
    for rec in records:
        for ch in str(rec.seq):
            counts[ch] += 1
            comp_total += 1
            if ch.islower():
                lowercase += 1

    seq_type = detect_type(counts, comp_total)
    fasta_bytes = Path(fasta).stat().st_size
    source_duplicates = name_duplicate_extras([r.id for r in all_records])
    # Index-visible extras (--no-dedup path); empty sequences are omitted from indexes.
    index_duplicates = name_duplicate_extras([r.id for r in records])
    out: dict = {
        "dedup_index": dedup,
        "source_duplicates": source_duplicates,
        "index_duplicates": index_duplicates,
        "fasta_bytes": fasta_bytes,
        "fasta_size": format_size(fasta_bytes),
        "fasta_name": Path(fasta).name,
        "num_seqs": num_seqs,
        "total_bases": total_bases,
        "shortest_len": min(lengths) if lengths else 0,
        "shortest_name": records[shortest_idx].id,
        "longest_len": max(lengths) if lengths else 0,
        "longest_name": records[longest_idx].id,
        "mean": mean,
        "median": median,
        "median_half_up": median_half_up,
        "n50": n50,
        "l50": l50,
        "n90": n90,
        "l90": l90,
        "au": au,
        "seq_type": seq_type,
    }

    soft_pct = lowercase / comp_total * 100.0 if comp_total else 0.0
    if seq_type == "nucleotide":
        a = counts.get("A", 0) + counts.get("a", 0)
        c = counts.get("C", 0) + counts.get("c", 0)
        g = counts.get("G", 0) + counts.get("g", 0)
        t = counts.get("T", 0) + counts.get("t", 0)
        n = counts.get("N", 0) + counts.get("n", 0)
        acgt = a + c + g + t
        other = comp_total - a - c - g - t - n
        out.update(
            {
                "a_pct": a / comp_total * 100.0 if comp_total else 0.0,
                "c_pct": c / comp_total * 100.0 if comp_total else 0.0,
                "g_pct": g / comp_total * 100.0 if comp_total else 0.0,
                "t_pct": t / comp_total * 100.0 if comp_total else 0.0,
                "n_pct": n / comp_total * 100.0 if comp_total else 0.0,
                "n_count": n,
                "other_pct": other / comp_total * 100.0 if comp_total and other > 0 else 0.0,
                "has_other": other > 0,
                "gc_pct": (g + c) / acgt * 100.0 if acgt else 0.0,
                "gc_excl_n": (n / comp_total * 100.0 if comp_total else 0.0) > 1.0,
                "gc_skew": (g - c) / (g + c) if (g + c) else 0.0,
                "has_gc_skew": (g + c) > 0,
                "soft_pct": soft_pct,
            }
        )
    else:
        out["soft_pct"] = soft_pct
        out["top_aa"] = top_amino_acids(counts, comp_total)
    return out


def strip_commas(s: str) -> int:
    return int(s.replace(",", ""))


def parse_zfasta(text: str) -> dict:
    out: dict = {"top_aa": []}
    for line in text.splitlines():
        raw = line
        line = line.strip()

        m = re.match(r"File:\s+(.+?)\s+\((.+?) on disk\)", line)
        if m:
            out["file_path"] = m.group(1).strip()
            out["file_size"] = m.group(2).strip()
            continue

        m = re.match(r"Index:\s+(.+)", line)
        if m:
            out["index_path"] = m.group(1).strip()
            continue

        m = re.match(r"^\s*([A-Z]):\s+([\d.]+)%\s+\((.+)\)", line)
        if m:
            out["top_aa"].append({"code": m.group(1), "pct": float(m.group(2)), "name": m.group(3)})
            continue

        if line == "(20 amino acids total)":
            out["protein_footer"] = True
            continue

        for pat, key, conv in (
            (r"Sequences:\s+(.+)", "num_seqs", strip_commas),
            (r"Total bases:\s+(.+)", "total_bases", strip_commas),
            (r"Shortest:\s+([\d,]+)\s+\((.+)\)", None, None),
            (r"Longest:\s+([\d,]+)\s+\((.+)\)", None, None),
            (r"Mean:\s+(.+)", "mean", strip_commas),
            (r"Median:\s+(.+)", "median", strip_commas),
            (r"N50:\s+(.+)", "n50", strip_commas),
            (r"L50:\s+(.+)", "l50", strip_commas),
            (r"N90:\s+(.+)", "n90", strip_commas),
            (r"L90:\s+(.+)", "l90", strip_commas),
            (r"AU:\s+(.+)", "au", strip_commas),
            (r"Type:\s+(.+)", "type", str.strip),
            (r"Duplicates:\s+(\d+)\s*$", "duplicates", int),
            (r"Duplicates:\s+(n/a(?:\s.*)?)\s*$", "duplicates_na", str.strip),
            (r"^\s*A:\s+([\d.]+)%", "a_pct", float),
            (r"^\s*C:\s+([\d.]+)%", "c_pct", float),
            (r"^\s*G:\s+([\d.]+)%", "g_pct", float),
            (r"^\s*T:\s+([\d.]+)%", "t_pct", float),
            (r"^\s*N:\s+([\d.]+)%", "n_pct", float),
            (r"Other:\s+([\d.]+)%", "other_pct", float),
            (r"GC:\s+([\d.]+)%", "gc_pct", float),
            (r"GC skew:\s+([+-]?[\d.]+)", "gc_skew", float),
            (r"Soft-masked:\s+([\d.]+)%", "soft_pct", float),
            (r"Lowercase:\s+([\d.]+)%", "soft_pct", float),
            (r"N content:\s+(.+)", "n_content", strip_commas),
        ):
            m = re.search(pat, raw)
            if not m:
                continue
            if key is None:
                if "Shortest" in pat:
                    out["shortest_len"] = strip_commas(m.group(1))
                    out["shortest_name"] = m.group(2)
                else:
                    out["longest_len"] = strip_commas(m.group(1))
                    out["longest_name"] = m.group(2)
            else:
                out[key] = conv(m.group(1))

        if "GC:" in line and "(excl. N)" in line:
            out["gc_excl_n"] = True
    return out


def parse_wrapper(text: str) -> dict:
    out: dict = {"top_aa": []}
    for line in text.splitlines():
        key, _, val = line.partition("\t")
        if not key:
            continue
        if key in (
            "sequences",
            "total_bases",
            "shortest_len",
            "longest_len",
            "mean",
            "median",
            "n50",
            "l50",
            "n90",
            "l90",
            "au",
            "n_content",
        ):
            out[key] = int(val)
        elif key in (
            "a_pct",
            "c_pct",
            "g_pct",
            "t_pct",
            "n_pct",
            "other_pct",
            "gc_pct",
            "gc_skew",
            "soft_pct",
        ):
            out[key] = float(val)
        elif key in ("shortest_name", "longest_name", "type"):
            out[key] = val
        elif key.startswith("top_aa_"):
            # code:pct:name
            parts = val.split(":", 2)
            if len(parts) == 3:
                out["top_aa"].append(
                    {"code": parts[0], "pct": float(parts[1]), "name": parts[2]}
                )
    return out


def compare_wrapper_to_expected(tool: str, expected: dict, got: dict, errors: list[str]) -> None:
    """Compare wrapper TSV stats to BioPython expected (z-fasta field formulas).

    Wrappers are clean-FASTA comparison peers only. They re-parse the FASTA with
    noodles/rust-bio; they do not strip messy whitespace or use side tables.
    Require the expanded field set so a stale 2-line binary fails loudly.
    """
    required = (
        "sequences",
        "total_bases",
        "shortest_len",
        "longest_len",
        "mean",
        "median",
        "n50",
        "l50",
        "n90",
        "l90",
        "au",
        "type",
        "soft_pct",
    )
    missing = [k for k in required if k not in got]
    if missing:
        errors.append(
            f"{tool}: stale or incomplete stats output (missing {', '.join(missing)}); "
            f"rebuild tools/{tool}_wrapper --target-dir ./target"
        )
        return

    int_ok(errors, f"{tool}.sequences", expected["num_seqs"], got.get("sequences", -2))
    int_ok(errors, f"{tool}.total_bases", expected["total_bases"], got.get("total_bases", -2))
    for key in (
        "shortest_len",
        "longest_len",
        "mean",
        "median",
        "n50",
        "l50",
        "n90",
        "l90",
        "au",
    ):
        int_ok(errors, f"{tool}.{key}", expected[key], got.get(key, -2))
    str_ok(errors, f"{tool}.shortest_name", expected["shortest_name"], got.get("shortest_name", "?"))
    str_ok(errors, f"{tool}.longest_name", expected["longest_name"], got.get("longest_name", "?"))
    want_type = "nucleotide" if expected["seq_type"] == "nucleotide" else "protein"
    str_ok(errors, f"{tool}.type", want_type, got.get("type", "?"))
    float_ok(errors, f"{tool}.soft_pct", expected["soft_pct"], got.get("soft_pct", -1.0))

    if expected["seq_type"] == "nucleotide":
        for key in ("a_pct", "c_pct", "g_pct", "t_pct", "n_pct", "gc_pct"):
            if key not in got:
                errors.append(f"{tool}.{key}: missing from wrapper output")
                continue
            float_ok(errors, f"{tool}.{key}", expected[key], got.get(key, -1.0))
        if expected.get("has_other"):
            float_ok(errors, f"{tool}.other_pct", expected["other_pct"], got.get("other_pct", -1.0))
        if expected.get("has_gc_skew"):
            float_ok(errors, f"{tool}.gc_skew", expected["gc_skew"], got.get("gc_skew", 9.0), tol=0.002)
        int_ok(errors, f"{tool}.n_content", expected["n_count"], got.get("n_content", -1))
    else:
        got_top = got.get("top_aa", [])
        if len(got_top) < 3:
            errors.append(f"{tool}.top_aa: expected 3 entries, got {len(got_top)}")
        for i, want in enumerate(expected.get("top_aa", [])):
            if i >= len(got_top):
                errors.append(f"{tool}.top_aa[{i}]: missing")
                continue
            float_ok(errors, f"{tool}.top_aa[{i}].pct", want["pct"], got_top[i].get("pct", -1.0))
            str_ok(errors, f"{tool}.top_aa[{i}].code", want["code"], got_top[i].get("code", "?"))
            str_ok(errors, f"{tool}.top_aa[{i}].name", want["name"], got_top[i].get("name", "?"))


def parse_int_field(s: str) -> int:
    return int(s.replace(",", ""))


def parse_seqkit_assembly(text: str) -> dict:
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("file")]
    if not lines:
        return {}
    # Prefer TSV (-T): header + data. Fall back to whitespace table.
    raw_lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(raw_lines) >= 2 and "\t" in raw_lines[0]:
        headers = raw_lines[0].split("\t")
        cols = raw_lines[-1].split("\t")
        if len(cols) < len(headers):
            return {}
        by = dict(zip(headers, cols))
        out: dict = {
            "num_seqs": parse_int_field(by.get("num_seqs", "-1")),
            "total_bases": parse_int_field(by.get("sum_len", "-1")),
            "n50": parse_int_field(by.get("N50", "-1")),
            "gc_pct": float(by.get("GC(%)", "-2")),
        }
        if "Q2" in by:
            out["q2"] = parse_int_field(by["Q2"])
        return out

    cols = lines[-1].split()
    if len(cols) < 18:
        return {}
    return {
        "num_seqs": parse_int_field(cols[3]),
        "total_bases": parse_int_field(cols[4]),
        "n50": parse_int_field(cols[12]),
        "gc_pct": float(cols[17]),
        # Whitespace table: Q2 is column 10 (0-based index 9) when -a is used.
        "q2": parse_int_field(cols[9]) if len(cols) > 9 else -1,
    }


def seqkit_q2_vs_median(expected: dict, q2: int, errors: list[str]) -> None:
    """Compare seqkit Q2 to z-fasta median.

    Odd N: Q2 must equal median.
    Even N: both average the two central lengths; z-fasta floors, seqkit rounds half up.
    Example proteome lengths 20 and 51: midpoint 35.5 -> z-fasta 35, seqkit Q2 36.
    """
    median = expected["median"]
    half_up = expected.get("median_half_up", median)
    n = expected["num_seqs"]
    if n <= 0:
        errors.append("seqkit.q2: expected num_seqs > 0")
        return
    if n % 2 == 1:
        int_ok(errors, "seqkit.q2_vs_median", median, q2)
        return
    int_ok(errors, "seqkit.q2_vs_median_half_up", half_up, q2)
    if half_up != median and q2 == median:
        errors.append(
            f"seqkit.q2_vs_median: got floor {median} but seqkit should half-up to {half_up}"
        )


def parse_seqtk_comp(text: str, indexed_names: list[str]) -> dict:
    want = set(indexed_names) if indexed_names else None
    seen: set[str] = set()
    a = c = g = t = n = total = 0
    for line in text.splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        name = cols[0]
        if want is not None:
            if name not in want or name in seen:
                continue
            seen.add(name)
        if len(cols) < 6:
            continue
        length = int(cols[1])
        if length == 0:
            continue
        a += int(cols[2])
        c += int(cols[3])
        g += int(cols[4])
        t += int(cols[5])
        n += int(cols[6]) if len(cols) > 6 else 0
        total += length
    acgt = a + c + g + t
    return {
        "total_bases": total,
        "a_pct": a / total * 100.0 if total else 0.0,
        "c_pct": c / total * 100.0 if total else 0.0,
        "g_pct": g / total * 100.0 if total else 0.0,
        "t_pct": t / total * 100.0 if total else 0.0,
        "n_pct": n / total * 100.0 if total else 0.0,
        "n_count": n,
        "gc_pct": (g + c) / acgt * 100.0 if acgt else 0.0,
    }


def fail(errors: list[str]) -> None:
    for line in errors:
        print(line, file=sys.stderr)
    sys.exit(1)


def int_ok(errors: list[str], label: str, expected: int, got: int) -> None:
    if expected != got:
        errors.append(f"{label}: expected {expected}, got {got}")


def str_ok(errors: list[str], label: str, expected: str, got: str) -> None:
    if expected != got:
        errors.append(f"{label}: expected {expected!r}, got {got!r}")


def float_ok(errors: list[str], label: str, expected: float, got: float, tol: float = 0.015) -> None:
    if abs(expected - got) > tol:
        errors.append(f"{label}: expected {expected:.4f}, got {got:.4f}")


def check_header(expected: dict, got: dict, errors: list[str], index_tag: str) -> None:
    want_suffix = ".zfi" if index_tag == "zfi" else ".fai"
    index_path = got.get("index_path", "")
    if not index_path.endswith(want_suffix):
        errors.append(f"header.index: expected path ending in {want_suffix!r}, got {index_path!r}")
    if expected["fasta_name"] not in got.get("file_path", ""):
        errors.append("header.file: fasta name missing from File line")
    str_ok(errors, "header.file_size", expected["fasta_size"], got.get("file_size", "?"))


def compare_duplicates(expected: dict, got: dict, errors: list[str], *, index_only: bool) -> None:
    """Product policy: full stats report source extras; index-only never fabricates 0."""
    if "source_duplicates" not in expected:
        return
    if index_only:
        # A number is knowable only when the index retains repeated names.
        # Unique index names → n/a (dedup may have dropped repeats, or source
        # may only duplicate empty records the indexer omits).
        if expected.get("index_duplicates", 0) > 0:
            int_ok(errors, "duplicates", expected["index_duplicates"], got.get("duplicates", -1))
        elif "duplicates_na" not in got:
            errors.append(
                f"duplicates: expected n/a under index-only with no retained repeats, got {got.get('duplicates', '?')!r}"
            )
        return
    int_ok(errors, "duplicates", expected["source_duplicates"], got.get("duplicates", -1))


def compare_index(expected: dict, got: dict, errors: list[str], *, index_only: bool = False) -> None:
    for key in ("num_seqs", "total_bases", "shortest_len", "longest_len", "mean", "median", "n50", "l50", "n90", "l90", "au"):
        int_ok(errors, key, expected[key], got.get(key, -1))
    str_ok(errors, "shortest_name", expected["shortest_name"], got.get("shortest_name", "?"))
    str_ok(errors, "longest_name", expected["longest_name"], got.get("longest_name", "?"))
    compare_duplicates(expected, got, errors, index_only=index_only)


def compare_full(expected: dict, got: dict, errors: list[str]) -> None:
    compare_index(expected, got, errors, index_only=False)
    want_type = "Nucleotide" if expected["seq_type"] == "nucleotide" else "Protein"
    str_ok(errors, "type", want_type, got.get("type", "?"))

    if expected["seq_type"] == "nucleotide":
        for key in ("a_pct", "c_pct", "g_pct", "t_pct", "n_pct", "gc_pct", "soft_pct"):
            float_ok(errors, key, expected[key], got.get(key, -1.0))
        if expected.get("has_other"):
            float_ok(errors, "other_pct", expected["other_pct"], got.get("other_pct", -1.0))
        elif got.get("other_pct") is not None:
            errors.append("other_pct: expected absent, got value")
        if expected.get("gc_excl_n"):
            if not got.get("gc_excl_n"):
                errors.append("gc_excl_n: expected '(excl. N)' suffix on GC line")
        if expected.get("has_gc_skew"):
            float_ok(errors, "gc_skew", expected["gc_skew"], got.get("gc_skew", 9.0), tol=0.002)
        elif got.get("gc_skew") is not None:
            errors.append("gc_skew: expected absent when G+C is zero")
        int_ok(errors, "n_content", expected["n_count"], got.get("n_content", -1))
    else:
        float_ok(errors, "soft_pct", expected["soft_pct"], got.get("soft_pct", -1.0))
        if not got.get("protein_footer"):
            errors.append("protein_footer: expected '(20 amino acids total)' line")
        got_top = got.get("top_aa", [])
        for i, want in enumerate(expected.get("top_aa", [])):
            if i >= len(got_top):
                errors.append(f"top_aa[{i}]: missing")
                continue
            float_ok(errors, f"top_aa[{i}].pct", want["pct"], got_top[i].get("pct", -1.0))
            str_ok(errors, f"top_aa[{i}].code", want["code"], got_top[i].get("code", "?"))
            str_ok(errors, f"top_aa[{i}].name", want["name"], got_top[i].get("name", "?"))


def cross_duplicates_ok(a: dict, b: dict, errors: list[str]) -> None:
    if "duplicates_na" in a or "duplicates_na" in b:
        if "duplicates_na" not in a or "duplicates_na" not in b:
            errors.append(
                f"cross.duplicates: n/a mismatch ({a.get('duplicates_na', a.get('duplicates'))!r} vs {b.get('duplicates_na', b.get('duplicates'))!r})"
            )
        return
    int_ok(errors, "cross.duplicates", a.get("duplicates", -1), b.get("duplicates", -2))


def compare_parsed_full(a: dict, b: dict, errors: list[str]) -> None:
    compare_index(
        {k: a.get(k, -1) for k in ("num_seqs", "total_bases", "shortest_len", "longest_len", "mean", "median", "n50", "l50", "n90", "l90", "au")}
        | {"shortest_name": a.get("shortest_name", "?"), "longest_name": a.get("longest_name", "?")},
        b,
        errors,
    )
    cross_duplicates_ok(a, b, errors)
    str_ok(errors, "cross.type", a.get("type", "?"), b.get("type", "!"))
    if a.get("type") == "Nucleotide":
        for key in ("a_pct", "c_pct", "g_pct", "t_pct", "n_pct", "gc_pct", "soft_pct", "other_pct"):
            if key in a or key in b:
                float_ok(errors, f"cross.{key}", a.get(key, -1.0), b.get(key, -2.0))
        if "gc_skew" in a or "gc_skew" in b:
            float_ok(errors, "cross.gc_skew", a.get("gc_skew", 0.0), b.get("gc_skew", 9.0), tol=0.002)
        int_ok(errors, "cross.n_content", a.get("n_content", -1), b.get("n_content", -2))
    else:
        float_ok(errors, "cross.soft_pct", a.get("soft_pct", -1.0), b.get("soft_pct", -2.0))
        for i in range(min(len(a.get("top_aa", [])), len(b.get("top_aa", [])))):
            float_ok(errors, f"cross.top_aa[{i}].pct", a["top_aa"][i]["pct"], b["top_aa"][i].get("pct", -2.0))


def cmd_expected(argv: list[str]) -> None:
    fasta = argv[2]
    dedup = not (len(argv) > 3 and argv[3] == "--no-dedup")
    print(json.dumps(compute_expected(fasta, dedup=dedup)))


def cmd_check(argv: list[str]) -> None:
    mode, index_tag, fasta_path, exp_path, stats_path = argv[2], argv[3], argv[4], argv[5], argv[6]
    expected = json.loads(Path(exp_path).read_text())
    text = Path(stats_path).read_text()
    got = parse_zfasta(text)
    errors: list[str] = []

    check_header(expected, got, errors, index_tag)

    if mode == "index-only":
        if "Composition:" in text:
            errors.append("index-only: composition section present")
        if "run without --index-only" not in got.get("type", ""):
            errors.append(f"index-only: bad type placeholder ({got.get('type', '?')!r})")
        compare_index(expected, got, errors, index_only=True)
    else:
        if "Composition:" not in text:
            errors.append("full: composition section missing")
        compare_full(expected, got, errors)

    if errors:
        fail(errors)


def cmd_verify_mode(argv: list[str]) -> None:
    mode, fasta_path, exp_path, zfi_path, fai_path = argv[2], argv[3], argv[4], argv[5], argv[6]
    expected = json.loads(Path(exp_path).read_text())
    errors: list[str] = []
    index_only = mode == "index-only"

    for tag, path in (("zfi", zfi_path), ("fai", fai_path)):
        text = Path(path).read_text()
        got = parse_zfasta(text)
        check_header(expected, got, errors, tag)
        if index_only:
            if "Composition:" in text:
                errors.append(f"{tag}: composition section present")
            if "run without --index-only" not in got.get("type", ""):
                errors.append(f"{tag}: bad index-only type placeholder")
            compare_index(expected, got, errors, index_only=True)
        else:
            if "Composition:" not in text:
                errors.append(f"{tag}: composition section missing")
            compare_full(expected, got, errors)

    zfi = parse_zfasta(Path(zfi_path).read_text())
    fai = parse_zfasta(Path(fai_path).read_text())
    if index_only:
        compare_index(
            {k: zfi.get(k, -1) for k in ("num_seqs", "total_bases", "shortest_len", "longest_len", "mean", "median", "n50", "l50", "n90", "l90", "au")}
            | {"shortest_name": zfi.get("shortest_name", "?"), "longest_name": zfi.get("longest_name", "?")},
            fai,
            errors,
        )
        cross_duplicates_ok(zfi, fai, errors)
    else:
        compare_parsed_full(zfi, fai, errors)

    if errors:
        fail(errors)


def cmd_parity(argv: list[str]) -> None:
    fasta, exp_path, zfi_idx, zfi_full = argv[2], argv[3], argv[4], argv[5]
    expected = json.loads(Path(exp_path).read_text())
    errors: list[str] = []
    pos = 6
    tools_validated = 0
    while pos < len(argv):
        tool, path = argv[pos], argv[pos + 1]
        pos += 2
        if tool in ("noodles", "rustbio"):
            got = parse_wrapper(Path(path).read_text())
            compare_wrapper_to_expected(tool, expected, got, errors)
            tools_validated += 1
        elif tool == "seqkit":
            zf = parse_zfasta(Path(zfi_full).read_text())
            sk = parse_seqkit_assembly(Path(path).read_text())
            int_ok(errors, "seqkit.num_seqs", expected["num_seqs"], sk.get("num_seqs", -1))
            int_ok(errors, "seqkit.total_bases", expected["total_bases"], sk.get("total_bases", -1))
            int_ok(errors, "seqkit.n50", expected["n50"], sk.get("n50", -1))
            if "q2" not in sk:
                errors.append("seqkit.q2: missing (run seqkit stats -a)")
            else:
                seqkit_q2_vs_median(expected, sk["q2"], errors)
            # Cross-check z-fasta printed median against BioPython expected.
            int_ok(errors, "seqkit.zfasta_median", expected["median"], zf.get("median", -1))
            if expected["seq_type"] == "nucleotide" and "gc_pct" in zf:
                float_ok(errors, "seqkit.gc_pct", zf["gc_pct"], sk.get("gc_pct", -2.0))
            tools_validated += 1
        elif tool == "seqtk":
            if expected["seq_type"] != "nucleotide":
                continue
            zf = parse_zfasta(Path(zfi_full).read_text())
            comp = parse_seqtk_comp(Path(path).read_text(), load_indexed_names(fasta))
            for key in ("a_pct", "c_pct", "g_pct", "t_pct", "n_pct", "gc_pct"):
                float_ok(errors, f"seqtk.{key}", expected[key], comp.get(key, -1.0))
            int_ok(errors, "seqtk.n_content", expected["n_count"], comp.get("n_count", -1))
            float_ok(errors, "seqtk.gc_pct_vs_zfasta", zf.get("gc_pct", -1.0), comp.get("gc_pct", -2.0))
            tools_validated += 1
        else:
            errors.append(f"parity: unknown tool {tool!r}")
    if tools_validated == 0:
        fail(["parity: no external tools validated"])
    if errors:
        fail(errors)


def cmd_same(argv: list[str]) -> None:
    mode, a_path, b_path = argv[2], argv[3], argv[4]
    a = parse_zfasta(Path(a_path).read_text())
    b = parse_zfasta(Path(b_path).read_text())
    errors: list[str] = []

    if mode == "index-only":
        compare_index(
            {k: a.get(k, -1) for k in ("num_seqs", "total_bases", "shortest_len", "longest_len", "mean", "median", "n50", "l50", "n90", "l90", "au")}
            | {"shortest_name": a.get("shortest_name", "?"), "longest_name": a.get("longest_name", "?")},
            b,
            errors,
        )
        cross_duplicates_ok(a, b, errors)
    else:
        compare_parsed_full(a, b, errors)

    if errors:
        fail(errors)


COMMANDS: dict[str, Callable[[list[str]], None]] = {
    "expected": cmd_expected,
    "check": cmd_check,
    "verify-mode": cmd_verify_mode,
    "same": cmd_same,
    "parity": cmd_parity,
}


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print("usage: oracle.py <expected|check|verify-mode|same|parity> ...", file=sys.stderr)
        sys.exit(2)
    COMMANDS[sys.argv[1]](sys.argv)


if __name__ == "__main__":
    main()
