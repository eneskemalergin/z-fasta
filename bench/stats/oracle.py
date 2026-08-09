#!/usr/bin/env python3
"""Independent exact oracle for the temporary stats correctness gate."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

IUPAC = b"ACGTURYSWKMBDHVN"
AMBIGUOUS = "RYSWKMBDHV"
PROTEIN_FIELDS = (
    ("a_alanine", "A"), ("r_arginine", "R"), ("n_asparagine", "N"),
    ("d_aspartate", "D"), ("c_cysteine", "C"), ("e_glutamate", "E"),
    ("q_glutamine", "Q"), ("g_glycine", "G"), ("h_histidine", "H"),
    ("i_isoleucine", "I"), ("l_leucine", "L"), ("k_lysine", "K"),
    ("m_methionine", "M"), ("f_phenylalanine", "F"), ("p_proline", "P"),
    ("s_serine", "S"), ("t_threonine", "T"), ("w_tryptophan", "W"),
    ("y_tyrosine", "Y"), ("v_valine", "V"), ("b_asx", "B"),
    ("z_glx", "Z"), ("j_xle", "J"), ("x_unknown", "X"),
    ("u_selenocysteine", "U"), ("o_pyrrolysine", "O"),
)


@dataclass(frozen=True)
class Record:
    name: str
    sequence: bytes


def read_fasta(path: Path) -> list[Record]:
    records: list[Record] = []
    name: str | None = None
    chunks: list[bytes] = []
    for raw in path.read_bytes().splitlines():
        if raw.startswith(b">"):
            if name is not None:
                records.append(Record(name, b"".join(chunks)))
            fields = raw[1:].split()
            name = fields[0].decode("utf-8", errors="replace") if fields else ""
            chunks = []
        elif name is not None:
            chunks.append(bytes(byte for byte in raw if byte > 32))
    if name is not None:
        records.append(Record(name, b"".join(chunks)))
    return records


def retained_records(path: Path, deduplicate: bool = True) -> list[Record]:
    retained: list[Record] = []
    seen: set[str] = set()
    for record in read_fasta(path):
        if not record.sequence or (deduplicate and record.name in seen):
            continue
        seen.add(record.name)
        retained.append(record)
    return retained


def median(values: list[int]) -> int:
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else values[middle - 1] + (values[middle] - values[middle - 1]) // 2


def fixed(numerator: int, denominator: int, places: int) -> str:
    negative = numerator < 0
    magnitude = abs(numerator)
    scale = 10**places
    whole, remainder = divmod(magnitude, denominator)
    fraction = (remainder * scale + denominator // 2) // denominator
    if fraction == scale:
        whole += 1
        fraction = 0
    sign = "-" if negative and (whole or fraction) else ""
    return f"{sign}{whole}.{fraction:0{places}d}"


def letter(counts: Counter[int], code: str) -> int:
    return counts[ord(code)] + counts[ord(code.lower())]


def calculate(path: Path, deduplicate: bool = True) -> dict[str, object]:
    records = retained_records(path, deduplicate)
    if not records:
        raise ValueError(f"no retained records in {path}")
    lengths = [len(record.sequence) for record in records]
    ascending = sorted(lengths)
    total = sum(lengths)
    half = len(lengths) // 2
    q1 = ascending[0] if len(lengths) == 1 else median(ascending[:half])
    q3 = ascending[0] if len(lengths) == 1 else median(ascending[-half:])
    cumulative = 0
    square_sum = 0
    n50 = l50 = n90 = l90 = 0
    for rank, length in enumerate(reversed(ascending), start=1):
        cumulative += length
        square_sum += length * length
        if not l50 and cumulative >= (total * 50 + 99) // 100:
            n50, l50 = length, rank
        if not l90 and cumulative >= (total * 90 + 99) // 100:
            n90, l90 = length, rank
    shortest = min(range(len(records)), key=lambda index: lengths[index])
    longest = max(range(len(records)), key=lambda index: lengths[index])
    counts = Counter(byte for record in records for byte in record.sequence)
    lowercase = sum(value for byte, value in counts.items() if ord("a") <= byte <= ord("z"))
    nucleotide = sum(letter(counts, chr(code)) for code in IUPAC)
    return {
        "records": records, "counts": counts, "count": len(records), "total": total,
        "shortest_length": lengths[shortest], "shortest_name": records[shortest].name,
        "longest_length": lengths[longest], "longest_name": records[longest].name,
        "mean": total // len(records), "q1": q1, "median": median(ascending), "q3": q3,
        "range": lengths[longest] - lengths[shortest], "n50": n50, "l50": l50,
        "n90": n90, "l90": l90, "aun_numerator": square_sum,
        "nucleotide": nucleotide * 10 > total * 9, "lowercase": lowercase,
    }


def composition_lines(values: dict[str, object]) -> list[str]:
    counts = values["counts"]
    assert isinstance(counts, Counter)
    total = int(values["total"])
    lowercase = int(values["lowercase"])
    lines: list[str] = []
    if values["nucleotide"]:
        bases = {code: letter(counts, code) for code in "ACGTUN"}
        ambiguous = sum(letter(counts, code) for code in AMBIGUOUS)
        invalid = total - sum(bases.values()) - ambiguous
        if bases["T"] and bases["U"]:
            type_name = "nucleotide_mixed_tu"
        elif bases["T"]:
            type_name = "nucleotide_t"
        elif bases["U"]:
            type_name = "nucleotide_u"
        else:
            type_name = "nucleotide"
        lines.extend((f"type\t{type_name}", "percent_denominator\ttotal_symbols"))
        for key in "acgtun":
            value = bases[key.upper()]
            lines.append(f"{key}\t{value}\t{fixed(value * 100, total, 2)}")
        lines.append(f"iupac_ambiguous\t{ambiguous}\t{fixed(ambiguous * 100, total, 2)}")
        lines.append(f"invalid\t{invalid}\t{fixed(invalid * 100, total, 2)}")
        canonical = sum(bases[code] for code in "ACGTU")
        lines.append(f"gc\t{fixed((bases['G'] + bases['C']) * 100, canonical, 2) if canonical else 'n/a'}")
        gc_total = bases["G"] + bases["C"]
        lines.append(f"gc_skew\t{fixed(bases['G'] - bases['C'], gc_total, 3) if gc_total else 'n/a'}")
    else:
        lines.extend(("type\tprotein", "percent_denominator\ttotal_symbols"))
        assigned = 0
        for key, code in PROTEIN_FIELDS:
            value = letter(counts, code)
            assigned += value
            lines.append(f"{key}\t{value}\t{fixed(value * 100, total, 2)}")
        stop = counts[ord("*")]
        lines.append(f"stop\t{stop}\t{fixed(stop * 100, total, 2)}")
        invalid = total - assigned - stop
        lines.append(f"invalid\t{invalid}\t{fixed(invalid * 100, total, 2)}")
    lines.append(f"lowercase\t{lowercase}\t{fixed(lowercase * 100, total, 2)}")
    return lines


def shared_lines(values: dict[str, object]) -> list[str]:
    return [
        f"indexed_records\t{values['count']}", f"total_symbols\t{values['total']}",
        f"shortest_length\t{values['shortest_length']}", f"shortest_name\t{values['shortest_name']}",
        f"longest_length\t{values['longest_length']}", f"longest_name\t{values['longest_name']}",
        f"mean\t{values['mean']}", f"q1\t{values['q1']}", f"median\t{values['median']}",
        f"q3\t{values['q3']}", f"range\t{values['range']}", f"n50\t{values['n50']}",
        f"l50\t{values['l50']}", f"n90\t{values['n90']}", f"l90\t{values['l90']}",
        f"aun\t{fixed(int(values['aun_numerator']), int(values['total']), 2)}",
        *composition_lines(values),
    ]


def expected_peer(path: Path) -> str:
    return "\n".join(shared_lines(calculate(path))) + "\n"


def expected_zfasta(path: Path, index_suffix: str, deduplicate: bool = True) -> str:
    values = calculate(path, deduplicate)
    shared = shared_lines(values)
    fields = dict(line.split("\t", 1) for line in shared[:16])
    lines = [
        "File:", f"  path: {path}", f"  index: {path}{index_suffix}", f"  size_bytes: {path.stat().st_size}", "",
        "Lengths:", *[f"  {key}: {fields[key]}" for key in (
            "indexed_records", "total_symbols", "shortest_length", "shortest_name", "longest_length",
            "longest_name", "mean", "q1", "median", "q3", "range")], "", "Nx:",
        *[f"  {key}: {fields[key]}" for key in ("n50", "l50", "n90", "l90", "aun")], "", "Composition:",
    ]
    for line in shared[16:]:
        parts = line.split("\t")
        if len(parts) == 3:
            lines.append(f"  {parts[0]}: {parts[1]} {parts[2]}%")
        elif parts[0] == "gc" and parts[1] != "n/a":
            lines.append(f"  gc: {parts[1]}%")
        else:
            lines.append(f"  {parts[0]}: {parts[1]}")
    return "\n".join(lines) + "\n"


def assert_equal(label: str, expected: str, actual: str) -> None:
    if actual == expected:
        return
    wanted = expected.splitlines()
    got = actual.splitlines()
    for index in range(max(len(wanted), len(got))):
        left = wanted[index] if index < len(wanted) else "<end>"
        right = got[index] if index < len(got) else "<end>"
        if left != right:
            raise SystemExit(f"{label}: line {index + 1}: expected {left!r}, got {right!r}")
    raise SystemExit(f"{label}: output differs")


def assert_fields(label: str, actual: dict[str, object], expected: dict[str, object]) -> None:
    for key, wanted in expected.items():
        got = actual[key]
        if got != wanted:
            raise AssertionError(f"{label} {key}: expected {wanted!r}, got {got!r}")


def composition_fields(values: dict[str, object]) -> dict[str, str]:
    fields: dict[str, str] = {}
    for line in composition_lines(values):
        key, value = line.split("\t", 1)
        fields[key] = value
    return fields


def verify(args: argparse.Namespace) -> None:
    fasta = Path(args.fasta)
    zfi_fasta = Path(args.zfi_fasta) if args.zfi_fasta else fasta
    if zfi_fasta.read_bytes() != fasta.read_bytes():
        raise SystemExit("z-fasta .zfi: FASTA view differs from oracle input")
    assert_equal("z-fasta .zfi", expected_zfasta(zfi_fasta, ".zfi", not args.no_dedup), Path(args.zfi).read_text())
    if args.fai:
        fai_fasta = Path(args.fai_fasta) if args.fai_fasta else fasta
        if fai_fasta.read_bytes() != fasta.read_bytes():
            raise SystemExit("z-fasta .fai: FASTA view differs from oracle input")
        assert_equal("z-fasta .fai", expected_zfasta(fai_fasta, ".fai", not args.no_dedup), Path(args.fai).read_text())
    peer = expected_peer(fasta)
    if args.noodles:
        assert_equal("noodles", peer, Path(args.noodles).read_text())
    if args.rustbio:
        assert_equal("rust-bio", peer, Path(args.rustbio).read_text())


def legacy_command(args: argparse.Namespace) -> None:
    """Serve the permanent runner's small command interface."""
    if args.command == "expected":
        values = calculate(Path(args.fasta), not args.no_dedup)
        print(json.dumps({"num_seqs": values["count"], "no_dedup": args.no_dedup}))
        return
    if args.command == "check":
        metadata = json.loads(Path(args.expected).read_text())
        expected = expected_zfasta(Path(args.fasta), f".{args.format}", not metadata["no_dedup"])
        assert_equal(f"z-fasta .{args.format}", expected, Path(args.output).read_text())
        return
    if args.command == "same":
        left = Path(args.left).read_text().split("Lengths:\n", 1)
        right = Path(args.right).read_text().split("Lengths:\n", 1)
        if len(left) != 2 or len(right) != 2:
            raise SystemExit("cross-format comparison: missing Lengths section")
        assert_equal("cross-format report", left[1], right[1])
        return
    if args.command == "parity":
        fasta = Path(args.fasta)
        if args.tool in ("noodles", "rustbio"):
            assert_equal(args.tool, expected_peer(fasta), Path(args.output).read_text())
        elif args.tool == "seqkit":
            verify_seqkit(fasta, Path(args.output).read_text())
        elif args.tool == "seqtk":
            verify_seqtk(fasta, Path(args.output).read_text())
        else:
            raise SystemExit(f"unsupported peer: {args.tool}")
        return
    raise AssertionError(args.command)


def verify_seqkit(fasta: Path, text: str) -> None:
    lines = [line for line in text.splitlines() if line]
    if len(lines) < 2:
        raise SystemExit("seqkit: missing TSV result")
    fields = dict(zip(lines[0].split("\t"), lines[-1].split("\t"), strict=True))
    values = calculate(fasta)
    records = values["records"]
    assert isinstance(records, list)
    lengths = sorted(len(record.sequence) for record in records)
    middle = len(lengths) // 2
    q2 = lengths[middle] if len(lengths) % 2 else (lengths[middle - 1] + lengths[middle] + 1) // 2
    expected = {
        "num_seqs": int(values["count"]), "sum_len": int(values["total"]),
        "min_len": int(values["shortest_length"]), "max_len": int(values["longest_length"]),
        "Q2": q2, "N50": int(values["n50"]), "N50_num": int(values["l50"]),
    }
    for key, wanted in expected.items():
        got = int(fields.get(key, "-1").replace(",", ""))
        if got != wanted:
            raise SystemExit(f"seqkit {key}: expected {wanted}, got {got}")


def verify_seqtk(fasta: Path, text: str) -> None:
    records = retained_records(fasta)
    wanted_names = {record.name for record in records}
    seen: set[str] = set()
    got = Counter()
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 7 or fields[0] not in wanted_names or fields[0] in seen:
            continue
        seen.add(fields[0])
        got["total"] += int(fields[1])
        for key, value in zip(("A", "C", "G", "T", "N"), fields[2:7], strict=True):
            got[key] += int(value)
    counts = Counter(byte for record in records for byte in record.sequence)
    expected = Counter({"total": sum(len(record.sequence) for record in records)})
    for key in "ACGTN":
        expected[key] = letter(counts, key)
    for key, wanted in expected.items():
        if got[key] != wanted:
            raise SystemExit(f"seqtk {key}: expected {wanted}, got {got[key]}")


def self_test(fixtures: Path, design_fixtures: Path) -> None:
    assert fixed(1 * 100, 32, 2) == "3.13"
    assert fixed(226, 16, 2) == "14.13"
    assert fixed(-30, 32, 3) == "-0.938"
    assembly = calculate(fixtures / "assembly.fasta")
    assert [assembly[key] for key in ("count", "total", "q1", "median", "q3", "n50", "l50", "n90", "l90")] == [4, 10, 1, 2, 3, 3, 2, 2, 3]
    assert calculate(fixtures / "threshold_90.fasta")["nucleotide"] is False
    assert calculate(fixtures / "threshold_above_90.fasta")["nucleotide"] is True

    length_answers = {
        "one_t.fasta": (1, 4, 4, "only", 4, "only", 4, 4, 4, 4, 0, 4, 1, 4, 1, 16),
        "two_u.fasta": (2, 6, 2, "short", 4, "long", 3, 2, 3, 4, 2, 4, 1, 2, 2, 20),
        "odd_ties_mixed.fasta": (
            5, 15, 1, "short_first", 5, "long_first", 3, 1, 3, 5, 4, 5, 2, 1, 4, 61,
        ),
        "even_lengths.fasta": (4, 10, 1, "a", 4, "d", 2, 1, 2, 3, 3, 3, 2, 2, 3, 30),
        "nucleotide_categories.fasta": (
            1, 17, 17, "categories", 17, "categories", 17, 17, 17, 17, 0, 17, 1, 17, 1, 289,
        ),
        "ambiguity_only.fasta": (
            1, 12, 12, "ambiguity", 12, "ambiguity", 12, 12, 12, 12, 0, 12, 1, 12, 1, 144,
        ),
        "gc_zero.fasta": (1, 6, 6, "zero", 6, "zero", 6, 6, 6, 6, 0, 6, 1, 6, 1, 36),
        "exact_90_protein.fasta": (
            1, 10, 10, "boundary", 10, "boundary", 10, 10, 10, 10, 0, 10, 1, 10, 1, 100,
        ),
        "protein_complete.fasta": (
            1, 30, 30, "protein", 30, "protein", 30, 30, 30, 30, 0, 30, 1, 30, 1, 900,
        ),
    }
    length_keys = (
        "count", "total", "shortest_length", "shortest_name", "longest_length",
        "longest_name", "mean", "q1", "median", "q3", "range", "n50", "l50",
        "n90", "l90", "aun_numerator",
    )
    design_values: dict[str, dict[str, object]] = {}
    for name, answer in length_answers.items():
        values = calculate(design_fixtures / name)
        design_values[name] = values
        assert_fields(name, values, dict(zip(length_keys, answer, strict=True)))

    composition_answers = {
        "one_t.fasta": {
            "type": "nucleotide_t", "a": "1\t25.00", "c": "1\t25.00",
            "g": "1\t25.00", "t": "1\t25.00", "u": "0\t0.00", "n": "0\t0.00",
            "iupac_ambiguous": "0\t0.00", "invalid": "0\t0.00", "gc": "50.00",
            "gc_skew": "0.000", "lowercase": "4\t100.00",
        },
        "two_u.fasta": {
            "type": "nucleotide_u", "a": "1\t16.67", "c": "2\t33.33",
            "g": "0\t0.00", "t": "0\t0.00", "u": "3\t50.00", "n": "0\t0.00",
            "iupac_ambiguous": "0\t0.00", "invalid": "0\t0.00", "gc": "33.33",
            "gc_skew": "-1.000", "lowercase": "2\t33.33",
        },
        "odd_ties_mixed.fasta": {
            "type": "nucleotide_mixed_tu", "a": "2\t13.33", "c": "2\t13.33",
            "g": "2\t13.33", "t": "2\t13.33", "u": "2\t13.33", "n": "5\t33.33",
            "iupac_ambiguous": "0\t0.00", "invalid": "0\t0.00", "gc": "40.00",
            "gc_skew": "0.000", "lowercase": "0\t0.00",
        },
        "even_lengths.fasta": {
            "type": "nucleotide_t", "a": "4\t40.00", "c": "3\t30.00",
            "g": "2\t20.00", "t": "1\t10.00", "u": "0\t0.00", "n": "0\t0.00",
            "iupac_ambiguous": "0\t0.00", "invalid": "0\t0.00", "gc": "50.00",
            "gc_skew": "-0.200", "lowercase": "0\t0.00",
        },
        "nucleotide_categories.fasta": {
            "type": "nucleotide_mixed_tu", "a": "1\t5.88", "c": "1\t5.88",
            "g": "1\t5.88", "t": "1\t5.88", "u": "1\t5.88", "n": "1\t5.88",
            "iupac_ambiguous": "10\t58.82", "invalid": "1\t5.88", "gc": "40.00",
            "gc_skew": "0.000", "lowercase": "5\t29.41",
        },
        "ambiguity_only.fasta": {
            "type": "nucleotide", "a": "0\t0.00", "c": "0\t0.00", "g": "0\t0.00",
            "t": "0\t0.00", "u": "0\t0.00", "n": "2\t16.67",
            "iupac_ambiguous": "10\t83.33", "invalid": "0\t0.00", "gc": "n/a",
            "gc_skew": "n/a", "lowercase": "0\t0.00",
        },
        "gc_zero.fasta": {
            "type": "nucleotide_t", "a": "3\t50.00", "c": "0\t0.00",
            "g": "0\t0.00", "t": "3\t50.00", "u": "0\t0.00", "n": "0\t0.00",
            "iupac_ambiguous": "0\t0.00", "invalid": "0\t0.00", "gc": "0.00",
            "gc_skew": "n/a", "lowercase": "6\t100.00",
        },
    }
    for name, answer in composition_answers.items():
        assert_fields(name, composition_fields(design_values[name]), answer)

    exact_90 = composition_fields(design_values["exact_90_protein.fasta"])
    assert_fields("exact_90_protein.fasta", exact_90, {
        "type": "protein", "a_alanine": "9\t90.00", "l_leucine": "1\t10.00",
        "invalid": "0\t0.00", "lowercase": "0\t0.00",
    })
    protein = composition_fields(design_values["protein_complete.fasta"])
    assert protein["type"] == "protein"
    assert protein["a_alanine"] == "2\t6.67"
    assert protein["r_arginine"] == "2\t6.67"
    for key, _ in PROTEIN_FIELDS[2:]:
        assert protein[key] == "1\t3.33", (key, protein[key])
    assert_fields("protein_complete.fasta", protein, {
        "stop": "1\t3.33", "invalid": "1\t3.33", "lowercase": "2\t6.67",
    })
    assert_equal("self-test equality", "a\n", "a\n")
    for changed in ("a\nextra\n", "", "b\n"):
        try:
            assert_equal("self-test rejection", "a\n", changed)
        except SystemExit:
            pass
        else:
            raise AssertionError("strict comparison accepted changed output")


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    self_parser = commands.add_parser("self-test")
    self_parser.add_argument("--fixtures", required=True)
    self_parser.add_argument("--design-fixtures", required=True)
    verify_parser = commands.add_parser("verify")
    verify_parser.add_argument("--fasta", required=True)
    verify_parser.add_argument("--zfi", required=True)
    verify_parser.add_argument("--zfi-fasta")
    verify_parser.add_argument("--fai")
    verify_parser.add_argument("--fai-fasta")
    verify_parser.add_argument("--noodles")
    verify_parser.add_argument("--rustbio")
    verify_parser.add_argument("--no-dedup", action="store_true")
    expected_parser = commands.add_parser("expected")
    expected_parser.add_argument("fasta")
    expected_parser.add_argument("--no-dedup", action="store_true")
    check_parser = commands.add_parser("check")
    check_parser.add_argument("format", choices=("zfi", "fai"))
    check_parser.add_argument("fasta")
    check_parser.add_argument("expected")
    check_parser.add_argument("output")
    same_parser = commands.add_parser("same")
    same_parser.add_argument("left")
    same_parser.add_argument("right")
    parity_parser = commands.add_parser("parity")
    parity_parser.add_argument("fasta")
    parity_parser.add_argument("expected")
    parity_parser.add_argument("stats")
    parity_parser.add_argument("tool")
    parity_parser.add_argument("output")
    args = parser.parse_args()
    if args.command == "self-test":
        self_test(Path(args.fixtures), Path(args.design_fixtures))
    elif args.command == "verify":
        verify(args)
    else:
        legacy_command(args)


if __name__ == "__main__":
    main()
