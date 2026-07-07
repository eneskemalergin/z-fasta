#!/usr/bin/env python3
"""Messy-FASTA oracles for bench/get/verify.sh only. Extend in place; no sibling scripts."""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path

RC = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def load_seqs(fasta: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    for line in Path(fasta).read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append("".join(c for c in line if not c.isspace()))
    return {name: "".join(parts) for name, parts in seqs.items()}


def slice_region(seqs: dict[str, str], name: str, start: int, end: int) -> str:
    return seqs[name][start - 1 : end]


def write_region(path: str, name: str, start: int, end: int, seq: str) -> None:
    Path(path).write_text(f">{name}:{start}-{end}\n{seq}\n")


def iter_bed_rows(bed: str):
    for line in Path(bed).read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            continue
        yield line.split("\t")


def cmd_region(argv: list[str], transform: Callable[[str], str] | None = None) -> None:
    fasta, name, start, end, out = argv[2], argv[3], int(argv[4]), int(argv[5]), argv[6]
    seq = slice_region(load_seqs(fasta), name, start, end)
    if transform is not None:
        seq = transform(seq)
    write_region(out, name, start, end, seq)


def cmd_bed(argv: list[str], stranded_rc: bool = False) -> None:
    fasta, bed, out = argv[2], argv[3], argv[4]
    seqs = load_seqs(fasta)
    parts: list[str] = []
    for fields in iter_bed_rows(bed):
        chrom, s0, e0 = fields[0], int(fields[1]), int(fields[2])
        seq = seqs[chrom][s0:e0]
        if stranded_rc and (fields[5] if len(fields) >= 6 else ".") != "-":
            seq = seq.translate(RC)[::-1]
        parts.append(f">{chrom}:{s0 + 1}-{e0}\n{seq}\n")
    Path(out).write_text("".join(parts))


COMMANDS = {
    "region": lambda argv: cmd_region(argv),
    "rc": lambda argv: cmd_region(argv, lambda s: s.translate(RC)[::-1]),
    "rev": lambda argv: cmd_region(argv, lambda s: s[::-1]),
    "comp": lambda argv: cmd_region(argv, lambda s: s.translate(RC)),
    "bed": lambda argv: cmd_bed(argv),
    "bed-stranded-rc": lambda argv: cmd_bed(argv, stranded_rc=True),
}


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print("usage: oracle.py <region|rc|rev|comp|bed|bed-stranded-rc> ...", file=sys.stderr)
        sys.exit(2)
    COMMANDS[sys.argv[1]](sys.argv)


if __name__ == "__main__":
    main()
