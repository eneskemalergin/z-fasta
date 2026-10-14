use std::io::{self, Write};

const IUPAC: &[u8] = b"ACGTURYSWKMBDHVN";
const AMBIGUOUS: &[u8] = b"RYSWKMBDHV";
const PROTEIN_FIELDS: &[(&str, u8)] = &[
    ("a_alanine", b'A'),
    ("r_arginine", b'R'),
    ("n_asparagine", b'N'),
    ("d_aspartate", b'D'),
    ("c_cysteine", b'C'),
    ("e_glutamate", b'E'),
    ("q_glutamine", b'Q'),
    ("g_glycine", b'G'),
    ("h_histidine", b'H'),
    ("i_isoleucine", b'I'),
    ("l_leucine", b'L'),
    ("k_lysine", b'K'),
    ("m_methionine", b'M'),
    ("f_phenylalanine", b'F'),
    ("p_proline", b'P'),
    ("s_serine", b'S'),
    ("t_threonine", b'T'),
    ("w_tryptophan", b'W'),
    ("y_tyrosine", b'Y'),
    ("v_valine", b'V'),
    ("b_asx", b'B'),
    ("z_glx", b'Z'),
    ("j_xle", b'J'),
    ("x_unknown", b'X'),
    ("u_selenocysteine", b'U'),
    ("o_pyrrolysine", b'O'),
];

fn write_fixed_unsigned<W: Write>(
    out: &mut W,
    numerator: u128,
    denominator: u128,
    places: u32,
) -> io::Result<()> {
    let scale = 10u128.pow(places);
    let mut whole = numerator / denominator;
    let mut fraction = ((numerator % denominator) * scale + denominator / 2) / denominator;
    if fraction == scale {
        whole += 1;
        fraction = 0;
    }
    write!(out, "{whole}.{fraction:0width$}", width = places as usize)
}

fn write_fixed_signed<W: Write>(
    out: &mut W,
    numerator: i128,
    denominator: u128,
    places: u32,
) -> io::Result<()> {
    let magnitude = numerator.unsigned_abs();
    let scale = 10u128.pow(places);
    let rounded_fraction = ((magnitude % denominator) * scale + denominator / 2) / denominator;
    if numerator < 0 && (magnitude / denominator != 0 || rounded_fraction != 0) {
        out.write_all(b"-")?;
    }
    write_fixed_unsigned(out, magnitude, denominator, places)
}

fn letter(counts: &[u64; 256], code: u8) -> u64 {
    counts[code as usize] + counts[code.to_ascii_lowercase() as usize]
}

fn median(values: &[u64]) -> u64 {
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        values[middle - 1] + (values[middle] - values[middle - 1]) / 2
    }
}

fn count_percent<W: Write>(out: &mut W, key: &str, count: u64, total: u64) -> io::Result<()> {
    write!(out, "{key}\t{count}\t")?;
    write_fixed_unsigned(out, count as u128 * 100, total as u128, 2)?;
    out.write_all(b"\n")
}

pub fn write_report<W: Write>(
    mut out: W,
    names: &[String],
    lengths: &[u64],
    counts: &[u64; 256],
    lowercase: u64,
) -> io::Result<()> {
    if lengths.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no retained records",
        ));
    }
    let total = lengths
        .iter()
        .try_fold(0u64, |sum, value| sum.checked_add(*value))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "sequence length overflow"))?;
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable();
    let half = sorted.len() / 2;
    let q1 = if sorted.len() == 1 {
        sorted[0]
    } else {
        median(&sorted[..half])
    };
    let q3 = if sorted.len() == 1 {
        sorted[0]
    } else {
        median(&sorted[sorted.len() - half..])
    };
    let mut shortest = 0usize;
    let mut longest = 0usize;
    for index in 1..lengths.len() {
        if lengths[index] < lengths[shortest] {
            shortest = index;
        }
        if lengths[index] > lengths[longest] {
            longest = index;
        }
    }
    let threshold_50 = (total as u128 * 50).div_ceil(100);
    let threshold_90 = (total as u128 * 90).div_ceil(100);
    let mut seen = 0u128;
    let mut square_sum = 0u128;
    let (mut n50, mut l50, mut n90, mut l90) = (0u64, 0usize, 0u64, 0usize);
    for (index, length) in sorted.iter().rev().enumerate() {
        seen += *length as u128;
        square_sum += *length as u128 * *length as u128;
        if l50 == 0 && seen >= threshold_50 {
            n50 = *length;
            l50 = index + 1;
        }
        if l90 == 0 && seen >= threshold_90 {
            n90 = *length;
            l90 = index + 1;
        }
    }

    writeln!(out, "indexed_records\t{}", lengths.len())?;
    writeln!(out, "total_symbols\t{total}")?;
    writeln!(out, "shortest_length\t{}", lengths[shortest])?;
    writeln!(out, "shortest_name\t{}", names[shortest])?;
    writeln!(out, "longest_length\t{}", lengths[longest])?;
    writeln!(out, "longest_name\t{}", names[longest])?;
    writeln!(out, "mean\t{}", total / lengths.len() as u64)?;
    writeln!(out, "q1\t{q1}")?;
    writeln!(out, "median\t{}", median(&sorted))?;
    writeln!(out, "q3\t{q3}")?;
    writeln!(out, "range\t{}", lengths[longest] - lengths[shortest])?;
    writeln!(out, "n50\t{n50}")?;
    writeln!(out, "l50\t{l50}")?;
    writeln!(out, "n90\t{n90}")?;
    writeln!(out, "l90\t{l90}")?;
    out.write_all(b"aun\t")?;
    write_fixed_unsigned(&mut out, square_sum, total as u128, 2)?;
    out.write_all(b"\n")?;

    let nucleotide = IUPAC.iter().map(|code| letter(counts, *code)).sum::<u64>() as u128 * 10
        > total as u128 * 9;
    if nucleotide {
        let a = letter(counts, b'A');
        let c = letter(counts, b'C');
        let g = letter(counts, b'G');
        let t = letter(counts, b'T');
        let u = letter(counts, b'U');
        let n = letter(counts, b'N');
        let ambiguous = AMBIGUOUS
            .iter()
            .map(|code| letter(counts, *code))
            .sum::<u64>();
        let invalid = total - a - c - g - t - u - n - ambiguous;
        let type_name = if t > 0 && u > 0 {
            "nucleotide_mixed_tu"
        } else if t > 0 {
            "nucleotide_t"
        } else if u > 0 {
            "nucleotide_u"
        } else {
            "nucleotide"
        };
        writeln!(out, "type\t{type_name}")?;
        writeln!(out, "percent_denominator\ttotal_symbols")?;
        for (key, count) in [
            ("a", a),
            ("c", c),
            ("g", g),
            ("t", t),
            ("u", u),
            ("n", n),
            ("iupac_ambiguous", ambiguous),
            ("invalid", invalid),
        ] {
            count_percent(&mut out, key, count, total)?;
        }
        let canonical = a + c + g + t + u;
        if canonical == 0 {
            writeln!(out, "gc\tn/a")?;
        } else {
            out.write_all(b"gc\t")?;
            write_fixed_unsigned(&mut out, (g + c) as u128 * 100, canonical as u128, 2)?;
            out.write_all(b"\n")?;
        }
        let gc_total = g + c;
        if gc_total == 0 {
            writeln!(out, "gc_skew\tn/a")?;
        } else {
            out.write_all(b"gc_skew\t")?;
            write_fixed_signed(&mut out, g as i128 - c as i128, gc_total as u128, 3)?;
            out.write_all(b"\n")?;
        }
    } else {
        writeln!(out, "type\tprotein")?;
        writeln!(out, "percent_denominator\ttotal_symbols")?;
        let mut assigned = 0u64;
        for (key, code) in PROTEIN_FIELDS {
            let count = letter(counts, *code);
            assigned += count;
            count_percent(&mut out, key, count, total)?;
        }
        let stop = counts[b'*' as usize];
        count_percent(&mut out, "stop", stop, total)?;
        count_percent(&mut out, "invalid", total - assigned - stop, total)?;
    }
    count_percent(&mut out, "lowercase", lowercase, total)?;
    out.flush()
}
