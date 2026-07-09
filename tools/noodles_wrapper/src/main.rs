// noodles_wrapper.rs
// Extended CLI wrapper around noodles-fasta for benchmarking get/ operations.
// Subcommands: index, get, stats
//
// index  uses noodles' own fai indexer
// get    supports single/multi positional regions, full-name lookups, BED batch,
//        names-file batch, and --rc / --honor-strand reverse complement.
// stats  assembly/composition TSV for clean-FASTA bench parity only.
//        Re-parses the FASTA with noodles (no index, no messy/side-table support).
//        Field set mirrors z-fasta formulas so we have something to compare beyond
//        sequences/total_bases. Do not treat this as a messy-FASTA peer.
//
// Usage:
//   noodles_wrapper index <fasta>                    Build .fai index
//   noodles_wrapper get <fasta> <region>...          Extract region(s)
//   noodles_wrapper get <fasta> --bed <file>         Extract from BED file
//   noodles_wrapper get <fasta> --names <file>       Extract from names file
//   noodles_wrapper stats <fasta>                    Clean-FASTA stats (TSV)
//
// Flags for get:
//   --rc             Reverse complement all output
//   --honor-strand   In BED mode, RC only for reverse-strand entries

use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};

use noodles_core::Region;
use noodles_fasta::{fai, fs, io::{IndexedReader, Reader}};

fn usage() -> ! {
    eprintln!("Usage: noodles_wrapper <index|get|stats> <args...>");
    eprintln!("  index <fasta>                    Build .fai index");
    eprintln!("  get   <fasta> <region>...        Extract region(s), e.g. chr1:1000-2000");
    eprintln!("  get   <fasta> --bed <file>       Extract from BED file");
    eprintln!("  get   <fasta> --names <file>     Extract from names file");
    eprintln!("  stats <fasta>                    Clean-FASTA stats TSV (bench peer; not messy)");
    eprintln!("Flags: --rc, --honor-strand");
    std::process::exit(1);
}

// ════════════════════════════════════════════════════════════════════
//  Argument parsing
// ════════════════════════════════════════════════════════════════════

struct GetConfig<'a> {
    fasta: &'a str,
    mode: GetMode<'a>,
    rc: bool,
    honor_strand: bool,
}

enum GetMode<'a> {
    Regions(Vec<&'a str>),
    Bed(&'a str),
    Names(&'a str),
}

fn parse_get_args<'a>(args: &'a [String]) -> Result<GetConfig<'a>, String> {
    if args.is_empty() {
        return Err("get requires at least a FASTA path".into());
    }

    let fasta = &args[0];
    let mut rc = false;
    let mut honor_strand = false;
    let mut regions: Vec<&'a str> = Vec::new();
    let mut bed_file: Option<&'a str> = None;
    let mut names_file: Option<&'a str> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--rc" => rc = true,
            "--honor-strand" => honor_strand = true,
            "--bed" => {
                i += 1;
                bed_file = args.get(i).map(|s| s.as_str());
            }
            "--names" => {
                i += 1;
                names_file = args.get(i).map(|s| s.as_str());
            }
            reg => regions.push(reg),
        }
        i += 1;
    }

    let mode = if let Some(bed) = bed_file {
        GetMode::Bed(bed)
    } else if let Some(names) = names_file {
        GetMode::Names(names)
    } else if !regions.is_empty() {
        GetMode::Regions(regions)
    } else {
        return Err("no regions, --bed, or --names specified".into());
    };

    Ok(GetConfig { fasta, mode, rc, honor_strand })
}

// ════════════════════════════════════════════════════════════════════
//  Helpers
// ════════════════════════════════════════════════════════════════════

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
        b'U' => b'A', b'W' => b'W', b'S' => b'S', b'M' => b'K',
        b'K' => b'M', b'R' => b'Y', b'Y' => b'R', b'B' => b'V',
        b'D' => b'H', b'H' => b'D', b'V' => b'B', b'N' => b'N',
        b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a',
        b'u' => b'a', b'w' => b'w', b's' => b's', b'm' => b'k',
        b'k' => b'm', b'r' => b'y', b'y' => b'r', b'b' => b'v',
        b'd' => b'h', b'h' => b'd', b'v' => b'b', b'n' => b'n',
        _ => b'N',
    }).collect()
}

fn open_reader(fasta_path: &str) -> io::Result<IndexedReader<BufReader<File>>> {
    let fai_path = format!("{fasta_path}.fai");
    let index = fai::fs::read(&fai_path)?;
    let file = File::open(fasta_path)?;
    Ok(IndexedReader::new(BufReader::new(file), index))
}

fn write_fasta(
    out: &mut BufWriter<io::StdoutLock<'_>>,
    header: &str,
    seq: &[u8],
    do_rc: bool,
) -> io::Result<()> {
    writeln!(out, ">{header}")?;
    let data = if do_rc { revcomp(seq) } else { seq.to_vec() };
    for chunk in data.chunks(60) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

// ════════════════════════════════════════════════════════════════════
//  cmd_index
// ════════════════════════════════════════════════════════════════════

fn cmd_index(fasta: &str) -> io::Result<()> {
    let index = fs::index(fasta)?;
    let fai_path = format!("{fasta}.fai");
    fai::fs::write(&fai_path, &index)
}

// ════════════════════════════════════════════════════════════════════
//  cmd_get: single/multi-region (positional or name-only)
// ════════════════════════════════════════════════════════════════════

fn cmd_get_regions(config: &GetConfig) -> io::Result<()> {
    let regions = match &config.mode {
        GetMode::Regions(r) => r,
        _ => unreachable!(),
    };

    let mut reader = open_reader(config.fasta)?;
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    for region_str in regions {
        let region: Region = region_str.parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("Invalid region: {e}"))
        })?;
        let record = reader.query(&region)?;
        write_fasta(&mut out, region_str, record.sequence().as_ref(), config.rc)?;
    }

    out.flush()?;
    Ok(())
}

// ════════════════════════════════════════════════════════════════════
//  cmd_get_bed: BED batch extraction
// ════════════════════════════════════════════════════════════════════

fn cmd_get_bed(config: &GetConfig) -> io::Result<()> {
    let bed_path = match &config.mode {
        GetMode::Bed(p) => p,
        _ => unreachable!(),
    };

    let mut reader = open_reader(config.fasta)?;
    let bed_file = File::open(bed_path)?;
    let bed_reader = BufReader::new(bed_file);
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    for line in bed_reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chrom = fields[0];
        let start: u64 = fields[1].parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid BED start: {e}"))
        })?;
        let end: u64 = fields[2].parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid BED end: {e}"))
        })?;

        // BED: 0-based start, half-open end -> noodles: 1-based inclusive
        let region_str = format!("{}:{}-{}", chrom, start + 1, end);
        let region: Region = region_str.parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("Invalid region: {e}"))
        })?;

        let record = reader.query(&region)?;

        // BED name (col 4) as header if present and non-empty
        let name = fields.get(3).filter(|s| !s.is_empty() && **s != ".");
        let header = match name {
            Some(n) => n.to_string(),
            None => region_str,
        };

        // Strand-aware RC: col 6 (BED6+)
        let strand = fields.get(5).copied();
        let rc_this = config.rc && (!config.honor_strand || strand == Some("-"));

        write_fasta(&mut out, &header, record.sequence().as_ref(), rc_this)?;
    }

    out.flush()?;
    Ok(())
}

// ════════════════════════════════════════════════════════════════════
//  cmd_get_names: names-file extraction
// ════════════════════════════════════════════════════════════════════

fn cmd_get_names(config: &GetConfig) -> io::Result<()> {
    let names_path = match &config.mode {
        GetMode::Names(p) => p,
        _ => unreachable!(),
    };

    let mut reader = open_reader(config.fasta)?;
    let names_file = File::open(names_path)?;
    let names_reader = BufReader::new(names_file);
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    for line in names_reader.lines() {
        let name = line?;
        let trimmed = name.trim().to_string();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let region: Region = trimmed.parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidInput, format!("Invalid region: {e}"))
        })?;
        let record = reader.query(&region)?;
        write_fasta(&mut out, &trimmed, record.sequence().as_ref(), config.rc)?;
    }

    out.flush()?;
    Ok(())
}

// ════════════════════════════════════════════════════════════════════
//  cmd_stats
//  Clean-FASTA comparison peer only. Emits the same metric formulas as
//  z-fasta so benches/verify have fields beyond sequences/total_bases.
//  Does not strip messy whitespace or use side tables.
// ════════════════════════════════════════════════════════════════════

const IUPAC_NUC: &[u8] = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn";
const AA_CODES: &[u8] = b"ARNDCEQGHILKMFPSTWYV";
const AA_NAMES: &[&str] = &[
    "Alanine",
    "Arginine",
    "Asparagine",
    "Aspartate",
    "Cysteine",
    "Glutamate",
    "Glutamine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Leucine",
    "Lysine",
    "Methionine",
    "Phenylalanine",
    "Proline",
    "Serine",
    "Threonine",
    "Tryptophan",
    "Tyrosine",
    "Valine",
];

fn detect_type(counts: &[u64; 256], total: u64) -> &'static str {
    if total == 0 {
        return "nucleotide";
    }
    let mut nuc = 0u64;
    for &b in IUPAC_NUC {
        nuc += counts[b as usize];
    }
    if nuc * 10 > total * 9 {
        "nucleotide"
    } else {
        "protein"
    }
}

fn cmd_stats(fasta_path: &str) -> io::Result<()> {
    let file = File::open(fasta_path)?;
    let mut reader = Reader::new(BufReader::new(file));
    let mut seen: HashSet<String> = HashSet::new();
    let mut names: Vec<String> = Vec::new();
    let mut lengths: Vec<u64> = Vec::new();
    let mut counts = [0u64; 256];
    let mut lowercase = 0u64;
    let mut comp_total = 0u64;

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let seq = record.sequence();
        if seq.is_empty() || !seen.insert(name.clone()) {
            continue;
        }
        names.push(name);
        lengths.push(seq.len() as u64);
        for &b in seq.as_ref() {
            counts[b as usize] += 1;
            comp_total += 1;
            if b.is_ascii_lowercase() {
                lowercase += 1;
            }
        }
    }

    let num_seqs = lengths.len() as u64;
    let total_bases: u64 = lengths.iter().sum();
    let mut sorted = lengths.clone();
    sorted.sort_by(|a, b| b.cmp(a));

    let mean = if num_seqs > 0 { total_bases / num_seqs } else { 0 };
    let median = if num_seqs == 0 {
        0
    } else if num_seqs % 2 == 1 {
        sorted[(num_seqs / 2) as usize]
    } else {
        (sorted[(num_seqs / 2 - 1) as usize] + sorted[(num_seqs / 2) as usize]) / 2
    };

    let threshold_50 = (total_bases + 1) / 2;
    let threshold_90 = (total_bases * 9 + 9) / 10;
    let mut bases_seen = 0u64;
    let mut au_sum: u128 = 0;
    let mut n50 = 0u64;
    let mut l50 = 0u64;
    let mut n90 = 0u64;
    let mut l90 = 0u64;
    let mut found_n50 = false;
    let mut found_n90 = false;
    for (i, &len) in sorted.iter().enumerate() {
        bases_seen += len;
        au_sum += (len as u128) * (len as u128);
        if !found_n50 && bases_seen >= threshold_50 {
            n50 = len;
            l50 = (i + 1) as u64;
            found_n50 = true;
        }
        if !found_n90 && bases_seen >= threshold_90 {
            n90 = len;
            l90 = (i + 1) as u64;
            found_n90 = true;
        }
    }
    let au = if total_bases > 0 {
        (au_sum / total_bases as u128) as u64
    } else {
        0
    };

    let (shortest_idx, _) = lengths
        .iter()
        .enumerate()
        .min_by_key(|(_, &l)| l)
        .unwrap_or((0, &0));
    let (longest_idx, _) = lengths
        .iter()
        .enumerate()
        .max_by_key(|(_, &l)| l)
        .unwrap_or((0, &0));
    let shortest_len = lengths.get(shortest_idx).copied().unwrap_or(0);
    let longest_len = lengths.get(longest_idx).copied().unwrap_or(0);
    let shortest_name = names.get(shortest_idx).map(|s| s.as_str()).unwrap_or("");
    let longest_name = names.get(longest_idx).map(|s| s.as_str()).unwrap_or("");

    let seq_type = detect_type(&counts, comp_total);
    let soft_pct = if comp_total > 0 {
        lowercase as f64 / comp_total as f64 * 100.0
    } else {
        0.0
    };

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    writeln!(out, "sequences\t{num_seqs}")?;
    writeln!(out, "total_bases\t{total_bases}")?;
    writeln!(out, "shortest_len\t{shortest_len}")?;
    writeln!(out, "shortest_name\t{shortest_name}")?;
    writeln!(out, "longest_len\t{longest_len}")?;
    writeln!(out, "longest_name\t{longest_name}")?;
    writeln!(out, "mean\t{mean}")?;
    writeln!(out, "median\t{median}")?;
    writeln!(out, "n50\t{n50}")?;
    writeln!(out, "l50\t{l50}")?;
    writeln!(out, "n90\t{n90}")?;
    writeln!(out, "l90\t{l90}")?;
    writeln!(out, "au\t{au}")?;
    writeln!(out, "type\t{seq_type}")?;
    writeln!(out, "soft_pct\t{soft_pct:.4}")?;

    if seq_type == "nucleotide" {
        let a = counts[b'A' as usize] + counts[b'a' as usize];
        let c = counts[b'C' as usize] + counts[b'c' as usize];
        let g = counts[b'G' as usize] + counts[b'g' as usize];
        let t = counts[b'T' as usize] + counts[b't' as usize];
        let n = counts[b'N' as usize] + counts[b'n' as usize];
        let acgt = a + c + g + t;
        let other = comp_total.saturating_sub(a + c + g + t + n);
        let f_total = comp_total as f64;
        let pct = |x: u64| if f_total > 0.0 { x as f64 / f_total * 100.0 } else { 0.0 };
        writeln!(out, "a_pct\t{:.4}", pct(a))?;
        writeln!(out, "c_pct\t{:.4}", pct(c))?;
        writeln!(out, "g_pct\t{:.4}", pct(g))?;
        writeln!(out, "t_pct\t{:.4}", pct(t))?;
        writeln!(out, "n_pct\t{:.4}", pct(n))?;
        writeln!(out, "other_pct\t{:.4}", pct(other))?;
        let gc_pct = if acgt > 0 {
            (g + c) as f64 / acgt as f64 * 100.0
        } else {
            0.0
        };
        writeln!(out, "gc_pct\t{gc_pct:.4}")?;
        let gc_sum = g + c;
        if gc_sum > 0 {
            let skew = (g as f64 - c as f64) / gc_sum as f64;
            writeln!(out, "gc_skew\t{skew:.6}")?;
        }
        writeln!(out, "n_content\t{n}")?;
    } else {
        let mut aa: Vec<(u8, u64, usize)> = AA_CODES
            .iter()
            .enumerate()
            .map(|(i, &code)| {
                let lower = code.to_ascii_lowercase();
                let cnt = counts[code as usize] + counts[lower as usize];
                (code, cnt, i)
            })
            .collect();
        aa.sort_by(|a, b| b.1.cmp(&a.1).then(a.2.cmp(&b.2)));
        for (rank, (code, cnt, idx)) in aa.iter().take(3).enumerate() {
            let pct = if comp_total > 0 {
                *cnt as f64 / comp_total as f64 * 100.0
            } else {
                0.0
            };
            writeln!(
                out,
                "top_aa_{}\t{}:{:.4}:{}",
                rank + 1,
                *code as char,
                pct,
                AA_NAMES[*idx]
            )?;
        }
    }

    Ok(())
}

// ════════════════════════════════════════════════════════════════════
//  Main dispatch
// ════════════════════════════════════════════════════════════════════

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        usage();
    }

    let result = match args[1].as_str() {
        "index" => {
            if args.len() < 3 {
                eprintln!("Error: index requires a FASTA path");
                std::process::exit(1);
            }
            cmd_index(&args[2])
        }
        "get" => {
            if args.len() < 4 {
                eprintln!("Error: get requires <fasta> <region|--bed|--names>");
                std::process::exit(1);
            }
            let config = parse_get_args(&args[2..]).unwrap_or_else(|e| {
                eprintln!("Error: {e}");
                std::process::exit(1);
            });
            match config.mode {
                GetMode::Regions(_) => cmd_get_regions(&config),
                GetMode::Bed(_) => cmd_get_bed(&config),
                GetMode::Names(_) => cmd_get_names(&config),
            }
        }
        "stats" => {
            if args.len() < 3 {
                eprintln!("Error: stats requires a FASTA path");
                std::process::exit(1);
            }
            cmd_stats(&args[2])
        }
        _ => {
            eprintln!("Unknown command: {}", args[1]);
            usage();
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
