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

#[path = "../../stats_peer.rs"]
mod stats_peer;

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

fn cmd_stats(fasta_path: &str) -> io::Result<()> {
    let file = File::open(fasta_path)?;
    let mut reader = Reader::new(BufReader::new(file));
    let mut seen: HashSet<String> = HashSet::new();
    let mut names: Vec<String> = Vec::new();
    let mut lengths: Vec<u64> = Vec::new();
    let mut counts = [0u64; 256];
    let mut lowercase = 0u64;

    for result in reader.records() {
        let record = result?;
        let name = String::from_utf8_lossy(record.name()).into_owned();
        let sequence = record.sequence();
        if sequence.is_empty() || !seen.insert(name.clone()) {
            continue;
        }
        names.push(name);
        lengths.push(sequence.len() as u64);
        for &byte in sequence.as_ref() {
            counts[byte as usize] += 1;
            if byte.is_ascii_lowercase() {
                lowercase += 1;
            }
        }
    }

    stats_peer::write_report(io::stdout().lock(), &names, &lengths, &counts, lowercase)
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
