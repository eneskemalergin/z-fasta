// rustbio_wrapper.rs
// Extended CLI wrapper around rust-bio for benchmarking get/ operations.
// Subcommands: index, get, stats
//
// index  uses a strict FAI-style scanner (uniform line widths, no blank lines).
// get    supports single/multi positional regions, full-name lookups, BED batch,
//        names-file batch, and --rc / --honor-strand reverse complement.
// stats  assembly/composition TSV for clean-FASTA bench parity only.
//        Re-parses the FASTA with rust-bio (no index, no messy/side-table support).
//        Field set mirrors z-fasta formulas so we have something to compare beyond
//        sequences/total_bases. Do not treat this as a messy-FASTA peer.
//
// Usage:
//   rustbio_wrapper index <fasta>                    Build .fai index
//   rustbio_wrapper get <fasta> <region>...          Extract region(s)
//   rustbio_wrapper get <fasta> --bed <file>         Extract from BED file
//   rustbio_wrapper get <fasta> --names <file>       Extract from names file
//   rustbio_wrapper stats <fasta>                    Clean-FASTA stats (TSV)
//
// Flags for get:
//   --rc             Reverse complement all output
//   --honor-strand   In BED mode with --rc, RC only for reverse-strand entries

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};
use std::path::PathBuf;

use bio::io::fasta;

#[path = "../../stats_peer.rs"]
mod stats_peer;

const PEER_VERSION: &str = "4.0.1";

fn print_version() {
    println!("rust-bio {PEER_VERSION} (custom indexer wrapper)");
}

fn usage() -> ! {
    eprintln!("Usage: rustbio_wrapper <index|get|stats> <args...>");
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
    /// One or more region strings (positional or name-only).
    Regions(Vec<&'a str>),
    /// BED file path.
    Bed(&'a str),
    /// Names file path.
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

/// IUPAC reverse complement (also handles U/u uracil).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', b'U' => b'A',
        b'W' => b'W', b'S' => b'S', b'M' => b'K', b'K' => b'M',
        b'R' => b'Y', b'Y' => b'R', b'B' => b'V', b'D' => b'H',
        b'H' => b'D', b'V' => b'B', b'N' => b'N',
        b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a', b'u' => b'a',
        b'w' => b'w', b's' => b's', b'm' => b'k', b'k' => b'm',
        b'r' => b'y', b'y' => b'r', b'b' => b'v', b'd' => b'h',
        b'h' => b'd', b'v' => b'b', b'n' => b'n',
        _ => b'N',
    }).collect()
}

/// Parse a region string into chrom and optional (1-based inclusive) coordinates.
fn parse_region(s: &str) -> (&str, Option<(u64, u64)>) {
    if let Some(colon) = s.find(':') {
        let chrom = &s[..colon];
        let coords = &s[colon + 1..];
        if let Some(dash) = coords.find('-') {
            let start: u64 = coords[..dash].parse().unwrap_or(0);
            let end: u64 = coords[dash + 1..].parse().unwrap_or(0);
            (chrom, Some((start, end)))
        } else {
            (s, None)
        }
    } else {
        (s, None)
    }
}

/// Read .fai to build sequence name -> length map.
fn read_fai_lengths(fasta_path: &str) -> io::Result<HashMap<String, u64>> {
    let fai_path = format!("{fasta_path}.fai");
    let file = File::open(&fai_path)?;
    let reader = BufReader::new(file);
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            if let Ok(len) = parts[1].parse::<u64>() {
                map.insert(parts[0].to_string(), len);
            }
        }
    }
    Ok(map)
}

/// Open an indexed FASTA reader.
fn open_reader(fasta_path: &str) -> io::Result<fasta::IndexedReader<File>> {
    let path = PathBuf::from(fasta_path);
    fasta::IndexedReader::from_file(&path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::Other,
            format!("failed to open indexed FASTA: {e}"),
        )
    })
}

/// Fetch a region (0-based) into buf; clears buf first.
fn fetch_into(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
    buf: &mut Vec<u8>,
) -> io::Result<()> {
    reader.fetch(chrom, start, end)?;
    buf.clear();
    reader.read(buf)?;
    Ok(())
}

/// Write a FASTA record with 60-char wrapping, optionally RC'd.
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
//  cmd_index  (unchanged strict-FAI scanner)
// ════════════════════════════════════════════════════════════════════

fn fail_index(writer: &mut BufWriter<File>, err: io::Error) -> io::Error {
    let _ = writer.get_mut().set_len(0);
    err
}

fn cmd_index(fasta: &str) -> io::Result<()> {
    let src = File::open(fasta)?;
    let reader = BufReader::new(src);
    let fai_path = format!("{fasta}.fai");
    let mut writer = BufWriter::new(File::create(&fai_path)?);

    let mut file_offset: u64 = 0;
    let mut seq_len: u64 = 0;
    let mut line_bases: u64 = 0;
    let mut line_len: u64 = 0;
    let mut data_offset: u64 = 0;
    let mut current_name: Option<String> = None;
    let mut is_first_seq = true;
    let mut had_short_line = false;
    let mut record_count: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        let raw_len = line.len() as u64;
        let line_bytes = raw_len + 1;

        if line.is_empty() {
            return Err(fail_index(
                &mut writer,
                io::Error::new(io::ErrorKind::InvalidData, "blank line in FASTA input"),
            ));
        }

        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                if seq_len == 0 {
                    return Err(fail_index(
                        &mut writer,
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("sequence '{name}' has no bases"),
                        ),
                    ));
                }
                writeln!(
                    writer,
                    "{name}\t{seq_len}\t{data_offset}\t{line_bases}\t{line_len}"
                )?;
                record_count += 1;
            }

            let header = &line[1..];
            let name = header.split_whitespace().next().unwrap_or(header).to_string();
            if name.is_empty() {
                return Err(fail_index(
                    &mut writer,
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "empty sequence name in FASTA header",
                    ),
                ));
            }

            current_name = Some(name);
            file_offset += line_bytes;
            data_offset = file_offset;
            seq_len = 0;
            line_bases = 0;
            line_len = 0;
            is_first_seq = true;
            had_short_line = false;
        } else {
            let name = current_name.as_ref().ok_or_else(|| {
                fail_index(
                    &mut writer,
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "sequence data before first FASTA header",
                    ),
                )
            })?;

            if is_first_seq {
                line_bases = raw_len;
                line_len = line_bytes;
                is_first_seq = false;
            } else if raw_len == line_bases {
                if had_short_line {
                    return Err(fail_index(
                        &mut writer,
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("non-uniform line length in sequence '{name}'"),
                        ),
                    ));
                }
            } else if raw_len > line_bases {
                return Err(fail_index(
                    &mut writer,
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("non-uniform line length in sequence '{name}'"),
                    ),
                ));
            } else if had_short_line {
                return Err(fail_index(
                    &mut writer,
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("non-uniform line length in sequence '{name}'"),
                    ),
                ));
            } else {
                had_short_line = true;
            }

            seq_len += raw_len;
            file_offset += line_bytes;
        }
    }

    if let Some(name) = current_name {
        if seq_len == 0 {
            return Err(fail_index(
                &mut writer,
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("sequence '{name}' has no bases"),
                ),
            ));
        }
        writeln!(
            writer,
            "{name}\t{seq_len}\t{data_offset}\t{line_bases}\t{line_len}"
        )?;
        record_count += 1;
    }

    if record_count == 0 {
        return Err(fail_index(
            &mut writer,
            io::Error::new(io::ErrorKind::InvalidData, "no FASTA records found"),
        ));
    }

    Ok(())
}

// ════════════════════════════════════════════════════════════════════
//  cmd_get: single/multi-region (positional or name-only)
// ════════════════════════════════════════════════════════════════════

fn cmd_get_regions(config: &GetConfig) -> io::Result<()> {
    let regions = match &config.mode {
        GetMode::Regions(r) => r,
        _ => unreachable!(),
    };

    let lengths = read_fai_lengths(config.fasta)?;
    let mut reader = open_reader(config.fasta)?;
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let mut seq_buf = Vec::new();

    for region_str in regions {
        let (chrom, coords) = parse_region(region_str);
        match coords {
            Some((start, end)) => {
                // samtools-style 1-based inclusive -> 0-based half-open
                fetch_into(&mut reader, chrom, start - 1, end, &mut seq_buf)?;
            }
            None => {
                let len = lengths.get(chrom).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("sequence '{chrom}' not found in index"),
                    )
                })?;
                fetch_into(&mut reader, chrom, 0, *len, &mut seq_buf)?;
            }
        }
        write_fasta(&mut out, region_str, &seq_buf, config.rc)?;
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
    let mut seq_buf = Vec::new();

    for line in bed_reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid BED start: {e}"))
        })?;
        let end: u64 = fields[2].parse().map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid BED end: {e}"))
        })?;

        let name = fields.get(3).filter(|s| !s.is_empty() && **s != ".");
        let header = match name {
            Some(n) => n.to_string(),
            None => format!("{}:{}-{}", chrom, start + 1, end),
        };

        let strand = fields.get(5).copied();
        let rc_this = config.rc && (!config.honor_strand || strand == Some("-"));

        fetch_into(&mut reader, &chrom, start, end, &mut seq_buf)?;
        write_fasta(&mut out, &header, &seq_buf, rc_this)?;
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

    let lengths = read_fai_lengths(config.fasta)?;
    let mut reader = open_reader(config.fasta)?;
    let names_file = File::open(names_path)?;
    let names_reader = BufReader::new(names_file);
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let mut seq_buf = Vec::new();

    for line in names_reader.lines() {
        let name = line?;
        let trimmed = name.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let len = lengths.get(trimmed).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("sequence '{trimmed}' not found in index"),
            )
        })?;
        fetch_into(&mut reader, trimmed, 0, *len, &mut seq_buf)?;
        write_fasta(&mut out, trimmed, &seq_buf, config.rc)?;
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
    let source = File::open(fasta_path)?;
    let reader = fasta::Reader::new(BufReader::new(source));
    let mut seen: HashSet<String> = HashSet::new();
    let mut names: Vec<String> = Vec::new();
    let mut lengths: Vec<u64> = Vec::new();
    let mut counts = [0u64; 256];
    let mut lowercase = 0u64;

    for result in reader.records() {
        let record = result?;
        let name = record.id().to_string();
        let sequence = record.seq();
        if sequence.is_empty() || !seen.insert(name.clone()) {
            continue;
        }
        names.push(name);
        lengths.push(sequence.len() as u64);
        for &byte in sequence {
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
    if args.get(1).map(String::as_str) == Some("--version") {
        print_version();
        return;
    }
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
