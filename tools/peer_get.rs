pub struct GetConfig<'a> {
    pub fasta: &'a str,
    pub mode: GetMode<'a>,
    pub rc: bool,
    pub honor_strand: bool,
}

pub enum GetMode<'a> {
    Regions(Vec<&'a str>),
    Bed(&'a str),
    Names(&'a str),
}

pub fn parse_get_args<'a>(args: &'a [String]) -> Result<GetConfig<'a>, String> {
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

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
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
