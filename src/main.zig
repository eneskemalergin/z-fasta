//! CLI entry and subcommand dispatcher for `index`, `get`, `stats`, and `validate`.

const std = @import("std");
const builtin = @import("builtin");

pub const index_format = @import("index_format.zig");
pub const indexer = @import("indexer.zig");
pub const getter = @import("getter.zig");
pub const bed_parser = @import("bed_parser.zig");
pub const stats = @import("stats.zig");
pub const validator = @import("validator.zig");

const printErrorAndExit = index_format.printErrorAndExit;

const VERSION = "0.3.4";
const STRAND_AWARE_FLAG = "--strand-aware";
const STRAND_AWARE_ALIAS = "--honor-strand";
const RC_FLAG = "--rc";
const COMPLEMENT_ONLY_FLAG = "--complement-only";
const REVERSE_ONLY_FLAG = "--reverse-only";
const ANNOTATE_RC_FLAG = "--annotate-rc";

const USAGE =
    \\usage: z-fasta <command> [options]
    \\
    \\Commands:
    \\  index    Build a FASTA index (.zfi binary or .fai text)
    \\  get      Extract a sequence or region from a FASTA file
    \\  stats    Show statistics for a FASTA file
    \\  validate Check FASTA structure, alphabets, headers, and fix safe format issues
    \\
    \\General options:
    \\  --help       Show this help message
    \\  --version    Print version
    \\
    \\Index options:
    \\  --emit-fai   Output FAI to stdout when every record has fixed line geometry;
    \\               otherwise fails and directs you to default .zfi indexing
    \\  --no-dedup   Keep duplicate sequence names in the index (default: first wins
    \\               at index time). get resolves duplicate names to the first record.
    \\
    \\Get usage:
    \\  z-fasta get <file.fasta> [options] <region> [region ...]
    \\  z-fasta get <file.fasta> [options] --names file.txt|-
    \\  z-fasta get <file.fasta> [options] --bed file.bed|- [--strand-aware]
    \\  Choose exactly one request source: positional regions, --names, or --bed.
    \\  Region formats: NAME, NAME:START-END, NAME:START-
    \\  --strand-aware  With --bed, respect column 6; '-' emits reverse-complement output
    \\                  (alias: --honor-strand)
    \\  --rc            Reverse-complement extracted sequence output
    \\  --complement-only  Complement extracted sequence output without reversing
    \\  --reverse-only  Reverse extracted sequence output without complementing
    \\  --annotate-rc   Annotate the final output orientation when it is transformed
    \\  Positional regions: max 1024 per invocation. Names and BED input stream.
    \\
    \\Validate usage:
    \\  z-fasta validate [options] <file.fasta>
    \\  --strict                 Treat warnings as errors
    \\  --json                   Emit JSON Lines instead of text
    \\  --summary                With --json, emit one summary object
    \\  --fix -o <file.fasta>    Write a fixed FASTA to a separate output path
    \\  --fix-format-only        With --fix, proceed despite alphabet errors
    \\  --schema NAME            Header schema: uniprot or refseq
    \\  --custom-alphabet CHARS  Override nucleotide/protein alphabet checks
    \\  --max-header-len N       Warn on headers longer than N bytes (default: 1024)
    \\  Event list stops at 10000 issues (error without --fix). --fix still streams to -o.
    \\
    \\Examples:
    \\  z-fasta index genome.fa                  Create .zfi binary index
    \\  z-fasta index --emit-fai genome.fa       Output FAI to stdout
    \\  z-fasta get genome.fa chr1:1000-2000     Extract region
    \\  z-fasta get genome.fa chr1               Extract full sequence
    \\  z-fasta get genome.fa chr1:1000-2000 --rc
    \\  z-fasta get genome.fa --bed regions.bed  Extract BED regions
    \\  z-fasta get genome.fa --names ids.txt    Extract whole sequences from a file
    \\  z-fasta get genome.fa --bed regions.bed --strand-aware --summary
    \\  z-fasta stats genome.fa                  Show stats with composition
    \\  z-fasta validate genome.fa               Check FASTA validity
    \\  z-fasta validate --json --summary genome.fa
    \\
;

const GET_USAGE_ERROR =
    "error: usage: z-fasta get <file.fasta> (<region>... | --names file.txt|- | --bed file.bed|-)\n";

const ParseTransformFlagsError = error{ConflictingTransformFlags};

const ParsedTransformFlags = struct {
    orientation: getter.Orientation,
    annotate_transform: bool,
};

fn parseTransformFlags(
    rc: bool,
    complement_only: bool,
    reverse_only: bool,
    annotate_transform: bool,
) ParseTransformFlagsError!ParsedTransformFlags {
    var selected: usize = 0;
    if (rc) selected += 1;
    if (complement_only) selected += 1;
    if (reverse_only) selected += 1;
    if (selected > 1) return error.ConflictingTransformFlags;

    const orientation = if (rc)
        getter.Orientation.reverseComplement()
    else if (complement_only)
        getter.Orientation.complementOnly()
    else if (reverse_only)
        getter.Orientation.reverseOnly()
    else
        getter.Orientation{};

    return .{
        .orientation = orientation,
        .annotate_transform = annotate_transform,
    };
}

fn printUsageAndExit() noreturn {
    std.debug.print("{s}", .{USAGE});
    std.process.exit(1);
}

// Positional argv tokens must not start with `-`. Known flags are matched earlier.
fn rejectUnknownOption(arg: []const u8) void {
    if (arg.len > 0 and arg[0] == '-') {
        printErrorAndExit("error: unknown option: {s}\n", .{arg});
    }
}

fn printHelpAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, USAGE) catch {};
    std.process.exit(0);
}

fn printVersionAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, "z-fasta " ++ VERSION ++ "\n") catch {};
    std.process.exit(0);
}

fn writeStdoutAll(bytes: []const u8) void {
    if (comptime builtin.os.tag == .linux) {
        var written: usize = 0;
        while (written < bytes.len) {
            const result = std.os.linux.write(1, bytes[written..].ptr, bytes.len - written);
            if (std.os.linux.errno(result) != .SUCCESS) return;
            if (result == 0) return;
            written += result;
        }
    } else {
        var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
        defer threaded.deinit();
        std.Io.File.writeStreamingAll(.stdout(), threaded.io(), bytes) catch {};
    }
}

fn printHelpFastAndExit() noreturn {
    writeStdoutAll(USAGE);
    std.process.exit(0);
}

fn printVersionFastAndExit() noreturn {
    writeStdoutAll("z-fasta " ++ VERSION ++ "\n");
    std.process.exit(0);
}

pub fn main(init: std.process.Init.Minimal) void {
    var args = std.process.Args.Iterator.initAllocator(init.args, std.heap.page_allocator) catch {
        std.debug.print("error: out of memory\n", .{});
        std.process.exit(1);
    };
    defer args.deinit();
    _ = args.skip();

    const cmd = args.next() orelse {
        printUsageAndExit();
    };

    if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) {
        printHelpFastAndExit();
    }
    if (std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V")) {
        printVersionFastAndExit();
    }

    var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    if (std.mem.eql(u8, cmd, "index")) {
        runIndex(io, init.environ, &args);
    } else if (std.mem.eql(u8, cmd, "get")) {
        runGetCmd(io, &args);
    } else if (std.mem.eql(u8, cmd, "stats")) {
        runStatsCmd(io, &args);
    } else if (std.mem.eql(u8, cmd, "validate")) {
        runValidateCmd(io, &args);
    } else {
        printUsageAndExit();
    }
}

// --- Index command ---

fn runIndex(io: std.Io, environ: std.process.Environ, args: *std.process.Args.Iterator) void {
    var emit_fai = false;
    var enable_dedup = true;
    var fasta_path: ?[]const u8 = null;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")) {
            printVersionAndExit(io);
        } else if (std.mem.eql(u8, arg, "--emit-fai")) {
            emit_fai = true;
        } else if (std.mem.eql(u8, arg, "--no-dedup")) {
            enable_dedup = false;
        } else if (std.mem.eql(u8, arg, "--dedup")) {
            enable_dedup = true;
        } else {
            rejectUnknownOption(arg);
            if (fasta_path != null) {
                printErrorAndExit("error: index accepts exactly one FASTA path\n", .{});
            }
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printUsageAndExit();
    };

    var failure: index_format.diagnostics.Failure = undefined;
    const result = indexer.runIndex(io, environ, path, .{
        .emit_fai = emit_fai,
        .enable_dedup = enable_dedup,
    }, &failure) catch {
        printIndexFailure(path, &failure);
    };
    if (!emit_fai) std.debug.print("wrote {s}.zfi ({d} sequences)\n", .{ path, result.record_count });
}

fn printIndexFailure(path: []const u8, failure: *const index_format.diagnostics.Failure) noreturn {
    const primary = failure.primary() orelse printErrorAndExit("error: processing failed\n", .{});
    switch (primary.code()) {
        .@"io.file_not_found" => printErrorAndExit("error: file not found: {s}\n", .{path}),
        .@"io.access_denied" => printErrorAndExit("error: access denied: {s}\n", .{path}),
        .@"io.open_failed" => if (primary.at.source == .fasta) printErrorAndExit("error: failed to open file: {s}\n", .{path}) else printErrorAndExit("error: write failed\n", .{}),
        .@"fasta.no_records" => printErrorAndExit("error: file is empty: {s}\n", .{path}),
        .@"index.initial_header_required" => printErrorAndExit("error: not a FASTA file: {s}\n", .{path}),
        .@"index.identifier_too_long" => printErrorAndExit(
            "error: sequence name exceeds {d} bytes: {s}\n",
            .{ indexer.MAX_INDEX_NAME_LEN, path },
        ),
        .@"index.sequence_line_too_long" => printErrorAndExit(
            "error: sequence line exceeds the index format limit: {s}\n",
            .{path},
        ),
        .@"index.record_count_limit_exceeded" => printErrorAndExit(
            "error: sequence count exceeds {d}: {s}\n",
            .{ std.math.maxInt(u32), path },
        ),
        .@"index.nonuniform_fai" => printErrorAndExit(
            "error: cannot emit .fai for non-uniform sequence layout; run 'z-fasta index' (default) to write .zfi\n",
            .{},
        ),
        .@"index.no_indexable_records" => printErrorAndExit(
            "error: no valid sequences found in: {s}\n",
            .{path},
        ),
        .@"io.temp_unavailable" => if (primary.at.source != .fai_spool) printErrorAndExit("error: write failed\n", .{}) else printErrorAndExit(
            "error: no usable temporary directory for FAI spool\n",
            .{},
        ),
        .@"input.changed" => printErrorAndExit(
            "error: source changed while indexing: {s}\n",
            .{path},
        ),
        .@"index.source_timestamp_out_of_range" => printErrorAndExit(
            "error: source modification time predates the Unix epoch: {s}\n",
            .{path},
        ),
        .@"io.invalid_path" => printErrorAndExit("error: path too long\n", .{}),
        .@"io.publish_failed" => printErrorAndExit(
            "error: failed to finalize index: {s}.zfi\n",
            .{path},
        ),
        .@"io.stat_failed" => if (primary.at.source == .fasta)
            printErrorAndExit("error: failed to stat file: {s}\n", .{path})
        else
            printErrorAndExit("error: failed to read temporary FAI spool\n", .{}),
        .@"io.read_failed", .@"io.unexpected_eof" => if (primary.at.source == .fasta)
            printErrorAndExit("error: failed to read file: {s}\n", .{path})
        else if (primary.at.source == .fai_spool)
            printErrorAndExit("error: failed to read temporary FAI spool\n", .{})
        else
            printErrorAndExit("error: write failed\n", .{}),
        .@"io.write_failed", .@"io.flush_failed" => if (primary.at.source == .fai_spool)
            printErrorAndExit("error: failed to write temporary FAI spool\n", .{})
        else if (primary.at.source == .stdout and primary.stage == .flush)
            printErrorAndExit("error: failed to flush FAI stdout\n", .{})
        else if (primary.at.source == .stdout)
            printErrorAndExit("error: failed to replay FAI stdout\n", .{})
        else
            printErrorAndExit("error: write failed\n", .{}),
        .@"io.create_failed" => printErrorAndExit("error: write failed\n", .{}),
        else => printErrorAndExit("error: processing failed\n", .{}),
    }
}

// --- Get command ---

fn runGetCmd(io: std.Io, args: *std.process.Args.Iterator) void {
    var fasta_path: ?[]const u8 = null;
    var bed_path: ?[]const u8 = null;
    var names_path: ?[]const u8 = null;
    var honor_strand = false;
    var summary = false;
    var rc = false;
    var complement_only = false;
    var reverse_only = false;
    var annotate_transform = false;
    // Static buffer: up to 1024 region strings without heap allocation.
    var region_buf: [1024][]const u8 = undefined;
    var region_count: usize = 0;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--bed")) {
            if (bed_path != null) {
                printErrorAndExit("error: --bed may be specified only once\n", .{});
            }
            bed_path = args.next() orelse {
                printErrorAndExit("error: --bed requires a path or '-'\n", .{});
            };
        } else if (std.mem.eql(u8, arg, "--names")) {
            if (names_path != null) {
                printErrorAndExit("error: --names may be specified only once\n", .{});
            }
            names_path = args.next() orelse {
                printErrorAndExit("error: --names requires a path\n", .{});
            };
        } else if (std.mem.eql(u8, arg, STRAND_AWARE_FLAG) or std.mem.eql(u8, arg, STRAND_AWARE_ALIAS)) {
            honor_strand = true;
        } else if (std.mem.eql(u8, arg, "--summary")) {
            summary = true;
        } else if (std.mem.eql(u8, arg, RC_FLAG)) {
            rc = true;
        } else if (std.mem.eql(u8, arg, COMPLEMENT_ONLY_FLAG)) {
            complement_only = true;
        } else if (std.mem.eql(u8, arg, REVERSE_ONLY_FLAG)) {
            reverse_only = true;
        } else if (std.mem.eql(u8, arg, ANNOTATE_RC_FLAG)) {
            annotate_transform = true;
        } else if (fasta_path == null) {
            rejectUnknownOption(arg);
            fasta_path = arg;
        } else {
            rejectUnknownOption(arg);
            if (region_count >= region_buf.len) {
                printErrorAndExit("error: too many regions (max 1024)\n", .{});
            }
            region_buf[region_count] = arg;
            region_count += 1;
        }
    }

    const path = fasta_path orelse {
        printErrorAndExit(GET_USAGE_ERROR, .{});
    };
    const source_count = @as(usize, @intFromBool(region_count > 0)) +
        @as(usize, @intFromBool(bed_path != null)) +
        @as(usize, @intFromBool(names_path != null));
    if (source_count == 0) {
        printErrorAndExit(GET_USAGE_ERROR, .{});
    }
    if (source_count > 1) {
        printErrorAndExit("error: choose one request source: positional regions, --names, or --bed\n", .{});
    }
    if (honor_strand and bed_path == null) {
        printErrorAndExit("error: --strand-aware requires --bed\n", .{});
    }

    const transform = parseTransformFlags(rc, complement_only, reverse_only, annotate_transform) catch {
        printErrorAndExit("error: --rc, --complement-only, and --reverse-only are mutually exclusive\n", .{});
    };
    if (transform.annotate_transform and transform.orientation.isIdentity() and !honor_strand) {
        printErrorAndExit("error: --annotate-rc requires a transform or --strand-aware BED input\n", .{});
    }

    const source: getter.RequestSource = if (bed_path) |bed|
        .{ .bed = bed }
    else if (names_path) |names|
        .{ .names = names }
    else
        .{ .positional = region_buf[0..region_count] };

    const options = getter.GetOptions{
        .source = source,
        .honor_strand = honor_strand,
        .measure_elapsed = summary,
        .orientation = transform.orientation,
        .annotate_transform = transform.annotate_transform,
    };
    var failure: index_format.diagnostics.Failure = undefined;
    const result = getter.runGetWithOptions(std.heap.page_allocator, io, path, options, &failure) catch {
        printGetFailure(path, options, &failure);
    };
    if (summary) {
        failure.sources.set(.stderr, .{ .label = .{ .borrowed = "stderr" }, .kind = .stream });
        var buffer: [512]u8 = undefined;
        var writer = std.Io.File.Writer.initStreaming(.stderr(), io, &buffer);
        failure.outputs.stderr.state = .streaming;
        failure.outputs.stderr.may_be_partial = true;
        writeGetSummary(&writer.interface, result) catch |err| {
            const captured = failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .at = .{ .source = .stderr }, .cause = .{ .zig = writer.err orelse err } });
            switch (captured) {
                else => {},
            }
            printGetFailure(path, options, &failure);
        };
        writer.flush() catch |err| {
            failure.outputs.stderr.final_flush = .failed;
            const captured = failure.fail(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .stderr }, .cause = .{ .zig = writer.err orelse err } });
            switch (captured) {
                else => {},
            }
            printGetFailure(path, options, &failure);
        };
        failure.outputs.stderr.final_flush = .succeeded;
    }
}

fn writeGetSummary(writer: *std.Io.Writer, result: getter.GetResult) std.Io.Writer.Error!void {
    const seconds = if (result.elapsed_ns == 0) 0.0 else @as(f64, @floatFromInt(result.elapsed_ns)) / @as(f64, @floatFromInt(std.time.ns_per_s));
    const regions_per_second = if (seconds == 0.0) 0.0 else @as(f64, @floatFromInt(result.region_count)) / seconds;
    try writer.print("summary: regions={d} total_bases={d} elapsed_s={d:.6} regions_per_s={d:.1}\n", .{ result.region_count, result.total_bases, seconds, regions_per_second });
}

fn getFailureIdentifier(options: getter.GetOptions, primary: *const index_format.diagnostics.Diagnostic, failure: *const index_format.diagnostics.Failure, preview: ?index_format.diagnostics.Preview) struct { bytes: []const u8, original_bytes: u64, truncated: bool } {
    if (options.source == .positional) {
        const requests = options.source.positional;
        if (primary.at.request) |ordinal| {
            if (ordinal > 0 and ordinal <= requests.len) {
                const name = getter.parseRegion(requests[@intCast(ordinal - 1)]).name;
                return .{ .bytes = name, .original_bytes = name.len, .truncated = false };
            }
        }
    }
    if (preview) |reference| return .{ .bytes = failure.preview(reference), .original_bytes = reference.original_bytes, .truncated = reference.truncated };
    return .{ .bytes = "", .original_bytes = 0, .truncated = false };
}

// Temporary legacy messages until the shared CLI renderer conversion.
fn printGetFailure(path: []const u8, options: getter.GetOptions, failure: *const index_format.diagnostics.Failure) noreturn {
    const primary = failure.primary() orelse printErrorAndExit("error: processing failed\n", .{});
    const source = primary.at.source;
    const label = if (source) |id| failure.sourceLabel(id) orelse path else path;
    const line = primary.at.location.line orelse 0;
    switch (primary.details) {
        .@"request.sequence_not_found" => |detail| {
            const identifier = getFailureIdentifier(options, primary, failure, detail.identifier);
            if (identifier.truncated) printErrorAndExit("error: sequence not found: {s} (identifier truncated; {d} bytes)\n", .{ identifier.bytes, identifier.original_bytes });
            printErrorAndExit("error: sequence not found: {s}\n", .{identifier.bytes});
        },
        .@"request.start_before_one" => printErrorAndExit("error: start position must be >= 1\n", .{}),
        .@"request.start_after_sequence" => |detail| printErrorAndExit("error: start position {d} exceeds sequence length {d}\n", .{ detail.start, detail.sequence_length }),
        .@"request.end_before_start" => printErrorAndExit("error: end position must be >= start position\n", .{}),
        .@"request.no_regions" => printErrorAndExit("error: no regions provided\n", .{}),
        .@"request.identifier_too_long" => |detail| if (options.source == .bed)
            printErrorAndExit("error: invalid BED line {d}: chrom exceeds {d} bytes\n", .{ line, detail.limit_bytes })
        else
            printErrorAndExit("error: name at line {d} exceeds {d} bytes\n", .{ line, detail.limit_bytes }),
        .@"request.bed.field_missing" => |detail| printErrorAndExit("error: invalid BED line {d}: missing {s}\n", .{ line, @tagName(detail.field) }),
        .@"request.bed.invalid_coordinate" => |detail| printErrorAndExit("error: invalid BED line {d}: invalid {s}\n", .{ line, @tagName(detail.field) }),
        .@"request.bed.coordinate_out_of_range" => |detail| printErrorAndExit("error: invalid BED line {d}: invalid {s}\n", .{ line, @tagName(detail.field) }),
        .@"request.bed.invalid_interval" => printErrorAndExit("error: invalid BED line {d}: end must be greater than start\n", .{line}),
        .@"request.bed.invalid_strand" => printErrorAndExit("error: invalid BED line {d}: invalid strand\n", .{line}),
        .@"transform.complement_on_protein" => |detail| {
            const identifier = getFailureIdentifier(options, primary, failure, detail.identifier);
            if (identifier.truncated) printErrorAndExit("error: reverse complement is not defined for protein sequences: {s} (identifier truncated; {d} bytes) (classified from up to {d} bases)\n", .{ identifier.bytes, identifier.original_bytes, stats.GET_TYPE_SAMPLE_BASES });
            printErrorAndExit("error: reverse complement is not defined for protein sequences: {s} (classified from up to {d} bases)\n", .{ identifier.bytes, stats.GET_TYPE_SAMPLE_BASES });
        },
        .@"resource.line_buffer_exceeded" => |detail| if (options.source == .bed)
            printErrorAndExit("error: BED line {d} exceeds {d}-byte reader buffer (no newline within limit)\n", .{ line, detail.buffer_bytes })
        else
            printErrorAndExit("error: name at line {d} exceeds {d} bytes\n", .{ line, std.math.maxInt(u16) }),
        .@"resource.out_of_memory" => if (source == .fasta or source == .zfi or source == .fai)
            printErrorAndExit("error: out of memory loading index\n", .{})
        else
            printErrorAndExit("error: out of memory\n", .{}),
        .@"io.file_not_found" => printErrorAndExit("error: file not found: {s}\n", .{label}),
        .@"io.access_denied" => printErrorAndExit("error: access denied: {s}\n", .{label}),
        .@"io.invalid_path" => printErrorAndExit("error: path too long\n", .{}),
        .@"io.open_failed" => if (source == .fasta or source == .names or source == .bed)
            printErrorAndExit("error: failed to open file: {s}\n", .{label})
        else
            printErrorAndExit("error: failed to load index for: {s}\n", .{path}),
        .@"io.map_failed" => printErrorAndExit("error: failed to mmap file: {s}\n", .{path}),
        .@"fasta.no_records" => printErrorAndExit("error: file is empty: {s}\n", .{path}),
        .@"index.not_found" => printErrorAndExit("error: no index found for {s}. Run 'z-fasta index {s}' first.\n", .{ path, path }),
        .@"index.older_than_source", .@"index.source_size_mismatch", .@"index.source_mtime_mismatch" => printErrorAndExit("error: index is stale (FASTA size or mtime does not match the index). Re-run 'z-fasta index <path>'.\n", .{}),
        .@"io.write_failed", .@"io.flush_failed" => printErrorAndExit("error: write failed\n", .{}),
        .@"io.read_failed" => if (source == .names or source == .bed or source == .stdin)
            if (options.source == .bed) printErrorAndExit("error: failed to read BED input\n", .{}) else printErrorAndExit("error: failed to read names input\n", .{})
        else if (source == .fasta and primary.stage == .read and primary.related != null)
            printErrorAndExit("error: failed to read FASTA\n", .{})
        else
            printErrorAndExit("error: failed to load index for: {s}\n", .{path}),
        .@"io.unexpected_eof", .@"pair.region_length_mismatch", .@"index.sequence_span_out_of_bounds" => if (primary.related != null)
            printErrorAndExit("error: read past end of FASTA\n", .{})
        else
            printErrorAndExit("error: corrupt index file for: {s}\n", .{path}),
        .@"input.changed" => printErrorAndExit("error: request input changed while extracting: {s}\n", .{label}),
        else => if (index_format.diagnostics.metadata(primary.code()).category == .index_integrity or index_format.diagnostics.metadata(primary.code()).category == .index_compatibility)
            printErrorAndExit("error: corrupt index file for: {s}\n", .{path})
        else
            printErrorAndExit("error: processing failed\n", .{}),
    }
}

// --- Stats command ---

fn runStatsCmd(io: std.Io, args: *std.process.Args.Iterator) void {
    var fasta_path: ?[]const u8 = null;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else {
            rejectUnknownOption(arg);
            if (fasta_path != null) {
                printErrorAndExit("error: stats accepts exactly one FASTA path\n", .{});
            }
            fasta_path = arg;
        }
    }

    const path = fasta_path orelse {
        printErrorAndExit("error: usage: z-fasta stats <file.fasta>\n", .{});
    };

    var failure: index_format.diagnostics.Failure = undefined;
    stats.runStats(std.heap.page_allocator, io, path, &failure) catch printStatsFailure(path, &failure);
}

// Temporary legacy messages; common rendering and status move in commit 10.
fn printStatsFailure(path: []const u8, failure: *const index_format.diagnostics.Failure) noreturn {
    const primary = failure.primary().?;
    switch (primary.code()) {
        .@"io.file_not_found" => printErrorAndExit("error: file not found: {s}\n", .{path}),
        .@"io.access_denied" => printErrorAndExit("error: access denied: {s}\n", .{path}),
        .@"io.invalid_path" => printErrorAndExit("error: path too long\n", .{}),
        .@"fasta.no_records" => printErrorAndExit("error: file is empty: {s}\n", .{path}),
        .@"io.map_failed" => printErrorAndExit("error: failed to mmap file: {s}\n", .{path}),
        .@"index.older_than_source", .@"index.source_size_mismatch", .@"index.source_mtime_mismatch" => printErrorAndExit("error: index is stale (FASTA size or mtime does not match the index). Re-run 'z-fasta index <path>'.\n", .{}),
        .@"index.not_found" => printErrorAndExit("error: no index found for {s}. Run 'z-fasta index {s}' first.\n", .{ path, path }),
        .@"resource.out_of_memory" => if (primary.stage == .allocate and (primary.details.@"resource.out_of_memory".allocation_role == .lengths or primary.details.@"resource.out_of_memory".allocation_role == .traversal_order))
            printErrorAndExit("error: out of memory\n", .{})
        else
            printErrorAndExit("error: out of memory loading index\n", .{}),
        .@"resource.count_overflow" => printErrorAndExit("error: sequence length overflow\n", .{}),
        .@"io.write_failed", .@"io.flush_failed" => printErrorAndExit("error: write failed\n", .{}),
        .@"pair.total_length_mismatch" => printErrorAndExit("error: index sequence lengths do not match FASTA data\n", .{}),
        .@"io.read_failed" => if (primary.stage == .read and primary.at.source == .fasta)
            printErrorAndExit("error: failed to read FASTA\n", .{})
        else
            printErrorAndExit("error: failed to load index for: {s}\n", .{path}),
        .@"io.unexpected_eof" => if (primary.at.source == .fasta)
            printErrorAndExit("error: index sequence lengths do not match FASTA data\n", .{})
        else
            printErrorAndExit("error: corrupt index file for: {s}\n", .{path}),
        .@"internal.consistency_failed" => switch (primary.details.@"internal.consistency_failed".check) {
            .retained_shortest_name => printErrorAndExit("error: index does not retain the shortest record name\n", .{}),
            .retained_longest_name => printErrorAndExit("error: index does not retain the longest record name\n", .{}),
            .composition_accounting => printErrorAndExit("error: composition does not match indexed symbols\n", .{}),
            else => printErrorAndExit("error: index sequence lengths do not match FASTA data\n", .{}),
        },
        .@"io.open_failed", .@"io.stat_failed" => printErrorAndExit("error: failed to load index for: {s}\n", .{path}),
        .@"resource.capacity_limit_exceeded" => if (primary.details.@"resource.capacity_limit_exceeded".capacity_role == .map_bytes)
            printErrorAndExit("error: failed to mmap file: {s}\n", .{path})
        else
            printErrorAndExit("error: out of memory loading index\n", .{}),
        .@"index.no_records" => if (primary.stage == .scan)
            printErrorAndExit("error: no sequences in index\n", .{})
        else
            printErrorAndExit("error: corrupt index file for: {s}\n", .{path}),
        else => printErrorAndExit("error: corrupt index file for: {s}\n", .{path}),
    }
}

// --- Validate command ---

fn runValidateCmd(io: std.Io, args: *std.process.Args.Iterator) void {
    var fasta_path: ?[]const u8 = null;
    var options = validator.Options{};
    var saw_json = false;
    var saw_summary = false;

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        } else if (std.mem.eql(u8, arg, "--strict")) {
            options.strict = true;
        } else if (std.mem.eql(u8, arg, "--json")) {
            saw_json = true;
            options.output_mode = if (saw_summary) .json_summary else .json_lines;
        } else if (std.mem.eql(u8, arg, "--summary")) {
            saw_summary = true;
            options.output_mode = .json_summary;
        } else if (std.mem.eql(u8, arg, "--fix")) {
            options.fix = true;
        } else if (std.mem.eql(u8, arg, "--fix-format-only")) {
            options.fix_format_only = true;
        } else if (std.mem.eql(u8, arg, "-o")) {
            options.output_path = args.next() orelse {
                printErrorAndExit("error: -o requires an output path\n", .{});
            };
        } else if (std.mem.eql(u8, arg, "--schema")) {
            const raw = args.next() orelse {
                printErrorAndExit("error: --schema requires uniprot or refseq\n", .{});
            };
            options.schema = parseValidateSchema(raw);
        } else if (std.mem.eql(u8, arg, "--custom-alphabet")) {
            options.custom_alphabet = args.next() orelse {
                printErrorAndExit("error: --custom-alphabet requires characters\n", .{});
            };
        } else if (std.mem.eql(u8, arg, "--max-header-len")) {
            const raw = args.next() orelse {
                printErrorAndExit("error: --max-header-len requires a positive integer\n", .{});
            };
            options.max_header_len = parsePositiveUsize(raw, "--max-header-len");
        } else if (fasta_path == null) {
            rejectUnknownOption(arg);
            fasta_path = arg;
        } else {
            rejectUnknownOption(arg);
            printErrorAndExit("error: validate accepts exactly one FASTA path\n", .{});
        }
    }

    if (saw_summary and !saw_json) {
        printErrorAndExit("error: validate --summary requires --json\n", .{});
    }
    if (options.fix_format_only and !options.fix) {
        printErrorAndExit("error: validate --fix-format-only requires --fix\n", .{});
    }
    if (options.output_path != null and !options.fix) {
        printErrorAndExit("error: validate -o requires --fix\n", .{});
    }

    const path = fasta_path orelse {
        printErrorAndExit("error: usage: z-fasta validate [options] <file.fasta>\n", .{});
    };

    var failure: index_format.diagnostics.Failure = undefined;
    const result = validator.runValidate(std.heap.page_allocator, io, path, options, &failure) catch printValidateFailure(path, options, &failure);
    if (result.truncated) {
        var buffer: [512]u8 = undefined;
        var stderr_fw = std.Io.File.Writer.initStreaming(.stderr(), io, &buffer);
        failure.sources.set(.stderr, .{ .label = .{ .borrowed = "stderr" }, .kind = .stream });
        failure.outputs.stderr.state = .streaming;
        failure.outputs.stderr.may_be_partial = true;
        writeValidateCapNotice(&stderr_fw.interface, result.fix_written) catch |err| {
            switch (failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .report, .at = .{ .source = .stderr }, .cause = .{ .zig = stderr_fw.err orelse err } })) {
                else => {},
            }
            printValidateFailure(path, options, &failure);
        };
        stderr_fw.flush() catch |err| {
            failure.outputs.stderr.final_flush = .failed;
            switch (failure.fail(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .stderr }, .cause = .{ .zig = stderr_fw.err orelse err } })) {
                else => {},
            }
            printValidateFailure(path, options, &failure);
        };
        failure.outputs.stderr.final_flush = .succeeded;
    }
    std.process.exit(result.status);
}

fn writeValidateCapNotice(writer: *std.Io.Writer, fixed: bool) std.Io.Writer.Error!void {
    if (fixed) try writer.print("warning: validate event list truncated at {d}; --fix output was still written\n", .{validator.MAX_VALIDATE_EVENTS}) else try writer.print("error: validate retained {d} events; scanning and counting completed; fix reported issues and re-run\n", .{validator.MAX_VALIDATE_EVENTS});
}

// Temporary legacy messages; C10 selects common rendering/status after cleanup.
fn printValidateFailure(path: []const u8, options: validator.Options, failure: *const index_format.diagnostics.Failure) noreturn {
    const primary = failure.primary().?;
    switch (primary.code()) {
        .@"io.file_not_found" => printErrorAndExit("error: file not found: {s}\n", .{path}),
        .@"io.access_denied" => printErrorAndExit("error: access denied: {s}\n", .{path}),
        .@"io.open_failed", .@"io.invalid_path" => printErrorAndExit("error: failed to open file: {s}\n", .{path}),
        .@"io.stat_failed" => printErrorAndExit("error: failed to stat file: {s}\n", .{path}),
        .@"io.map_failed", .@"resource.capacity_limit_exceeded" => printErrorAndExit("error: failed to mmap file: {s}\n", .{path}),
        .@"resource.out_of_memory" => printErrorAndExit("error: out of memory\n", .{}),
        .@"index.sequence_line_too_long" => printErrorAndExit("error: sequence line exceeds the index format limit: {s}\n", .{path}),
        .@"repair.output_is_input" => printErrorAndExit("error: validate --fix will not overwrite input: {s}\n", .{path}),
        .@"repair.output_identity_unknown" => printErrorAndExit("error: validate --fix cannot establish output identity: {s}\n", .{options.output_path orelse ""}),
        .@"repair.blocked" => switch (primary.details.@"repair.blocked".blocker) {
            .content_before_header => printErrorAndExit("error: validate --fix refuses content before the first header at line {d}\n", .{primary.at.location.line orelse 1}),
            .invalid_symbol => printErrorAndExit("error: validate --fix refuses character-level errors; use --fix-format-only to keep them unchanged\n", .{}),
            .duplicate_identifier => printErrorAndExit("error: validate --fix cannot repair duplicate_name\n", .{}),
            .no_records => printErrorAndExit("error: validate --fix cannot repair no_sequences\n", .{}),
            .header_control_byte, .sequence_control_byte => printErrorAndExit("error: validate --fix cannot repair null_byte\n", .{}),
        },
        .@"usage.option_requires" => printErrorAndExit("error: validate --fix requires -o <output.fa>\n", .{}),
        .@"io.create_failed", .@"io.publish_failed", .@"input.changed" => printErrorAndExit("error: failed to write fixed FASTA: {s}\n", .{options.output_path orelse ""}),
        .@"io.write_failed", .@"io.flush_failed" => if (primary.at.source == .output)
            printErrorAndExit("error: failed to write fixed FASTA: {s}\n", .{options.output_path orelse ""})
        else
            printErrorAndExit("error: write failed\n", .{}),
        else => printErrorAndExit("error: validation failed: {s}\n", .{path}),
    }
}

fn parseValidateSchema(raw: []const u8) validator.Schema {
    if (std.mem.eql(u8, raw, "uniprot")) return .uniprot;
    if (std.mem.eql(u8, raw, "refseq")) return .refseq;
    printErrorAndExit("error: --schema must be uniprot or refseq\n", .{});
}

fn parsePositiveUsize(raw: []const u8, comptime flag: []const u8) usize {
    const parsed = std.fmt.parseInt(usize, raw, 10) catch {
        printErrorAndExit("error: {s} requires a positive integer\n", .{flag});
    };
    if (parsed == 0) {
        printErrorAndExit("error: {s} requires a positive integer\n", .{flag});
    }
    return parsed;
}

test "[unit] - [get transform flags]: selects the requested orientation" {
    const rc = try parseTransformFlags(true, false, false, false);
    const complement_only = try parseTransformFlags(false, true, false, true);
    const reverse_only = try parseTransformFlags(false, false, true, false);

    try std.testing.expect(rc.orientation.reverse);
    try std.testing.expect(rc.orientation.complement);
    try std.testing.expect(!complement_only.orientation.reverse);
    try std.testing.expect(complement_only.orientation.complement);
    try std.testing.expect(complement_only.annotate_transform);
    try std.testing.expect(reverse_only.orientation.reverse);
    try std.testing.expect(!reverse_only.orientation.complement);
}

test "[failure] - [get transform flags]: rejects conflicting modes" {
    const conflicts = [_][3]bool{
        .{ true, true, false },
        .{ true, false, true },
        .{ false, true, true },
        .{ true, true, true },
    };

    for (conflicts) |flags| {
        try std.testing.expectError(
            error.ConflictingTransformFlags,
            parseTransformFlags(flags[0], flags[1], flags[2], false),
        );
    }
}

test "[cli] - [dispatch]: writes help and rejects missing or unknown commands" {
    const zfasta_bin = if (builtin.os.tag == .windows) "zig-out\\bin\\z-fasta.exe" else "zig-out/bin/z-fasta";
    const allocator = std.testing.allocator;
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    const help = try std.process.run(allocator, io, .{
        .argv = &.{ zfasta_bin, "--help" },
        .stdout_limit = .limited(64 * 1024),
        .stderr_limit = .limited(64 * 1024),
    });
    defer allocator.free(help.stdout);
    defer allocator.free(help.stderr);
    switch (help.term) {
        .exited => |code| try std.testing.expectEqual(@as(u8, 0), code),
        else => return error.ChildProcessFailed,
    }
    try std.testing.expectEqual(@as(usize, 0), help.stderr.len);
    try std.testing.expectEqualStrings(USAGE, help.stdout);

    const cases = [_][]const []const u8{
        &.{zfasta_bin},
        &.{ zfasta_bin, "index" },
        &.{ zfasta_bin, "not-a-command" },
    };
    for (cases) |argv| {
        const result = try std.process.run(allocator, io, .{
            .argv = argv,
            .stdout_limit = .limited(64 * 1024),
            .stderr_limit = .limited(64 * 1024),
        });
        defer allocator.free(result.stdout);
        defer allocator.free(result.stderr);

        switch (result.term) {
            .exited => |code| try std.testing.expectEqual(@as(u8, 1), code),
            else => return error.ChildProcessFailed,
        }
        try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
        try std.testing.expectEqualStrings(help.stdout, result.stderr);
    }
}

test "[failure] - [get summary]: propagates a full destination" {
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try std.testing.expectError(error.WriteFailed, writeGetSummary(&writer, .{ .region_count = 1, .total_bases = 1, .elapsed_ns = 1 }));
}
