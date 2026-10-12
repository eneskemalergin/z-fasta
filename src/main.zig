//! CLI entry and subcommand dispatcher for `index`, `get`, `stats`, and `validate`.

const std = @import("std");
const builtin = @import("builtin");

pub const index_format = @import("index_format.zig");
pub const indexer = @import("indexer.zig");
pub const getter = @import("getter.zig");
pub const bed_parser = @import("bed_parser.zig");
pub const stats = @import("stats.zig");
pub const validator = @import("validator.zig");

const diagnostics = @import("diagnostics");
pub const reporting = @import("reporting.zig");

pub fn main(init: std.process.Init.Minimal) void {
    const status = runCli(init);
    if (status != 0) std.process.exit(status);
}

const TEST_BINARY_PATH = if (builtin.is_test) blk: {
    const options = @import("test_options");
    break :blk if (std.fs.path.isAbsolute(options.binary_path)) options.binary_path else options.build_root ++ std.fs.path.sep_str ++ options.binary_path;
} else "";

const VERSION = "0.3.5";
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
    \\  --error-format=text|json  Select stderr diagnostics (default: text)
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
    \\  Retain 10000 examples; scan/count all input (error without --fix).
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

const ErrorFormat = enum { text, json };

fn selectErrorFormat(args: std.process.Args.Iterator, operation: diagnostics.Operation, initial: ErrorFormat) error{InvalidFormat}!ErrorFormat {
    var probe = args;
    var format = initial;
    while (probe.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h") or (operation == .index and (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")))) break;
        if ((std.mem.eql(u8, arg, "--error-format") or std.mem.startsWith(u8, arg, "--error-format="))) {
            format = try parseErrorFormat(arg);
        } else if (takesValue(operation, arg)) {
            _ = probe.next();
        }
    }
    return format;
}

fn takesValue(operation: diagnostics.Operation, arg: []const u8) bool {
    return switch (operation) {
        .get => std.mem.eql(u8, arg, "--bed") or std.mem.eql(u8, arg, "--names"),
        .validate => std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--schema") or std.mem.eql(u8, arg, "--custom-alphabet") or std.mem.eql(u8, arg, "--max-header-len"),
        else => false,
    };
}

fn parseErrorFormat(arg: []const u8) error{InvalidFormat}!ErrorFormat {
    if (std.mem.eql(u8, arg, "--error-format=text")) return .text;
    if (std.mem.eql(u8, arg, "--error-format=json")) return .json;
    return error.InvalidFormat;
}

const Cli = struct {
    io: std.Io,
    args: *std.process.Args.Iterator,
    format: ErrorFormat,
    failure: *diagnostics.Failure,

    fn next(self: *Cli) ?[]const u8 {
        while (self.args.next()) |arg| {
            if (std.mem.startsWith(u8, arg, "--error-format=")) continue;
            return arg;
        }
        return null;
    }

    fn fail(self: *Cli, details: diagnostics.Details) diagnostics.Error {
        return self.failure.fail(.{ .details = details, .stage = .arguments });
    }

    fn value(self: *Cli, option: diagnostics.Argument, expected: diagnostics.Expected) diagnostics.Error![]const u8 {
        return self.args.next() orelse self.fail(.{ .@"usage.missing_option_value" = .{ .option = option, .expected = expected } });
    }

    fn positional(self: *Cli, arg: []const u8) diagnostics.Error!void {
        if (arg.len != 0 and arg[0] == '-') return self.fail(.{ .@"usage.unknown_option" = .{ .option = self.failure.retain(arg, 32) } });
    }

    fn stdout(self: *Cli, bytes: []const u8) diagnostics.Error!u8 {
        self.failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
        self.failure.outputs.stdout.state = .streaming;
        self.failure.outputs.stdout.may_be_partial = true;
        std.Io.File.writeStreamingAll(.stdout(), self.io, bytes) catch |err| return self.failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .at = .{ .source = .stdout }, .cause = .{ .zig = err } });
        self.failure.succeed();
        return 0;
    }
};

fn reportFailure(io: std.Io, format: ErrorFormat, failure: *diagnostics.Failure) u8 {
    var buffer: [4096]u8 = undefined;
    var fw = std.Io.File.Writer.initStreaming(.stderr(), io, &buffer);
    if (failure.sources.get(.stderr) == null) failure.sources.set(.stderr, .{ .label = .{ .borrowed = "stderr" }, .kind = .stream });
    failure.outputs.stderr.state = .streaming;
    failure.outputs.stderr.may_be_partial = true;
    const rendered = switch (format) {
        .text => reporting.writeText(&fw.interface, failure),
        .json => reporting.writeJson(&fw.interface, failure, .{ .version = VERSION }),
    };
    rendered catch |err| {
        fw.interface.end = 0;
        switch (failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .report, .at = .{ .source = .stderr }, .cause = .{ .zig = fw.err orelse err } })) {
            else => {},
        }
        return 1;
    };
    fw.flush() catch |err| {
        fw.interface.end = 0;
        failure.outputs.stderr.final_flush = .failed;
        switch (failure.fail(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .stderr }, .cause = .{ .zig = fw.err orelse err } })) {
            else => {},
        }
        return 1;
    };
    failure.outputs.stderr.final_flush = .succeeded;
    return 1;
}

fn bootstrap(reason: reporting.Bootstrap) u8 {
    var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
    defer threaded.deinit();
    var buffer: [512]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    reporting.writeBootstrap(&writer, reason) catch return 1;
    std.Io.File.writeStreamingAll(.stderr(), threaded.io(), writer.buffered()) catch {};
    return 1;
}

const DirectWriteFailure = struct { cause: ?diagnostics.Cause };

fn writeStdoutDirect(bytes: []const u8) ?DirectWriteFailure {
    if (comptime builtin.os.tag != .linux) {
        var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
        defer threaded.deinit();
        std.Io.File.writeStreamingAll(.stdout(), threaded.io(), bytes) catch |err| return .{ .cause = .{ .zig = err } };
        return null;
    }
    var written: usize = 0;
    while (written < bytes.len) {
        const result = std.os.linux.write(1, bytes[written..].ptr, bytes.len - written);
        const err = std.os.linux.errno(result);
        switch (err) {
            .SUCCESS => if (result == 0) return .{ .cause = null },
            .INTR => continue,
            else => return .{ .cause = .{ .native = .{ .domain = .posix_errno, .number = @intFromEnum(err) } } },
        }
        written += result;
    }
    return null;
}

fn runCli(init: std.process.Init.Minimal) u8 {
    var args = std.process.Args.Iterator.initAllocator(init.args, std.heap.page_allocator) catch return bootstrap(.argv_oom);
    defer args.deinit();
    _ = args.skip();
    var format: ErrorFormat = .text;
    var command = args.next();
    while (command) |arg| {
        if (!(std.mem.eql(u8, arg, "--error-format") or std.mem.startsWith(u8, arg, "--error-format="))) break;
        format = parseErrorFormat(arg) catch return bootstrap(.invalid_format);
        command = args.next();
    }
    const terminal_command = if (command) |cmd| std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h") or std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V") else false;
    if (command) |cmd| {
        if (format == .text and terminal_command) {
            const bytes = if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) USAGE else "z-fasta " ++ VERSION ++ "\n";
            const write_failure = writeStdoutDirect(bytes) orelse return 0;
            var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
            defer threaded.deinit();
            var failure: diagnostics.Failure = undefined;
            failure.reset(.cli);
            failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
            failure.outputs.stdout.state = .streaming;
            failure.outputs.stdout.may_be_partial = true;
            switch (failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .at = .{ .source = .stdout }, .cause = write_failure.cause })) {
                else => {},
            }
            return reportFailure(threaded.io(), .text, &failure);
        }
    }
    const operation: diagnostics.Operation = if (command) |cmd| std.meta.stringToEnum(diagnostics.Operation, cmd) orelse .cli else .cli;
    if (!terminal_command) format = selectErrorFormat(args, operation, format) catch return bootstrap(.invalid_format);
    var threaded = std.Io.Threaded.init(std.heap.page_allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    var failure: diagnostics.Failure = undefined;
    failure.reset(operation);
    var cli = Cli{ .io = io, .args = &args, .format = format, .failure = &failure };
    const status = dispatch(&cli, init.environ, command) catch {
        if (failure.primary().?.at.source == .stderr) return 1;
        return reportFailure(io, format, &failure);
    };
    return status;
}

fn dispatch(cli: *Cli, environ: std.process.Environ, command: ?[]const u8) diagnostics.Error!u8 {
    const cmd = command orelse return cli.fail(.{ .@"usage.command_required" = {} });
    if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) return cli.stdout(USAGE);
    if (std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V")) return cli.stdout("z-fasta " ++ VERSION ++ "\n");
    return switch (cli.failure.operation) {
        .index => runIndex(cli, environ),
        .get => runGetCmd(cli),
        .stats => runStatsCmd(cli),
        .validate => runValidateCmd(cli),
        .cli => cli.fail(.{ .@"usage.unknown_command" = .{ .command = cli.failure.retain(cmd, 32) } }),
    };
}

fn runIndex(cli: *Cli, environ: std.process.Environ) diagnostics.Error!u8 {
    var emit_fai = false;
    var enable_dedup = true;
    var path: ?[]const u8 = null;
    while (cli.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) return cli.stdout(USAGE);
        if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")) return cli.stdout("z-fasta " ++ VERSION ++ "\n");
        if (std.mem.eql(u8, arg, "--emit-fai")) emit_fai = true else if (std.mem.eql(u8, arg, "--no-dedup")) enable_dedup = false else if (std.mem.eql(u8, arg, "--dedup")) enable_dedup = true else {
            try cli.positional(arg);
            if (path != null) return cli.fail(.{ .@"usage.extra_argument" = .{ .argument = .fasta_path, .allowed_count = 1 } });
            path = arg;
        }
    }
    const fasta = path orelse return cli.fail(.{ .@"usage.missing_argument" = .{ .argument = .fasta_path } });
    const result = try indexer.runIndex(std.heap.page_allocator, cli.io, environ, fasta, .{ .emit_fai = emit_fai, .enable_dedup = enable_dedup }, cli.failure);
    if (!emit_fai and cli.format == .text) try writeNotice(cli, .{ .index = .{ .path = fasta, .count = result.record_count } });
    return 0;
}

fn runGetCmd(cli: *Cli) diagnostics.Error!u8 {
    var fasta_path: ?[]const u8 = null;
    var bed_path: ?[]const u8 = null;
    var names_path: ?[]const u8 = null;
    var honor_strand = false;
    var summary = false;
    var rc = false;
    var complement_only = false;
    var reverse_only = false;
    var annotate_transform = false;
    var regions: [1024][]const u8 = undefined;
    var count: usize = 0;
    while (cli.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) return cli.stdout(USAGE);
        if (std.mem.eql(u8, arg, "--bed")) {
            if (bed_path != null) return cli.fail(.{ .@"usage.repeated_option" = .{ .option = .bed } });
            bed_path = try cli.value(.bed, .path);
        } else if (std.mem.eql(u8, arg, "--names")) {
            if (names_path != null) return cli.fail(.{ .@"usage.repeated_option" = .{ .option = .names } });
            names_path = try cli.value(.names, .path);
        } else if (std.mem.eql(u8, arg, STRAND_AWARE_FLAG) or std.mem.eql(u8, arg, STRAND_AWARE_ALIAS)) honor_strand = true else if (std.mem.eql(u8, arg, "--summary")) summary = true else if (std.mem.eql(u8, arg, RC_FLAG)) rc = true else if (std.mem.eql(u8, arg, COMPLEMENT_ONLY_FLAG)) complement_only = true else if (std.mem.eql(u8, arg, REVERSE_ONLY_FLAG)) reverse_only = true else if (std.mem.eql(u8, arg, ANNOTATE_RC_FLAG)) annotate_transform = true else {
            try cli.positional(arg);
            if (fasta_path == null) fasta_path = arg else {
                if (count == regions.len) return cli.fail(.{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = .positional_requests, .limit = regions.len, .required = .{ .exact = count + 1 } } });
                regions[count] = arg;
                count += 1;
            }
        }
    }
    const path = fasta_path orelse return cli.fail(.{ .@"usage.missing_argument" = .{ .argument = .fasta_path } });
    const source_count = @as(usize, @intFromBool(count != 0)) + @as(usize, @intFromBool(bed_path != null)) + @as(usize, @intFromBool(names_path != null));
    if (source_count == 0) return cli.fail(.{ .@"usage.missing_argument" = .{ .argument = .request_source } });
    if (source_count > 1) {
        var arguments = diagnostics.Arguments.initEmpty();
        if (count != 0) arguments.insert(.positional);
        if (bed_path != null) arguments.insert(.bed);
        if (names_path != null) arguments.insert(.names);
        return cli.fail(.{ .@"usage.conflicting_options" = .{ .arguments = arguments } });
    }
    if (honor_strand and bed_path == null) return cli.fail(.{ .@"usage.option_requires" = .{ .option = .strand_aware, .required = .bed } });
    const transform = parseTransformFlags(rc, complement_only, reverse_only, annotate_transform) catch {
        var arguments = diagnostics.Arguments.initEmpty();
        if (rc) arguments.insert(.rc);
        if (complement_only) arguments.insert(.complement_only);
        if (reverse_only) arguments.insert(.reverse_only);
        return cli.fail(.{ .@"usage.conflicting_options" = .{ .arguments = arguments } });
    };
    if (transform.annotate_transform and transform.orientation.isIdentity() and !honor_strand) return cli.fail(.{ .@"usage.option_requires" = .{ .option = .annotate_rc, .required = .transform_or_strand } });
    const source: getter.RequestSource = if (bed_path) |bed| .{ .bed = bed } else if (names_path) |names| .{ .names = names } else .{ .positional = regions[0..count] };
    const result = try getter.runGet(std.heap.page_allocator, cli.io, path, .{ .source = source, .honor_strand = honor_strand, .measure_elapsed = summary, .orientation = transform.orientation, .annotate_transform = transform.annotate_transform }, cli.failure);
    if (summary) try writeNotice(cli, .{ .get = result });
    return 0;
}

fn runStatsCmd(cli: *Cli) diagnostics.Error!u8 {
    var path: ?[]const u8 = null;
    while (cli.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) return cli.stdout(USAGE);
        try cli.positional(arg);
        if (path != null) return cli.fail(.{ .@"usage.extra_argument" = .{ .argument = .fasta_path, .allowed_count = 1 } });
        path = arg;
    }
    const fasta = path orelse return cli.fail(.{ .@"usage.missing_argument" = .{ .argument = .fasta_path } });
    try stats.runStats(std.heap.page_allocator, cli.io, fasta, cli.failure);
    return 0;
}

fn runValidateCmd(cli: *Cli) diagnostics.Error!u8 {
    var path: ?[]const u8 = null;
    var options = validator.Options{};
    var saw_json = false;
    var saw_summary = false;
    while (cli.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) return cli.stdout(USAGE);
        if (std.mem.eql(u8, arg, "--strict")) options.strict = true else if (std.mem.eql(u8, arg, "--json")) {
            saw_json = true;
            options.output_mode = if (saw_summary) .json_summary else .json_lines;
        } else if (std.mem.eql(u8, arg, "--summary")) {
            saw_summary = true;
            options.output_mode = .json_summary;
        } else if (std.mem.eql(u8, arg, "--fix")) options.fix = true else if (std.mem.eql(u8, arg, "--fix-format-only")) options.fix_format_only = true else if (std.mem.eql(u8, arg, "-o")) options.output_path = try cli.value(.output, .path) else if (std.mem.eql(u8, arg, "--custom-alphabet")) options.custom_alphabet = try cli.value(.custom_alphabet, .alphabet) else if (std.mem.eql(u8, arg, "--schema")) {
            const raw = try cli.value(.schema, .schema);
            options.schema = if (std.mem.eql(u8, raw, "uniprot")) .uniprot else if (std.mem.eql(u8, raw, "refseq")) .refseq else return cli.fail(.{ .@"usage.invalid_option_value" = .{ .option = .schema, .expected = .schema, .value = cli.failure.retain(raw, 32) } });
        } else if (std.mem.eql(u8, arg, "--max-header-len")) {
            const raw = try cli.value(.max_header_len, .positive_usize);
            const width = std.fmt.parseInt(usize, raw, 10) catch return invalidWidth(cli, raw);
            if (width == 0) return invalidWidth(cli, raw);
            options.max_header_len = width;
        } else {
            try cli.positional(arg);
            if (path != null) return cli.fail(.{ .@"usage.extra_argument" = .{ .argument = .fasta_path, .allowed_count = 1 } });
            path = arg;
        }
    }
    if (saw_summary and !saw_json) return cli.fail(.{ .@"usage.option_requires" = .{ .option = .summary, .required = .json } });
    if (options.fix_format_only and !options.fix) return cli.fail(.{ .@"usage.option_requires" = .{ .option = .fix_format_only, .required = .fix } });
    if (options.output_path != null and !options.fix) return cli.fail(.{ .@"usage.option_requires" = .{ .option = .output, .required = .fix } });
    const fasta = path orelse return cli.fail(.{ .@"usage.missing_argument" = .{ .argument = .fasta_path } });
    const result = try validator.runValidate(std.heap.page_allocator, cli.io, fasta, options, cli.failure);
    if (result.truncated and cli.format == .text) try writeNotice(cli, .{ .cap = result.fix_written });
    return result.status;
}

fn invalidWidth(cli: *Cli, raw: []const u8) diagnostics.Error {
    return cli.fail(.{ .@"usage.invalid_option_value" = .{ .option = .max_header_len, .expected = .positive_usize, .value = cli.failure.retain(raw, 32) } });
}

const Notice = union(enum) { index: struct { path: []const u8, count: u64 }, get: getter.Result, cap: bool };

fn writeNotice(cli: *Cli, notice: Notice) diagnostics.Error!void {
    const failure = cli.failure;
    if (failure.sources.get(.stderr) == null) failure.sources.set(.stderr, .{ .label = .{ .borrowed = "stderr" }, .kind = .stream });
    failure.outputs.stderr.state = .streaming;
    failure.outputs.stderr.may_be_partial = true;
    var buffer: [4096]u8 = undefined;
    var fw = std.Io.File.Writer.initStreaming(.stderr(), cli.io, &buffer);
    errdefer fw.interface.end = 0;
    writeNoticeBytes(&fw.interface, cli.format, failure, notice) catch |err| return failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .report, .at = .{ .source = .stderr }, .cause = .{ .zig = fw.err orelse err } });
    fw.flush() catch |err| {
        failure.outputs.stderr.final_flush = .failed;
        return failure.fail(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .stderr }, .cause = .{ .zig = fw.err orelse err } });
    };
    failure.outputs.stderr.final_flush = .succeeded;
}

fn writeNoticeBytes(writer: *std.Io.Writer, format: ErrorFormat, failure: *const diagnostics.Failure, notice: Notice) std.Io.Writer.Error!void {
    switch (notice) {
        .index => |value| {
            try writer.writeAll("wrote ");
            try reporting.writeEscapedBytes(writer, value.path);
            try writer.print(".zfi ({d} sequences)\n", .{value.count});
        },
        .get => |value| if (format == .json) try reporting.writeJson(writer, failure, .{ .version = VERSION, .accepted = true, .metrics = .{ .region_count = value.region_count, .total_bases = value.total_bases, .elapsed_ns = value.elapsed_ns } }) else try writeGetSummary(writer, value),
        .cap => |fixed| if (fixed) try writer.print("warning: validate event list truncated at {d}; --fix output was still written\n", .{validator.MAX_VALIDATE_EVENTS}) else try writer.print("error: validate retained {d} events; scanning and counting completed; fix reported issues and re-run\n", .{validator.MAX_VALIDATE_EVENTS}),
    }
}

fn writeGetSummary(writer: *std.Io.Writer, result: getter.Result) std.Io.Writer.Error!void {
    const seconds = if (result.elapsed_ns == 0) 0.0 else @as(f64, @floatFromInt(result.elapsed_ns)) / @as(f64, @floatFromInt(std.time.ns_per_s));
    const regions_per_second = if (seconds == 0.0) 0.0 else @as(f64, @floatFromInt(result.region_count)) / seconds;
    try writer.print("summary: regions={d} total_bases={d} elapsed_s={d:.6} regions_per_s={d:.1}\n", .{ result.region_count, result.total_bases, seconds, regions_per_second });
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
    const zfasta_bin = TEST_BINARY_PATH;
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
    const errors = [_][]const u8{
        "error [usage.command_required] cli: A command is required.\nhint: Correct the command arguments; use --help for accepted options.\n",
        "error [usage.missing_argument] index: A required argument is missing. argument=fasta_path\nhint: Correct the command arguments; use --help for accepted options.\n",
        "error [usage.unknown_command] cli: Unknown command. command=\"not-a-command\"\nhint: Correct the command arguments; use --help for accepted options.\n",
    };
    for (cases, errors) |argv, expected| {
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
        try std.testing.expectEqualStrings(expected, result.stderr);
    }
}

test "[failure] - [get summary]: propagates a full destination" {
    var buffer: [1]u8 = undefined;
    var writer = std.Io.Writer.fixed(&buffer);
    try std.testing.expectError(error.WriteFailed, writeGetSummary(&writer, .{ .region_count = 1, .total_bases = 1, .elapsed_ns = 1 }));
}

test "[cli] - [diagnostic format]: respects placement, values, repetition and help priority" {
    const binary = TEST_BINARY_PATH;
    var threaded = std.Io.Threaded.init(std.testing.allocator, .{});
    defer threaded.deinit();
    const Case = struct { argv: []const []const u8, code: diagnostics.Code, json: bool = true };
    const cases = [_]Case{
        .{ .argv = &.{ binary, "--error-format=json" }, .code = .@"usage.command_required" },
        .{ .argv = &.{ binary, "--error-format=json", "unknown" }, .code = .@"usage.unknown_command" },
        .{ .argv = &.{ binary, "stats", "--bad", "--error-format=json" }, .code = .@"usage.unknown_option" },
        .{ .argv = &.{ binary, "get", "--error-format=text", "--error-format=json" }, .code = .@"usage.missing_argument" },
        .{ .argv = &.{ binary, "get", "--error-format=json", "--error-format=text" }, .code = .@"usage.missing_argument", .json = false },
        .{ .argv = &.{ binary, "validate", "--custom-alphabet", "--error-format=json", "--schema", "bad" }, .code = .@"usage.invalid_option_value", .json = false },
        .{ .argv = &.{ binary, "get", "tests/data/simple.fasta", "--names", "--error-format=json" }, .code = .@"io.file_not_found", .json = false },
        .{ .argv = &.{ binary, "stats", "--error-format=json", "--error-format=bad" }, .code = .@"usage.invalid_option_value", .json = false },
        .{ .argv = &.{ binary, "stats", "--error-format=json", "--error-formatting" }, .code = .@"usage.unknown_option" },
    };
    for (cases) |case| {
        const result = try std.process.run(std.testing.allocator, threaded.io(), .{ .argv = case.argv });
        defer std.testing.allocator.free(result.stdout);
        defer std.testing.allocator.free(result.stderr);
        try std.testing.expectEqual(@as(u8, 1), result.term.exited);
        try std.testing.expectEqualStrings("", result.stdout);
        if (case.json) {
            const parsed = try std.json.parseFromSlice(std.json.Value, std.testing.allocator, result.stderr, .{});
            defer parsed.deinit();
            try std.testing.expectEqualStrings(@tagName(case.code), parsed.value.object.get("failure").?.object.get("primary").?.object.get("code").?.string);
            try std.testing.expect(!parsed.value.object.get("accepted").?.bool);
        } else {
            var prefix: [128]u8 = undefined;
            try std.testing.expect(std.mem.startsWith(u8, result.stderr, try std.fmt.bufPrint(&prefix, "error [{s}]", .{@tagName(case.code)})));
        }
    }
    const help = try std.process.run(std.testing.allocator, threaded.io(), .{ .argv = &.{ binary, "get", "--error-format=json", "--help", "--error-format=bad" } });
    defer std.testing.allocator.free(help.stdout);
    defer std.testing.allocator.free(help.stderr);
    try std.testing.expectEqual(@as(u8, 0), help.term.exited);
    try std.testing.expectEqualStrings(USAGE, help.stdout);
    try std.testing.expectEqualStrings("", help.stderr);
    for ([_][]const u8{ "--error-format=text", "--error-format=json" }) |format| {
        for ([_][]const u8{ "--help", "-h", "--version", "-V" }) |option| {
            const result = try std.process.run(std.testing.allocator, threaded.io(), .{
                .argv = &.{ binary, format, option, "--error-format=bad" },
            });
            defer std.testing.allocator.free(result.stdout);
            defer std.testing.allocator.free(result.stderr);
            try std.testing.expectEqual(@as(u8, 0), result.term.exited);
            const expected = if (std.mem.eql(u8, option, "--help") or std.mem.eql(u8, option, "-h")) USAGE else "z-fasta " ++ VERSION ++ "\n";
            try std.testing.expectEqualStrings(expected, result.stdout);
            try std.testing.expectEqualStrings("", result.stderr);
        }
    }
}

test "[failure] - [CLI delivery]: retains primary after failed stderr and observes stdout failures" {
    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => .{ .file_write_streaming = error.BrokenPipe },
                else => std.testing.io.vtable.operate(userdata, operation),
            };
        }
    };
    var vtable = std.testing.io.vtable.*;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = std.testing.io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    failure.reset(.cli);
    switch (failure.fail(.{ .details = .{ .@"usage.command_required" = {} }, .stage = .arguments })) {
        else => {},
    }
    try std.testing.expectEqual(@as(u8, 1), reportFailure(fault_io, .json, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"usage.command_required", failure.primary().?.code());
    try std.testing.expectEqual(@as(u8, 1), failure.secondary_len);
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.secondary[0].code());
    try std.testing.expectEqual(error.BrokenPipe, failure.secondary[0].cause.?.zig);
    try std.testing.expectEqual(.failed, failure.outputs.stderr.final_flush);
    failure.reset(.cli);
    var cli = Cli{ .io = fault_io, .args = undefined, .format = .json, .failure = &failure };
    try std.testing.expectError(error.Io, cli.stdout("help\n"));
    try std.testing.expectEqual(diagnostics.Code.@"io.write_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    failure.reset(.index);
    failure.succeed();
    try std.testing.expectError(error.Io, writeNotice(&cli, .{ .index = .{ .path = "input.fa", .count = 1 } }));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expect(failure.work_complete and !failure.complete);
}

test "[cli] - [Linux delivery]: help, version and required stderr failures return failure" {
    if (comptime builtin.os.tag != .linux) return error.SkipZigTest;
    const binary = TEST_BINARY_PATH;
    var threaded = std.Io.Threaded.init(std.testing.allocator, .{});
    defer threaded.deinit();
    for ([_][]const u8{ "--help", "--version" }) |option| {
        for ([_]bool{ false, true }) |json| {
            const script = if (json) "exec \"$1\" --error-format=json \"$2\" --error-format=bad >/dev/full" else "exec \"$1\" \"$2\" --error-format=bad >/dev/full";
            const result = try std.process.run(std.testing.allocator, threaded.io(), .{ .argv = &.{ "/bin/sh", "-c", script, "sh", binary, option } });
            defer std.testing.allocator.free(result.stdout);
            defer std.testing.allocator.free(result.stderr);
            try std.testing.expectEqual(@as(u8, 1), result.term.exited);
            try std.testing.expectEqualStrings("", result.stdout);
            if (json) {
                const parsed = try std.json.parseFromSlice(std.json.Value, std.testing.allocator, result.stderr, .{});
                defer parsed.deinit();
                const primary = parsed.value.object.get("failure").?.object.get("primary").?.object;
                try std.testing.expectEqualStrings("io.write_failed", primary.get("code").?.string);
                try std.testing.expectEqualStrings("stdout", primary.get("source").?.string);
            } else {
                try std.testing.expect(std.mem.startsWith(u8, result.stderr, "error [io.write_failed]"));
                try std.testing.expect(std.mem.indexOf(u8, result.stderr, "cause=posix_errno:28") != null);
            }
        }
    }
    const result = try std.process.run(std.testing.allocator, threaded.io(), .{ .argv = &.{ "/bin/sh", "-c", "exec \"$1\" get tests/data/simple.fasta seq1:1-4 --summary 2>/dev/full", "sh", binary } });
    defer std.testing.allocator.free(result.stdout);
    defer std.testing.allocator.free(result.stderr);
    try std.testing.expectEqual(@as(u8, 1), result.term.exited);
    try std.testing.expectEqualStrings(">seq1:1-4\nACGT\n", result.stdout);
    try std.testing.expectEqualStrings("", result.stderr);
}
