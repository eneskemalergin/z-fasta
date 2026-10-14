//! Validates FASTA for `z-fasta validate` in one streaming pass through a fixed buffer; `--fix`
//! writes the repaired copy during the same pass. Rules follow `fasta.zig`: a symbol is a byte
//! above space, an identifier ends at space, tab, CR, or LF, and a line ends at LF or CRLF.
//! Each finding code keeps an exact count, its affected records, and its first occurrence only.
//! Memory is the buffer, one identifier, and the duplicate-name set: on a file, a fingerprint and
//! an offset per name. The alphabet comes from the first 100,000 symbols unless set. Equally
//! wrapped lines are checked in vector blocks; the first irregular byte drops to the per-line path.
//! `--fix` wraps symbols at a fixed width, refuses content it cannot keep, and publishes through
//! atomic replacement only when the source is unchanged. `--index` runs `index`'s scanner on the
//! same windows and compares each record it would write with the loaded index.

const std = @import("std");

const fasta = @import("fasta.zig");
const diagnostics = @import("diagnostics");
const index_format = @import("index_format.zig");
const reporting = @import("reporting.zig");

pub const TYPE_SAMPLE_SYMBOLS: u64 = 100_000;
pub const READ_BUFFER_BYTES = 256 * 1024;
pub const DEFAULT_MAX_HEADER_BYTES: u64 = 1024;
pub const DEFAULT_FIX_WIDTH: u64 = 60;

pub const Schema = enum { none, uniprot, refseq };

/// `auto` chooses nucleotide or protein from the first 100,000 symbols; `any` accepts every symbol.
pub const Alphabet = enum { auto, dna, rna, protein, any };

/// What the scan checks. `custom_alphabet` is an exact, case-sensitive byte set and cannot be
/// combined with an `alphabet` other than `auto`.
pub const Rules = struct {
    schema: Schema = .none,
    alphabet: Alphabet = .auto,
    custom_alphabet: ?[]const u8 = null,
    max_header_bytes: u64 = DEFAULT_MAX_HEADER_BYTES,
    /// Bytes the duplicate check may hold; past it the check stops and the report is `incomplete`.
    max_name_bytes: ?u64 = null,
};

/// Not validated here: `output_path` and `width` apply only with `fix`, and `width` is at least 1.
pub const Options = struct {
    rules: Rules = .{},
    strict: bool = false,
    json: bool = false,
    fix: bool = false,
    output_path: ?[]const u8 = null,
    width: u64 = DEFAULT_FIX_WIDTH,
    /// The index to compare the FASTA with: `.automatic` is the one `get` would load. Needs a
    /// file, not stdin.
    index: ?index_format.IndexSelection = null,
    /// Producer version in a JSON report: trusted UTF-8 of at most 32 bytes.
    producer_version: ?[]const u8 = null,
};

/// Where the duplicate check re-reads names from; without one, names are copied into memory.
pub const NameSource = struct { io: std.Io, file: std.Io.File, size: u64 };

pub const Result = struct { status: u8, fix_written: bool };

/// One input's counts and first occurrences; `validation` returns a view that borrows it.
pub const Report = struct {
    completion: reporting.Validation.Completion = .interrupted,
    findings: [REPORT_CODES.len]reporting.Finding = blk: {
        var findings: [REPORT_CODES.len]reporting.Finding = undefined;
        for (REPORT_CODES, &findings) |code, *entry| entry.* = .{ .code = code };
        break :blk findings;
    },
    examined: reporting.Validation.Examined = .{},
    alphabet: diagnostics.AlphabetPolicy = .{ .alphabet = .nucleotide, .selection = .automatic },
    records_with_errors: u64 = 0,
    records_with_warnings: u64 = 0,
    invalid_symbols: [256]u64 = @splat(0),
    name_budget_stop: ?reporting.Validation.BudgetStop = null,
    index: ?reporting.Validation.IndexComparison = null,

    pub fn validation(self: *const Report, rules: Rules, strict: bool) reporting.Validation {
        return .{
            .completion = self.completion,
            .strict = strict,
            .alphabet = self.alphabet,
            .max_header_bytes = rules.max_header_bytes,
            .convention = convention(rules.schema),
            .examined = self.examined,
            .records_with_errors = self.records_with_errors,
            .records_with_warnings = self.records_with_warnings,
            .findings = &self.findings,
            .invalid_symbols = &self.invalid_symbols,
            .name_budget_stop = self.name_budget_stop,
            .index = self.index,
        };
    }

    pub fn finding(self: *const Report, comptime code: diagnostics.Code) *const reporting.Finding {
        return &self.findings[comptime reportSlot(code)];
    }

    /// 0 clean, 2 warnings only, 1 for errors, warnings under `strict`, or an unfinished report.
    /// `repaired` leaves out the format warnings a published `--fix` copy no longer has.
    pub fn status(self: *const Report, strict: bool, repaired: bool) u8 {
        if (self.completion != .complete) return 1;
        var warnings = false;
        for (self.findings) |entry| {
            if (entry.occurrences == 0) continue;
            const metadata = diagnostics.metadata(entry.code);
            if (metadata.severity == .@"error") return 1;
            if (!(repaired and metadata.category == .format)) warnings = true;
        }
        if (warnings) return if (strict) 1 else 2;
        return 0;
    }
};

pub fn runValidate(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!Result {
    failure.reset(.validate);
    const from_stdin = std.mem.eql(u8, fasta_path, "-");
    failure.sources.set(.fasta, if (from_stdin) .{ .label = .{ .borrowed = "stdin" }, .kind = .stream } else .{ .label = .{ .borrowed = fasta_path }, .kind = .path });
    failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
    if (options.output_path) |path| failure.sources.set(.output, .{ .label = .{ .borrowed = path }, .kind = .path });
    const file = if (from_stdin) std.Io.File.stdin() else std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| return failure.fail(.{ .details = diagnostics.openFailure(err), .stage = .open, .at = .{ .source = .fasta }, .cause = .{ .zig = err } });
    defer if (!from_stdin) file.close(io);
    const before: ?std.Io.File.Stat = if (from_stdin) null else file.stat(io) catch |err| return failure.fail(.{ .details = .{ .@"io.stat_failed" = {} }, .stage = .stat, .at = .{ .source = .fasta }, .cause = .{ .zig = err } });

    var pair: Pair = undefined;
    if (options.index) |selection| {
        if (from_stdin) return failure.fail(.{ .details = .{ .@"usage.option_requires" = .{ .option = .index, .required = .fasta_file } }, .stage = .arguments, .at = .{ .source = .fasta } });
        try pair.open(allocator, io, fasta_path, selection, failure);
    }
    defer if (options.index != null) pair.deinit();

    var staged: Staged = undefined;
    var staged_open = false;
    defer if (staged_open) staged.cleanup(io, failure);
    var rewrite: ?Rewrite = null;
    if (options.fix) {
        const output_path = options.output_path orelse return failure.fail(.{ .details = .{ .@"usage.option_requires" = .{ .option = .fix, .required = .output } }, .stage = .arguments, .at = .{ .source = .fasta } });
        if (!from_stdin) try checkRewriteDestination(io, fasta_path, output_path, failure);
        try staged.create(allocator, io, output_path, failure);
        staged_open = true;
        rewrite = .{ .writer = &staged.writer.interface, .width = options.width, .file = &staged.writer };
    }

    const report = allocator.create(Report) catch return failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .validation_report, .requested_bytes = @sizeOf(Report) } }, .stage = .allocate, .at = .{ .source = .fasta }, .cause = .{ .zig = error.OutOfMemory } });
    defer allocator.destroy(report);
    const read_buf = allocator.alloc(u8, READ_BUFFER_BYTES) catch return failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .reader_buffer, .requested_bytes = READ_BUFFER_BYTES } }, .stage = .allocate, .at = .{ .source = .fasta }, .cause = .{ .zig = error.OutOfMemory } });
    defer allocator.free(read_buf);
    var io_buf: [FILE_IO_BUFFER_BYTES]u8 = undefined;
    var reader = if (from_stdin) file.readerStreaming(io, &io_buf) else file.reader(io, &io_buf);
    const names: ?NameSource = if (before) |stat| .{ .io = io, .file = file, .size = stat.size } else null;
    const compared: ?*Pair = if (options.index != null) &pair else null;
    const scanned = if (rewrite) |*sink|
        repairInto(allocator, &reader.interface, read_buf, names, options.rules, sink, compared, report, failure)
    else
        scan(allocator, &reader.interface, read_buf, names, options.rules, null, compared, report, failure);
    scanned catch |err| {
        diagnostics.retainReaderCause(failure, &reader);
        writeReportOutput(io, report, options, false, failure) catch {};
        return err;
    };
    if (options.fix) try staged.publish(io, if (before) |stat| .{ .file = file, .before = stat } else null, failure);
    failure.succeed();
    try writeReportOutput(io, report, options, options.fix, failure);
    return .{ .status = report.status(options.strict, options.fix), .fix_written = options.fix };
}

/// Fills `report` from `source`, read through `buffer`. After a failure, `report` covers what was
/// examined and is `interrupted`.
pub fn validate(allocator: std.mem.Allocator, source: *std.Io.Reader, buffer: []u8, names: ?NameSource, rules: Rules, report: *Report, failure: *diagnostics.Failure) diagnostics.Error!void {
    failure.reset(.validate);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "validated input" }, .kind = .stream });
    try scan(allocator, source, buffer, names, rules, null, null, report, failure);
    failure.succeed();
}

/// Validates as `validate` does and writes the repaired copy to `output` in the same pass:
/// headers as they are, symbols wrapped at `width`, LF line ends. On any error, including
/// `repair.blocked`, the caller discards what `output` received.
pub fn repair(allocator: std.mem.Allocator, source: *std.Io.Reader, buffer: []u8, names: ?NameSource, rules: Rules, width: u64, output: *std.Io.Writer, report: *Report, failure: *diagnostics.Failure) diagnostics.Error!void {
    failure.reset(.validate);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "validated input" }, .kind = .stream });
    failure.sources.set(.output, .{ .label = .{ .borrowed = "repaired output" }, .kind = .stream });
    if (width == 0) return failure.fail(.{ .details = .{ .@"usage.invalid_option_value" = .{ .option = .width, .expected = .positive_usize } }, .stage = .arguments });
    var sink = Rewrite{ .writer = output, .width = width };
    try repairInto(allocator, source, buffer, names, rules, &sink, null, report, failure);
    failure.succeed();
}

const OUTPUT_BUFFER_BYTES = 64 * 1024;
// A report is a few KiB; a longer one takes more than one write.
const REPORT_BUFFER_BYTES = 16 * 1024;
const FILE_IO_BUFFER_BYTES = 8 * 1024;
const MAX_LINK_HOPS = 32;
const UTF8_BOM = "\xEF\xBB\xBF";

const REPORT_CODES = [_]diagnostics.Code{
    .@"fasta.no_records",
    .@"fasta.content_before_header",
    .@"fasta.header_control_byte",
    .@"fasta.sequence_control_byte",
    .@"fasta.empty_record",
    .@"fasta.empty_identifier",
    .@"fasta.duplicate_identifier",
    .@"fasta.identifier_too_long",
    .@"fasta.long_header",
    .@"header.convention_mismatch",
    .@"sequence.invalid_symbol",
    .@"format.sequence_whitespace",
    .@"format.blank_sequence_line",
    .@"format.inconsistent_wrapping",
    .@"format.leading_bom",
    .@"format.missing_terminal_newline",
    .@"format.mixed_line_endings",
    .@"pair.record_missing_from_index",
    .@"pair.record_missing_from_fasta",
    .@"pair.sequence_length_mismatch",
    .@"pair.sequence_offset_mismatch",
    .@"pair.line_geometry_mismatch",
    .@"pair.side_table_mismatch",
};

const REPORT_SLOTS = blk: {
    var slots: [std.enums.values(diagnostics.Code).len]?u8 = @splat(null);
    for (REPORT_CODES, 0..) |code, slot| slots[@intFromEnum(code)] = slot;
    break :blk slots;
};

fn reportSlot(comptime code: diagnostics.Code) usize {
    return REPORT_SLOTS[@intFromEnum(code)] orelse @compileError("not a validate code: " ++ @tagName(code));
}

// Bit `n` of a byte's class marks a symbol outside candidate `n` (nucleotide, protein, or custom twice).
const WHITESPACE: u8 = 1 << 2;
const CONTROL: u8 = 1 << 3;

fn outsideBit(candidate: u1) u8 {
    return @as(u8, 1) << candidate;
}

fn convention(schema: Schema) ?diagnostics.Convention {
    return switch (schema) {
        .none => null,
        .uniprot => .uniprot,
        .refseq => .refseq,
    };
}

fn alphabetPolicy(rules: Rules, candidate: u1, sampled: u64) diagnostics.AlphabetPolicy {
    if (rules.custom_alphabet) |custom| return .{ .alphabet = .custom, .selection = .custom, .custom_symbols = customSymbols(custom) };
    switch (rules.alphabet) {
        .auto => {},
        .dna => return .{ .alphabet = .dna, .selection = .explicit },
        .rna => return .{ .alphabet = .rna, .selection = .explicit },
        .protein => return .{ .alphabet = .protein, .selection = .explicit },
        .any => return .{ .alphabet = .any, .selection = .explicit },
    }
    return .{ .alphabet = if (candidate == 0) .nucleotide else .protein, .selection = .automatic, .sample = .{ .symbols = @intCast(sampled), .rule = .nucleotide_above_90_percent } };
}

fn buildAlphabet(chars: []const u8) [256]bool {
    var table = [_]bool{false} ** 256;
    for (chars) |byte| table[byte] = true;
    return table;
}

fn customSymbols(custom: []const u8) [32]u8 {
    var bits: [32]u8 = @splat(0);
    for (custom) |byte| bits[byte / 8] |= @as(u8, 1) << @intCast(byte % 8);
    return bits;
}

// `repaired`: a `--fix` copy was published, as `Report.status` takes it.
fn writeReportOutput(io: std.Io, report: *const Report, options: Options, repaired: bool, failure: *diagnostics.Failure) diagnostics.Error!void {
    var buffer: [REPORT_BUFFER_BYTES]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &buffer);
    errdefer stdout_fw.interface.end = 0;
    failure.outputs.stdout.state = .streaming;
    failure.outputs.stdout.may_be_partial = true;
    const view = report.validation(options.rules, options.strict);
    const written = if (options.json)
        reporting.writeJson(&stdout_fw.interface, failure, .{ .version = options.producer_version, .accepted = report.status(options.strict, repaired) != 1, .validation = &view })
    else
        reporting.writeValidationText(&stdout_fw.interface, failure, &view);
    written catch |err| return failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .report, .at = .{ .source = .stdout }, .cause = .{ .zig = stdout_fw.err orelse err } });
    if (diagnostics.flushOutput(&stdout_fw, &failure.outputs.stdout)) |cause| return failure.fail(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .stdout }, .cause = .{ .zig = cause } });
}

const Staged = struct {
    atomic: std.Io.File.Atomic,
    writer: std.Io.File.Writer,
    allocator: std.mem.Allocator,
    buffer: []u8,

    // The buffer is allocated first, so running out of memory leaves no staged file.
    fn create(self: *Staged, allocator: std.mem.Allocator, io: std.Io, path: []const u8, failure: *diagnostics.Failure) diagnostics.Error!void {
        self.allocator = allocator;
        self.buffer = allocator.alloc(u8, OUTPUT_BUFFER_BYTES) catch return failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .output_buffer, .requested_bytes = OUTPUT_BUFFER_BYTES } }, .stage = .allocate, .at = .{ .source = .output }, .cause = .{ .zig = error.OutOfMemory } });
        errdefer allocator.free(self.buffer);
        self.atomic = std.Io.Dir.cwd().createFileAtomic(io, path, .{ .replace = true }) catch |err| return failure.fail(.{ .details = .{ .@"io.create_failed" = {} }, .stage = .create, .at = .{ .source = .output }, .cause = .{ .zig = err } });
        failure.outputs.file.state = .staged;
        self.writer = self.atomic.file.writer(io, self.buffer);
    }

    // `source` is null for stdin, which has no metadata to compare.
    fn publish(self: *Staged, io: std.Io, source: ?struct { file: std.Io.File, before: std.Io.File.Stat }, failure: *diagnostics.Failure) diagnostics.Error!void {
        if (diagnostics.flushOutput(&self.writer, &failure.outputs.file)) |cause| return failure.fail(.{ .details = .{ .@"io.flush_failed" = {} }, .stage = .flush, .at = .{ .source = .output }, .cause = .{ .zig = cause } });
        if (source) |input| try checkRewriteSource(io, input.file, input.before, failure);
        self.atomic.replace(io) catch |err| {
            failure.outputs.file.state = .unknown;
            return failure.fail(.{ .details = .{ .@"io.publish_failed" = {} }, .stage = .publish, .at = .{ .source = .output }, .cause = .{ .zig = err } });
        };
        failure.outputs.file.state = .published;
    }

    fn cleanup(self: *Staged, io: std.Io, failure: *diagnostics.Failure) void {
        self.writer.interface.end = 0;
        cleanupFixed(&self.atomic, io, failure);
        self.allocator.free(self.buffer);
    }
};

const FixBlocker = struct { blocker: diagnostics.Blocker, at: diagnostics.Reference };

// Duplicates and symbols outside the alphabet do not block: the rewrite copies them.
fn fixBlocker(report: *const Report) ?FixBlocker {
    if (report.finding(.@"fasta.no_records").occurrences != 0) return .{ .blocker = .{ .no_records = {} }, .at = .{ .source = .fasta } };
    var best: ?FixBlocker = null;
    if (report.finding(.@"fasta.content_before_header").first) |first| keepEarliest(&best, .{ .content_before_header = first.details.@"fasta.content_before_header" }, first.at);
    if (report.finding(.@"fasta.header_control_byte").first) |first| keepEarliest(&best, .{ .header_control_byte = first.details.@"fasta.header_control_byte" }, first.at);
    if (report.finding(.@"fasta.sequence_control_byte").first) |first| keepEarliest(&best, .{ .sequence_control_byte = first.details.@"fasta.sequence_control_byte" }, first.at);
    return best;
}

fn keepEarliest(best: *?FixBlocker, blocker: diagnostics.Blocker, at: diagnostics.Reference) void {
    if (best.*) |current| {
        if (current.at.location.offset.? <= at.location.offset.?) return;
    }
    best.* = .{ .blocker = blocker, .at = at };
}

const Rewrite = struct {
    writer: *std.Io.Writer,
    width: u64,
    // Its `err` holds the cause when a write through `writer` fails.
    file: ?*std.Io.File.Writer = null,
    column: u64 = 0,

    fn startHeader(self: *Rewrite) std.Io.Writer.Error!void {
        if (self.column != 0) try self.put("\n");
        self.column = 0;
        try self.put(">");
    }

    fn header(self: *Rewrite, bytes: []const u8) std.Io.Writer.Error!void {
        try self.put(bytes);
    }

    fn endHeader(self: *Rewrite) std.Io.Writer.Error!void {
        try self.put("\n");
    }

    fn symbols(self: *Rewrite, bytes: []const u8) std.Io.Writer.Error!void {
        var rest = bytes;
        while (rest.len != 0) {
            if (self.column == self.width) {
                try self.put("\n");
                self.column = 0;
            }
            const take: usize = @intCast(@min(@as(u64, rest.len), self.width - self.column));
            try self.put(rest[0..take]);
            self.column += take;
            rest = rest[take..];
        }
    }

    // `count` lines of `width` symbols, each followed by its ending, `line_bytes` apart. Lines already
    // at the output width with LF endings are copied as they are.
    fn lines(self: *Rewrite, data: []const u8, width: usize, line_bytes: usize, count: usize) std.Io.Writer.Error!void {
        if (width != self.width or line_bytes != width + 1 or (self.column != 0 and self.column != width)) {
            for (0..count) |k| try self.symbols(data[k * line_bytes ..][0..width]);
            return;
        }
        if (self.column == width) try self.put("\n");
        try self.put(data[0 .. count * line_bytes - 1]);
        self.column = width;
    }

    fn finish(self: *Rewrite) std.Io.Writer.Error!void {
        if (self.column != 0) try self.put("\n");
        self.column = 0;
    }

    // Fills the buffer before each flush, so the file is written in whole buffers: a write that
    // ends inside a page costs a filesystem about twice as much as one that ends on a boundary.
    fn put(self: *Rewrite, bytes: []const u8) std.Io.Writer.Error!void {
        var rest = bytes;
        while (rest.len != 0) {
            var room = self.writer.unusedCapacitySlice();
            if (room.len == 0) {
                try self.writer.flush();
                room = self.writer.unusedCapacitySlice();
                if (room.len == 0) return self.writer.writeAll(rest);
            }
            const take = @min(room.len, rest.len);
            @memcpy(room[0..take], rest[0..take]);
            self.writer.advance(take);
            rest = rest[take..];
        }
    }
};

fn repairInto(allocator: std.mem.Allocator, source: *std.Io.Reader, buffer: []u8, names: ?NameSource, rules: Rules, sink: *Rewrite, pair: ?*Pair, report: *Report, failure: *diagnostics.Failure) diagnostics.Error!void {
    try scan(allocator, source, buffer, names, rules, sink, pair, report, failure);
    if (fixBlocker(report)) |blocker| return failure.fail(.{ .details = .{ .@"repair.blocked" = .{ .blocker = blocker.blocker, .repair_policy = .preserve_symbols } }, .stage = .scan, .at = blocker.at });
    sink.finish() catch |err| return failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .at = .{ .source = .output }, .cause = .{ .zig = if (sink.file) |file| file.err orelse err else err } });
}

fn scan(allocator: std.mem.Allocator, source: *std.Io.Reader, buffer: []u8, names: ?NameSource, rules: Rules, rewrite: ?*Rewrite, pair: ?*Pair, report: *Report, failure: *diagnostics.Failure) diagnostics.Error!void {
    report.* = .{ .alphabet = alphabetPolicy(rules, 0, 0) };
    if (buffer.len == 0) return failure.fail(.{ .details = .{ .@"usage.invalid_option_value" = .{ .option = .read_buffer, .expected = .nonempty_buffer } }, .stage = .scan, .at = .{ .source = .fasta } });
    if (rules.custom_alphabet != null and rules.alphabet != .auto) {
        var arguments = diagnostics.Arguments.initEmpty();
        arguments.insert(.custom_alphabet);
        arguments.insert(.sequence_type);
        return failure.fail(.{ .details = .{ .@"usage.conflicting_options" = .{ .arguments = arguments } }, .stage = .arguments });
    }
    const state = allocator.create(Scan) catch return failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .validation_report, .requested_bytes = @sizeOf(Scan) } }, .stage = .allocate, .at = .{ .source = .fasta }, .cause = .{ .zig = error.OutOfMemory } });
    defer allocator.destroy(state);
    state.* = .init(allocator, rules, report, rewrite, failure, names);
    state.pair = pair;
    state.startNames();
    defer state.deinit();
    while (true) {
        const n = source.readSliceShort(buffer) catch |err| {
            state.finish(.interrupted);
            return failure.fail(.{ .details = diagnostics.readFailure(buffer.len, null), .stage = .read, .at = .{ .source = .fasta, .location = .{ .offset = state.offset } }, .cause = .{ .zig = err } });
        };
        if (n == 0) break;
        state.feed(buffer[0..n]) catch |err| {
            state.finish(.interrupted);
            return err;
        };
        if (n < buffer.len) break;
    }
    state.end() catch |err| {
        state.finish(.interrupted);
        return err;
    };
    state.finish(.complete);
}

const LineEnding = enum { lf, crlf };

// Where a line or header starts. Plain numbers, so the per-line state is copied as words; `Scan.located`
// makes the reference when a finding needs one.
const Position = struct { offset: u64, line: u64, column: u64 };

const SymbolLine = struct {
    symbols: u64,
    start: Position,
    // Every byte was in the fast set, so the next lines may be checked in blocks.
    dense: bool,
    ending: ?LineEnding,
};

const Candidate = struct {
    records: u64 = 0,
    records_with_errors: u64 = 0,
    first: ?reporting.Occurrence = null,
};

const Record = struct {
    ordinal: u64 = 0,
    // `Scan.name_buf` holds the identifier, up to `MAX_INDEX_NAME_BYTES`, until the next header.
    named: bool = false,
    identifier_bytes: u64 = 0,
    header: Position = .{ .offset = 0, .line = 0, .column = 0 },
    symbols: u64 = 0,
    width: ?u64 = null,
    last_line: ?SymbolLine = null,
    wrapping_reported: bool = false,
    classes: u8 = 0,
    has_error: bool = false,
    has_warning: bool = false,
    duplicate: bool = false,
};

const Line = struct {
    kind: enum { start, preamble, header, sequence } = .start,
    start: u64 = 0,
    len: u64 = 0,
    content_at: ?u64 = null,
    in_identifier: bool = false,
    /// The whole identifier's length; `Scan.name_buf` keeps at most `MAX_INDEX_NAME_BYTES` of it.
    identifier_bytes: u64 = 0,
    controls: u64 = 0,
    first_control: ?struct { byte: u8, index: u64 } = null,
    symbols: u64 = 0,
    whitespace: u64 = 0,
    first_whitespace: u64 = 0,
    dense: bool = true,
};

const FastSet = enum { none, plain_nucleotides, plain_dna, plain_rna, protein, any };

const ScanDiag = struct {
    failure: ?*diagnostics.Failure = null,
    at: diagnostics.Reference = .{ .source = .fasta },
    budget: ?*const NameBudget = null,
    format: diagnostics.Format = .zfi,

    pub fn fail(self: @This(), details: diagnostics.Details, stage: diagnostics.Stage, source: diagnostics.SourceId, cause: ?anyerror, legacy: anyerror) anyerror {
        @branchHint(.cold);
        const failure = self.failure orelse return legacy;
        return failure.fail(.{ .details = details, .stage = stage, .at = if (self.at.source == source) self.at else .{ .source = source }, .cause = if (cause) |err| .{ .zig = err } else null });
    }

    pub fn oom(self: @This(), role: diagnostics.AllocationRole, bytes: ?u64) error{OutOfMemory} {
        @branchHint(.cold);
        if (self.budget) |budget| {
            if (budget.exhausted) return error.OutOfMemory;
        }
        if (self.failure) |failure| switch (failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = role, .requested_bytes = bytes } }, .stage = .allocate, .at = self.at, .cause = .{ .zig = error.OutOfMemory } })) {
            else => {},
        };
        return error.OutOfMemory;
    }

    pub fn capacity(self: @This(), role: diagnostics.CapacityRole, limit: u128, required: u128, legacy: anyerror) anyerror {
        return self.fail(.{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = role, .limit = limit, .required = .{ .exact = required } } }, .allocate, .fasta, null, legacy);
    }

    pub fn consistency(self: @This(), check: diagnostics.Check, legacy: anyerror) anyerror {
        return self.fail(.{ .@"internal.consistency_failed" = .{ .component = .validator, .check = check } }, .scan, .fasta, null, legacy);
    }

    pub fn recordLimit(self: @This(), count: u64) anyerror {
        return self.capacity(.records, std.math.maxInt(u32), count, error.OutOfMemory);
    }

    pub fn lineLimit(self: @This(), field: diagnostics.IndexField, observed: u64, width_role: diagnostics.LineWidthRole) anyerror {
        return self.fail(.{ .@"index.sequence_line_too_long" = .{ .field = field, .observed = .{ .exact = observed }, .limit = std.math.maxInt(u32), .width_role = width_role, .format = self.format } }, .scan, .fasta, null, error.SequenceLineTooLong);
    }

    pub fn identifierTooLong(self: @This(), observed: u128, limit: usize) anyerror {
        return self.fail(.{ .@"index.identifier_too_long" = .{ .observed_bytes = .{ .at_least = observed }, .limit_bytes = limit, .format = self.format } }, .scan, .fasta, null, error.HeaderTooLong);
    }
};

// Refuses an allocation that would take the duplicate check past `limit` live bytes.
const NameBudget = struct {
    child: std.mem.Allocator,
    limit: u64,
    live: u64 = 0,
    exhausted: bool = false,

    fn allocator(self: *NameBudget) std.mem.Allocator {
        return .{ .ptr = self, .vtable = &.{ .alloc = alloc, .resize = resize, .remap = remap, .free = free } };
    }

    fn admit(self: *NameBudget, more: usize) bool {
        if (self.live + more <= self.limit) return true;
        self.exhausted = true;
        return false;
    }

    fn alloc(context: *anyopaque, len: usize, alignment: std.mem.Alignment, ret_addr: usize) ?[*]u8 {
        const self: *NameBudget = @ptrCast(@alignCast(context));
        if (!self.admit(len)) return null;
        const memory = self.child.rawAlloc(len, alignment, ret_addr) orelse return null;
        self.live += len;
        return memory;
    }

    fn resize(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ret_addr: usize) bool {
        const self: *NameBudget = @ptrCast(@alignCast(context));
        if (new_len > memory.len and !self.admit(new_len - memory.len)) return false;
        if (!self.child.rawResize(memory, alignment, new_len, ret_addr)) return false;
        self.live = self.live - memory.len + new_len;
        return true;
    }

    fn remap(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ret_addr: usize) ?[*]u8 {
        const self: *NameBudget = @ptrCast(@alignCast(context));
        if (new_len > memory.len and !self.admit(new_len - memory.len)) return null;
        const moved = self.child.rawRemap(memory, alignment, new_len, ret_addr) orelse return null;
        self.live = self.live - memory.len + new_len;
        return moved;
    }

    fn free(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, ret_addr: usize) void {
        const self: *NameBudget = @ptrCast(@alignCast(context));
        self.child.rawFree(memory, alignment, ret_addr);
        self.live -= memory.len;
    }
};

const SourceNames = fasta.SourceDedupWith(ScanDiag);
const PackedNames = fasta.PackedNamesWith(ScanDiag);

const MAX_INDEX_NAME_BYTES = std.math.maxInt(@FieldType(index_format.IndexRecord, "name_len"));
comptime {
    std.debug.assert(MAX_INDEX_NAME_BYTES <= fasta.MAX_NAME_BYTES);
}
const NO_ENTRY = std.math.maxInt(u32);

// The `--index` comparison: `index`'s scanner runs on the same windows, and each record it would
// write is matched with an index entry. Entries stream in file order through fixed buffers while
// names agree; the first record out of order loads the index's catalog and a name lookup over it.
const Pair = struct {
    allocator: std.mem.Allocator,
    io: std.Io,
    failure: *diagnostics.Failure,
    fasta_path: []const u8,
    selection: index_format.IndexSelection,
    stream: *index_format.EntryStream,
    // The next entry in file order, read but not yet taken.
    head: ?*const index_format.EntryStream.Entry = null,
    catalog: ?index_format.LoadedIndex = null,
    lookup: ?Lookup = null,
    scanner: fasta.RecordScannerWith(ScanDiag),
    fed: u64 = 0,
    taken_in_order: usize = 0,
    matched: u64 = 0,
    omitted: u64 = 0,

    // `heads` maps a name to its first entry not yet taken; `next` chains entries by name.
    const Lookup = struct {
        heads: std.StringHashMapUnmanaged(u32),
        next: []u32,
        taken: std.DynamicBitSetUnmanaged,
    };

    const Association = union(enum) { entry: index_format.IndexRecord, omitted, missing };

    // `self` must not move afterwards: the scanner's diagnostics and the stream stay in place.
    fn open(self: *Pair, allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, selection: index_format.IndexSelection, failure: *diagnostics.Failure) diagnostics.Error!void {
        const stream = allocator.create(index_format.EntryStream) catch return failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .reader_buffer, .requested_bytes = @sizeOf(index_format.EntryStream) } }, .stage = .allocate, .cause = .{ .zig = error.OutOfMemory } });
        errdefer allocator.destroy(stream);
        try stream.open(allocator, io, fasta_path, selection, failure);
        const format: diagnostics.Format = if (stream.source == .zfi) .zfi else .fai;
        self.* = .{
            .allocator = allocator,
            .io = io,
            .failure = failure,
            .fasta_path = fasta_path,
            .selection = selection,
            .stream = stream,
            .scanner = .{ .allocator = allocator, .collect_side_tables = format == .zfi, .diag = .{ .failure = failure, .format = format } },
        };
    }

    fn deinit(self: *Pair) void {
        self.stream.close();
        self.allocator.destroy(self.stream);
        self.scanner.name.deinit(self.allocator);
        self.scanner.side_table.deinit(self.allocator);
        if (self.lookup) |*lookup| {
            lookup.heads.deinit(self.allocator);
            self.allocator.free(lookup.next);
            lookup.taken.deinit(self.allocator);
        }
        if (self.catalog) |*catalog| catalog.deinit();
    }

    // Entries read so far: every entry once the scan completes, fewer when it stops early.
    fn entriesRead(self: *const Pair) u64 {
        return if (self.catalog) |*catalog| catalog.records.len() else self.stream.ordinal;
    }

    fn source(self: *const Pair) diagnostics.SourceId {
        return if (self.stream.source == .zfi) .zfi else .fai;
    }

    // A later duplicate the index left out is `omitted`; so is any name whose entries are taken.
    fn associate(self: *Pair, name: []const u8, duplicate: bool) diagnostics.Error!Association {
        if (self.lookup == null) {
            if (self.head == null) self.head = try self.stream.next();
            if (self.head) |head| {
                if (std.mem.eql(u8, head.name, name)) {
                    self.head = null;
                    self.taken_in_order += 1;
                    return .{ .entry = head.record };
                }
            }
            if (duplicate) return .omitted;
            try self.buildLookup();
        }
        const catalog = &self.catalog.?;
        const lookup = &self.lookup.?;
        const head = lookup.heads.getPtr(name) orelse return .missing;
        while (head.* != NO_ENTRY and lookup.taken.isSet(head.*)) head.* = lookup.next[head.*];
        if (head.* == NO_ENTRY) return .omitted;
        const entry = head.*;
        lookup.taken.set(entry);
        head.* = lookup.next[entry];
        return .{ .entry = catalog.records.get(entry) };
    }

    fn buildLookup(self: *Pair) diagnostics.Error!void {
        self.catalog = try index_format.loadIndex(self.allocator, self.io, self.fasta_path, self.selection, .catalog, self.failure);
        const catalog = &self.catalog.?;
        const entries = catalog.records.len();
        if (entries >= NO_ENTRY) return self.failure.fail(.{ .details = .{ .@"resource.capacity_limit_exceeded" = .{ .capacity_role = .name_map, .limit = NO_ENTRY - 1, .required = .{ .exact = entries } } }, .stage = .allocate, .at = .{ .source = self.source() } });
        const oom = diagnostics.Diagnostic{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = .name_map } }, .stage = .allocate, .at = .{ .source = self.source() }, .cause = .{ .zig = error.OutOfMemory } };
        var lookup: Lookup = .{ .heads = .empty, .next = self.allocator.alloc(u32, entries) catch return self.failure.fail(oom), .taken = undefined };
        errdefer self.allocator.free(lookup.next);
        lookup.taken = std.DynamicBitSetUnmanaged.initEmpty(self.allocator, entries) catch return self.failure.fail(oom);
        errdefer lookup.taken.deinit(self.allocator);
        lookup.heads.ensureTotalCapacity(self.allocator, @intCast(entries)) catch return self.failure.fail(oom);
        var entry = entries;
        while (entry > 0) {
            entry -= 1;
            const slot = lookup.heads.getOrPutAssumeCapacity(catalog.recordName(entry).?);
            lookup.next[entry] = if (slot.found_existing) slot.value_ptr.* else NO_ENTRY;
            slot.value_ptr.* = @intCast(entry);
        }
        lookup.taken.setRangeValue(.{ .start = 0, .end = self.taken_in_order }, true);
        self.lookup = lookup;
    }

    // The first difference between the indexed side table of `record` and the scanned one.
    fn sideTableDifference(self: *Pair, record: index_format.IndexRecord, observed: []const u8) diagnostics.Error!?diagnostics.SideTableMismatch {
        const lines = (observed.len - @sizeOf(u64)) / @sizeOf(fasta.LineEntry);
        const mapped: ?[]const index_format.SideTableLine = if (self.catalog) |*catalog| catalog.sideTableLines(record) else null;
        const indexed: u64 = if (mapped) |table| table.len else try self.stream.sideLineCount(record);
        for (0..@min(lines, indexed)) |k| {
            const line = std.mem.bytesToValue(fasta.LineEntry, observed[@sizeOf(u64) + k * @sizeOf(fasta.LineEntry) ..][0..@sizeOf(fasta.LineEntry)]);
            const entry = if (mapped) |table| table[k] else try self.stream.sideLine(record, k);
            if (line.base_start != entry.base_start or line.byte_offset != entry.byte_offset or line.line_bytes != entry.line_bytes or line.line_bases != entry.line_bases) {
                return .{ .indexed_lines = indexed, .observed_lines = lines, .line = k + 1 };
            }
        }
        if (lines == indexed) return null;
        return .{ .indexed_lines = indexed, .observed_lines = lines, .line = @min(lines, indexed) + 1 };
    }
};

const Scan = struct {
    allocator: std.mem.Allocator,
    rules: Rules,
    report: *Report,
    rewrite: ?*Rewrite,
    failure: *diagnostics.Failure,
    names: union(enum) { source: SourceNames, copied: PackedNames, stopped },
    budget: ?NameBudget = null,
    name_buf: std.ArrayList(u8) = .empty,
    first_duplicate: ?[]u8 = null,
    // The fingerprint of the current record's identifier while its duplicate check waits.
    pending_name: ?u32 = null,
    file: ?NameSource,
    classes: [256]u8,
    fast_after_choice: FastSet,
    fast: FastSet = .none,
    chosen: ?u1 = null,
    sample_counts: [256]u64 = @splat(0),
    sample_total: u64 = 0,
    // Counts every byte of a line holding a symbol outside the alphabet; only outside bytes are read.
    byte_counts: [256]u64 = @splat(0),
    candidates: [2]Candidate = .{ .{}, .{} },
    last_records: [REPORT_CODES.len]u64 = @splat(0),
    record: Record = .{},
    line: Line = .{},
    offset: u64 = 0,
    line_number: u64 = 1,
    origin: u64 = 0,
    first_ending: ?LineEnding = null,
    last_byte: ?u8 = null,
    pending_cr: bool = false,
    bom_matched: u2 = 0,
    bom_open: bool = true,
    pair: ?*Pair = null,
    window: []const u8 = &.{},
    window_base: u64 = 0,

    fn init(allocator: std.mem.Allocator, rules: Rules, report: *Report, rewrite: ?*Rewrite, failure: *diagnostics.Failure, file: ?NameSource) Scan {
        const alphabets: [2][256]bool = if (rules.custom_alphabet) |custom|
            .{ buildAlphabet(custom), buildAlphabet(custom) }
        else switch (rules.alphabet) {
            .auto => .{ buildAlphabet(fasta.NUCLEOTIDE_SYMBOLS), buildAlphabet(fasta.PROTEIN_SYMBOLS) },
            .dna => .{ buildAlphabet(fasta.DNA_SYMBOLS), buildAlphabet(fasta.DNA_SYMBOLS) },
            .rna => .{ buildAlphabet(fasta.RNA_SYMBOLS), buildAlphabet(fasta.RNA_SYMBOLS) },
            .protein => .{ buildAlphabet(fasta.PROTEIN_SYMBOLS), buildAlphabet(fasta.PROTEIN_SYMBOLS) },
            .any => .{ @splat(true), @splat(true) },
        };
        var classes: [256]u8 = undefined;
        for (&classes, 0..) |*class, index| {
            const byte: u8 = @intCast(index);
            class.* = if (byte == ' ' or byte == '\t') WHITESPACE else if (byte < ' ') CONTROL else @as(u8, @intFromBool(!alphabets[0][index])) | @as(u8, @intFromBool(!alphabets[1][index])) << 1;
        }
        var self = Scan{
            .allocator = allocator,
            .rules = rules,
            .report = report,
            .rewrite = rewrite,
            .failure = failure,
            .file = file,
            .classes = classes,
            .fast_after_choice = if (rules.custom_alphabet != null) .none else switch (rules.alphabet) {
                .auto => .plain_nucleotides,
                .dna => .plain_dna,
                .rna => .plain_rna,
                .protein => .protein,
                .any => .any,
            },
            .names = .stopped,
            .budget = if (rules.max_name_bytes) |limit| .{ .child = allocator, .limit = limit } else null,
        };
        if (rules.custom_alphabet != null or rules.alphabet != .auto) self.choose(0);
        return self;
    }

    // The name set's allocator may point into `self`, so this runs once `self` has its final address.
    fn startNames(self: *Scan) void {
        const names_allocator = if (self.budget) |*budget| budget.allocator() else self.allocator;
        self.names = if (self.file) |source| .{ .source = SourceNames.init(names_allocator, source.io, source.file, source.size) } else .{ .copied = PackedNames.init(names_allocator) };
        const diag = ScanDiag{ .failure = self.failure, .budget = if (self.budget) |*budget| budget else null };
        switch (self.names) {
            .stopped => unreachable,
            inline else => |*set| set.diag = diag,
        }
    }

    fn deinit(self: *Scan) void {
        self.stopNames();
        self.name_buf.deinit(self.allocator);
        if (self.first_duplicate) |name| self.allocator.free(name);
    }

    fn feed(self: *Scan, data: []const u8) diagnostics.Error!void {
        const base = self.offset;
        self.offset = std.math.add(u64, base, data.len) catch return self.failure.fail(.{ .details = .{ .@"resource.count_overflow" = .{ .counter = .source_bytes, .max_value = std.math.maxInt(u64), .current = base, .increment = data.len } }, .stage = .scan, .at = .{ .source = .fasta } });
        self.window = data;
        self.window_base = base;
        var rest = data;
        var at = base;
        while (self.bom_open and rest.len != 0) {
            if (rest[0] == UTF8_BOM[self.bom_matched]) {
                rest = rest[1..];
                at += 1;
                if (self.bom_matched == 2) {
                    self.bom_open = false;
                    self.count(.{ .@"format.leading_bom" = {} }, fileReference(0, 1, 1));
                } else self.bom_matched += 1;
            } else {
                self.bom_open = false;
                if (self.bom_matched != 0) try self.process(UTF8_BOM[0..self.bom_matched], 0);
            }
        }
        if (rest.len != 0) try self.process(rest, at);
        try self.advancePair(self.offset);
        self.last_byte = data[data.len - 1];
        self.report.examined.bytes = self.offset;
    }

    fn end(self: *Scan) diagnostics.Error!void {
        if (self.bom_open and self.bom_matched != 0) {
            self.bom_open = false;
            try self.process(UTF8_BOM[0..self.bom_matched], 0);
        }
        if (self.pending_cr) {
            self.pending_cr = false;
            try self.content("\r", self.offset - 1);
        }
        if (self.line.kind != .start) try self.endLine(null, self.offset);
        try self.checkPendingName();
        try self.finishPair();
        self.finishRecord();
        if (self.report.examined.records == 0) self.count(.{ .@"fasta.no_records" = {} }, .{ .source = .fasta });
        if (self.last_byte) |byte| {
            if (byte != '\n') self.count(.{ .@"format.missing_terminal_newline" = {} }, self.where(self.offset));
        }
        try self.resolveFirstDefinition();
    }

    fn finish(self: *Scan, completion: reporting.Validation.Completion) void {
        const report = self.report;
        report.completion = if (completion == .complete and report.name_budget_stop != null) .incomplete else completion;
        if (self.chosen == null) self.choose(@intFromBool(fasta.detectType(&self.sample_counts, self.sample_total) == .protein));
        const chosen = self.chosen.?;
        report.alphabet = alphabetPolicy(self.rules, chosen, self.sample_total);
        var occurrences: u64 = 0;
        for (self.byte_counts, &report.invalid_symbols, self.classes) |seen, *invalid, class| {
            invalid.* = if (class & outsideBit(chosen) != 0) seen else 0;
            occurrences += invalid.*;
        }
        const invalid = &report.findings[comptime reportSlot(.@"sequence.invalid_symbol")];
        invalid.occurrences = occurrences;
        invalid.records = self.candidates[chosen].records;
        invalid.first = self.candidates[chosen].first;
        if (invalid.first) |*first| first.details.@"sequence.invalid_symbol".policy = report.alphabet;
        report.records_with_errors = self.candidates[chosen].records_with_errors;
        if (self.pair) |pair| report.index = .{ .format = pair.scanner.diag.format, .entries = pair.entriesRead(), .matched = pair.matched, .omitted_duplicates = pair.omitted };
    }

    fn choose(self: *Scan, candidate: u1) void {
        self.chosen = candidate;
        self.fast = if (self.fast_after_choice == .plain_nucleotides and candidate == 1) .protein else self.fast_after_choice;
    }

    fn process(self: *Scan, data: []const u8, base: u64) diagnostics.Error!void {
        switch (self.names) {
            .source => |*set| set.setChunk(data, base),
            else => {},
        }
        var i: usize = 0;
        if (self.pending_cr) {
            self.pending_cr = false;
            if (data[0] == '\n') {
                try self.endLine(.crlf, base - 1);
                i = 1;
            } else try self.content("\r", base - 1);
        }
        while (i < data.len) {
            if (self.line.kind == .start) {
                const at = base + i;
                if (data[i] == '>') {
                    try self.startHeader(at);
                    i += 1;
                    continue;
                }
                if (self.record.ordinal == 0) {
                    self.line = .{ .kind = .preamble, .start = at };
                } else {
                    const lines = self.denseLines(data, i);
                    if (lines != 0) {
                        try self.takeDenseLines(data[i..], at, lines);
                        i += lines * @as(usize, @intCast(self.record.width.? + self.endingBytes()));
                        continue;
                    }
                    self.line = .{ .kind = .sequence, .start = at };
                }
            }
            const lf = std.mem.indexOfScalarPos(u8, data, i, '\n') orelse data.len;
            if (lf == data.len) {
                var stop = lf;
                if (data[stop - 1] == '\r') {
                    stop -= 1;
                    self.pending_cr = true;
                }
                try self.content(data[i..stop], base + i);
                return;
            }
            const crlf = lf > i and data[lf - 1] == '\r';
            const stop = lf - @intFromBool(crlf);
            try self.content(data[i..stop], base + i);
            try self.endLine(if (crlf) .crlf else .lf, base + stop);
            i = lf + 1;
        }
    }

    fn endingBytes(self: *const Scan) u64 {
        return switch (self.record.last_line.?.ending.?) {
            .lf => 1,
            .crlf => 2,
        };
    }

    fn denseLines(self: *const Scan, data: []const u8, i: usize) usize {
        if (self.fast == .none) return 0;
        const last = self.record.last_line orelse return 0;
        if (!last.dense or last.ending == null) return 0;
        const width = self.record.width orelse last.symbols;
        if (last.symbols != width or width == 0) return 0;
        const line_bytes = width + self.endingBytes();
        if (line_bytes > std.math.maxInt(u32)) return 0;
        return switch (self.fast) {
            .none => 0,
            .plain_nucleotides => fasta.countDenseLines(fasta.PlainNucleotides, data, i, @intCast(width), @intCast(line_bytes)),
            .plain_dna => fasta.countDenseLines(fasta.PlainDna, data, i, @intCast(width), @intCast(line_bytes)),
            .plain_rna => fasta.countDenseLines(fasta.PlainRna, data, i, @intCast(width), @intCast(line_bytes)),
            .protein => fasta.countDenseLines(fasta.ProteinSymbols, data, i, @intCast(width), @intCast(line_bytes)),
            .any => fasta.countDenseLines(fasta.Bases, data, i, @intCast(width), @intCast(line_bytes)),
        };
    }

    fn takeDenseLines(self: *Scan, data: []const u8, at: u64, lines: usize) diagnostics.Error!void {
        const width = self.record.width orelse self.record.last_line.?.symbols;
        self.record.width = width;
        const symbols_per_line: usize = @intCast(width);
        const line_bytes: usize = symbols_per_line + @as(usize, @intCast(self.endingBytes()));
        if (self.rewrite) |sink| sink.lines(data, symbols_per_line, line_bytes, lines) catch |err| return self.writeFailed(err);
        const symbols = width * lines;
        self.record.symbols += symbols;
        self.report.examined.symbols += symbols;
        self.report.examined.lines += lines;
        const last_start = at + (lines - 1) * line_bytes;
        self.line_number += lines - 1;
        self.record.last_line = .{
            .symbols = width,
            .start = .{ .offset = last_start, .line = self.line_number, .column = 1 },
            .dense = true,
            .ending = self.record.last_line.?.ending,
        };
        self.line_number += 1;
        self.origin = at + lines * line_bytes;
        self.report.examined.bytes = self.origin;
    }

    fn startHeader(self: *Scan, at: u64) diagnostics.Error!void {
        try self.checkPendingName();
        try self.advancePair(at + 1);
        self.finishRecord();
        self.report.examined.records += 1;
        self.record = .{ .ordinal = self.report.examined.records };
        self.record.header = self.position(at);
        self.line = .{ .kind = .header, .start = at, .in_identifier = true };
        self.name_buf.clearRetainingCapacity();
        if (self.rewrite) |sink| sink.startHeader() catch |err| return self.writeFailed(err);
    }

    fn content(self: *Scan, bytes: []const u8, at: u64) diagnostics.Error!void {
        switch (self.line.kind) {
            .start => unreachable,
            .preamble => {
                if (self.line.content_at == null) {
                    for (bytes, 0..) |byte, index| {
                        if (byte != ' ' and byte != '\t' and byte != '\r') {
                            self.line.content_at = self.line.len + index;
                            break;
                        }
                    }
                }
            },
            .header => try self.headerContent(bytes, at),
            .sequence => try self.sequenceContent(bytes, at),
        }
        self.line.len += bytes.len;
    }

    fn headerContent(self: *Scan, bytes: []const u8, at: u64) diagnostics.Error!void {
        if (self.line.in_identifier) {
            const stop = fasta.findNameDelimiter(bytes, 0);
            const kept = @min(stop, MAX_INDEX_NAME_BYTES - self.name_buf.items.len);
            self.name_buf.appendSlice(self.allocator, bytes[0..kept]) catch return self.outOfMemory(.validation_identifiers, at);
            self.line.identifier_bytes += stop;
            if (stop < bytes.len) self.line.in_identifier = false;
        }
        for (bytes, 0..) |byte, index| {
            if (byte < ' ' and byte != '\t') {
                if (self.line.first_control == null) self.line.first_control = .{ .byte = byte, .index = self.line.len + index };
                self.line.controls += 1;
            }
        }
        if (self.rewrite) |sink| sink.header(bytes) catch |err| return self.writeFailed(err);
    }

    fn sequenceContent(self: *Scan, bytes: []const u8, at: u64) diagnostics.Error!void {
        if (bytes.len == 0) return;
        const all_fast = switch (self.fast) {
            .none => false,
            .plain_nucleotides => !fasta.PlainNucleotides.anyOutside(DENSE_VECTOR_LEN, bytes),
            .plain_dna => !fasta.PlainDna.anyOutside(DENSE_VECTOR_LEN, bytes),
            .plain_rna => !fasta.PlainRna.anyOutside(DENSE_VECTOR_LEN, bytes),
            .protein => !fasta.ProteinSymbols.anyOutside(DENSE_VECTOR_LEN, bytes),
            .any => !fasta.Bases.anyOutside(DENSE_VECTOR_LEN, bytes),
        };
        if (all_fast) {
            self.line.symbols += bytes.len;
            self.report.examined.symbols += bytes.len;
            if (self.rewrite) |sink| sink.symbols(bytes) catch |err| return self.writeFailed(err);
            return;
        }
        self.line.dense = false;
        const outside: u8 = if (self.chosen) |chosen| outsideBit(chosen) else outsideBit(0) | outsideBit(1);
        var classes: u8 = 0;
        for (bytes) |byte| classes |= self.classes[byte];
        if (classes & (WHITESPACE | CONTROL) != 0) {
            try self.irregularContent(bytes, at, outside);
        } else {
            if (self.chosen == null) self.addSample(bytes);
            const found = classes & outside;
            if (found != 0) {
                for (bytes) |byte| self.byte_counts[byte] += 1;
                self.record.classes |= found;
                inline for (0..2) |candidate| {
                    if (found & outsideBit(candidate) != 0 and self.candidates[candidate].first == null) {
                        const index = for (bytes, 0..) |byte, index| {
                            if (self.classes[byte] & outsideBit(candidate) != 0) break index;
                        } else unreachable;
                        self.keepFirstOutside(candidate, bytes[index], at + index);
                    }
                }
            }
            self.line.symbols += bytes.len;
            self.report.examined.symbols += bytes.len;
            if (self.rewrite) |sink| sink.symbols(bytes) catch |err| return self.writeFailed(err);
        }
        if (self.chosen == null and self.sample_total == TYPE_SAMPLE_SYMBOLS) self.choose(@intFromBool(fasta.detectType(&self.sample_counts, self.sample_total) == .protein));
    }

    fn irregularContent(self: *Scan, bytes: []const u8, at: u64, outside: u8) diagnostics.Error!void {
        const mask = WHITESPACE | CONTROL | outside;
        var run_start: usize = 0;
        var symbols: u64 = 0;
        for (bytes, 0..) |byte, index| {
            const class = self.classes[byte] & mask;
            if (class & (WHITESPACE | CONTROL) != 0) {
                if (self.rewrite) |sink| sink.symbols(bytes[run_start..index]) catch |err| return self.writeFailed(err);
                run_start = index + 1;
                if (class == WHITESPACE) {
                    if (self.line.whitespace == 0) self.line.first_whitespace = self.line.len + index;
                    self.line.whitespace += 1;
                } else {
                    self.count(.{ .@"fasta.sequence_control_byte" = .{ .byte = byte, .grammar = .controls_except_tab } }, self.reference(at + index));
                }
                continue;
            }
            symbols += 1;
            if (self.chosen == null) self.addSample(bytes[index..][0..1]);
            if (class != 0) {
                self.byte_counts[byte] += 1;
                self.record.classes |= class;
                inline for (0..2) |candidate| {
                    if (class & outsideBit(candidate) != 0 and self.candidates[candidate].first == null) self.keepFirstOutside(candidate, byte, at + index);
                }
            }
        }
        if (self.rewrite) |sink| sink.symbols(bytes[run_start..]) catch |err| return self.writeFailed(err);
        self.line.symbols += symbols;
        self.report.examined.symbols += symbols;
    }

    fn addSample(self: *Scan, symbols: []const u8) void {
        const take = @min(symbols.len, TYPE_SAMPLE_SYMBOLS - self.sample_total);
        for (symbols[0..take]) |byte| self.sample_counts[byte] += 1;
        self.sample_total += take;
    }

    fn keepFirstOutside(self: *Scan, comptime candidate: u1, byte: u8, offset: u64) void {
        self.candidates[candidate].first = .{
            .details = .{ .@"sequence.invalid_symbol" = .{ .byte = byte, .policy = self.report.alphabet } },
            .at = self.reference(offset),
            .identifier = self.identifierPreview(),
        };
    }

    fn endLine(self: *Scan, ending: ?LineEnding, at: u64) diagnostics.Error!void {
        if (ending) |kind| {
            if (self.first_ending) |first| {
                if (kind != first and self.report.finding(.@"format.mixed_line_endings").occurrences == 0) self.count(.{ .@"format.mixed_line_endings" = {} }, self.where(at));
            } else self.first_ending = kind;
        }
        switch (self.line.kind) {
            .start => unreachable,
            .preamble => if (self.line.content_at) |index| {
                self.count(.{ .@"fasta.content_before_header" = .{ .content_bytes = .{ .exact = self.line.len } } }, self.where(self.line.start + index));
            },
            .header => try self.endHeader(),
            .sequence => self.endSequenceLine(ending),
        }
        self.report.examined.lines += 1;
        if (ending) |kind| {
            self.line_number += 1;
            self.origin = at + @as(u64, if (kind == .crlf) 2 else 1);
            self.report.examined.bytes = self.origin;
        } else self.report.examined.bytes = at;
        self.line = .{};
    }

    fn endHeader(self: *Scan) diagnostics.Error!void {
        const identifier = self.name_buf.items;
        const start = self.line.start;
        self.record.named = true;
        self.record.identifier_bytes = self.line.identifier_bytes;
        if (identifier.len == 0) {
            self.count(.{ .@"fasta.empty_identifier" = {} }, self.reference(start + 1));
        } else if (self.line.identifier_bytes > MAX_INDEX_NAME_BYTES) {
            // No index can hold it, so it takes no part in the duplicate check.
            self.count(.{ .@"fasta.identifier_too_long" = .{ .observed_bytes = self.line.identifier_bytes, .limit_bytes = MAX_INDEX_NAME_BYTES } }, self.reference(start + 1));
        } else switch (self.names) {
            .stopped => {},
            .source => |*set| self.pending_name = set.prefetch(identifier),
            .copied => self.pending_name = 0,
        }
        if (self.line.first_control) |control| {
            self.countMany(.{ .@"fasta.header_control_byte" = .{ .byte = control.byte, .grammar = .controls_except_tab } }, self.reference(start + 1 + control.index), self.line.controls);
        }
        if (self.line.len > self.rules.max_header_bytes) self.count(.{ .@"fasta.long_header" = .{ .observed_bytes = self.line.len, .limit_bytes = self.rules.max_header_bytes } }, self.located(self.record.header));
        if (convention(self.rules.schema)) |selected| {
            if (!matchesSchema(self.rules.schema, identifier)) self.count(.{ .@"header.convention_mismatch" = .{ .convention = selected } }, self.located(self.record.header));
        }
        if (self.rewrite) |sink| sink.endHeader() catch |err| return self.writeFailed(err);
    }

    fn stopNames(self: *Scan) void {
        switch (self.names) {
            .stopped => {},
            inline else => |*set| set.deinit(),
        }
        self.names = .stopped;
    }

    // The duplicate check `endHeader` deferred, so the table slots it prefetched have arrived. It
    // runs before anything reads `record.duplicate`: the next header, or the end of input.
    fn checkPendingName(self: *Scan) diagnostics.Error!void {
        const fingerprint = self.pending_name orelse return;
        self.pending_name = null;
        const identifier = self.name_buf.items;
        const header = self.record.header;
        const at: diagnostics.Reference = .{ .source = .fasta, .location = .{ .offset = header.offset + 1, .line = header.line, .column = header.column + 1 }, .record = self.record.ordinal };
        if (!try self.seenBefore(identifier, fingerprint, at)) return;
        self.record.duplicate = true;
        if (self.first_duplicate == null) self.first_duplicate = self.allocator.dupe(u8, identifier) catch return self.outOfMemoryAt(.validation_identifiers, self.located(header));
        self.count(.{ .@"fasta.duplicate_identifier" = .{} }, at);
    }

    fn seenBefore(self: *Scan, identifier: []const u8, fingerprint: u32, at: diagnostics.Reference) diagnostics.Error!bool {
        const offset = at.location.offset.?;
        const found = switch (self.names) {
            .stopped => return false,
            .source => |*set| set.observeFingerprint(identifier, offset, fingerprint),
            .copied => |*set| set.observeAt(identifier, offset),
        } catch |err| {
            if (self.budget != null and self.budget.?.exhausted) {
                self.report.name_budget_stop = .{ .budget_bytes = self.budget.?.limit, .stopped_at_record = self.record.ordinal };
                self.stopNames();
                return false;
            }
            if (self.failure.primary()) |primary| return diagnostics.metadata(primary.code()).returned_error;
            return self.failure.fail(.{ .details = .{ .@"internal.consistency_failed" = .{ .component = .validator, .check = .retained_name } }, .stage = .scan, .at = at, .cause = .{ .zig = err } });
        };
        return found;
    }

    fn endSequenceLine(self: *Scan, ending: ?LineEnding) void {
        const line = self.line;
        if (line.len == 0) return self.count(.{ .@"format.blank_sequence_line" = {} }, self.reference(line.start));
        if (line.whitespace != 0) self.count(.{ .@"format.sequence_whitespace" = .{ .whitespace_bytes = line.whitespace } }, self.reference(line.start + line.first_whitespace));
        if (line.symbols == 0) return;
        self.record.symbols += line.symbols;
        if (self.record.last_line) |previous| {
            const width = self.record.width orelse previous.symbols;
            self.record.width = width;
            if (previous.symbols != width) self.reportWrapping(width, previous);
        }
        self.record.last_line = .{ .symbols = line.symbols, .start = self.position(line.start), .dense = line.dense and line.symbols == line.len, .ending = ending };
    }

    fn reportWrapping(self: *Scan, width: u64, line: SymbolLine) void {
        if (self.record.wrapping_reported) return;
        self.record.wrapping_reported = true;
        self.count(.{ .@"format.inconsistent_wrapping" = .{ .expected_width = width, .actual_width = line.symbols } }, self.located(line.start));
    }

    fn finishRecord(self: *Scan) void {
        const record = &self.record;
        if (record.ordinal == 0) return;
        if (record.last_line) |last| {
            if (record.width) |width| {
                if (last.symbols > width) self.reportWrapping(width, last);
            }
        }
        if (record.symbols == 0) self.count(.{ .@"fasta.empty_record" = {} }, self.located(record.header));
        inline for (0..2) |candidate| {
            const outside = record.classes & outsideBit(candidate) != 0;
            if (outside) self.candidates[candidate].records += 1;
            if (outside or record.has_error) self.candidates[candidate].records_with_errors += 1;
        }
        if (record.has_warning) self.report.records_with_warnings += 1;
        self.report.examined.records_completed += 1;
        record.ordinal = 0;
    }

    fn count(self: *Scan, details: diagnostics.Details, at: diagnostics.Reference) void {
        self.countMany(details, at, 1);
    }

    fn countMany(self: *Scan, details: diagnostics.Details, at: diagnostics.Reference, occurrences: u64) void {
        const code = std.meta.activeTag(details);
        const slot = REPORT_SLOTS[@intFromEnum(code)].?;
        const entry = &self.report.findings[slot];
        entry.occurrences += occurrences;
        if (at.record) |record| {
            if (self.last_records[slot] != record) {
                self.last_records[slot] = record;
                entry.records += 1;
            }
            switch (diagnostics.metadata(code).severity) {
                .@"error" => self.record.has_error = true,
                .warning => self.record.has_warning = true,
            }
        }
        if (entry.first == null) entry.first = .{ .details = details, .at = at, .identifier = if (at.record != null) self.identifierPreview() else null };
    }

    // Re-reading the source up to the first duplicate costs nothing on files without one.
    fn resolveFirstDefinition(self: *Scan) diagnostics.Error!void {
        const name = self.first_duplicate orelse return;
        const source = self.file orelse return;
        const first = &(self.report.findings[comptime reportSlot(.@"fasta.duplicate_identifier")].first orelse return);
        const stop = first.at.location.offset.?;
        const buffer = self.allocator.alloc(u8, @max(FIRST_DEFINITION_BUFFER_BYTES, name.len + UTF8_BOM.len + 2)) catch return self.outOfMemory(.validation_identifiers, stop);
        defer self.allocator.free(buffer);
        var offset: u64 = 0;
        var line: u64 = 1;
        var record: u64 = 0;
        var line_start: u64 = 0;
        var at_line_start = true;
        while (offset < stop) {
            const got = source.file.readPositionalAll(source.io, buffer, offset) catch |err| return self.failure.fail(.{ .details = diagnostics.readFailure(buffer.len, null), .stage = .read, .at = .{ .source = .fasta, .location = .{ .offset = offset } }, .cause = .{ .zig = err } });
            const data = buffer[0..got];
            const full = got == buffer.len;
            var i: usize = 0;
            while (i < data.len) {
                if (at_line_start) {
                    if (offset + i >= stop) return;
                    const skip: usize = if (offset + i == 0 and std.mem.startsWith(u8, data, UTF8_BOM)) UTF8_BOM.len else 0;
                    const marker = i + skip;
                    if (marker < data.len and data[marker] == '>') {
                        const id_at = marker + 1;
                        if (data.len - id_at <= name.len and full and i != 0) break;
                        record += 1;
                        if (data.len - id_at > name.len and std.mem.eql(u8, data[id_at..][0..name.len], name) and fasta.isNameDelimiter(data[id_at + name.len])) {
                            const id_offset = offset + id_at;
                            first.related = .{ .source = .fasta, .location = .{ .offset = id_offset, .line = line, .column = id_offset - line_start + 1 }, .record = record };
                            return;
                        }
                    }
                    at_line_start = false;
                }
                const lf = std.mem.indexOfScalarPos(u8, data, i, '\n') orelse {
                    i = data.len;
                    break;
                };
                line += 1;
                i = lf + 1;
                line_start = offset + i;
                at_line_start = true;
            }
            if (!full) return;
            offset += i;
        }
    }

    // Feeds `index`'s scanner the current window up to `to`. Fed through a header's `>`, it emits
    // the record before that header while `record` still holds it.
    fn advancePair(self: *Scan, to: u64) diagnostics.Error!void {
        const pair = self.pair orelse return;
        if (to <= pair.fed) return;
        const bytes = self.window[@intCast(pair.fed - self.window_base)..@intCast(to - self.window_base)];
        const from = pair.fed;
        pair.fed = to;
        pair.scanner.processChunk(bytes, from, MAX_INDEX_NAME_BYTES, self, compareRecord) catch |err| return self.pairFailed(err);
    }

    fn finishPair(self: *Scan) diagnostics.Error!void {
        const pair = self.pair orelse return;
        if (pair.scanner.active) pair.scanner.emitRecord(self, compareRecord) catch |err| return self.pairFailed(err);
        if (pair.catalog) |*catalog| {
            const lookup = &pair.lookup.?;
            for (0..catalog.records.len()) |entry| {
                if (lookup.taken.isSet(entry)) continue;
                const offset = if (pair.source() == .zfi) @sizeOf(index_format.ZfiHeader) + entry * index_format.ZFI_INDEX_RECORD_BYTES else catalog.records.get(entry).name_offset;
                self.countMissingEntry(offset, entry + 1, catalog.recordName(entry).?);
            }
            return;
        }
        if (pair.head) |head| self.countMissingEntry(head.offset, head.ordinal, head.name);
        pair.head = null;
        while (try pair.stream.next()) |entry| self.countMissingEntry(entry.offset, entry.ordinal, entry.name);
    }

    fn countMissingEntry(self: *Scan, offset: u64, ordinal: u64, name: []const u8) void {
        const finding = &self.report.findings[comptime reportSlot(.@"pair.record_missing_from_fasta")];
        finding.occurrences += 1;
        finding.records += 1;
        if (finding.first == null) finding.first = .{ .details = .{ .@"pair.record_missing_from_fasta" = {} }, .at = .{ .source = self.pair.?.source(), .location = .{ .offset = offset }, .record = ordinal }, .identifier = .init(name) };
    }

    fn compareRecord(self: *Scan, record: fasta.ScannedRecord) anyerror!void {
        const pair = self.pair.?;
        const current = self.record.ordinal != 0 and self.record.header.offset == record.header_offset;
        const at: diagnostics.Reference = if (current) self.located(self.record.header) else .{ .source = .fasta, .location = .{ .offset = record.header_offset } };
        const entry = switch (try pair.associate(record.name, current and self.record.duplicate)) {
            .entry => |entry| entry,
            .omitted => {
                pair.omitted += 1;
                return;
            },
            .missing => return self.count(.{ .@"pair.record_missing_from_index" = {} }, at),
        };
        pair.matched += 1;
        if (entry.seq_len != record.seq_len) self.count(.{ .@"pair.sequence_length_mismatch" = .{ .indexed_length = entry.seq_len, .observed_length = record.seq_len } }, at);
        if (entry.seq_offset != record.seq_offset) self.count(.{ .@"pair.sequence_offset_mismatch" = .{ .indexed_offset = entry.seq_offset, .observed_offset = record.seq_offset } }, at);
        const indexed_layout: diagnostics.LineLayout = if (entry.isUniformWidth()) .uniform else .side_table;
        const observed_layout: diagnostics.LineLayout = if (record.uses_uniform_formula) .uniform else .side_table;
        const geometry_differs = indexed_layout == .uniform and (entry.line_bases != record.line_bases or entry.line_bytes != record.line_bytes);
        if (indexed_layout != observed_layout or geometry_differs) {
            self.count(.{ .@"pair.line_geometry_mismatch" = .{ .indexed_line_bases = entry.line_bases, .indexed_line_bytes = entry.line_bytes, .indexed_layout = indexed_layout, .observed_line_bases = record.line_bases, .observed_line_bytes = record.line_bytes, .observed_layout = observed_layout } }, at);
        } else if (indexed_layout == .side_table) {
            if (try pair.sideTableDifference(entry, record.side_table)) |mismatch| self.count(.{ .@"pair.side_table_mismatch" = mismatch }, at);
        }
    }

    fn pairFailed(self: *Scan, err: anyerror) diagnostics.Error {
        @branchHint(.cold);
        if (self.failure.primary()) |primary| return diagnostics.metadata(primary.code()).returned_error;
        return self.failure.fail(.{ .details = .{ .@"internal.consistency_failed" = .{ .component = .validator, .check = .side_table_state } }, .stage = .scan, .at = .{ .source = .fasta, .location = .{ .offset = self.offset } }, .cause = .{ .zig = err } });
    }

    fn position(self: *const Scan, offset: u64) Position {
        return .{ .offset = offset, .line = self.line_number, .column = offset - self.origin + 1 };
    }

    // A `position` in the current record.
    fn located(self: *const Scan, start: Position) diagnostics.Reference {
        return .{ .source = .fasta, .location = .{ .offset = start.offset, .line = start.line, .column = start.column }, .record = self.record.ordinal };
    }

    fn identifierPreview(self: *const Scan) ?reporting.Identifier {
        if (!self.record.named) return null;
        var preview = reporting.Identifier.init(self.name_buf.items);
        preview.original_bytes = self.record.identifier_bytes;
        return preview;
    }

    fn reference(self: *const Scan, offset: u64) diagnostics.Reference {
        return .{ .source = .fasta, .location = .{ .offset = offset, .line = self.line_number, .column = offset - self.origin + 1 }, .record = if (self.record.ordinal == 0) null else self.record.ordinal };
    }

    fn where(self: *const Scan, offset: u64) diagnostics.Reference {
        return fileReference(offset, self.line_number, offset - self.origin + 1);
    }

    fn outOfMemory(self: *Scan, role: diagnostics.AllocationRole, offset: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.outOfMemoryAt(role, self.reference(offset));
    }

    fn outOfMemoryAt(self: *Scan, role: diagnostics.AllocationRole, at: diagnostics.Reference) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = .{ .@"resource.out_of_memory" = .{ .allocation_role = role } }, .stage = .allocate, .at = at, .cause = .{ .zig = error.OutOfMemory } });
    }

    fn writeFailed(self: *Scan, err: anyerror) diagnostics.Error {
        @branchHint(.cold);
        const file = if (self.rewrite) |sink| sink.file else null;
        const cause = if (file) |writer| writer.err orelse err else err;
        return self.failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .at = .{ .source = .output }, .cause = .{ .zig = cause } });
    }
};

const DENSE_VECTOR_LEN = std.simd.suggestVectorLength(u8) orelse 1;
const FIRST_DEFINITION_BUFFER_BYTES = 64 * 1024;

fn fileReference(offset: u64, line: u64, column: u64) diagnostics.Reference {
    return .{ .source = .fasta, .location = .{ .offset = offset, .line = line, .column = column } };
}

fn matchesSchema(schema: Schema, identifier: []const u8) bool {
    return switch (schema) {
        .none => true,
        .uniprot => matchesUniprot(identifier),
        .refseq => matchesRefseq(identifier),
    };
}

fn matchesUniprot(header: []const u8) bool {
    if (!(std.mem.startsWith(u8, header, "sp|") or
        std.mem.startsWith(u8, header, "tr|") or
        std.mem.startsWith(u8, header, "db|")))
    {
        return false;
    }
    const rest = header[3..];
    const first_pipe = std.mem.indexOfScalar(u8, rest, '|') orelse return false;
    return first_pipe > 0 and first_pipe + 1 < rest.len;
}
fn matchesRefseq(header: []const u8) bool {
    const prefixes = [_][]const u8{
        "NC_", "NM_", "NR_", "NP_", "NW_", "NG_", "NT_", "NZ_", "XM_", "XR_", "XP_",
    };
    for (prefixes) |prefix| {
        if (std.mem.startsWith(u8, header, prefix)) return true;
    }
    return false;
}
fn cleanupFixed(file: *std.Io.File.Atomic, io: std.Io, failure: *diagnostics.Failure) void {
    if (file.file_open) {
        file.file.close(io);
        file.file_open = false;
    }
    if (file.file_exists) {
        var removed = true;
        const basename = std.fmt.hex(file.file_basename_hex);
        file.dir.deleteFile(io, &basename) catch |err| {
            removed = err == error.FileNotFound;
            if (!removed) {
                failure.sources.set(.staged_output, .{ .label = .{ .retained = failure.retain(&basename, 32) }, .kind = .path });
                switch (failure.fail(.{ .details = .{ .@"io.cleanup_failed" = .{ .resource_role = .staged_output } }, .stage = .cleanup, .at = .{ .source = .staged_output }, .cause = .{ .zig = err } })) {
                    else => {},
                }
            }
        };
        if (failure.outputs.file.state == .staged) failure.outputs.file.state = if (removed) .aborted else .unknown;
        file.file_exists = false;
    }
    file.deinit(io);
}
fn unknownOutputIdentity(failure: *diagnostics.Failure, source: diagnostics.SourceId, cause: anyerror) diagnostics.Error {
    return failure.fail(.{
        .details = .{ .@"repair.output_identity_unknown" = .{ .input = .fasta, .output = .output } },
        .stage = .check_identity,
        .at = .{ .source = source },
        .cause = .{ .zig = cause },
    });
}
fn checkRewriteDestination(io: std.Io, input_path: []const u8, output_path: []const u8, failure: *diagnostics.Failure) diagnostics.Error!void {
    if (std.mem.eql(u8, input_path, output_path)) return failure.fail(.{
        .details = .{ .@"repair.output_is_input" = .{ .input = .fasta, .output = .output, .comparison = .lexical_path } },
        .stage = .check_identity,
        .at = .{ .source = .output },
    });
    const output_file = std.Io.Dir.cwd().openFile(io, output_path, .{}) catch |err| {
        if (err == error.FileNotFound) return;
        return unknownOutputIdentity(failure, .output, err);
    };
    output_file.close(io);
    var input_buf: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const input_entry = canonicalEntry(io, input_path, &input_buf) catch |err| return unknownOutputIdentity(failure, .fasta, err);
    var output_buf: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const output_entry = canonicalEntry(io, output_path, &output_buf) catch |err| return unknownOutputIdentity(failure, .output, err);
    if (std.mem.eql(u8, input_entry, output_entry)) return failure.fail(.{
        .details = .{ .@"repair.output_is_input" = .{ .input = .fasta, .output = .output, .comparison = .canonical_path } },
        .stage = .check_identity,
        .at = .{ .source = .output },
    });
}
fn canonicalEntry(io: std.Io, path: []const u8, out: *[std.Io.Dir.max_path_bytes]u8) ![]const u8 {
    var hop_bufs: [2][std.Io.Dir.max_path_bytes]u8 = undefined;
    var target_buf: [std.Io.Dir.max_path_bytes]u8 = undefined;
    var current = path;
    var hops: usize = 0;
    while (true) : (hops += 1) {
        if (hops == MAX_LINK_HOPS) return error.SymLinkLoop;
        const target_len = std.Io.Dir.cwd().readLink(io, current, &target_buf) catch |err| switch (err) {
            error.NotLink => break,
            else => return err,
        };
        const target = target_buf[0..target_len];
        const next = &hop_bufs[hops % 2];
        const parent = std.fs.path.dirname(current);
        current = if (std.fs.path.isAbsolute(target) or parent == null)
            std.fmt.bufPrint(next, "{s}", .{target}) catch return error.NameTooLong
        else
            std.fmt.bufPrint(next, "{s}{c}{s}", .{ parent.?, std.fs.path.sep, target }) catch return error.NameTooLong;
    }
    const name = std.fs.path.basename(current);
    var dir = try std.Io.Dir.cwd().openDir(io, std.fs.path.dirname(current) orelse ".", .{});
    defer dir.close(io);
    const dir_len = try dir.realPath(io, out);
    if (dir_len + 1 + name.len > out.len) return error.NameTooLong;
    out[dir_len] = std.fs.path.sep;
    @memcpy(out[dir_len + 1 ..][0..name.len], name);
    return out[0 .. dir_len + 1 + name.len];
}
fn checkRewriteSource(io: std.Io, input: std.Io.File, before: std.Io.File.Stat, failure: *diagnostics.Failure) diagnostics.Error!void {
    const after = input.stat(io) catch |err| return failure.fail(.{
        .details = .{ .@"io.stat_failed" = {} },
        .stage = .compare_source,
        .at = .{ .source = .fasta },
        .cause = .{ .zig = err },
    });
    if (before.size != after.size or before.mtime.nanoseconds != after.mtime.nanoseconds) return failure.fail(.{
        .details = .{ .@"input.changed" = .{
            .changed_fields = .{ .size = before.size != after.size, .mtime = before.mtime.nanoseconds != after.mtime.nanoseconds },
            .before = .{ .bytes = before.size, .mtime_ns = before.mtime.nanoseconds },
            .after = .{ .bytes = after.size, .mtime_ns = after.mtime.nanoseconds },
        } },
        .stage = .compare_source,
        .at = .{ .source = .fasta },
    });
}

fn writeStaged(io: std.Io, data: []const u8, path: []const u8, input: std.Io.File, before: std.Io.File.Stat, failure: *diagnostics.Failure) diagnostics.Error!void {
    var staged: Staged = undefined;
    try staged.create(std.testing.allocator, io, path, failure);
    defer staged.cleanup(io, failure);
    staged.writer.interface.writeAll(data) catch |err| return failure.fail(.{ .details = .{ .@"io.write_failed" = {} }, .stage = .write, .at = .{ .source = .output }, .cause = .{ .zig = staged.writer.err orelse err } });
    try staged.publish(io, .{ .file = input, .before = before }, failure);
}

test "[failure] - [validate fix]: an unresolved output identity keeps its cause" {
    const Fault = struct {
        fn realPath(_: ?*anyopaque, _: std.Io.Dir, _: []u8) std.Io.Dir.RealPathError!usize {
            return error.AccessDenied;
        }
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const input = try tmp.dir.createFile(std.testing.io, "input.fa", .{});
    input.close(std.testing.io);
    const output = try tmp.dir.createFile(std.testing.io, "output.fa", .{});
    output.close(std.testing.io);
    const input_path = try tmp.dir.realPathFileAlloc(std.testing.io, "input.fa", std.testing.allocator);
    defer std.testing.allocator.free(input_path);
    const output_path = try tmp.dir.realPathFileAlloc(std.testing.io, "output.fa", std.testing.allocator);
    defer std.testing.allocator.free(output_path);
    var vtable = std.testing.io.vtable.*;
    vtable.dirRealPath = Fault.realPath;
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    try std.testing.expectError(error.RepairRefused, checkRewriteDestination(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, input_path, output_path, &failure));
    try std.testing.expectEqual(diagnostics.SourceId.fasta, failure.primary().?.at.source.?);
    try std.testing.expectEqual(diagnostics.Code.@"repair.output_identity_unknown", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
}

test "[failure] - [validate fix]: an observed source change prevents publication" {
    const Fault = struct {
        fn stat(userdata: ?*anyopaque, file: std.Io.File) std.Io.File.StatError!std.Io.File.Stat {
            var value = try std.testing.io.vtable.fileStat(userdata, file);
            value.size += 1;
            return value;
        }
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const input = try tmp.dir.createFile(std.testing.io, "input.fa", .{ .read = true });
    defer input.close(std.testing.io);
    const data = ">seq1\nACGT";
    try std.Io.File.writeStreamingAll(input, std.testing.io, data);
    const before = try input.stat(std.testing.io);
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "output.fa", .data = "previous destination\n" });
    const path = try tmp.dir.realPathFileAlloc(std.testing.io, "output.fa", std.testing.allocator);
    defer std.testing.allocator.free(path);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    var vtable = std.testing.io.vtable.*;
    vtable.fileStat = Fault.stat;
    try std.testing.expectError(error.InputChanged, writeStaged(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, ">seq1\nACGT\n", path, input, before, &failure));
    const detail = failure.primary().?.details.@"input.changed";
    try std.testing.expect(detail.changed_fields.size and !detail.changed_fields.mtime);
    try std.testing.expectEqual(@as(?u64, data.len), detail.before.bytes);
    try std.testing.expectEqual(@as(?u64, data.len + 1), detail.after.?.bytes);
    try std.testing.expectEqual(.aborted, failure.outputs.file.state);
    const after = try tmp.dir.readFileAlloc(std.testing.io, "output.fa", std.testing.allocator, .limited(1024));
    defer std.testing.allocator.free(after);
    try std.testing.expectEqualStrings("previous destination\n", after);
}

test "[failure] - [validate fix]: a failed publication stays unknown after staged cleanup" {
    const Fault = struct {
        var scenario: enum { published, refused, cleanup_failed } = .published;

        fn rename(userdata: ?*anyopaque, old_dir: std.Io.Dir, old_path: []const u8, new_dir: std.Io.Dir, new_path: []const u8) std.Io.Dir.RenameError!void {
            if (scenario == .published) try std.testing.io.vtable.dirRename(userdata, old_dir, old_path, new_dir, new_path);
            return error.AccessDenied;
        }
        fn remove(_: ?*anyopaque, _: std.Io.Dir, _: []const u8) std.Io.Dir.DeleteFileError!void {
            return error.AccessDenied;
        }
    };
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const input = try tmp.dir.createFile(io, "in.fa", .{ .read = true });
    defer input.close(io);
    const data = ">seq\nACGT\n";
    try std.Io.File.writeStreamingAll(input, io, data);
    const before = try input.stat(io);
    try tmp.dir.writeFile(io, .{ .sub_path = "out.fa", .data = "old output\n" });
    const out_path = try tmp.dir.realPathFileAlloc(io, "out.fa", std.testing.allocator);
    defer std.testing.allocator.free(out_path);
    for ([_]@TypeOf(Fault.scenario){ .published, .refused, .cleanup_failed }) |scenario| {
        Fault.scenario = scenario;
        try tmp.dir.writeFile(io, .{ .sub_path = "out.fa", .data = "old output\n" });
        var vtable = io.vtable.*;
        vtable.dirRename = Fault.rename;
        if (scenario == .cleanup_failed) vtable.dirDeleteFile = Fault.remove;
        var failure: diagnostics.Failure = undefined;
        failure.reset(.validate);
        try std.testing.expectError(error.Io, writeStaged(.{ .userdata = io.userdata, .vtable = &vtable }, data, out_path, input, before, &failure));
        const after = try tmp.dir.readFileAlloc(io, "out.fa", std.testing.allocator, .limited(1024));
        defer std.testing.allocator.free(after);
        try std.testing.expectEqualStrings(if (scenario == .published) data else "old output\n", after);
        const primary = failure.primary().?;
        try std.testing.expectEqual(diagnostics.Code.@"io.publish_failed", primary.code());
        try std.testing.expectEqual(error.AccessDenied, primary.cause.?.zig);
        try std.testing.expectEqual(diagnostics.SourceId.output, primary.at.source);
        try std.testing.expectEqual(.unknown, failure.outputs.file.state);
        try std.testing.expectEqual(.succeeded, failure.outputs.file.final_flush);
        try std.testing.expectEqual(@as(usize, if (scenario == .cleanup_failed) 1 else 0), failure.secondaryItems().len);
        if (scenario == .cleanup_failed) {
            try std.testing.expectEqual(diagnostics.Code.@"io.cleanup_failed", failure.secondary[0].code());
            try std.testing.expectEqual(error.AccessDenied, failure.secondary[0].cause.?.zig);
        }
    }
}

test "[failure] - [validate fix]: a failed flush aborts staging and keeps the destination" {
    const Fault = struct {
        fn write(_: ?*anyopaque, _: std.Io.File, _: []const u8, _: []const []const u8, _: usize, _: u64) std.Io.File.WritePositionalError!usize {
            return error.InputOutput;
        }
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const input = try tmp.dir.createFile(std.testing.io, "in.fa", .{ .read = true });
    defer input.close(std.testing.io);
    const data = ">seq\nACGT\n";
    try std.Io.File.writeStreamingAll(input, std.testing.io, data);
    const before = try input.stat(std.testing.io);
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "out.fa", .data = "old output\n" });
    const path = try tmp.dir.realPathFileAlloc(std.testing.io, "out.fa", std.testing.allocator);
    defer std.testing.allocator.free(path);
    var vtable = std.testing.io.vtable.*;
    vtable.fileWritePositional = Fault.write;
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    try std.testing.expectError(error.Io, writeStaged(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, data, path, input, before, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.InputOutput, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(.failed, failure.outputs.file.final_flush);
    try std.testing.expectEqual(.aborted, failure.outputs.file.state);
    try std.testing.expectEqual(@as(u8, 0), failure.secondary_len);
    const after = try tmp.dir.readFileAlloc(std.testing.io, "out.fa", std.testing.allocator, .limited(1024));
    defer std.testing.allocator.free(after);
    try std.testing.expectEqualStrings("old output\n", after);
}

fn checkPackedNameAllocations(allocator: std.mem.Allocator) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    const diag = ScanDiag{ .failure = &failure };
    errdefer |err| {
        if (err == error.OutOfMemory) {
            std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code()) catch @panic("missing allocation diagnostic");
            std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage) catch @panic("wrong allocation stage");
        }
    }

    var names = PackedNames.init(allocator);
    names.diag = diag;
    defer names.deinit();
    try names.ensureCapacity(16);
    var name = [_]u8{'A'} ** MAX_INDEX_NAME_BYTES;
    for (0..3) |i| {
        name[0] = @intCast('A' + i);
        try std.testing.expect(!try names.observeAt(&name, 0));
    }
    for (0..300) |i| {
        var buffer: [32]u8 = undefined;
        const short = try std.fmt.bufPrint(&buffer, "name-{d}", .{i});
        try std.testing.expect(!try names.observeAt(short, 0));
    }
    for (0..3) |i| {
        name[0] = @intCast('A' + i);
        try std.testing.expect(try names.observeAt(&name, 0));
    }
}

fn checkSourceNameAllocations(allocator: std.mem.Allocator, data: []const u8) !void {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    const diag = ScanDiag{ .failure = &failure };
    errdefer |err| {
        if (err == error.OutOfMemory) {
            std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code()) catch @panic("missing allocation diagnostic");
            std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage) catch @panic("wrong allocation stage");
        }
    }

    var names = SourceNames.init(allocator, std.testing.io, undefined, data.len);
    names.diag = diag;
    defer names.deinit();
    names.setChunk(data, 0);
    for ([_]bool{ false, true }) |duplicate| {
        var lines = std.mem.tokenizeScalar(u8, data, '\n');
        while (lines.next()) |name| {
            const offset = @intFromPtr(name.ptr) - @intFromPtr(data.ptr);
            try std.testing.expectEqual(duplicate, try names.observe(name, offset));
        }
    }
    try std.testing.expectEqual(@as(u32, 5000), names.tables.narrow.count());
}

test "[edge] - [name catalog]: full source tables retain duplicate lookups" {
    inline for (.{ u32, u64 }) |Offset| {
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        var source = SourceNames.init(failing.allocator(), std.testing.io, undefined, if (Offset == u32) 0 else std.math.maxInt(u32));
        defer source.deinit();
        source.cache_full = true;
        source.setChunk("alpha\nbeta\n", 0);
        const set = if (Offset == u32) &source.tables.narrow else &source.tables.wide;
        try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
        set.len = std.math.maxInt(u32);
        failing.fail_index = failing.alloc_index;
        try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", 0, 0));
        try std.testing.expectError(error.OutOfMemory, set.observeWithFingerprint(&source, "beta", 6, 0));
        try std.testing.expect(!failing.has_induced_failure);
        try std.testing.expectEqual(std.math.maxInt(u32), set.len);
    }
}

test "[property] - [name catalog]: complete names resolve fingerprint collisions" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const long = "a-much-longer-different-name";
    const data = long ++ "\nalpha\nalphabet\nbeta\nalpha\n";
    try file.writePositionalAll(io, data, 0);
    var source = SourceNames.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    source.cache_full = true;
    const set = &source.tables.narrow;
    const start = long.len + 1;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", start, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", start + 6, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "beta", start + 15, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", start + 20, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, long, 0, 0));
}

test "[unit] - [name catalog]: exact comparisons use the current input chunk" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const data = "alpha\nalphabet\nbeta\nalpha\n";
    var source = SourceNames.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    source.cache_full = true;
    source.setChunk(data, 0);
    const set = &source.tables.narrow;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", 6, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", 20, 0));
}

test "[edge] - [name catalog]: retain offsets above four GiB" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const offset: u64 = @as(u64, std.math.maxInt(u32)) + 64;
    const data = "alpha\nalphabet\nbeta\nalpha\n";
    var source = SourceNames.init(std.testing.allocator, io, file, offset + data.len);
    defer source.deinit();
    source.cache_full = true;
    source.setChunk(data, offset);
    const set = &source.tables.wide;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", offset, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", offset + 6, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", offset + 20, 0));
    try std.testing.expectEqual(offset, set.keyAt(0).offset);
    try std.testing.expectEqual(@as(u64, 0), (try file.stat(io)).size);
}

test "[property] - [name catalog]: cached names resolve fingerprint collisions" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    const data = "alpha\nalphabet\nbeta\nalpha\n";
    try file.writePositionalAll(io, data, 0);
    var source = SourceNames.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    const set = &source.tables.narrow;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alphabet", 6, 0));
    try std.testing.expect(!try set.observeWithFingerprint(&source, "beta", 15, 0));
    try std.testing.expect(try set.observeWithFingerprint(&source, "alpha", 20, 0));
}

test "[failure] - [name catalog]: read errors leave the insertion unfinished" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{});
    defer file.close(io);
    const data = "alpha\nbeta\n";
    try file.writePositionalAll(io, data, 0);
    var source = SourceNames.init(std.testing.allocator, io, file, data.len);
    defer source.deinit();
    source.cache_full = true;
    const set = &source.tables.narrow;
    try std.testing.expect(!try set.observeWithFingerprint(&source, "alpha", 0, 0));
    try std.testing.expectError(error.SourceReadFailed, set.observeWithFingerprint(&source, "beta", 6, 0));
    try std.testing.expectEqual(@as(u32, 1), set.count());
}

test "[failure] - [name catalog]: block growth cleans every allocation failure" {
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkPackedNameAllocations, .{});
}

test "[failure] - [name catalog]: table and entry growth clean every allocation failure" {
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(std.testing.allocator);
    for (0..5000) |i| {
        var buffer: [32]u8 = undefined;
        try data.appendSlice(std.testing.allocator, try std.fmt.bufPrint(&buffer, "lookup-{d:0>8}\n", .{i}));
    }
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkSourceNameAllocations, .{data.items});
}

test "[edge] - [name catalog]: maximum names cross copied and retained cache blocks" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    var name: [MAX_INDEX_NAME_BYTES]u8 = undefined;
    const stride = MAX_INDEX_NAME_BYTES + 1;
    for (0..6) |i| {
        @memset(&name, @intCast('A' + i));
        try file.writePositionalAll(io, &name, i * stride);
        try file.writePositionalAll(io, "\n", i * stride + name.len);
    }
    var names = SourceNames.init(std.testing.allocator, io, file, 6 * stride);
    defer names.deinit();
    for ([_]bool{ false, true, true }) |duplicate| {
        for (0..6) |i| {
            @memset(&name, @intCast('A' + i));
            try std.testing.expectEqual(duplicate, try names.observe(&name, i * stride));
        }
    }
    try std.testing.expect(names.cache_full);
    try std.testing.expect(names.hot_blocks.items.len >= 2);
}

test "[failure] - [name catalog]: short reads distinguish collisions from truncation" {
    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    const io = std.testing.io;
    const file = try temporary.dir.createFile(io, "names.fa", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, "alpha\n", 0);
    var names = SourceNames.init(std.testing.allocator, io, file, 6);
    defer names.deinit();
    names.cache_full = true;
    const flag = SourceNames.Tables(u32).CACHE_FLAG;
    try std.testing.expect(!names.sameName("alphabet", 0, flag));
    try std.testing.expect(!names.read_failed);
    try file.setLength(io, 3);
    try std.testing.expect(!names.sameName("alpha", 0, flag));
    try std.testing.expect(names.read_failed);
}

test "[failure] - [name catalog]: source offsets cannot overlap cache flags" {
    var narrow = SourceNames.init(std.testing.allocator, std.testing.io, undefined, 0);
    defer narrow.deinit();
    try std.testing.expectError(error.SourceChanged, narrow.observe("name", SourceNames.Tables(u32).CACHE_FLAG));
    try std.testing.expectEqual(@as(u32, 0), narrow.tables.narrow.count());
    var wide = SourceNames.init(std.testing.allocator, std.testing.io, undefined, std.math.maxInt(u32));
    defer wide.deinit();
    try std.testing.expectError(error.SourceChanged, wide.observe("name", SourceNames.Tables(u64).CACHE_FLAG));
    try std.testing.expectEqual(@as(u32, 0), wide.tables.wide.count());
}

test "[failure] - [name catalog]: oversized offsets roll back the insertion" {
    if (comptime @bitSizeOf(usize) <= 32) return error.SkipZigTest;
    var names = PackedNames.init(std.testing.allocator);
    defer names.deinit();
    names.position = @as(usize, std.math.maxInt(u32)) + 1;
    try std.testing.expectError(error.OutOfMemory, names.observeAt("name", 0));
    try std.testing.expectEqual(@as(u32, 0), names.map.count());
}

test "[failure] - [name catalog]: cache representation limits are capacity failures" {
    var failure: diagnostics.Failure = undefined;
    inline for (.{ u32, u64 }) |Offset| {
        failure.reset(.validate);
        var names = SourceNames.init(std.testing.allocator, std.testing.io, undefined, if (Offset == u32) 0 else std.math.maxInt(u32));
        names.diag = .{ .failure = &failure };
        defer names.deinit();
        const limit = SourceNames.Tables(Offset).CACHE_FLAG;
        try std.testing.expectError(error.CapacityExceeded, names.observe("name", limit));
        const details = failure.primary().?.details.@"resource.capacity_limit_exceeded";
        try std.testing.expectEqual(diagnostics.CapacityRole.name_cache_offset, details.capacity_role);
        try std.testing.expectEqual(@as(u128, limit - 1), details.limit);
        try std.testing.expectEqual(@as(u128, limit), details.required.exact);
    }
    if (comptime @bitSizeOf(usize) > 32) {
        failure.reset(.validate);
        var names = PackedNames.init(std.testing.allocator);
        names.diag = .{ .failure = &failure };
        defer names.deinit();
        names.position = @as(usize, std.math.maxInt(u32)) + 1;
        try std.testing.expectError(error.CapacityExceeded, names.observeAt("name", 0));
        try std.testing.expectEqual(diagnostics.CapacityRole.name_offset, failure.primary().?.details.@"resource.capacity_limit_exceeded".capacity_role);
        try std.testing.expectEqual(@as(u32, 0), names.map.count());
    }
}

test "[failure] - [validate fix]: an output buffer allocation failure leaves the destination as it was" {
    var tmp = std.testing.tmpDir(.{ .iterate = true });
    defer tmp.cleanup();
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "input.fa", .data = ">seq1\nACGT\n" });
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "output.fa", .data = "previous destination\n" });
    const input_path = try tmp.dir.realPathFileAlloc(std.testing.io, "input.fa", std.testing.allocator);
    defer std.testing.allocator.free(input_path);
    const output_path = try tmp.dir.realPathFileAlloc(std.testing.io, "output.fa", std.testing.allocator);
    defer std.testing.allocator.free(output_path);
    var fail_index: usize = 0;
    const failure = while (fail_index < 16) : (fail_index += 1) {
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = fail_index });
        var failure: diagnostics.Failure = undefined;
        _ = runValidate(failing.allocator(), std.testing.io, input_path, .{ .fix = true, .output_path = output_path }, &failure) catch {
            const details = failure.primary().?.details;
            if (details == .@"resource.out_of_memory" and details.@"resource.out_of_memory".allocation_role == .output_buffer) break failure;
            continue;
        };
        return error.OutputBufferNeverFailed;
    } else return error.OutputBufferNeverFailed;
    try std.testing.expectEqual(diagnostics.Stage.allocate, failure.primary().?.stage);
    try std.testing.expectEqual(@as(?u64, OUTPUT_BUFFER_BYTES), failure.primary().?.details.@"resource.out_of_memory".requested_bytes);
    var buffer: [64]u8 = undefined;
    try std.testing.expectEqualStrings("previous destination\n", try tmp.dir.readFile(std.testing.io, "output.fa", &buffer));
    var entries: usize = 0;
    var iterator = tmp.dir.iterate();
    while (try iterator.next(std.testing.io)) |_| entries += 1;
    try std.testing.expectEqual(@as(usize, 2), entries);
}
