//! FASTA structure, alphabet, and header validation with optional fix streaming.
//!
//! Events are capped (`MAX_VALIDATE_EVENTS`). JSON output escapes arbitrary header bytes.
//! Shared sequence-type sampling uses `stats.detectType`.

const std = @import("std");

const platform = @import("platform.zig");
const stats = @import("stats.zig");
const diagnostics = @import("diagnostics");
const reporting = @import("reporting.zig");

/// Hard cap on retained validate events. Further issues set `Summary.truncated`
/// instead of growing an unbounded list.
pub const MAX_VALIDATE_EVENTS: usize = 10_000;

/// Selects the stable text, JSON-lines, or aggregate JSON report.
pub const OutputMode = enum {
    text,
    json_lines,
    json_summary,
};

/// Optional header convention checked in addition to ordinary FASTA validation.
pub const Schema = enum {
    none,
    uniprot,
    refseq,
};

/// Controls validation, reporting, and format-only rewriting.
///
/// `output_path` and `custom_alphabet` are borrowed for the call. Direct CLI-equivalent
/// callers must pair `fix_format_only` and `output_path` with `fix`.
pub const Options = struct {
    strict: bool = false,
    output_mode: OutputMode = .text,
    fix: bool = false,
    fix_format_only: bool = false,
    output_path: ?[]const u8 = null,
    schema: Schema = .none,
    custom_alphabet: ?[]const u8 = null,
    max_header_len: usize = DEFAULT_MAX_HEADER_LEN,
};

/// Severity used for report events and exit-status accounting.
pub const Level = enum {
    error_level,
    warning,

    fn text(self: Level) []const u8 {
        return switch (self) {
            .error_level => "error",
            .warning => "warning",
        };
    }

    fn label(self: Level) []const u8 {
        return switch (self) {
            .error_level => "ERROR",
            .warning => "WARNING",
        };
    }
};

/// Stable validation issue families used by text and JSON reports.
pub const Kind = enum {
    no_sequences,
    duplicate_name,
    invalid_character,
    null_byte,
    utf8_bom,
    inconsistent_line_widths,
    trailing_whitespace,
    empty_sequence,
    missing_terminal_newline,
    mixed_line_endings,
    long_header,
    schema_violation,

    fn text(self: Kind) []const u8 {
        return switch (self) {
            .no_sequences => "no_sequences",
            .duplicate_name => "duplicate_name",
            .invalid_character => "invalid_character",
            .null_byte => "null_byte",
            .utf8_bom => "utf8_bom",
            .inconsistent_line_widths => "inconsistent_line_widths",
            .trailing_whitespace => "trailing_whitespace",
            .empty_sequence => "empty_sequence",
            .missing_terminal_newline => "missing_terminal_newline",
            .mixed_line_endings => "mixed_line_endings",
            .long_header => "long_header",
            .schema_violation => "schema_violation",
        };
    }
};

/// One retained issue. `name` borrows the validated input; other optional fields apply only
/// to the issue kinds that expose them in the report schema.
pub const Event = struct {
    level: Level,
    kind: Kind,
    line: usize,
    name: []const u8 = "",
    first_line: usize = 0,
    byte: u8 = 0,
    expected_width: u32 = 0,
    actual_width: u32 = 0,
    limit: usize = 0,
    offset: ?u64 = null,
    column: ?u64 = null,
    record: ?u64 = null,
    control_place: enum { header, sequence, preamble } = .sequence,
};

/// Owns retained events and rewrite widths while borrowing event names from the input.
/// Counts remain complete when `truncated` is true, although `events` does not.
pub const Summary = struct {
    events: std.ArrayList(Event),
    record_widths: std.ArrayList(u32),
    sequence_type: stats.SequenceType = .nucleotide,
    /// Bases fed into `detectType` (capped at `stats.VALIDATE_TYPE_SAMPLE_BASES`).
    type_bases_sampled: u64 = 0,
    sequence_count: usize = 0,
    header_count: usize = 0,
    error_count: usize = 0,
    warning_count: usize = 0,
    format_warning_count: usize = 0,
    kind_counts: [KIND_COUNT]usize = .{0} ** KIND_COUNT,
    /// True when more issues were found after `MAX_VALIDATE_EVENTS` were retained.
    truncated: bool = false,
    /// Fix blocker found after event retention reached its cap.
    truncated_fix_rejection: ?Kind = null,
    first_blocker: ?Event = null,
    first_invalid: ?Event = null,
    preamble: ?PreamblePosition = null,

    /// Releases owned storage with the allocator passed to `validateData`.
    pub fn deinit(self: *Summary, allocator: std.mem.Allocator) void {
        self.events.deinit(allocator);
        self.record_widths.deinit(allocator);
    }
};

/// Scalar validation outcome after complete report delivery. No input names borrow
/// from the released mapping; findings and their legacy counts stay in stdout.
pub const Result = struct { status: u8, truncated: bool, fix_written: bool };

/// Validates, optionally rewrites, and delivers the legacy report. Resets once
/// and returns after cleanup. Source path labels borrow caller strings through
/// handling/reset; the input mapping remains alive while report names borrow it.
pub fn runValidate(allocator: std.mem.Allocator, io: std.Io, fasta_path: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!Result {
    failure.reset(.validate);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = fasta_path }, .kind = .path });
    failure.sources.set(.stdout, .{ .label = .{ .borrowed = "stdout" }, .kind = .stream });
    if (options.output_path) |path| failure.sources.set(.output, .{ .label = .{ .borrowed = path }, .kind = .path });
    var diag = Diagnostics{ .failure = failure };
    const file = std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| return diag.fail(switch (err) {
        error.FileNotFound => .{ .@"io.file_not_found" = {} },
        error.AccessDenied => .{ .@"io.access_denied" = {} },
        error.NameTooLong => .{ .@"io.invalid_path" = .{ .reason = .too_long } },
        else => .{ .@"io.open_failed" = {} },
    }, .open, .fasta, err);
    defer file.close(io);
    const stat = file.stat(io) catch |err| return diag.fail(.{ .@"io.stat_failed" = {} }, .stat, .fasta, err);
    var fasta_map: ?std.Io.File.MemoryMap = null;
    defer if (fasta_map) |*map| map.destroy(io);
    const data: []const u8 = if (stat.size == 0) &.{} else blk: {
        const len = std.math.cast(usize, stat.size) orelse return diag.fail(.{ .@"resource.capacity_limit_exceeded" = .{
            .capacity_role = .map_bytes,
            .limit = std.math.maxInt(usize),
            .required = .{ .exact = stat.size },
        } }, .allocate, .fasta, null);
        fasta_map = platform.mapFileReadOnly(io, file, len) catch |err| return diag.fail(.{ .@"io.map_failed" = .{ .map_bytes = len } }, .map, .fasta, err);
        break :blk fasta_map.?.memory;
    };
    if (fasta_map) |*map| platform.adviseSequential(map.memory);
    var summary = try validateDataInternal(allocator, data, options, &diag);
    defer summary.deinit(allocator);
    failure.work_complete = true;
    diag.at = .{ .source = .fasta };
    if (options.fix) {
        const output_path = options.output_path orelse return diag.fail(.{ .@"usage.option_requires" = .{ .option = .fix, .required = .output } }, .arguments, .fasta, null);
        try checkRewriteDestination(io, file, fasta_path, output_path, failure);
        try fixRejectionChecked(&summary, options, failure);
        try writeFixed(io, data, summary.record_widths.items, output_path, file, stat, failure);
    }
    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    errdefer stdout_fw.interface.end = 0;
    failure.outputs.stdout.state = .streaming;
    failure.outputs.stdout.may_be_partial = true;
    writeSummary(&stdout_fw.interface, &summary, options) catch |err| return diag.fail(.{ .@"io.write_failed" = {} }, .report, .stdout, stdout_fw.err orelse err);
    stdout_fw.flush() catch |err| {
        failure.outputs.stdout.final_flush = .failed;
        return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .stdout, stdout_fw.err orelse err);
    };
    failure.outputs.stdout.final_flush = .succeeded;
    const result = Result{ .status = if (summary.truncated and !options.fix) 1 else exitCodeForOptions(&summary, options), .truncated = summary.truncated, .fix_written = options.fix };
    failure.succeed();
    return result;
}

/// Scans all input bytes and returns an owning summary whose event names borrow `data`.
/// Retained events are capped, but counts, fix blockers, and rewrite widths remain complete.
/// Compatibility projection: native callers use validateDataChecked. A sequence
/// line above u32 returns SequenceLineTooLong; allocation/count failures keep tags.
pub fn validateData(allocator: std.mem.Allocator, data: []const u8, options: Options) error{ OutOfMemory, SequenceLineTooLong, CountOverflow, Internal }!Summary {
    var failure: diagnostics.Failure = undefined;
    return validateDataChecked(allocator, data, options, &failure) catch |err| {
        if (failure.primary().?.code() == .@"index.sequence_line_too_long") return error.SequenceLineTooLong;
        return switch (err) {
            error.OutOfMemory => error.OutOfMemory,
            error.CountOverflow => error.CountOverflow,
            else => error.Internal,
        };
    };
}

/// Scans borrowed bytes with native failure facts. Resets once; returned Summary
/// owns lists but event/blocker names borrow data until Summary.deinit.
pub fn validateDataChecked(allocator: std.mem.Allocator, data: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!Summary {
    failure.reset(.validate);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "validated input" }, .kind = .memory });
    var diag = Diagnostics{ .failure = failure };
    const result = try validateDataInternal(allocator, data, options, &diag);
    failure.succeed();
    return result;
}

/// Returns the validate status for complete severity counts and strictness.
pub fn exitCode(error_count: usize, warning_count: usize, strict: bool) u8 {
    if (error_count > 0) return 1;
    if (warning_count > 0) return if (strict) 1 else 2;
    return 0;
}

/// Returns validate status after excluding warnings repaired by a requested format rewrite.
pub fn exitCodeForOptions(summary: *const Summary, options: Options) u8 {
    if (!options.fix) return exitCode(summary.error_count, summary.warning_count, options.strict);

    return exitCode(
        summary.error_count,
        summary.warning_count - summary.format_warning_count,
        options.strict,
    );
}

/// Returns the first legacy Kind that blocks a rewrite, including past-cap events.
/// Use fixRejectionChecked to include preamble safety and retain original facts.
pub fn fixRejection(summary: *const Summary, options: Options) ?Kind {
    for (summary.events.items) |event| {
        switch (event.kind) {
            .invalid_character => if (!options.fix_format_only) return .invalid_character,
            .no_sequences, .duplicate_name, .null_byte => return event.kind,
            else => {},
        }
    }
    if (summary.truncated_fix_rejection) |kind| {
        if (kind != .invalid_character or !options.fix_format_only) return kind;
    }
    return null;
}

/// Records the original blocker for a requested fix without resetting Failure.
/// Summary event names must still borrow live input; retained previews own bytes.
/// Use the options that produced Summary, with the requested fix policy.
pub fn fixRejectionChecked(summary: *const Summary, options: Options, failure: *diagnostics.Failure) diagnostics.Error!void {
    if (summary.preamble) |position| return failure.fail(.{
        .details = .{ .@"repair.blocked" = .{ .blocker = .{ .content_before_header = .{ .content_bytes = .{ .at_least = 1 } } }, .repair_policy = if (options.fix_format_only) .legacy_format_only else .legacy_safe } },
        .stage = .scan,
        .at = .{ .source = .fasta, .location = .{ .offset = position.offset, .line = position.line, .column = position.column } },
    });
    const event = summary.first_blocker orelse if (!options.fix_format_only) summary.first_invalid orelse return else return;
    const blocker: diagnostics.Blocker = switch (event.kind) {
        .no_sequences => .{ .no_records = {} },
        .duplicate_name => .{ .duplicate_identifier = .{ .identifier = failure.retain(event.name, 32) } },
        .null_byte => if (event.control_place == .header) .{ .header_control_byte = .{ .byte = 0, .grammar = .legacy_nul_only } } else if (event.control_place == .sequence) .{ .sequence_control_byte = .{ .byte = 0, .grammar = .legacy_nul_only } } else .{ .content_before_header = .{ .content_bytes = .{ .at_least = 1 } } },
        .invalid_character => blk: {
            var policy = diagnostics.AlphabetPolicy{ .alphabet = if (options.custom_alphabet != null) .custom else if (summary.sequence_type == .nucleotide) .nucleotide else .protein, .selection = if (options.custom_alphabet != null) .custom else .automatic, .sample = if (options.custom_alphabet == null) .{ .symbols = @intCast(summary.type_bases_sampled), .rule = .nucleotide_above_90_percent } else null };
            if (options.custom_alphabet) |custom| {
                var bits: [32]u8 = .{0} ** 32;
                for (custom) |byte| bits[byte / 8] |= @as(u8, 1) << @intCast(byte % 8);
                policy.custom_symbols = bits;
            }
            break :blk .{ .invalid_symbol = .{ .byte = event.byte, .policy = policy } };
        },
        else => return failure.fail(.{ .details = .{ .@"internal.consistency_failed" = .{ .component = .validator, .check = .codec_state } }, .stage = .scan, .at = eventReference(event) }),
    };
    return failure.fail(.{ .details = .{ .@"repair.blocked" = .{ .blocker = blocker, .repair_policy = if (options.fix_format_only) .legacy_format_only else .legacy_safe } }, .stage = .scan, .at = eventReference(event), .related = if (event.kind == .duplicate_name) .{ .source = .fasta, .location = .{ .line = event.first_line } } else null });
}

/// Rewrite format-level issues to normalized LF FASTA bytes. Caller must run validateData first.
/// Allocates the full result for tests and library callers; CLI `--fix` streams via `writeFixed`.
/// Allocator-backed writer failures are reported as `error.OutOfMemory`.
/// Non-format content before the first header returns `error.RepairRefused`.
/// Returns normalized bytes allocated with `allocator`; release with `allocator.free`.
/// Input and widths borrow only during this call.
/// Compatibility helper; native callers use fixDataChecked.
pub fn fixData(allocator: std.mem.Allocator, data: []const u8, record_widths: []const u32) diagnostics.Error![]u8 {
    var failure: diagnostics.Failure = undefined;
    return fixDataChecked(allocator, data, record_widths, &failure);
}

/// Returns normalized bytes allocated with `allocator`; release with `allocator.free`.
/// Input and widths borrow only during this call. Captures native failure facts;
/// resets once, rejects unsafe preamble content and never publishes a destination.
pub fn fixDataChecked(allocator: std.mem.Allocator, data: []const u8, record_widths: []const u32, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
    failure.reset(.validate);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "rewrite input" }, .kind = .memory });
    var diag = Diagnostics{ .failure = failure };
    try checkPreamble(data, .{}, failure);
    var aw: std.Io.Writer.Allocating = .init(allocator);
    errdefer aw.deinit();
    writeFixedContent(&aw.writer, data, record_widths) catch |err| return diag.oom(.rewrite_bytes, err);
    const result = aw.toOwnedSlice() catch |err| return diag.oom(.rewrite_bytes, err);
    failure.succeed();
    return result;
}

/// Returns one JSON event object allocated with `allocator`, without a trailing newline.
/// Release with `allocator.free`; event names borrow only during this call.
/// Allocator-backed writer failures are reported as `error.OutOfMemory`.
/// Compatibility projection; native callers use renderJsonEventChecked.
pub fn renderJsonEvent(allocator: std.mem.Allocator, event: Event) error{OutOfMemory}![]u8 {
    var failure: diagnostics.Failure = undefined;
    return renderJsonEventChecked(allocator, event, &failure) catch return error.OutOfMemory;
}

/// Returns a legacy v1 event allocated with `allocator`; release with `allocator.free`.
/// Captures native allocation facts and resets once; names borrow only during this
/// call. Invalid UTF-8 escaping is not byte-reversible.
pub fn renderJsonEventChecked(allocator: std.mem.Allocator, event: Event, failure: *diagnostics.Failure) diagnostics.Error![]u8 {
    failure.reset(.validate);
    failure.sources.set(.fasta, .{ .label = .{ .borrowed = "validation event" }, .kind = .memory });
    var diag = Diagnostics{ .failure = failure, .at = eventReference(event) };
    var aw: std.Io.Writer.Allocating = .init(allocator);
    errdefer aw.deinit();
    writeJsonEventObject(&aw.writer, event) catch |err| return diag.oom(.serialized_event, err);
    const result = aw.toOwnedSlice() catch |err| return diag.oom(.serialized_event, err);
    failure.succeed();
    return result;
}

const Diagnostics = struct {
    failure: *diagnostics.Failure,
    at: diagnostics.Reference = .{ .source = .fasta },

    fn fail(self: *@This(), details: diagnostics.Details, stage: diagnostics.Stage, source: diagnostics.SourceId, cause: ?anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{ .details = details, .stage = stage, .at = if (self.at.source == source) self.at else .{ .source = source }, .cause = if (cause) |err| .{ .zig = err } else null });
    }

    fn oom(self: *@This(), role: diagnostics.AllocationRole, cause: anyerror) diagnostics.Error {
        return self.fail(.{ .@"resource.out_of_memory" = .{ .allocation_role = role } }, .allocate, .fasta, cause);
    }

    fn add(self: *@This(), counter: diagnostics.Counter, current: usize, increment: usize) diagnostics.Error!usize {
        return std.math.add(usize, current, increment) catch return self.overflow(counter, current, increment);
    }

    fn overflow(self: *@This(), counter: diagnostics.Counter, current: usize, increment: usize) diagnostics.Error {
        @branchHint(.cold);
        return self.fail(.{ .@"resource.count_overflow" = .{ .counter = counter, .max_value = std.math.maxInt(usize), .current = current, .increment = increment } }, .scan, .fasta, null);
    }
};

const UTF8_BOM = "\xEF\xBB\xBF";

const DEFAULT_MAX_HEADER_LEN = 1024;

const DEFAULT_FIX_WIDTH = 60;

const SIMD_CHUNK_SIZE = 32;

const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);

comptime {
    // Keep the validate help line in src/main.zig synchronized with this value.
    if (MAX_VALIDATE_EVENTS != 10_000) {
        @compileError("update validate help text for MAX_VALIDATE_EVENTS");
    }
}

const KIND_COUNT = @typeInfo(Kind).@"enum".fields.len;

const WidthCount = struct {
    width: u32,
    count: usize,
};

const RecordState = struct {
    active: bool = false,
    name: []const u8 = "",
    header_line: usize = 0,
    has_bases: bool = false,
    expected_width: ?u32 = null,
    pending_width: ?u32 = null,
    pending_line: usize = 0,
    width_warning_emitted: bool = false,
    widths: std.ArrayList(WidthCount) = .empty,

    fn reset(self: *RecordState) void {
        self.active = false;
        self.name = "";
        self.header_line = 0;
        self.has_bases = false;
        self.expected_width = null;
        self.pending_width = null;
        self.pending_line = 0;
        self.width_warning_emitted = false;
        self.widths.clearRetainingCapacity();
    }
};

fn validateDataInternal(allocator: std.mem.Allocator, data: []const u8, options: Options, diag: *Diagnostics) diagnostics.Error!Summary {
    var summary = Summary{
        .events = .empty,
        .record_widths = .empty,
    };
    errdefer summary.deinit(allocator);

    var seen_names: std.StringHashMapUnmanaged(usize) = .empty;
    defer seen_names.deinit(allocator);

    var type_counts: [256]u64 = .{0} ** 256;
    var type_total: u64 = 0;
    var byte_counts: [256]u64 = .{0} ** 256;
    var first_byte_at: [256]diagnostics.Reference = .{diagnostics.Reference{}} ** 256;

    var state = RecordState{};
    defer state.widths.deinit(allocator);

    if (std.mem.startsWith(u8, data, UTF8_BOM)) {
        try appendEvent(allocator, &summary, diag, .{
            .level = .warning,
            .kind = .utf8_bom,
            .line = 1,
        });
    }

    if (data.len == 0 or data[data.len - 1] != '\n') {
        try appendEvent(allocator, &summary, diag, .{
            .level = .warning,
            .kind = .missing_terminal_newline,
            .line = try countLines(data, diag),
        });
    }

    var saw_crlf = false;
    var saw_lf = false;
    var first_crlf_line: usize = 0;
    var first_lf_line: usize = 0;

    var pos: usize = 0;
    var line_number: usize = 1;
    if (std.mem.startsWith(u8, data, UTF8_BOM)) pos = UTF8_BOM.len;

    while (pos <= data.len) {
        const line_start = pos;
        const first_column: usize = if (line_number == 1) line_start + 1 else 1;
        const line_end = findNextNewline(data, pos);

        if (line_start == data.len and line_end == data.len) break;

        var content_end = line_end;
        const has_lf = line_end < data.len;
        const has_crlf = has_lf and content_end > line_start and data[content_end - 1] == '\r';
        if (has_crlf) content_end -= 1;

        if (has_crlf) {
            if (!saw_crlf) first_crlf_line = line_number;
            saw_crlf = true;
        } else if (has_lf) {
            if (!saw_lf) first_lf_line = line_number;
            saw_lf = true;
        }

        const line = data[line_start..content_end];
        diag.at = .{ .source = .fasta, .location = .{ .offset = line_start, .line = line_number, .column = first_column }, .record = if (state.active) summary.header_count else null };
        if (line.len > 0 and line[0] == '>') {
            try finalizeRecord(allocator, &summary, &state, diag);
            diag.at = .{ .source = .fasta, .location = .{ .offset = line_start, .line = line_number, .column = first_column } };
            try startRecord(allocator, &summary, &seen_names, &state, line, line_number, options, diag);
        } else if (state.active) {
            try scanSequenceLine(allocator, &summary, &state, line, line_number, &type_counts, &type_total, &byte_counts, &first_byte_at, diag);
        } else {
            try scanNullsOutsideSequence(allocator, &summary, line, line_number, diag);
            if (summary.preamble == null) {
                for (line, 0..) |byte, column| {
                    if (byte != ' ' and byte != '\t' and byte != '\r') {
                        summary.preamble = .{ .offset = line_start + column, .line = line_number, .column = first_column + column };
                        break;
                    }
                }
            }
        }

        if (!has_lf) break;
        pos = line_end + 1;
        line_number = try diag.add(.physical_lines, line_number, 1);
    }

    try finalizeRecord(allocator, &summary, &state, diag);

    if (summary.header_count == 0) {
        try appendEvent(allocator, &summary, diag, .{
            .level = .error_level,
            .kind = .no_sequences,
            .line = 1,
        });
    }

    if (saw_crlf and saw_lf) {
        try appendEvent(allocator, &summary, diag, .{
            .level = .warning,
            .kind = .mixed_line_endings,
            .line = @max(first_crlf_line, first_lf_line),
        });
    }

    summary.sequence_type = stats.detectType(&type_counts, type_total);
    summary.type_bases_sampled = type_total;
    try appendInvalidCharacterEvents(allocator, &summary, options, summary.sequence_type, byte_counts, first_byte_at, diag);
    return summary;
}

fn isFixedByFormatRewrite(kind: Kind) bool {
    return switch (kind) {
        .utf8_bom,
        .inconsistent_line_widths,
        .trailing_whitespace,
        .missing_terminal_newline,
        .mixed_line_endings,
        => true,
        else => false,
    };
}

fn eventReference(event: Event) diagnostics.Reference {
    return .{ .source = .fasta, .location = .{ .offset = event.offset, .line = event.line, .column = event.column }, .record = event.record };
}

fn appendEvent(allocator: std.mem.Allocator, summary: *Summary, diag: *Diagnostics, original: Event) diagnostics.Error!void {
    const scan_at = diag.at;
    defer diag.at = scan_at;
    var event = original;
    if (diag.at.location.line == event.line) {
        if (event.offset == null) event.offset = diag.at.location.offset;
        if (event.column == null) event.column = diag.at.location.column;
        if (event.record == null) event.record = diag.at.record;
    }
    diag.at = eventReference(event);
    switch (event.kind) {
        .duplicate_name, .null_byte, .no_sequences => if (summary.first_blocker == null) {
            summary.first_blocker = event;
        },
        .invalid_character => if (summary.first_invalid == null) {
            summary.first_invalid = event;
        },
        else => {},
    }
    switch (event.level) {
        .error_level => summary.error_count = try diag.add(.events, summary.error_count, 1),
        .warning => summary.warning_count = try diag.add(.events, summary.warning_count, 1),
    }
    const kind_index = @intFromEnum(event.kind);
    summary.kind_counts[kind_index] = try diag.add(.events, summary.kind_counts[kind_index], 1);
    if (isFixedByFormatRewrite(event.kind)) summary.format_warning_count = try diag.add(.events, summary.format_warning_count, 1);

    if (summary.events.items.len >= MAX_VALIDATE_EVENTS) {
        summary.truncated = true;
        noteTruncatedFixRejection(summary, event.kind);
        return;
    }
    summary.events.append(allocator, event) catch |err| return diag.oom(.validation_events, err);
}

fn noteTruncatedFixRejection(summary: *Summary, kind: Kind) void {
    switch (kind) {
        .invalid_character => {
            if (summary.truncated_fix_rejection == null) {
                summary.truncated_fix_rejection = .invalid_character;
            }
        },
        .no_sequences, .duplicate_name, .null_byte => {
            if (summary.truncated_fix_rejection == null or
                summary.truncated_fix_rejection == .invalid_character)
            {
                summary.truncated_fix_rejection = kind;
            }
        },
        else => {},
    }
}

fn startRecord(
    allocator: std.mem.Allocator,
    summary: *Summary,
    seen_names: *std.StringHashMapUnmanaged(usize),
    state: *RecordState,
    line: []const u8,
    line_number: usize,
    options: Options,
    diag: *Diagnostics,
) diagnostics.Error!void {
    const line_offset = diag.at.location.offset.?;
    const first_column = diag.at.location.column.?;
    state.reset();
    state.active = true;
    state.header_line = line_number;
    summary.header_count = try diag.add(.records, summary.header_count, 1);
    diag.at.record = summary.header_count;

    var name_end: usize = 1;
    while (name_end < line.len and line[name_end] != ' ' and line[name_end] != '\t') : (name_end += 1) {}
    state.name = line[1..name_end];

    const gop = seen_names.getOrPut(allocator, state.name) catch |err| return diag.oom(.validation_identifiers, err);
    if (gop.found_existing) {
        try appendEvent(allocator, summary, diag, .{
            .level = .error_level,
            .kind = .duplicate_name,
            .line = line_number,
            .name = state.name,
            .first_line = gop.value_ptr.*,
        });
    } else {
        gop.value_ptr.* = line_number;
    }

    if (line.len - 1 > options.max_header_len) {
        try appendEvent(allocator, summary, diag, .{
            .level = .warning,
            .kind = .long_header,
            .line = line_number,
            .limit = options.max_header_len,
        });
    }

    if (options.schema != .none and !matchesSchema(options.schema, line)) {
        try appendEvent(allocator, summary, diag, .{
            .level = .warning,
            .kind = .schema_violation,
            .line = line_number,
            .name = state.name,
        });
    }

    for (line, 0..) |byte, column| {
        if (byte == 0) {
            try appendEvent(allocator, summary, diag, .{
                .level = .error_level,
                .kind = .null_byte,
                .offset = line_offset + column,
                .column = first_column + column,
                .control_place = .header,
                .line = line_number,
            });
            break;
        }
    }
}

fn finalizeRecord(allocator: std.mem.Allocator, summary: *Summary, state: *RecordState, diag: *Diagnostics) diagnostics.Error!void {
    if (!state.active) return;

    if (!state.width_warning_emitted) {
        if (state.expected_width) |expected| {
            if (state.pending_width) |final_width| {
                if (final_width > expected) {
                    try appendEvent(allocator, summary, diag, .{
                        .level = .warning,
                        .kind = .inconsistent_line_widths,
                        .line = state.pending_line,
                        .name = state.name,
                        .expected_width = expected,
                        .actual_width = final_width,
                    });
                }
            }
        }
    }

    if (!state.has_bases) {
        try appendEvent(allocator, summary, diag, .{
            .level = .warning,
            .kind = .empty_sequence,
            .line = state.header_line,
            .name = state.name,
        });
        summary.record_widths.append(allocator, DEFAULT_FIX_WIDTH) catch |err| return diag.oom(.line_widths, err);
    } else {
        summary.sequence_count += 1;
        summary.record_widths.append(allocator, modalWidth(state)) catch |err| return diag.oom(.line_widths, err);
    }

    state.reset();
}

fn scanSequenceLine(
    allocator: std.mem.Allocator,
    summary: *Summary,
    state: *RecordState,
    line: []const u8,
    line_number: usize,
    type_counts: *[256]u64,
    type_total: *u64,
    byte_counts: *[256]u64,
    first_byte_at: *[256]diagnostics.Reference,
    diag: *Diagnostics,
) diagnostics.Error!void {
    const line_offset = diag.at.location.offset.?;
    const first_column = diag.at.location.column.?;
    for (line, 0..) |byte, column| {
        if (byte == 0) {
            try appendEvent(allocator, summary, diag, .{
                .level = .error_level,
                .kind = .null_byte,
                .offset = line_offset + column,
                .column = first_column + column,
                .control_place = .sequence,
                .line = line_number,
            });
            break;
        }
    }

    const trimmed = trimSequenceRight(line);
    if (trimmed.len != line.len) {
        try appendEvent(allocator, summary, diag, .{
            .level = .warning,
            .kind = .trailing_whitespace,
            .line = line_number,
            .name = state.name,
        });
    }

    const width = try checkedWidth(trimmed.len, diag);
    if (width > 0) {
        state.has_bases = true;
        try updateWidth(allocator, state, width, diag);
        if (state.pending_width) |pending_width| {
            const expected = state.expected_width orelse pending_width;
            state.expected_width = expected;
            if (!state.width_warning_emitted and pending_width != expected) {
                try appendEvent(allocator, summary, diag, .{
                    .level = .warning,
                    .kind = .inconsistent_line_widths,
                    .line = state.pending_line,
                    .name = state.name,
                    .expected_width = expected,
                    .actual_width = pending_width,
                });
                state.width_warning_emitted = true;
            }
        }
        state.pending_width = width;
        state.pending_line = line_number;
    }

    for (trimmed, 0..) |byte, column| {
        byte_counts[byte] += 1;
        if (first_byte_at[byte].location.line == null) first_byte_at[byte] = .{ .source = .fasta, .location = .{ .line = line_number, .offset = line_offset + column, .column = first_column + column }, .record = summary.header_count };

        if (type_total.* < stats.VALIDATE_TYPE_SAMPLE_BASES) {
            type_counts[byte] += 1;
            type_total.* += 1;
        }
    }
}

fn scanNullsOutsideSequence(
    allocator: std.mem.Allocator,
    summary: *Summary,
    line: []const u8,
    line_number: usize,
    diag: *Diagnostics,
) diagnostics.Error!void {
    const line_offset = diag.at.location.offset.?;
    const first_column = diag.at.location.column.?;
    for (line, 0..) |byte, column| {
        if (byte == 0) {
            try appendEvent(allocator, summary, diag, .{
                .level = .error_level,
                .kind = .null_byte,
                .offset = line_offset + column,
                .column = first_column + column,
                .control_place = .preamble,
                .line = line_number,
            });
            return;
        }
    }
}

fn updateWidth(allocator: std.mem.Allocator, state: *RecordState, width: u32, diag: *Diagnostics) diagnostics.Error!void {
    for (state.widths.items) |*entry| {
        if (entry.width == width) {
            entry.count += 1;
            return;
        }
    }

    state.widths.append(allocator, .{
        .width = width,
        .count = 1,
    }) catch |err| return diag.oom(.line_widths, err);
}

fn modalWidth(state: *const RecordState) u32 {
    var best_width: u32 = DEFAULT_FIX_WIDTH;
    var best_count: usize = 0;
    for (state.widths.items) |entry| {
        if (entry.count > best_count) {
            best_width = entry.width;
            best_count = entry.count;
        }
    }
    return if (best_width == 0) DEFAULT_FIX_WIDTH else best_width;
}

fn appendInvalidCharacterEvents(
    allocator: std.mem.Allocator,
    summary: *Summary,
    options: Options,
    sequence_type: stats.SequenceType,
    byte_counts: [256]u64,
    first_byte_at: [256]diagnostics.Reference,
    diag: *Diagnostics,
) diagnostics.Error!void {
    const alphabet = if (options.custom_alphabet) |custom|
        buildAlphabet(custom)
    else switch (sequence_type) {
        .nucleotide => buildAlphabet("ACGTUNRYWSMKHBVDacgtunrywsmkhbvd"),
        .protein => buildAlphabet("ACDEFGHIKLMNPQRSTVWYBZJXUO*-acdefghiklmnpqrstvwybzjxuo"),
    };

    for (byte_counts, 0..) |count, i| {
        if (count == 0) continue;
        const byte: u8 = @intCast(i);
        if (byte == 0) continue;
        if (!alphabet[byte]) {
            try appendEvent(allocator, summary, diag, .{
                .level = .error_level,
                .kind = .invalid_character,
                .line = @intCast(first_byte_at[byte].location.line.?),
                .offset = first_byte_at[byte].location.offset,
                .column = first_byte_at[byte].location.column,
                .record = first_byte_at[byte].record,
                .byte = byte,
            });
        }
    }
}

fn buildAlphabet(chars: []const u8) [256]bool {
    var table = [_]bool{false} ** 256;
    for (chars) |byte| table[byte] = true;
    return table;
}

fn trimSequenceRight(line: []const u8) []const u8 {
    var end = line.len;
    while (end > 0 and (line[end - 1] == ' ' or line[end - 1] == '\t')) {
        end -= 1;
    }
    return line[0..end];
}

fn countLines(data: []const u8, diag: *Diagnostics) diagnostics.Error!usize {
    if (data.len == 0) return 1;
    var lines: usize = 1;

    var pos: usize = 0;
    while (data.len - pos >= SIMD_CHUNK_SIZE) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const mask = chunk == @as(SimdVec, @splat('\n'));
        lines = try diag.add(.physical_lines, lines, @popCount(@as(u32, @bitCast(mask))));
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) : (pos += 1) {
        if (data[pos] == '\n') lines = try diag.add(.physical_lines, lines, 1);
    }
    return lines;
}

fn findNextNewline(data: []const u8, start: usize) usize {
    var pos = start;
    while (data.len - pos >= SIMD_CHUNK_SIZE) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const mask = chunk == @as(SimdVec, @splat('\n'));
        if (@reduce(.Or, mask)) {
            inline for (0..SIMD_CHUNK_SIZE) |i| {
                if (mask[i]) return pos + i;
            }
        }
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) : (pos += 1) {
        if (data[pos] == '\n') return pos;
    }
    return data.len;
}

fn matchesSchema(schema: Schema, header_line: []const u8) bool {
    if (header_line.len == 0 or header_line[0] != '>') return false;
    const header = header_line[1..];
    return switch (schema) {
        .none => true,
        .uniprot => matchesUniprot(header),
        .refseq => matchesRefseq(header),
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

fn checkedWidth(width: usize, diag: *Diagnostics) diagnostics.Error!u32 {
    return std.math.cast(u32, width) orelse diag.fail(.{ .@"index.sequence_line_too_long" = .{
        .field = .line_bytes,
        .observed = .{ .exact = width },
        .limit = std.math.maxInt(u32),
        .width_role = .validator_width,
    } }, .scan, .fasta, null);
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

fn writeFixed(io: std.Io, data: []const u8, record_widths: []const u32, output_path: []const u8, input: std.Io.File, before: std.Io.File.Stat, failure: *diagnostics.Failure) diagnostics.Error!void {
    var diag = Diagnostics{ .failure = failure };
    var atomic_file = std.Io.Dir.cwd().createFileAtomic(io, output_path, .{ .replace = true }) catch |err| return diag.fail(.{ .@"io.create_failed" = {} }, .create, .output, err);
    defer cleanupFixed(&atomic_file, io, failure);
    failure.outputs.file.state = .staged;
    var out_buf: [65536]u8 = undefined;
    var file_fw = atomic_file.file.writer(io, &out_buf);
    errdefer file_fw.interface.end = 0;
    writeFixedContent(&file_fw.interface, data, record_widths) catch |err| return diag.fail(.{ .@"io.write_failed" = {} }, .write, .output, file_fw.err orelse err);
    file_fw.flush() catch |err| {
        failure.outputs.file.final_flush = .failed;
        return diag.fail(.{ .@"io.flush_failed" = {} }, .flush, .output, file_fw.err orelse err);
    };
    failure.outputs.file.final_flush = .succeeded;
    try checkRewriteSource(io, input, before, failure);
    atomic_file.replace(io) catch |err| return diag.fail(.{ .@"io.publish_failed" = {} }, .publish, .output, err);
    failure.outputs.file.state = .published;
}

const PreamblePosition = struct { offset: usize, line: usize, column: usize };

fn firstPreamble(data: []const u8) ?PreamblePosition {
    var pos: usize = if (std.mem.startsWith(u8, data, UTF8_BOM)) UTF8_BOM.len else 0;
    var line_number: usize = 1;
    while (pos < data.len) {
        const end = findNextNewline(data, pos);
        const line = data[pos..end];
        const first_column: usize = if (line_number == 1) pos + 1 else 1;
        if (line.len > 0 and line[0] == '>') return null;
        for (line, 0..) |byte, column| {
            if (byte != ' ' and byte != '\t' and byte != '\r')
                return .{ .offset = pos + column, .line = line_number, .column = first_column + column };
        }
        if (end == data.len) break;
        pos = end + 1;
        line_number += 1;
    }
    return null;
}

fn checkPreamble(data: []const u8, options: Options, failure: *diagnostics.Failure) diagnostics.Error!void {
    const at = firstPreamble(data) orelse return;
    return failure.fail(.{
        .details = .{ .@"repair.blocked" = .{
            .blocker = .{ .content_before_header = .{ .content_bytes = .{ .at_least = 1 } } },
            .repair_policy = if (options.fix_format_only) .legacy_format_only else .legacy_safe,
        } },
        .stage = .scan,
        .at = .{ .source = .fasta, .location = .{ .offset = at.offset, .line = at.line, .column = at.column } },
    });
}

fn unknownOutputIdentity(failure: *diagnostics.Failure, source: diagnostics.SourceId, cause: anyerror) diagnostics.Error {
    return failure.fail(.{
        .details = .{ .@"repair.output_identity_unknown" = .{ .input = .fasta, .output = .output } },
        .stage = .check_identity,
        .at = .{ .source = source },
        .cause = .{ .zig = cause },
    });
}

fn checkRewriteDestination(io: std.Io, input_file: std.Io.File, input_path: []const u8, output_path: []const u8, failure: *diagnostics.Failure) diagnostics.Error!void {
    if (std.mem.eql(u8, input_path, output_path)) return failure.fail(.{
        .details = .{ .@"repair.output_is_input" = .{ .input = .fasta, .output = .output, .comparison = .lexical_path } },
        .stage = .check_identity,
        .at = .{ .source = .output },
    });
    const output_file = std.Io.Dir.cwd().openFile(io, output_path, .{}) catch |err| {
        if (err == error.FileNotFound) return;
        return unknownOutputIdentity(failure, .output, err);
    };
    defer output_file.close(io);
    var input_buf: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const input_len = input_file.realPath(io, &input_buf) catch |err| return unknownOutputIdentity(failure, .fasta, err);
    var output_buf: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const output_len = output_file.realPath(io, &output_buf) catch |err| return unknownOutputIdentity(failure, .output, err);
    if (std.mem.eql(u8, input_buf[0..input_len], output_buf[0..output_len])) return failure.fail(.{
        .details = .{ .@"repair.output_is_input" = .{ .input = .fasta, .output = .output, .comparison = .canonical_path } },
        .stage = .check_identity,
        .at = .{ .source = .output },
    });
    // Atomic replacement replaces the output directory entry. A distinct hard
    // link is safe: replacing that entry leaves the input entry and inode intact.
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

fn writeFixedContent(writer: *std.Io.Writer, data: []const u8, record_widths: []const u32) !void {
    var pos: usize = if (std.mem.startsWith(u8, data, UTF8_BOM)) UTF8_BOM.len else 0;
    var record_index: usize = 0;
    var current_width: u32 = DEFAULT_FIX_WIDTH;
    var line_pos: u32 = 0;
    var in_record = false;

    while (pos <= data.len) {
        const line_start = pos;
        const line_end = findNextNewline(data, pos);
        if (line_start == data.len and line_end == data.len) break;

        var content_end = line_end;
        if (line_end < data.len and content_end > line_start and data[content_end - 1] == '\r') {
            content_end -= 1;
        }
        const line = data[line_start..content_end];

        if (line.len > 0 and line[0] == '>') {
            if (in_record and line_pos > 0) {
                try writer.writeByte('\n');
            }
            try writer.writeAll(line);
            try writer.writeByte('\n');

            current_width = if (record_index < record_widths.len and record_widths[record_index] > 0)
                record_widths[record_index]
            else
                DEFAULT_FIX_WIDTH;
            record_index += 1;
            line_pos = 0;
            in_record = true;
        } else if (in_record) {
            const trimmed = trimSequenceRight(line);
            for (trimmed) |byte| {
                if (line_pos == current_width) {
                    try writer.writeByte('\n');
                    line_pos = 0;
                }
                try writer.writeByte(byte);
                line_pos += 1;
            }
        }

        if (line_end >= data.len) break;
        pos = line_end + 1;
    }

    if (in_record and line_pos > 0) {
        try writer.writeByte('\n');
    }
}

fn writeSummary(writer: anytype, summary: *const Summary, options: Options) !void {
    switch (options.output_mode) {
        .text => try writeText(writer, summary),
        .json_lines => try writeJsonLines(writer, summary),
        .json_summary => try writeJsonSummary(writer, summary),
    }
}

fn writeText(writer: anytype, summary: *const Summary) !void {
    if (summary.events.items.len == 0) {
        try writer.writeAll("OK: no issues found\n");
        return;
    }

    for (summary.events.items) |event| {
        try writer.print("{s}: line {d}: ", .{ event.level.label(), event.line });
        try writeMessage(writer, event, false);
        try writer.writeByte('\n');
    }
}

fn writeMessage(writer: anytype, event: Event, comptime escape_names: bool) !void {
    switch (event.kind) {
        .no_sequences => try writer.writeAll("no sequences found"),
        .duplicate_name => {
            try writer.writeAll("duplicate name '");
            try writeName(writer, event.name, escape_names);
            try writer.print("' (first seen line {d})", .{event.first_line});
        },
        .invalid_character => try writer.print("invalid sequence character 0x{x:0>2}", .{event.byte}),
        .null_byte => try writer.writeAll("null byte found"),
        .utf8_bom => try writer.writeAll("UTF-8 BOM at start of file"),
        .inconsistent_line_widths => {
            try writer.writeAll("inconsistent line width in '");
            try writeName(writer, event.name, escape_names);
            try writer.print("' (expected {d}, found {d})", .{ event.expected_width, event.actual_width });
        },
        .trailing_whitespace => {
            try writer.writeAll("trailing whitespace on sequence line in '");
            try writeName(writer, event.name, escape_names);
            try writer.writeByte('\'');
        },
        .empty_sequence => {
            try writer.writeAll("empty sequence '");
            try writeName(writer, event.name, escape_names);
            try writer.writeByte('\'');
        },
        .missing_terminal_newline => try writer.writeAll("missing terminal newline"),
        .mixed_line_endings => try writer.writeAll("mixed CRLF and LF line endings"),
        .long_header => try writer.print("header exceeds {d} bytes", .{event.limit}),
        .schema_violation => {
            try writer.writeAll("header for '");
            try writeName(writer, event.name, escape_names);
            try writer.writeAll("' does not match schema");
        },
    }
}

fn writeName(writer: anytype, name: []const u8, comptime escape: bool) !void {
    if (escape) {
        try writeJsonStringBytes(writer, name);
    } else {
        try reporting.writeEscapedBytes(writer, name);
    }
}

fn writeJsonLines(writer: anytype, summary: *const Summary) !void {
    for (summary.events.items) |event| {
        try writeJsonEventObject(writer, event);
        try writer.writeByte('\n');
    }
}

fn writeJsonSummary(writer: anytype, summary: *const Summary) !void {
    const type_str: []const u8 = switch (summary.sequence_type) {
        .nucleotide => "nucleotide",
        .protein => "protein",
    };
    try writer.print(
        "{{\"schema_version\":\"v1\",\"truncated\":{s},\"sequence_type\":\"{s}\",\"type_bases_sampled\":{d},\"type_sample_cap\":{d},\"counts\":{{",
        .{
            if (summary.truncated) "true" else "false",
            type_str,
            summary.type_bases_sampled,
            stats.VALIDATE_TYPE_SAMPLE_BASES,
        },
    );
    for (allKinds(), 0..) |kind, i| {
        if (i > 0) try writer.writeByte(',');
        try writer.print("\"{s}\":{d}", .{ kind.text(), countKind(summary, kind) });
    }
    try writer.writeAll("},\"first_examples\":{");
    var wrote_example = false;
    for (allKinds()) |kind| {
        if (firstKind(summary, kind)) |event| {
            if (wrote_example) try writer.writeByte(',');
            try writer.print("\"{s}\":", .{kind.text()});
            try writeJsonEventObject(writer, event);
            wrote_example = true;
        }
    }
    try writer.writeAll("}}\n");
}

fn writeJsonEventObject(writer: anytype, event: Event) !void {
    try writer.print(
        "{{\"schema_version\":\"v1\",\"level\":\"{s}\",\"line\":{d},\"kind\":\"{s}\",\"message\":\"",
        .{ event.level.text(), event.line, event.kind.text() },
    );
    try writeMessage(writer, event, true);
    try writer.writeByte('"');
    if (event.name.len > 0) {
        try writer.writeAll(",\"name\":\"");
        try writeJsonStringBytes(writer, event.name);
        try writer.writeByte('"');
    }
    if (event.first_line > 0) try writer.print(",\"first_line\":{d}", .{event.first_line});
    if (event.kind == .invalid_character) try writer.print(",\"byte\":{d}", .{event.byte});
    try writer.writeByte('}');
}

// Invalid UTF-8 bytes become `\u00XX` in legacy v1. The result is valid JSON,
// but it cannot distinguish those bytes from the corresponding Unicode text.
fn writeJsonStringBytes(writer: anytype, text: []const u8) !void {
    var i: usize = 0;
    while (i < text.len) {
        const byte = text[i];
        if (byte < 0x80) {
            switch (byte) {
                '"' => try writer.writeAll("\\\""),
                '\\' => try writer.writeAll("\\\\"),
                '\n' => try writer.writeAll("\\n"),
                '\r' => try writer.writeAll("\\r"),
                '\t' => try writer.writeAll("\\t"),
                0...8, 11...12, 14...0x1f => try writer.print("\\u{x:0>4}", .{byte}),
                else => try writer.writeByte(byte),
            }
            i += 1;
            continue;
        }

        const seq_len = std.unicode.utf8ByteSequenceLength(byte) catch {
            try writer.print("\\u{x:0>4}", .{byte});
            i += 1;
            continue;
        };
        if (i + seq_len > text.len) {
            try writer.print("\\u{x:0>4}", .{byte});
            i += 1;
            continue;
        }
        _ = std.unicode.utf8Decode(text[i..][0..seq_len]) catch {
            try writer.print("\\u{x:0>4}", .{byte});
            i += 1;
            continue;
        };
        try writer.writeAll(text[i .. i + seq_len]);
        i += seq_len;
    }
}

fn allKinds() []const Kind {
    return std.meta.tags(Kind);
}

fn countKind(summary: *const Summary, kind: Kind) usize {
    return summary.kind_counts[@intFromEnum(kind)];
}

fn firstKind(summary: *const Summary, kind: Kind) ?Event {
    for (summary.events.items) |event| {
        if (event.kind == kind) return event;
    }
    return null;
}

test "[failure] - [rewrite safety]: keeps original preamble location and identity cause" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    const data = UTF8_BOM ++ "\t \n \r\n  junk\n>seq1\nACGT\n";
    try std.testing.expectError(error.RepairRefused, checkPreamble(data, .{ .fix_format_only = true }, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(@as(?u64, 3), primary.at.location.line);
    try std.testing.expectEqual(@as(?u64, 3), primary.at.location.column);
    try std.testing.expectEqual(@as(?u64, std.mem.indexOf(u8, data, "junk").?), primary.at.location.offset);
    try std.testing.expectEqual(diagnostics.RepairPolicy.legacy_format_only, primary.details.@"repair.blocked".repair_policy);
    failure.reset(.validate);
    try checkPreamble(UTF8_BOM ++ " \t\r\n>seq1\nACGT\n", .{}, &failure);
    try std.testing.expect(failure.primary() == null);

    const Fault = struct {
        fn realPath(_: ?*anyopaque, _: std.Io.File, _: []u8) std.Io.File.RealPathError!usize {
            return error.AccessDenied;
        }
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const input = try tmp.dir.createFile(std.testing.io, "input.fa", .{ .read = true });
    defer input.close(std.testing.io);
    const output = try tmp.dir.createFile(std.testing.io, "output.fa", .{});
    output.close(std.testing.io);
    const output_path = try tmp.dir.realPathFileAlloc(std.testing.io, "output.fa", std.testing.allocator);
    defer std.testing.allocator.free(output_path);
    var vtable = std.testing.io.vtable.*;
    vtable.fileRealPath = Fault.realPath;
    failure.reset(.validate);
    try std.testing.expectError(error.RepairRefused, checkRewriteDestination(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, input, "input.fa", output_path, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"repair.output_identity_unknown", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
}

test "[failure] - [rewrite safety]: observed source change prevents atomic publication" {
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
    const output = try tmp.dir.createFile(std.testing.io, "output.fa", .{});
    try std.Io.File.writeStreamingAll(output, std.testing.io, "previous destination\n");
    output.close(std.testing.io);
    const path = try tmp.dir.realPathFileAlloc(std.testing.io, "output.fa", std.testing.allocator);
    defer std.testing.allocator.free(path);
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    var vtable = std.testing.io.vtable.*;
    vtable.fileStat = Fault.stat;
    try std.testing.expectError(error.InputChanged, writeFixed(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, data, &.{4}, path, input, before, &failure));
    const detail = failure.primary().?.details.@"input.changed";
    try std.testing.expect(detail.changed_fields.size and !detail.changed_fields.mtime);
    try std.testing.expectEqual(@as(?u64, data.len), detail.before.bytes);
    try std.testing.expectEqual(@as(?u64, data.len + 1), detail.after.?.bytes);
    const after = try tmp.dir.readFileAlloc(std.testing.io, "output.fa", std.testing.allocator, .limited(1024));
    defer std.testing.allocator.free(after);
    try std.testing.expectEqualStrings("previous destination\n", after);
}

test "[failure] - [validator counters]: captures limits before arithmetic and allocation" {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    var diag = Diagnostics{ .failure = &failure, .at = .{ .source = .fasta, .location = .{ .offset = 5, .line = 2, .column = 1 }, .record = 1 } };
    var summary = Summary{ .events = .empty, .record_widths = .empty, .error_count = std.math.maxInt(usize) };
    defer summary.deinit(std.testing.allocator);
    try std.testing.expectError(error.CountOverflow, appendEvent(std.testing.allocator, &summary, &diag, .{ .level = .error_level, .kind = .null_byte, .line = 2, .offset = 9, .column = 5 }));
    const detail = failure.primary().?.details.@"resource.count_overflow";
    try std.testing.expectEqual(@as(?u128, std.math.maxInt(usize)), detail.current);
    try std.testing.expectEqual(@as(?u128, 1), detail.increment);
    try std.testing.expectEqual(@as(?u64, 9), failure.primary().?.at.location.offset);
    try std.testing.expectEqual(@as(?u64, 5), diag.at.location.offset);
    try std.testing.expectEqual(@as(usize, 0), summary.events.items.len);
    failure.reset(.validate);
    try std.testing.expectError(error.CountOverflow, diag.add(.physical_lines, std.math.maxInt(usize), 1));
    try std.testing.expectEqual(diagnostics.Counter.physical_lines, failure.primary().?.details.@"resource.count_overflow".counter);
    if (comptime @sizeOf(usize) > @sizeOf(u32)) {
        failure.reset(.validate);
        const width: usize = @as(usize, std.math.maxInt(u32)) + 1;
        try std.testing.expectError(error.InvalidIndex, checkedWidth(width, &diag));
        const limit = failure.primary().?.details.@"index.sequence_line_too_long";
        try std.testing.expectEqual(@as(u64, std.math.maxInt(u32)), limit.limit);
        try std.testing.expectEqual(@as(u64, width), limit.observed.exact);
        try std.testing.expectEqual(diagnostics.LineWidthRole.validator_width, limit.width_role);
    }
}

test "[failure] - [validator output]: preserves publication primary and records failed staged cleanup" {
    const Fault = struct {
        fn rename(_: ?*anyopaque, _: std.Io.Dir, _: []const u8, _: std.Io.Dir, _: []const u8) std.Io.Dir.RenameError!void {
            return error.AccessDenied;
        }
        fn remove(_: ?*anyopaque, _: std.Io.Dir, _: []const u8) std.Io.Dir.DeleteFileError!void {
            return error.AccessDenied;
        }
    };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const input = try tmp.dir.createFile(std.testing.io, "in.fa", .{ .read = true });
    defer input.close(std.testing.io);
    const data = ">seq\nACGT\n";
    try std.Io.File.writeStreamingAll(input, std.testing.io, data);
    const before = try input.stat(std.testing.io);
    const output = try tmp.dir.createFile(std.testing.io, "out.fa", .{});
    try std.Io.File.writeStreamingAll(output, std.testing.io, "old output\n");
    output.close(std.testing.io);
    const path = try tmp.dir.realPathFileAlloc(std.testing.io, "out.fa", std.testing.allocator);
    defer std.testing.allocator.free(path);
    var vtable = std.testing.io.vtable.*;
    vtable.dirRename = Fault.rename;
    vtable.dirDeleteFile = Fault.remove;
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    try std.testing.expectError(error.Io, writeFixed(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, data, &.{4}, path, input, before, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.publish_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(@as(u8, 1), failure.secondary_len);
    try std.testing.expectEqual(diagnostics.Code.@"io.cleanup_failed", failure.secondary[0].code());
    try std.testing.expectEqual(.unknown, failure.outputs.file.state);
    const after = try tmp.dir.readFileAlloc(std.testing.io, "out.fa", std.testing.allocator, .limited(1024));
    defer std.testing.allocator.free(after);
    try std.testing.expectEqualStrings("old output\n", after);
}

test "[failure] - [validator rewrite]: failed flush aborts staging and preserves the destination" {
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
    const output = try tmp.dir.createFile(std.testing.io, "out.fa", .{});
    try std.Io.File.writeStreamingAll(output, std.testing.io, "old output\n");
    output.close(std.testing.io);
    const path = try tmp.dir.realPathFileAlloc(std.testing.io, "out.fa", std.testing.allocator);
    defer std.testing.allocator.free(path);
    var vtable = std.testing.io.vtable.*;
    vtable.fileWritePositional = Fault.write;
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    try std.testing.expectError(error.Io, writeFixed(.{ .userdata = std.testing.io.userdata, .vtable = &vtable }, data, &.{4}, path, input, before, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.InputOutput, failure.primary().?.cause.?.zig);
    try std.testing.expectEqual(diagnostics.SourceId.output, failure.primary().?.at.source);
    try std.testing.expectEqual(.failed, failure.outputs.file.final_flush);
    try std.testing.expectEqual(.aborted, failure.outputs.file.state);
    try std.testing.expectEqual(@as(u8, 0), failure.secondary_len);
    const after = try tmp.dir.readFileAlloc(std.testing.io, "out.fa", std.testing.allocator, .limited(1024));
    defer std.testing.allocator.free(after);
    try std.testing.expectEqualStrings("old output\n", after);
}
