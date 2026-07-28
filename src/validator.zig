//! FASTA structure, alphabet, and header validation with optional fix streaming.
//!
//! Events are capped (`max_validate_events`). JSON output escapes arbitrary header bytes.
//! Shared sequence-type sampling uses `stats.detectType`.

const std = @import("std");
const index_format = @import("index_format.zig");
const platform = @import("platform.zig");
const stats = @import("stats.zig");

const printErrorAndExit = index_format.printErrorAndExit;

const UTF8_BOM = "\xEF\xBB\xBF";
const DEFAULT_MAX_HEADER_LEN = 1024;
const DEFAULT_FIX_WIDTH = 60;
const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);

/// Hard cap on retained validate events. Further issues set `Summary.truncated`
/// instead of growing an unbounded list.
pub const max_validate_events: usize = 10_000;

comptime {
    // Keep the validate help line in src/main.zig synchronized with this value.
    if (max_validate_events != 10_000) {
        @compileError("update validate help text for max_validate_events");
    }
}

pub const OutputMode = enum {
    text,
    json_lines,
    json_summary,
};

pub const Schema = enum {
    none,
    uniprot,
    refseq,
};

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

pub const ValidateEvent = struct {
    level: Level,
    kind: Kind,
    line: usize,
    name: []const u8 = "",
    first_line: usize = 0,
    byte: u8 = 0,
    expected_width: u32 = 0,
    actual_width: u32 = 0,
    limit: usize = 0,
};

pub const Summary = struct {
    events: std.ArrayList(ValidateEvent),
    record_widths: std.ArrayList(u32),
    sequence_type: stats.SequenceType = .nucleotide,
    /// Bases fed into `detectType` (capped at `stats.validate_type_sample_bases`).
    type_bases_sampled: u64 = 0,
    sequence_count: usize = 0,
    header_count: usize = 0,
    error_count: usize = 0,
    warning_count: usize = 0,
    /// True when more issues were found after `max_validate_events` were retained.
    truncated: bool = false,

    pub fn deinit(self: *Summary, allocator: std.mem.Allocator) void {
        self.events.deinit(allocator);
        self.record_widths.deinit(allocator);
    }
};

const WidthCount = struct {
    width: u32,
    count: u32,
    first_order: usize,
};

const RecordState = struct {
    active: bool = false,
    name: []const u8 = "",
    header_line: usize = 0,
    bases: u64 = 0,
    expected_width: ?u32 = null,
    pending_width: ?u32 = null,
    pending_line: usize = 0,
    width_warning_emitted: bool = false,
    widths: std.ArrayList(WidthCount) = .empty,
    width_order: usize = 0,

    fn reset(self: *RecordState) void {
        self.active = false;
        self.name = "";
        self.header_line = 0;
        self.bases = 0;
        self.expected_width = null;
        self.pending_width = null;
        self.pending_line = 0;
        self.width_warning_emitted = false;
        self.widths.clearRetainingCapacity();
        self.width_order = 0;
    }
};

pub fn runValidate(io: std.Io, fasta_path: []const u8, options: Options) void {
    // Mapping ownership matches LoadedIndex: one FileView (or empty), destroyed here.
    // Validation does not load an index; GET/stats use LoadedIndex.deinit(io) instead.
    const file = std.Io.Dir.cwd().openFile(io, fasta_path, .{}) catch |err| switch (err) {
        error.FileNotFound => printErrorAndExit("error: file not found: {s}\n", .{fasta_path}),
        error.AccessDenied => printErrorAndExit("error: access denied: {s}\n", .{fasta_path}),
        else => printErrorAndExit("error: failed to open file: {s}\n", .{fasta_path}),
    };
    defer file.close(io);

    const stat = file.stat(io) catch {
        printErrorAndExit("error: failed to stat file: {s}\n", .{fasta_path});
    };

    var fasta_view: ?platform.FileView = null;
    defer if (fasta_view) |*view| view.destroy(io);

    const data: []const u8 = if (stat.size == 0)
        &[_]u8{}
    else blk: {
        fasta_view = platform.FileView.mapFile(io, file, @intCast(stat.size)) catch {
            printErrorAndExit("error: failed to mmap file: {s}\n", .{fasta_path});
        };
        break :blk fasta_view.?.bytes();
    };

    if (data.len > 0) {
        platform.advise(data, .sequential);
    }

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var summary = validateData(allocator, data, options) catch |err| switch (err) {
        error.OutOfMemory => printErrorAndExit("error: out of memory\n", .{}),
    };
    defer summary.deinit(allocator);

    if (options.fix) {
        const output_path = options.output_path orelse {
            printErrorAndExit("error: validate --fix requires -o <output.fa>\n", .{});
        };
        if (std.mem.eql(u8, fasta_path, output_path)) {
            printErrorAndExit("error: validate --fix will not overwrite input: {s}\n", .{fasta_path});
        }
        ensureFixAllowed(&summary, options);
        writeFixed(io, data, summary.record_widths.items, output_path) catch |err| switch (err) {
            else => printErrorAndExit("error: failed to write fixed FASTA: {s}\n", .{output_path}),
        };
    }

    var out_buf: [65536]u8 = undefined;
    var stdout_fw = std.Io.File.Writer.initStreaming(.stdout(), io, &out_buf);
    writeSummary(allocator, &stdout_fw.interface, &summary, options) catch {
        stdout_fw.flush() catch {};
        printErrorAndExit("error: write failed\n", .{});
    };
    stdout_fw.flush() catch {
        printErrorAndExit("error: write failed\n", .{});
    };

    if (summary.truncated) {
        printErrorAndExit(
            "error: validate stopped after {d} events; fix reported issues and re-run\n",
            .{max_validate_events},
        );
    }

    std.process.exit(exitCodeForOptions(&summary, options));
}

pub fn validateData(allocator: std.mem.Allocator, data: []const u8, options: Options) !Summary {
    var summary = Summary{
        .events = .empty,
        .record_widths = .empty,
    };
    errdefer summary.deinit(allocator);

    var seen_names = std.StringHashMap(usize).init(allocator);
    defer seen_names.deinit();

    var type_counts: [256]u64 = .{0} ** 256;
    var type_total: u64 = 0;
    var byte_counts: [256]u64 = .{0} ** 256;
    var first_byte_line: [256]usize = .{0} ** 256;

    var state = RecordState{};
    defer state.widths.deinit(allocator);

    if (std.mem.startsWith(u8, data, UTF8_BOM)) {
        try appendEvent(allocator, &summary, .{
            .level = .warning,
            .kind = .utf8_bom,
            .line = 1,
        });
    }

    if (data.len == 0 or data[data.len - 1] != '\n') {
        try appendEvent(allocator, &summary, .{
            .level = .warning,
            .kind = .missing_terminal_newline,
            .line = countLines(data),
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
        if (line.len > 0 and line[0] == '>') {
            try finalizeRecord(allocator, &summary, &state);
            try startRecord(allocator, &summary, &seen_names, &state, line, line_number, options);
        } else if (state.active) {
            try scanSequenceLine(allocator, &summary, &state, line, line_number, &type_counts, &type_total, &byte_counts, &first_byte_line);
        } else {
            try scanNullsOutsideSequence(allocator, &summary, line, line_number);
        }

        if (!has_lf) break;
        pos = line_end + 1;
        line_number += 1;
    }

    try finalizeRecord(allocator, &summary, &state);

    if (summary.header_count == 0) {
        try appendEvent(allocator, &summary, .{
            .level = .error_level,
            .kind = .no_sequences,
            .line = 1,
        });
    }

    if (saw_crlf and saw_lf) {
        try appendEvent(allocator, &summary, .{
            .level = .warning,
            .kind = .mixed_line_endings,
            .line = @max(first_crlf_line, first_lf_line),
        });
    }

    summary.sequence_type = stats.detectType(&type_counts, type_total);
    summary.type_bases_sampled = type_total;
    try appendInvalidCharacterEvents(allocator, &summary, options, summary.sequence_type, byte_counts, first_byte_line);
    return summary;
}

pub fn exitCode(error_count: usize, warning_count: usize, strict: bool) u8 {
    if (error_count > 0) return 1;
    if (warning_count > 0) return if (strict) 1 else 2;
    return 0;
}

pub fn exitCodeForOptions(summary: *const Summary, options: Options) u8 {
    if (!options.fix) return exitCode(summary.error_count, summary.warning_count, options.strict);

    var remaining_errors: usize = 0;
    var remaining_warnings: usize = 0;
    for (summary.events.items) |event| {
        if (isFixedByFormatRewrite(event.kind)) continue;
        switch (event.level) {
            .error_level => remaining_errors += 1,
            .warning => remaining_warnings += 1,
        }
    }
    return exitCode(remaining_errors, remaining_warnings, options.strict);
}

fn isFixedByFormatRewrite(kind: Kind) bool {
    // --fix rewrites only these kinds. Everything else is left unchanged or blocks fix.
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

pub fn fixRejection(summary: *const Summary, options: Options) ?Kind {
    for (summary.events.items) |event| {
        switch (event.kind) {
            .invalid_character => if (!options.fix_format_only) return .invalid_character,
            .no_sequences, .duplicate_name, .null_byte => return event.kind,
            else => {},
        }
    }
    return null;
}

fn appendEvent(allocator: std.mem.Allocator, summary: *Summary, event: ValidateEvent) !void {
    if (summary.events.items.len >= max_validate_events) {
        summary.truncated = true;
        return;
    }
    try summary.events.append(allocator, event);
    switch (event.level) {
        .error_level => summary.error_count += 1,
        .warning => summary.warning_count += 1,
    }
}

fn startRecord(
    allocator: std.mem.Allocator,
    summary: *Summary,
    seen_names: *std.StringHashMap(usize),
    state: *RecordState,
    line: []const u8,
    line_number: usize,
    options: Options,
) !void {
    state.reset();
    state.active = true;
    state.header_line = line_number;
    summary.header_count += 1;

    var name_end: usize = 1;
    while (name_end < line.len and line[name_end] != ' ' and line[name_end] != '\t') : (name_end += 1) {}
    state.name = line[1..name_end];

    const gop = try seen_names.getOrPut(state.name);
    if (gop.found_existing) {
        try appendEvent(allocator, summary, .{
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
        try appendEvent(allocator, summary, .{
            .level = .warning,
            .kind = .long_header,
            .line = line_number,
            .limit = options.max_header_len,
        });
    }

    if (options.schema != .none and !matchesSchema(options.schema, line)) {
        try appendEvent(allocator, summary, .{
            .level = .warning,
            .kind = .schema_violation,
            .line = line_number,
            .name = state.name,
        });
    }

    for (line) |byte| {
        if (byte == 0) {
            try appendEvent(allocator, summary, .{
                .level = .error_level,
                .kind = .null_byte,
                .line = line_number,
            });
            break;
        }
    }
}

fn finalizeRecord(allocator: std.mem.Allocator, summary: *Summary, state: *RecordState) !void {
    if (!state.active) return;

    if (state.bases == 0) {
        try appendEvent(allocator, summary, .{
            .level = .warning,
            .kind = .empty_sequence,
            .line = state.header_line,
            .name = state.name,
        });
        try summary.record_widths.append(allocator, DEFAULT_FIX_WIDTH);
    } else {
        summary.sequence_count += 1;
        try summary.record_widths.append(allocator, modalWidth(state));
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
    first_byte_line: *[256]usize,
) !void {
    for (line) |byte| {
        if (byte == 0) {
            try appendEvent(allocator, summary, .{
                .level = .error_level,
                .kind = .null_byte,
                .line = line_number,
            });
            break;
        }
    }

    const trimmed = trimSequenceRight(line);
    if (trimmed.len != line.len) {
        try appendEvent(allocator, summary, .{
            .level = .warning,
            .kind = .trailing_whitespace,
            .line = line_number,
            .name = state.name,
        });
    }

    const width: u32 = @intCast(trimmed.len);
    if (width > 0) {
        try updateWidth(allocator, state, width);
        if (state.pending_width) |pending_width| {
            const expected = state.expected_width orelse pending_width;
            state.expected_width = expected;
            if (!state.width_warning_emitted and pending_width != expected) {
                try appendEvent(allocator, summary, .{
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

    state.bases += width;
    for (trimmed) |byte| {
        byte_counts[byte] += 1;
        if (first_byte_line[byte] == 0) first_byte_line[byte] = line_number;

        if (type_total.* < stats.validate_type_sample_bases) {
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
) !void {
    for (line) |byte| {
        if (byte == 0) {
            try appendEvent(allocator, summary, .{
                .level = .error_level,
                .kind = .null_byte,
                .line = line_number,
            });
            return;
        }
    }
}

fn updateWidth(allocator: std.mem.Allocator, state: *RecordState, width: u32) !void {
    for (state.widths.items) |*entry| {
        if (entry.width == width) {
            entry.count += 1;
            return;
        }
    }

    try state.widths.append(allocator, .{
        .width = width,
        .count = 1,
        .first_order = state.width_order,
    });
    state.width_order += 1;
}

fn modalWidth(state: *const RecordState) u32 {
    var best_width: u32 = DEFAULT_FIX_WIDTH;
    var best_count: u32 = 0;
    var best_order: usize = std.math.maxInt(usize);
    for (state.widths.items) |entry| {
        if (entry.count > best_count or (entry.count == best_count and entry.first_order < best_order)) {
            best_width = entry.width;
            best_count = entry.count;
            best_order = entry.first_order;
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
    first_byte_line: [256]usize,
) !void {
    const alphabet = if (options.custom_alphabet) |custom|
        buildAlphabet(custom)
    else switch (sequence_type) {
        .nucleotide => buildAlphabet("ACGTUNRYWSMKHBVDacgtunrywsmkhbvd"),
        .protein => buildAlphabet("ACDEFGHIKLMNPQRSTVWYUOX*-acdefghiklmnpqrstvwyuox"),
    };

    for (byte_counts, 0..) |count, i| {
        if (count == 0) continue;
        const byte: u8 = @intCast(i);
        if (byte == 0) continue;
        if (!alphabet[byte]) {
            try appendEvent(allocator, summary, .{
                .level = .error_level,
                .kind = .invalid_character,
                .line = first_byte_line[byte],
                .byte = byte,
            });
        }
    }
}

fn buildAlphabet(chars: ?[]const u8) [256]bool {
    var table = [_]bool{false} ** 256;
    if (chars) |bytes| {
        for (bytes) |byte| table[byte] = true;
    }
    return table;
}

fn trimSequenceRight(line: []const u8) []const u8 {
    var end = line.len;
    while (end > 0 and (line[end - 1] == ' ' or line[end - 1] == '\t')) {
        end -= 1;
    }
    return line[0..end];
}

fn countLines(data: []const u8) usize {
    if (data.len == 0) return 1;
    var lines: usize = 1;

    var pos: usize = 0;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const mask = chunk == @as(SimdVec, @splat('\n'));
        lines += @popCount(@as(u32, @bitCast(mask)));
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) : (pos += 1) {
        if (data[pos] == '\n') lines += 1;
    }
    return lines;
}

fn findNextNewline(data: []const u8, start: usize) usize {
    var pos = start;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
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

fn ensureFixAllowed(summary: *const Summary, options: Options) void {
    const rejection = fixRejection(summary, options) orelse return;
    switch (rejection) {
        .invalid_character => {
            printErrorAndExit("error: validate --fix refuses character-level errors; use --fix-format-only to keep them unchanged\n", .{});
        },
        else => printErrorAndExit("error: validate --fix cannot repair {s}\n", .{rejection.text()}),
    }
}

/// Rewrite format-level issues to normalized LF FASTA bytes. Caller must run validateData first.
/// Allocates the full result for tests and library callers; CLI `--fix` streams via `writeFixed`.
pub fn fixData(allocator: std.mem.Allocator, data: []const u8, record_widths: []const u32) ![]u8 {
    var aw: std.Io.Writer.Allocating = .init(allocator);
    errdefer aw.deinit();
    try writeFixedContent(&aw.writer, data, record_widths);
    return try aw.toOwnedSlice();
}

fn writeFixed(io: std.Io, data: []const u8, record_widths: []const u32, output_path: []const u8) !void {
    const cwd = std.Io.Dir.cwd();
    const out_file = try cwd.createFile(io, output_path, .{ .truncate = true });
    errdefer cwd.deleteFile(io, output_path) catch {};
    defer out_file.close(io);

    var out_buf: [65536]u8 = undefined;
    var file_fw = out_file.writer(io, &out_buf);
    try writeFixedContent(&file_fw.interface, data, record_widths);
    try file_fw.flush();
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
    } else if (in_record and (data.len == 0 or data[data.len - 1] != '\n')) {
        try writer.writeByte('\n');
    }
}

fn writeSummary(allocator: std.mem.Allocator, writer: anytype, summary: *const Summary, options: Options) !void {
    switch (options.output_mode) {
        .text => try writeText(writer, summary),
        .json_lines => try writeJsonLines(allocator, writer, summary),
        .json_summary => try writeJsonSummary(allocator, writer, summary),
    }
}

fn writeText(writer: anytype, summary: *const Summary) !void {
    if (summary.events.items.len == 0) {
        try writer.writeAll("OK: no issues found\n");
        return;
    }

    for (summary.events.items) |event| {
        try writer.print("{s}: line {d}: ", .{ event.level.label(), event.line });
        try writeHumanMessage(writer, event);
        try writer.writeByte('\n');
    }
}

fn writeHumanMessage(writer: anytype, event: ValidateEvent) !void {
    switch (event.kind) {
        .no_sequences => try writer.writeAll("no sequences found"),
        .duplicate_name => try writer.print("duplicate name '{s}' (first seen line {d})", .{ event.name, event.first_line }),
        .invalid_character => try writer.print("invalid sequence character 0x{x:0>2}", .{event.byte}),
        .null_byte => try writer.writeAll("null byte found"),
        .utf8_bom => try writer.writeAll("UTF-8 BOM at start of file"),
        .inconsistent_line_widths => try writer.print("inconsistent line width in '{s}' (expected {d}, found {d})", .{ event.name, event.expected_width, event.actual_width }),
        .trailing_whitespace => try writer.print("trailing whitespace on sequence line in '{s}'", .{event.name}),
        .empty_sequence => try writer.print("empty sequence '{s}'", .{event.name}),
        .missing_terminal_newline => try writer.writeAll("missing terminal newline"),
        .mixed_line_endings => try writer.writeAll("mixed CRLF and LF line endings"),
        .long_header => try writer.print("header exceeds {d} bytes", .{event.limit}),
        .schema_violation => try writer.print("header for '{s}' does not match schema", .{event.name}),
    }
}

fn writeJsonLines(allocator: std.mem.Allocator, writer: anytype, summary: *const Summary) !void {
    for (summary.events.items) |event| {
        try writeJsonEventObject(allocator, writer, event);
        try writer.writeByte('\n');
    }
}

fn writeJsonSummary(allocator: std.mem.Allocator, writer: anytype, summary: *const Summary) !void {
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
            stats.validate_type_sample_bases,
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
            try writeJsonEventObject(allocator, writer, event);
            wrote_example = true;
        }
    }
    try writer.writeAll("}}\n");
}

fn writeJsonEventObject(allocator: std.mem.Allocator, writer: anytype, event: ValidateEvent) !void {
    try writer.print(
        "{{\"schema_version\":\"v1\",\"level\":\"{s}\",\"line\":{d},\"kind\":\"{s}\",\"message\":\"",
        .{ event.level.text(), event.line, event.kind.text() },
    );
    try writeJsonMessage(allocator, writer, event);
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

fn writeJsonMessage(allocator: std.mem.Allocator, writer: anytype, event: ValidateEvent) !void {
    var aw: std.Io.Writer.Allocating = .init(allocator);
    defer aw.deinit();
    try writeHumanMessage(&aw.writer, event);
    try writeJsonStringBytes(writer, aw.written());
}

/// Write `text` inside a JSON string. ASCII controls and quotes are escaped.
/// Invalid UTF-8 bytes become `\u00XX` so the overall document stays valid UTF-8 JSON.
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

/// Render one JSON event object (no trailing newline). For tests and tooling.
pub fn renderJsonEvent(allocator: std.mem.Allocator, event: ValidateEvent) ![]u8 {
    var aw: std.Io.Writer.Allocating = .init(allocator);
    errdefer aw.deinit();
    try writeJsonEventObject(allocator, &aw.writer, event);
    return try aw.toOwnedSlice();
}

fn allKinds() []const Kind {
    return &.{
        .no_sequences,
        .duplicate_name,
        .invalid_character,
        .null_byte,
        .utf8_bom,
        .inconsistent_line_widths,
        .trailing_whitespace,
        .empty_sequence,
        .missing_terminal_newline,
        .mixed_line_endings,
        .long_header,
        .schema_violation,
    };
}

fn countKind(summary: *const Summary, kind: Kind) usize {
    var count: usize = 0;
    for (summary.events.items) |event| {
        if (event.kind == kind) count += 1;
    }
    return count;
}

fn firstKind(summary: *const Summary, kind: Kind) ?ValidateEvent {
    for (summary.events.items) |event| {
        if (event.kind == kind) return event;
    }
    return null;
}
