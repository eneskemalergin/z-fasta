//! Validate tests: finding codes and counts, read splits, the vector path, `--fix`, and the CLI.

const std = @import("std");
const main = @import("main");
const utility = @import("utility.zig");

const io = std.testing.io;
const validator = main.validator;
const diagnostics = main.index_format.diagnostics;
const reporting = main.reporting;

const ZFASTA_BIN = utility.ZFASTA_BIN;
const expectCliResult = utility.expectCliResult;
const readTestFile = utility.readRequiredFile;
const writeFastaArtifact = utility.writeFastaArtifact;
const writeZfi = utility.writeZfi;
const captureExtractRegion = utility.captureExtractRegion;

const Expected = struct {
    code: diagnostics.Code,
    occurrences: u64,
    records: u64,
    offset: ?u64,
    line: ?u64,
    column: ?u64,
    record: ?u64,
};

fn reportWith(allocator: std.mem.Allocator, data: []const u8, buffer_len: usize, rules: validator.Rules) !*validator.Report {
    const report = try allocator.create(validator.Report);
    errdefer allocator.destroy(report);
    const buffer = try allocator.alloc(u8, buffer_len);
    defer allocator.free(buffer);
    var reader = std.Io.Reader.fixed(data);
    var failure: diagnostics.Failure = undefined;
    try validator.validate(allocator, &reader, buffer, null, rules, report, &failure);
    return report;
}

fn reportBytes(allocator: std.mem.Allocator, data: []const u8, rules: validator.Rules) !*validator.Report {
    return reportWith(allocator, data, validator.READ_BUFFER_BYTES, rules);
}

fn repairBytes(allocator: std.mem.Allocator, data: []const u8, width: u64) ![]u8 {
    const report = try allocator.create(validator.Report);
    defer allocator.destroy(report);
    const buffer = try allocator.alloc(u8, 4096);
    defer allocator.free(buffer);
    var reader = std.Io.Reader.fixed(data);
    var output: std.Io.Writer.Allocating = .init(allocator);
    errdefer output.deinit();
    var failure: diagnostics.Failure = undefined;
    try validator.repair(allocator, &reader, buffer, null, .{}, width, &output.writer, report, &failure);
    return output.toOwnedSlice();
}

fn renderJson(allocator: std.mem.Allocator, report: *const validator.Report, rules: validator.Rules) ![]u8 {
    var failure: diagnostics.Failure = undefined;
    failure.reset(.validate);
    failure.succeed();
    const view = report.validation(rules, false);
    var output: std.Io.Writer.Allocating = .init(allocator);
    errdefer output.deinit();
    try reporting.writeJson(&output.writer, &failure, .{ .validation = &view });
    return output.toOwnedSlice();
}

fn expectFindings(report: *const validator.Report, expected: []const Expected) !void {
    var occurred: usize = 0;
    for (report.findings) |finding| {
        if (finding.occurrences != 0) occurred += 1;
    }
    try std.testing.expectEqual(expected.len, occurred);
    for (expected) |want| {
        const finding = for (report.findings) |finding| {
            if (finding.code == want.code) break finding;
        } else return error.TestUnexpectedResult;
        errdefer std.debug.print("code {s}\n", .{@tagName(want.code)});
        try std.testing.expectEqual(want.occurrences, finding.occurrences);
        try std.testing.expectEqual(want.records, finding.records);
        const first = finding.first.?;
        try std.testing.expectEqual(want.code, std.meta.activeTag(first.details));
        try std.testing.expectEqual(want.offset, first.at.location.offset);
        try std.testing.expectEqual(want.line, first.at.location.line);
        try std.testing.expectEqual(want.column, first.at.location.column);
        try std.testing.expectEqual(want.record, first.at.record);
    }
}

// Alphabet policies may differ.
fn expectSameFindings(a: *const validator.Report, b: *const validator.Report) !void {
    try std.testing.expectEqual(a.completion, b.completion);
    try std.testing.expectEqual(a.examined, b.examined);
    try std.testing.expectEqual(a.records_with_errors, b.records_with_errors);
    try std.testing.expectEqual(a.records_with_warnings, b.records_with_warnings);
    try std.testing.expectEqualSlices(u64, &a.invalid_symbols, &b.invalid_symbols);
    for (a.findings, b.findings) |left, right| {
        errdefer std.debug.print("code {s}\n", .{@tagName(left.code)});
        try std.testing.expectEqual(left.occurrences, right.occurrences);
        try std.testing.expectEqual(left.records, right.records);
        try std.testing.expectEqual(left.first == null, right.first == null);
        const first = left.first orelse continue;
        const other = right.first.?;
        try std.testing.expectEqual(first.at, other.at);
        try std.testing.expectEqual(first.related, other.related);
        try std.testing.expectEqual(first.identifier == null, other.identifier == null);
        if (first.identifier) |identifier| try std.testing.expectEqualSlices(u8, identifier.slice(), other.identifier.?.slice());
        switch (first.details) {
            .@"sequence.invalid_symbol" => |facts| try std.testing.expectEqual(facts.byte, other.details.@"sequence.invalid_symbol".byte),
            else => try std.testing.expectEqualDeep(first.details, other.details),
        }
    }
}

// The clean lead exceeds the alphabet sample, so the defects after it go through the vector path.
fn largeInput(allocator: std.mem.Allocator, random: std.Random, protein: bool) ![]u8 {
    const alphabet = if (protein) "ACDEFGHIKLMNPQRSTVWY" else "ACGTNacgtn";
    const odd = "RYKMSWBDHVU?@1.*-\x00\x01\x7f\xff \t\r";
    const eol: []const u8 = if (random.boolean()) "\n" else "\r\n";
    const width = random.intRangeAtMost(usize, 50, 80);
    var out: std.ArrayList(u8) = .empty;
    errdefer out.deinit(allocator);
    try out.appendSlice(allocator, ">lead clean");
    try out.appendSlice(allocator, eol);
    for (0..110_000 / width + 1) |_| {
        for (0..width) |_| try out.append(allocator, alphabet[random.uintLessThan(usize, alphabet.len)]);
        try out.appendSlice(allocator, eol);
    }
    for (0..random.intRangeAtMost(usize, 30, 60)) |_| {
        try out.print(allocator, ">{s}r{d}{s}", .{ if (random.uintLessThan(u8, 20) == 0) " " else "", random.uintLessThan(u8, 20), if (random.uintLessThan(u8, 20) == 0) "\x01" else " desc" });
        try out.appendSlice(allocator, if (random.uintLessThan(u8, 30) == 0) "\r\n" else eol);
        const lines = random.uintLessThan(usize, 90);
        const record_width = if (random.uintLessThan(u8, 4) == 0) random.intRangeAtMost(usize, 1, 12) else width;
        for (0..lines) |k| {
            const len = if (k + 1 < lines) record_width else random.intRangeAtMost(usize, 1, record_width + 2);
            if (random.uintLessThan(u8, 100) == 0) {
                try out.appendSlice(allocator, eol);
                continue;
            }
            const start = out.items.len;
            for (0..len) |_| try out.append(allocator, alphabet[random.uintLessThan(usize, alphabet.len)]);
            if (random.uintLessThan(u8, 25) == 0) out.items[start + random.uintLessThan(usize, len)] = odd[random.uintLessThan(usize, odd.len)];
            try out.appendSlice(allocator, if (random.uintLessThan(u8, 100) == 0) "\r\n" else eol);
        }
    }
    if (random.uintLessThan(u8, 5) == 0) _ = out.pop();
    return out.toOwnedSlice(allocator);
}

const FailingSource = struct {
    data: []const u8,
    fail_at: usize,
    position: usize = 0,
    reader: std.Io.Reader,

    fn init(data: []const u8, fail_at: usize) FailingSource {
        return .{ .data = data, .fail_at = fail_at, .reader = .{ .vtable = &.{ .stream = stream }, .buffer = &.{}, .seek = 0, .end = 0 } };
    }

    fn stream(reader: *std.Io.Reader, writer: *std.Io.Writer, limit: std.Io.Limit) std.Io.Reader.StreamError!usize {
        const self: *FailingSource = @fieldParentPtr("reader", reader);
        if (self.position == self.fail_at) return error.ReadFailed;
        const stop = @min(self.fail_at, self.data.len);
        if (self.position == stop) return error.EndOfStream;
        const rest = self.data[self.position..stop];
        const n = try writer.write(rest[0..limit.minInt(rest.len)]);
        self.position += n;
        return n;
    }
};

fn fuzzValidator(_: void, smith: *std.testing.Smith) !void {
    var storage: [512]u8 = undefined;
    try checkValidatorProperties(storage[0..smith.slice(&storage)]);
}

fn checkValidatorProperties(data: []const u8) !void {
    const allocator = std.testing.allocator;
    const report = try reportBytes(allocator, data, .{});
    defer allocator.destroy(report);
    try std.testing.expectEqual(reporting.Validation.Completion.complete, report.completion);
    try std.testing.expectEqual(@as(u64, data.len), report.examined.bytes);
    try std.testing.expectEqual(report.examined.records, report.examined.records_completed);
    const fixed = repairBytes(allocator, data, 7) catch |err| switch (err) {
        error.RepairRefused => return,
        else => return err,
    };
    defer allocator.free(fixed);
    const again = try reportBytes(allocator, fixed, .{});
    defer allocator.destroy(again);
    for (again.findings) |finding| {
        if (finding.occurrences != 0) try std.testing.expect(diagnostics.metadata(finding.code).category != .format);
    }
    const twice = try repairBytes(allocator, fixed, 7);
    defer allocator.free(twice);
    try std.testing.expectEqualStrings(fixed, twice);
}

fn checkReportAllocations(allocator: std.mem.Allocator, data: []const u8) !void {
    const report = try std.testing.allocator.create(validator.Report);
    defer std.testing.allocator.destroy(report);
    var buffer: [64]u8 = undefined;
    var reader = std.Io.Reader.fixed(data);
    var failure: diagnostics.Failure = undefined;
    validator.validate(allocator, &reader, &buffer, null, .{}, report, &failure) catch |err| {
        try std.testing.expectEqual(error.OutOfMemory, err);
        try std.testing.expectEqual(reporting.Validation.Completion.interrupted, report.completion);
        try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
        try std.testing.expect(report.examined.bytes <= data.len);
        try std.testing.expect(report.examined.records_completed <= report.examined.records);
        return error.OutOfMemory;
    };
    try std.testing.expectEqual(reporting.Validation.Completion.complete, report.completion);
}

fn scanZfi(allocator: std.mem.Allocator, data: []const u8) !main.indexer.ZfiIndex {
    return main.indexer.scanZfiData(allocator, data);
}

fn expectFixPreservesRegion(allocator: std.mem.Allocator, stem: []const u8, broken: []const u8, region: []const u8, expected: []const u8) !void {
    var broken_index = try scanZfi(allocator, broken);
    defer broken_index.deinit(allocator);
    try std.testing.expect(broken_index.side_tables.items.len > 0);
    const fixed = try repairBytes(allocator, broken, 4);
    defer allocator.free(fixed);
    var fixed_index = try scanZfi(allocator, fixed);
    defer fixed_index.deinit(allocator);
    for (fixed_index.records.items) |record| try std.testing.expect(record.isUniformWidth());
    try std.testing.expectEqual(@as(usize, 0), fixed_index.side_tables.items.len);

    const broken_stem = try std.fmt.allocPrint(allocator, "{s}-broken", .{stem});
    const fixed_stem = try std.fmt.allocPrint(allocator, "{s}-fixed", .{stem});
    const broken_path = try writeFastaArtifact(allocator, broken_stem, broken);
    defer std.Io.Dir.cwd().deleteFile(io, broken_path) catch {};
    const fixed_path = try writeFastaArtifact(allocator, fixed_stem, fixed);
    defer std.Io.Dir.cwd().deleteFile(io, fixed_path) catch {};
    const broken_zfi = try std.fmt.allocPrint(allocator, "{s}.zfi", .{broken_path});
    defer std.Io.Dir.cwd().deleteFile(io, broken_zfi) catch {};
    const fixed_zfi = try std.fmt.allocPrint(allocator, "{s}.zfi", .{fixed_path});
    defer std.Io.Dir.cwd().deleteFile(io, fixed_zfi) catch {};
    try writeZfi(allocator, broken_path, broken);
    try writeZfi(allocator, fixed_path, fixed);
    const broken_output = try captureExtractRegion(allocator, broken_path, region);
    const fixed_output = try captureExtractRegion(allocator, fixed_path, region);
    try std.testing.expectEqualStrings(broken_output, fixed_output);
    try std.testing.expectEqualStrings(expected, fixed_output);
}

// --- Counts and first occurrences ---

test "[unit] - [validate report]: counts every invalid byte" {
    const report = try reportBytes(std.testing.allocator, ">s\nACGT????\n", .{});
    defer std.testing.allocator.destroy(report);
    try expectFindings(report, &.{
        .{ .code = .@"sequence.invalid_symbol", .occurrences = 4, .records = 1, .offset = 7, .line = 2, .column = 5, .record = 1 },
    });
    try std.testing.expectEqual(@as(u64, 4), report.invalid_symbols['?']);
    try std.testing.expectEqual(diagnostics.Alphabet.protein, report.alphabet.alphabet);
    try std.testing.expectEqual(@as(u8, 1), report.status(false, false));
}

test "[unit] - [validate report]: counts structure and identifier codes in their units" {
    const data = "\xEF\xBB\xBFjunk\n>a desc\nACGT\n>a\x01x\r\n>a\nAC\n>\n>NC_1 a long description\nACGT";
    const rules: validator.Rules = .{ .max_header_bytes = 20, .schema = .refseq };
    const report = try reportBytes(std.testing.allocator, data, rules);
    defer std.testing.allocator.destroy(report);
    try expectFindings(report, &.{
        .{ .code = .@"fasta.content_before_header", .occurrences = 1, .records = 0, .offset = 3, .line = 1, .column = 4, .record = null },
        .{ .code = .@"fasta.header_control_byte", .occurrences = 1, .records = 1, .offset = 23, .line = 4, .column = 3, .record = 2 },
        .{ .code = .@"fasta.empty_record", .occurrences = 2, .records = 2, .offset = 21, .line = 4, .column = 1, .record = 2 },
        .{ .code = .@"fasta.empty_identifier", .occurrences = 1, .records = 1, .offset = 34, .line = 7, .column = 2, .record = 4 },
        .{ .code = .@"fasta.duplicate_identifier", .occurrences = 1, .records = 1, .offset = 28, .line = 5, .column = 2, .record = 3 },
        .{ .code = .@"fasta.long_header", .occurrences = 1, .records = 1, .offset = 35, .line = 8, .column = 1, .record = 5 },
        .{ .code = .@"header.convention_mismatch", .occurrences = 4, .records = 4, .offset = 8, .line = 2, .column = 1, .record = 1 },
        .{ .code = .@"format.leading_bom", .occurrences = 1, .records = 0, .offset = 0, .line = 1, .column = 1, .record = null },
        .{ .code = .@"format.missing_terminal_newline", .occurrences = 1, .records = 0, .offset = 64, .line = 9, .column = 5, .record = null },
        .{ .code = .@"format.mixed_line_endings", .occurrences = 1, .records = 0, .offset = 25, .line = 4, .column = 5, .record = null },
    });
    const control = report.finding(.@"fasta.header_control_byte").first.?.details.@"fasta.header_control_byte";
    try std.testing.expectEqual(@as(u8, 1), control.byte);
    try std.testing.expectEqual(diagnostics.Grammar.controls_except_tab, control.grammar);
    const long = report.finding(.@"fasta.long_header").first.?.details.@"fasta.long_header";
    try std.testing.expectEqual(@as(u64, 23), long.observed_bytes);
    try std.testing.expectEqual(@as(u64, 20), long.limit_bytes);
    try std.testing.expectEqual(@as(u64, 4), report.finding(.@"fasta.content_before_header").first.?.details.@"fasta.content_before_header".content_bytes.exact);
    try std.testing.expectEqual(reporting.Validation.Examined{ .bytes = 64, .lines = 9, .records = 5, .records_completed = 5, .symbols = 10 }, report.examined);
    try std.testing.expectEqual(@as(u64, 3), report.records_with_errors);
    try std.testing.expectEqual(@as(u64, 5), report.records_with_warnings);
    try std.testing.expectEqual(reporting.Validation.Completion.complete, report.completion);
}

test "[unit] - [validate report]: counts sequence and format codes in their units" {
    const data = ">s1\nACGTACGT\nAC GT\tAC\n\nAC?G\x00T\n>s2\nACGTACGT\nACGTACGTAC\n  \nA*\x7f\n";
    const report = try reportBytes(std.testing.allocator, data, .{});
    defer std.testing.allocator.destroy(report);
    try expectFindings(report, &.{
        .{ .code = .@"fasta.sequence_control_byte", .occurrences = 1, .records = 1, .offset = 27, .line = 5, .column = 5, .record = 1 },
        .{ .code = .@"sequence.invalid_symbol", .occurrences = 3, .records = 2, .offset = 25, .line = 5, .column = 3, .record = 1 },
        .{ .code = .@"format.sequence_whitespace", .occurrences = 2, .records = 2, .offset = 15, .line = 3, .column = 3, .record = 1 },
        .{ .code = .@"format.blank_sequence_line", .occurrences = 1, .records = 1, .offset = 22, .line = 4, .column = 1, .record = 1 },
        .{ .code = .@"format.inconsistent_wrapping", .occurrences = 2, .records = 2, .offset = 13, .line = 3, .column = 1, .record = 1 },
    });
    try std.testing.expectEqual(@as(u64, 2), report.finding(.@"format.sequence_whitespace").first.?.details.@"format.sequence_whitespace".whitespace_bytes);
    const wrapping = report.finding(.@"format.inconsistent_wrapping").first.?.details.@"format.inconsistent_wrapping";
    try std.testing.expectEqual(@as(u64, 8), wrapping.expected_width);
    try std.testing.expectEqual(@as(u64, 6), wrapping.actual_width);
    const invalid = report.finding(.@"sequence.invalid_symbol").first.?.details.@"sequence.invalid_symbol";
    try std.testing.expectEqual(@as(u8, '?'), invalid.byte);
    try std.testing.expectEqual(diagnostics.Alphabet.nucleotide, invalid.policy.alphabet);
    try std.testing.expectEqual(@as(u32, 40), invalid.policy.sample.?.symbols);
    for (report.invalid_symbols, 0..) |occurrences, byte| {
        try std.testing.expectEqual(@as(u64, if (byte == '*' or byte == '?' or byte == 0x7f) 1 else 0), occurrences);
    }
    try std.testing.expectEqual(reporting.Validation.Examined{ .bytes = 61, .lines = 10, .records = 2, .records_completed = 2, .symbols = 40 }, report.examined);
    try std.testing.expectEqual(@as(u64, 2), report.records_with_errors);
    try std.testing.expectEqual(@as(u64, 2), report.records_with_warnings);
}

test "[edge] - [validate report]: repeated findings in one record count the record once" {
    const report = try reportBytes(std.testing.allocator, ">a\nA\x01C\x02\nA C\n AC\n>b\nACGT\nACGTAC\n", .{});
    defer std.testing.allocator.destroy(report);
    try expectFindings(report, &.{
        .{ .code = .@"fasta.sequence_control_byte", .occurrences = 2, .records = 1, .offset = 4, .line = 2, .column = 2, .record = 1 },
        .{ .code = .@"format.sequence_whitespace", .occurrences = 2, .records = 1, .offset = 9, .line = 3, .column = 2, .record = 1 },
        .{ .code = .@"format.inconsistent_wrapping", .occurrences = 1, .records = 1, .offset = 24, .line = 7, .column = 1, .record = 2 },
    });
    try std.testing.expectEqual(@as(u64, 1), report.records_with_errors);
    try std.testing.expectEqual(@as(u64, 2), report.records_with_warnings);
}

test "[edge] - [validate report]: empty and whitespace-only input report only missing records" {
    for ([_][]const u8{ "", "  \n\t\n", "\r\n" }) |data| {
        const report = try reportBytes(std.testing.allocator, data, .{});
        defer std.testing.allocator.destroy(report);
        try expectFindings(report, &.{
            .{ .code = .@"fasta.no_records", .occurrences = 1, .records = 0, .offset = null, .line = null, .column = null, .record = null },
        });
        try std.testing.expectEqual(@as(u64, data.len), report.examined.bytes);
        try std.testing.expectEqual(@as(u8, 1), report.status(false, false));
    }
}

test "[edge] - [validate report]: a duplicate points at its first definition" {
    const data = "\xEF\xBB\xBF>dup one\nAC\n>other\nAC\n>dup two\nAC\n>dup three\n";
    const report = try reportBytes(std.testing.allocator, data, .{});
    defer std.testing.allocator.destroy(report);
    const duplicate = report.finding(.@"fasta.duplicate_identifier");
    try std.testing.expectEqual(@as(u64, 2), duplicate.occurrences);
    try std.testing.expectEqual(@as(?diagnostics.Reference, null), duplicate.first.?.related);

    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    try temporary.dir.writeFile(io, .{ .sub_path = "dup.fa", .data = data });
    const file = try temporary.dir.openFile(io, "dup.fa", .{});
    defer file.close(io);
    var buffer: [8]u8 = undefined;
    var file_reader = file.reader(io, &.{});
    const from_file = try std.testing.allocator.create(validator.Report);
    defer std.testing.allocator.destroy(from_file);
    var failure: diagnostics.Failure = undefined;
    try validator.validate(std.testing.allocator, &file_reader.interface, &buffer, .{ .io = io, .file = file, .size = data.len }, .{}, from_file, &failure);
    const related = from_file.finding(.@"fasta.duplicate_identifier").first.?.related.?;
    try std.testing.expectEqual(@as(?u64, 4), related.location.offset);
    try std.testing.expectEqual(@as(?u64, 1), related.location.line);
    try std.testing.expectEqual(@as(?u64, 5), related.location.column);
    try std.testing.expectEqual(@as(?u64, 1), related.record);
}

test "[edge] - [validate report]: compares identifiers up to 65,535 bytes and reports longer ones" {
    const allocator = std.testing.allocator;
    const at_limit = "N" ** 65_535;
    const over = "N" ** 70_000;
    const data = ">" ++ at_limit ++ "\nAC\n>b\nAC\n>" ++ at_limit ++ "\nAC\n>" ++ over ++ "\nAC\n>" ++ over ++ "x\nAC\n>" ++ at_limit ++ "N\nAC\n";
    const copied = try reportBytes(allocator, data, .{});
    defer allocator.destroy(copied);

    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    try temporary.dir.writeFile(io, .{ .sub_path = "long.fa", .data = data });
    const file = try temporary.dir.openFile(io, "long.fa", .{});
    defer file.close(io);
    var buffer: [4096]u8 = undefined;
    var file_reader = file.reader(io, &.{});
    const from_file = try allocator.create(validator.Report);
    defer allocator.destroy(from_file);
    var failure: diagnostics.Failure = undefined;
    try validator.validate(allocator, &file_reader.interface, &buffer, .{ .io = io, .file = file, .size = data.len }, .{}, from_file, &failure);

    for ([_]*const validator.Report{ copied, from_file }) |report| {
        const duplicate = report.finding(.@"fasta.duplicate_identifier");
        try std.testing.expectEqual(@as(u64, 1), duplicate.occurrences);
        try std.testing.expectEqual(@as(?u64, 3), duplicate.first.?.at.record);
        const long = report.finding(.@"fasta.identifier_too_long");
        try std.testing.expectEqual(@as(u64, 3), long.occurrences);
        try std.testing.expectEqual(@as(?u64, 4), long.first.?.at.record);
        const facts = long.first.?.details.@"fasta.identifier_too_long";
        try std.testing.expectEqual(@as(u64, 70_000), facts.observed_bytes);
        try std.testing.expectEqual(@as(u64, 65_535), facts.limit_bytes);
        try std.testing.expectEqual(@as(u64, 70_000), long.first.?.identifier.?.original_bytes);
        try std.testing.expectEqual(@as(u8, 1), report.status(false, false));
    }
}

test "[edge] - [validate report]: layout does not change symbol and alphabet counts" {
    const layouts = [_][]const u8{
        ">a\nACGTNACGTN\nACG*\n>b\nMKV?\n",
        ">a\r\nACGTN\r\nACGTN\r\nACG*\r\n>b\r\nMK\r\nV?\r\n",
        ">a\nAC GT\tNACGTNACG*\n>b\nM K V ?\n",
    };
    for (layouts) |data| {
        const report = try reportBytes(std.testing.allocator, data, .{});
        defer std.testing.allocator.destroy(report);
        try std.testing.expectEqual(@as(u64, 18), report.examined.symbols);
        try std.testing.expectEqual(@as(u64, 2), report.examined.records);
        try std.testing.expectEqual(diagnostics.Alphabet.protein, report.alphabet.alphabet);
        const invalid = report.finding(.@"sequence.invalid_symbol");
        try std.testing.expectEqual(@as(u64, 1), invalid.occurrences);
        try std.testing.expectEqual(@as(u64, 1), invalid.records);
        try std.testing.expectEqual(@as(u8, '?'), invalid.first.?.details.@"sequence.invalid_symbol".byte);
        try std.testing.expectEqual(@as(?u64, 2), invalid.first.?.at.record);
    }
}

test "[edge] - [validate report]: more than 10,000 findings keep exact counts and acceptance" {
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(std.testing.allocator);
    for (0..10_001) |index| try data.print(std.testing.allocator, ">r{d}\n", .{index});
    const report = try reportBytes(std.testing.allocator, data.items, .{});
    defer std.testing.allocator.destroy(report);
    try std.testing.expectEqual(@as(u64, 10_001), report.finding(.@"fasta.empty_record").occurrences);
    try std.testing.expectEqual(reporting.Validation.Completion.complete, report.completion);
    try std.testing.expectEqual(@as(u8, 2), report.status(false, false));
    try std.testing.expectEqual(@as(u8, 1), report.status(true, false));
}

test "[edge] - [validate report]: a name budget stops only the duplicate check" {
    const allocator = std.testing.allocator;
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(allocator);
    for (0..40_000) |index| try data.print(allocator, ">name{d:0>6}\nAC\n", .{index % 30_000});

    const whole = try reportBytes(allocator, data.items, .{});
    defer allocator.destroy(whole);
    try std.testing.expectEqual(@as(u64, 10_000), whole.finding(.@"fasta.duplicate_identifier").occurrences);
    try std.testing.expectEqual(reporting.Validation.Completion.complete, whole.completion);

    var temporary = std.testing.tmpDir(.{});
    defer temporary.cleanup();
    try temporary.dir.writeFile(io, .{ .sub_path = "names.fa", .data = data.items });
    const file = try temporary.dir.openFile(io, "names.fa", .{});
    defer file.close(io);
    for ([_]bool{ false, true }) |from_file| {
        errdefer std.debug.print("from file: {}\n", .{from_file});
        const report = try allocator.create(validator.Report);
        defer allocator.destroy(report);
        var buffer: [4096]u8 = undefined;
        var fixed = std.Io.Reader.fixed(data.items);
        var file_reader = file.reader(io, &.{});
        const source: *std.Io.Reader = if (from_file) &file_reader.interface else &fixed;
        const names: ?validator.NameSource = if (from_file) .{ .io = io, .file = file, .size = data.items.len } else null;
        var failure: diagnostics.Failure = undefined;
        // A file-backed set reserves its name cache first, so it needs about 1 MiB before any name.
        const budget: u64 = if (from_file) 1024 * 1024 else 700 * 1024;
        try validator.validate(allocator, source, &buffer, names, .{ .max_name_bytes = budget }, report, &failure);
        try std.testing.expect(failure.complete);
        const stop = report.name_budget_stop.?;
        try std.testing.expectEqual(budget, stop.budget_bytes);
        try std.testing.expect(stop.stopped_at_record > 1 and stop.stopped_at_record <= 30_000);
        try std.testing.expectEqual(@as(u64, 0), report.finding(.@"fasta.duplicate_identifier").occurrences);
        try std.testing.expectEqual(reporting.Validation.Completion.incomplete, report.completion);
        try std.testing.expectEqual(@as(u64, data.items.len), report.examined.bytes);
        try std.testing.expectEqual(@as(u64, 40_000), report.examined.records_completed);
        const view = report.validation(.{}, false);
        try std.testing.expectEqual(reporting.Validation.CheckState.incomplete, view.check(.identifiers));
        try std.testing.expectEqual(reporting.Validation.CheckState.passed, view.check(.structure));
        try std.testing.expectEqual(reporting.Validation.Outcome.incomplete, view.outcome());
        try std.testing.expectEqual(@as(u8, 1), report.status(false, false));
    }
}

// --- Reading ---

test "[property] - [validate scan]: read splits give identical reports" {
    const allocator = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(std.testing.random_seed);
    var inputs: [5][]const u8 = undefined;
    inputs[0] = try allocator.dupe(u8, "\xEF\xBB\xBFjunk\n>a desc\nACGT\n>a\x01x\r\n>a\nAC\n>\n>NC_1 a long description\nACGT");
    inputs[1] = try allocator.dupe(u8, ">s1\nACGTACGT\nAC GT\tAC\n\nAC?G\x00T\r\n>s2\nACGTACGT\nACGTACGTAC\n  \nA*\x7f\r");
    inputs[2] = try allocator.dupe(u8, "\xEF\xBB");
    inputs[3] = try largeInput(allocator, prng.random(), false);
    inputs[4] = try largeInput(allocator, prng.random(), true);
    defer for (inputs) |input| allocator.free(input);
    for (inputs) |data| {
        const whole = try reportBytes(allocator, data, .{});
        defer allocator.destroy(whole);
        const expected = try renderJson(allocator, whole, .{});
        defer allocator.free(expected);
        for ([_]usize{ 1, 2, 3, 31, 64, 4096 }) |size| {
            errdefer std.debug.print("buffer {d}, input {d} bytes\n", .{ size, data.len });
            const split = try reportWith(allocator, data, size, .{});
            defer allocator.destroy(split);
            const actual = try renderJson(allocator, split, .{});
            defer allocator.free(actual);
            try std.testing.expectEqualStrings(expected, actual);
        }
    }
}

test "[property] - [validate scan]: vector checks agree with the per-byte path" {
    const allocator = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(std.testing.random_seed);
    for (0..6) |round| {
        const protein = round % 2 == 1;
        const data = try largeInput(allocator, prng.random(), protein);
        defer allocator.free(data);
        const vector = try reportBytes(allocator, data, .{});
        defer allocator.destroy(vector);
        try std.testing.expectEqual(if (protein) diagnostics.Alphabet.protein else diagnostics.Alphabet.nucleotide, vector.alphabet.alphabet);
        const per_byte = try reportBytes(allocator, data, .{ .custom_alphabet = if (protein) main.fasta.PROTEIN_SYMBOLS else main.fasta.NUCLEOTIDE_SYMBOLS });
        defer allocator.destroy(per_byte);
        try expectSameFindings(vector, per_byte);
    }
}

test "[failure] - [validate report]: a read failure keeps what was examined" {
    const data = ">a\nACGT\n>b\nACGT\n>c\nACGT\n";
    var source = FailingSource.init(data, 13);
    const report = try std.testing.allocator.create(validator.Report);
    defer std.testing.allocator.destroy(report);
    var buffer: [4]u8 = undefined;
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.validate(std.testing.allocator, &source.reader, &buffer, null, .{}, report, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.read_failed", failure.primary().?.code());
    try std.testing.expectEqual(@as(?u64, 12), failure.primary().?.at.location.offset);
    try std.testing.expectEqual(reporting.Validation.Completion.interrupted, report.completion);
    try std.testing.expectEqual(@as(u64, 12), report.examined.bytes);
    try std.testing.expectEqual(@as(u64, 3), report.examined.lines);
    try std.testing.expectEqual(@as(u64, 2), report.examined.records);
    try std.testing.expectEqual(@as(u64, 1), report.examined.records_completed);
    try std.testing.expectEqual(@as(u64, 0), report.finding(.@"fasta.no_records").occurrences);
    try std.testing.expectEqual(@as(u8, 1), report.status(false, false));
}

test "[failure] - [validate report]: running out of memory interrupts the report" {
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(std.testing.allocator);
    for (0..64) |index| try data.print(std.testing.allocator, ">r{d:0>2} {s}\nACGT\n", .{ index % 40, "description" });
    try std.testing.checkAllAllocationFailures(std.testing.allocator, checkReportAllocations, .{data.items});
}

test "[fuzz] - [validate scan]: completes on any bytes and repairs idempotently" {
    try std.testing.fuzz({}, fuzzValidator, .{});
}

// --- Rules ---

test "[property] - [validate alphabet]: accepts every symbol of the chosen alphabet" {
    const nucleotide = try std.fmt.allocPrint(std.testing.allocator, ">n\n{s}\n", .{main.fasta.NUCLEOTIDE_SYMBOLS});
    defer std.testing.allocator.free(nucleotide);
    const protein = try std.fmt.allocPrint(std.testing.allocator, ">p\n{s}\n", .{main.fasta.PROTEIN_SYMBOLS});
    defer std.testing.allocator.free(protein);
    for ([_][]const u8{ nucleotide, protein }) |data| {
        const report = try reportBytes(std.testing.allocator, data, .{});
        defer std.testing.allocator.destroy(report);
        try std.testing.expectEqual(@as(u64, 0), report.finding(.@"sequence.invalid_symbol").occurrences);
    }
}

test "[unit] - [validate alphabet]: a custom alphabet is exact and case-sensitive" {
    const report = try reportBytes(std.testing.allocator, ">s\nACGTacgt-\n", .{ .custom_alphabet = "ACGT-" });
    defer std.testing.allocator.destroy(report);
    try std.testing.expectEqual(@as(u64, 4), report.finding(.@"sequence.invalid_symbol").occurrences);
    try std.testing.expectEqual(diagnostics.Alphabet.custom, report.alphabet.alphabet);
    try std.testing.expectEqual(@as(?diagnostics.Sample, null), report.alphabet.sample);
}

test "[unit] - [validate alphabet]: explicit alphabets check every symbol" {
    const data = ">a\nACGTU\n>b\nacgun\n>c\nMKV*-\n";
    const Case = struct { alphabet: validator.Alphabet, occurrences: u64, records: u64, offset: ?u64, column: ?u64 };
    const cases = [_]Case{
        .{ .alphabet = .dna, .occurrences = 4, .records = 3, .offset = 7, .column = 5 },
        .{ .alphabet = .rna, .occurrences = 3, .records = 2, .offset = 6, .column = 4 },
        .{ .alphabet = .protein, .occurrences = 0, .records = 0, .offset = null, .column = null },
        .{ .alphabet = .any, .occurrences = 0, .records = 0, .offset = null, .column = null },
    };
    for (cases) |case| {
        errdefer std.debug.print("alphabet {s}\n", .{@tagName(case.alphabet)});
        const report = try reportBytes(std.testing.allocator, data, .{ .alphabet = case.alphabet });
        defer std.testing.allocator.destroy(report);
        const invalid = report.finding(.@"sequence.invalid_symbol");
        try std.testing.expectEqual(case.occurrences, invalid.occurrences);
        try std.testing.expectEqual(case.records, invalid.records);
        try std.testing.expectEqual(case.offset, if (invalid.first) |first| first.at.location.offset else null);
        try std.testing.expectEqual(case.column, if (invalid.first) |first| first.at.location.column else null);
        try std.testing.expectEqual(diagnostics.AlphabetSelection.explicit, report.alphabet.selection);
        try std.testing.expectEqual(@as(?diagnostics.Sample, null), report.alphabet.sample);
    }
    const dna = try reportBytes(std.testing.allocator, data, .{ .alphabet = .dna });
    defer std.testing.allocator.destroy(dna);
    for ([_]u8{ 'U', 'u', '*', '-' }) |byte| try std.testing.expectEqual(@as(u64, 1), dna.invalid_symbols[byte]);
    const any = try reportBytes(std.testing.allocator, ">a\nA?\x01 C\n", .{ .alphabet = .any });
    defer std.testing.allocator.destroy(any);
    try std.testing.expectEqual(@as(u64, 1), any.finding(.@"fasta.sequence_control_byte").occurrences);
    try std.testing.expectEqual(@as(u64, 1), any.finding(.@"format.sequence_whitespace").occurrences);
    try std.testing.expectEqual(@as(u64, 3), any.examined.symbols);
}

test "[property] - [validate alphabet]: an explicit alphabet does not depend on record order" {
    const allocator = std.testing.allocator;
    const first = try reportBytes(allocator, ">p\nMKVLAQ\n>n\nACGTACGT\n", .{ .alphabet = .dna });
    defer allocator.destroy(first);
    const second = try reportBytes(allocator, ">n\nACGTACGT\n>p\nMKVLAQ\n", .{ .alphabet = .dna });
    defer allocator.destroy(second);
    try std.testing.expectEqual(first.finding(.@"sequence.invalid_symbol").occurrences, second.finding(.@"sequence.invalid_symbol").occurrences);
    try std.testing.expectEqual(@as(u64, 2), first.finding(.@"sequence.invalid_symbol").occurrences);
    try std.testing.expectEqualSlices(u64, &first.invalid_symbols, &second.invalid_symbols);
}

test "[property] - [validate scan]: explicit vector checks agree with the per-byte path" {
    const allocator = std.testing.allocator;
    var everything: [223]u8 = undefined;
    for (&everything, 0x21..) |*byte, value| byte.* = @intCast(value);
    var prng = std.Random.DefaultPrng.init(std.testing.random_seed);
    const Case = struct { protein: bool, alphabet: validator.Alphabet, symbols: []const u8 };
    const cases = [_]Case{
        .{ .protein = false, .alphabet = .dna, .symbols = main.fasta.DNA_SYMBOLS },
        .{ .protein = false, .alphabet = .rna, .symbols = main.fasta.RNA_SYMBOLS },
        .{ .protein = true, .alphabet = .protein, .symbols = main.fasta.PROTEIN_SYMBOLS },
        .{ .protein = false, .alphabet = .any, .symbols = &everything },
    };
    for (cases) |case| {
        errdefer std.debug.print("alphabet {s}\n", .{@tagName(case.alphabet)});
        const data = try largeInput(allocator, prng.random(), case.protein);
        defer allocator.free(data);
        const vector = try reportBytes(allocator, data, .{ .alphabet = case.alphabet });
        defer allocator.destroy(vector);
        const per_byte = try reportBytes(allocator, data, .{ .custom_alphabet = case.symbols });
        defer allocator.destroy(per_byte);
        try expectSameFindings(vector, per_byte);
    }
}

test "[unit] - [validate schemas]: checks UniProt and RefSeq identifiers" {
    const data = ">sp|P12345|NAME_HUMAN desc\nM\n>tr|Q1|X\nM\n>db|A|B\nM\n>sp|P1\nM\n>NC_000001.11 chr1\nA\n>chr2\nA\n";
    const uniprot = try reportBytes(std.testing.allocator, data, .{ .schema = .uniprot });
    defer std.testing.allocator.destroy(uniprot);
    try std.testing.expectEqual(@as(u64, 3), uniprot.finding(.@"header.convention_mismatch").occurrences);
    try std.testing.expectEqual(@as(?u64, 4), uniprot.finding(.@"header.convention_mismatch").first.?.at.record);
    const refseq = try reportBytes(std.testing.allocator, data, .{ .schema = .refseq });
    defer std.testing.allocator.destroy(refseq);
    try std.testing.expectEqual(@as(u64, 5), refseq.finding(.@"header.convention_mismatch").occurrences);
}

test "[edge] - [validate scan]: warns only above the header byte limit" {
    const report = try reportBytes(std.testing.allocator, ">abcd\nA\n>abcde\nA\n>ab\r\nA\r\n", .{ .max_header_bytes = 4 });
    defer std.testing.allocator.destroy(report);
    const long = report.finding(.@"fasta.long_header");
    try std.testing.expectEqual(@as(u64, 1), long.occurrences);
    try std.testing.expectEqual(@as(?u64, 2), long.first.?.at.record);
}

// --- Identifiers in JSON ---

test "[unit] - [validate json]: identifiers keep their bytes" {
    const allocator = std.testing.allocator;
    var long_name: [300]u8 = undefined;
    @memset(&long_name, 'n');
    const data = ">bad\xff\n>bad\xc3\xbf\nAC\n>bad\xc3\xbf\n>c\\q\"\x01\nAC\n";
    const report = try reportBytes(allocator, data, .{ .max_header_bytes = 4 });
    defer allocator.destroy(report);
    const json = try renderJson(allocator, report, .{ .max_header_bytes = 4 });
    defer allocator.free(json);
    var parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    const codes = parsed.value.object.get("validation").?.object.get("codes").?.array.items;
    const Want = struct { code: []const u8, encoding: []const u8, bytes: []const u8, original: []const u8 };
    const wants = [_]Want{
        .{ .code = "fasta.header_control_byte", .encoding = "utf8", .bytes = "c\\q\"\x01", .original = "5" },
        .{ .code = "fasta.empty_record", .encoding = "base64", .bytes = "bad\xff", .original = "4" },
        .{ .code = "fasta.duplicate_identifier", .encoding = "utf8", .bytes = "bad\xc3\xbf", .original = "5" },
    };
    var buffer: [diagnostics.PREVIEW_BYTES]u8 = undefined;
    for (wants) |want| {
        const code = for (codes) |code| {
            if (std.mem.eql(u8, code.object.get("code").?.string, want.code)) break code;
        } else return error.TestUnexpectedResult;
        const identifier = code.object.get("first").?.object.get("identifier").?.object;
        try std.testing.expectEqualStrings(want.encoding, identifier.get("encoding").?.string);
        const value = identifier.get("value").?.string;
        const bytes = if (std.mem.eql(u8, want.encoding, "base64")) blk: {
            const len = try std.base64.standard.Decoder.calcSizeForSlice(value);
            try std.base64.standard.Decoder.decode(buffer[0..len], value);
            break :blk buffer[0..len];
        } else value;
        try std.testing.expectEqualSlices(u8, want.bytes, bytes);
        try std.testing.expectEqualStrings(want.original, identifier.get("original_bytes").?.string);
        try std.testing.expect(!identifier.get("truncated").?.bool);
    }

    const long_input = ">" ++ long_name ++ "\n";
    const long_report = try reportBytes(allocator, long_input, .{});
    defer allocator.destroy(long_report);
    const long_json = try renderJson(allocator, long_report, .{});
    defer allocator.free(long_json);
    var long_parsed = try std.json.parseFromSlice(std.json.Value, allocator, long_json, .{});
    defer long_parsed.deinit();
    const preview = long_parsed.value.object.get("validation").?.object.get("codes").?.array.items[0].object.get("first").?.object.get("identifier").?.object;
    try std.testing.expect(preview.get("truncated").?.bool);
    try std.testing.expectEqualStrings("300", preview.get("original_bytes").?.string);
    try std.testing.expectEqualStrings(long_name[0..diagnostics.PREVIEW_BYTES], preview.get("value").?.string);
}

test "[unit] - [validate json]: states the outcome, checks, and totals" {
    const allocator = std.testing.allocator;
    const report = try reportBytes(allocator, ">a\nAC?\n>a\n", .{});
    defer allocator.destroy(report);
    const json = try renderJson(allocator, report, .{});
    defer allocator.free(json);
    var parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    const validation = parsed.value.object.get("validation").?.object;
    try std.testing.expectEqualStrings("invalid", validation.get("outcome").?.string);
    const checks = validation.get("checks").?.object;
    try std.testing.expectEqualStrings("warnings", checks.get("structure").?.string);
    try std.testing.expectEqualStrings("failed", checks.get("identifiers").?.string);
    try std.testing.expectEqualStrings("failed", checks.get("alphabet").?.string);
    try std.testing.expectEqualStrings("passed", checks.get("format").?.string);
    const findings = validation.get("findings").?.object;
    try std.testing.expectEqualStrings("2", findings.get("errors").?.string);
    try std.testing.expectEqualStrings("1", findings.get("warnings").?.string);
    try std.testing.expectEqualStrings("2", findings.get("error_codes").?.string);
    const categories = validation.get("categories").?.object;
    try std.testing.expectEqualStrings("1", categories.get("identifiers").?.object.get("errors").?.string);
    try std.testing.expectEqualStrings("1", categories.get("structure").?.object.get("warnings").?.string);
    try std.testing.expectEqualStrings("2", validation.get("examined").?.object.get("records_completed").?.string);
}

// --- Repair ---

test "[unit] - [validate fix]: rewrites layout and keeps headers and symbols" {
    const allocator = std.testing.allocator;
    const fixed = try repairBytes(allocator, "\xEF\xBB\xBF>a desc\r\nAC GT\r\nACG\n\nTTTTTT\t\n>b\n>c\nAC?", 4);
    defer allocator.free(fixed);
    try std.testing.expectEqualStrings(">a desc\nACGT\nACGT\nTTTT\nT\n>b\n>c\nAC?\n", fixed);
    const clean = try repairBytes(allocator, ">x\nACGTACGT\nACG\n", 8);
    defer allocator.free(clean);
    try std.testing.expectEqualStrings(">x\nACGTACGT\nACG\n", clean);
}

test "[failure] - [validate fix]: refuses what it cannot keep" {
    for ([_][]const u8{ "junk\n>a\nAC\n", ">a\x01\nAC\n", ">a\nA\x00C\n", "", "  \n" }) |data| {
        try std.testing.expectError(error.RepairRefused, repairBytes(std.testing.allocator, data, 60));
    }
}

test "[property] - [validate fix]: repairs once and keeps duplicates and invalid symbols" {
    const inputs = [_][]const u8{
        ">dup\nAC??\n>dup\nAAAA  \nCCCC\nGGGGTT\n",
        ">mixed\r\nAAAACCCCGGGG\nTTTTAA\r\nAACCCCGGGGTT\nTT",
        "\xEF\xBB\xBF>e\n>f\n\nACGT\n",
    };
    for (inputs) |data| try checkValidatorProperties(data);
}

test "[integration] - [validate fix]: preserves indexed retrieval across layout rewrites" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const cases = [_]struct { stem: []const u8, broken: []const u8, region: []const u8, expected: []const u8 }{
        .{ .stem = "validate-fix-crafted", .broken = ">messy_seq widths and trailing ws\nAAAA    \nCCCC\nGGGGTT\n", .region = "messy_seq:1-14", .expected = ">messy_seq:1-14\nAAAACCCCGGGGTT\n" },
        .{ .stem = "validate-fix-mixed", .broken = ">mixed_widths internal line widths vary\nAAAACCCCGGGG\nTTTTAA\nAACCCCGGGGTT\nTT\n", .region = "mixed_widths:3-24", .expected = ">mixed_widths:3-24\nAACCCCGGGGTTTTAAAACCCC\n" },
        .{ .stem = "validate-fix-trailing", .broken = ">trailing_whitespace spaces and tabs after sequence bytes\nAAAACCCC    \nGGGGTTTT\t\nCCCCAAAA\n", .region = "trailing_whitespace:1-16", .expected = ">trailing_whitespace:1-16\nAAAACCCCGGGGTTTT\n" },
    };
    for (cases) |case| try expectFixPreservesRegion(allocator, case.stem, case.broken, case.region, case.expected);
}

test "[integration] - [validate]: agrees with index on records and symbols" {
    const allocator = std.testing.allocator;
    const data = try readTestFile(allocator, "tests/data/validator_indexer_agreement.fasta");
    defer allocator.free(data);
    const report = try reportBytes(allocator, data, .{});
    defer allocator.destroy(report);
    try std.testing.expectEqual(@as(u64, 6), report.examined.records);
    try std.testing.expectEqual(@as(u64, 1), report.finding(.@"fasta.empty_record").occurrences);
    for ([_]diagnostics.Code{ .@"format.inconsistent_wrapping", .@"format.sequence_whitespace", .@"format.mixed_line_endings", .@"format.missing_terminal_newline", .@"fasta.long_header" }) |code| {
        const found = for (report.findings) |finding| {
            if (finding.code == code) break finding.occurrences;
        } else 0;
        try std.testing.expect(found >= 1);
    }
    var index = try scanZfi(allocator, data);
    defer index.deinit(allocator);
    try std.testing.expectEqual(@as(usize, 5), index.records.items.len);
    var indexed_symbols: u64 = 0;
    for (index.records.items) |record| indexed_symbols += record.seq_len;
    try std.testing.expectEqual(indexed_symbols, report.examined.symbols);
}

// --- CLI ---

test "[cli] - [validate]: rejects unknown options regardless of position" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const fasta = "tests/data/simple.fasta";
    const stderr = "error [usage.unknown_option] validate: Unknown option. option=\"--not-a-flag\"\nhint: Correct the command arguments; use --help for accepted options.\n";
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--not-a-flag", fasta }, 1, "", stderr);
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", fasta, "--not-a-flag" }, 1, "", stderr);
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--summary", fasta }, 1, "", "error [usage.unknown_option] validate: Unknown option. option=\"--summary\"\nhint: Correct the command arguments; use --help for accepted options.\n");
}

test "[cli] - [validate]: rejects invalid arguments with exact diagnostics" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const hint = "\nhint: Correct the command arguments; use --help for accepted options.\n";
    const cases = [_]struct { argv: []const []const u8, stderr: []const u8 }{
        .{ .argv = &.{ ZFASTA_BIN, "validate", "-o", "zig-cache/test-artifacts/unused.fasta", "tests/data/simple.fasta" }, .stderr = "error [usage.option_requires] validate: An option prerequisite is missing. option=output required=fix" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--width", "70", "tests/data/simple.fasta" }, .stderr = "error [usage.option_requires] validate: An option prerequisite is missing. option=width required=fix" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "--width", "0", "tests/data/simple.fasta" }, .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=width expected=positive_usize value=\"0\"" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--width" }, .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=width expected=positive_usize" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate" }, .stderr = "error [usage.missing_argument] validate: A required argument is missing. argument=fasta_path" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "tests/data/definitely-missing-cli-failure.fasta" }, .stderr = "error [io.file_not_found] validate: \"tests/data/definitely-mi\" (47 bytes): Required file was not found.; cause=\"FileNotFound\"\nhint: Check the named input/output and retained I/O cause.\n" },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "tests/data/simple.fasta/child" }, .stderr = "error [io.open_failed] validate: \"tests/data/simple.fasta/\" (29 bytes): Could not open the required input.; cause=\"NotDir\"\nhint: Check the named input/output and retained I/O cause.\n" },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "tests/data/simple.fasta", "tests/data/single.fasta" }, .stderr = "error [usage.extra_argument] validate: Too many arguments. argument=fasta_path allowed_count=1" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "-o" }, .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=output expected=path" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--schema" }, .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=schema expected=schema" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--schema", "ena", "tests/data/simple.fasta" }, .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=schema expected=schema value=\"ena\"" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--custom-alphabet" }, .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=custom_alphabet expected=alphabet" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--max-header-len" }, .stderr = "error [usage.missing_option_value] validate: A required option value is missing. option=max_header_len expected=positive_usize" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--max-header-len", "0", "tests/data/simple.fasta" }, .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=max_header_len expected=positive_usize value=\"0\"" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "tests/data/simple.fasta" }, .stderr = "error [usage.option_requires] validate: \"tests/data/simple.fasta\": An option prerequisite is missing. option=fix required=output" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--max-name-memory", "0", "tests/data/simple.fasta" }, .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=max_name_memory expected=positive_usize value=\"0\"" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--type", "dnaa", "tests/data/simple.fasta" }, .stderr = "error [usage.invalid_option_value] validate: Invalid option value. option=sequence_type expected=sequence_type value=\"dnaa\"" ++ hint },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--type", "dna", "--custom-alphabet", "ACGT", "tests/data/simple.fasta" }, .stderr = "error [usage.conflicting_options] validate: Selected options conflict. arguments=custom_alphabet,sequence_type" ++ hint },
    };
    for (cases) |case| try expectCliResult(allocator, case.argv, 1, "", case.stderr);
}

test "[cli] - [validate]: returns exact status and text for each outcome" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const warn = try writeFastaArtifact(allocator, "cli-validate-warn", ">empty_rec\n");
    defer std.Io.Dir.cwd().deleteFile(io, warn) catch {};
    const invalid = try writeFastaArtifact(allocator, "cli-validate-error", "not FASTA\n");
    defer std.Io.Dir.cwd().deleteFile(io, invalid) catch {};

    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "tests/data/simple.fasta" }, 0,
        \\OK tests/data/simple.fasta: no findings
        \\complete: examined 82 bytes, 5 lines, 2 records, 36 symbols; nucleotide alphabet chosen from 36 symbols
        \\
    , "");
    const warning_text = try std.fmt.allocPrint(allocator,
        \\WARNINGS {s}: 0 errors in 0 codes, 1 warning in 1 code; 0 of 1 record with errors, 1 with warnings
        \\warning [fasta.empty_record] 1 record; first at line 1 column 1, record 1 "empty_rec": FASTA record contains no sequence symbols.
        \\complete: examined 11 bytes, 1 line, 1 record, 0 symbols; nucleotide alphabet chosen from 0 symbols
        \\
    , .{warn});
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", warn }, 2, warning_text, "");
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--strict", warn }, 1, warning_text, "");
    const error_text = try std.fmt.allocPrint(allocator,
        \\INVALID {s}: 2 errors in 2 codes, 0 warnings in 0 codes; 0 of 0 records with errors, 0 with warnings
        \\error [fasta.no_records]: No FASTA records were found.
        \\error [fasta.content_before_header] 1 line; first at line 1 column 1: Content appears before the first FASTA header. content_bytes=9
        \\complete: examined 10 bytes, 1 line, 0 records, 0 symbols; nucleotide alphabet chosen from 0 symbols
        \\
    , .{invalid});
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", invalid }, 1, error_text, "");
}

test "[cli] - [validate]: each option reaches its rule" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "cli-validate-options", ">ab x\nACGT\n>cd\nAC\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const cases = [_]struct { options: []const []const u8, code: []const u8, occurrences: []const u8 }{
        .{ .options = &.{ "--type", "rna" }, .code = "sequence.invalid_symbol", .occurrences = "1" },
        .{ .options = &.{ "--custom-alphabet", "ACG" }, .code = "sequence.invalid_symbol", .occurrences = "1" },
        .{ .options = &.{ "--schema", "uniprot" }, .code = "header.convention_mismatch", .occurrences = "2" },
        .{ .options = &.{ "--max-header-len", "2" }, .code = "fasta.long_header", .occurrences = "1" },
    };
    for (cases) |case| {
        const argv = try std.mem.concat(allocator, []const u8, &.{ &.{ ZFASTA_BIN, "validate", "--json" }, case.options, &.{path} });
        const result = try std.process.run(allocator, threaded.io(), .{ .argv = argv });
        const validation = (try std.json.parseFromSlice(std.json.Value, allocator, result.stdout, .{})).value.object.get("validation").?.object;
        const codes = validation.get("codes").?.array.items;
        try std.testing.expectEqual(@as(usize, 1), codes.len);
        try std.testing.expectEqualStrings(case.code, codes[0].object.get("code").?.string);
        try std.testing.expectEqualStrings(case.occurrences, codes[0].object.get("occurrences").?.string);
    }
    const budget = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", "--max-name-memory", "1", path } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, budget.term);
    const report = (try std.json.parseFromSlice(std.json.Value, allocator, budget.stdout, .{})).value.object.get("validation").?.object;
    try std.testing.expectEqualStrings("incomplete", report.get("completion").?.string);
    try std.testing.expectEqualStrings("1", report.get("name_budget").?.object.get("stopped_at_record").?.string);
}

test "[cli] - [validate json]: writes one report object with the status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "cli-validate-json", ">empty_rec\n");
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    for ([_][]const []const u8{
        &.{ ZFASTA_BIN, "validate", "--json", path },
        &.{ ZFASTA_BIN, "validate", "--json", "--error-format=json", path },
    }) |argv| {
        const result = try std.process.run(allocator, threaded.io(), .{ .argv = argv });
        switch (result.term) {
            .exited => |code| try std.testing.expectEqual(@as(u8, 2), code),
            else => return error.ChildProcessFailed,
        }
        try std.testing.expectEqualStrings("", result.stderr);
        try std.testing.expect(std.mem.endsWith(u8, result.stdout, "}\n"));
        try std.testing.expectEqual(@as(usize, 1), std.mem.count(u8, result.stdout, "\n"));
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, result.stdout, .{});
        const object = parsed.value.object;
        try std.testing.expectEqualStrings("summary", object.get("kind").?.string);
        try std.testing.expect(object.get("accepted").?.bool);
        try std.testing.expect(object.get("failure").? == .null);
        const validation = object.get("validation").?.object;
        try std.testing.expectEqualStrings("warnings", validation.get("outcome").?.string);
        try std.testing.expectEqualStrings("complete", validation.get("completion").?.string);
        const code = validation.get("codes").?.array.items[0].object;
        try std.testing.expectEqualStrings("fasta.empty_record", code.get("code").?.string);
        try std.testing.expectEqualStrings("empty_rec", code.get("first").?.object.get("identifier").?.object.get("value").?.string);
    }
}

test "[cli] - [validate fix]: writes the repaired copy and reports the input" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input = try writeFastaArtifact(allocator, "validate-fix-input", ">a\r\nACGT ACGT\r\nAC\r\n>a\nAC?\n");
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const output = try utility.uniqueArtifactPath(allocator, "validate-fix-output", "fa");
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "--width", "3", "-o", output, input } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, result.term);
    try std.testing.expect(std.mem.startsWith(u8, result.stdout, "INVALID "));
    try std.testing.expectEqualStrings("", result.stderr);
    try std.testing.expectEqualStrings(">a\nACG\nTAC\nGTA\nC\n>a\nAC?\n", try readTestFile(allocator, output));

    const warnings = try writeFastaArtifact(allocator, "validate-fix-warnings", ">b\nACGT  \nAC\n");
    defer std.Io.Dir.cwd().deleteFile(io, warnings) catch {};
    const repaired = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", "--fix", "-o", output, warnings } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 0 }, repaired.term);
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, repaired.stdout, .{});
    try std.testing.expectEqualStrings("warnings", parsed.value.object.get("validation").?.object.get("outcome").?.string);
    try std.testing.expectEqualStrings(">b\nACGTAC\n", try readTestFile(allocator, output));

    // Under --strict, the published copy's repaired warnings no longer fail the run, in the status and in `accepted`.
    const strict = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", "--strict", "--fix", "-o", output, warnings } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 0 }, strict.term);
    try std.testing.expect((try std.json.parseFromSlice(std.json.Value, allocator, strict.stdout, .{})).value.object.get("accepted").?.bool);
    const unrepaired = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", "--strict", warnings } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, unrepaired.term);
    try std.testing.expect(!(try std.json.parseFromSlice(std.json.Value, allocator, unrepaired.stdout, .{})).value.object.get("accepted").?.bool);
}

test "[cli] - [validate fix]: refuses blocked rewrites and keeps the destination" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const output = try writeFastaArtifact(allocator, "validate-blocked-output", "previous destination\n");
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    for ([_][]const u8{ "junk\n>seq1\nACGT\n", ">a\x01\nAC\n", ">a\nA\x00C\n", "" }) |data| {
        const input = try writeFastaArtifact(allocator, "validate-blocked-input", data);
        defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
        const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--error-format=json", "--fix", "-o", output, input } });
        try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, result.term);
        try std.testing.expect(std.mem.startsWith(u8, result.stdout, "INVALID "));
        try utility.expectDiagnostic(allocator, result.stderr, .@"repair.blocked");
        try std.testing.expectEqualStrings("previous destination\n", try readTestFile(allocator, output));
        try std.testing.expectEqualStrings(data, try readTestFile(allocator, input));
    }
}

test "[cli] - [validate fix]: rejects an aliased input path without modifying it" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const broken = ">seq\nAAAA  \n";
    const in_path = try writeFastaArtifact(allocator, "validate-fix-input-alias", broken);
    defer std.Io.Dir.cwd().deleteFile(io, in_path) catch {};
    const output_alias = try std.fmt.allocPrint(allocator, "./{s}", .{in_path});
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--error-format=json", "--fix", "-o", output_alias, in_path } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, result.term);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try utility.expectDiagnostic(allocator, result.stderr, .@"repair.output_is_input");
    try std.testing.expectEqualStrings(broken, try readTestFile(allocator, in_path));
}

test "[cli] - [validate fix]: rejects a symlink alias and safely replaces a separate hard link" {
    const allocator = std.testing.allocator;
    const broken = ">seq1\nACGT";
    const input = try writeFastaArtifact(allocator, "validate-link-input", broken);
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const hard = try utility.uniqueArtifactPath(allocator, "validate-hard-output", "fa");
    defer allocator.free(hard);
    defer std.Io.Dir.cwd().deleteFile(io, hard) catch {};
    const sym = try utility.uniqueArtifactPath(allocator, "validate-sym-output", "fa");
    defer allocator.free(sym);
    defer std.Io.Dir.cwd().deleteFile(io, sym) catch {};
    const absolute_input = try std.Io.Dir.cwd().realPathFileAlloc(io, input, allocator);
    defer allocator.free(absolute_input);
    try std.Io.Dir.cwd().symLink(io, absolute_input, sym, .{});
    try std.Io.Dir.cwd().hardLink(input, std.Io.Dir.cwd(), hard, io, .{});
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const rejected = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", sym, input } });
    defer allocator.free(rejected.stdout);
    defer allocator.free(rejected.stderr);
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, rejected.term);
    const accepted = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", hard, input } });
    defer allocator.free(accepted.stdout);
    defer allocator.free(accepted.stderr);
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 0 }, accepted.term);
    const input_after = try readTestFile(allocator, input);
    defer allocator.free(input_after);
    const output_after = try readTestFile(allocator, hard);
    defer allocator.free(output_after);
    try std.testing.expectEqualStrings(broken, input_after);
    try std.testing.expectEqualStrings(">seq1\nACGT\n", output_after);
}

// --- Operational failures ---

test "[failure] - [validate fix]: refuses an output open failure before staging" {
    const allocator = std.testing.allocator;
    const broken = ">seq1\nACGT";
    const input = try writeFastaArtifact(allocator, "validate-identity-input", broken);
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const output = try utility.uniqueArtifactPath(allocator, "validate-identity-loop", "fa");
    defer allocator.free(output);
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    const directory = try std.Io.Dir.cwd().realPathFileAlloc(io, ".", allocator);
    defer allocator.free(directory);
    const absolute_output = try std.fs.path.join(allocator, &.{ directory, output });
    defer allocator.free(absolute_output);
    try std.Io.Dir.cwd().symLink(io, absolute_output, output, .{});
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.RepairRefused, validator.runValidate(allocator, io, input, .{ .fix = true, .output_path = output }, &failure));
    const primary = failure.primary().?;
    try std.testing.expectEqual(diagnostics.Code.@"repair.output_identity_unknown", primary.code());
    try std.testing.expectEqual(diagnostics.Stage.check_identity, primary.stage);
    try std.testing.expectEqual(diagnostics.SourceId.output, primary.at.source.?);
    try std.testing.expectEqual(error.SymLinkLoop, primary.cause.?.zig);
    try std.testing.expectEqual(.not_started, failure.outputs.file.state);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
    try std.testing.expect(!failure.work_complete and !failure.complete);
    const after = try readTestFile(allocator, input);
    defer allocator.free(after);
    try std.testing.expectEqualStrings(broken, after);
}

test "[failure] - [validate diagnostics]: retains a stat cause before scanning" {
    const Fault = struct {
        fn stat(_: ?*anyopaque, _: std.Io.File) std.Io.File.StatError!std.Io.File.Stat {
            return error.AccessDenied;
        }
    };
    var vtable = io.vtable.*;
    vtable.fileStat = Fault.stat;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{}, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.stat_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.AccessDenied, failure.primary().?.cause.?.zig);
    try std.testing.expect(!failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.not_started, failure.outputs.stdout.state);
}

test "[failure] - [validate diagnostics]: preserves published output after report flush failure" {
    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => |write| if (write.file.handle == std.Io.File.stdout().handle)
                    .{ .file_write_streaming = error.BrokenPipe }
                else
                    io.vtable.operate(userdata, operation),
                else => io.vtable.operate(userdata, operation),
            };
        }
    };
    const allocator = std.testing.allocator;
    const input = try writeFastaArtifact(allocator, "validate-published-input", ">seq\nACGT");
    defer allocator.free(input);
    defer std.Io.Dir.cwd().deleteFile(io, input) catch {};
    const output = try writeFastaArtifact(allocator, "validate-published-output", "old destination\n");
    defer allocator.free(output);
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    var vtable = io.vtable.*;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.runValidate(allocator, fault_io, input, .{ .fix = true, .output_path = output }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(diagnostics.SourceId.stdout, failure.primary().?.at.source);
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.published, failure.outputs.file.state);
    try std.testing.expectEqual(.succeeded, failure.outputs.file.final_flush);
    try std.testing.expectEqual(.failed, failure.outputs.stdout.final_flush);
    const written = try readTestFile(allocator, output);
    defer allocator.free(written);
    try std.testing.expectEqualStrings(">seq\nACGT\n", written);
}

test "[failure] - [validate diagnostics]: retains the final flush cause and resets for reuse" {
    const Fault = struct {
        fn operate(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => .{ .file_write_streaming = error.BrokenPipe },
                else => io.vtable.operate(userdata, operation),
            };
        }
        fn discard(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            return switch (operation) {
                .file_write_streaming => |write| blk: {
                    var size = write.header.len;
                    for (write.data, 0..) |data, i| size += data.len * (if (i + 1 == write.data.len) write.splat else 1);
                    break :blk .{ .file_write_streaming = size };
                },
                else => io.vtable.operate(userdata, operation),
            };
        }
    };
    var vtable = io.vtable.*;
    vtable.operate = Fault.operate;
    const fault_io = std.Io{ .userdata = io.userdata, .vtable = &vtable };
    var failure: diagnostics.Failure = undefined;
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{ .json = true }, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.flush_failed", failure.primary().?.code());
    try std.testing.expectEqual(error.BrokenPipe, failure.primary().?.cause.?.zig);
    try std.testing.expect(failure.work_complete and !failure.complete);
    try std.testing.expectEqual(.failed, failure.outputs.stdout.final_flush);
    try std.testing.expectError(error.Io, validator.runValidate(std.testing.allocator, fault_io, "tests/data/validate-file-does-not-exist.fa", .{}, &failure));
    try std.testing.expectEqual(diagnostics.Code.@"io.file_not_found", failure.primary().?.code());
    try std.testing.expect(!failure.work_complete and !failure.complete);
    vtable.operate = Fault.discard;
    const result = try validator.runValidate(std.testing.allocator, fault_io, "tests/data/simple.fasta", .{ .json = true }, &failure);
    try std.testing.expectEqual(@as(u8, 0), result.status);
    try std.testing.expect(failure.work_complete and failure.complete);
    try std.testing.expect(failure.primary() == null);
    try std.testing.expectEqual(.succeeded, failure.outputs.stdout.final_flush);
}

test "[cli] - [validate]: reads standard input" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const child_io = threaded.io();
    const data = ">a\nAC GT\n>a\nAC\n";
    const output = try utility.uniqueArtifactPath(allocator, "validate-stdin-output", "fa");
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    const Run = struct { argv: []const []const u8, status: u8, stdout: []const u8 };
    const runs = [_]Run{
        .{ .argv = &.{ ZFASTA_BIN, "validate", "-" }, .status = 1, .stdout =
        \\INVALID stdin: 1 error in 1 code, 1 warning in 1 code; 1 of 2 records with errors, 1 with warnings
        \\error [fasta.duplicate_identifier] 1 record; first at line 3 column 2, record 2 "a": Identifier duplicates an earlier record.
        \\warning [format.sequence_whitespace] 1 line in 1 record; first at line 2 column 3, record 1 "a": Sequence line contains whitespace. whitespace_bytes=1
        \\complete: examined 15 bytes, 4 lines, 2 records, 6 symbols; nucleotide alphabet chosen from 6 symbols
        \\
        },
        .{ .argv = &.{ ZFASTA_BIN, "validate", "--fix", "-o", output, "-" }, .status = 1, .stdout = "" },
    };
    for (runs) |run| {
        var child = try std.process.spawn(child_io, .{ .argv = run.argv, .stdin = .pipe, .stdout = .pipe, .stderr = .pipe });
        try child.stdin.?.writeStreamingAll(child_io, data);
        child.stdin.?.close(child_io);
        child.stdin = null;
        var stdout_buffer: [4096]u8 = undefined;
        var stdout_reader = child.stdout.?.readerStreaming(child_io, &stdout_buffer);
        const stdout = try stdout_reader.interface.allocRemaining(allocator, .limited(64 * 1024));
        var stderr_buffer: [4096]u8 = undefined;
        var stderr_reader = child.stderr.?.readerStreaming(child_io, &stderr_buffer);
        const stderr = try stderr_reader.interface.allocRemaining(allocator, .limited(64 * 1024));
        const term = try child.wait(child_io);
        try std.testing.expectEqual(std.process.Child.Term{ .exited = run.status }, term);
        try std.testing.expectEqualStrings("", stderr);
        if (run.stdout.len != 0) try std.testing.expectEqualStrings(run.stdout, stdout) else try std.testing.expect(std.mem.startsWith(u8, stdout, "INVALID stdin: "));
    }
    try std.testing.expectEqualStrings(">a\nACGT\n>a\nAC\n", try readTestFile(allocator, output));
}

// Writes `fasta`, then a `.zfi` built from `indexed`, which has the same size, with the file's identity.
fn discardValidateStdout(userdata: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
    return switch (operation) {
        .file_write_streaming => |write| blk: {
            var size = write.header.len;
            for (write.data, 0..) |data, i| size += data.len * (if (i + 1 == write.data.len) write.splat else 1);
            break :blk .{ .file_write_streaming = size };
        },
        else => io.vtable.operate(userdata, operation),
    };
}

fn runValidateForAllocationCheck(allocator: std.mem.Allocator, path: []const u8, options: validator.Options) !void {
    var vtable = io.vtable.*;
    vtable.operate = discardValidateStdout;
    var failure: diagnostics.Failure = undefined;
    _ = validator.runValidate(allocator, .{ .userdata = io.userdata, .vtable = &vtable }, path, options, &failure) catch |err| {
        if (err == error.OutOfMemory) try std.testing.expectEqual(diagnostics.Code.@"resource.out_of_memory", failure.primary().?.code());
        return err;
    };
}

fn writeIndexedPair(allocator: std.mem.Allocator, stem: []const u8, fasta: []const u8, indexed: []const u8) ![]u8 {
    try std.testing.expectEqual(fasta.len, indexed.len);
    const path = try writeFastaArtifact(allocator, stem, fasta);
    try writeZfi(allocator, path, indexed);
    return path;
}

fn deletePair(allocator: std.mem.Allocator, path: []const u8, suffix: []const u8) void {
    const sidecar = std.fmt.allocPrint(allocator, "{s}{s}", .{ path, suffix }) catch return;
    std.Io.Dir.cwd().deleteFile(io, sidecar) catch {};
    std.Io.Dir.cwd().deleteFile(io, path) catch {};
}

test "[failure] - [validate]: every allocation failure in a run is reported and released" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const data = ">a desc\nACGT  \nAC\n>b\nGG\nGGGG\n>a\nTT\n";
    const path = try writeFastaArtifact(allocator, "validate-allocation", data);
    defer deletePair(allocator, path, ".zfi");
    try writeZfi(allocator, path, data);
    const output = try utility.uniqueArtifactPath(allocator, "validate-allocation-output", "fa");
    defer std.Io.Dir.cwd().deleteFile(io, output) catch {};
    for ([_]validator.Options{
        .{},
        .{ .fix = true, .output_path = output },
        .{ .index = .automatic },
    }) |options| {
        try std.testing.checkAllAllocationFailures(std.testing.allocator, runValidateForAllocationCheck, .{ path, options });
    }
}

test "[cli] - [validate index]: matches the records index writes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const data = ">e\n>b desc\nACGT\nAC\n>b\nAAAA\n>c\nGG\nG\nGGG\n";
    {
        const path = try writeIndexedPair(allocator, "validate-index-clean", data, data);
        defer deletePair(allocator, path, ".zfi");
        const expected = try std.fmt.allocPrint(allocator,
            \\INVALID {s}: 1 error in 1 code, 2 warnings in 2 codes; 1 of 4 records with errors, 2 with warnings
            \\warning [fasta.empty_record] 1 record; first at line 1 column 1, record 1 "e": FASTA record contains no sequence symbols.
            \\error [fasta.duplicate_identifier] 1 record; first at line 5 column 2, record 3 "b": Identifier duplicates an earlier record. First defined at line 2 column 2.
            \\warning [format.inconsistent_wrapping] 1 record; first at line 9 column 1, record 4 "c": Sequence line width differs from the record's wrap width. expected_width=2 actual_width=1
            \\complete: examined 39 bytes, 10 lines, 4 records, 16 symbols; nucleotide alphabet chosen from 16 symbols; {s}
            \\
        , .{ path, "zfi index: 3 of 3 entries matched" });
        try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index", path }, 1, expected, "");
    }

    const uniform = try writeFastaArtifact(allocator, "validate-index-fai", ">a\nACGT\nAC\n>b\nGG\n");
    defer deletePair(allocator, uniform, ".fai");
    {
        const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{uniform});
        const fai = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai.close(io);
        try std.Io.File.writeStreamingAll(fai, io, "a\t6\t3\t4\t5\nb\t2\t14\t2\t3\n");
    }
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index", uniform }, 0, try std.fmt.allocPrint(allocator,
        \\OK {s}: no findings
        \\complete: examined 17 bytes, 5 lines, 2 records, 8 symbols; nucleotide alphabet chosen from 8 symbols; fai index: 2 of 2 entries matched
        \\
    , .{uniform}), "");

    // With both sidecars, `--index` takes the `.zfi` as `get` does, and `--index=PATH` the named one.
    try writeZfi(allocator, uniform, ">a\nACGT\nAC\n>b\nGG\n");
    defer deletePair(allocator, uniform, ".zfi");
    const explicit = try std.fmt.allocPrint(allocator, "--index={s}.fai", .{uniform});
    for ([_]struct { option: []const u8, format: []const u8 }{ .{ .option = "--index", .format = "zfi" }, .{ .option = explicit, .format = "fai" } }) |case| {
        try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", case.option, uniform }, 0, try std.fmt.allocPrint(allocator,
            \\OK {s}: no findings
            \\complete: examined 17 bytes, 5 lines, 2 records, 8 symbols; nucleotide alphabet chosen from 8 symbols; {s} index: 2 of 2 entries matched
            \\
        , .{ uniform, case.format }), "");
    }

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--json", explicit, uniform } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 0 }, result.term);
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, result.stdout, .{});
    const validation = parsed.value.object.get("validation").?.object;
    try std.testing.expectEqualStrings("passed", validation.get("checks").?.object.get("index_match").?.string);
    const index = validation.get("index").?.object;
    try std.testing.expectEqualStrings("fai", index.get("format").?.string);
    try std.testing.expectEqualStrings("2", index.get("entries").?.string);
    try std.testing.expectEqualStrings("2", index.get("matched").?.string);
    try std.testing.expectEqualStrings("0", index.get("omitted_duplicates").?.string);
}

test "[cli] - [validate index]: an index that leaves out later duplicates still matches" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try writeFastaArtifact(allocator, "validate-index-omitted", ">a\nAC\n>b\nGT\n>a\nTT\n");
    defer deletePair(allocator, path, ".fai");
    {
        const fai_path = try std.fmt.allocPrint(allocator, "{s}.fai", .{path});
        const fai = try std.Io.Dir.cwd().createFile(io, fai_path, .{});
        defer fai.close(io);
        try std.Io.File.writeStreamingAll(fai, io, "a\t2\t3\t2\t3\nb\t2\t9\t2\t3\n");
    }
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index", path }, 1, try std.fmt.allocPrint(allocator,
        \\INVALID {s}: 1 error in 1 code, 0 warnings in 0 codes; 1 of 3 records with errors, 0 with warnings
        \\error [fasta.duplicate_identifier] 1 record; first at line 5 column 2, record 3 "a": Identifier duplicates an earlier record. First defined at line 1 column 2.
        \\complete: examined 18 bytes, 6 lines, 3 records, 6 symbols; nucleotide alphabet chosen from 6 symbols; fai index: 2 of 2 entries matched, 1 later duplicate not indexed
        \\
    , .{path}), "");
}

test "[cli] - [validate index]: reports each difference from the index" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const cases = [_]struct { indexed: []const u8, fasta: []const u8, findings: []const u8, summary: []const u8 }{
        .{
            .indexed = ">a\nAC\n>b\nGT\n>a\nTT\n",
            .fasta = ">b\nAC\n>a\nGT\n>a\nTT\n",
            .findings =
            \\INVALID {s}: 3 errors in 2 codes, 0 warnings in 0 codes; 3 of 3 records with errors, 0 with warnings
            \\error [fasta.duplicate_identifier] 1 record; first at line 5 column 2, record 3 "a": Identifier duplicates an earlier record. First defined at line 3 column 2.
            \\error [pair.sequence_offset_mismatch] 2 records; first at line 1 column 1, record 1 "b": Indexed sequence offset differs from the FASTA. indexed_offset=9 observed_offset=3
            \\
            ,
            .summary = "examined 18 bytes, 6 lines, 3 records, 6 symbols; nucleotide alphabet chosen from 6 symbols; zfi index: 3 of 3 entries matched",
        },
        .{
            .indexed = ">a\nAC\n>b\nGT\n",
            .fasta = ">a\nAC\n>c\nGT\n",
            .findings =
            \\INVALID {s}: 2 errors in 2 codes, 0 warnings in 0 codes; 1 of 2 records with errors, 0 with warnings
            \\error [pair.record_missing_from_index] 1 record; first at line 3 column 1, record 2 "c": FASTA record has no index entry.
            \\error [pair.record_missing_from_fasta] 1 record; first at zfi entry 2 "b": Index entry has no FASTA record.
            \\
            ,
            .summary = "examined 12 bytes, 4 lines, 2 records, 4 symbols; nucleotide alphabet chosen from 4 symbols; zfi index: 1 of 2 entries matched",
        },
        .{
            .indexed = ">a\nACGT\n>b\nG\n",
            .fasta = ">a\nACG\n>b\nGT\n",
            .findings =
            \\INVALID {s}: 5 errors in 3 codes, 0 warnings in 0 codes; 2 of 2 records with errors, 0 with warnings
            \\error [pair.sequence_length_mismatch] 2 records; first at line 1 column 1, record 1 "a": Indexed sequence length differs from the FASTA. indexed_length=4 observed_length=3
            \\error [pair.sequence_offset_mismatch] 1 record; first at line 3 column 1, record 2 "b": Indexed sequence offset differs from the FASTA. indexed_offset=11 observed_offset=10
            \\error [pair.line_geometry_mismatch] 2 records; first at line 1 column 1, record 1 "a": Indexed line layout differs from the FASTA. indexed_line_bases=4 indexed_line_bytes=5 indexed_layout=uniform observed_line_bases=3 observed_line_bytes=4 observed_layout=uniform
            \\
            ,
            .summary = "examined 13 bytes, 4 lines, 2 records, 5 symbols; nucleotide alphabet chosen from 5 symbols; zfi index: 2 of 2 entries matched",
        },
        .{
            .indexed = ">a\nACGT\nAC\n",
            .fasta = ">a\nAC\nACGT\n",
            .findings =
            \\INVALID {s}: 1 error in 1 code, 1 warning in 1 code; 1 of 1 record with errors, 1 with warnings
            \\warning [format.inconsistent_wrapping] 1 record; first at line 3 column 1, record 1 "a": Sequence line width differs from the record's wrap width. expected_width=2 actual_width=4
            \\error [pair.line_geometry_mismatch] 1 record; first at line 1 column 1, record 1 "a": Indexed line layout differs from the FASTA. indexed_line_bases=4 indexed_line_bytes=5 indexed_layout=uniform observed_line_bases=2 observed_line_bytes=3 observed_layout=side_table
            \\
            ,
            .summary = "examined 11 bytes, 3 lines, 1 record, 6 symbols; nucleotide alphabet chosen from 6 symbols; zfi index: 1 of 1 entry matched",
        },
        .{
            .indexed = ">a\nAC\nACGT\nA\n",
            .fasta = ">a\nAC\nACG\nAA\n",
            .findings =
            \\INVALID {s}: 1 error in 1 code, 1 warning in 1 code; 1 of 1 record with errors, 1 with warnings
            \\warning [format.inconsistent_wrapping] 1 record; first at line 3 column 1, record 1 "a": Sequence line width differs from the record's wrap width. expected_width=2 actual_width=3
            \\error [pair.side_table_mismatch] 1 record; first at line 1 column 1, record 1 "a": Indexed side table differs from the FASTA lines. indexed_lines=3 observed_lines=3 line=2
            \\
            ,
            .summary = "examined 13 bytes, 4 lines, 1 record, 7 symbols; nucleotide alphabet chosen from 7 symbols; zfi index: 1 of 1 entry matched",
        },
    };
    for (cases) |case| {
        const path = try writeIndexedPair(allocator, "validate-index-mismatch", case.fasta, case.indexed);
        defer deletePair(allocator, path, ".zfi");
        const expected = try std.fmt.allocPrint(allocator, "{s}complete: {s}\n", .{ try std.mem.replaceOwned(u8, allocator, case.findings, "{s}", path), case.summary });
        try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index", path }, 1, expected, "");
        try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", try std.fmt.allocPrint(allocator, "--index={s}.zfi", .{path}), path }, 1, expected, "");
    }
}

test "[cli] - [validate index]: needs a file and a loadable index" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index", "-" }, 1, "", "error [usage.option_requires] validate: An option prerequisite is missing. option=index required=fasta_file\nhint: Correct the command arguments; use --help for accepted options.\n");

    const unindexed = try writeFastaArtifact(allocator, "validate-index-none", ">a\nAC\n");
    defer std.Io.Dir.cwd().deleteFile(io, unindexed) catch {};
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "validate", "--index", unindexed }, .@"index.not_found");

    const stale = try writeFastaArtifact(allocator, "validate-index-stale", ">a\nACGT\n");
    defer deletePair(allocator, stale, ".zfi");
    try writeZfi(allocator, stale, ">a\nAC\n");
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "validate", "--index", stale }, .@"index.source_size_mismatch");

    // A damaged entry stops the scan where the index stream reaches it.
    const damaged = try writeFastaArtifact(allocator, "validate-index-damaged", ">a\nACGT\n>b\nGG\n");
    defer deletePair(allocator, damaged, ".zfi");
    try writeZfi(allocator, damaged, ">a\nACGT\n>b\nGG\n");
    {
        const zfi_path = try std.fmt.allocPrint(allocator, "{s}.zfi", .{damaged});
        const file = try std.Io.Dir.cwd().openFile(io, zfi_path, .{ .mode = .read_write });
        defer file.close(io);
        try file.writePositionalAll(io, &std.mem.toBytes(@as(u32, 0)), @sizeOf(main.index_format.ZfiHeader) + main.index_format.ZFI_INDEX_RECORD_BYTES + 19);
    }
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const result = try std.process.run(allocator, threaded.io(), .{ .argv = &.{ ZFASTA_BIN, "validate", "--index", "--error-format=json", damaged } });
    try std.testing.expectEqual(std.process.Child.Term{ .exited = 1 }, result.term);
    try std.testing.expectEqualStrings(try std.fmt.allocPrint(allocator,
        \\INTERRUPTED {s}: no findings
        \\interrupted: examined 14 bytes, 4 lines, 2 records (1 completed), 6 symbols; nucleotide alphabet chosen from 6 symbols; zfi index: 1 of 1 entry matched
        \\
    , .{damaged}), result.stdout);
    try utility.expectDiagnostic(allocator, result.stderr, .@"index.invalid_line_geometry");
    try utility.expectCliFailureCode(allocator, &.{ ZFASTA_BIN, "validate", "--index=zig-cache/test-artifacts/definitely-missing.zfi", stale }, .@"io.file_not_found");

    const hint = "\nhint: Correct the command arguments; use --help for accepted options.\n";
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index=sample.idx", stale }, 1, "", "error [usage.invalid_option_value] validate: Invalid option value. option=index expected=index_path value=\"sample.idx\"" ++ hint);
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index=", stale }, 1, "", "error [usage.invalid_option_value] validate: Invalid option value. option=index expected=index_path value=\"\"" ++ hint);
    try expectCliResult(allocator, &.{ ZFASTA_BIN, "validate", "--index", "--index=sample.fai", stale }, 1, "", "error [usage.repeated_option] validate: This option cannot be repeated. option=index" ++ hint);
}
