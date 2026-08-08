//! BED line parser: 0-based half-open intervals to 1-based inclusive regions plus strand.
//!
//! Used by GET `--bed` batching (`BedRegion`, strand-aware extract).

const std = @import("std");

pub const BedStrand = enum {
    plus,
    minus,
    none,
    invalid,
};

pub const BedRegion = struct {
    chrom: []const u8,
    start_0based: u64,
    end_0based: u64,
    strand: BedStrand,
    line_number: usize,

    pub fn start1Based(self: BedRegion) u64 {
        return self.start_0based + 1;
    }

    pub fn end1BasedInclusive(self: BedRegion) u64 {
        return self.end_0based;
    }
};

pub const ParseError = error{
    MissingChrom,
    MissingStart,
    MissingEnd,
    InvalidStart,
    InvalidEnd,
    EmptyInterval,
};

pub const ParseResult = union(enum) {
    skip,
    region: BedRegion,
};

fn parseStrand(field: []const u8) BedStrand {
    if (field.len == 0) return .none;
    if (field.len != 1) return .invalid;

    return switch (field[0]) {
        '+' => .plus,
        '-' => .minus,
        '.' => .none,
        else => .invalid,
    };
}

/// Parses one BED line. Returned region slices borrow from `line`.
pub fn parseBedLine(line: []const u8, line_number: usize) ParseError!ParseResult {
    const trimmed = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;

    if (trimmed.len == 0) return .skip;
    if (trimmed[0] == '#') return .skip;
    if (std.mem.eql(u8, trimmed, "track") or std.mem.startsWith(u8, trimmed, "track ")) return .skip;
    if (std.mem.eql(u8, trimmed, "browser") or std.mem.startsWith(u8, trimmed, "browser ")) return .skip;

    var fields = std.mem.splitScalar(u8, trimmed, '\t');

    const chrom = fields.next() orelse return error.MissingChrom;
    if (chrom.len == 0) return error.MissingChrom;

    const start_field = fields.next() orelse return error.MissingStart;
    if (start_field.len == 0) return error.MissingStart;

    const end_field = fields.next() orelse return error.MissingEnd;
    if (end_field.len == 0) return error.MissingEnd;

    const start_0based = std.fmt.parseInt(u64, start_field, 10) catch return error.InvalidStart;
    const end_0based = std.fmt.parseInt(u64, end_field, 10) catch return error.InvalidEnd;

    if (end_0based <= start_0based) return error.EmptyInterval;

    _ = fields.next(); // column 4 name, ignored
    _ = fields.next(); // column 5 score, ignored

    const strand_field = fields.next();
    const strand = if (strand_field) |field| parseStrand(field) else .none;

    return .{ .region = .{
        .chrom = chrom,
        .start_0based = start_0based,
        .end_0based = end_0based,
        .strand = strand,
        .line_number = line_number,
    } };
}

test "parseBedLine parses basic three-column BED" {
    const parsed = try parseBedLine("chr1\t100\t200", 4);
    switch (parsed) {
        .skip => return error.UnexpectedResult,
        .region => |region| {
            try std.testing.expectEqualStrings("chr1", region.chrom);
            try std.testing.expectEqual(@as(u64, 100), region.start_0based);
            try std.testing.expectEqual(@as(u64, 200), region.end_0based);
            try std.testing.expectEqual(BedStrand.none, region.strand);
            try std.testing.expectEqual(@as(u64, 101), region.start1Based());
            try std.testing.expectEqual(@as(u64, 200), region.end1BasedInclusive());
            try std.testing.expectEqual(@as(usize, 4), region.line_number);
        },
    }
}

test "parseBedLine parses sixth-column strand" {
    const plus = try parseBedLine("chr1\t0\t10\tname\t0\t+", 1);
    const minus = try parseBedLine("chr1\t0\t10\tname\t0\t-", 2);
    const dot = try parseBedLine("chr1\t0\t10\tname\t0\t.", 3);
    const invalid = try parseBedLine("chr1\t0\t10\tname\t0\t?", 4);

    try std.testing.expectEqual(BedStrand.plus, plus.region.strand);
    try std.testing.expectEqual(BedStrand.minus, minus.region.strand);
    try std.testing.expectEqual(BedStrand.none, dot.region.strand);
    try std.testing.expectEqual(BedStrand.invalid, invalid.region.strand);
}

test "parseBedLine skips comments and empty lines" {
    try std.testing.expectEqual(ParseResult.skip, try parseBedLine("", 1));
    try std.testing.expectEqual(ParseResult.skip, try parseBedLine("# comment", 2));
    try std.testing.expectEqual(ParseResult.skip, try parseBedLine("track name=foo", 3));
    try std.testing.expectEqual(ParseResult.skip, try parseBedLine("browser position chr1:1-10", 4));
}

test "parseBedLine preserves chromosome names that resemble directives" {
    const cases = [_][]const u8{
        "track1\t0\t10",
        "browser1\t0\t10",
        "track\t0\t10",
        "browser\t0\t10",
    };

    for (cases) |line| {
        const parsed = try parseBedLine(line, 1);

        try std.testing.expectEqualStrings(line[0..std.mem.indexOfScalar(u8, line, '\t').?], parsed.region.chrom);
    }
}

test "parseBedLine trims CRLF line endings" {
    const parsed = try parseBedLine("chr2\t5\t9\r", 7);
    try std.testing.expectEqualStrings("chr2", parsed.region.chrom);
    try std.testing.expectEqual(@as(u64, 6), parsed.region.start1Based());
    try std.testing.expectEqual(@as(u64, 9), parsed.region.end1BasedInclusive());
}

test "parseBedLine rejects malformed coordinates" {
    try std.testing.expectError(error.InvalidStart, parseBedLine("chr1\t-1\t10", 1));
    try std.testing.expectError(error.InvalidEnd, parseBedLine("chr1\t0\tfoo", 1));
    try std.testing.expectError(error.EmptyInterval, parseBedLine("chr1\t10\t10", 1));
    try std.testing.expectError(error.EmptyInterval, parseBedLine("chr1\t10\t9", 1));
}

test "parseBedLine rejects missing required columns" {
    try std.testing.expectError(error.MissingStart, parseBedLine("chr1", 1));
    try std.testing.expectError(error.MissingEnd, parseBedLine("chr1\t0", 1));
    try std.testing.expectError(error.MissingChrom, parseBedLine("\t0\t10", 1));
}
