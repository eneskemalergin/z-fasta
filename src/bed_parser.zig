//! BED line parser: 0-based half-open intervals to 1-based inclusive requests.

const std = @import("std");

pub const BedStrand = enum {
    forward,
    reverse,
    invalid,
};

pub const BedRegion = struct {
    chrom: []const u8,
    start_1based: u64,
    end_1based: u64,
    strand: BedStrand,
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
    if (field.len == 0) return .forward;
    if (field.len != 1) return .invalid;

    return switch (field[0]) {
        '+', '.' => .forward,
        '-' => .reverse,
        else => .invalid,
    };
}

fn parseCoordinate(field: []const u8) ?u64 {
    if (std.mem.findScalar(u8, field, '_') != null) return null;
    return std.fmt.parseUnsigned(u64, field, 10) catch null;
}

/// Parses one BED row into a 1-based inclusive request.
///
/// Empty, comment, track, and browser lines return `.skip`. Coordinates must
/// contain only decimal digits and satisfy `end > start`. `chrom` borrows
/// from `line`.
pub fn parseBedLine(line: []const u8) ParseError!ParseResult {
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

    const start_0based = parseCoordinate(start_field) orelse return error.InvalidStart;
    const end_0based = parseCoordinate(end_field) orelse return error.InvalidEnd;

    if (end_0based <= start_0based) return error.EmptyInterval;

    _ = fields.next();
    _ = fields.next();

    const strand_field = fields.next();
    const strand = if (strand_field) |field| parseStrand(field) else .forward;

    return .{ .region = .{
        .chrom = chrom,
        .start_1based = start_0based + 1,
        .end_1based = end_0based,
        .strand = strand,
    } };
}
