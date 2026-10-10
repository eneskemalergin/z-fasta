//! BED line parser: 0-based half-open intervals to 1-based inclusive requests.

const std = @import("std");
const diagnostics = @import("diagnostics");

pub const Strand = enum {
    forward,
    reverse,
    invalid,
};

pub const Region = struct {
    chrom: []const u8,
    start_1based: u64,
    end_1based: u64,
    strand: Strand,
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
    region: Region,
};

/// Parses a BED row and records precise field failures in caller-owned storage.
/// `at` identifies the original row; offsets and columns are byte locations.
/// The caller initializes `failure` once for the operation. This helper does
/// not reset it. Returned chromosome bytes borrow `line`; failure previews own
/// bounded copies. Invalid strands are rejected only when `honor_strand` is set.
pub fn parseBedLineChecked(
    line: []const u8,
    honor_strand: bool,
    at: diagnostics.Reference,
    failure: *diagnostics.Failure,
) diagnostics.Error!ParseResult {
    const diag = BedDiagnostics{ .failure = failure, .at = at };
    const parsed = parseLine(line, diag) catch {
        return diagnostics.metadata(failure.primary().?.code()).returned_error;
    };
    if (honor_strand and parsed == .region and parsed.region.strand == .invalid) {
        var fields = std.mem.splitScalar(u8, requestLine(line).?, '\t');
        for (0..5) |_| _ = fields.next();
        const field = fields.next().?;
        diag.capture(.{ .@"request.bed.invalid_strand" = .{ .strand = diag.preview(field) } }, @intFromPtr(field.ptr) - @intFromPtr(line.ptr));
        return diagnostics.metadata(failure.primary().?.code()).returned_error;
    }
    return parsed;
}

/// Parses one BED row into a 1-based inclusive request.
///
/// Empty, comment, track, and browser lines return `.skip`. Coordinates must
/// contain only decimal digits and satisfy `end > start`. `chrom` borrows
/// from `line`. Invalid strands are retained for the caller to accept or reject.
pub fn parseBedLine(line: []const u8) ParseError!ParseResult {
    return parseLine(line, .{});
}

/// Borrows a chromosome for lookup collection. Output callers still parse the complete row.
pub fn requestChrom(line: []const u8) ?[]const u8 {
    const trimmed = requestLine(line) orelse return null;
    const end = std.mem.findScalar(u8, trimmed, '\t') orelse trimmed.len;
    if (end == 0) return null;
    return trimmed[0..end];
}

fn parseStrand(field: []const u8) Strand {
    if (field.len == 0) return .forward;
    if (field.len != 1) return .invalid;

    return switch (field[0]) {
        '+', '.' => .forward,
        '-' => .reverse,
        else => .invalid,
    };
}

const BedDiagnostics = struct {
    failure: ?*diagnostics.Failure = null,
    at: diagnostics.Reference = .{},

    fn fail(self: @This(), details: diagnostics.Details, field_start: usize, legacy: ParseError) ParseError {
        self.capture(details, field_start);
        return legacy;
    }

    fn capture(self: @This(), details: diagnostics.Details, field_start: usize) void {
        if (self.failure) |failure| {
            var at = self.at;
            if (at.location.offset) |offset| {
                at.location.offset = std.math.add(u64, offset, field_start) catch null;
            }
            at.location.column = std.math.add(u64, at.location.column orelse 1, field_start) catch null;
            switch (failure.fail(.{ .details = details, .stage = .parse, .at = at })) {
                else => {},
            }
        }
    }

    fn preview(self: @This(), bytes: []const u8) ?diagnostics.Preview {
        return if (self.failure) |failure| failure.retain(bytes, 32) else null;
    }
};

const Coordinate = union(enum) { value: u64, invalid, out_of_range };

fn parseCoordinate(field: []const u8) Coordinate {
    if (field.len == 0) return .invalid;
    var value: u64 = 0;
    for (field, 0..) |byte, i| {
        if (byte < '0' or byte > '9') return .invalid;
        value = std.math.mul(u64, value, 10) catch return rejectedCoordinate(field[i..]);
        value = std.math.add(u64, value, byte - '0') catch return rejectedCoordinate(field[i..]);
    }
    return .{ .value = value };
}

fn rejectedCoordinate(remaining: []const u8) Coordinate {
    for (remaining) |byte| if (byte < '0' or byte > '9') return .invalid;
    return .out_of_range;
}

fn checkedCoordinate(field: []const u8, kind: diagnostics.BedCoordinate, field_start: usize, diag: BedDiagnostics) ParseError!u64 {
    const legacy: ParseError = if (kind == .start) error.InvalidStart else error.InvalidEnd;
    return switch (parseCoordinate(field)) {
        .value => |value| value,
        .invalid => diag.fail(.{ .@"request.bed.invalid_coordinate" = .{
            .field = kind,
            .coordinate_system = .zero_based_half_open,
            .value = diag.preview(field),
        } }, field_start, legacy),
        .out_of_range => diag.fail(.{ .@"request.bed.coordinate_out_of_range" = .{
            .field = kind,
            .max_value = std.math.maxInt(u64),
            .coordinate_system = .zero_based_half_open,
            .value = diag.preview(field),
        } }, field_start, legacy),
    };
}

fn parseLine(line: []const u8, diag: BedDiagnostics) ParseError!ParseResult {
    const trimmed = requestLine(line) orelse return .skip;
    var fields = std.mem.splitScalar(u8, trimmed, '\t');

    const chrom = fields.next() orelse "";
    if (chrom.len == 0) return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .chrom,
        .column_number = 1,
    } }, 0, error.MissingChrom);

    const start_field = fields.next() orelse return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .start,
        .column_number = 2,
    } }, trimmed.len, error.MissingStart);
    const start_pos = @intFromPtr(start_field.ptr) - @intFromPtr(trimmed.ptr);
    if (start_field.len == 0) return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .start,
        .column_number = 2,
    } }, start_pos, error.MissingStart);

    const end_field = fields.next() orelse return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .end,
        .column_number = 3,
    } }, trimmed.len, error.MissingEnd);
    const end_pos = @intFromPtr(end_field.ptr) - @intFromPtr(trimmed.ptr);
    if (end_field.len == 0) return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .end,
        .column_number = 3,
    } }, end_pos, error.MissingEnd);

    const start_0based = try checkedCoordinate(start_field, .start, start_pos, diag);
    const end_0based = try checkedCoordinate(end_field, .end, end_pos, diag);
    if (end_0based <= start_0based) return diag.fail(.{ .@"request.bed.invalid_interval" = .{
        .start = start_0based,
        .end = end_0based,
        .coordinate_system = .zero_based_half_open,
    } }, start_pos, error.EmptyInterval);

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

fn requestLine(line: []const u8) ?[]const u8 {
    const trimmed = if (line.len > 0 and line[line.len - 1] == '\r') line[0 .. line.len - 1] else line;
    if (trimmed.len == 0 or trimmed[0] == '#') return null;
    if (std.mem.eql(u8, trimmed, "track") or std.mem.startsWith(u8, trimmed, "track ")) return null;
    if (std.mem.eql(u8, trimmed, "browser") or std.mem.startsWith(u8, trimmed, "browser ")) return null;
    return trimmed;
}
