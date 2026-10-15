//! Parses BED rows for `get --bed`: zero-based half-open in, one-based inclusive out.
//! Comment, `track`, and `browser` rows are skipped. The chrom slice borrows the input line.

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

pub const ParseResult = union(enum) {
    skip,
    region: Region,
};

pub fn parseBedLine(
    line: []const u8,
    honor_strand: bool,
    at: diagnostics.Reference,
    failure: *diagnostics.Failure,
) diagnostics.Error!ParseResult {
    const diag = BedDiagnostics{ .failure = failure, .at = at };
    const parsed = try parseLine(line, diag);
    if (honor_strand and parsed == .region and parsed.region.strand == .invalid) {
        var fields = std.mem.splitScalar(u8, requestLine(line).?, '\t');
        for (0..5) |_| _ = fields.next();
        const field = fields.next().?;
        return diag.fail(.{ .@"request.bed.invalid_strand" = .{ .strand = diag.preview(field) } }, @intFromPtr(field.ptr) - @intFromPtr(line.ptr));
    }
    return parsed;
}

/// Reads only the chrom column; the rest of the row is not checked.
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
    failure: *diagnostics.Failure,
    at: diagnostics.Reference = .{},

    fn fail(self: @This(), details: diagnostics.Details, field_start: usize) diagnostics.Error {
        @branchHint(.cold);
        var at = self.at;
        if (at.location.offset) |offset| {
            at.location.offset = std.math.add(u64, offset, field_start) catch null;
        }
        at.location.column = std.math.add(u64, at.location.column orelse 1, field_start) catch null;
        return self.failure.fail(.{ .details = details, .stage = .parse, .at = at });
    }

    fn preview(self: @This(), bytes: []const u8) ?diagnostics.Preview {
        return self.failure.retain(bytes, 32);
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

fn checkedCoordinate(field: []const u8, kind: diagnostics.BedCoordinate, field_start: usize, diag: BedDiagnostics) diagnostics.Error!u64 {
    return switch (parseCoordinate(field)) {
        .value => |value| value,
        .invalid => diag.fail(.{ .@"request.bed.invalid_coordinate" = .{
            .field = kind,
            .coordinate_system = .zero_based_half_open,
            .value = diag.preview(field),
        } }, field_start),
        .out_of_range => diag.fail(.{ .@"request.bed.coordinate_out_of_range" = .{
            .field = kind,
            .max_value = std.math.maxInt(u64),
            .coordinate_system = .zero_based_half_open,
            .value = diag.preview(field),
        } }, field_start),
    };
}

fn parseLine(line: []const u8, diag: BedDiagnostics) diagnostics.Error!ParseResult {
    const trimmed = requestLine(line) orelse return .skip;
    var fields = std.mem.splitScalar(u8, trimmed, '\t');

    const chrom = fields.next() orelse "";
    if (chrom.len == 0) return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .chrom,
        .column_number = 1,
    } }, 0);

    const start_field = fields.next() orelse return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .start,
        .column_number = 2,
    } }, trimmed.len);
    const start_pos = @intFromPtr(start_field.ptr) - @intFromPtr(trimmed.ptr);
    if (start_field.len == 0) return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .start,
        .column_number = 2,
    } }, start_pos);

    const end_field = fields.next() orelse return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .end,
        .column_number = 3,
    } }, trimmed.len);
    const end_pos = @intFromPtr(end_field.ptr) - @intFromPtr(trimmed.ptr);
    if (end_field.len == 0) return diag.fail(.{ .@"request.bed.field_missing" = .{
        .field = .end,
        .column_number = 3,
    } }, end_pos);

    const start_0based = try checkedCoordinate(start_field, .start, start_pos, diag);
    const end_0based = try checkedCoordinate(end_field, .end, end_pos, diag);
    if (end_0based <= start_0based) return diag.fail(.{ .@"request.bed.invalid_interval" = .{
        .start = start_0based,
        .end = end_0based,
        .coordinate_system = .zero_based_half_open,
    } }, start_pos);

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
