const std = @import("std");

fn setPair(table: *[256]u8, left: u8, right: u8) void {
    table[left] = right;
    table[right] = left;

    table[std.ascii.toLower(left)] = std.ascii.toLower(right);
    table[std.ascii.toLower(right)] = std.ascii.toLower(left);
}

fn buildComplementTable() [256]u8 {
    var table: [256]u8 = undefined;
    for (&table, 0..) |*slot, i| {
        slot.* = @intCast(i);
    }

    setPair(&table, 'A', 'T');
    setPair(&table, 'C', 'G');
    setPair(&table, 'R', 'Y');
    setPair(&table, 'W', 'W');
    setPair(&table, 'S', 'S');
    setPair(&table, 'K', 'M');
    setPair(&table, 'B', 'V');
    setPair(&table, 'D', 'H');
    setPair(&table, 'N', 'N');

    table['U'] = 'A';
    table['u'] = 'a';

    return table;
}

pub const iupac_complement_table: [256]u8 = buildComplementTable();

/// Return the IUPAC complement for one base.
/// Unknown bytes pass through unchanged so headers and non-sequence bytes are not mangled.
pub fn complement(byte: u8) u8 {
    return iupac_complement_table[byte];
}

/// Write the IUPAC complement of `src` into `dst` (same length).
pub fn complementInto(dst: []u8, src: []const u8) void {
    std.debug.assert(dst.len == src.len);
    var i: usize = 0;
    // Process 32-byte chunks with a table lookup per lane (no branch on base).
    while (i + 32 <= src.len) : (i += 32) {
        inline for (0..32) |j| {
            dst[i + j] = iupac_complement_table[src[i + j]];
        }
    }
    while (i < src.len) : (i += 1) {
        dst[i] = iupac_complement_table[src[i]];
    }
}

/// Write the reverse complement of `src` into `dst`.
/// Caller must provide a same-length destination buffer.
pub fn reverseComplementInto(dst: []u8, src: []const u8) void {
    std.debug.assert(dst.len == src.len);

    for (src, 0..) |byte, i| {
        dst[src.len - 1 - i] = complement(byte);
    }
}

/// Read exactly `len` bytes, reverse-complement them, and write the result.
/// Allocation-backed helper for generic readers that cannot seek backward.
///
/// Not used by the shipped `get` command. Production RC/reverse paths walk the
/// mmapped FASTA backward in `getter.zig` (zero extra copy of the region span).
/// See `bench/get/RC_STRATEGY.md`. Kept exported for tests and future streaming callers.
pub fn reverseComplementStream(allocator: std.mem.Allocator, reader: anytype, writer: anytype, len: usize) !void {
    const src = try allocator.alloc(u8, len);
    defer allocator.free(src);

    const dst = try allocator.alloc(u8, len);
    defer allocator.free(dst);

    try reader.readSliceAll(src);
    reverseComplementInto(dst, src);
    try writer.writeAll(dst);
}

test "complement handles standard bases" {
    try std.testing.expectEqual(@as(u8, 'T'), complement('A'));
    try std.testing.expectEqual(@as(u8, 'A'), complement('T'));
    try std.testing.expectEqual(@as(u8, 'G'), complement('C'));
    try std.testing.expectEqual(@as(u8, 'C'), complement('G'));
}

test "complement handles lowercase bases" {
    try std.testing.expectEqual(@as(u8, 't'), complement('a'));
    try std.testing.expectEqual(@as(u8, 'a'), complement('t'));
    try std.testing.expectEqual(@as(u8, 'g'), complement('c'));
    try std.testing.expectEqual(@as(u8, 'c'), complement('g'));
}

test "complement handles IUPAC ambiguity codes" {
    const pairs = [_][2]u8{
        .{ 'R', 'Y' },
        .{ 'W', 'W' },
        .{ 'S', 'S' },
        .{ 'K', 'M' },
        .{ 'B', 'V' },
        .{ 'D', 'H' },
        .{ 'N', 'N' },
    };

    for (pairs) |pair| {
        try std.testing.expectEqual(pair[1], complement(pair[0]));
        try std.testing.expectEqual(pair[0], complement(pair[1]));
        try std.testing.expectEqual(std.ascii.toLower(pair[1]), complement(std.ascii.toLower(pair[0])));
        try std.testing.expectEqual(std.ascii.toLower(pair[0]), complement(std.ascii.toLower(pair[1])));
    }
}

test "complement maps uracil to adenine" {
    try std.testing.expectEqual(@as(u8, 'A'), complement('U'));
    try std.testing.expectEqual(@as(u8, 'a'), complement('u'));
}

test "complement leaves unknown bytes unchanged" {
    try std.testing.expectEqual(@as(u8, '-'), complement('-'));
    try std.testing.expectEqual(@as(u8, 'X'), complement('X'));
    try std.testing.expectEqual(@as(u8, 'x'), complement('x'));
}

test "reverseComplementInto handles mixed case and ambiguity" {
    var dst: [8]u8 = undefined;
    reverseComplementInto(&dst, "AaCGnRyU");
    try std.testing.expectEqualStrings("ArYnCGtT", &dst);
}

test "reverseComplementStream reads exact length and writes reversed complement" {
    var reader = std.Io.Reader.fixed("ACGTN");
    var out_buf: [5]u8 = undefined;
    var writer = std.Io.Writer.fixed(&out_buf);

    try reverseComplementStream(std.testing.allocator, &reader, &writer, 5);

    try std.testing.expectEqualStrings("NACGT", writer.buffered());
}
