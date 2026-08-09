//! IUPAC byte complementation for GET transforms.

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

const IUPAC_COMPLEMENT_TABLE: [256]u8 = buildComplementTable();

/// Returns the IUPAC complement for one byte.
/// Unknown bytes pass through unchanged.
pub fn complement(byte: u8) u8 {
    return IUPAC_COMPLEMENT_TABLE[byte];
}

/// Write the IUPAC complement of `src` into `dst`; asserts equal lengths.
pub fn complementInto(dst: []u8, src: []const u8) void {
    std.debug.assert(dst.len == src.len);
    var i: usize = 0;
    while (src.len - i >= 32) : (i += 32) {
        inline for (0..32) |j| {
            dst[i + j] = IUPAC_COMPLEMENT_TABLE[src[i + j]];
        }
    }
    while (i < src.len) : (i += 1) {
        dst[i] = IUPAC_COMPLEMENT_TABLE[src[i]];
    }
}

test "[unit] - [complement]: maps IUPAC symbols and preserves other bytes" {
    const pairs = [_][2]u8{
        .{ 'A', 'T' },
        .{ 'C', 'G' },
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
    try std.testing.expectEqual(@as(u8, 'A'), complement('U'));
    try std.testing.expectEqual(@as(u8, 'a'), complement('u'));

    for ([_]u8{ 0, '-', 'X', 'x', 0xff }) |byte| {
        try std.testing.expectEqual(byte, complement(byte));
    }
}

test "[property] - [complementInto]: matches scalar mapping across chunk boundaries" {
    const alphabet = "ACGTRYSWKMBDHVNacgtryswkmbdhvnUu-Xx";
    var src: [97]u8 = undefined;
    for (&src, 0..) |*byte, i| byte.* = alphabet[i % alphabet.len];

    for (0..src.len + 1) |len| {
        var dst: [src.len]u8 = @splat(0xff);
        complementInto(dst[0..len], src[0..len]);

        for (src[0..len], dst[0..len]) |input, actual| {
            try std.testing.expectEqual(complement(input), actual);
        }
        if (len < dst.len) try std.testing.expectEqual(@as(u8, 0xff), dst[len]);
    }
}
