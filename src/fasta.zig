//! What FASTA bytes mean: bases, identifiers, the nucleotide or protein alphabet, and the
//! complement. A base is any byte above ASCII space; bytes at or below it count as whitespace,
//! which `index`, `get`, and `stats` skip in sequence lines. An identifier ends at the first
//! space, tab, CR, or LF. The complement keeps case, maps U to A one way, and passes every other
//! byte through.

const std = @import("std");

pub const SequenceType = enum { nucleotide, protein };

pub fn isBase(byte: u8) bool {
    return byte > ' ';
}

pub fn isNameDelimiter(byte: u8) bool {
    return switch (byte) {
        ' ', '\t', '\r', '\n' => true,
        else => false,
    };
}

pub fn countBases(data: []const u8) u64 {
    var count: u64 = 0;
    var pos: usize = 0;
    while (pos + SIMD_CHUNK_SIZE <= data.len) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        count += std.simd.countTrues(chunk > space_char);
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < data.len) {
        if (isBase(data[pos])) count += 1;
        pos += 1;
    }
    return count;
}

pub fn countBasesInLineSlice(data: []const u8, start: usize, end: usize) u32 {
    if (start >= end) return 0;
    if (!lineSliceHasWhitespace(data, start, end)) {
        return @intCast(end - start);
    }
    return @intCast(countBases(data[start..end]));
}

pub fn countWhitespace(comptime vector_len: usize, data: []const u8) usize {
    if (data.len < 2048) return countWhitespaceDirect(vector_len, data);
    const Vec = @Vector(vector_len, u8);
    const Wide = @Vector(vector_len, u16);
    const group = 128;
    var count: usize = 0;
    var pos: usize = 0;
    while (data.len - pos >= vector_len * group) {
        var lanes: Vec = @splat(0);
        for (0..group) |_| {
            const bytes: Vec = data[pos..][0..vector_len].*;
            lanes += @select(u8, bytes <= @as(Vec, @splat(' ')), @as(Vec, @splat(1)), @as(Vec, @splat(0)));
            pos += vector_len;
        }
        count += @reduce(.Add, @as(Wide, lanes));
    }
    var lanes: Vec = @splat(0);
    while (data.len - pos >= vector_len) : (pos += vector_len) {
        const bytes: Vec = data[pos..][0..vector_len].*;
        lanes += @select(u8, bytes <= @as(Vec, @splat(' ')), @as(Vec, @splat(1)), @as(Vec, @splat(0)));
    }
    count += @reduce(.Add, @as(Wide, lanes));
    while (pos < data.len) : (pos += 1) {
        count += @intFromBool(!isBase(data[pos]));
    }
    return count;
}

pub fn findNameDelimiter(data: []const u8, start: usize) usize {
    const width = 16;
    const Vec = @Vector(width, u8);
    var position = start;
    while (data.len - position >= width) : (position += width) {
        const bytes: Vec = data[position..][0..width].*;
        const separators = (bytes == @as(Vec, @splat(' '))) |
            (bytes == @as(Vec, @splat('\t'))) |
            (bytes == @as(Vec, @splat('\r'))) |
            (bytes == @as(Vec, @splat('\n')));
        if (std.simd.firstTrue(separators)) |offset| return position + offset;
    }
    while (position < data.len) : (position += 1) {
        if (isNameDelimiter(data[position])) return position;
    }
    return position;
}

/// Nucleotide when IUPAC nucleotide letters are more than 90% of `total`; empty is nucleotide.
pub fn detectType(counts: *const [256]u64, total: u64) SequenceType {
    if (total == 0) return .nucleotide;

    const nuc_count = counts['A'] + counts['a'] +
        counts['C'] + counts['c'] +
        counts['G'] + counts['g'] +
        counts['T'] + counts['t'] +
        counts['R'] + counts['r'] +
        counts['Y'] + counts['y'] +
        counts['S'] + counts['s'] +
        counts['W'] + counts['w'] +
        counts['K'] + counts['k'] +
        counts['M'] + counts['m'] +
        counts['B'] + counts['b'] +
        counts['D'] + counts['d'] +
        counts['H'] + counts['h'] +
        counts['V'] + counts['v'] +
        counts['N'] + counts['n'] +
        counts['U'] + counts['u'];

    if (@as(u128, nuc_count) * 10 > @as(u128, total) * 9) {
        return .nucleotide;
    }
    return .protein;
}

pub fn complement(byte: u8) u8 {
    return IUPAC_COMPLEMENT_TABLE[byte];
}

/// Caller assumes `dst.len == src.len`; only Debug checks it.
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

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);

fn lineSliceHasWhitespace(data: []const u8, start: usize, end: usize) bool {
    var pos = start;
    while (pos + SIMD_CHUNK_SIZE <= end) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        const space_char: SimdVec = @splat(' ');
        const mask = chunk <= space_char;
        if (@reduce(.Or, mask)) return true;
        pos += SIMD_CHUNK_SIZE;
    }
    while (pos < end) : (pos += 1) {
        if (!isBase(data[pos])) return true;
    }
    return false;
}

fn countWhitespaceDirect(comptime vector_len: usize, data: []const u8) usize {
    const Vec = @Vector(vector_len, u8);
    var count: usize = 0;
    var pos: usize = 0;
    while (pos + vector_len <= data.len) : (pos += vector_len) {
        const chunk: Vec = data[pos..][0..vector_len].*;
        count += std.simd.countTrues(chunk <= @as(Vec, @splat(' ')));
    }
    while (pos < data.len) : (pos += 1) {
        count += @intFromBool(!isBase(data[pos]));
    }
    return count;
}

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

test "[property] - [portable scanner]: whitespace counts match scalar counts" {
    var data: [12289]u8 = undefined;
    const lengths = [_]usize{ 0, 15, 16, 31, 32, 127, 128, 2047, 2048, 2049, 4095, 4096, 4097, 8191, 8192, 8193, data.len };
    for ([_]?u8{ null, ' ', 'X' }) |fill| {
        for (&data, 0..) |*byte, i| byte.* = fill orelse @as(u8, @truncate(i * 17 + i / 31));
        for (lengths) |len| {
            var expected: usize = 0;
            for (data[0..len]) |byte| expected += @intFromBool(byte <= ' ');
            inline for (.{ 1, 8, 16, 32, 64 }) |width| {
                try std.testing.expectEqual(expected, countWhitespace(width, data[0..len]));
            }
        }
    }
}

test "[property] - [portable scanner]: name delimiters match scalar scanning" {
    var data = [_]u8{'A'} ** 65;
    for ([_]u8{ ' ', '\t', '\r', '\n' }) |separator| {
        for ([_]usize{ 0, 1, 15, 16, 17, 31, 32, 63, 64 }) |at| {
            data[at] = separator;
            for ([_]usize{ 0, 1, 15, 16, 32, data.len }) |start| {
                try std.testing.expectEqual(if (start <= at) at else data.len, findNameDelimiter(&data, start));
            }
            data[at] = 'A';
        }
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
