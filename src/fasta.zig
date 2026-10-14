//! What FASTA bytes mean: bases, identifiers, the nucleotide or protein alphabet, and the
//! complement. A base is any byte above ASCII space; bytes at or below it count as whitespace,
//! which `index`, `get`, and `stats` skip in sequence lines. An identifier ends at the first
//! space, tab, CR, or LF. The complement keeps case, maps U to A one way, and passes every other
//! byte through.
//! `RecordReader` streams whole records under the same rules into caller storage. As in `index`,
//! the first byte must be `>`, and a record without bases is skipped.
//! Identifier sets compare names by exact bytes: `SourceDedupWith` keeps a fingerprint and an offset
//! per name and re-reads the source on a fingerprint match; `PackedNamesWith` copies names.
//! `RecordScannerWith` is `index`'s scanner: one `ScannedRecord` per record, with its line layout.

const std = @import("std");

const diagnostics = @import("diagnostics");

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
    return data.len - countWhitespace(SIMD_CHUNK_SIZE, data);
}

pub fn countBasesInLineSlice(data: []const u8, start: usize, end: usize) u32 {
    if (start >= end) return 0;
    if (!lineSliceHasWhitespace(data, start, end)) {
        return @intCast(end - start);
    }
    return @intCast(countBases(data[start..end]));
}

// Aligned so the hot loop's position in a cache line does not depend on where the linker puts it.
pub fn countWhitespace(comptime vector_len: usize, data: []const u8) align(64) usize {
    return countVectorized(vector_len, data, isWhitespace);
}

fn isWhitespace(comptime Vec: type, bytes: Vec) @Vector(@typeInfo(Vec).vector.len, bool) {
    return bytes <= @as(Vec, @splat(' '));
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

/// The rule `index` uses: a symbol is any byte above space.
pub const Bases = struct {
    pub const countOutside = countWhitespace;

    pub fn anyOutside(comptime vector_len: usize, data: []const u8) bool {
        return anyVectorized(vector_len, data, isWhitespace);
    }
};

/// IUPAC nucleotide symbols with both T and U, in either case.
pub const NUCLEOTIDE_SYMBOLS = "ACGTUNRYWSMKHBVDacgtunrywsmkhbvd";
/// IUPAC DNA symbols, without U, in either case.
pub const DNA_SYMBOLS = "ACGTNRYWSMKHBVDacgtnrywsmkhbvd";
/// IUPAC RNA symbols, U instead of T, in either case.
pub const RNA_SYMBOLS = "ACGUNRYWSMKHBVDacgunrywsmkhbvd";
/// The 26 letters, `*` for a stop, and `-` for a gap, in either case.
pub const PROTEIN_SYMBOLS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ*-abcdefghijklmnopqrstuvwxyz";

/// The common symbols of a nucleotide alphabet, which a vector check accepts: `letters` (lowercase)
/// in either case. Other IUPAC symbols count as outside, so callers check them through the full
/// alphabet.
pub fn PlainLetters(comptime letters: []const u8) type {
    return struct {
        pub fn countOutside(comptime vector_len: usize, data: []const u8) usize {
            return countOutsideWith(vector_len, data, isPlain);
        }

        pub fn anyOutside(comptime vector_len: usize, data: []const u8) bool {
            return anyOutsideWith(vector_len, data, isPlain);
        }

        fn isPlain(comptime Vec: type, bytes: Vec) @Vector(@typeInfo(Vec).vector.len, bool) {
            const folded = bytes | @as(Vec, @splat(0x20));
            var inside: @Vector(@typeInfo(Vec).vector.len, bool) = @splat(false);
            inline for (letters) |letter| inside = inside | (folded == @as(Vec, @splat(letter)));
            return inside;
        }
    };
}

pub const PlainNucleotides = PlainLetters("acgtun");
pub const PlainDna = PlainLetters("acgtn");
pub const PlainRna = PlainLetters("acgun");

/// Exactly `PROTEIN_SYMBOLS`.
pub const ProteinSymbols = struct {
    pub fn countOutside(comptime vector_len: usize, data: []const u8) usize {
        return countOutsideWith(vector_len, data, isProteinSymbol);
    }

    pub fn anyOutside(comptime vector_len: usize, data: []const u8) bool {
        return anyOutsideWith(vector_len, data, isProteinSymbol);
    }

    fn isProteinSymbol(comptime Vec: type, bytes: Vec) @Vector(@typeInfo(Vec).vector.len, bool) {
        const folded = bytes | @as(Vec, @splat(0x20));
        const letter = (folded >= @as(Vec, @splat('a'))) & (folded <= @as(Vec, @splat('z')));
        return letter | (bytes == @as(Vec, @splat('*'))) | (bytes == @as(Vec, @splat('-')));
    }
};

/// True when a line of `line_bases` symbols plus an LF or CRLF fills exactly `line_bytes`.
pub fn isDenseLineGeometry(line_bases: u32, line_bytes: u32) bool {
    if (line_bases == 0 or line_bytes <= line_bases) return false;
    const sep_len = line_bytes - line_bases;
    return (sep_len == 1 or sep_len == 2) and line_bases + sep_len == line_bytes;
}

/// Counts whole lines from `start` that repeat the geometry: `line_bases` bytes that `Symbols`
/// accepts, then the same LF or CRLF. Stops at a header, a short line, or any other byte.
pub fn countDenseLines(comptime Symbols: type, data: []const u8, start: usize, line_bases: u32, line_bytes: u32) usize {
    return countDenseLinesWith(DENSE_VECTOR_LEN, Symbols, data, start, line_bases, line_bytes);
}

/// Nucleotide when IUPAC nucleotide letters are more than 90% of `total`; empty is nucleotide.
pub fn detectType(counts: *const [256]u64, total: u64) SequenceType {
    var nucleotides: u64 = 0;
    for (IUPAC_NUCLEOTIDE_LETTERS) |letter| nucleotides += counts[letter] + counts[std.ascii.toLower(letter)];
    return typeFromCounts(nucleotides, total);
}

/// `detectType` from the number of IUPAC nucleotide letters among `total` symbols.
pub fn typeFromCounts(nucleotides: u64, total: u64) SequenceType {
    if (total == 0) return .nucleotide;
    return if (@as(u128, nucleotides) * 10 > @as(u128, total) * 9) .nucleotide else .protein;
}

pub fn isIupacNucleotide(byte: u8) bool {
    return IUPAC_NUCLEOTIDE_TABLE[byte];
}

const IUPAC_NUCLEOTIDE_LETTERS = "ACGTRYSWKMBDHVNU";
const IUPAC_NUCLEOTIDE_TABLE: [256]bool = blk: {
    var table = [_]bool{false} ** 256;
    for (IUPAC_NUCLEOTIDE_LETTERS) |letter| {
        table[letter] = true;
        table[std.ascii.toLower(letter)] = true;
    }
    break :blk table;
};

pub fn complement(byte: u8) u8 {
    return IUPAC_COMPLEMENT_TABLE[byte];
}

const Block = @Vector(32, u8);

// `complement` on 32 bytes: letters fold to upper case, swap with their partner, and get their
// case back; every other byte is left as it is.
fn complementBlock(block: Block) Block {
    const upper = block & @as(Block, @splat(0xDF));
    // Blocks of A, C, G, T, and N in either case, nearly every block of a genome, swap by XOR.
    const at = (upper == @as(Block, @splat('A'))) | (upper == @as(Block, @splat('T')));
    const cg = (upper == @as(Block, @splat('C'))) | (upper == @as(Block, @splat('G')));
    const plain = at | cg | (upper == @as(Block, @splat('N'))) | (block < @as(Block, @splat('A')));
    if (@reduce(.And, plain)) {
        const flip = @select(u8, at, @as(Block, @splat('A' ^ 'T')), @select(u8, cg, @as(Block, @splat('C' ^ 'G')), @as(Block, @splat(0))));
        return block ^ flip;
    }
    var swapped = upper;
    inline for (.{ "AT", "TA", "CG", "GC", "RY", "YR", "KM", "MK", "BV", "VB", "DH", "HD", "UA" }) |pair| {
        swapped = @select(u8, upper == @as(Block, @splat(pair[0])), @as(Block, @splat(pair[1])), swapped);
    }
    const letter = (upper >= @as(Block, @splat('A'))) & (upper <= @as(Block, @splat('Z')));
    return @select(u8, letter, swapped | (block & @as(Block, @splat(0x20))), block);
}

fn reverseBlock(block: Block) Block {
    const order = comptime blk: {
        var indices: [32]i32 = undefined;
        for (&indices, 0..) |*index, i| index.* = 31 - @as(i32, @intCast(i));
        break :blk indices;
    };
    return @shuffle(u8, block, undefined, order);
}

/// Caller assumes `dst.len == src.len`; only Debug checks it.
pub fn complementInto(dst: []u8, src: []const u8) void {
    std.debug.assert(dst.len == src.len);
    var i: usize = 0;
    while (src.len - i >= 32) : (i += 32) {
        dst[i..][0..32].* = complementBlock(src[i..][0..32].*);
    }
    while (i < src.len) : (i += 1) {
        dst[i] = IUPAC_COMPLEMENT_TABLE[src[i]];
    }
}

/// `src` written backwards into `dst`, complemented when `do_complement`; they must not overlap.
pub fn reverseInto(dst: []u8, src: []const u8, do_complement: bool) void {
    std.debug.assert(dst.len == src.len);
    var i: usize = 0;
    while (src.len - i >= 32) : (i += 32) {
        const block = reverseBlock(src[src.len - i - 32 ..][0..32].*);
        dst[i..][0..32].* = if (do_complement) complementBlock(block) else block;
    }
    while (i < src.len) : (i += 1) {
        const byte = src[src.len - 1 - i];
        dst[i] = if (do_complement) IUPAC_COMPLEMENT_TABLE[byte] else byte;
    }
}

pub const RecordReader = struct {
    source: *std.Io.Reader,
    storage: []u8,
    failure: *diagnostics.Failure,
    state: enum { start, header, done } = .start,
    failed: ?diagnostics.Error = null,
    offset: u64 = 0,
    header_offset: u64 = 0,
    headers: u64 = 0,
    byte: [1]u8 = undefined,
    byte_held: bool = false,

    pub const Record = struct {
        id: []const u8,
        description: []const u8,
        sequence: []const u8,
    };

    /// `storage` holds one record: the bytes of its header line before LF, then its bases.
    pub fn init(source: *std.Io.Reader, storage: []u8, failure: *diagnostics.Failure) RecordReader {
        failure.reset(.index);
        return .{ .source = source, .storage = storage, .failure = failure };
    }

    /// The record borrows `storage` until the next call; after a failure, each call returns that error.
    pub fn next(self: *RecordReader) diagnostics.Error!?Record {
        if (self.failed) |err| return err;
        return self.readRecord() catch |err| {
            self.failed = err;
            return err;
        };
    }

    fn readRecord(self: *RecordReader) diagnostics.Error!?Record {
        while (true) {
            switch (self.state) {
                .done => return null,
                .header => {},
                .start => {
                    const bytes = try self.window();
                    if (bytes.len == 0) {
                        self.finish();
                        return null;
                    }
                    if (bytes[0] != '>') return self.failure.fail(.{
                        .details = .{ .@"fasta.content_before_header" = .{ .content_bytes = .{ .at_least = 1 } } },
                        .stage = .scan,
                        .at = .{ .source = .fasta, .location = .{ .offset = 0, .line = 1, .column = 1 } },
                    });
                    self.startHeader();
                },
            }
            const header_len = try self.readHeader();
            const end = try self.readSequence(header_len);
            if (end > header_len) return makeRecord(self.storage[0..header_len], self.storage[header_len..end]);
        }
    }

    fn readHeader(self: *RecordReader) diagnostics.Error!usize {
        var len: usize = 0;
        while (true) {
            const bytes = try self.window();
            if (bytes.len == 0) return len;
            const line_end = std.mem.indexOfScalar(u8, bytes, '\n');
            const piece = bytes[0 .. line_end orelse bytes.len];
            len = try self.append(len, piece);
            self.consume(piece.len + @intFromBool(line_end != null));
            if (line_end != null) return len;
        }
    }

    fn readSequence(self: *RecordReader, start: usize) diagnostics.Error!usize {
        var len = start;
        var line_start = true;
        while (true) {
            const bytes = try self.window();
            if (bytes.len == 0) {
                self.finish();
                return len;
            }
            const header = findHeader(bytes, line_start);
            if (header == 0) {
                self.startHeader();
                return len;
            }
            const segment = bytes[0 .. header orelse bytes.len];
            len = try self.appendBases(len, segment);
            line_start = segment[segment.len - 1] == '\n';
            self.consume(segment.len);
        }
    }

    fn findHeader(bytes: []const u8, line_start: bool) ?usize {
        if (line_start and bytes[0] == '>') return 0;
        var from: usize = 0;
        while (std.mem.indexOfScalarPos(u8, bytes, from, '>')) |at| : (from = at + 1) {
            if (at > 0 and bytes[at - 1] == '\n') return at;
        }
        return null;
    }

    fn appendBases(self: *RecordReader, start: usize, segment: []const u8) diagnostics.Error!usize {
        var out = start;
        var in: usize = 0;
        while (segment.len - in >= 2 * SIMD_CHUNK_SIZE and self.storage.len - out >= 2 * SIMD_CHUNK_SIZE) : (in += SIMD_CHUNK_SIZE) {
            const chunk: SimdVec = segment[in..][0..SIMD_CHUNK_SIZE].*;
            const lanes: @Vector(SIMD_CHUNK_SIZE, u1) = @bitCast(chunk > @as(SimdVec, @splat(' ')));
            var bases: u64 = @as(u32, @bitCast(lanes));
            while (bases != 0) {
                const run_start = @ctz(bases);
                const run_len = @ctz(~(bases >> @intCast(run_start)));
                // The loop bounds keep this full-width store in storage; bytes past the run are rewritten.
                self.storage[out..][0..SIMD_CHUNK_SIZE].* = segment[in + run_start ..][0..SIMD_CHUNK_SIZE].*;
                out += run_len;
                bases &= ~((@as(u64, 1) << @intCast(run_start + run_len)) - 1);
            }
        }
        for (segment[in..], in..) |byte, at| {
            if (!isBase(byte)) continue;
            if (out == self.storage.len) return self.tooLarge(@as(u64, out) + countBases(segment[at..]));
            self.storage[out] = byte;
            out += 1;
        }
        return out;
    }

    fn append(self: *RecordReader, len: usize, bytes: []const u8) diagnostics.Error!usize {
        if (bytes.len > self.storage.len - len) return self.tooLarge(@as(u64, len) + bytes.len);
        @memcpy(self.storage[len..][0..bytes.len], bytes);
        return len + bytes.len;
    }

    fn window(self: *RecordReader) diagnostics.Error![]const u8 {
        if (self.source.buffer.len != 0) {
            return self.source.peekGreedy(1) catch |err| switch (err) {
                error.EndOfStream => return &.{},
                error.ReadFailed => return self.readFailed(err),
            };
        }
        if (!self.byte_held) {
            const count = self.source.readSliceShort(&self.byte) catch |err| return self.readFailed(err);
            if (count == 0) return &.{};
            self.byte_held = true;
        }
        return &self.byte;
    }

    fn consume(self: *RecordReader, count: usize) void {
        if (self.source.buffer.len != 0) {
            self.source.toss(count);
        } else {
            std.debug.assert(count == 1);
            self.byte_held = false;
        }
        self.offset +|= count;
    }

    fn startHeader(self: *RecordReader) void {
        self.header_offset = self.offset;
        self.headers +|= 1;
        self.consume(1);
        self.state = .header;
    }

    fn finish(self: *RecordReader) void {
        self.state = .done;
        self.failure.succeed();
    }

    fn reference(self: *const RecordReader, offset: u64) diagnostics.Reference {
        return .{ .source = .fasta, .location = .{ .offset = offset }, .record = if (self.headers == 0) null else self.headers };
    }

    fn tooLarge(self: *RecordReader, required: u64) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{
            .details = .{ .@"resource.line_buffer_exceeded" = .{ .buffer_bytes = self.storage.len, .observed_bytes = .{ .at_least = required } } },
            .stage = .read,
            .at = self.reference(self.header_offset),
        });
    }

    fn readFailed(self: *RecordReader, err: anyerror) diagnostics.Error {
        @branchHint(.cold);
        return self.failure.fail(.{
            .details = diagnostics.readFailure(1, null),
            .stage = .read,
            .at = self.reference(self.offset),
            .cause = .{ .zig = err },
        });
    }

    fn makeRecord(header: []const u8, sequence: []const u8) Record {
        const id_end = findNameDelimiter(header, 0);
        return .{
            .id = header[0..id_end],
            .description = std.mem.trim(u8, header[id_end..], " \t\r"),
            .sequence = sequence,
        };
    }
};

// --- Record layout ---

/// The layout of one line of a record that does not follow the uniform formula. The same 32 bytes
/// as `index_format.SideTableLine`, which `index` checks at compile time.
pub const LineEntry = extern struct {
    base_start: u64,
    byte_offset: u64,
    line_bytes: u32,
    line_bases: u32,
};

/// A record with bases, as `RecordScannerWith` found it. Slices borrow the scanner until the next
/// record.
pub const ScannedRecord = struct {
    name: []const u8,
    header_offset: u64,
    reference: *const diagnostics.Reference,
    seq_offset: u64,
    seq_len: u64,
    line_bases: u32,
    line_bytes: u32,
    uses_uniform_formula: bool,
    /// A `u64` line count, then one `LineEntry` per line; empty under the uniform formula.
    side_table: []const u8,
    nonuniform_reason: diagnostics.NonuniformReason,
};

/// `index`'s record scanner: windows in, one `ScannedRecord` per record with bases out. Line
/// geometry and side-table lines use `u32` widths. `Diag` provides `at`, `oom`, `lineLimit`,
/// `consistency`, and `identifierTooLong`.
pub fn RecordScannerWith(comptime Diag: type) type {
    return struct {
        const Self = @This();

        diag: Diag = .{},
        allocator: std.mem.Allocator,
        collect_side_tables: bool = false,
        active: bool = false,
        in_header: bool = false,
        parsing_name: bool = false,
        at_line_start: bool = true,
        header_start_offset: u64 = 0,
        source_record: u64 = 0,
        source_reference: diagnostics.Reference = .{ .source = .fasta },
        name: std.ArrayList(u8) = .empty,
        seq_offset: u64 = 0,
        seq_offset_set: bool = false,
        line_builder: LineMetricsBuilder = .{},
        in_pending_line: bool = false,
        pending_line_bases: u32 = 0,
        pending_line_bytes: u32 = 0,
        pending_pre_lf_bytes: u32 = 0,
        pending_pre_lf_ends_with_cr: bool = false,
        pending_line_file_offset: u64 = 0,
        side_table: std.ArrayList(u8) = .empty,
        side_base_start: u64 = 0,
        side_table_active: bool = false,
        deferred_side_line: ?struct {
            bases: u32,
            actual_bytes: u32,
            line_file_offset: u64,
        } = null,
        deferred_short_tail: ?struct {
            bases: u32,
            actual_bytes: u32,
            line_file_offset: u64,
        } = null,

        fn recordDiagnostics(self: *const Self) Diag {
            var diag = self.diag;
            diag.at.location.offset = self.header_start_offset;
            diag.at.record = if (self.source_record == 0) null else self.source_record;
            return diag;
        }

        fn ensureSideTableHeader(self: *Self) !void {
            if (self.side_table.items.len >= @sizeOf(u64)) return;
            var line_count: u64 = 0;
            self.side_table.appendSlice(self.allocator, std.mem.asBytes(&line_count)) catch return self.recordDiagnostics().oom(.side_tables, @sizeOf(u64));
        }

        fn activateSideTable(self: *Self, lines_to_materialize: u64) !void {
            if (self.side_table_active) return;
            self.side_table_active = true;
            try self.ensureSideTableHeader();
            if (lines_to_materialize == 0) {
                self.deferred_side_line = null;
                return;
            }
            const first = self.deferred_side_line orelse return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
            var i: u64 = 0;
            while (i < lines_to_materialize) : (i += 1) {
                const off = first.line_file_offset + i * @as(u64, first.actual_bytes);
                try self.appendSideTableLine(first.bases, first.actual_bytes, off);
            }
            self.deferred_side_line = null;
        }

        fn startHeader(self: *Self, header_start_offset: u64) !void {
            self.active = true;
            self.in_header = true;
            self.parsing_name = true;
            self.at_line_start = false;
            self.header_start_offset = header_start_offset;
            self.source_record += 1;
            self.source_reference.location.offset = header_start_offset;
            self.source_reference.record = self.source_record;
            self.name.clearRetainingCapacity();
            self.seq_offset = 0;
            self.seq_offset_set = false;
            self.line_builder = .{};
            self.in_pending_line = false;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_pre_lf_bytes = 0;
            self.pending_pre_lf_ends_with_cr = false;
            self.pending_line_file_offset = 0;
            self.side_base_start = 0;
            self.side_table.clearRetainingCapacity();
            self.side_table_active = false;
            self.deferred_side_line = null;
            self.deferred_short_tail = null;
        }

        fn appendSideTableLine(
            self: *Self,
            bases: u32,
            actual_bytes: u32,
            line_file_offset: u64,
        ) !void {
            const entry = LineEntry{
                .base_start = self.side_base_start,
                .byte_offset = line_file_offset,
                .line_bytes = actual_bytes,
                .line_bases = bases,
            };
            self.side_table.appendSlice(self.allocator, std.mem.asBytes(&entry)) catch return self.recordDiagnostics().oom(.side_tables, @sizeOf(LineEntry));
            self.side_base_start += bases;
        }

        fn flushDeferredShortTailThen(
            self: *Self,
            metrics_line_count: u64,
            bases: u32,
            actual_bytes: u32,
            line_file_offset: u64,
        ) !void {
            const short = self.deferred_short_tail orelse return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
            self.deferred_short_tail = null;
            if (metrics_line_count < 2) return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
            try self.activateSideTable(metrics_line_count - 2);
            try self.appendSideTableLine(short.bases, short.actual_bytes, short.line_file_offset);
            try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
        }

        fn recordSequenceLine(
            self: *Self,
            bases: u32,
            actual_bytes: u32,
            has_lf: bool,
            content_len: u32,
            line_file_offset: u64,
        ) !void {
            self.line_builder.ingestLine(bases, actual_bytes, has_lf, content_len);

            if (!self.collect_side_tables) return;

            if (bases == 0) return;

            const metrics = self.line_builder.metrics;
            const matches_stride = metrics.line_bases > 0 and
                bases == metrics.line_bases and
                actual_bytes == self.line_builder.first_actual_bytes;
            const dense = metrics.first_content_dense and isDenseLineGeometry(metrics.line_bases, metrics.line_bytes);
            const can_use_formula = metrics.is_uniform_width and matches_stride and dense;

            if (can_use_formula) {
                if (metrics.line_count == 1) {
                    self.deferred_side_line = .{
                        .bases = bases,
                        .actual_bytes = actual_bytes,
                        .line_file_offset = line_file_offset,
                    };
                }
                return;
            }

            if (self.deferred_short_tail != null) {
                try self.flushDeferredShortTailThen(metrics.line_count, bases, actual_bytes, line_file_offset);
                return;
            }

            const short_ok_if_final = metrics.is_uniform_width and
                metrics.line_bases > 0 and
                dense and
                bases <= metrics.line_bases and
                content_len <= self.line_builder.first_content_len and
                actual_bytes <= self.line_builder.first_actual_bytes + 1;
            if (short_ok_if_final) {
                self.deferred_short_tail = .{
                    .bases = bases,
                    .actual_bytes = actual_bytes,
                    .line_file_offset = line_file_offset,
                };
                return;
            }

            try self.activateSideTable(metrics.line_count - 1);
            try self.appendSideTableLine(bases, actual_bytes, line_file_offset);
        }

        fn noteSeqOffset(self: *Self, line_file_offset: u64) void {
            if (self.seq_offset_set) return;
            self.seq_offset = line_file_offset;
            self.seq_offset_set = true;
        }

        fn strideEligible(self: *const Self) bool {
            if (!self.active or self.in_header or self.in_pending_line) return false;
            if (self.side_table_active or self.deferred_short_tail != null) return false;
            if (self.line_builder.blank_after_bases) return false;
            const m = self.line_builder.metrics;
            if (m.line_count == 0 or !m.is_uniform_width) return false;
            if (!m.first_content_dense or !isDenseLineGeometry(m.line_bases, m.line_bytes)) return false;
            if (self.line_builder.first_actual_bytes != m.line_bytes) return false;
            if (self.line_builder.have_pending and
                (self.line_builder.pending_bases != m.line_bases or
                    self.line_builder.pending_actual_bytes != self.line_builder.first_actual_bytes))
            {
                return false;
            }
            return true;
        }

        fn ingestStrideRun(self: *Self, n_lines: u64) void {
            std.debug.assert(n_lines >= 1);
            std.debug.assert(self.line_builder.metrics.line_count >= 1);
            const bases = self.line_builder.metrics.line_bases;
            const line_bytes = self.line_builder.metrics.line_bytes;
            self.line_builder.metrics.seq_len += n_lines * @as(u64, bases);
            self.line_builder.metrics.line_count += n_lines;
            self.line_builder.pending_bases = bases;
            self.line_builder.pending_actual_bytes = line_bytes;
            self.line_builder.pending_content_len = self.line_builder.first_content_len;
            self.line_builder.have_pending = true;
            self.at_line_start = true;
        }

        fn ingestSequenceLine(
            self: *Self,
            data: []const u8,
            line_start: usize,
            line_end: usize,
            has_lf: bool,
            line_file_offset: u64,
        ) !void {
            var content_end = line_end;
            if (content_end > line_start and data[content_end - 1] == '\r') {
                content_end -= 1;
            }

            const actual_bytes: u32 = @intCast((if (has_lf) line_end + 1 else line_end) - line_start);
            const content_len: u32 = @intCast(content_end - line_start);
            const bases = countBasesInLineSlice(data, line_start, content_end);
            if (bases > 0) self.noteSeqOffset(line_file_offset);

            try self.recordSequenceLine(bases, actual_bytes, has_lf, content_len, line_file_offset);
            self.at_line_start = has_lf;
            self.in_pending_line = false;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_pre_lf_bytes = 0;
            self.pending_pre_lf_ends_with_cr = false;
        }

        fn appendPendingFragment(self: *Self, fragment: []const u8, fragment_file_offset: u64) !void {
            if (!self.in_pending_line) {
                self.in_pending_line = true;
                self.pending_line_bases = 0;
                self.pending_line_bytes = 0;
                self.pending_pre_lf_bytes = 0;
                self.pending_pre_lf_ends_with_cr = false;
                self.pending_line_file_offset = fragment_file_offset;
            }
            const fragment_bases = countBasesInLineSlice(fragment, 0, fragment.len);
            var diag = self.recordDiagnostics();
            diag.at.location.offset = self.pending_line_file_offset;
            const next_bases = std.math.add(u32, self.pending_line_bases, fragment_bases) catch
                return diag.lineLimit(.line_bases, @as(u64, self.pending_line_bases) + fragment_bases, .index_geometry);
            const next_bytes = std.math.add(u32, self.pending_line_bytes, @intCast(fragment.len)) catch
                return diag.lineLimit(.line_bytes, @as(u64, self.pending_line_bytes) + fragment.len, .index_geometry);
            self.pending_line_bases = next_bases;
            self.pending_line_bytes = next_bytes;
            self.pending_pre_lf_bytes = next_bytes;
            if (fragment.len > 0) {
                self.pending_pre_lf_ends_with_cr = fragment[fragment.len - 1] == '\r';
            }
        }

        fn commitPendingLine(self: *Self, has_lf: bool) !void {
            if (!self.in_pending_line and self.pending_line_bytes == 0) return;
            var diag = self.recordDiagnostics();
            diag.at.location.offset = self.pending_line_file_offset;
            if (!has_lf and self.pending_line_bytes == std.math.maxInt(u32)) {
                return diag.lineLimit(.line_bytes, @as(u64, self.pending_line_bytes) + 1, .synthetic_separator);
            }
            if (self.pending_line_bases > 0) self.noteSeqOffset(self.pending_line_file_offset);
            var content_len = if (has_lf) self.pending_pre_lf_bytes else self.pending_line_bytes;
            if (content_len > 0 and self.pending_pre_lf_ends_with_cr) {
                content_len -= 1;
            }
            try self.recordSequenceLine(
                self.pending_line_bases,
                self.pending_line_bytes,
                has_lf,
                content_len,
                self.pending_line_file_offset,
            );
            self.at_line_start = has_lf;
            self.in_pending_line = false;
            self.pending_line_bases = 0;
            self.pending_line_bytes = 0;
            self.pending_pre_lf_bytes = 0;
            self.pending_pre_lf_ends_with_cr = false;
        }

        fn finalizeHeader(self: *Self) void {
            self.in_header = false;
            self.at_line_start = true;
        }

        /// Emits the record in progress when it has bases; empty records are skipped.
        pub fn emitRecord(self: *Self, ctx: anytype, comptime emit: fn (@TypeOf(ctx), ScannedRecord) anyerror!void) !void {
            try self.commitPendingLine(false);
            const metrics = self.line_builder.finish();
            const seq_info = sequenceLengthInfoFromMetrics(metrics);
            if (seq_info.seq_len == 0) {
                self.deferred_short_tail = null;
                return;
            }

            var side_table_slice: []const u8 = &.{};
            if (seq_info.uses_uniform_formula) {
                self.deferred_short_tail = null;
            } else if (self.collect_side_tables) {
                if (self.deferred_short_tail) |short| {
                    if (!self.side_table_active) {
                        if (metrics.line_count < 1) return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
                        try self.activateSideTable(metrics.line_count - 1);
                        try self.appendSideTableLine(short.bases, short.actual_bytes, short.line_file_offset);
                    }
                    self.deferred_short_tail = null;
                } else if (!self.side_table_active) {
                    try self.activateSideTable(metrics.line_count);
                }
                if (!self.side_table_active) {
                    return self.recordDiagnostics().consistency(.missing_side_table, error.MissingSideTable);
                }
                @memcpy(self.side_table.items[0..@sizeOf(u64)], std.mem.asBytes(&metrics.line_count));
                side_table_slice = self.side_table.items;
            }
            try emit(ctx, .{
                .name = self.name.items,
                .header_offset = self.header_start_offset,
                .reference = &self.source_reference,
                .seq_offset = self.seq_offset,
                .seq_len = seq_info.seq_len,
                .line_bases = metrics.line_bases,
                .line_bytes = metrics.line_bytes,
                .uses_uniform_formula = seq_info.uses_uniform_formula,
                .side_table = side_table_slice,
                .nonuniform_reason = if (!metrics.first_content_dense) .sequence_whitespace else self.line_builder.nonuniform_reason,
            });
        }

        /// Reads the next window, which starts at `file_offset`. `max_name_len` limits identifiers.
        pub fn processChunk(
            state: *Self,
            data: []const u8,
            file_offset: u64,
            max_name_len: ?usize,
            ctx: anytype,
            comptime emit: fn (@TypeOf(ctx), ScannedRecord) anyerror!void,
        ) !void {
            var i: usize = 0;

            if (state.in_pending_line) {
                const nl = findLineFeed(data, 0, data.len);
                if (nl >= data.len) {
                    try state.appendPendingFragment(data, file_offset);
                    return;
                }

                try state.appendPendingFragment(data[0..nl], file_offset);
                state.pending_line_bytes = std.math.add(u32, state.pending_line_bytes, 1) catch {
                    var diag = state.diag;
                    diag.at.location.offset = state.pending_line_file_offset;
                    return diag.lineLimit(.line_bytes, @as(u64, state.pending_line_bytes) + 1, .index_geometry);
                };
                try state.commitPendingLine(true);
                i = nl + 1;
            }

            while (i < data.len) {
                const byte_offset = file_offset + @as(u64, @intCast(i));
                const byte = data[i];

                if (state.at_line_start and byte == '>') {
                    if (state.active) {
                        try state.emitRecord(ctx, emit);
                    }
                    try state.startHeader(byte_offset);
                    i += 1;
                    continue;
                }

                if (state.in_header) {
                    if (state.parsing_name) {
                        const name_end = findNameDelimiter(data, i);

                        const fragment = data[i..name_end];
                        if (max_name_len) |limit| {
                            if (state.name.items.len > limit or fragment.len > limit - state.name.items.len) {
                                return state.recordDiagnostics().identifierTooLong(@as(u128, state.name.items.len) + fragment.len, limit);
                            }
                        }
                        state.name.appendSlice(state.allocator, fragment) catch return state.diag.oom(.names, fragment.len);
                        i = name_end;
                        if (i >= data.len) return;
                        if (data[i] == '\n') {
                            state.finalizeHeader();
                            i += 1;
                            continue;
                        }
                        state.parsing_name = false;
                    }

                    const header_end = findLineFeed(data, i, data.len);
                    if (header_end >= data.len) return;
                    state.finalizeHeader();
                    i = header_end + 1;
                    continue;
                }

                if (!state.active) {
                    state.at_line_start = byte == '\n';
                    i += 1;
                    continue;
                }

                if (state.at_line_start and state.strideEligible()) {
                    const line_bases = state.line_builder.metrics.line_bases;
                    const line_bytes = state.line_builder.metrics.line_bytes;
                    const run_lines = countDenseLines(Bases, data, i, line_bases, line_bytes);
                    if (run_lines > 0) {
                        state.ingestStrideRun(run_lines);
                        i += run_lines * @as(usize, line_bytes);
                    }
                    if (i >= data.len) return;
                    if (data[i] == '>') continue;
                    if (data.len - i < line_bytes) {
                        if (findLineFeed(data, i, data.len) >= data.len) {
                            try state.appendPendingFragment(
                                data[i..],
                                file_offset + @as(u64, @intCast(i)),
                            );
                            return;
                        }
                    }
                }

                const line_start = i;
                const nl = findLineFeed(data, i, data.len);
                if (nl >= data.len) {
                    try state.appendPendingFragment(
                        data[line_start..],
                        file_offset + @as(u64, @intCast(line_start)),
                    );
                    return;
                }

                try state.ingestSequenceLine(
                    data,
                    line_start,
                    nl,
                    true,
                    file_offset + @as(u64, @intCast(line_start)),
                );
                i = nl + 1;
            }
        }
    };
}

const LineMetrics = struct {
    seq_len: u64 = 0,
    line_bases: u32 = 0,
    line_bytes: u32 = 0,
    is_uniform_width: bool = true,
    line_count: u64 = 0,
    first_content_dense: bool = false,
};

const SequenceLengthInfo = struct {
    seq_len: u64,
    uses_uniform_formula: bool,
};

fn findLineFeed(data: []const u8, start: usize, end: usize) usize {
    const relative = std.mem.findScalar(u8, data[start..end], '\n') orelse return end;
    return start + relative;
}

const LineMetricsBuilder = struct {
    metrics: LineMetrics = .{},
    nonuniform_reason: diagnostics.NonuniformReason = .nonuniform_layout,
    first_actual_bytes: u32 = 0,
    first_content_len: u32 = 0,
    pending_bases: u32 = 0,
    pending_actual_bytes: u32 = 0,
    pending_content_len: u32 = 0,
    have_pending: bool = false,
    blank_after_bases: bool = false,

    fn ingestLine(self: *LineMetricsBuilder, bases: u32, actual_bytes: u32, has_lf: bool, content_len: u32) void {
        if (bases == 0) {
            if (self.metrics.line_count > 0) {
                self.blank_after_bases = true;
            }
            return;
        }
        if (self.blank_after_bases) {
            self.metrics.is_uniform_width = false;
            self.nonuniform_reason = .interior_blank_line;
            self.blank_after_bases = false;
        }
        self.metrics.seq_len += bases;
        self.metrics.line_count += 1;
        if (bases != content_len) {
            if (self.metrics.is_uniform_width) self.nonuniform_reason = .sequence_whitespace;
            self.metrics.is_uniform_width = false;
        }

        if (self.metrics.line_count == 1) {
            self.first_actual_bytes = actual_bytes;
            self.first_content_len = content_len;
            self.metrics.line_bases = bases;
            self.metrics.line_bytes = if (has_lf) actual_bytes else actual_bytes + 1;
            const sep_len = if (has_lf) actual_bytes -| content_len else @as(u32, 1);
            self.metrics.first_content_dense = bases == content_len and (sep_len == 1 or sep_len == 2);
        } else if (self.have_pending) {
            if (self.pending_bases != self.metrics.line_bases or
                self.pending_actual_bytes != self.first_actual_bytes)
            {
                if (self.metrics.is_uniform_width) self.nonuniform_reason = if (self.pending_bases < self.metrics.line_bases) .interior_short_line else if (self.pending_bases > self.metrics.line_bases) .wider_line else .mixed_line_endings;
                self.metrics.is_uniform_width = false;
            }
        }

        self.pending_bases = bases;
        self.pending_actual_bytes = actual_bytes;
        self.pending_content_len = content_len;
        self.have_pending = true;
    }

    fn finish(self: *LineMetricsBuilder) LineMetrics {
        self.blank_after_bases = false;
        if (self.metrics.line_count > 1 and self.have_pending) {
            if (self.pending_bases > self.metrics.line_bases or
                self.pending_content_len > self.first_content_len)
            {
                if (self.metrics.is_uniform_width) self.nonuniform_reason = .wider_line;
                self.metrics.is_uniform_width = false;
            }
        }
        return self.metrics;
    }
};

fn sequenceLengthInfoFromMetrics(metrics: LineMetrics) SequenceLengthInfo {
    if (metrics.seq_len == 0) return .{ .seq_len = 0, .uses_uniform_formula = false };
    if (!metrics.is_uniform_width or metrics.line_bases == 0 or !metrics.first_content_dense) {
        return .{ .seq_len = metrics.seq_len, .uses_uniform_formula = false };
    }
    return .{
        .seq_len = metrics.seq_len,
        .uses_uniform_formula = isDenseLineGeometry(metrics.line_bases, metrics.line_bytes),
    };
}

// --- Identifier sets ---

/// Exact identifier set that copies names into 128 KiB blocks. `Diag` reports failures in the
/// caller's codes: it provides `at` and the `oom`, `capacity`, `fail`, `consistency`, and
/// `recordLimit` methods.
pub fn PackedNamesWith(comptime Diag: type) type {
    return struct {
        diag: Diag = .{},
        allocator: std.mem.Allocator,
        blocks: std.ArrayList([]u8) = .empty,
        position: usize = 0,
        map: Map,
        const BLOCK_BYTES = 128 * 1024;
        const Ref = struct { offset: u32, len: u32 };
        const Context = struct {
            owner: *PackedNames,
            pub fn hash(self: @This(), ref: Ref) u64 {
                return std.hash_map.hashString(self.owner.bytes(ref));
            }
            pub fn eql(self: @This(), a: Ref, b: Ref) bool {
                return a.len == b.len and std.mem.eql(u8, self.owner.bytes(a), self.owner.bytes(b));
            }
        };
        const Lookup = struct {
            owner: *PackedNames,
            pub fn hash(_: @This(), name: []const u8) u64 {
                return std.hash_map.hashString(name);
            }
            pub fn eql(self: @This(), name: []const u8, ref: Ref) bool {
                return name.len == ref.len and std.mem.eql(u8, name, self.owner.bytes(ref));
            }
        };
        const Map = std.HashMap(Ref, void, Context, std.hash_map.default_max_load_percentage);

        pub fn init(allocator: std.mem.Allocator) @This() {
            return .{ .allocator = allocator, .map = Map.initContext(allocator, .{ .owner = undefined }) };
        }

        pub fn ensureCapacity(self: *@This(), capacity: u32) !void {
            self.map.ctx.owner = self;
            self.map.ensureTotalCapacity(capacity) catch return self.diag.oom(.dedup_table, null);
        }

        pub fn deinit(self: *@This()) void {
            self.map.deinit();
            for (self.blocks.items) |block| self.allocator.free(block);
            self.blocks.deinit(self.allocator);
        }

        fn bytes(self: *@This(), ref: Ref) []const u8 {
            if (ref.len == 0) return &.{};
            return self.blocks.items[ref.offset / BLOCK_BYTES][ref.offset % BLOCK_BYTES ..][0..ref.len];
        }

        pub fn observeAt(self: *@This(), name: []const u8, source_offset: u64) !bool {
            _ = source_offset;
            self.map.ctx.owner = self;
            const entry = self.map.getOrPutAdapted(name, Lookup{ .owner = self }) catch return self.diag.oom(.dedup_table, null);
            if (entry.found_existing) return true;
            errdefer self.map.removeByPtr(entry.key_ptr);
            const remaining = BLOCK_BYTES - self.position % BLOCK_BYTES;
            if (name.len > remaining) {
                self.position = std.math.add(usize, self.position, remaining) catch return self.diag.capacity(.name_offset, std.math.maxInt(usize), @as(u128, self.position) + remaining, error.OutOfMemory);
            }
            const offset = std.math.cast(u32, self.position) orelse return self.diag.capacity(.name_offset, std.math.maxInt(u32), self.position, error.OutOfMemory);
            const next_position = std.math.add(usize, self.position, name.len) catch return self.diag.capacity(.name_offset, std.math.maxInt(usize), @as(u128, self.position) + name.len, error.OutOfMemory);
            if (name.len != 0 and self.position / BLOCK_BYTES == self.blocks.items.len) {
                const block = self.allocator.alloc(u8, BLOCK_BYTES) catch return self.diag.oom(.names, BLOCK_BYTES);
                self.blocks.append(self.allocator, block) catch {
                    self.allocator.free(block);
                    return self.diag.oom(.names, null);
                };
            }
            const reference = Ref{ .offset = offset, .len = @intCast(name.len) };
            if (name.len != 0) @memcpy(self.blocks.items[self.position / BLOCK_BYTES][self.position % BLOCK_BYTES ..][0..name.len], name);
            self.position = next_position;
            entry.key_ptr.* = reference;
            return false;
        }

        const PackedNames = @This();
    };
}

const CACHE_BYTES = 512 * 1024;
const CACHE_BLOCK_BYTES: usize = 128 * 1024;

/// Exact identifier set over a seekable source: a fingerprint and a source offset per name, with
/// early names cached. Names are re-read from the source only on a fingerprint match. `Diag` as
/// for `PackedNamesWith`.
pub fn SourceDedupWith(comptime Diag: type) type {
    return struct {
        diag: Diag = .{},
        allocator: std.mem.Allocator,
        io: std.Io,
        file: std.Io.File,
        tables: union(enum) { narrow: Tables(u32), wide: Tables(u64) },
        chunk: []const u8 = &.{},
        chunk_offset: u64 = 0,
        cache: std.ArrayList(u8) = .empty,
        scratch: std.ArrayList(u8) = .empty,
        cache_full: bool = false,
        allocation_failed: bool = false,
        source_match: bool = false,
        hot_blocks: std.ArrayList([]u8) = .empty,
        hot_len: usize = 0,
        read_failed: bool = false,

        pub fn init(allocator: std.mem.Allocator, io: std.Io, file: std.Io.File, size: u64) @This() {
            return .{
                .allocator = allocator,
                .io = io,
                .file = file,
                .tables = if (size <= std.math.maxInt(u31) - CACHE_BYTES) .{ .narrow = Tables(u32).init(allocator) } else .{ .wide = Tables(u64).init(allocator) },
            };
        }

        pub fn deinit(self: *@This()) void {
            switch (self.tables) {
                .narrow => |*set| set.deinit(),
                .wide => |*set| set.deinit(),
            }
            self.cache.deinit(self.allocator);
            self.scratch.deinit(self.allocator);
            for (self.hot_blocks.items) |block| self.allocator.free(block);
            self.hot_blocks.deinit(self.allocator);
        }

        pub fn setChunk(self: *@This(), data: []const u8, offset: u64) void {
            self.chunk = data;
            self.chunk_offset = offset;
        }

        pub fn observe(self: *@This(), name: []const u8, offset: u64) !bool {
            return self.observeFingerprint(name, offset, fingerprintOf(name));
        }

        /// The fingerprint `observeFingerprint` takes for `name`. The table slots its probe starts
        /// at are prefetched, so a caller that observes later finds them in cache.
        pub fn prefetch(self: *@This(), name: []const u8) u32 {
            const fingerprint = fingerprintOf(name);
            switch (self.tables) {
                inline else => |*set| set.prefetch(fingerprint),
            }
            return fingerprint;
        }

        pub fn observeFingerprint(self: *@This(), name: []const u8, offset: u64, fingerprint: u32) !bool {
            return switch (self.tables) {
                .narrow => |*set| if (offset >= Tables(u32).CACHE_FLAG)
                    self.diag.capacity(.name_cache_offset, Tables(u32).CACHE_FLAG - 1, offset, error.SourceChanged)
                else
                    set.observeWithFingerprint(self, name, @intCast(offset), fingerprint),
                .wide => |*set| if (offset >= Tables(u64).CACHE_FLAG)
                    self.diag.capacity(.name_cache_offset, Tables(u64).CACHE_FLAG - 1, offset, error.SourceChanged)
                else
                    set.observeWithFingerprint(self, name, offset, fingerprint),
            };
        }

        fn fingerprintOf(name: []const u8) u32 {
            return @truncate(std.hash_map.hashString(name));
        }

        const CacheEncoding = struct { header: [7]u8, size: usize, middle: []const u8 };

        fn encodeCacheName(name: []const u8, prefix: usize, suffix: usize) CacheEncoding {
            const middle = name[prefix .. name.len - suffix];
            const short = prefix <= 255 and suffix <= 255 and middle.len <= 255;
            const header_size: usize = if (short) 4 else 7;
            var header: [7]u8 = undefined;
            if (name.len <= 253 and name.len + 1 <= header_size + middle.len) {
                header[0] = @intCast(name.len);
                return .{ .header = header, .size = 1, .middle = name };
            }
            if (short) {
                header[0] = 254;
                header[1] = @intCast(prefix);
                header[2] = @intCast(suffix);
                header[3] = @intCast(middle.len);
            } else {
                header[0] = 255;
                std.mem.writeInt(u16, header[1..3], @intCast(prefix), .little);
                std.mem.writeInt(u16, header[3..5], @intCast(suffix), .little);
                std.mem.writeInt(u16, header[5..7], @intCast(middle.len), .little);
            }
            return .{ .header = header, .size = header_size, .middle = middle };
        }

        fn cacheAnchor(self: *@This(), base: usize) []const u8 {
            return switch (self.cache.items[base]) {
                0...253 => |len| self.cache.items[base + 1 ..][0..len],
                254 => self.cache.items[base + 4 ..][0..self.cache.items[base + 3]],
                255 => self.cache.items[base + 7 ..][0..std.mem.readInt(u16, self.cache.items[base + 5 ..][0..2], .little)],
            };
        }

        fn remember(self: *@This(), name: []const u8, offset: u64, flag: u64) !u64 {
            if (self.cache_full or self.cache.items.len >= CACHE_BYTES) return offset;
            var prefix: usize = 0;
            var suffix: usize = 0;
            if (name.len >= 16 and self.cache.items.len % CACHE_BLOCK_BYTES != 0) {
                const base = (self.cache.items.len - 1) & ~(CACHE_BLOCK_BYTES - 1);
                const anchor = self.cacheAnchor(base);
                prefix = commonPrefix(name, anchor);
                suffix = commonSuffix(name, anchor, prefix);
            }
            var encoded = encodeCacheName(name, prefix, suffix);
            var start = self.cache.items.len;
            const remaining = CACHE_BLOCK_BYTES - (start % CACHE_BLOCK_BYTES);
            if (encoded.size + encoded.middle.len > remaining) {
                start += remaining;
                encoded = encodeCacheName(name, 0, 0);
            }
            if (start + encoded.size + encoded.middle.len > CACHE_BYTES) {
                self.cache_full = true;
                return offset;
            }
            if (self.cache.capacity == 0) self.cache.ensureTotalCapacity(self.allocator, CACHE_BYTES) catch return self.diag.oom(.name_cache, CACHE_BYTES);
            self.cache.items.len = start;
            self.cache.appendSliceAssumeCapacity(encoded.header[0..encoded.size]);
            self.cache.appendSliceAssumeCapacity(encoded.middle);
            return flag | start;
        }

        fn sameCached(self: *@This(), name: []const u8, start: usize) bool {
            const tag = self.cache.items[start];
            if (tag <= 253) {
                return name.len == tag and std.mem.eql(u8, name, self.cache.items[start + 1 ..][0..tag]);
            }
            const prefix: usize = if (tag == 254) self.cache.items[start + 1] else std.mem.readInt(u16, self.cache.items[start + 1 ..][0..2], .little);
            const suffix: usize = if (tag == 254) self.cache.items[start + 2] else std.mem.readInt(u16, self.cache.items[start + 3 ..][0..2], .little);
            const middle_len: usize = if (tag == 254) self.cache.items[start + 3] else std.mem.readInt(u16, self.cache.items[start + 5 ..][0..2], .little);
            const header_len: usize = if (tag == 254) 4 else 7;
            if (name.len != prefix + suffix + middle_len) return false;
            const anchor = self.cacheAnchor(start & ~(CACHE_BLOCK_BYTES - 1));
            return std.mem.eql(u8, name[0..prefix], anchor[0..prefix]) and
                std.mem.eql(u8, name[prefix..][0..middle_len], self.cache.items[start + header_len ..][0..middle_len]) and
                std.mem.eql(u8, name[name.len - suffix ..], anchor[anchor.len - suffix ..]);
        }

        fn appendHot(self: *@This(), input: []const u8) !void {
            var remaining = input;
            while (remaining.len != 0) {
                const block_index = self.hot_len / CACHE_BLOCK_BYTES;
                const within = self.hot_len % CACHE_BLOCK_BYTES;
                if (block_index == self.hot_blocks.items.len) {
                    const block = self.allocator.alloc(u8, CACHE_BLOCK_BYTES) catch return self.diag.oom(.retained_names, CACHE_BLOCK_BYTES);
                    errdefer self.allocator.free(block);
                    self.hot_blocks.append(self.allocator, block) catch return self.diag.oom(.retained_names, null);
                }
                const n = @min(remaining.len, CACHE_BLOCK_BYTES - within);
                @memcpy(self.hot_blocks.items[block_index][within..][0..n], remaining[0..n]);
                self.hot_len += n;
                remaining = remaining[n..];
            }
        }

        fn rememberHot(self: *@This(), name: []const u8, flag: u64) !u64 {
            if (self.hot_len % CACHE_BLOCK_BYTES == CACHE_BLOCK_BYTES - 1) try self.appendHot(&.{0});
            const start = self.hot_len;
            if (@as(u128, CACHE_BYTES) + start >= flag) return self.diag.capacity(.name_cache_offset, flag - 1, @as(u128, CACHE_BYTES) + start, error.OutOfMemory);
            var header: [2]u8 = undefined;
            std.mem.writeInt(u16, &header, @intCast(name.len), .little);
            try self.appendHot(&header);
            try self.appendHot(name);
            return flag | (CACHE_BYTES + start);
        }

        fn sameHot(self: *@This(), name: []const u8, start: usize) bool {
            const header_block = self.hot_blocks.items[start / CACHE_BLOCK_BYTES];
            const len = std.mem.readInt(u16, header_block[start % CACHE_BLOCK_BYTES ..][0..2], .little);
            if (len != name.len) return false;
            if (name.len == 0) return true;
            const bytes_start = start + 2;
            const within = bytes_start % CACHE_BLOCK_BYTES;
            const n = @min(name.len, CACHE_BLOCK_BYTES - within);
            const block_index = bytes_start / CACHE_BLOCK_BYTES;
            if (!std.mem.eql(u8, name[0..n], self.hot_blocks.items[block_index][within..][0..n])) return false;
            return n == name.len or std.mem.eql(u8, name[n..], self.hot_blocks.items[block_index + 1][0 .. name.len - n]);
        }

        pub fn sameName(self: *@This(), name: []const u8, offset: u64, flag: u64) bool {
            if (offset & flag != 0) {
                const start: usize = @intCast(offset & ~flag);
                return if (start < CACHE_BYTES) self.sameCached(name, start) else self.sameHot(name, start - CACHE_BYTES);
            }
            const needed = name.len + 1;
            if (offset >= self.chunk_offset and offset - self.chunk_offset <= self.chunk.len) {
                const start: usize = @intCast(offset - self.chunk_offset);
                if (needed <= self.chunk.len - start) return sameBytes(name, self.chunk[start..][0..needed]);
            }
            self.scratch.resize(self.allocator, needed) catch {
                switch (self.diag.oom(.retained_names, needed)) {
                    else => {},
                }
                self.allocation_failed = true;
                return false;
            };
            const read = self.file.readPositionalAll(self.io, self.scratch.items[0..needed], offset) catch |err| {
                var diag = self.diag;
                diag.at.location.offset = offset;
                switch (diag.fail(diagnostics.readFailure(needed, null), .read, .fasta, err, error.SourceReadFailed)) {
                    else => {},
                }
                self.read_failed = true;
                return false;
            };
            if (read != needed) {
                for (self.scratch.items[0..read]) |byte| {
                    if (isNameDelimiter(byte)) return false;
                }
                var diag = self.diag;
                diag.at.location.offset = offset;
                switch (diag.fail(diagnostics.readFailure(needed, read), .read, .fasta, null, error.SourceReadFailed)) {
                    else => {},
                }
                self.read_failed = true;
                return false;
            }
            const matches = sameBytes(name, self.scratch.items[0..needed]);
            if (matches) self.source_match = true;
            return matches;
        }

        fn sameBytes(name: []const u8, bytes: []const u8) bool {
            return isNameDelimiter(bytes[name.len]) and std.mem.eql(u8, name, bytes[0..name.len]);
        }

        const SourceDedup = @This();

        pub fn Tables(comptime Offset: type) type {
            return struct {
                pub const CACHE_FLAG: u64 = @as(u64, 1) << (@bitSizeOf(Offset) - 1);
                const Key = struct { fingerprint: u32, offset: Offset };
                const BLOCK_ENTRIES = 4096;
                const Vec = @Vector(16, u8);
                allocator: std.mem.Allocator,
                buckets: []u32 = &.{},
                metadata: []u8 = &.{},
                blocks: std.ArrayList([]Key) = .empty,
                len: u32 = 0,

                pub fn init(allocator: std.mem.Allocator) @This() {
                    return .{ .allocator = allocator };
                }

                pub fn deinit(self: *@This()) void {
                    self.allocator.free(self.buckets);
                    self.allocator.free(self.metadata);
                    for (self.blocks.items) |block| self.allocator.free(block);
                    self.blocks.deinit(self.allocator);
                }

                pub fn count(self: *const @This()) u32 {
                    return self.len;
                }

                pub fn keyAt(self: *@This(), index: u32) *Key {
                    return &self.blocks.items[index / BLOCK_ENTRIES][index % BLOCK_ENTRIES];
                }

                fn tag(fingerprint: u32) u8 {
                    return 0x80 | @as(u8, @truncate(fingerprint >> 25));
                }

                fn emptySlot(metadata: []const u8, fingerprint: u32) usize {
                    const mask = metadata.len - 1;
                    var group = spread(fingerprint) & mask & ~@as(usize, 15);
                    while (true) : (group = (group + 16) & mask) {
                        const bytes: Vec = metadata[group..][0..16].*;
                        if (std.simd.firstTrue(bytes == @as(Vec, @splat(0)))) |index| return group + index;
                    }
                }

                fn grow(self: *@This(), diag: Diag) !void {
                    var required: u128 = if (self.buckets.len == 0) 512 else @as(u128, self.buckets.len) * 2;
                    while ((@as(u128, self.len) + 1) * 5 > required * 4) required *= 2;
                    const capacity = std.math.cast(usize, required) orelse return diag.capacity(.name_map, std.math.maxInt(usize), required, error.OutOfMemory);
                    if (required * @sizeOf(u32) > std.math.maxInt(usize)) return diag.capacity(.name_map, std.math.maxInt(usize), required * @sizeOf(u32), error.OutOfMemory);
                    self.allocator.free(self.buckets);
                    self.allocator.free(self.metadata);
                    self.buckets = &.{};
                    self.metadata = &.{};
                    const metadata = self.allocator.alloc(u8, capacity) catch return diag.oom(.dedup_table, capacity);
                    errdefer self.allocator.free(metadata);
                    @memset(metadata, 0);
                    const buckets = self.allocator.alloc(u32, capacity) catch return diag.oom(.dedup_table, capacity * @sizeOf(u32));
                    for (0..self.len) |entry_index| {
                        const key = self.keyAt(@intCast(entry_index));
                        const position = emptySlot(metadata, key.fingerprint);
                        buckets[position] = @intCast(entry_index + 1);
                        metadata[position] = tag(key.fingerprint);
                    }
                    self.buckets = buckets;
                    self.metadata = metadata;
                }

                pub fn observe(self: *@This(), source: *SourceDedup, name: []const u8, offset: Offset) !bool {
                    return self.observeWithFingerprint(source, name, offset, @truncate(std.hash_map.hashString(name)));
                }

                pub fn observeWithFingerprint(self: *@This(), source: *SourceDedup, name: []const u8, offset: Offset, fingerprint: u32) !bool {
                    if (self.buckets.len == 0) try self.grow(source.diag);
                    const mask = self.buckets.len - 1;
                    var group = spread(fingerprint) & mask & ~@as(usize, 15);
                    const wanted: Vec = @splat(tag(fingerprint));
                    source.source_match = false;
                    var position: usize = undefined;
                    while (true) : (group = (group + 16) & mask) {
                        const bytes: Vec = self.metadata[group..][0..16].*;
                        var matches: u16 = @bitCast(bytes == wanted);
                        while (matches != 0) {
                            const candidate = group + @ctz(matches);
                            matches &= matches - 1;
                            const key = self.keyAt(self.buckets[candidate] - 1);
                            if (key.fingerprint != fingerprint) continue;
                            if (source.sameName(name, key.offset, CACHE_FLAG)) {
                                if (source.source_match) key.offset = @intCast(try source.rememberHot(name, CACHE_FLAG));
                                return true;
                            }
                            if (source.read_failed) return error.SourceReadFailed;
                            if (source.allocation_failed) return error.OutOfMemory;
                        }
                        if (std.simd.firstTrue(bytes == @as(Vec, @splat(0)))) |empty| {
                            position = group + empty;
                            break;
                        }
                    }
                    if (self.len == std.math.maxInt(u32)) {
                        @branchHint(.cold);
                        return source.diag.recordLimit(@as(u64, self.len) + 1);
                    }
                    if ((@as(u64, self.len) + 1) * 5 > @as(u64, self.buckets.len) * 4) {
                        try self.grow(source.diag);
                        position = emptySlot(self.metadata, fingerprint);
                    }
                    if (self.len == self.blocks.items.len * BLOCK_ENTRIES) {
                        const block = self.allocator.alloc(Key, BLOCK_ENTRIES) catch return source.diag.oom(.dedup_table, BLOCK_ENTRIES * @sizeOf(Key));
                        errdefer self.allocator.free(block);
                        self.blocks.append(self.allocator, block) catch return source.diag.oom(.dedup_table, null);
                    }
                    const reference = try source.remember(name, offset, CACHE_FLAG);
                    self.keyAt(self.len).* = .{ .fingerprint = fingerprint, .offset = @intCast(reference) };
                    self.buckets[position] = self.len + 1;
                    self.metadata[position] = tag(fingerprint);
                    self.len += 1;
                    return false;
                }

                fn prefetch(self: *const @This(), fingerprint: u32) void {
                    if (self.buckets.len == 0) return;
                    const group = spread(fingerprint) & (self.buckets.len - 1) & ~@as(usize, 15);
                    @prefetch(self.metadata.ptr + group, .{});
                    @prefetch(self.buckets.ptr + group, .{});
                }

                fn spread(fingerprint: u32) usize {
                    return @truncate(@as(u64, fingerprint) *% 0x9e3779b97f4a7c15);
                }
            };
        }
    };
}

pub fn commonPrefix(a: []const u8, b: []const u8) usize {
    const Vec = @Vector(16, u8);
    const limit = @min(a.len, b.len);
    var position: usize = 0;
    while (limit - position >= 16) : (position += 16) {
        const av: Vec = a[position..][0..16].*;
        const bv: Vec = b[position..][0..16].*;
        if (std.simd.firstTrue(av != bv)) |different| return position + different;
    }
    while (position < limit and a[position] == b[position]) : (position += 1) {}
    return position;
}

pub fn commonSuffix(a: []const u8, b: []const u8, prefix: usize) usize {
    const Vec = @Vector(16, u8);
    const limit = @min(a.len, b.len) - prefix;
    var matched: usize = 0;
    while (limit - matched >= 16) : (matched += 16) {
        const av: Vec = a[a.len - matched - 16 ..][0..16].*;
        const bv: Vec = b[b.len - matched - 16 ..][0..16].*;
        if (std.simd.lastTrue(av != bv)) |different| return matched + 15 - different;
    }
    while (matched < limit and a[a.len - matched - 1] == b[b.len - matched - 1]) : (matched += 1) {}
    return matched;
}

const SIMD_CHUNK_SIZE = 32;
const SimdVec = @Vector(SIMD_CHUNK_SIZE, u8);

const DENSE_VECTOR_LEN = std.simd.suggestVectorLength(u8) orelse 1;
const DENSE_BLOCK_LINES = 256;

fn countDenseLinesWith(comptime vector_len: usize, comptime Symbols: type, data: []const u8, start: usize, line_bases_u32: u32, line_bytes_u32: u32) usize {
    const line_bases: usize = line_bases_u32;
    const line_bytes: usize = line_bytes_u32;
    if (!isDenseLineGeometry(line_bases_u32, line_bytes_u32)) return 0;
    var cursor = start;
    var accepted_lines: usize = 0;
    while (data.len - cursor >= line_bytes) {
        const block_lines = @min((data.len - cursor) / line_bytes, DENSE_BLOCK_LINES);
        const valid_lines = denseBlockLines(vector_len, Symbols, data, cursor, block_lines, line_bases, line_bytes);
        accepted_lines += valid_lines;
        if (valid_lines != block_lines) return accepted_lines;
        cursor += block_lines * line_bytes;
    }
    return accepted_lines;
}

fn denseBlockLines(comptime vector_len: usize, comptime Symbols: type, data: []const u8, start: usize, line_count: usize, line_bases: usize, line_bytes: usize) usize {
    const sep_len = line_bytes - line_bases;
    const candidate_lines = switch (sep_len) {
        1 => separatedLines("\n", data[start..][0 .. line_count * line_bytes], line_bases, line_bytes),
        2 => separatedLines("\r\n", data[start..][0 .. line_count * line_bytes], line_bases, line_bytes),
        else => unreachable,
    };
    if (candidate_lines == 0) return 0;
    const block_end = start + candidate_lines * line_bytes;
    if (Symbols.countOutside(vector_len, data[start..block_end]) == candidate_lines * sep_len) return candidate_lines;
    var line_index: usize = 0;
    while (line_index < candidate_lines) : (line_index += 1) {
        const line_start = start + line_index * line_bytes;
        if (Symbols.countOutside(vector_len, data[line_start .. line_start + line_bases]) != 0) return line_index;
    }
    return 0;
}

// Lines of `block` that start with no `>` and end in `separator`; `block` holds whole lines.
fn separatedLines(comptime separator: []const u8, block: []const u8, line_bases: usize, line_bytes: usize) usize {
    var line_start: usize = 0;
    while (line_start != block.len) : (line_start += line_bytes) {
        if (block[line_start] == '>') break;
        if (!std.mem.eql(u8, block[line_start + line_bases ..][0..separator.len], separator)) break;
    }
    return line_start / line_bytes;
}

fn lineSliceHasWhitespace(data: []const u8, start: usize, end: usize) bool {
    const space: SimdVec = @splat(' ');
    if (end - start < SIMD_CHUNK_SIZE) {
        for (data[start..end]) |byte| if (byte <= ' ') return true;
        return false;
    }
    var pos = start;
    while (end - pos > SIMD_CHUNK_SIZE) : (pos += SIMD_CHUNK_SIZE) {
        const chunk: SimdVec = data[pos..][0..SIMD_CHUNK_SIZE].*;
        if (@reduce(.Or, chunk <= space)) return true;
    }
    const last: SimdVec = data[end - SIMD_CHUNK_SIZE ..][0..SIMD_CHUNK_SIZE].*;
    return @reduce(.Or, last <= space);
}

fn countOutsideWith(comptime vector_len: usize, data: []const u8, comptime inside: anytype) usize {
    return countVectorized(vector_len, data, struct {
        fn outside(comptime Vec: type, bytes: Vec) @Vector(@typeInfo(Vec).vector.len, bool) {
            return !inside(Vec, bytes);
        }
    }.outside);
}

fn anyOutsideWith(comptime vector_len: usize, data: []const u8, comptime inside: anytype) bool {
    return anyVectorized(vector_len, data, struct {
        fn outside(comptime Vec: type, bytes: Vec) @Vector(@typeInfo(Vec).vector.len, bool) {
            return !inside(Vec, bytes);
        }
    }.outside);
}

// Whether `matches` selects any byte: masks ORed across vectors, the last partial vector read as the
// slice's final full vector.
fn anyVectorized(comptime vector_len: usize, data: []const u8, comptime matches: anytype) bool {
    const Vec = @Vector(vector_len, u8);
    if (data.len < vector_len) {
        for (data) |byte| {
            if (matches(@Vector(1, u8), .{byte})[0]) return true;
        }
        return false;
    }
    var found: @Vector(vector_len, bool) = @splat(false);
    var pos: usize = 0;
    while (data.len - pos > vector_len) : (pos += vector_len) found = found | matches(Vec, data[pos..][0..vector_len].*);
    found = found | matches(Vec, data[data.len - vector_len ..][0..vector_len].*);
    return @reduce(.Or, found);
}

// Counts the bytes `matches` selects: u8 lanes summed per group, so no lane overflows, and a last
// partial vector read as the slice's final full vector with its already-counted lanes masked off.
fn countVectorized(comptime vector_len: usize, data: []const u8, comptime matches: anytype) usize {
    const Vec = @Vector(vector_len, u8);
    const Wide = @Vector(vector_len, u16);
    const one: Vec = @splat(1);
    const zero: Vec = @splat(0);
    if (data.len < vector_len) {
        var count: usize = 0;
        for (data) |byte| count += @intFromBool(matches(@Vector(1, u8), .{byte})[0]);
        return count;
    }
    const group = 128;
    const step = 4 * vector_len;
    var count: usize = 0;
    var pos: usize = 0;
    while (data.len - pos >= vector_len * group) {
        var lanes: Vec = zero;
        for (0..group / 4) |_| {
            inline for (0..4) |k| lanes += @select(u8, matches(Vec, data[pos + k * vector_len ..][0..vector_len].*), one, zero);
            pos += step;
        }
        count += @reduce(.Add, @as(Wide, lanes));
    }
    var lanes: Vec = zero;
    while (data.len - pos >= step) : (pos += step) {
        inline for (0..4) |k| lanes += @select(u8, matches(Vec, data[pos + k * vector_len ..][0..vector_len].*), one, zero);
    }
    while (data.len - pos >= vector_len) : (pos += vector_len) {
        lanes += @select(u8, matches(Vec, data[pos..][0..vector_len].*), one, zero);
    }
    const rest = data.len - pos;
    if (rest != 0) {
        const fresh = std.simd.iota(u8, vector_len) >= @as(Vec, @splat(@intCast(vector_len - rest)));
        lanes += @select(u8, matches(Vec, data[data.len - vector_len ..][0..vector_len].*), one, zero) & @select(u8, fresh, one, zero);
    }
    return count + @reduce(.Add, @as(Wide, lanes));
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

test "[property] - [fasta scanner]: whitespace counts match scalar counts" {
    var data: [12289]u8 = undefined;
    const lengths = [_]usize{ 0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 2047, 2048, 2049, 4095, 4096, 4097, 8191, 8192, 8193, data.len };
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
    var line = [_]u8{'A'} ** 80;
    for ([_]usize{ 0, 1, 7 }) |start| {
        for (start..line.len + 1) |end| {
            try std.testing.expect(!lineSliceHasWhitespace(&line, start, end));
            for (start..end) |at| {
                line[at] = ' ';
                try std.testing.expect(lineSliceHasWhitespace(&line, start, end));
                inline for (.{ 1, 16, 32, 64 }) |width| {
                    try std.testing.expectEqual(@as(usize, 1), countWhitespace(width, line[start..end]));
                }
                line[at] = 'A';
            }
        }
    }
}

test "[property] - [fasta scanner]: name delimiters match scalar scanning" {
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

test "[property] - [fasta scanner]: matches scalar structural boundaries" {
    const line_bases: u32 = 4;
    const line_count = 520;
    const target_line = 256;
    const cases = [_][]const u8{
        "AAAA\n",
        "AAAA\r\n",
    };
    const mutations = [_]u8{ ' ', '\t', '\r', '\n', 0, '>' };

    var storage: [line_count * 6]u8 = undefined;
    for (cases) |line| {
        const line_bytes: u32 = @intCast(line.len);
        const block = storage[0 .. line_count * line.len];
        for (0..line_count) |line_index| {
            @memcpy(block[line_index * line.len ..][0..line.len], line);
        }

        for (0..line.len + 1) |line_offset| {
            const mutation_offset = target_line * line.len + line_offset;
            const original = block[mutation_offset];
            for (mutations) |mutation| {
                block[mutation_offset] = mutation;

                const native = countDenseLines(Bases, block, 0, line_bases, line_bytes);
                const scalar = reference: {
                    for (0..line_count) |line_index| {
                        const row = block[line_index * line.len ..][0..line.len];
                        if (row[0] == '>') break :reference line_index;
                        for (row[0..line_bases]) |byte| {
                            if (byte <= ' ') break :reference line_index;
                        }
                        if (!std.mem.eql(u8, row[line_bases..], line[line_bases..])) break :reference line_index;
                    }
                    break :reference line_count;
                };

                try std.testing.expectEqual(scalar, native);
                const at_next_line = line_offset == line.len;
                const invalid = if (line_offset < line_bases)
                    mutation <= ' ' or (line_offset == 0 and mutation == '>')
                else if (!at_next_line)
                    mutation != line[line_offset]
                else
                    mutation <= ' ' or mutation == '>';
                const stop_line = target_line + @as(usize, @intFromBool(at_next_line));
                try std.testing.expectEqual(if (invalid) stop_line else line_count, native);
            }
            block[mutation_offset] = original;
        }
    }
}

test "[regression] - [fasta scanner]: checks whitespace before a later structural stop" {
    const cases = [_][]const u8{
        "AAAA\nAA A\n>b x\n",
        "AAAA\nAA A\nAAAAx",
        "AAAA\r\nAA A\r\n>b x\r\n",
        "AAAA\r\nAA A\r\nAAAAx\n",
    };
    for (cases) |data| {
        const line_bytes: u32 = if (data[4] == '\r') 6 else 5;
        try std.testing.expectEqual(@as(usize, 1), countDenseLinesWith(1, Bases, data, 0, 4, line_bytes));
        try std.testing.expectEqual(@as(usize, 1), countDenseLines(Bases, data, 0, 4, line_bytes));
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

test "[property] - [complement]: vector blocks match the table for every byte" {
    const plain = "ACGTNacgtn-\n*" ** 8;
    var plain_out: [plain.len]u8 = undefined;
    complementInto(&plain_out, plain);
    for (plain, plain_out) |input, actual| try std.testing.expectEqual(complement(input), actual);
    var src: [256 + 31]u8 = undefined;
    for (&src, 0..) |*byte, i| byte.* = @truncate(i);
    var dst: [src.len]u8 = undefined;
    complementInto(&dst, &src);
    for (src, dst) |input, actual| try std.testing.expectEqual(complement(input), actual);
    for (0..src.len + 1) |len| {
        for ([_]bool{ false, true }) |do_complement| {
            reverseInto(dst[0..len], src[0..len], do_complement);
            for (dst[0..len], 0..) |actual, i| {
                const input = src[len - 1 - i];
                try std.testing.expectEqual(if (do_complement) complement(input) else input, actual);
            }
        }
    }
}

test "[property] - [complement]: matches scalar mapping across chunk boundaries" {
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

test "[property] - [symbol sets]: anyOutside agrees with countOutside" {
    var prng = std.Random.DefaultPrng.init(0xa11);
    const random = prng.random();
    var data: [200]u8 = undefined;
    const alphabet = "ACGTNacgtnU*-Xx\n \t\r";
    for (0..4000) |_| {
        const len = random.uintAtMost(usize, data.len);
        const stray = random.uintLessThan(u8, 4) == 0;
        for (data[0..len]) |*byte| byte.* = if (stray) alphabet[random.uintLessThan(usize, alphabet.len)] else "ACGTacgt"[random.uintLessThan(usize, 8)];
        const bytes = data[0..len];
        inline for (.{ PlainNucleotides, PlainDna, PlainRna, ProteinSymbols, Bases }) |Symbols| {
            inline for (.{ 1, 16, 32 }) |vector_len| {
                try std.testing.expectEqual(Symbols.countOutside(vector_len, bytes) != 0, Symbols.anyOutside(vector_len, bytes));
            }
        }
    }
}
