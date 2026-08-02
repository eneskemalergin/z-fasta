//! Product build: z-fasta exe (ReleaseFast + strip by default) and unit/integration tests.

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    // Ship path is ReleaseFast + strip. ReleaseSmall slows .text; reject it.
    if (optimize == .ReleaseSmall) {
        std.debug.print(
            "error: ReleaseSmall is unsupported; use -Doptimize=ReleaseFast (strips by default)\n",
            .{},
        );
        std.process.exit(1);
    }
    // Zig auto-strips only ReleaseSmall; RF needs an explicit strip flag.
    const strip = b.option(bool, "strip", "Strip debug info") orelse (optimize != .Debug);

    const exe = b.addExecutable(.{
        .name = "z-fasta",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .strip = strip,
        }),
    });

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run z-fasta");
    run_step.dependOn(&run_cmd.step);

    // Test modules stay unstripped; ship strip applies to the installed exe only.
    const main_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    const test_index_module = b.createModule(.{
        .root_source_file = b.path("tests/test_index.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "main", .module = main_module }},
    });
    // Required for this suite on current Zig/OS (same as test_validate).
    test_index_module.link_libc = true;
    const run_test_index = b.addRunArtifact(b.addTest(.{ .root_module = test_index_module }));

    const run_test_get = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_get.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    }));
    // CLI subprocess tests spawn zig-out/bin/z-fasta; keep it current with this suite.
    run_test_get.step.dependOn(b.getInstallStep());

    const run_test_stats = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("tests/test_stats.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{.{ .name = "main", .module = main_module }},
        }),
    }));
    run_test_stats.step.dependOn(b.getInstallStep());

    const run_test_complement = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/complement.zig"),
            .target = target,
            .optimize = optimize,
        }),
    }));

    const run_test_bed_parser = b.addRunArtifact(b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/bed_parser.zig"),
            .target = target,
            .optimize = optimize,
        }),
    }));

    // Messy correctness FASTAs live in gitignored cache; materialize before tests that read them.
    const gen_messy_fixtures = b.addSystemCommand(&.{
        "python3",
        "bench/shared/generate_messy.py",
        "--fixtures",
    });

    // GET and stats subprocess tests load sidecars from tests/data. Generate them
    // with the selected build mode so `zig build test` works on a clean checkout.
    const test_fasta_paths = [_][]const u8{
        "tests/data/simple.fasta",
        "tests/data/proteome.fasta",
        "tests/data/single.fasta",
        "tests/data/edge_cases.fasta",
        "tests/data/mixed_widths.fasta",
    };
    for (test_fasta_paths) |fasta_path| {
        const gen_test_index = b.addRunArtifact(exe);
        gen_test_index.addArgs(&.{ "index", fasta_path });
        run_test_get.step.dependOn(&gen_test_index.step);
        run_test_stats.step.dependOn(&gen_test_index.step);
    }

    const test_validate_module = b.createModule(.{
        .root_source_file = b.path("tests/test_validate.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "main", .module = main_module }},
    });
    // Required for this suite on current Zig/OS (same as test_index).
    test_validate_module.link_libc = true;
    const run_test_validate = b.addRunArtifact(b.addTest(.{ .root_module = test_validate_module }));

    run_test_index.step.dependOn(&gen_messy_fixtures.step);
    run_test_validate.step.dependOn(&gen_messy_fixtures.step);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_test_index.step);
    test_step.dependOn(&run_test_get.step);
    test_step.dependOn(&run_test_stats.step);
    test_step.dependOn(&run_test_complement.step);
    test_step.dependOn(&run_test_bed_parser.step);
    test_step.dependOn(&run_test_validate.step);
}
