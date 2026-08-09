# Installation

If you just want to use z-fasta, download a release archive. You get one static executable with no runtime dependencies, plus the MIT license.

## Download a release

Open the [z-fasta Releases page](https://github.com/eneskemalergin/z-fasta/releases) and choose the archive that matches your operating system and CPU:

- `z-fasta_0.3.2_linux_x86_64.tar.gz`
- `z-fasta_0.3.2_linux_arm64.tar.gz`
- `z-fasta_0.3.2_macos_x86_64.tar.gz`
- `z-fasta_0.3.2_macos_arm64.tar.gz`

Extract the archive and place `z-fasta` in a directory on your `PATH`.

### Linux and macOS

```bash
tar -xzf z-fasta_0.3.2_linux_x86_64.tar.gz
chmod +x z-fasta
./z-fasta --version
```

To call it from anywhere, move it to a user-owned binary directory:

```bash
mkdir -p "$HOME/.local/bin"
mv z-fasta "$HOME/.local/bin/z-fasta"
```

Make sure `$HOME/.local/bin` is on your `PATH`.

### Windows through WSL

I would like to support Windows properly, but keeping the native builds reliable currently takes more time than I can justify. I would rather be honest about that than publish an executable I cannot support well. I am sorry for the extra step.

For now, install [Windows Subsystem for Linux](https://learn.microsoft.com/windows/wsl/install), open its Linux shell, and use the Linux archive matching the WSL architecture. Run the Linux installation and verification commands above inside WSL.

> [!IMPORTANT]
> Confirm that `z-fasta --version` prints the release you intended to install before building indexes used by a pipeline.

## Build from source

Build from source when you are developing z-fasta or need a commit that has not reached a release. Download Zig 0.16.0 from the [official Zig download page](https://ziglang.org/download/), add its executable to your `PATH`, and confirm the version before building.

```bash
zig version
git clone https://github.com/eneskemalergin/z-fasta.git
cd z-fasta
zig build -Doptimize=ReleaseFast
./zig-out/bin/z-fasta --version
```

The first command must print `0.16.0`. Other Zig versions are not part of the supported source-build contract.

The executable is written to `zig-out/bin/z-fasta`.

Run the test suite when you are building a development checkout:

```bash
zig build test --summary all
zig build test -Doptimize=ReleaseFast --summary all
```

## Verify the installation

```bash
z-fasta --help
z-fasta --version
```

Expected version output:

```text
z-fasta 0.3.2
```

Next: [Getting started](Getting-Started).
