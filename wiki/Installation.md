# Installation

If you just want to use z-fasta, download a release archive. You get one static executable with no runtime dependencies, plus the MIT license.

## Download a release

Open the [z-fasta Releases page](https://github.com/eneskemalergin/z-fasta/releases) and choose the archive that matches your operating system and CPU:

- `z-fasta_0.3.1_linux_x86_64.tar.gz`
- `z-fasta_0.3.1_linux_arm64.tar.gz`
- `z-fasta_0.3.1_macos_x86_64.tar.gz`
- `z-fasta_0.3.1_macos_arm64.tar.gz`
- `z-fasta_0.3.1_windows_x86_64.zip`
- `z-fasta_0.3.1_windows_arm64.zip`

Extract the archive and place `z-fasta` or `z-fasta.exe` in a directory on your `PATH`.

### Linux and macOS

```bash
tar -xzf z-fasta_0.3.1_linux_x86_64.tar.gz
chmod +x z-fasta
./z-fasta --version
```

To call it from anywhere, move it to a user-owned binary directory:

```bash
mkdir -p "$HOME/.local/bin"
mv z-fasta "$HOME/.local/bin/z-fasta"
```

Make sure `$HOME/.local/bin` is on your `PATH`.

### Windows PowerShell

```powershell
Expand-Archive .\z-fasta_0.3.1_windows_x86_64.zip -DestinationPath .\z-fasta
.\z-fasta\z-fasta.exe --version
```

Add the extracted directory to your user `PATH` to call `z-fasta.exe` without its full path.

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

The executable is written to `zig-out/bin/z-fasta` or `zig-out/bin/z-fasta.exe`. On Windows PowerShell, verify the source build with:

```powershell
.\zig-out\bin\z-fasta.exe --version
```

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
z-fasta 0.3.1
```

Next: [Getting started](Getting-Started).
