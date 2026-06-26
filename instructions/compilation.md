# Building MC/DC

Clone the repository:

```bash
git clone https://github.com/jonhrafe/MCDC_Simulator_public.git
cd MCDC_Simulator_public
```

The only external dependency (the Eigen linear-algebra library) is bundled in `src/Eigen`, so no
package installation is required beyond a C++ compiler.

## Requirements

- A C++ compiler with **C++17** support (g++ ≥ 7, clang ≥ 5, or equivalent).
- **CMake ≥ 3.16** (recommended build path).

Supported platforms: Linux and macOS. On Windows, use the
[Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) and follow
the Linux instructions.

## Recommended: build with CMake

From the repository root:

```bash
cmake -B build
cmake --build build -j
```

This produces two self-contained executables in the repository root:

- **`MC-DC_Simulator`** — the simulator (takes a `.conf` file).
- **`dataSynth`** — a standalone signal-synthesis helper.

CMake defaults to an optimized `Release` build (`-O3`, matching `compile.sh`).

### Run the test suite (optional)

```bash
ctest --test-dir build
```

The CTest suite validates the physics (free/restricted diffusion, multi-compartment T2/D,
permeability, gamma packings, PGSE and waveform sequences, …) against golden references.

## Alternative: quick build

A single shell script builds the binaries with `g++` (no CMake):

```bash
bash compile.sh
```

It compiles `MC-DC_Simulator` and `dataSynth` into the repository root using
`-O3 -std=c++17 -march=native`. Note that `-march=native` produces a non-portable binary (best for
local runs, not for distributing to other machines).

## Getting started

To test your build, run the first simulation in the
[Getting started page](GettingStarted.md):

```bash
./MC-DC_Simulator docs/conf_file_examples/freeDiffusion.conf
```
