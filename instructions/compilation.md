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
- *(Optional)* a **CUDA toolkit** (`nvcc`) to build the GPU backend. It is auto-detected; without it
  the full simulator still builds and runs on the CPU.

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

CMake defaults to an optimized `Release` build. If a CUDA toolkit is found, the GPU backend is
compiled and the `gpu` configuration flag becomes available; otherwise the build is CPU-only. You
will see one of these lines during configuration:

```
-- CUDA found (...); GPU path ENABLED, arch 86.
-- CUDA not found; GPU path DISABLED (CPU-only build).
```

To target a specific GPU architecture (default is `86`, Ampere/RTX 30xx):

```bash
cmake -B build -DCMAKE_CUDA_ARCHITECTURES=75   # e.g. Turing
```

### Run the test suite (optional)

```bash
ctest --test-dir build
```

The CTest suite validates the physics (free/restricted diffusion, multi-compartment T2/D,
permeability, PGSE and waveform sequences, …) against golden references. GPU tests automatically
skip (pass) on machines without CUDA.

## Alternative: quick CPU-only build

A single shell script builds a CPU-only binary with `g++` (no CMake, no CUDA):

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
