# Getting Started

## Basic usage

```bash
./MC-DC_Simulator configuration_file.conf
```

- `MC-DC_Simulator`: the application (see [building](compilation.md)).
- `configuration_file.conf`: a text file with **all** the simulation parameters, listed one per line.

## Test your installation

The folder `docs/conf_file_examples/` contains the commented examples used in the tutorials. From
the repository root, run:

```bash
./MC-DC_Simulator docs/conf_file_examples/freeDiffusion.conf
```

If everything is set up correctly, the output is written next to the configured `exp_prefix`
(`instructions/demos/output/` for this example).

## A first free-diffusion simulation

Below is the body of `docs/conf_file_examples/freeDiffusion.conf` (the file itself carries a
detailed units header as comments):

```
N 1000
T 1000
duration 0.100
diffusivity 0.6e-9

exp_prefix instructions/demos/output/free_diffusion_test

scheme_file docs/scheme_files/PGSE_sample_scheme.scheme

write_txt 1
write_bin 1
write_traj_file 1

<voxel>
0 0 0
0.5e-3 0.5e-3 0.5e-3
</voxel>

num_process 1

<END>
```

Parameters are read **until the `<END>` tag**. Each parameter is one line: the name, a space, then
the value. Any path must be either an absolute system path or relative to the directory the
simulator is launched from.

> **Units.** By default the `.conf` (and its scheme file) is in **SI units** — metres, seconds,
> Tesla — scaled internally to mm/ms. Add `use_mm_ms 1` to declare a file already in the internal
> units instead. The legacy `scale_from_stu 1` flag is still accepted but deprecated. See the
> [Units section of the README](../README.md#units).

### Parameters

- **`N`** [int]: number of spin particles to diffuse.
- **`T`** [int]: number of time steps over the experiment duration.
- **`duration`** [float]: total diffusion time in **seconds** (SI) — at least as long as the longest
  echo time (TE) in the scheme.
- **`diffusivity`** [float]: diffusion coefficient of the medium in **m²/s** (SI).
- **`exp_prefix`** [string]: output path and filename prefix for the experiment.
- **`scheme_file`** [path]: the acquisition protocol (PGSE or general waveform).
- **`write_txt`** [0,1]: write text output (reduced numerical precision).
- **`write_bin`** [0,1]: write binary output (full float32 precision, **recommended**).
- **`write_traj_file`** [0,1]: write the per-particle trajectories (**warning: very large files**),
  in text or binary depending on the flags above.
- **`<voxel> … </voxel>`**: the simulation voxel, in **metres** (SI). Three numbers for the minimum
  corner (x_min, y_min, z_min), then three for the maximum corner (x_max, y_max, z_max).
- **`num_process`** [int]: number of CPU threads to use.

### Output files

Outputs share the `exp_prefix` name with an appended suffix:

- real part of the DW-MRI signal (`_DWI`),
- imaginary part of the DW-MRI signal (`_DWI_img`),
- a simulation info file (`_simulation_info`),
- and, if enabled, the trajectory file.

### What this run does

This example diffuses 1,000 spins for 100 ms in 1,000 steps (one step every 0.1 ms). The signal is
computed with the PGSE scheme `PGSE_sample_scheme.scheme`. Walkers are seeded in an arbitrary cubic
voxel from (0, 0, 0) to (0.5, 0.5, 0.5) mm (written as `0.5e-3` m in SI), and the run uses a single
CPU thread.
