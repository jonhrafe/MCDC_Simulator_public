# Automatic initialization of spins in the internal or external compartment

In this example we show several options for the automatic initialization of spins in two geometries (obstacles):
- [Gamma distributed cylinders](#gamma-distributed-cylinders)
- [PLY models](#ply-meshes)

Before starting, build the latest version: [build instructions](compilation.md).

> **Units.** All `.conf` and scheme files below are in **SI units** (metres, seconds, Tesla) by
> default; values are scaled internally to mm/ms. The legacy `scale_from_stu` flag is deprecated
> (use `use_mm_ms 1` to declare a file already in internal units). See the
> [Units section of the README](../README.md#units).

# Gamma distributed cylinders

__The setup:__ a toy substrate of cylinders with radii drawn from a Gamma distribution. The `.conf`
file is in [`docs/conf_file_examples/`](../docs/conf_file_examples/).

```
N 1000
T 1000
duration 0.050
diffusivity 0.6e-9

exp_prefix instructions/demos/output/cylinder_gamma_packing_test

scheme_file docs/scheme_files/PGSE_sample_scheme.scheme

write_txt 1
write_bin 0
write_traj_file 1

<obstacle>
<cylinder_gamma_packing>
alpha 1.5
beta 0.5e-6
icvf 0.70
num_cylinders 1000
</cylinder_gamma_packing>
</obstacle>

ini_walkers_pos extra
num_process 1
<END>
```

### Basic parameters

- `N 1000` — number of spins.
- `T 1000` — number of time steps.
- `duration 0.050` — total diffusion time, in **seconds** (≥ the scheme's TE of 45 ms).
- `diffusivity 0.6e-9` — diffusion coefficient of the medium, in **m²/s** (near ex-vivo diffusion).
- `exp_prefix [string]` — output path and prefix.
- `scheme_file [string]` — the PGSE scheme file provided in the repository.
- `write_txt [0,1]` — text output (reduced precision; **warning**, large files that scale with N·T).
- `write_bin [0,1]` — binary output (full float32 precision, recommended).
- `write_traj_file [0,1]` — write per-particle trajectories (warning: very large files).
- `num_process [int]` — number of CPU threads. (With a trajectory file, the `.traj` is split per thread.)
- `<END>` — the `.conf` must end here; anything after is ignored.

### Defining the obstacle

An obstacle is declared between `<obstacle> … </obstacle>`. Inside, a gamma cylinder packing is
defined with `<cylinder_gamma_packing> … </cylinder_gamma_packing>` and four parameters:

- `alpha 1.5` — **shape** parameter of the Gamma distribution (dimensionless). (`shape` is an alias.)
- `beta 0.5e-6` — **scale** parameter, a length in **metres** (SI). (`scale` is an alias.)
- `icvf 0.70` — target intra-cylinder volume fraction (fraction of the voxel filled by cylinders).
- `num_cylinders 1000` — number of cylinders to pack; the voxel size is adjusted automatically.

With `alpha = 1.5` and `beta = 0.5e-6 m`, the Gamma distribution has a mean **radius** of
`alpha·beta = 0.75 µm` (mean diameter ≈ 1.5 µm). The packer aims for the target ICVF as closely as
possible; ICVFs above ~0.8 are nearly unfeasible.

### Spin initialization

To seed spins in one compartment only, use `ini_walkers_pos`:
`ini_walkers_pos extra` initializes spins **outside** the cylinders (extra-axonal space).

## Running and outputs

```bash
./MC-DC_Simulator docs/conf_file_examples/gammaDistributedCylinders.conf
```

Outputs are stored under `instructions/demos/output/`:

- `cylinder_gamma_packing_test_DWI.txt` — real part of the (un-normalized) DW signal.
- `cylinder_gamma_packing_test_DWI_img.txt` — imaginary part (≈ 0 here).
- `cylinder_gamma_packing_test_gamma_distributed_cylinder_list.txt` — positions and radii of the packed cylinders.
- `cylinder_gamma_packing_test_simulation_info.txt` — simulation info and any warnings.
- `cylinder_gamma_packing_test_0.traj.txt` (+ `.hdr`) — particle trajectories (for visualization).

### Visualizing the trajectories

The trajectory file holds all particle positions over time: `N · (T+1) · 3` values, ordered as
(x, y, z) for each step `0…T`, for each particle `1…N`. Visualizing them gives:

![simulation](https://user-images.githubusercontent.com/4105920/88835884-31fb9800-d1d6-11ea-8dcb-5210ae50793e.gif)

# PLY Meshes

The same automatic initialization can restrict spins to the **inside** of closed triangulated meshes
of arbitrary shape (e.g. the intra-axonal space).

![simulation](https://user-images.githubusercontent.com/4105920/88836514-18a71b80-d1d7-11ea-9e2a-a43df479a889.gif)

__The setup:__ 100 steps over 0.05 s with a low diffusivity of 0.6e-10 m²/s. The example mesh
(`instructions/meshes/PorusMedia.ply`) is stored in **micrometres**.

```
N 1000
T 100
duration 0.050
diffusivity 0.6e-10
exp_prefix instructions/demos/output/mesh_initialization_test

scheme_file docs/scheme_files/PGSE_sample_scheme.scheme

write_txt 1
write_bin 0
write_traj_file 1

<obstacle>
ply instructions/meshes/PorusMedia.ply
ply_scale 1e-6
</obstacle>

<voxels>
-1.0e-5 -1.0e-5 -5.0e-5
 1.0e-5  1.0e-5  5.0e-5
</voxels>

<spawning_area>
-5.0e-6 -5.0e-6 -5.0e-6
 5.0e-6  5.0e-6  5.0e-6
</spawning_area>

ini_walkers_pos intra

num_process 1

<END>
```

```bash
./MC-DC_Simulator docs/conf_file_examples/meshIntraInitialization.conf
```

### New parameters

- `ply [string]` — a closed, fully **triangulated** PLY mesh (vertices and faces only). The vertex
  units are whatever the file uses; `ply_scale` maps them to metres.
- `ply_scale [float]` — scale factor applied to the mesh vertices, in **metres per file unit** (SI).
  PorusMedia.ply is stored in micrometres, so `ply_scale 1e-6` converts micrometre vertices to
  metres. (Under the old mm-based convention this value was `0.001`; with SI units it is `1e-6`.)
- `<voxels> … </voxels>` — the voxel limits (where the signal is synthesized), six numbers in
  **metres** (SI): the lower corner then the upper corner.
- `<spawning_area> … </spawning_area>` — an optional custom region (in **metres**) where spins are
  seeded uniformly; useful when targeting a specific compartment. If omitted, spins fill the voxel.
- `ini_walkers_pos intra` — seed spins **only inside** the mesh: spins are placed uniformly in the
  spawning area, then any that land outside the mesh are discarded.

### Visualizing the outputs

As before, visualize the trajectory file to check the custom initialization:

![mesh](https://user-images.githubusercontent.com/4105920/88842439-bef71f00-d1df-11ea-9616-ff607cc2b6fe.gif)
