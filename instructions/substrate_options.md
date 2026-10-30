# Substrate and run-control options

This page describes the obstacle list files, the PLY list files, the voxel boundary conditions and
the options that control a run (`max_sim_time`, `verbatim`, `use_mm_ms`). Every option here is checked
by `tests/test_substrate_options.py` (ctest `substrate_options`) against a known answer; the checks
are listed next to each option.

## Units: SI or `use_mm_ms 1`

By default every value in the `.conf` and in the files it points to is in SI units (metres, seconds,
Tesla). `use_mm_ms 1` declares that the `.conf`, the scheme file and the list files are already in
the internal units (mm, ms, T/mm); nothing is then converted.

| quantity | SI (default) | with `use_mm_ms 1` |
|---|---|---|
| lengths (`<voxels>`, `<spawning_area>`, `<delta>`, positions file) | m | mm |
| scale factor (first line of a sphere/cylinder list, `ply_scale`, PLY list scale) | m per file unit | mm per file unit |
| `duration`, T2 (`.conf` and list files) | s | ms |
| diffusivity | m²/s | mm²/ms |
| permeability | m/s | mm/ms (the same number) |
| scheme file: G / Δ, δ, TE | T/m / s | T/mm / ms |

Output files are always in the internal units (mm, ms); see [output_options.md](output_options.md).

*Checked:* the same substrate written both ways gives the same signals (total, intra, extra) to within
3e-6 of the b = 0 signal, for a `spheres_list` (permeable, with a per-sphere T2), a `cylinders_list`
(per-cylinder T2), a `ply_file_list_scale_permeability` (permeable mesh), and for start positions
from an `ini_walkers_file`.

## Sphere and cylinder lists

Inside the `<obstacle>` block:

```
<obstacle>
spheres_list   path/to/spheres.list
cylinders_list path/to/cylinders.list
</obstacle>
```

The first line of the file is the scale factor alone (metres per file unit, e.g. `1e-6` for files
written in micrometres). Each following line is one obstacle:

| file | columns |
|---|---|
| spheres | `x y z r permeability T2` |
| cylinders | `x y z r` or `x y z r permeability T2` (the same number of columns on every line) |

Positions and radii are in file units (multiplied by the scale), permeability in m/s (0 =
impermeable), T2 in seconds (mm/ms files: ms). Every cylinder is infinite and along z, through the
point (x, y, z).

### `oriented_cylinders_list`

Reads the same file format as `cylinders_list`; the cylinders are along z. Versions up to 1.5 read a
7-column format (`x1 y1 z1 x2 y2 z2 r`, a cylinder through two points); this format is no longer
supported and is rejected at start-up with an explanation, so an old file is never misread.

*Checked:* a 6-column file gives the same signal through `oriented_cylinders_list` and
`cylinders_list`; a 7-column file is rejected.

## PLY meshes and PLY lists

There are four ways to give meshes, all inside the `<obstacle>` block:

| option | what it sets |
|---|---|
| `ply file.ply` + `ply_scale s` | one mesh; scale s (m per file unit); permeability from `obstacle_permeability` |
| `ply_file_list list.txt` | one mesh per line: `<ply> <scale>` |
| `ply_file_list_scale_permeability list.txt` | one mesh per line: `<ply> <scale> <permeability>` |
| `ply_extended_file_list list.txt` | one mesh per line: `<ply> <scale> <d_intra> <T2> <permeability>` (see [tutorial_multi_compartment.md](tutorial_multi_compartment.md)) |

In a list file, lines starting with `#` are comments. Scale is in metres per file unit (mm per file
unit with `use_mm_ms 1`), permeability in m/s (0 = impermeable).

A global `obstacle_permeability` greater than 0 **overrides** every per-obstacle permeability, from a
sphere/cylinder list or a PLY list alike.

*Checked:* `ply + ply_scale`, `ply_file_list` and `ply_file_list_scale_permeability` with permeability
0 give identical signals for the same mesh; a per-mesh permeability of 3e-5 m/s gives the identical
signal as `ply_file_list` with `obstacle_permeability 3e-5`, and a different one from the
impermeable mesh.

Example: [`docs/conf_file_examples/plyListSymmetric.conf`](../docs/conf_file_examples/plyListSymmetric.conf)
with [`three_spheres_permeable.list`](../docs/conf_file_examples/three_spheres_permeable.list) (three
spheres, three permeabilities, symmetric boundaries).

## `boundary_conditions periodic | symmetric`

What happens when a walker reaches a face of the voxel:

- `periodic` (default): it leaves through the face and comes back through the opposite face, as in an
  infinite tiling of copies of the voxel.
- `symmetric` (also `mirror`): it is reflected back into the voxel, as in an infinite tiling of
  mirrored copies of the voxel.

In both cases the signal is computed from the walker's real (unfolded) position, so free diffusion
gives exp(-b D) whatever the voxel size.

*Checked:* free diffusion in a 5 µm box, where each walker crosses the faces many times: the signal
matches exp(-b D) within Monte Carlo noise for both (123 measurements, max |z| < 5, mean z within ±1).
This check found a bug, now fixed: with `symmetric`, the part of a step up to the face moved the real
position the wrong way after an odd number of reflections, which lowered the signal (up to 0.06 below
exp(-b D) in this test).

## `max_sim_time`

```
max_sim_time 3600
```

Stops each process after the given **wall-clock** time in seconds (not a simulated time, so it is not
converted by `use_mm_ms`). The walkers not yet simulated are dropped: the signal is the sum over the
walkers simulated so far, and `simulation_info.txt` reports their number ("Number of simulated
particles"). Divide by that number, not by N. The clock has 1-second resolution.

*Checked:* a run of 3 000 000 walkers with `max_sim_time 2` stops after about 1.5 s; the b = 0 signal
equals the reported number of simulated walkers.

## `verbatim`

A flag (no value). Prints more information while running and writes
`<prefix>_perm_counters.txt` (counts of hits and crossings per permeable obstacle, a diagnostic of the
permeability model) when there are permeable obstacles.

*Checked:* the file is written with `verbatim` and not without.

## `coll_sphere_dist`

No longer has an effect; it is accepted with a warning so that old `.conf` files still run.

## Common mistakes

- **A list scale in the old units.** The first line of a sphere/cylinder list and the PLY list scale
  are metres per file unit in an SI `.conf`: `1e-6` for a file in micrometres. `1e-3` (mm per unit)
  makes everything 1000 times too large.
- **Reusing a written `*_gamma_distributed_*_list.txt` as an input list in an SI `.conf`.** Its
  first line is `0.001` (the coordinates are in µm and the scale is mm per µm); change it to `1e-6`
  before using it as a `spheres_list` / `cylinders_list` in an SI `.conf`, or use it with
  `use_mm_ms 1` (then also give T2 in ms).
- **Setting `obstacle_permeability` and per-obstacle permeabilities.** The global value replaces
  all of them; leave it out (or 0) to keep the per-obstacle values.
- **T2 in a list file in ms with an SI `.conf`.** T2 in a list is in seconds like every other time.
- **A 7-column cylinder list (two points per cylinder).** Rejected: cylinders are along z.
- **Dividing a `max_sim_time` signal by N.** Divide by the number of simulated walkers reported in
  `simulation_info.txt`.
- **Writing a value after `verbatim`.** It is a flag; see [output_options.md](output_options.md).
