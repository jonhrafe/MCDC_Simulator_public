# Walker initialisation options

This page describes where the walkers start and the options that control it. Every option here is
checked by `tests/test_init_options.py` (ctest `init_options`) against a known answer, using an SI
`.conf`; the checks are listed next to each option.

All positions in the `.conf` and in a positions file are in **metres** (SI), unless `use_mm_ms 1`
is set (then mm).

## Which option decides the start

Only one start mode applies. From highest to lowest priority:

| mode | how to request it |
|---|---|
| positions file | `ini_walkers_file <file>` |
| single point | `ini_walkers_pos delta` + `<delta>` block |
| inside obstacles | `ini_walkers_pos intra` |
| outside obstacles | `ini_walkers_pos extra` |
| regular grid | `<regular_sampling>` |
| uniform random (default) | nothing |

Combinations that would silently ignore one of the options are rejected at start-up (see
[Common mistakes](#common-mistakes)).

## Uniform start (default) and `<spawning_area>`

Without any option, walkers start uniformly at random in the voxel. A `<spawning_area>` restricts the
start to a box (in metres), which also applies to `ini_walkers_pos intra` / `extra` and to
`<regular_sampling>`:

```
<spawning_area>
-5e-6 -5e-6 -5e-6
 5e-6  5e-6  5e-6
</spawning_area>
```

The walkers still diffuse in the whole voxel; only the start is restricted.

*Checked:* with spheres in the voxel, every start (uniform, intra and extra) is inside the box.

## `ini_walkers_pos intra` / `extra`

Starts every walker inside an obstacle (`intra`) or outside all obstacles (`extra`), uniformly within
that compartment (and within the `<spawning_area>` if given). A random position is drawn and redrawn
until it falls in the requested compartment; the rejected draws are not counted as walkers. If no
position is found after 100000 draws (`intra`) or 10000 (`extra`) the run stops with an error: the
sampling area does not reach the requested compartment (or reaches very little of it; make the
`<spawning_area>` fit the obstacles).

- For spheres and cylinders the inside test is exact.
- For PLY meshes the inside test counts how many times a ray from the point crosses the mesh
  (odd = inside), so it works for concave meshes too. Positions very close to a wall are accepted,
  as a uniform start requires.

*Checked:* among spheres, every intra start is inside a sphere and every extra start outside all of
them; on a concave mesh (a torus), every intra / extra start is on the right side by an independent
ray-parity test, and a uniform start fills the mesh in proportion to its volume.

## `ini_walkers_pos delta` + `<delta>`

Every walker starts at the same point (metres):

```
ini_walkers_pos delta
<delta>
1e-6 -2e-6 0
</delta>
```

*Checked:* all walkers start at the point (2 processes); for free diffusion the displacement variance
after t is 2 D t.

Example: [`docs/conf_file_examples/deltaInitialization.conf`](../docs/conf_file_examples/deltaInitialization.conf).

## `ini_walkers_file`

Reads the start positions from a text file with one `x y z` per line, in **metres** (mm with
`use_mm_ms 1`). Walker number g, counted over all processes, starts at row `g mod K` (K = number of
rows), so with more walkers than rows the file is reused from the top. Every position must be inside
the voxel, otherwise the run stops with an error.

*Checked:* with 7 walkers, 3 rows and 2 processes, walker g starts at row g mod 3, converted from
metres; a position outside the voxel is rejected.

Example: [`docs/conf_file_examples/walkerPositionsFile.conf`](../docs/conf_file_examples/walkerPositionsFile.conf)
with [`walker_positions.txt`](../docs/conf_file_examples/walker_positions.txt).

A second example, [`walkerPositionsSpheres.conf`](../docs/conf_file_examples/walkerPositionsSpheres.conf)
with [`walker_positions_spheres.txt`](../docs/conf_file_examples/walker_positions_spheres.txt), puts
the walkers at the centres of three impermeable spheres with different T2 and at one point outside;
its comments give the exact b = 0 intra and extra signals to expect (checked by the test).

## `<regular_sampling>`

A flag (no value). Places the walkers on a regular grid instead of at random: the voxel (or the
`<spawning_area>`) is divided into n × n × n cells, n = ceil(cbrt(N)), and each walker sits at the
centre of one cell. If N is a cube (1000, 8000, ...) every cell holds exactly one walker; otherwise
the N walkers are spread evenly over the n³ cells. Walkers of different processes never share a cell.

*Checked:* N = 1000 on 2 processes gives exactly the 10³ cell centres, each once; N = 2000 spreads
over 13 layers per axis with balanced counts (153–154) and a centred mean.

Example: [`docs/conf_file_examples/regularSampling.conf`](../docs/conf_file_examples/regularSampling.conf).

## Checking the start positions yourself

Write the trajectory at step 0 only, and read it (positions are in mm, see
[output_options.md](output_options.md)):

```
write_bin 1
write_traj_file 1
<log>
<positions>
0
</positions>
</log>
```

```python
import numpy as np
pre = "instructions/demos/output/delta_initialization"
hdr = np.fromfile(pre + "_0.bhdr", "<f4")                       # duration [ms], walkers, positions
start = np.fromfile(pre + "_0.traj", "<f4").reshape(int(hdr[1]), int(hdr[2]), 3)[:, 0]   # [mm]
```

There is one `_<p>.traj` file per process.

## Common mistakes

- **`ini_walkers_pos delta` without a `<delta>` block.** Rejected: the walkers would otherwise start
  uniformly without notice.
- **`ini_walkers_pos intra` / `extra` together with `<delta>` or `ini_walkers_file`.** Rejected: the
  fixed positions would win and the intra/extra request would be ignored.
- **`<regular_sampling>` together with `ini_walkers_pos`.** Rejected, for the same reason.
- **More than one of `<delta>`, `ini_walkers_file` and `<regular_sampling>`.** Rejected: only one of
  them can set the start positions.
- **A positions file in mm with an SI `.conf`.** The file is in metres like every other quantity
  (`1e-6` = 1 µm); values in mm would put the walkers 1000 times too far out (rejected if outside the
  voxel) or in the wrong place.
- **Expecting the `<spawning_area>` to confine the walkers.** It only sets where they start; they then
  diffuse in the whole voxel.
- **A `<spawning_area>` (or voxel) that misses the obstacles with `ini_walkers_pos intra`, or lies
  entirely inside them with `extra`.** The run stops with "Max. number of tries to find an
  intra-celular / extra-celular compartment reached" (earlier release builds looped forever
  printing it).
- **Reading `compute_volume` with `intra` / `extra` starts.** The ICVF estimate is meaningful only
  with a uniform start (see [output_options.md](output_options.md)).
