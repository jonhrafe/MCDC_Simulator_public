# Output and logging options

This page describes the options that control what the simulator writes besides the signal itself,
the files each one produces, their layout and their units.

Every option here is checked by the test `tests/test_log_outputs.py` (ctest `log_outputs`) against a
known answer, using an SI `.conf`. The checks are listed next to each option.

## Units of the output files

The `.conf` is read in SI units (m, s, m²/s, T/m), but **these output files are written in the
engine's internal units**: lengths in **mm**, times in **ms**, squared lengths in **mm²**, volumes in
**mm³**. This is the same for every output on this page; the signal files are unit-free (sums of
cos(phase) over walkers).

## Signal files

| option | what it does |
|---|---|
| `write_txt 1` | writes the signal as text: `<prefix>_DWI.txt`, one value per scheme row. Text values have **6 significant digits**. |
| `write_bin 1` | writes the same values as binary: `<prefix>_DWI.bfloat`, little-endian float32, one value per scheme row (full float32 precision). Every other output on this page that has a `.txt` file also gets a `.bfloat` twin. |
| `separate_signals` | also writes the signal split by compartment, `<prefix>_DWI_intra.txt` and `<prefix>_DWI_extra.txt` (and `.bfloat`). It is a flag: its presence turns it on, it takes no value. |

The signal value of a row is the sum over walkers of cos(phase), i.e. `S/S0 = value / N` at b = 0
with no T2 weighting.

*Checked:* each `.bfloat` equals its `.txt` file; free diffusion gives `_DWI_intra = 0`,
`_DWI_extra = _DWI`, and `_DWI / N = exp(-bD)` within Monte Carlo noise.

## Trajectories: `write_traj_file`

`write_traj_file 1` writes every walker's positions, one file per process `p`:

| file | content |
|---|---|
| `<prefix>_<p>.traj` (with `write_bin 1`) | float32 `x y z` per position, walker after walker: all positions of walker 0, then walker 1, ... |
| `<prefix>_<p>.bhdr` | float32 header: `duration [ms]`, `number of walkers in this file`, `number of positions` (see below) |
| `<prefix>_<p>.traj.txt`, `<prefix>_<p>.hdr.txt` (with `write_txt 1`) | the same, as text, one number per line |

- Positions are in **mm** and are the walker's real position (not wrapped back into the voxel by
  periodic boundaries), so displacements can be computed directly.
- Without a `<positions>` list (below) the file holds **T+1 positions per walker** (the start and the
  T steps) while the header's third value is **T**.
- Files are large (N × positions × 12 bytes).

*Checked:* `.traj` equals `.traj.txt`; the walkers of all processes add up to N; each step of a full
trajectory has length exactly `sqrt(6 D dt)`; the per-axis displacement variance is `2 D t`.

## `<log>` block

```
<log>
<positions>
0 100 t
</positions>
<propagator>
directions dirs.txt
100 t
</propagator>
</log>
```

Times are **step indices** (0 to T); `t` means the last step, T.

### `<positions>`

Restricts the trajectory file to the listed steps (it has an effect only with `write_traj_file 1`).
The header's third value is then the number of listed steps, and the file holds exactly that many
positions per walker.

*Checked:* only the listed steps are written, step 0 is the start position inside the voxel, and the
displacement variance at each listed step is `2 D t`.

### `<propagator>`

Despite its name, this records the **mean squared displacement** along given directions, not a full
propagator histogram.

- `directions <file>`: a plain text file with one direction `x y z` per line and **no header line**
  (a scheme file cannot be used: its first line stops the reader and no direction is loaded).
  Directions are normalised.
- the listed step indices: the times at which the displacement from the start is evaluated.
- Output `<prefix>_propagator.txt`: one row per time (sorted, duplicates removed), one column per
  direction; each value is the mean over all walkers of (displacement · direction)², in **mm²**.
  Binary twin: `<prefix>_propagator.txt.bfloat` (float32, row after row).
- It does not need `write_traj_file`.

*Checked:* every value equals the mean squared displacement recomputed from the trajectory file
(to float precision), and equals `2 D t` for free diffusion.

## `log_phase_shift`

A flag. Writes the distribution of the walkers' accumulated phase at the echo for every scheme row:
`<prefix>_phase_shift.txt` (and `.bfloat`), one row per scheme row, **3600 columns**: a histogram of
|phase| over [0, 2π) with bins of 0.1° (π/1800 rad). Each row holds N counts. Because cos is even,
the mean of cos(bin centre) over the histogram reproduces the signal `_DWI / N` to within the bin
width.

*Checked:* each row sums to N, and the histogram's mean cos(phase) equals `_DWI / N` (difference
< 2·10⁻⁵).

## `log_opp`

A flag. Also writes the signal of each process separately: `<prefix>_<p>_DWI.txt`, the imaginary
part `<prefix>_<p>_DWI_img.txt` and, with `log_phase_shift`, `<prefix>_<p>_phase_shift.txt`.

*Checked:* the per-process `_DWI.txt` files add up to the merged `_DWI.txt`.

## `compute_volume`

A flag. Adds a volume estimate to `<prefix>_simulation_info.txt`:

- `Estimated ICVF from sampling`: the fraction of walkers that end the walk inside an obstacle;
- `Total sampling volume (mm^3)`: the volume of `<spawning_area>`, or of the voxel if there is none;
- `Estimated Intra-axonal volume from sampling (mm^3)`: ICVF × total sampling volume.

The estimate is meaningful only with a **uniform** initialisation (no `ini_walkers_pos`): with
`ini_walkers_pos intra` every walker starts inside and the ICVF reads ~1.

*Checked:* a 4 µm sphere in a 10 µm voxel gives ICVF = (4/3 π R³)/V = 0.268 within binomial noise;
the total volume is the voxel volume; volume = ICVF × total.

## Example

A complete, runnable example with every option on this page is
[`docs/conf_file_examples/outputOptions.conf`](../docs/conf_file_examples/outputOptions.conf), with
the directions file
[`docs/conf_file_examples/propagator_directions.txt`](../docs/conf_file_examples/propagator_directions.txt).
Run it from the repository root:

```
./MC-DC_Simulator --conf docs/conf_file_examples/outputOptions.conf
```

It writes to `instructions/demos/output/output_options_*`. A quick check of the result: the two rows
of `output_options_propagator.txt` (steps 500 and 1000, i.e. 25 ms and 50 ms) are close to
2 D t = 3.0e-5 and 6.0e-5 mm² in every direction.

Reading the outputs in Python:

```python
import numpy as np
pre = "instructions/demos/output/output_options"
dwi   = np.fromfile(pre + "_DWI.bfloat", dtype="<f4")                  # one value per scheme row
hdr   = np.fromfile(pre + "_0.bhdr", dtype="<f4")                      # duration [ms], walkers, positions
traj  = np.fromfile(pre + "_0.traj", dtype="<f4").reshape(int(hdr[1]), int(hdr[2]), 3)   # [mm]
# without a <positions> list the file holds hdr[2] + 1 positions per walker (start + T steps)
prop  = np.loadtxt(pre + "_propagator.txt", ndmin=2)                   # rows = times, cols = directions [mm^2]
phase = np.loadtxt(pre + "_phase_shift.txt")                           # rows = scheme rows, 3600 bins of 0.1 deg
```

## Common mistakes

- **Giving a flag a value.** `separate_signals`, `log_phase_shift`, `log_opp`, `compute_volume`,
  `verbatim` and `debug` are switched on by writing the word alone. A value after them is read as
  the next parameter: `separate_signals 1` prints `Parameter: 1 Unknown`, and `separate_signals 0`
  still turns the option on.
- **Using a scheme file as the propagator `directions` file.** The directions file must contain only
  `x y z` lines. A scheme file starts with a header line, the reader stops there, no direction is
  loaded, and `_propagator.txt` comes out empty.
- **Giving times in seconds.** The times in `<positions>` and `<propagator>` are step indices
  (0 to T, `t` = T), not seconds.
- **Reading a full `.traj` with the header's count.** Without `<positions>`, the header's third value
  is T but the file holds T+1 positions per walker.
- **Expecting SI units in the outputs.** Positions are in mm, the trajectory duration in ms, the
  propagator in mm², volumes in mm³, even when the `.conf` is in SI units.
- **Expecting a trajectory without `write_traj_file 1`.** `<positions>` only selects which steps go
  into the trajectory file; on its own it writes nothing.
- **Reading `compute_volume` with a non-uniform start.** With `ini_walkers_pos intra` (or `extra`)
  the ICVF estimate is ~1 (or ~0) by construction.
- **Spheres list with fewer than 6 columns.** Each row must be `x y z radius permeability T2`; a
  shorter row is rejected by the format check.
