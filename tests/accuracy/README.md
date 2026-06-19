# Accuracy golden-master regression tests

Fixed experiments whose diffusion signals are frozen as baselines, so any change
elsewhere in the code can be checked for *"does it still produce the same
physics?"*. These are also the reference a future GPU / re-implementation is
validated against.

All fixtures are **self-contained and tracked** under `tests/accuracy/` (meshes
in `meshes/`, the PGSE scheme, the configs and list) — nothing depends on the
gitignored `debug/` directory, so the suite works on a clean checkout.

## Experiment 1 — `sphere_impermeable.conf`

- **Substrate:** a simple sphere mesh (`meshes/unitMesh.ply`) inside a voxel.
- **Walkers:** `N=2000`, `T=1000`, seeded **uniformly** in the voxel (no
  `ini_walkers_pos`), so the ensemble is a natural **intra + extra** mix.
- **Impermeable:** no `permeability` key. The deliberate invariant — the
  permeability work (and anything else not meant to change the signal) must
  leave this run untouched.
- **Deterministic:** fixed `seed 12345` and `num_process 2`.

## Experiment 2 — `two_meshes_impermeable.conf` (stress)

- **Substrate:** TWO impermeable meshes in one simulation (`meshes/unitMesh.ply`
  + a displaced deformed sphere `meshes/Mesh_O200.ply`), loaded via
  `ply_extended_file_list` (`two_meshes.list`) with **distinct `d_intra` and `T2`
  per mesh**, plus **distinct extra-cellular `d_extra` / `t2_extra`**.
- **Walkers:** `N=2000`, `T=1000`, seeded **uniformly** in a voxel that tightly
  bounds the two meshes (combined PLY bbox + ~0.1 µm), so the intra compartment is
  heavily sampled (~28%: ~440 in sphere1, ~120 in sphere2) while extra is still
  populated — this stresses that each particle picks up the Di/T2 of the
  compartment it starts in.
- **Impermeable:** permeability `0` for both meshes. `seed 12345`, `num_process 2`.

To validate the per-compartment assignment directly, add the **`debug`** flag
(presence-only — just the word `debug` on its own line; there is no `debug 1`,
and `debug 0` would NOT turn it off — remove the line to disable). The simulator
then writes a per-process `*_debug_trace.txt` logging, for every walker and step,
`walker step x y z location compartment Di T2` (position in mm, Di in mm²/ms, T2
in ms; the model has no T1). Each row's `Di`/`T2` should match the compartment in
the `compartment` column (e.g. `ply0`, `ply1`, `extra`). The file is large, so it
is only written when `debug` is present.

### PLY list options

- `ply_file_list <file>` — each line `<ply_file> <scale>`; permeability, T2 and
  d_intra inherited from the global params.
- `ply_extended_file_list <file>` — each line `<ply_file> <scale> <d_intra> <t2>
  <permeability>`, **per mesh**, in standard units (d_intra m²/s, t2 s,
  permeability m/s; scaled internally like the rest of the .conf; `scale` is the
  geometry file→mm factor and is never unit-scaled).
- `ply_file_list_scale_permeability <file>` — legacy `<ply_file> <scale>
  <permeability>` (T2/d_intra global). Kept for back-compatibility.

The geometry `scale` is mandatory in every form (it has no global default).

Captured outputs (the actual scientific signal, 270 PGSE measurements each):
`*_DWI.txt` (real), `*_DWI_intra.txt`, `*_DWI_extra.txt`, stored in `golden/`.

## Running

```bash
# all accuracy + reproducibility tests
ctest --test-dir build

# experiment 1 (default conf) against the committed golden
python3 tests/test_accuracy_golden.py

# experiment 2 (two meshes)
python3 tests/test_accuracy_golden.py --conf tests/accuracy/two_meshes_impermeable.conf

# loosen tolerance, e.g. when comparing a GPU port or a different compiler
python3 tests/test_accuracy_golden.py --rtol 1e-4

# regenerate a baseline after an INTENTIONAL, reviewed behaviour change
python3 tests/test_accuracy_golden.py [--conf <conf>] --update
```

Comparison is `|run - golden| <= atol + rtol*|golden|` per value (defaults
`rtol=1e-9`, `atol=1e-12` — effectively bit-identical on the same build). The
test reports the worst absolute and relative deviation per signal.

## When the baseline legitimately changes

`num_process 2` makes the run reproducible only for that process count (see the
seed warning in `SimErrno::checkSimulationParameters`). Any *intentional* physics
change that alters this impermeable signal — regenerate with `--update` and call
it out explicitly in the commit. A change you did **not** expect to move the
signal failing here is the test doing its job.
