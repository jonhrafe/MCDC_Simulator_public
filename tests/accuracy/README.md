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

## Experiment 2 — `two_meshes_impermeable.conf`

- **Substrate:** TWO impermeable meshes in one simulation (`meshes/unitMesh.ply`
  + a displaced deformed sphere `meshes/Mesh_O200.ply`), loaded via a simple
  `ply_file_list` (`two_meshes.list`, each line `<ply_file> <scale>`); permeability,
  T2 and d_intra come from the global params.
- **Walkers:** `N=2000`, `T=1000`, **initialised intra** (inside the meshes) in a
  larger voxel — exercises multi-mesh restricted diffusion. With impermeable
  membranes the **extra signal must stay ~0** (a built-in leak detector).
- **Impermeable:** no global permeability; a finite `T2` (0.080 s) set globally.
  `seed 12345`, `num_process 2`.

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
