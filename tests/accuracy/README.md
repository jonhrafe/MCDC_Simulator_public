# Accuracy golden-master regression test

A single fixed experiment whose diffusion signal is frozen as a baseline, so any
change elsewhere in the code can be checked for *"does it still produce the same
physics?"*. This is also the reference a future GPU / re-implementation is
validated against.

## The experiment — `sphere_impermeable.conf`

- **Substrate:** a simple sphere mesh (`debug/unitMesh.ply`) inside a voxel.
- **Walkers:** `N=2000`, `T=1000`, seeded **uniformly** in the voxel (no
  `ini_walkers_pos`), so the ensemble is a natural **intra + extra** mix.
- **Impermeable:** no `permeability` key. This is the deliberate invariant — the
  permeability work (and anything else not meant to change the signal) must
  leave this run untouched.
- **Deterministic:** fixed `seed 12345` and `num_process 2`.

Captured outputs (the actual scientific signal, 270 PGSE measurements each):
`*_DWI.txt` (real), `*_DWI_intra.txt`, `*_DWI_extra.txt`, stored in `golden/`.

## Running

```bash
# check the current build against the committed golden (also: ctest -R accuracy_golden)
python3 tests/test_accuracy_golden.py

# loosen tolerance, e.g. when comparing a GPU port or a different compiler
python3 tests/test_accuracy_golden.py --rtol 1e-4

# regenerate the baseline after an INTENTIONAL, reviewed behaviour change
python3 tests/test_accuracy_golden.py --update
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
