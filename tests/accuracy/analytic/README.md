# Frozen analytic reference curves

Deterministic closed-form / dmipy signals, frozen as committed `.npy` arrays so
the test-time suite needs only **numpy** (no dmipy at run time).

Each curve is produced once by [`tests/tools/gen_analytic.py`](../../tools/gen_analytic.py)
(which needs `pip install .[analytic]` → dmipy) and ships with a
`<name>.provenance.json` recording the dmipy version, model, physical parameters
(radius, D, scheme sha256) and the date — so the curve is auditable and
regeneratable.

Files:

| file | oracle | validity |
|------|--------|----------|
| `cylinder_r5um_d0p6.npy` | dmipy `C4CylinderGaussianPhaseApproximation` (van Gelderen, infinite cylinder) | **impermeable, small b (≲ 5000 s/mm²)**; perpendicular restricted diffusion only |

Loaded via `mcdc_testkit.oracles.load_analytic(name)`. A test that needs a curve
not yet frozen **self-skips** with a pointer to the generator (a fresh checkout,
or an env without dmipy, is not a failure).
