# MC/DC test suite

The suite runs the simulator end to end and inspects the DW-MRI signal it produces. The runner is
CTest (see the repo `CMakeLists.txt`); every test is a standalone Python script that needs only
`numpy`. Shared plumbing lives in [`mcdc_testkit/`](mcdc_testkit/).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
ctest --test-dir build --output-on-failure            # the whole suite (~23 min with -j 3 on 16 cores)
ctest --test-dir build -L validation_analytic         # one family (labels below)
python3 tests/validation/test_free_exp.py --n 20000   # one script, with its own options
```

CTest passes the built binary to each script as `--bin`; a script run by hand uses `$MCDC_BIN` or
`./MC-DC_Simulator` at the repo root.

## Three families

| family | question | tests (CTest names) | reference |
|---|---|---|---|
| **Regression** | does the code still produce the signal it produced yesterday? | `accuracy_golden*` (sphere mesh, two meshes, hex packings, gamma packings), `rng_reproducibility`, `population_history`, `obstacle_properties`, `accuracy_permeability`, `permeability_list`, `permeability_equilibrium*`, `permeability_coarse_step`, `list_units`, `waveform_protocol` | same-build golden signals in [`accuracy/golden/`](accuracy/golden/) (bit-exact, fixed `seed` and `num_process`), or an invariant (crossing counters, volume fractions) |
| **Option checks** | does every `.conf` option do what the documentation says? | `format_checks`, `log_outputs`, `init_options`, `substrate_options` | a known answer per option; the documented behaviour is in [`../instructions/`](../instructions/) |
| **Analytic validation** | is the physics correct, not merely unchanged? | `validation_*` ([`validation/`](validation/)) | closed-form signals: free exp(-bD), van Gelderen cylinder and Neuman sphere (Gaussian phase approximation), Kärger exchange, exp(-TE/T2), volume-weighted GPA over a packing's radii; plus conservation and monotonicity checks where no closed form exists |
| **Cross-version** | does the current build still match the historical releases? | `crossversion_free_diffusion_v1.5`, `_v2.0` ([`crossversion/`](crossversion/)) | the v1.5 and v2.0 releases, compiled from this repository's own history on first use (`tools/build_old_version.py`, cached under the gitignored `crossversion/old_binaries/`) and run on the same experiment |

Labels: `validation_analytic` (the `validation_*` tests), `crossversion`, `validation` (option
checks), `accuracy`, `permeability`, `reproducibility`.

## Reading a result

- **Golden tests** compare bit-exactly (`|run − gold| ≤ atol + rtol·|gold|`, tight). A difference
  means the physics or the RNG stream changed. After an intentional, reviewed change regenerate
  the golden with `python3 tests/test_accuracy_golden.py --conf <conf> --update` and say so in the
  commit ([`accuracy/README.md`](accuracy/README.md)).
- **Analytic tests** compare within a statistical band at a converged N (the band is printed with
  the result). Each script's docstring and [`validation/README.md`](validation/README.md) give the
  oracle and its validity limits (the GPA curves hold for impermeable obstacles at small b).
- **Cross-version tests** compare statistically, never bit-exactly: the releases use different
  RNGs. Per scheme row the Monte Carlo standard error of the difference is known, so the test
  reports the largest |z| over the rows and a z per b-shell ([`crossversion/README.md`](crossversion/README.md)).

A run is reproducible only for a fixed positive `seed` **and** a fixed `num_process` (the
per-thread RNG streams depend on the thread count).

## Layout

```
tests/
├── test_*.py              regression and option checks (one script per CTest entry)
├── test_crossversion.py   current build vs an old release on one case
├── validation/            analytic validation scripts + README (oracles, validity)
├── crossversion/          old-release registry (binaries.toml), per-version .conf files, comparator
├── mcdc_testkit/          io (conf rewriting, signal loading), run (error / skip detection),
│                          scheme, metrics, oracles, provenance, paths
├── tools/                 generators: build_old_version.py, gen_analytic.py (needs dmipy),
│                          gen_reference.py, make_meshes.py, ...
└── accuracy/              fixtures: .conf files, lists, meshes, schemes, golden/ and reference_xv/
```

`test_gpu_*.py` are the scripts of the GPU line; they are not registered here. A `.conf` with
`gpu 1` stops this build with "built WITHOUT CUDA"; `gpu 0` is accepted, so the GPU line's `.conf`
files run unchanged on the CPU.

Outputs go to `tests/accuracy/output/` and `tests/output/` (both gitignored).
