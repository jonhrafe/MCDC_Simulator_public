# MC-DC Simulator — Phase 0 Development Plan (CPU Correctness Foundation)

**Branch:** `mcdc2_dev` (all Phase 0 work stays here; GPU work will be a separate
branch later). **Goal:** a *correct, reproducible, validated* CPU simulator
before any GPU work begins. GPU readiness is a deliberate non-goal of Phase 0,
but Phase 0 fixes (especially the RNG redesign) are what make a clean GPU port
possible later.

File paths below are relative to the simulator repo
(`.../Simulators/Claude/MCDC_Simulator_public/`), `src/` unless noted.

---

## Why Phase 0 comes before everything

The simulator's physics *model* is largely sound (proper bounce-resolution
collision engine, AABB-grid acceleration, Powles-style permeable membranes,
compartment T2, gamma substrates). But it is **not reproducible** and has **no
validation harness**. You cannot certify a permeability improvement or a future
GPU port without (a) deterministic runs and (b) analytic ground truth. The RNG
redesign that fixes reproducibility is also the single biggest enabler of a SIMT
GPU port — one fix unblocks all three research lines.

---

## The reproducibility problem (4 independent RNG sites)

All four must be routed through one seeded engine:

1. **Per-thread engines share the same seed** — `parallelmcsimulation.cpp:180`
   (seed copied, not offset). Parallel fixed-seed runs produce *duplicated*
   ensembles, not independent samples → statistically wrong.
2. **Walker initial placement** constructs a fresh `random_device`-seeded
   `mt19937` per call — `walker.cpp:197-199`, `dynamicsSimulation.cpp:503-505`,
   `:550-552`. Even single-threaded fixed-seed runs are nondeterministic.
3. **Gamma substrate generation** ignores the seed —
   `cylindergammadistribution.cpp:71-72,81` (and the sphere equivalent). Geometry
   differs every run.
4. **Membrane-crossing draw uses C `rand()/RAND_MAX`** — `cylinder.cpp:136`,
   `sphere.cpp:112`, `plyobstacle.cpp:383`. Unseeded, not thread-safe, ~32767-level
   quantization. Most damaging for permeability studies specifically.

**Target design:** a counter-based RNG (Philox/Threefry) keyed on
`(walker_id, step, purpose)`. Deterministic regardless of thread count, and the
RNG you'll want on the GPU anyway. Kill every `rand()` and every per-call
`random_device`.

---

## Phase 0 task list (priority order)

### P0.0 — Build & tooling (done / in progress)
- [x] `CMakeLists.txt` mirroring `compile.sh`, exporting `compile_commands.json`
  for the clangd LSP plugin. Binaries land at repo root as before.
- [ ] `.gitignore` for `build/` and the `compile_commands.json` symlink.

### P0.1 — RNG unification (the keystone task)
- [x] Introduce a single seeded RNG abstraction (`src/rng.h`, `RandomEngine`):
  `std::mt19937_64` + splitmix64 key-mixing `seedFrom(base, keys...)`, ready to
  become counter-based (Philox/Threefry) keyed on `(walker, step, purpose)`.
- [x] Route per-thread engine seeding via `seed + worker_index`
  (`parallelmcsimulation.cpp`, both sub-sim loops). Workers now draw
  independent ensembles instead of duplicating one.
- [x] Route walker placement RNG through the seeded engine
  (`walker.cpp` `setRandomInitialPosition`, `dynamicsSimulation.cpp`
  `getAnIntra/ExtraCellularPosition`). Per-sim seed resolved once in
  `initBaseSeed()`.
- [x] Route gamma substrate RNG through the seed (`cylindergammadistribution.cpp`,
  `spheregammadistribution.cpp`; `seed` member set from `params.seed`).
- [x] Replace `rand()` in the percolation draws with a seeded, per-walker engine
  (`Walker::rng`, used in `cylinder.cpp`, `sphere.cpp`, `plyobstacle.cpp`).
- [x] Acceptance test passes: same seed → bit-identical trajectory (1 and 2
  procs verified; different seed differs). `ctest -R rng_reproducibility`.
- [ ] **Remaining:** bit-identical results *across different thread counts*
  (1 vs 2 vs N). Today placement/steps draw from a per-worker engine, so the
  output is reproducible for a *fixed* `num_process` but not invariant to it.
  True invariance needs placement+steps keyed on a *global* walker id (the
  per-walker `Walker::rng` already is, for the crossing draw). Deferrable.

### P0.2 — Validation suite + CI
- [ ] CMake `add_test` targets asserting analytic limits:
  free diffusion `S = exp(-bD)`; restricted cylinder (van Gelderen); restricted
  sphere (Neuman / Murday-Cotts); **two-compartment exchange (Kärger)** for a
  known permeability `kappa` (directly serves the permeability research line).
- [ ] GitHub Actions workflow building with CMake and running the tests.

### P0.3 — Permeability correctness (research line 1 foundations)
- [x] **Fix the `perm_crossed_flag` one-way latch** — was set false only at init
  and never reset (per-step reset was commented out), permanently disabling the
  numerical-leak sentinel after a walker's first legitimate crossing. Now cleared
  each step with `initial_location` rebaselined to the current compartment
  (`dynamicsSimulation.cpp` step loop). No-op for impermeable runs. (commit fc39fb7)
- [x] Use the per-obstacle `d_intra` in `updateStepLength` (was computed then
  discarded in favour of global `params.diff_intra`); now used when set (>0),
  else falls back to the global. (commit fc39fb7)
- [ ] Apply `sqrt(D_new/D_old)` rescaling to the remaining sub-step at a crossing
  event and update compartment membership immediately (not next step).
  **Changes permeable physics — needs a permeable baseline first.**
- [ ] Scale `obstacle_permeability` in `scale_from_stu` and document its units
  (`parameters.cpp`). **Changes permeable magnitude/units.**
- [ ] Implement permeability at mesh edges/vertices
  (`plyobstacle.cpp` currently forces reflection there). **Changes permeable physics.**

### P0.4 — Release-build correctness landmines
- [ ] Move real work out of `assert()`: `assert(fread(...))`
  (`trajectory.cpp:373`) and `assert(checkScheme/PLY)` (`simerrno.cpp:92,97`)
  vanish under `-DNDEBUG`.
- [ ] Replace the `print + assert(0)` validation idiom with real error
  paths/exceptions throughout `simerrno.cpp`.
- [ ] Fix `PLYObstacle` raw-pointer ownership (`plyobstacle.h:28-29`): no
  destructor/copy-ctor while stored in a `std::vector` → leak + double-free.
  Move to `std::vector<Vertex>`/`std::vector<Triangle>` or rule-of-five.
- [ ] Fix validation bugs: `checkSphereListFile` opens `cylinders_files`
  (`simerrno.cpp:611`); duplicated `diff_extra` check, `diff_intra` unchecked
  (`:70`).
- [ ] Close the APGSE null-deref (`mcsimulation.cpp:84`) — parses as valid but
  no class is constructed.

### P0.5 — Hygiene (opportunistic)
- [ ] Replace `std::endl` with `'\n'` in trajectory/IO hot loops
  (`trajectory.cpp:231,242`).
- [ ] Hoist magic `3600` (phase-histogram bins) and scattered tolerances into
  named constants.
- [ ] Optional `double`-precision output (currently all output downcast to float).

---

## Deferred to later phases (explicitly NOT Phase 0)
- **GPU port** — separate branch. Restructure obstacles to struct-of-arrays +
  type-id dispatch (no virtuals), keep `FixedGrid`/CSR as device accel structure,
  CUDA walker kernel with online phase accumulation. Depends on P0.1 + P0.2.
- Sequence factory/registry (OGSE, double-diffusion, STEAM as first-class types).
- Python bindings (pybind11), HDF5/NIfTI output.
- Physics extensions: T1 relaxation, susceptibility/off-resonance, more geometries.

---

## Immediate next step
P0.1 — the RNG unification — starting from the seeding sites listed above, with
the P0.2 reproducibility acceptance test written first so we can prove the fix.
