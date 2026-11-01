# Analytical validation tier

Per-feature validation against **analytical oracles** at Monte-Carlo-converged N.
Each page below states the physics, the oracle and *why it is valid* (including
the regime where it is **not**), the fixture, the tolerance and its
justification, and how to read a pass/fail.

These tests are `pytest`-based (`pip install .[test]`) and run in the
`validation-analytic` CTest tier. They are **slow on CPU, fast on GPU**; the
GPU variants **self-skip** without CUDA. The binary is located via
`$MCDC_BIN` (set by CTest) or the repo-root default.

> **Golden vs analytical.** The `accuracy/` golden tests pin the signal
> bit-exactly to *this build's* own output — they catch drift, not physics
> errors. The tests here compare to an **independent** analytical truth, so they
> catch physics errors a golden cannot.

---

## ✅ Free diffusion — `test_free_exp.py`

**Physics.** With no obstacles the PGSE signal is Gaussian: `S/S0 = exp(-b·D)`
for every measurement (`b` from the scheme, `D` from the config, both SI).

**Oracle.** Closed-form `exp(-b·D)` (`mcdc_testkit.oracles.free_exp`) — exact, no
dependencies. Free diffusion is **isotropic**, so the physical quantity is the
**per-shell mean** (every gradient direction in a shell shares the same
`exp(-b·D)`); per-measurement scatter is pure MC noise. The test asserts on the
per-shell mean and keeps a loose per-measurement bound only as a gross-outlier
guard.

**Tolerance.** A few percent: fixed-length MC steps approximate Gaussian
diffusion to `O(1/T)` (a small bias that grows with `b`), on top of MC noise
`~1/√N`. Verified on the current build: per-shell mean deviation **0.05 % at
N=50 000** (physics essentially exact — no collisions), and **2.5 % at the smoke
size N=5 000, T=1 000** (noise-dominated, deterministic at the fixed seed).

**Read it.** `worst per-shell mean dev` is the headline; `per-meas worst` is the
outlier guard. Both are printed. Runs on CPU and GPU (GPU self-skips without
CUDA). Also the fast **smoke** gate, via
`python3 tests/validation/test_free_exp.py --n 5000 --t 1000 --gpu 0 --rtol 0.05`.

---

## ✅ Restricted cylinder — `test_cylinder_vangelderen.py`  *(analytic curve must be frozen first)*

**Physics.** Diffusion inside an infinite cylinder is restricted perpendicular to
the axis and free along it, producing a characteristic anisotropic, b- and
angle-dependent attenuation.

**Oracle.** dmipy `C4CylinderGaussianPhaseApproximation` — the **van Gelderen
Gaussian-Phase-Approximation** for an infinite cylinder.

> **Validity limits — important.** The van Gelderen GPA cylinder is valid only
> for:
> - an **infinite, impermeable** cylinder (no exchange, no caps),
> - **small b-values (≲ 5000 s/mm²)** — the Gaussian-phase approximation degrades
>   at high b,
> - **perpendicular** restriction (the axis is the simulator's z-axis; native
>   cylinders are always z-oriented).
>
> So this validates restricted diffusion and its anisotropy — **not**
> permeability, high-b behaviour, or finite/capped geometry. Permeable and
> high-b cases are validated by other oracles (crossing-counter, cross-version).

**Fixture.** `accuracy/cylinder_restricted.conf` (single z-cylinder, r = 5 µm,
walkers initialised intra, T2 ≈ ∞). The frozen curve is produced **once** by
`tests/tools/gen_analytic.py` in a dmipy environment (`pip install .[analytic]`)
and committed to `accuracy/analytic/` with a provenance JSON; the run-time test
loads the `.npy` and needs only numpy.

**Tolerance.** Per-shell NRMSE within a statistical band (MC noise + finite-step
bias), evaluated **only on the small-b shells** where the GPA holds.

**Read it.** The test does two things: (1) an **always-on, dmipy-free anisotropy
sanity check** — restricted (perpendicular) shells must attenuate more than the
near-parallel directions — and (2) the **analytic comparison**, which
**self-skips** with a pointer to `gen_analytic.py` until the frozen curve exists
(a fresh checkout, or a curve not yet generated, is a skip, not a failure).

---

## ✅ Gamma packing generator — `test_gamma_packing.py`

Validates the **substrate generator**, not the walk: the simulator's built-in
gamma packers (`{cylinder,sphere}gammadistribution.cpp`) draw radii from a
Gamma(α, β) law, **auto-size** the voxel from the requested ICVF, pack
non-overlapping obstacles with periodic images, and write the realised packing to
`<prefix>_gamma_distributed_{cylinder,sphere}_list.txt`.

**Fixtures (tissue-scale, ~100 µm voxel — the calibration target):**
`accuracy/gamma_cylinder_100um.conf` (α1.5, β1e-6 m, 340 cyl → mean r **1.5 µm**)
and `accuracy/gamma_sphere_100um.conf` (α1.5, β2e-6 m, 900 sph → mean r **3.0 µm**).
Radii are **SI/metres**; **mean radius = α·β**; `min_radius` default 0 (no
truncation). The voxel is *derived* — `L=√(Σπr²/icvf)` (cyl) / `(Σ⁴⁄₃πr³/icvf)^⅓`
(sph) — so you hit ~100 µm by tuning `num_obstacles`.

**Oracles (geometric / statistical invariants — no analytic DWI needed):**
achieved ICVF (recomputed from the list == requested == the engine's report);
radius distribution (sample mean ≈ α·β pins the scale, CV ≈ 1/√α pins the shape);
non-overlap of all obstacles incl. periodic images; the auto-sized voxel matches
the closed form; and **determinism** (same `seed`+`num_process` → byte-identical
list, different seed → different — guards the P0.1 seeding fix). Numpy-only, CPU,
`tier2_fast`. Frozen example packings live in `accuracy/gamma_lists/`.

---

## ✅ Packing restricted diffusion — `test_packing_restricted.py`  *(experiments 1 & 2)*

Runs the calibrated tissue-scale packings intra-only + impermeable and compares the
Monte-Carlo signal to the **volume-weighted GPA** over the packing's own radii:
`S_ana(b) = Σ_i V_i · S_GPA(b; r_i) / Σ_i V_i` (spheres V∝r³ Neuman; cylinders V∝r²
van Gelderen, per gradient direction). Exact in expectation because uniform intra
seeding puts a walker count in each obstacle ∝ its volume.

**Oracle.** Pure-numpy GPA (`oracles.sphere_neuman`, `oracles.cylinder_vangelderen`,
`oracles.packing_signal`) — the SAME model as dmipy `S4Sphere…`/`C4Cylinder…GaussianPhaseApproximation`
(see `../tools/gen_analytic.py`),
re-implemented so it needs no dmipy and is evaluated live from each packing's radii
(the curve is specific to the exact packing, so there is nothing to freeze). Validated
against the single-obstacle MC (≤0.006) and the packings (≤0.0015 sphere / ≤0.0008 cyl,
N=30000). **Validity:** GPA = impermeable, small b (≲5000 s/mm²) — the scheme tops out
at 2000. `--plot DIR` writes a per-shell overlay (markdown always; PNG if matplotlib).

## ✅ Packing exchange — `test_packing_exchange.py`  *(experiments 3 & 4)*

Same substrates, permeable / two-compartment. No closed form for a polydisperse
permeable packing, so the oracle is a combination: **κ→0 recovers the restricted GPA**
(the strong anchor); the extra fraction and total ADC **increase monotonically with κ**;
the `separate_signals` conservation identity `S_tot = f_i·S_intra + f_e·S_extra`; and the
extra compartment is **hindered** (`0 < ADC_extra < d_extra`). A Kärger overlay is drawn
as reference only (not pass/fail — restricted intra is non-Gaussian). Groups with exp 1/2
on one master plot per geometry. Measured: κ=0 intra vs GPA ≤0.004; exchange rises 0→50%
across κ∈{0,1e-5,1e-4}; conservation exact.

Note: `test_sphere_neuman.py` / `test_cylinder_vangelderen.py` (single obstacle) now fall
back to the **live pure-numpy GPA** when the frozen dmipy curve is absent, so they run in
a numpy-only env instead of skipping.

## ✅ Real-mesh consistency — `test_mesh_cell.py`

A small **real** triangulated cell (`accuracy/meshes/voronoi_cell.ply`, ~66k tris,
~42 µm span) imported from the tissue-forge astrocyte "voronoi" meshes. VTK exports
carry per-vertex normals + RGB colour, which the simulator's PLY reader (xyz-only,
`unitMesh`-style) misreads → "Non-triangle face" → segfault; `tests/tools/import_mesh.py`
strips a source PLY to **xyz + triangle faces only** (fan-triangulating any polygon)
so it ingests cleanly. `ply_scale 1e-6` (coords are µm); voxel auto-sized to the bbox.

**Oracle = consistency** (no closed form for a real cell): (1) CPU sanity — the
signal is monotone in b and the ADC is bounded by the free value (mild restriction —
the cell is large); (2) **CPU ≡ GPU parity** — both engines walk the identical mesh,
per-shell means agree within MC noise (self-skips without CUDA); (3) **cross-version**
— the same conf on the historical binaries in `../crossversion/binaries.toml`, per-shell
within a loose band (self-skips if absent). Mesh walking has no triangle-acceleration
structure (each step tests all ~66k tris), so the CPU is slow (~7 min at N=1500) while
the **GPU is ~1000× faster (~0.4 s)** — the whole point of the GPU path for meshes, and
why a tiny toy cell is the right fixture. The test runs one CPU baseline at a small N
(~2-3 min) and adds the GPU parity when CUDA is present; measured CPU≡GPU |Δ|=0.01.
Import more cells with `import_mesh.py --dir SRC --out-dir DST`.

## ✅ Time-dependence — `test_time_dependence.py`

The physics the rest of the suite validated only at a single diffusion time. Uses
`accuracy/nexi_time_dependence.scheme` — a 4-Δ subset (Δ ∈ {10,20,30,45} ms, TE=58 ms,
so the conf needs `duration ≥ 0.060`) of a NEXI protocol — on one 5 µm sphere as **both**
a native list and a PLY mesh (same geometry → one oracle).

**Restriction × Δ** (impermeable): the signal matches Neuman GPA(Δ) at b ≤ 3000 (GPA valid
only at small b; b=5000/10000 drive the parity + monotonicity checks) and the **ADC decreases
with Δ** — restriction becomes visible as the diffusion time grows.

**Permeability × Δ** (κ ∈ {0,1e-5,5e-5,2e-4}): the exchange effect E = S_impermeable −
S_permeable grows with **κ AND Δ** — permeability is nearly invisible at short Δ, clearly
visible at long Δ (the NEXI idea); f_extra grows with κ. **A quality test, not a sweep** —
a few meaningful κ/Δ levels, the one physically-coupled interaction. Native ≡ PLY and
CPU ≡ GPU throughout (the PLY parity flushed out the GPU mesh intra-seeding bug, now fixed).

Companion "meaningful-levels" checks: `test_free_exp.py` recovers ADC == D at a low and a
high diffusivity; `test_t2_decay.py` verifies exp(-TE/T2) at a low (40 ms) and high (150 ms) T2.

## ⏳ Planned (Phase B+)

| Feature | Oracle | Notes |
|---|---|---|
| Restricted **sphere** | Neuman / Murday-Cotts | small-b, impermeable; SPH master chain (analytic ≡ list ≡ PLY) |
| **T2** decay | `S(b0) ∝ exp(-TE/T2)`; multi-TE exponential | promote from golden to analytic |
| **Permeability** (exchange signal) | **Kärger** two-compartment | slow-exchange regime only; **no closed form** for general cyl/sphere exchange → otherwise crossing-counter + cross-version |
| **Cross-representation** | list ≡ ply ≡ hex-packing | consistency oracle (no analytic needed) |
| **CPU/GPU parity** | cross-hardware, `~k/√N` | parametrized over all GPU-eligible cases; collapses the hand-written `gpu_*` CTest entries |

See the top-level [`tests/README.md`](../README.md) for the test families and how to read
their results.
