# Membrane Permeability in Monte‑Carlo dMRI Simulators — Units, Model, and Calibration

**Status:** research note (no code changes yet). Prepared for the P0.3 permeability work.
**Audience:** computational physicist implementing/calibrating permeability in the MC‑DC simulator.
**Method:** deep multi‑source review with adversarial fact‑checking (23/25 extracted claims confirmed by ≥2/3 independent verifiers; 2 refuted claims listed at the end). Primary sources are 2010–2022 peer‑reviewed dMRI literature; see [References](#references).

---

## TL;DR (the answers)

1. **Permeability `κ` is a velocity** — units **length/time** (µm/ms ≡ mm/s, or m/s), defined `κ = D_m/l_m` (membrane diffusivity ÷ thickness, thin‑membrane limit). It enters a **flux / Robin boundary condition**, and the per‑encounter **crossing probability is *derived* from κ**. It is **not** itself a probability.
2. **Our code already implements the canonical Powles/RMS scheme** — the `2.0/3.0` is the 3‑D geometric constant `C_d`, and the `1 + 0.5·(…)` denominator is the exact two‑compartment normalization. **But the user input `percolation` is being used as κ₀ (a velocity), not a probability** — and the simulator currently treats it as the genuine κ, which is only valid when the per‑encounter probability `P ≪ 1` (small Δt).
3. **κ ↔ exchange** is governed by the **surface‑to‑volume ratio**: exchange rate `k = κ·(S/V)`, residence time `τ = (V/S)/κ`. For a sphere `τ = r/(3κ)`, for a cylinder `τ = r/(2κ) = d/(4κ)`. The Kärger two‑compartment exchange time is `t_ex = 1/(r_ie + r_ei)`.

Physiological scale: **κ ≈ 2–33 µm/s** (NEXI gray‑matter neurites) up to **≈ 50–125 µm/s** (red‑blood‑cell water, fast neurite exchange). So κ ~ 10⁻⁶–10⁻⁴ m/s.

---

## 1. The Powles / transit‑probability model: definition, units, formula, limitations

### 1.1 What "permeability" physically is

In the continuum picture (Novikov & Fieremans), a thin membrane is the limit `D_m → 0`, `l_m → 0` with the ratio held fixed:

```
κ = D_m / l_m            [length²/time ÷ length = length/time]   →  κ has units of VELOCITY
```

It enters the steady‑state **flux (Robin / "leaky") boundary condition** at the membrane surface `r_m` with normal `n`:

```
−n·J|_rm  =  D₀ (n·∇ψ)|_rm  =  κ · [ ψ(r_m⁺) − ψ(r_m⁻) ]
```

i.e. the diffusive flux through the wall is proportional to the concentration jump, with proportionality constant κ. This is the standard definition across the field (Novikov/Fieremans 2011; Alemany/Rose 2022: *"membrane permeability κ_i, with units of distance over time, effectively defines the rate of exchange"*). **Confidence: high (3‑0, 5 primary sources).**

> ⚠️ **Convention flag.** This κ (velocity) is distinct from a *dimensionless transmission probability* `P_T`. Some MC tools (notably **Camino**, `-p`) expose the **dimensionless step‑through probability directly**, which is *not* a permeability in m/s. Our simulator instead scales by `δ_s`, `D`, `C_d`, so the input plays the role of κ₀ (a velocity) — see §1.4.

### 1.2 The per‑encounter crossing probability (the Powles form)

A walker that reaches the membrane crosses with probability (Lee/RMS, citing Powles 1992; Szafer 1995; Fieremans & Lee 2018):

```
P  ≈  (κ₀ · δ_s / D₀) · C_d ,        δ_s = √(2 d D₀ Δt)   (RMS step length, d = spatial dim)

           ⎧ 1     (d = 1)
   C_d  =  ⎨ π/4   (d = 2)
           ⎩ 2/3   (d = 3)
```

Equivalent linearized form (Fieremans 2010): `P_T,i→e = 2 d_s κ / D_i`.
Exact (non‑linearized) form (Novikov/Fieremans App C.1; Alemany/Rose Eq. 7):

```
P = 2 κ δ_x / (D₀ + 2 κ δ_x)     ⇔     P/(1−P) = 2 κ δ_x / D₀
```

where `δ_x` is the walker‑to‑barrier distance. `P ≈ 2κδ_s/D` only when `P ≪ 1`. **Confidence: high (3‑0).**

**The 3‑D factor `2/3` in our code is exactly `C_d` for d = 3.** Our step length `√(6 D Δt) = √(2·3·D·Δt)` is the `d = 3` RMS step, consistent with `C_d = 2/3`.

### 1.3 Two compartments with different diffusivities — the normalization

For an interface between compartments with `D₁, D₂`, the asymmetric crossing probability is (Lee/RMS Eq. 3):

```
                 (κ δ_s1 / D₁) C_d
P_{1→2}  =  ────────────────────────────────────
            1 + (κ/2)(δ_s1/D₁ + δ_s2/D₂) C_d
```

**This is term‑for‑term our code.** Verified directly against `src/parallelmcsimulation.cpp` (≈ lines 670–712):

```cpp
dsi          = sqrt(6*time_step*d_intra);                 // δ_s,i = √(2·3·D_i·Δt)
prob_cross_i_e = percolation * dsi * 2.0/3.0 / d_intra;   // = κ₀ · δ_s,i · C_d / D_i
// …symmetric prob_cross_e_i with D_extra…
prob_cross_i_e_final = prob_cross_i_e
                     / (1 + 0.5*(prob_cross_i_e + prob_cross_e_i));
```

Factoring κ₀ out of the RMS denominator reproduces `1 + 0.5·(…)` exactly, with **`percolation` ≡ κ₀**. **Confidence: high (3‑0; one source is our own code).**

### 1.4 UNITS — the key conclusion for our simulator

Because `prob_cross = percolation · δ_s · C_d / D`, and `δ_s` has units of length while `D` has units of length²/time, **for `prob_cross` to be a (dimensionless) probability, `percolation` must have units of length/time — i.e. it is κ₀, a velocity, not a probability.** It must be supplied in **units consistent with `δ_s` and `D`** (after `scale_from_stu`, the code's `D` is in mm²/ms, lengths in mm, so κ₀ must be in **mm/ms**).

> This is the crux of the units confusion. The config word "permeability"/"percolation" is **κ₀ in length/time**, *not* the Camino‑style dimensionless `-p`. Item P0.4 ("scale `obstacle_permeability` in `scale_from_stu`") is therefore **necessary**: κ₀ has length/time units and must be unit‑converted alongside `D` and lengths, or the percolation probability will be wrong by the unit‑scale factor.

### 1.5 Documented limitations & pitfalls

| Pitfall | Statement | Source / confidence |
|---|---|---|
| **Small‑step validity** | The Powles form is valid only for `P ≪ 1`. The input `κ₀` equals the *genuine* permeability only in that limit. | RMS, Fieremans, Alemany/Rose — **high (3‑0)** |
| **Finite‑step bias** | Genuine `κ = κ₀/(1 − α) > κ₀` (α = mean permeation probability). The discrete scheme **under‑estimates** exchange; the bias grows with Δt and κ. | RMS Eq. C.7 — **high (3‑0)** |
| **Recommended bounds** | Keep max transit probability **< ~0.01** and step/radius **δ_s/R < ~0.1**; highly permeable membranes need small Δt. | Fieremans 2010; Alemany/Rose — **high (3‑0)** |
| **Time‑step dependence** | Camino's `-p` is dt‑dependent and must be recalibrated when dt changes — the generic symptom of the unscaled probabilistic approach. | Camino docs — **high (3‑0)** |
| **Interface reflection / detailed balance** | For `D₁ ≠ D₂`, correct physics requires `p_t(i→j)/p_t(j→i) = √(D_j/D_i)` even as κ → ∞ (from the diffusivity discontinuity, not the membrane). Naïve transit models violate it (ratio → 1). | Alemany/Rose Eq. 4 (cites Szafer/Zhong/Gore 1995) — **medium (2‑1, single source)** |
| **Mesh / S‑V discretization** | Encounters per unit time scale with `S/V` and step length; triangle‑size discretization can over/under‑count encounters and bias effective κ on PLY meshes. | Inferred from S/V law (see §2); **open question — needs a convergence test** |

---

## 2. A surface‑to‑volume‑scaled permeability model (mesh/Δt‑independent target)

### 2.1 The macroscopic link

The single robust, geometry‑independent quantity is the **macroscopic exchange rate**, set by κ and the surface‑to‑volume ratio (RMS, citing Slijkerman & Hofman 1998):

```
k = κ · (S/V)          [ (length/time)·(1/length) = 1/time ]          (exchange RATE)
τ = (V/S) / κ          (single‑compartment residence time)
```

For the canonical shapes:

| Geometry | S/V | Residence time τ = (V/S)/κ |
|---|---|---|
| Sphere, radius r | 3/r | `τ = r/(3κ)` |
| Cylinder, radius r (diam d) | 2/r = 4/d | `τ = r/(2κ) = d/(4κ)` |
| Slab/plane gap w | 1/w (one side) | `τ = w/κ` |

**Confidence: high (3‑0).** Example: cylinder r = 1 µm, τ = 4 ms ⇒ κ = r/(2τ) = 125 µm/s (Olesen).

### 2.2 Why this gives mesh/Δt independence

The number of membrane encounters per unit time in an MC walk scales with `S/V` and the step length, while each encounter crosses with `P ≈ 2κδ_s/D`. Their product reproduces `k = κ·(S/V)` **independent of how the surface is triangulated**, *provided* `P ≪ 1` (small Δt). So the calibration recipe is:

> **Fix the physical target `κ` (in length/time). Derive the per‑encounter probability `P = (κ δ_s/D)·C_d` from κ, δ_s(Δt), D. The simulated exchange rate then converges to `k = κ·S/V` as Δt → 0, regardless of mesh resolution.**

This is exactly the structure already in the code (§1.3); the missing pieces are (a) treating `percolation` explicitly as κ₀ with correct units (§1.4), and (b) a **finite‑step correction** so the *realized* κ matches the target:

```
κ_target = κ₀ / (1 − α)      ⇒  to hit κ_target, set  κ₀ = κ_target·(1 − α)
```

with α the mean per‑encounter probability (≈ P at the operating Δt). For `P ≪ 1`, κ₀ ≈ κ_target and no correction is needed — which is why the small‑Δt bound matters.

### 2.3 Practical calibration knobs

- Given a desired κ and a chosen Δt: `δ_s = √(6 D Δt)`, `P = (2/3)·κ·δ_s/D`. **Choose Δt so P < 0.01.**
- To make exchange reproducible across meshes, validate that the **measured** exchange rate (e.g. from intra/extra fraction decay, or an apparent exchange‑rate fit) equals `κ·S/V` and is **flat vs triangle count** — this is the convergence test to add.

---

## 3. Permeability ↔ exchange / residence time (Kärger and tissue‑exchange dMRI)

### 3.1 Kärger two‑compartment model

```
t_ex = 1 / r ,     r = r_ie + r_ei        (sum of inter‑compartment exchange rates)
detailed balance:   f · r_ie = (1 − f) · r_ei      (f = intra fraction)
```

With S/V scaling, the apparent (two‑compartment) exchange time carries an extracellular‑volume‑fraction prefactor (Alemany/Rose):

```
τ_ex = α · (V/S) / κ        (α = extracellular volume fraction)   ← FEXI/AXR "apparent" exchange time
```

versus the **bare single‑compartment** residence time `τ_i = (V/S)/κ` (no α). **Do not conflate these** — papers differ. **Confidence: high (3‑0), with this convention caveat.**

### 3.2 Neurite (NEXI/SMEX) form

For cylindrical neurites (Jelescu/NEXI):

```
P = d / [ 4 ( t_ex − d²/(32 D∥) ) ]   ≈   d / (4 t_ex)
```

i.e. the cylinder result `τ = d/(4κ)` inverted, with a small intra‑cylinder mixing correction. **Confidence: high (3‑0).**

### 3.3 Quantitative anchors (sanity targets for the simulator)

| Scenario | Geometry | Exchange/residence | Permeability κ |
|---|---|---|---|
| Healthy GM neurites (NEXI) | d ≈ 0.5–2 µm | t_ex ≈ 15–60 ms | **P ≈ 2–33 µm/s** |
| Fast neurite exchange (Olesen) | r = 1 µm | τ ≈ 4 ms | **κ ≈ 125 µm/s** |
| Red‑blood‑cell water (reference) | — | — | **49–112 µm/s @ 37 °C** |

**Confidence: high (3‑0).** These confirm κ lives at the **µm/s (≈ 10⁻⁶–10⁻⁴ m/s)** scale — a useful range‑check for whatever units the config ends up expecting.

---

## 4. Implications & recommendations for the MC‑DC simulator

1. **Document and treat `percolation`/`permeability` as κ₀ in length/time** (after `scale_from_stu`: **mm/ms**), not a probability. This resolves the units question and tells us P0.4 (scaling κ₀ in `scale_from_stu`) is *required for correctness*, not cosmetic.
2. **Add a small‑step guard / diagnostic.** Compute the per‑encounter `P = (2/3)κ₀δ_s/D` at setup and **warn if P ≳ 0.01** (and/or `δ_s/R ≳ 0.1`). Optionally report the *realized* κ and the finite‑step‑corrected `κ₀ = κ_target(1−α)`.
3. **Report exchange‑time / S–V diagnostics.** Given the substrate's S/V, print the implied `τ = (V/S)/κ` and `k = κ·S/V` so users can sanity‑check against the §3.3 anchors. Decide and document **which** exchange time is reported (bare `τ_i` vs α‑weighted `τ_ex`).
4. **Mesh‑independence test (new regression).** Run the *same* κ on coarse vs fine triangulations of the same sphere and confirm the measured exchange rate ≈ `κ·S/V` and is flat vs triangle count. This is the permeable analogue of the impermeable golden test.
5. **Check the interface‑reflection ratio** `p_t(i→j)/p_t(j→i) = √(D_j/D_i)` when `D_intra ≠ D_extra` (medium‑confidence requirement) — numerically test for spurious membrane accumulation; if violated, the Alemany/Rose *hybrid* model is the fix.
6. **Validation target:** a permeable two‑compartment run should reproduce the **Kärger** signal for the κ‑implied `t_ex` — the natural "permeability physics as intended" test the user asked for.

---

## 5. Open questions (carry into implementation)

- Does/should the code apply the finite‑step correction `κ = κ₀/(1−α)`, or expose a calibrated‑κ diagnostic? What is the realized κ for the user's typical Δt, D?
- On triangulated PLY meshes, does the per‑encounter probability reproduce `k = κ·S/V` independent of triangle size, or does surface discretization bias effective κ? (→ convergence test, rec. #4.)
- Is `p_t(i→j)/p_t(j→i) = √(D_j/D_i)` satisfied by the current symmetric `δ_s,i/δ_s,e` construction for `D_intra ≠ D_extra`?
- Which exchange‑time convention (bare `τ_i` vs α‑weighted `τ_ex`) should be reported for NEXI/FEXI validation, and is the intended κ the **diffusional (tracer)** or **osmotic** permeability? (NEXI's P is the *diffusional* one, explicitly < osmotic.)

---

## 6. Caveats & convention differences (read before coding)

1. The input `percolation` is the Powles **κ₀**; it equals the genuine κ only for `P ≪ 1`. For high κ the genuine `κ = κ₀/(1−α)` is larger ⇒ the simulator **under‑estimates exchange** unless Δt is small.
2. `percolation` carries **velocity units** once interpreted as κ₀ — unlike Camino's dimensionless `-p`. Verify the user value is length/time consistent with δ_s and D.
3. S/V residence time appears **with and without** the α (ECV) prefactor: `τ_i = (V/S)/κ` (bare) vs `τ_ex = α(V/S)/κ` (FEXI/AXR apparent). Don't conflate.
4. The interface‑reflection detailed‑balance condition rests on a **single source, 2‑1 vote (medium confidence)** — test numerically.
5. NEXI's P is the **diffusional (tracer)** water permeability, smaller than the **osmotic** permeability — match modeling intent.
6. `P_T = 2 d_s κ/D` is the **small‑P linearization** of the exact `P/(1−P) = 2κδ_s/D`; use the full form when P is not tiny.

*Time‑sensitivity:* all primary sources are 2010–2022 and current as of 2026. The Camino claim relied on a search cache of the docs page (direct fetch blocked by a self‑signed TLS cert) but multiple retrievals agreed.

---

## 7. Validation experiments (proposed)

These tests validate the *transport* behaviour of the permeability model before comparing DWI signals. The goal is to verify that the model preserves compartment volume, reproduces the target membrane permeability, and remains stable across time step, diffusivity contrast, and substrate geometry.

### 7.1 Compartment particle fractions

Let `i` denote the intra-compartment and `e` the extra-compartment.

- `N_i(t)`: number of walkers inside the intra-compartment at time `t`.
- `N_e(t)`: number of walkers inside the extra-compartment at time `t`.
- `N(t) = N_i(t) + N_e(t)`: total number of valid walkers.

The measured intra-compartment particle fraction is:

$$
f_i(t) = \frac{N_i(t)}{N_i(t) + N_e(t)}
$$

For a substrate with true intra-volume fraction `f_true`, volume preservation requires:

$$
f_i(t) \rightarrow f_{\mathrm{true}}
$$

after sufficient mixing time.

### 7.2 Crossing counts and membrane hits

Over a simulation interval `[0,T]`, define:

- `C_i→e`: number of successful intra-to-extra crossings.
- `C_e→i`: number of successful extra-to-intra crossings.
- `H_i`: number of membrane hits from the intra side.
- `H_e`: number of membrane hits from the extra side.

The empirical per-encounter crossing probabilities are:

$$
\hat{p}^{enc}_{i\to e} = \frac{C_{i\to e}}{H_i}
$$

$$
\hat{p}^{enc}_{e\to i} = \frac{C_{e\to i}}{H_e}
$$

### 7.3 Implementation prerequisites (not yet in the code)

The quantities above are **not currently produced** by the simulator. Before these
validation experiments can run, the following must be implemented. Listed with the
nearest existing hook so the gap is concrete. *(Notes only — no code written here.)*

1. **Per-direction crossing counters `C_{i→e}`, `C_{e→i}`.** Today each obstacle has a
   single `count_perc_crossings`, incremented on *any* successful crossing
   (`cylinder.cpp`, `sphere.cpp`, `plyobstacle.cpp`) — not split by direction and not
   aggregated/output per ensemble. The crossing site already knows the direction via
   `colision.col_location` (inside/outside), so it is a matter of two counters instead
   of one, plus aggregation across walkers/threads.

2. **Per-direction membrane-hit counters `H_i`, `H_e`.** Not tracked at all. A "hit" =
   a step whose collision lands on a permeable membrane (where the percolation draw is
   made), whether or not it crosses. Needs a counter incremented in the same
   `if(percolation>0)` block, by side, *before* the crossing decision. Required for the
   empirical per-encounter probabilities `p̂ᵉⁿᶜ` in §7.2.

3. **Time-resolved compartment populations `N_i(t)`, `N(t)`.** Partially present:
   `global_particle_positions[t]` already accumulates an intra count per step
   (`walker_particle_positions`, added recently), but it is marked "for debugging, todo:
   remove", its semantics (does it count only *valid* walkers?) need confirming, and it
   is not written to any output file. Need: confirm semantics and expose `N_i(t)`, `N(t)`
   so `f_i(t)` (§7.1) can be formed; `N_e(t) = N(t) − N_i(t)`.

4. **Substrate surface area and volume → `S/V`.** Volume is computed (volume
   subsampling / subdivisions); **surface area is not exposed** (PLY: sum of triangle
   areas; sphere/cylinder: analytic). Both `S` and `V` must be reported so the predicted
   `k = κ·S/V` and `τ = (V/S)/κ` can be compared against the measurements.

5. **A transport-diagnostics output.** The counters and populations above need to be
   written out (e.g. a per-step `*_transport.txt` and an end-of-run summary carrying
   `C_{i→e}, C_{e→i}, H_i, H_e, N_i(t), S, V`). Only DWI / trajectory / info outputs
   exist today.

6. **κ as a units-correct input (depends on P0.4).** These tests assume `percolation` ≡
   κ₀ in length/time, unit-scaled in `scale_from_stu` (see §1.4). Until P0.4 lands the
   "target κ" is ill-defined, so §7 results cannot be compared *quantitatively* to
   `κ·S/V`. A **calibrated-κ / small-step diagnostic** (compute `P = (2/3)κ₀δ_s/D` at
   setup, warn if `P ≳ 0.01`; §4 rec. #2) belongs here too.

7. **Prescribed crossing probabilities exposed.** `prob_cross_i_e`, `prob_cross_e_i` are
   computed at setup (`parallelmcsimulation.cpp`) but not output. Exposing them lets §7.2
   compare the empirical `p̂ᵉⁿᶜ` to the prescribed values and check the interface-reflection
   ratio `p̂_{i→e}/p̂_{e→i} ?= √(D_e/D_i)` for `D_intra ≠ D_extra`.

Already available (no work needed): an out-of-equilibrium start via `ini_walkers_pos intra`
(for the §7.1 mixing test), and `diff_intra`/`diff_extra` for the diffusivity-contrast /
interface-reflection check. The analytic **Kärger** reference for a signal-level test is a
*post-processing* (Python test-harness) dependency, not a simulator change.

## References

Verified primary sources (with the research angle they supported):

1. **Novikov & Fieremans**, *Relating extracellular diffusivity to cell size…* — arXiv:1004.2701 (Nature Physics 2011). κ = D_m/l_m; Robin BC; exact `P/(1−P) = 2κδ/D`.
2. **Lee, Fieremans, Novikov — Realistic Microstructure Simulator (RMS)** — arXiv:2012.06478. `P ≈ κ₀δ_s/D·C_d`, `C_d = {1, π/4, 2/3}`; two‑compartment normalization (Eq. 3); finite‑step bias `κ = κ₀/(1−α)` (Eq. C.7); `k = ρ·S/V`.
3. **Fieremans, Novikov, Jensen, Helpern**, *Monte‑Carlo study of a two‑compartment exchange model of diffusion* — NMR Biomed 2010, doi:10.1002/nbm.1577. `P_T = 2 d_s κ/D`; dt bound (transmission < 0.01, dr/R < 0.1).
4. **Alemany, Rose et al.**, hybrid permeability model — arXiv:2201.09609 / PMC9232609. Exact transit form (Eq. 7); interface reflection `p_t ratio = √(D_j/D_i)` (Eq. 4); `τ_ex = α(V/S)/κ`; small‑Δt requirement of the reference model.
5. **Jelescu et al. — NEXI** (anisotropic Kärger) — arXiv:2108.06121 / PMC10363376. `t_ex = 1/(r_ie+r_ei)`, detailed balance; `P = d/[4(t_ex − d²/32D∥)] ≈ d/(4t_ex)`; κ ≈ 2–33 µm/s.
6. **Olesen et al.**, diffusion exchange — PMC8961002. Cylinder `τ = r/(2κ)`; r = 1 µm, τ = 4 ms ⇒ κ = 125 µm/s; RBC water 49–112 µm/s.
7. **Camino MC simulator docs** — http://camino.cs.ucl.ac.uk/index.php?n=Tutorials.MCSimulator. `-p` = dimensionless step‑through probability, dt‑dependent, needs recalibration.
8. Additional corroborating primary sources: mrm.22782, mrm.28714, Nature s41598‑022‑14541‑y, PMC8026575, PMC2997614, ScienceDirect S1090780710001771, bioRxiv 2025.07.16.664944.

### Refuted / low‑confidence (excluded from conclusions)
- *"τ_i = R²/κ form"* via geometry from nbm.1577 — **refuted 0‑3** (the correct form is `τ = (V/S)/κ`, i.e. `r/(cκ)`, linear in r, not R²/κ).
- *"mean neurite residence time ≈ 4 ms in gray matter"* as a general value — **refuted 1‑2** (4 ms was a specific worked example, not a consensus GM value; consensus t_ex is ~10s of ms).
