# 3 permeable PLY spheres — CPU vs GPU validation + permeability equilibrium

**Substrate:** 3 separate PLY spheres (icospheres, 482 v / 960 f each, from the
`meshed_sphere` array), r = 5 µm, centres along z at 0 / 21 / 42 µm, in a
24×24×62 µm voxel. Driven through `ply_extended_file_list` with **DISTINCT per-sphere
Di, T2 and permeability** — fully supported by the mesh path on both CPU and GPU
(per-obstacle `obs_step` from Di, `obs_T2`, `prob_i_e/prob_e_i` from κ); no code change.

| sphere | Di [m²/s] | T2 [s] | κ [m/s] |
|---|---|---|---|
| A (z=0)  | 0.2e-9 | 0.040 | 1.0e-4 |
| B (z=21) | 0.5e-9 | 0.060 | 2.0e-4 |
| C (z=42) | 1.0e-9 | 0.080 | 5.0e-5 |

Extra (free water) D = 2.0e-9, T2 = ∞. Protocol: PGSE b-sweep 0..3000 s/mm²
(`PGSE_bsweep.scheme`), T=1000, duration 50 ms, periodic BC, uniform full-voxel seeding.

## CPU vs GPU validation — PASS

Distinct-{Di,T2,κ} config (`three_perm_spheres.conf`, N=500k):
GPU intra-frac **0.0423** vs CPU **0.0416** (|Δ|=0.0007); total 3.2 = 3.2; stuck=0.
→ All three distinct per-sphere properties are honored **identically** on CPU and GPU.

Per-sphere apparent ADC is distinct (0.885 / 1.401 / 0.766, units 1e-9 m²/s) — a
physically-correct convolution of intrinsic Di, restriction (all r=5 µm), and κ, so it
does **not** rank by Di alone (sphere C has high Di=1.0 but its diffusion length 12 µm >
the 10 µm diameter → strongly wall-restricted → low apparent ADC).

## Permeability equilibrium density (T2=∞ variant, `three_perm_spheres_eq.conf`)

To read density cleanly, a variant sets **all T2=∞** (so b0 = raw particle count),
keeping distinct Di + κ. Compared against an **impermeable** baseline (κ=0 → intra
fraction = placement = true volume fraction):

- impermeable placement (volume fraction): **4.32 %**
- permeable equilibrium, b0 **count** fraction: **4.64 %**  (CPU↔GPU agree)
- permeable, signal-summed over all b (NOT a density): 6.35 % — inflated because intra
  is restricted and retains high-b signal. Do not read this as a density.

So the Powles permeable membrane maintains an **approximately uniform** density
equilibrium (1.07× globally), with a **small residual diffusivity-dependent enrichment**
of the low-D intra compartments (per-sphere C_in/C_ex, measured from the 100³ volume):

| sphere | D_in/D_ex | C_in / C_ex |
|---|---|---|
| A | 0.10 | 1.16 |
| B | 0.25 | 1.14 |
| C | 0.50 | 1.04 |

Monotonic with the D-contrast (lowest-D sphere most enriched). This is the Powles
model's known D-dependent partition bias — relevant to the open Powles-vs-S/V question.
(Seeding is uniform, so this is the *maintained* equilibrium; to watch the *dynamics*
of equilibration, intra-seeding would be needed — currently CPU-only, GPU seeds uniform.)

## Visualise (GPU volumetric, 50M particles, 46 s; `tests/benchmark/output/`)

- `three_perm_spheres_vol_intra_b0.nii.gz` — ending points **inside** the 3 spheres (the compartments)
- `three_perm_spheres_vol_extra_b0.nii.gz` — ending points **outside**
- `three_perm_spheres_vol_total_b3000.nii.gz` — restriction contrast (spheres bright at high b)
- non-cubic sub-voxel: x,y = 0.24 µm, z = 0.62 µm (affine baked into the NIfTIs)

Note: GPU runs do not fill the CPU-side per-subdivision *density* counters
(`*_volFractions.txt` intra/extra are 0 on GPU); density is read from the separate-signal
b0 above. The empirical crossing counts in `*_perm_counters.txt` are CPU-only (the GPU
does its crossing in-kernel); the prescribed per-sphere `prob_i_e` (0.077/0.096/0.018) are
correctly distinct.
