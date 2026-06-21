# Performance & capability report — 2.0.3.beta vs 2.1.2.beta

Comparison of the MC-DC simulator at **`db833fb` (2.0.3.beta, 2025-05-09)** vs
**HEAD (2.1.2.beta, this session)**. Performance-focused (not bit-exact; different
RNG between versions, by design).

## Method
- Both built from source with **identical flags** via the shipped `compile.sh`:
  `g++ -O3 -std=c++17 -march=native -flto -w` (asserts on for both; no `-DNDEBUG`).
  2.0.3's `compile.sh` source list was stale (missing `AABBFixedGrid.cpp`), so both
  were compiled as `main.cpp` + all non-`main` sources.
- Same machine (16 cores, 62 GB), runs sequential (no contention), `/usr/bin/time -v`.
- **Same workload**: PorusMedia.ply (82,688 triangles), `N=20000`, `T=1500`,
  `duration=0.050 s`, `D=0.6e-9 m²/s`, the same 270-direction PGSE scheme,
  `num_process 10`, walkers initialised **intra**, trajectory output off, DWI on.
- Identical inputs (same mesh + scheme files). 2.0.3 needed the voxel supplied
  **pre-scaled ×1000** — see "config differences" below — to land on the same
  internal geometry; both then placed walkers intra and completed (exit 0).

## Performance result

| | 2.0.3.beta | 2.1.2.beta | change |
|---|---|---|---|
| wall clock | **52.63 s** | **21.35 s** | **2.46× faster** |
| sim time (avg/worker) | 48 s | 19 s | 2.5× faster |
| peak RSS | 125.9 MB | 128.2 MB | ~same |
| result | b0 intra placement OK | b0 = 19828.7 | — |

**~2.5× faster on the representative mesh workload, same memory footprint.** The
gain is dominated by the **B4 optimisation** (skip phase-shift timesteps where the
PGSE gradient is off): signal synthesis was ~78 % of runtime and is now a small
fraction. Smaller contributions: B1 (reuse the AABB grid-query buffer instead of a
per-step `unordered_set`), B3 (drop an unused per-walker VLA), B2 (trajectory I/O,
not exercised here). Single run each; expect ~5 % wall-time variance.

## Capability differences

| capability | 2.0.3.beta | 2.1.2.beta |
|---|---|---|
| **Load a gigantic mesh** (hp_outer_erode, 193M tris, 3.5 GB) | **OOM** — has the identical double-AABB-copy (`InitializeGrid: this->aabbs = aabbs`) and never frees the PLY AABBs → ~2× AABB memory (not run here, per request) | **Loads** — peak RSS 47.4 GB after the `fix(mem)` that drops the duplicate AABB copy + frees the transient (measured earlier this session) |
| Reproducible RNG (`seed`) | no unified seeding | **yes** (P0.1: one seeded engine; same seed+thread-count → identical output) |
| Config units | `scale_from_stu` scales **time & diffusivity only**, NOT the voxel/geometry | **SI by default** (scales time, diffusivity AND geometry); `use_mm_ms` opt-out; `scale_from_stu` kept as a deprecated alias |
| Input validation | weaker; `assert`-based (no-op under `-DNDEBUG`) | **enforced in release** (P0.4: clean `exit` on malformed config) |
| Per-mesh properties | global only | `ply_file_list` / `ply_extended_file_list` (per-mesh `d_intra`, `T2`, permeability) |
| Permeability | basic | validated Powles model + directional hit/cross counters |
| `#` config comments, documented units header | no | yes |
| Regression/validation tests | none in-tree | **8-target ctest** (RNG repro, mesh/sphere/cylinder goldens, permeability, per-obstacle props, format checks) |
| Per-step debug trace (`debug`) | no | yes (pos, compartment, Di, T2) |

## Errors / behaviour changes to be aware of
- **Config geometry scaling changed.** In 2.0.3, `scale_from_stu` does *not* scale
  the voxel, so a config written for current (SI metres) places the voxel 1000×
  too small on 2.0.3 → `"Cannot initialize intra-axonal walkers"` abort. Old 2.0.3
  configs likewise need their voxels divided by 1000 if reused on current. (This is
  why the 2.0.3 run needed the ×1000 voxel.)
- **Gigantic-mesh OOM is NOT new in 2.1** — it predates it (the double-AABB-copy is
  already in 2.0.3). The 2.1.2 `fix(mem)` is what removes it; so "v2.0 could load
  it" does not hold for 2.0.3.beta specifically (it has the same bug).
- Current aborts cleanly on bad configs where 2.0.3 might silently continue
  (asserts compiled out under release).

## Caveats
- Single run per version (no averaging); flags identical; same machine.
- Not bit-exact (different RNG); this is a timing/capability comparison only.
- 2.0.3 required the ×1000 voxel to match internal geometry; both completed with
  intra placement, so the per-step collision + synthesis workload is equivalent.
- The gigantic-mesh capability row for 2.0.3 is from code inspection (same memory
  bug), not a run (skipped by request); the 2.1.2 figure (47.4 GB) is measured.
