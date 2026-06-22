# GPU PLY stress test — tissue strand model

A full tissue model of **10,924 individual PLY strands** (~24.5M vertices,
**48.9M triangles**, geometry in µm). Validates the GPU PLY pipeline at scale and
compares it against the CPU on identical geometry.

## Setup
- Global vertex bbox (all PLYs): `[-494.7, -487.5, -505.9] .. [495.0, 488.8, 507.9]` µm
  → a ~1 mm cube centred near the origin. `ply_scale 0.001` (µm → internal mm).
- Voxel = the global bbox; walkers spawn in a central dense cube (`<spawning_area>`
  ±200 µm), uniform placement.
- Crossing sentinel OFF (`deportation 0`) — the strands are **not watertight**
  (the loader reports ~360 highly irregular triangles).
- `N=10000`, `T=1000`, `duration 0.050 s`, `D=0.6e-9 m²/s`, T2 = ∞, `num_process 1`.

## Reproduce
```
# 1. merge all strand PLYs into one binary PLY (offsets vertex indices)
python3 tests/benchmark/merge_plys.py <MESH_DIR> <OUT>/strands_merged.ply
# 2. point strands_merged.list at that file (one line: "<path> 0.001")
# 3. run, toggling gpu 0 / gpu 1
./MC-DC_Simulator --conf tests/benchmark/strands_merged.conf
# (run_strands_stress.sh runs the un-merged 10,924-PLY form directly on GPU)
```

## Results

### Loading: merged binary PLY vs 10,924 ASCII PLYs (GPU)
| input | wall time |
|---|---|
| 10,924 ASCII PLYs (`ply_file_list`) | 50.6 s |
| 1 merged binary PLY (929 MB) | **10.3 s** |

Merging + binary format cuts load time ~5×; identical signal. **Critical fix:** on
GPU runs the per-PLY CPU collision grid is skipped (`PLYObstacle::build_collision_grid
= !gpu`) — otherwise one grid per PLY costs ~61 GB for just 100 PLYs (the GPU builds
its own single device grid over all triangles).

### CPU vs GPU on the identical merged model (single global grid, 1 core)
| | wall time | RAM | mean S/N |
|---|---|---|---|
| GPU (float32) | **10.3 s** | 7.0 GB | 0.9269 |
| CPU (double)  | 72 s | 23.1 GB | 0.9543 |

- GPU ~**7× faster** at **3× less RAM** vs the single-core CPU; the GPU is nearly
  flat in N, so its lead grows with particle count.
- **Signal agreement ~2.9%** (GPU slightly *less* restricted). This is far larger
  than the **0.01%** the GPU matches the CPU on a *clean, watertight* mesh
  (`mesh_restricted`/unitMesh). The gap is attributed to the **non-watertight mesh**:
  float32 vs double at irregular/degenerate triangle edges, leakage through mesh
  gaps, and stuck-particle handling (CPU discards walkers exceeding max bounces; the
  GPU caps bounces and keeps them).

**Takeaway:** the GPU PLY pipeline scales to tens of millions of triangles and is
several× faster, but **mesh quality matters** — on a non-watertight substrate the
GPU/CPU signals agree only to a few %, vs ~0.01% on clean geometry. Decided
(2026-06): accept this as expected for imperfect meshes and document it here.

> Note: intra/extra placement is not yet classified per-mesh on the GPU (it would
> need a point-in-mesh test); the GPU mesh kernel places uniformly in the spawning
> area, so use uniform placement for apples-to-apples CPU/GPU mesh comparisons.
