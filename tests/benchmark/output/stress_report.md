# MC-DC stress test — GPU vs CPU (merged strands mesh, volumetric DWI)

**Substrate:** `strands_merged.ply` — one merged mesh, **48,873,064 triangles / 24,458,498 vertices** (886 MB).
**Protocol:** PGSE, **T = 1000** steps, **270 acquisitions**, duration 50 ms, D = 0.6e-9 m²/s.
**Volumetric output:** `subdivisions_number 10` → **10×10×10 = 1000 sub-signals**, each of length 270.
**Same voxel** (±~0.5 mm) and same merged mesh for both engines. Uniform seeding in the central ±0.2 mm.
**Hardware:** 1× NVIDIA RTX 3080 (10 GB) for GPU; 10 CPU cores for CPU. Host 62 GB.

## Headline

| | GPU (100M particles) | CPU (10M particles, 10 cores) |
|---|---|---|
| Wall-clock | **1:34:07** (5,647 s) | **1:13:22** (4,402 s) |
| Throughput | **17,710 walkers/s** | **2,272 walkers/s** |
| Peak device VRAM | **8.17 GB** (chunked) | — |
| Peak process RAM | 6.9 GB | **22.6 GB** |
| CPU utilisation | ~1 host core + GPU | 911% (~9.1 of 10 cores) |
| Dropped (stuck) | 0 | 0 |

**Throughput gain: ≈ 7.8× (1 RTX 3080 vs 10 CPU cores), ≈ 71× vs a single CPU core.**
Both ran the identical 48.9M-triangle mesh, T=1000, 270 acquisitions, 1000 subdivisions.

## Memory — the point of the chunking work

The GPU walks particles by materialising their full trajectory on-device
(`N·(T+1)·3` floats + `t2sum`). Un-chunked, **100M × T=1000 would need ≈ 1.60 TB of VRAM** —
impossible on a 10 GB card. The new **walker-chunking** sizes each chunk from the actually-free
VRAM (80% headroom, rounded to whole 256-thread blocks):

```
[GPU][mesh] N=100000000, T=1000: free=7.90/10.36 GB, 16020 B/walker
            -> chunk=394240 walkers x 254 chunk(s), block=256.
```

→ peak VRAM **8.17 GB, flat in N** (100M or 1B would use the same). The grid is built once and
stays resident; only the per-chunk trajectory is allocated/reused. Chunk size is **result-identical**
to a single chunk (seeding is by global walker id) — verified: single-chunk vs 5-chunk differ by
rel **1.95e-8** (float atomic-ordering only).

CPU memory (22.6 GB) is dominated by the single merged-mesh CSR grid + the 886 MB mesh; it is also
flat in N (particles are walked sequentially per thread).

## Correctness

- **Volumetric signal GPU vs CPU:** per-subdivision population fraction matches to **|Δ| ≤ 4e-4**;
  both populate **216 / 1000** subdivisions (the rest are outside the central seeding region).
  GPU sub-sum 97.2M ≈ CPU sub-sum × 10 (98.1M). ✓
- **Bug found & fixed by this test — float32 accumulator saturation.** At 100M the *global*
  DWI bin sums ~1 per walker over 1e8 walkers (~1e8); a float32 accumulator (24-bit mantissa,
  saturates near 2²⁴ ≈ 1.7e7) capped it at 33.5M instead of the true ~97.2M. The per-subdivision
  bins (~450k each) stayed accurate, which is why the *volume* was correct but the *global scalar*
  was not. Fixed by accumulating the on-device DWI in **double** (`atomicAdd(double*)`, sm_60+).
  Verified at N=20M: global b0 = 19.58M (≈0.98·N ✓), global-vs-subsum invariant **1.5e-6**
  (was ~190% before the fix).
  ⚠ The on-disk `stress_gpu_DWI.txt` (global vector) is from the pre-fix run and is float-saturated;
  **`stress_gpu_voxels_DWI.txt` (the volume) is correct.** Re-run the GPU config to regenerate a
  correct global vector if needed (the fix is now in the binary).

## Files (all under `tests/benchmark/output/`)

- `stress_gpu_voxels_DWI.txt`, `stress_cpu_voxels_DWI.txt` — volumetric signal: header `1000` then
  `270`, then 1000×270 values → reshape **(10,10,10,270)**.
- `stress_{gpu,cpu}_volFractions.txt` — per-subdivision density / intra / extra (n³×3).
- `stress_{gpu,cpu}_DWI.txt` — global signal (270). *(GPU one float-saturated, see above.)*
- `stress_{gpu,cpu}_simulation_info.txt` — params, timing, sentinel counts.
- `metrics/` — raw logs: `stress_*.log` (incl. `/usr/bin/time -v`), `vram_*.csv` (VRAM/RAM time series).

### Visualise
```
python3 tests/benchmark/voxels_to_nifti.py \
    tests/benchmark/output/stress_gpu_voxels_DWI.txt 0.1 tests/benchmark/output/stress_gpu
# -> stress_gpu_volume.nii.gz (10x10x10x270) and stress_gpu_b0.nii.gz (10x10x10)
```
Open in fsleyes / mrview / ITK-SNAP. `volFractions.txt` can be turned into a density/fraction NIfTI
with `src/volumeVolume.py` (n=10).
