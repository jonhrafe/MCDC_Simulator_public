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

---

## Addendum — GPU scaling stress (mirror BC, full voxel, 1,000,000 subdivisions)

Pure-GPU runs on the same merged 48.9M-triangle mesh, **reflective (mirror) BC**, **full-voxel**
seeding, **`subdivisions_number 100` → 100³ = 1,000,000 sub-signals** (volume 100×100×100×270).

| | 10M particles | **500M particles** |
|---|---|---|
| GPU compute | 79 s | **4,048 s (67 min)** |
| Total wall | 9:14 | **1:15:17** |
| Throughput | 126k walkers/s | **123.5k walkers/s** |
| **Peak VRAM** | **9.19 GB** | **9.19 GB** (flat in N) |
| Peak host RSS | 15.9 GB | 15.9 GB |
| Chunks | 178,688 × 56 | 178,688 × **2,799** |
| Output `_voxels_DWI.txt` | 2.0 GB | 2.0 GB (N-independent) |
| Correctness | global==sub-sum, 100% populated, no NaN | **global b0 = 4.798e8 == sub-sum (4.798e8)** |

**Key results**
- **VRAM is flat in N**: 500M particles use the *same* 9.19 GB as 10M — chunking bounds the
  trajectory, and the 4.3 GB DWI accumulator + ~2 GB grid are N-independent. The chunker sized
  178,688 walkers/chunk (2,799 chunks) to fit alongside the 1M-subdivision accumulator.
- **O(1) subdivision binning** (regular grid) makes 1,000,000 sub-voxels practical — a per-walker
  index computation instead of an O(n_sub) scan (which would be ~1e6 iters × 5e8 walkers).
- **Double accumulator is exact at 5e8**: global b0 = 4.798e8 matches the sum over 1M subdivisions
  to the digit (float32 would have saturated at ~1.7e7).
- **Live GPU ETA** (extrapolated from completed chunks) printed throughout, e.g.
  `[GPU][mesh] 1% | chunk 28/2799 | ETA 3920s (65m) | 126285 walkers/s`.
- Full-voxel seeding over the sparse strand mesh runs ~7× faster per walker than the earlier
  dense-middle seeding (fewer triangle collisions), so 500M finished in ~75 min, not hours.

**I/O note:** GPU compute was 67 min but writing the 2 GB text `_voxels_DWI.txt` (270M values,
single-threaded) adds ~8 min and is N-independent — the main remaining wall-time tax at 1M
subdivisions. A binary output path would largely remove it.

Caveat: mirror BC on the GPU is not validated against the CPU (only periodic is), so these are a
*systems* stress (memory/throughput/stability), not a signal-accuracy reference.
