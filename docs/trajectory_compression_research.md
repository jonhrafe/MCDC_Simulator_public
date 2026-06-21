# Trajectory (`.traj`) size reduction — design note

Status: **research / design only** (not scheduled for implementation). Goal: shrink
`.traj` files (can be hundreds of GB) for storage, visualization, and re-synthesis
of the signal with different parameters, with **no measurable effect on the DWI**.

## 1. Current format (the contract)

- MCDC writes positions as **binary float32**, 3 per position
  (`trajectory.cpp:writePositionBinary`), preceded by a small **text header**
  (`dyn_duration`, `N`, `T` / number of saved steps). Internally positions are
  `double`; they are already down-cast to `float` on write.
- CUDASynth (`CUDASynth.cu`) reads the `.traj` as a flat `float` array and
  **`seekg`s to a per-walker block** (`new float[blocks*block_size]`,
  `in_traj.seekg(file_pos*block_size)`), then recomputes the DWI on the GPU.
- Layout per walker: `(num_saved_steps) × 3` float32, contiguous. Random access is
  **by walker** (block), sequential within a walker.
- Size ≈ `N · S · 3 · 4` bytes (`S` = saved steps). E.g. `N=1e4, S=5e3` → ~600 GB;
  `N=1e5, S=1e4` → ~12 TB. Matches the "hundreds of GB" observation.

**Hard constraint (from the consumer):** any new format must support random access
*per walker* — so whole-file gzip/zstd (no seek) is out. Compression must be
**per-block (per-walker) with an offset index**.

## 1a. Chosen architecture (decided)

A **standalone, offline tool** — NOT part of the simulation, never run during it:

    [ MCDC sim ] --float32 .traj-->  (unchanged)
                                       │
                          [ trajzip CLI ]  reads a fully-computed .traj + header,
                                       │   compresses (delta+quantize+zstd, per
                                       │   walker, + offset index)
                                       ▼
                              archive  .traj.zc  on disk   (the space saving)
                                       │
                          [ CUDASynth ]  detects .traj.zc, reads the index, per
                                       │   walker: seek → read block → decompress
                                       ▼   → reconstruct float32 → existing path

Properties this buys us:
- **The simulator is untouched** (no per-step cost, no risk, no new deps in MCDC).
- Compression is a deliberate **archival** step for storing many runs cheaply.
- **CUDASynth owns decompression** and reads the compressed file directly — no temp
  uncompressed `.traj` on disk, no whole-file decompress pre-pass.
- The format is **self-describing** (header carries `N`, `S`, quantization grid,
  per-walker offset index) so the tool and CUDASynth only need to share the spec,
  not side metadata.

## 2. How much position precision does the signal actually need?

The DWI is `S = (1/N) Σ_k exp(i φ_k)`, with phase
`φ_k = γ Σ_t G(t)·x_k(t) Δt`. The sensitivity of phase to a position error is set by
the **q-value** `q = γ·G·δ` (δ = pulse width). A position error `Δx` perturbs the
phase by `≈ q·Δx`. To keep the per-step phase error well below the signal noise
floor (say < 0.01 rad), we need:

    Δx  <  0.01 / q

Worked numbers (γ = 2.675e8 rad/s/T):

| gradient G | δ | q = γGδ | 1/q | needed Δx (<0.01/q) |
|---|---|---|---|---|
| 0.3 T/m | 4 ms | 3.2e5 rad/m | 3.1 µm | ~31 nm |
| 1.0 T/m | 10 ms | 2.7e6 rad/m | 0.37 µm | ~3.7 nm |

So **~3–30 nm absolute position precision is plenty** for the signal, across
realistic gradients.

### Key insight: float32 is the *wrong* precision model
float32 keeps **relative** precision (~`x·2⁻²³` ≈ `x·1.2e-7`): ~12 nm at 100 µm,
but ~120 nm at 1 mm. So for strong gradients + large voxels, **float32 absolute
positions are already borderline-lossy** at the far edges, *and* waste bits on
fine precision near the origin. A **fixed-point integer grid** (uniform absolute
precision, e.g. 5 nm everywhere) is simultaneously **more signal-accurate and more
compressible** than float32. This is the real lever — not "float32 → float16"
(float16 ≈ 50 nm at 100 µm, too coarse at large magnitudes).

## 3. Options (smallest change → best ratio)

1. **Per-walker zstd of the float32 block + offset index.** Lossless, ~1.3–2× only
   (float32 mantissas are high-entropy). Minimal: keep producing float32, compress
   each walker block, store offsets. CUDASynth links zstd and decompresses per block.
2. **Delta encoding** (store per-step displacements, not absolute positions). A step
   is `~√(6 D Δt)` (sub-µm), tiny vs the voxel (tens–hundreds µm), so deltas have a
   far smaller dynamic range and entropy-code much better. Store walker's initial
   position at full precision; deltas thereafter.
3. **Fixed-point quantization** of positions/deltas to a signal-safe grid (e.g.
   1–5 nm). Quantized deltas fit in **int16** (±~160 µm at 5 nm) → 6 B/position vs
   12 B before, *before* entropy coding. Signal-lossless by §2.
4. **Combine (recommended target): per-walker block = [float32 x0,y0,z0] +
   [int16 quantized deltas] → zstd, with an offset index.** Estimated **~4–8×**
   overall, with no measurable DWI change. Random access preserved (seek to the
   walker's compressed block via the index, decompress, integrate deltas).

### Orthogonal lever (already half-built): store fewer steps
The signal integral only "sees" timesteps where the gradient is ON. For PGSE that's
the two δ-pulses — the **active timesteps** we already compute for perf-B4
(`PGSESequence::buildActiveTimesteps`). Saving only those steps is *lossless for a
fixed scheme* and can be ~10–25× smaller. **Caveat:** it breaks "re-synthesize with
*different* timing parameters" (different Δ/δ → different active steps). The existing
`steps_subset`/`pos_times` mechanism already lets the user pick a step subset; this
is a knob, not a default. For visualization, a coarse step subset is also fine.

## 4. Recommendation / opinion

- **Your instinct is right, with one refinement:** dropping position precision will
  not affect the signal — but do it as a **fixed-point grid (≈1–5 nm), on deltas**,
  not as a smaller float. Fixed-point is both more accurate for the signal (uniform
  absolute precision) and far more compressible. "Lossy" here is below the physical
  resolution, so it is effectively lossless *for the signal*.
- **Two components, one shared spec** (see §1a): an offline `trajzip` CLI (producer)
  and a CUDASynth read path (consumer). MCDC is not modified at all.
- **Format:** self-describing — header (`N`, `S`, quantization grid `Δq`, flags) +
  per-walker compressed blocks + an offset index. Per-walker blocks are the only
  shape compatible with CUDASynth's per-walker `seekg`.
- **Per-block codec:** `[float32 x0,y0,z0]` then int16 quantized deltas → zstd.
  Reconstruction is exact-to-the-grid (signal-lossless, §2).
- **Dependency:** zstd only (fast, excellent ratio, permissive, tiny) — added to the
  `trajzip` tool and to CUDASynth, **not** to MCDC. Delta+quantize is our own code.
- **Suggested staging when this becomes the job:**
  1. `trajzip` CLI: read `.traj` + header → write `.traj.zc` (delta+quantize+zstd
     per walker + index). Add a `--verify` that decompresses and checks max position
     error ≤ grid, and (ideally) re-runs the signal to confirm DWI matches the
     float32 baseline within noise.
  2. CUDASynth: detect `.traj.zc`, read header+index, decompress per walker block
     into the existing float buffer. Gate behind file-extension/magic so plain
     `.traj` still works.
  3. Tune the quantization grid from the real `Δt`/gradient ranges (§5).
- **Round-trip safety:** keep the original `.traj` until `trajzip --verify` passes;
  the compressed archive is only trustworthy once the decompressed signal matches.

## 5. Open questions for later
- Exact `Δt` and gradient ranges actually used → fixes the quantization grid.
- Is the `.traj` ever read by tools other than CUDASynth / the MCDC reader? (Defines
  how many readers must learn the new format.)
- Visualization path: does it tolerate a step subset / lower precision (almost
  certainly yes)?
