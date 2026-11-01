# Per-version configs (`tests/crossversion/confs/<version>/<case>.conf`)

The old binaries do **not** parse the current SI `tests/accuracy/*.conf` files
unmodified — see `tests/crossversion/README.md` and `tests/crossversion/binaries.toml`.
Each historical version needs its own config per case, kept here.

## Unit convention (v1.5, v2.0 — verified against `src/parameters.cpp` at both commits)

- `scale_from_stu 1` — `diffusivity`/`duration` are given in SI (`m^2/s`, `s`)
  and the parser auto-converts them to the engine's internal `mm^2/ms`, `ms`.
  Omit it and those two fields are assumed **already** in `mm^2/ms` / `ms`.
- `<voxel>` bounds are **always in millimetres**, regardless of `scale_from_stu`
  — that flag only touches `diffusivity` and `sim_duration`, not the voxel
  reader. A current-SI voxel of `-5e-3..5e-3` m is `-5..5` here.
- The scheme file columns (`G`, `Δ`, `δ`, `TE`) go through the **same**
  `scale_from_stu` flag (`src/scheme.cpp::readPGSE`) — when set, they're also
  read as SI (`T/m`, seconds) and auto-converted. The current SI scheme file
  (`tests/accuracy/PGSE_sample_scheme.scheme`) can therefore be **reused
  as-is** for v1.5/v2.0 confs; no separate per-version scheme file needed.
- No `gpu` line for v1.5/v2.0 — both predate the CUDA path. An unrecognized
  keyword only warns in these parsers, so leaving it in wouldn't break the run,
  but it's confusing to have it there.
- **No `#` comments.** v1.5/v2.0's `readParams` has no comment-skipping at
  all (`#` support was added later) — the reader tokenizes the whole file, and
  a comment word that happens to fuzzy-match a real keyword (`str_dist`
  allows edit distance up to 1-3 depending on the field) gets treated as that
  parameter, consuming the next token(s) as its value and silently shifting
  everything after it out of alignment. Concretely: an inline `#`-header on
  the first `free_diffusion.conf` draft here shifted `N` and blew straight
  through `checkSimulationParameters`'s `num_walkers > 1e9` guard. **Keep
  these confs comment-free** — put rationale in this README instead.

v1.5 and v2.0 share this convention (confirmed identical `scale_from_stu`
handling in `parameters.cpp` at both commits) — a new version's conf should
re-verify against its own `parameters.cpp`/`scheme.cpp` before assuming it.

## What's covered

Two cases have per-version confs today: **`free_diffusion`** and the **single
cylinder** (native list + capped PLY). A case earns a historical conf only when two
things line up:

1. **Units** (above) — mechanical, same fix for every case.
2. **Feature existence** — v1.5/v2.0 are `2016`–`2018`-era commits; plenty of current
   `tests/accuracy/*.conf` features postdate them entirely (gamma packing generators,
   per-obstacle `d_intra`/T2/permeability lists, subdivisions, mesh compartments, the
   modern permeability model). Those configs have **no** historical equivalent to run.

### Single cylinder

The scale question the earlier draft flagged is now resolved by the reference confs that
`script_comparing_3v.py` actually runs (`one_cylinder_v1p5.conf`, `one_cylinder_z_axis_1p5.txt`,
`one_cylinder_ply_v2p0.conf`):

- **Native list** (`cylinder_restricted.conf` → `single_cylinder.list`) — the list is a
  `scale` line then `x y z radius`. `scale_from_stu` does **not** touch list geometry, so
  the list is in the engine's native millimetres: scale `1`, row `0 0 0 0.005` → a 5 µm
  cylinder on z. (The current SI list is scale `1e-06`, row `0 0 0 5`; the old mm build
  wants `1` / `0.005`.) Covered for **v1.5 and v2.0**.
- **Capped PLY** (`cylinder_ply.conf`) — same mesh as the current SI conf
  (`tests/accuracy/meshes/cylinder_capped.ply`) but `ply_scale 1` (mm per file-unit = the
  current SI `0.001` × 1000) and a central `<spawning_area>` in mm. Covered for **v2.0**;
  v1.5 uses the native list only (no v1.5 PLY reference conf exists, and the native list
  already validates the same 5 µm cylinder there).
- **No `t2_extra`** in either — it cancels in the b0-normalised signal (constant TE) and
  does not exist at v1.5 at all, so it is simply omitted.
- **`permeability 0`** is the v1.5/v2.0 keyword (the current build uses
  `obstacle_permeability`); `<voxels>` bounds and the list/PLY scales are all in mm.

## v2.2

No v2.2 entry: that line is the current code (`VERSION_ID` went `2.1.1.beta` → `2.1.2.beta` →
`2.1.3.beta` → 2.3.x). If a specific commit should stand in for "v2.2", add it to `KNOWN_COMMITS` in
`tests/tools/build_old_version.py` and a `tests/crossversion/confs/v2.2/` directory following this same pattern.

## Adding a new case

1. Confirm the feature existed at that commit (`git show <commit>:src/*.cpp`).
2. Convert the current SI conf: apply `scale_from_stu 1`, convert voxel/obstacle
   geometry to mm, verify against that commit's `parameters.cpp` rather than
   assuming this file's rules hold unchanged.
3. Drop it at `tests/crossversion/confs/<version>/<case>.conf`.
4. Freeze it: `python3 tests/tools/gen_reference.py --version <version> --config
   tests/crossversion/confs/<version>/<case>.conf --seed <fixed> --num-process
   <fixed> --case <case>` (builds the binary automatically on first use).
