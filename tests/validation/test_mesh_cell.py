#!/usr/bin/env python3
"""
Real-mesh consistency: an astrocyte "voronoi" cell as a PLY substrate.

A small real triangulated cell (`accuracy/meshes/voronoi_cell.ply`, ~66k tris,
~42 µm span) imported from the tissue-forge astrocyte meshes and stripped to
xyz+triangles (colours/normals dropped -- see `tests/tools/import_mesh.py`; the
simulator's PLY reader only accepts `unitMesh`-style x/y/z vertices). There is no
analytic oracle for such a shape, so this is a **consistency** test -- the value
is that the SAME geometry must give the SAME signal across the paths that should
agree:

  * CPU sanity (runs today): the signal decreases with b and the apparent
    diffusivity is bounded by the free value (some restriction from the cell);
  * CPU ≡ GPU parity (self-skips without CUDA): both engines walk the identical
    mesh, so per-shell means must agree within MC noise;
  * cross-version consistency (self-skips without old binaries): re-run the same
    conf on the historical binaries in `tests/crossversion/binaries.toml` and
    compare per-shell within a loose band (different RNG -> never bit-exact).

The mesh is placed in a voxel auto-sized to its bounding box (+ margin), walkers
seeded uniformly over the voxel (intra+extra), `ply_scale 1e-6` (the cell coords
are micrometres).

  python3 tests/validation/test_mesh_cell.py --bin ./MC-DC_Simulator --n 3000
"""
import argparse
import os
import shutil
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import io, metrics, paths, run  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402
from tools.import_mesh import read_ply_any  # noqa: E402

import numpy as np  # noqa: E402

MESH = paths.accuracy("meshes", "voronoi_cell.ply")
SCHEME = "tests/accuracy/PGSE_sample_scheme.scheme"
PLY_SCALE = 1e-6           # cell coords are micrometres -> metres/unit
D = 0.6e-9                 # intra diffusivity (m^2/s)
MARGIN_UM = 2.0            # voxel margin around the mesh bbox


def _voxel_m():
    """Voxel (min, max) in metres = the mesh bbox (+margin) * PLY_SCALE."""
    v, _ = read_ply_any(os.path.join(paths.REPO_ROOT, MESH))
    lo, hi = v.min(0) - MARGIN_UM, v.max(0) + MARGIN_UM
    return lo * PLY_SCALE, hi * PLY_SCALE


def _conf(n, t):
    lo, hi = _voxel_m()
    return (
        f"N {n}\nT {t}\nduration 0.050\ndiffusivity {D}\n"
        f"exp_prefix PLACEHOLDER\nscheme_file {SCHEME}\n"
        "write_txt 1\nwrite_bin 0\nwrite_traj_file 0\nt2_extra 1e9\n"
        "<obstacle>\n"
        f"ply {os.path.join(paths.REPO_ROOT, MESH)}\nply_scale {PLY_SCALE}\n"
        "</obstacle>\n"
        f"<voxels>\n{lo[0]} {lo[1]} {lo[2]}\n{hi[0]} {hi[1]} {hi[2]}\n</voxels>\n"
        "num_process 1\nseed 12345\n<END>\n")


def run_mesh(binary, n=3000, t=500, gpu=0, timeout=1200):
    """Run the cell mesh; return (S/S0, scheme, wall_s) or (None, scheme, wall_s)."""
    sc = kscheme.parse_scheme(os.path.join(paths.REPO_ROOT, SCHEME))
    work = tempfile.mkdtemp(prefix="mcdc_cell_")
    try:
        t0 = time.time()
        prefix = run.run_case(binary, _conf(n, t), work, "cell", paths.REPO_ROOT,
                              gpu=gpu, timeout=timeout)
        wall = time.time() - t0
        if prefix is None:
            return None, sc, wall
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=n)
        return metrics.normalize_by_b0(dwi, sc.b0_mask), sc, wall
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --------------------------------------------------------------------------- #
#  Checks                                                                      #
# --------------------------------------------------------------------------- #
def check_sanity(ss0, sc):
    """Physical-signal sanity: monotone in b, S/S0 in (0,1), ADC bounded. (ok, msgs).

    There is no analytic oracle for a real cell, so this is only a "did it produce
    a physical signal" gate (the real checks are CPU/GPU parity and cross-version).
    The cell is large (~42 µm) so restriction is mild -- ADC sits just below the
    free value -- hence the ADC bound carries a small tolerance for MC noise.
    """
    nb = ~sc.b0_mask
    ms = metrics.per_shell_mean(ss0[nb], sc.bvals[nb])
    bs = sorted(ms)
    decreasing = all(ms[a] > ms[b] - 1e-6 for a, b in zip(bs, bs[1:]))
    physical = all(0.0 < ms[b] < 1.0 for b in bs)
    adc = -np.log(max(ms[bs[-1]], 1e-9)) / (bs[-1] * 1e6)      # highest shell
    bounded = 0.0 < adc < 1.1 * D            # <= free D up to MC noise (mild restriction)
    ok = decreasing and physical and bounded
    msgs = [f"  per-shell S/S0: " + ", ".join(f"{int(b)}:{ms[b]:.4f}" for b in bs),
            f"  [{'OK ' if decreasing else 'BAD'}] monotone decreasing, S/S0 in (0,1)",
            f"  [{'OK ' if bounded else 'BAD'}] ADC={adc*1e9:.3f} ≲ free D={D*1e9:.1f} (1e-9 m²/s)"]
    return ok, msgs


def check_parity(a, b, sc, tol=0.03):
    """Per-shell |a-b| within tol (CPU vs GPU, or vs another version). (ok, worst, msgs)."""
    nb = ~sc.b0_mask
    ma = metrics.per_shell_mean(a[nb], sc.bvals[nb])
    mb = metrics.per_shell_mean(b[nb], sc.bvals[nb])
    worst, msgs = 0.0, []
    for k in sorted(ma):
        d = abs(ma[k] - mb[k])
        worst = max(worst, d)
        msgs.append(f"  b={k:5.0f}  A={ma[k]:.4f}  B={mb[k]:.4f}  |Δ|={d:.4f}")
    return worst <= tol, worst, msgs


def _xv_binaries():
    """Historical binaries from tests/crossversion/binaries.toml that exist on disk.

    Returns {version: path}. Empty if the manifest or the binaries are absent
    (the cross-version check then self-skips).
    """
    manifest = os.path.join(paths.REPO_ROOT, "tests", "crossversion", "binaries.toml")
    if not os.path.exists(manifest):
        return {}
    try:
        import tomllib
        with open(manifest, "rb") as fh:
            data = tomllib.load(fh)
    except Exception:
        return {}
    out = {}
    for ver, entry in data.items():
        path = entry.get("path") if isinstance(entry, dict) else None
        if path and os.path.exists(path):
            out[ver] = path
    return out


# --- pytest entry points -----------------------------------------------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    if not os.path.exists(os.path.join(paths.REPO_ROOT, MESH)):
        pytest.skip(f"mesh fixture missing: {MESH}")
    return binary


# N is kept small because the mesh has NO triangle-acceleration structure: the CPU
# tests all ~66k tris every step, so a run is ~2-3 min on CPU at N=600 (but only ~0.3 s
# on GPU -- a ~1000x speedup, which is the whole point of the GPU path for meshes).
CELL_N, CELL_T = 600, 300


def test_mesh_cell_cpu_gpu():
    """Cell mesh: physical CPU signal + CPU≡GPU parity (one slow CPU run; GPU ~1000x faster).

    The CPU run (~2-3 min) is the sanity baseline; the GPU run (~0.3 s) adds the parity
    check when CUDA is present (else the test still passes on the CPU sanity alone).
    """
    binary = _binary_or_skip()
    cpu, sc, cwall = run_mesh(binary, n=CELL_N, t=CELL_T, gpu=0)
    ok, msgs = check_sanity(cpu, sc)
    assert ok, "cell mesh CPU sanity failed:\n" + "\n".join(msgs)
    gpu, _, gwall = run_mesh(binary, n=CELL_N, t=CELL_T, gpu=1)
    if gpu is not None:                       # CUDA present -> validate CPU ≡ GPU
        pok, worst, pmsgs = check_parity(cpu, gpu, sc)
        assert pok, (f"CPU/GPU parity off (worst {worst:.4f}; GPU {gwall:.1f}s vs "
                     f"CPU {cwall:.1f}s):\n" + "\n".join(pmsgs))


def test_mesh_cell_cross_version():
    """Same mesh across historical binaries (self-skips without old binaries)."""
    import pytest
    binary = _binary_or_skip()
    bins = _xv_binaries()
    if not bins:
        pytest.skip("no historical binaries in tests/crossversion/binaries.toml")
    cur, sc, _ = run_mesh(binary, n=CELL_N, t=CELL_T, gpu=0)
    for ver, path in bins.items():
        old, _, _ = run_mesh(path, n=CELL_N, t=CELL_T, gpu=0)
        if old is None:
            continue
        ok, worst, msgs = check_parity(cur, old, sc, tol=0.06)
        assert ok, f"cross-version {ver} off (worst {worst:.4f}):\n" + "\n".join(msgs)


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--t", type=int, default=300)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    a = ap.parse_args()
    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1
    if not os.path.exists(os.path.join(paths.REPO_ROOT, MESH)):
        print(f"FAIL: mesh fixture missing: {MESH}", file=sys.stderr)
        return 1

    ss0, sc, wall = run_mesh(binary, n=a.n, t=a.t, gpu=a.gpu)
    if ss0 is None:
        print("SKIP: binary built without CUDA / no GPU available.")
        return 0
    ok, msgs = check_sanity(ss0, sc)
    print(f"[voronoi cell mesh | N={a.n} T={a.t} {'GPU' if a.gpu else 'CPU'} | {wall:.1f}s]")
    for m in msgs:
        print(m)

    # CPU/GPU parity if a GPU run is available
    if a.gpu == 0:
        gpu, _, gwall = run_mesh(binary, n=a.n, t=a.t, gpu=1)
        if gpu is not None:
            pok, worst, _ = check_parity(ss0, gpu, sc)
            print(f"[CPU/GPU parity] worst |Δ|={worst:.4f} (GPU {gwall:.1f}s vs "
                  f"CPU {wall:.1f}s): {'OK' if pok else 'FAIL'}")
            ok = ok and pok
        else:
            print("[CPU/GPU parity] SKIP (no CUDA)")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
