"""Parsers for the gamma packer's generated substrate (its list dump + voxel).

The simulator writes the realised packing to
``<prefix>_gamma_distributed_{cylinder,sphere}_list.txt`` (a native list) and the
auto-sized voxel to ``<prefix>_simulation_info.txt``. Both the generator-fidelity
test and the MC-vs-analytic packing tests need to read those back, so the parsing
lives here once. File units are µm (the generator prints internal-mm × 1e3).
"""
import numpy as np


def voxel_side_um(sim_info_path):
    """Voxel side length in µm, parsed from a ``*_simulation_info.txt``.

    The 'Voxel limits:' header is followed by two ``( x y z )`` lines (min, max)
    in internal mm; the side is ``max(hi - lo) * 1e3``.
    """
    lines = open(sim_info_path).read().splitlines()
    for i, ln in enumerate(lines):
        if "Voxel limits" in ln:
            def vec(s):
                return np.array([float(x) for x in s.strip().strip("()").split()])
            lo, hi = vec(lines[i + 1]), vec(lines[i + 2])
            return float((hi - lo).max()) * 1e3
    raise ValueError(f"no 'Voxel limits' in {sim_info_path}")


def read_list(path):
    """Read a gamma list dump -> (centres Nx3, radii N) in µm.

    Row layout: ``x y z r percolation T2``. Line 0 is the scale factor (skipped;
    values are kept in the file's µm units). Includes periodic-image replicas.
    """
    rows = [ln.split() for ln in open(path).read().splitlines()[1:] if ln.strip()]
    a = np.array(rows, dtype=float)
    return a[:, 0:3], a[:, 3]


def primary(centres, radii, voxel_um, dims):
    """Keep one representative per obstacle: centre in the half-open [0, L).

    The list carries periodic replicas shifted by exactly ±L, so every obstacle
    has EXACTLY ONE representative in any half-open interval of length L; a
    half-open [0, L) window on the packed axes (``dims``: (0,1) for z-cylinders,
    (0,1,2) for spheres) de-dups cleanly (a closed [0, L] double-counts
    face-straddling obstacles). ``eps`` absorbs float dust at the 0 boundary.
    """
    eps = 1e-6 * voxel_um
    mask = np.ones(len(centres), dtype=bool)
    for d in dims:
        mask &= (centres[:, d] >= -eps) & (centres[:, d] < voxel_um - eps)
    return centres[mask], radii[mask]
