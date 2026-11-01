#!/usr/bin/env python3
"""
Generate small closed watertight PLY meshes for the bounding-box tests (numpy only).

Writes two shells that mirror ``tests/accuracy/meshes/unitMesh.ply``'s format
(binary little-endian; ``property float x/y/z`` + ``property list uchar uint
vertex_indices``), so the simulator ingests them exactly like the shipped mesh:

  * a CONVEX icosphere  (subdivided icosahedron)      -> convex_icosphere.ply
  * a NON-CONVEX dumbbell (icosphere pinched at a waist) -> concave_dumbbell.ply

The dumbbell is a smooth surface of revolution (no sharp reflex edges or sliver
triangles), so the simulator walks it as fast/stably as the icosphere while still
being genuinely non-convex: its silhouette radius is non-monotonic in z (a waist),
which is exactly the reentrant geometry that trips up the CPU nearest-normal
inside/outside test but not the GPU parity ray-cast.

Both are centred on the origin (like unitMesh), consistently outward-oriented
(positive signed volume), and sized in "file units" so that with the SI-default
``ply_scale 1`` (metres/file-unit, silently x1000 to the internal mm) they span a
few tens of micrometres -- the same order as unitMesh.

The reader / geometry helpers here are reused by
``tests/validation/test_bounding_box.py`` (single PLY implementation), which needs
the mesh AABB (to check the auto-computed voxel) and the volume/AABB ratio (the
geometric intra-fraction the b0 signal must reproduce for a convex shell).

    python3 tests/tools/make_meshes.py                     # write both defaults
    python3 tests/tools/make_meshes.py --convex-out a.ply --concave-out b.ply
"""
import argparse
import os
import struct
import sys

import numpy as np

# Default file-unit sizes. With SI ``ply_scale 1`` (x1000 -> mm) these give a
# ~20 um icosphere diameter and a similarly-sized dumbbell -- comparable to unitMesh.
DEF_SPHERE_R = 1.0e-5      # icosphere radius            (file units)
DEF_SPHERE_SUB = 3         # icosahedron subdivisions (20*4^sub faces)
DEF_DB_R = 1.0e-5          # dumbbell pole radius        (file units)
DEF_DB_SUB = 3             # icosahedron subdivisions for the dumbbell
DEF_DB_WAIST = 0.30        # waist radius / pole radius (<1 pinches -> non-convex)


# --------------------------------------------------------------------------- #
#  Geometry builders                                                          #
# --------------------------------------------------------------------------- #
def icosphere(radius=DEF_SPHERE_R, subdivisions=DEF_SPHERE_SUB):
    """A CONVEX icosphere: subdivide an icosahedron and project onto the sphere.

    Returns (verts Nx3, faces Mx3) with outward-facing, consistently wound
    triangles (20 * 4**subdivisions faces).
    """
    t = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array([
        (-1,  t,  0), (1,  t,  0), (-1, -t,  0), (1, -t,  0),
        (0, -1,  t), (0,  1,  t), (0, -1, -t), (0,  1, -t),
        (t,  0, -1), (t,  0,  1), (-t,  0, -1), (-t,  0,  1),
    ], dtype=np.float64)
    faces = np.array([
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ], dtype=np.int64)

    # Loop subdivision (split each triangle into 4), caching shared midpoints so
    # the mesh stays watertight (one vertex per shared edge).
    for _ in range(int(subdivisions)):
        vlist = [v for v in verts]
        midcache = {}

        def midpoint(a, b):
            key = (a, b) if a < b else (b, a)
            idx = midcache.get(key)
            if idx is None:
                m = (verts[a] + verts[b]) * 0.5
                idx = len(vlist)
                vlist.append(m)
                midcache[key] = idx
            return idx

        new_faces = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        verts = np.asarray(vlist, dtype=np.float64)
        faces = np.asarray(new_faces, dtype=np.int64)

    # Project every vertex onto the sphere of the requested radius.
    verts = verts / np.linalg.norm(verts, axis=1, keepdims=True) * float(radius)
    return _orient_outward(verts, faces)


def dumbbell(radius=DEF_DB_R, subdivisions=DEF_DB_SUB, waist=DEF_DB_WAIST):
    """A NON-CONVEX dumbbell: an icosphere with its radius pinched at the equator.

    Reuses the icosphere topology (so it stays watertight and well-triangulated,
    no slivers) and rescales each vertex radially by ``waist + (1-waist)*uz**2``:
    the radius is ``radius`` at the poles (uz=+-1) and only ``waist*radius`` at the
    equator (uz=0). With ``waist < 1`` the silhouette radius is non-monotonic in z
    -- a genuine reentrant waist -- so the body is non-convex. Centred on origin.
    Returns (verts Nx3, faces Mx3), outward-oriented.
    """
    u, faces = icosphere(1.0, subdivisions)      # unit sphere (|u| == 1)
    uz = u[:, 2]
    r = float(radius) * (float(waist) + (1.0 - float(waist)) * uz ** 2)
    verts = u * r[:, None]
    return _orient_outward(verts, faces)


def _orient_outward(verts, faces):
    """Flip all windings if the signed volume is negative (so normals point out)."""
    if signed_volume(verts, faces) < 0.0:
        faces = faces[:, ::-1].copy()
    return verts, faces


# --------------------------------------------------------------------------- #
#  Geometry measures (reused by the test)                                     #
# --------------------------------------------------------------------------- #
def signed_volume(verts, faces):
    """Signed volume of a closed triangle mesh (divergence theorem)."""
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(np.einsum("ij,ij->i", v0, np.cross(v1, v2)).sum() / 6.0)


def aabb(verts):
    """Axis-aligned bounding box -> (min xyz, max xyz) as float64 arrays."""
    return verts.min(axis=0), verts.max(axis=0)


def volume_fraction(verts, faces):
    """|signed volume| / AABB volume -- the geometric intra-fraction."""
    lo, hi = aabb(verts)
    return abs(signed_volume(verts, faces)) / float(np.prod(hi - lo))


# --------------------------------------------------------------------------- #
#  Binary little-endian PLY I/O (mirrors unitMesh.ply)                         #
# --------------------------------------------------------------------------- #
def write_ply_binary(path, verts, faces, comment="generated by make_meshes.py"):
    """Write a binary little-endian PLY matching unitMesh.ply's property layout."""
    verts = np.ascontiguousarray(verts, dtype="<f4")
    faces = np.ascontiguousarray(faces, dtype="<u4")
    header = (
        "ply\n"
        "format binary_little_endian 1.0\n"
        f"comment {comment}\n"
        f"element vertex {len(verts)}\n"
        "property float x\n"
        "property float y\n"
        "property float z\n"
        f"element face {len(faces)}\n"
        "property list uchar uint vertex_indices\n"
        "end_header\n"
    )
    with open(path, "wb") as fh:
        fh.write(header.encode("ascii"))
        fh.write(verts.tobytes())
        three = struct.pack("<B", 3)
        for tri in faces:
            fh.write(three)
            fh.write(tri.tobytes())


def read_ply(path):
    """Read a binary/ascii PLY (float x/y/z + uchar/uint face list) -> (verts, faces).

    Minimal reader for the meshes this tool writes (and unitMesh-style files).
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    hdr_end = raw.index(b"end_header\n") + len(b"end_header\n")
    header = raw[:hdr_end].decode("ascii", "replace")
    fmt = "ascii"
    nv = nf = 0
    for line in header.splitlines():
        if line.startswith("format"):
            fmt = "binary" if "binary" in line else "ascii"
        elif line.startswith("element vertex"):
            nv = int(line.split()[2])
        elif line.startswith("element face"):
            nf = int(line.split()[2])

    if fmt == "binary":
        verts = np.frombuffer(raw, dtype="<f4", count=nv * 3,
                              offset=hdr_end).reshape(nv, 3).astype(np.float64)
        p = hdr_end + nv * 3 * 4
        faces = np.empty((nf, 3), dtype=np.int64)
        for i in range(nf):
            cnt = raw[p]; p += 1
            faces[i] = struct.unpack_from("<%dI" % cnt, raw, p)[:3]
            p += 4 * cnt
    else:
        toks = raw[hdr_end:].split()
        verts = np.array(toks[:nv * 3], dtype=np.float64).reshape(nv, 3)
        p = nv * 3
        faces = np.empty((nf, 3), dtype=np.int64)
        for i in range(nf):
            cnt = int(toks[p]); p += 1
            faces[i] = [int(x) for x in toks[p:p + 3]]
            p += cnt
    return verts, faces


# --------------------------------------------------------------------------- #
#  CLI                                                                         #
# --------------------------------------------------------------------------- #
def _report(name, verts, faces):
    lo, hi = aabb(verts)
    print(f"  {name}: {len(verts)} verts, {len(faces)} tris  "
          f"AABB min=({lo[0]:.3e},{lo[1]:.3e},{lo[2]:.3e}) "
          f"max=({hi[0]:.3e},{hi[1]:.3e},{hi[2]:.3e})  "
          f"vol_frac={volume_fraction(verts, faces):.4f}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    default_mesh_dir = os.path.join(os.path.dirname(here), "accuracy", "meshes")
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--convex-out",
                    default=os.path.join(default_mesh_dir, "convex_icosphere.ply"))
    ap.add_argument("--concave-out",
                    default=os.path.join(default_mesh_dir, "concave_dumbbell.ply"))
    ap.add_argument("--sphere-radius", type=float, default=DEF_SPHERE_R)
    ap.add_argument("--sphere-subdivisions", type=int, default=DEF_SPHERE_SUB)
    ap.add_argument("--dumbbell-radius", type=float, default=DEF_DB_R)
    ap.add_argument("--dumbbell-subdivisions", type=int, default=DEF_DB_SUB)
    ap.add_argument("--dumbbell-waist", type=float, default=DEF_DB_WAIST)
    a = ap.parse_args()

    cv, cf = icosphere(a.sphere_radius, a.sphere_subdivisions)
    dv, df = dumbbell(a.dumbbell_radius, a.dumbbell_subdivisions, a.dumbbell_waist)

    write_ply_binary(a.convex_out, cv, cf,
                     comment="convex icosphere (make_meshes.py)")
    write_ply_binary(a.concave_out, dv, df,
                     comment="non-convex dumbbell (make_meshes.py)")

    print("Wrote meshes:")
    _report(os.path.relpath(a.convex_out), cv, cf)
    _report(os.path.relpath(a.concave_out), dv, df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
