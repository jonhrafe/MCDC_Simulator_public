#!/usr/bin/env python3
"""
Import an external PLY into the test suite: strip to xyz + triangles (binary LE).

Real meshes (VTK exports, the astrocyte "voronoi" cells, ...) carry per-vertex
NORMALS and RGB COLOUR. The simulator's PLY reader expects the ``unitMesh.ply``
layout -- vertices that are ONLY ``property float x/y/z`` -- so those extra vertex
properties misalign the byte cursor and the face section is read as garbage
("Non-triangle face detected" -> segfault). This tool copies a source PLY into a
clean mesh with ONLY ``property float x/y/z`` vertices and ``property list uchar
uint vertex_indices`` triangle faces (exactly unitMesh's format), dropping
colour/normals; the geometry is preserved and any non-triangle face is
fan-triangulated.

    python3 tests/tools/import_mesh.py SRC.ply DST.ply
    python3 tests/tools/import_mesh.py --dir SRC_DIR --out-dir DST_DIR   # batch *.ply
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_meshes import write_ply_binary  # noqa: E402

# PLY property type -> (numpy dtype string, byte size)
_NP = {
    "char": ("i1", 1), "int8": ("i1", 1), "uchar": ("u1", 1), "uint8": ("u1", 1),
    "short": ("<i2", 2), "int16": ("<i2", 2), "ushort": ("<u2", 2), "uint16": ("<u2", 2),
    "int": ("<i4", 4), "int32": ("<i4", 4), "uint": ("<u4", 4), "uint32": ("<u4", 4),
    "float": ("<f4", 4), "float32": ("<f4", 4), "double": ("<f8", 8), "float64": ("<f8", 8),
}


def read_ply_any(path):
    """Read a binary/ascii PLY with ARBITRARY vertex properties -> (verts, tris).

    Returns ``verts`` (M x 3 float64, the x/y/z columns only) and ``tris`` (T x 3
    int64), fan-triangulating any polygon face. Extra vertex properties (normals,
    colour) and any non-position columns are skipped.
    """
    raw = open(path, "rb").read()
    hdr_end = raw.index(b"end_header\n") + len(b"end_header\n")
    header = raw[:hdr_end].decode("ascii", "replace")

    fmt = "ascii"
    nv = nf = 0
    vprops = []           # [(name, np_dtype, size)] in order, for the vertex element
    face_count_t = face_idx_t = None
    cur = None
    for line in header.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "format":
            fmt = "binary" if "binary" in line else "ascii"
        elif t[0] == "element":
            cur = t[1]
            if t[1] == "vertex":
                nv = int(t[2])
            elif t[1] == "face":
                nf = int(t[2])
        elif t[0] == "property":
            if cur == "vertex":
                np_dt, sz = _NP[t[1]]
                vprops.append((t[2], np_dt, sz))
            elif cur == "face" and t[1] == "list":
                face_count_t = _NP[t[2]]      # (dtype, size) of the count
                face_idx_t = _NP[t[3]]        # (dtype, size) of each index

    if fmt != "binary":
        raise ValueError(f"{path}: only binary_little_endian PLY is supported here")
    if face_count_t is None or face_idx_t is None:
        raise ValueError(f"{path}: no 'property list' face definition in header")

    # vertices: build a structured dtype spanning ALL properties, keep x/y/z
    rec = np.dtype([(n, d) for n, d, _ in vprops])
    varr = np.frombuffer(raw, dtype=rec, count=nv, offset=hdr_end)
    verts = np.stack([varr["x"], varr["y"], varr["z"]], axis=1).astype(np.float64)

    # faces: variable-length list -> fan-triangulate
    p = hdr_end + nv * rec.itemsize
    csz = face_count_t[1]
    isz = face_idx_t[1]
    tris = []
    for _ in range(nf):
        cnt = int(np.frombuffer(raw, dtype=face_count_t[0], count=1, offset=p)[0])
        p += csz
        idx = np.frombuffer(raw, dtype=face_idx_t[0], count=cnt, offset=p).astype(np.int64)
        p += isz * cnt
        for k in range(1, cnt - 1):           # fan: (0, k, k+1)
            tris.append((idx[0], idx[k], idx[k + 1]))
    return verts, np.asarray(tris, dtype=np.int64)


def convert(src, dst):
    verts, tris = read_ply_any(src)
    write_ply_binary(dst, verts, tris, comment=f"stripped xyz+tris from {os.path.basename(src)}")
    lo, hi = verts.min(0), verts.max(0)
    print(f"  {os.path.basename(src)} -> {os.path.relpath(dst)}  "
          f"{len(verts)} verts, {len(tris)} tris  "
          f"span=({hi[0]-lo[0]:.3g},{hi[1]-lo[1]:.3g},{hi[2]-lo[2]:.3g})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", nargs="?", help="source .ply")
    ap.add_argument("dst", nargs="?", help="destination .ply")
    ap.add_argument("--dir", help="batch: convert every *.ply in this directory")
    ap.add_argument("--out-dir", help="batch: destination directory")
    a = ap.parse_args()

    if a.dir:
        if not a.out_dir:
            print("--dir requires --out-dir", file=sys.stderr)
            return 2
        os.makedirs(a.out_dir, exist_ok=True)
        for src in sorted(glob.glob(os.path.join(a.dir, "*.ply"))):
            convert(src, os.path.join(a.out_dir, os.path.basename(src)))
        return 0
    if not (a.src and a.dst):
        ap.print_help()
        return 2
    os.makedirs(os.path.dirname(os.path.abspath(a.dst)), exist_ok=True)
    convert(a.src, a.dst)
    return 0


if __name__ == "__main__":
    sys.exit(main())
