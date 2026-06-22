#!/usr/bin/env python3
"""
Merge a folder of (ASCII) PLY meshes into ONE binary PLY.

Concatenates all vertices and all faces, offsetting each mesh's face indices by the
running vertex count, and writes a single binary_little_endian PLY in the exact layout
the simulator's reader expects (float32 x/y/z vertices; faces as uint8 count + 3x uint32
indices). Merging into one file lets BOTH the CPU and the GPU load a single AABB grid
over all triangles -- the CPU otherwise builds one grid per PLY (O(n_PLY) per step).

Usage: merge_plys.py [MESH_DIR] [OUT.ply]
"""
import glob, os, sys, numpy as np

DIR = sys.argv[1] if len(sys.argv) > 1 else \
      "/home/jonathan/Documents/EPFL/Main_Studies/M_DISCO/decimated_holeless/decimated2"
OUT = sys.argv[2] if len(sys.argv) > 2 else \
      "/home/jonathan/Documents/EPFL/Main_Studies/M_DISCO/decimated_holeless/strands_merged.ply"

files = sorted(glob.glob(os.path.join(DIR, "*.ply")))
if not files:
    sys.exit(f"no PLYs in {DIR}")

# pass 1: total counts from headers
TV = TF = 0
hdr = []
for f in files:
    nv = nf = 0
    with open(f) as fh:
        for l in fh:
            if l.startswith("element vertex"): nv = int(l.split()[2])
            elif l.startswith("element face"): nf = int(l.split()[2])
            elif l.startswith("end_header"):   break
    hdr.append((nv, nf)); TV += nv; TF += nf
print(f"files={len(files)} total_verts={TV} total_faces={TF}", flush=True)

verts = np.empty((TV, 3), dtype=np.float32)
facev = np.empty((TF, 3), dtype=np.uint32)
vo = fo = vbase = 0
for fi, f in enumerate(files):
    nv, nf = hdr[fi]
    with open(f) as fh:
        for l in fh:
            if l.startswith("end_header"): break
        for _ in range(nv):
            p = fh.readline().split()
            verts[vo, 0] = p[0]; verts[vo, 1] = p[1]; verts[vo, 2] = p[2]; vo += 1
        for _ in range(nf):
            p = fh.readline().split()    # "3 i0 i1 i2"
            facev[fo, 0] = int(p[1]) + vbase
            facev[fo, 1] = int(p[2]) + vbase
            facev[fo, 2] = int(p[3]) + vbase; fo += 1
    vbase += nv
    if fi % 2000 == 0: print(f"  ...{fi}/{len(files)}", flush=True)

rec = np.empty(TF, dtype=np.dtype([('n', np.uint8), ('v', np.uint32, 3)]))  # packed, 13 B/face
assert rec.dtype.itemsize == 13, rec.dtype.itemsize
rec['n'] = 3; rec['v'] = facev
with open(OUT, "wb") as o:
    o.write(("ply\nformat binary_little_endian 1.0\n"
             f"element vertex {TV}\nproperty float x\nproperty float y\nproperty float z\n"
             f"element face {TF}\nproperty list uchar int vertex_indices\nend_header\n").encode("ascii"))
    verts.tofile(o); rec.tofile(o)
print(f"DONE {OUT}  {os.path.getsize(OUT)/1e6:.1f} MB", flush=True)
