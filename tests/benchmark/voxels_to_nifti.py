#!/usr/bin/env python3
"""
Convert an MC-DC volumetric DWI (from subdivisions_number) into a 4-D NIfTI you can
open in fsleyes / mrview / ITK-SNAP. Accepts either output format:
  *_voxels_DWI.txt    -- text: line1 = n_sub, line2 = num_rep, then n_sub*num_rep values.
  *_voxels_DWI.bfloat -- binary float32: [n_sub, num_rep, then n_sub*num_rep values]
                         (the fast `write_bin 1` output; ~20x quicker to write/read).
Values are in subdivision order (x outer, y middle, z inner -- matching
Parameters::addSubdivisions), so they reshape to (n, n, n, num_rep), n = round(n_sub**(1/3)).

Output:
  <prefix>_volume.nii.gz   -- the full (n,n,n,num_rep) signal (num_rep "DWI volumes")
  <prefix>_b0.nii.gz       -- just the first acquisition (n,n,n), handy as a density-ish map

Usage:
  python3 voxels_to_nifti.py <voxels_DWI.txt|.bfloat> [voxel_size_mm] [out_prefix]
    voxel_size_mm : physical size of ONE subdivision along each axis (default 1.0; only
                    affects the NIfTI header spacing, not the values).
"""
import sys, numpy as np
try:
    import nibabel as nib
except ImportError:
    sys.exit("nibabel not installed:  pip install nibabel")

def load_voxels(path):
    """Return (n_sub, num_rep, body float32 array). Handles big text files without
    building a giant Python list (np.fromstring parses the C way), and .bfloat directly."""
    if path.endswith(".bfloat"):
        raw = np.fromfile(path, dtype=np.float32)
        n_sub, num_rep = int(raw[0]), int(raw[1])
        return n_sub, num_rep, raw[2:]
    with open(path) as f:
        n_sub = int(f.readline()); num_rep = int(f.readline())
        body = np.fromstring(f.read(), sep="\n", dtype=np.float32)   # C-fast text parse
    return n_sub, num_rep, body

def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    vox  = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    out  = sys.argv[3] if len(sys.argv) > 3 else path.rsplit(".", 1)[0]

    n_sub, num_rep, body = load_voxels(path)
    assert body.size == n_sub * num_rep, f"expected {n_sub*num_rep} values, got {body.size}"
    n = round(n_sub ** (1/3))
    assert n*n*n == n_sub, f"n_sub={n_sub} is not a perfect cube"

    vol = body.reshape(n_sub, num_rep).reshape(n, n, n, num_rep)  # (x,y,z,acq)
    aff = np.diag([vox, vox, vox, 1.0])

    nib.save(nib.Nifti1Image(vol, aff), out + "_volume.nii.gz")
    nib.save(nib.Nifti1Image(np.ascontiguousarray(vol[..., 0]), aff), out + "_b0.nii.gz")
    nz = int(np.count_nonzero(vol[..., 0]))
    print(f"n={n} (={n_sub} subdivisions), num_rep={num_rep}")
    print(f"populated subdivisions (nonzero b0): {nz}/{n_sub}")
    print(f"wrote {out}_volume.nii.gz  (shape {vol.shape})")
    print(f"wrote {out}_b0.nii.gz      (shape {vol[...,0].shape})")

if __name__ == "__main__":
    main()
