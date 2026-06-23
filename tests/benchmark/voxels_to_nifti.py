#!/usr/bin/env python3
"""
Convert an MC-DC `*_voxels_DWI.txt` (volumetric DWI from subdivisions_number) into a
4-D NIfTI you can open in fsleyes / mrview / ITK-SNAP.

The text file is: line1 = n_sub, line2 = num_rep, then n_sub*num_rep signal values in
subdivision order (x outer, y middle, z inner -- matching Parameters::addSubdivisions),
so it reshapes to (n, n, n, num_rep) where n = round(n_sub**(1/3)).

Output:
  <prefix>_volume.nii.gz   -- the full (n,n,n,num_rep) signal (num_rep "DWI volumes")
  <prefix>_b0.nii.gz       -- just the first acquisition (n,n,n), handy as a density-ish map

Usage:
  python3 voxels_to_nifti.py <voxels_DWI.txt> [voxel_size_mm] [out_prefix]
    voxel_size_mm : physical size of ONE subdivision along each axis (default 1.0; only
                    affects the NIfTI header spacing, not the values).
"""
import sys, numpy as np
try:
    import nibabel as nib
except ImportError:
    sys.exit("nibabel not installed:  pip install nibabel")

def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    vox  = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    out  = sys.argv[3] if len(sys.argv) > 3 else path.rsplit(".", 1)[0]

    vals = [float(x) for x in open(path) if x.strip()]
    n_sub, num_rep = int(vals[0]), int(vals[1])
    body = np.asarray(vals[2:], dtype=np.float32)
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
