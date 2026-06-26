# Tutorial: Multi-compartment substrate (multi-D, multi-T2, multi-PLY)

This example puts the per-compartment physics together: **three separate PLY meshes** (spheres),
each given its **own intra diffusivity and T2**, walked simultaneously, with a low-resolution
**voxel-subdivision volume** so we can read the signal back per sphere.

| sphere (PLY) | centre z | radius | D_intra | T2 |
|---|---|---|---|---|
| `sphere_a` | 0 µm | 5 µm | 0.2×10⁻⁹ m²/s | 140 ms |
| `sphere_b` | 21 µm | 5 µm | 0.5×10⁻⁹ m²/s | 160 ms |
| `sphere_c` | 42 µm | 5 µm | 1.0×10⁻⁹ m²/s | 180 ms |

The per-mesh properties come from an **extended PLY list**
(`docs/conf_file_examples/three_spheres_multi.list`):

```
# <ply>  <scale [m/file-unit]>  <d_intra [m^2/s]>  <T2 [s]>  <permeability [m/s]>
instructions/meshes/three_spheres/sphere_a.ply 5e-06 0.2e-9 0.140 0
instructions/meshes/three_spheres/sphere_b.ply 5e-06 0.5e-9 0.160 0
instructions/meshes/three_spheres/sphere_c.ply 5e-06 1.0e-9 0.180 0
```

(permeability `0` → impermeable, so each sphere keeps its own spins and its own D/T2.)

## The configuration

`docs/conf_file_examples/threeSpheresMultiCompartment.conf`:

```
N 100000
T 500
duration 0.050
diffusivity 2.0e-9
t2_extra 1e9

exp_prefix instructions/demos/output/three_spheres_multi
scheme_file docs/scheme_files/PGSE_multiTE.scheme

write_txt 1
write_bin 0
write_traj_file 0
separate_signals 1
deportation 0

<obstacle>
ply_extended_file_list docs/conf_file_examples/three_spheres_multi.list
</obstacle>

<voxels>
-6.0e-6 -6.0e-6 -6.0e-6
 6.0e-6  6.0e-6  4.8e-5
</voxels>

subdivisions_number 20
num_process 4
seed 12345
<END>
```

- `subdivisions_number 20` → a 20×20×20 volumetric grid (low quality, fast — finished in ~10 s).
- `separate_signals 1` → intra and extra signals are written separately.
- The scheme `PGSE_multiTE.scheme` has two parts: a **b-sweep** (b = 0…3000 s/mm², TE = 50 ms — for
  **diffusivity**) and a series of **b0 acquisitions at increasing TE** (10→50 ms — for **T2**).
- Spins are seeded uniformly; ~20 % land inside the three spheres (intra), the rest are extra.

## Run it

```bash
./MC-DC_Simulator docs/conf_file_examples/threeSpheresMultiCompartment.conf
```

## Output

Alongside the global/separated signals (`..._DWI[_intra|_extra].txt`), the subdivision volume is
written as `..._voxels_DWI[_intra|_extra].txt` (a 2-value header `n_sub num_rep`, then one signal per
sub-voxel per acquisition) and the per-sub-voxel densities as `..._volFractions.txt`.

## Visualize

**1. The three spheres (multi-PLY).** Convert the densities to NIfTI with the bundled
`src/volumeVolume.py` (uses `nibabel`), or check the intra-density profile along z — three peaks at
the sphere centres:

```python
import numpy as np
n = 20
vf = np.loadtxt("instructions/demos/output/three_spheres_multi_volFractions.txt").reshape(n, n, n, 3)
intra = vf[..., 1]                       # columns: total, intra, extra
print("intra density per z-slice:", intra.sum(axis=(0, 1)).astype(int))   # 3 peaks
```

**2. Multi-T2 and multi-D, read back per sphere.** The volume lets us pull each sphere's signal out
of its z-band and recover its T2 (from the b0-vs-TE decay) and its restriction (from the b-sweep):

```python
import numpy as np, matplotlib.pyplot as plt

pref = "instructions/demos/output/three_spheres_multi"; n = 20
sch = np.loadtxt("docs/scheme_files/PGSE_multiTE.scheme", skiprows=1)
G, TE = sch[:, 3], sch[:, 6]
b0 = G == 0                              # the multi-TE b0 acquisitions

raw = np.loadtxt(pref + "_voxels_DWI_intra.txt")
V = raw[2:].reshape(n, n, n, len(G))     # skip the [n_sub, num_rep] header

bands = {"a (D=0.2, T2=140)": range(0, 5),     # z-slices of each sphere centre
         "b (D=0.5, T2=160)": range(7, 13),
         "c (D=1.0, T2=180)": range(15, 20)}

for name, zr in bands.items():
    sig = V[:, :, list(zr), :].sum(axis=(0, 1, 2))
    te = TE[b0]; s = sig[b0]; o = np.argsort(te)
    T2 = -1.0 / np.polyfit(te[o], np.log(s[o]), 1)[0]      # S0·exp(-TE/T2)
    plt.plot(te[o] * 1e3, s[o] / s[o].max(), "o-", label=f"{name} → T2≈{T2*1e3:.0f} ms")

plt.xlabel("TE [ms]"); plt.ylabel("b0 signal (norm.)"); plt.legend()
plt.tight_layout(); plt.savefig("three_spheres_T2.png", dpi=120)
```

The fitted T2s come back at **≈ 140 / 160 / 180 ms** — the per-sphere values. The b-sweep part
(`G != 0`) likewise shows the lowest-D sphere (`a`) attenuating the least at high b.

(Requires `numpy`, `matplotlib`; the NIfTI step needs `nibabel`.)

---
See also: [free diffusion](tutorial_free_diffusion.md) · [gamma cylinders](tutorial_gamma_cylinders.md) · [PLY meshes](tutorial_ply_meshes.md)
