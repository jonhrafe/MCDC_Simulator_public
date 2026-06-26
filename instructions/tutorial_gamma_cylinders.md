# Tutorial: Simulation in gamma-distributed cylinders

A restricted-diffusion substrate: a packing of parallel cylinders whose radii are drawn from a
Gamma distribution (a common white-matter / axon model). The packer fills the voxel to a target
intra-cylinder volume fraction (ICVF) and sizes the voxel automatically.

> Units are SI by default (metres, seconds, Tesla); the Gamma `beta`/`min_radius` are lengths in
> metres. See the parameter reference in [Getting started](GettingStarted.md).

## The configuration

`docs/conf_file_examples/gammaDistributedCylinders.conf` (body):

```
N 1000
T 1000
duration 0.050
diffusivity 0.6e-9

exp_prefix instructions/demos/output/cylinder_gamma_packing_test
scheme_file docs/scheme_files/PGSE_sample_scheme.scheme

write_txt 1
write_bin 0
write_traj_file 1

<obstacle>
<cylinder_gamma_packing>
alpha 1.5
beta 0.5e-6
icvf 0.70
num_cylinders 1000
</cylinder_gamma_packing>
</obstacle>

ini_walkers_pos extra
num_process 1
<END>
```

`alpha`/`beta` are the Gamma shape/scale (mean radius ≈ `alpha·beta` = 0.75 µm), `icvf 0.70` is the
target packing fraction, and `ini_walkers_pos extra` seeds spins **outside** the cylinders.

## Run it

```bash
./MC-DC_Simulator docs/conf_file_examples/gammaDistributedCylinders.conf
```

The log reports the achieved ICVF (close to the 0.70 target).

## Expected output (under `instructions/demos/output/`)

| file | contents |
|---|---|
| `cylinder_gamma_packing_test_DWI.txt` | real part of the DW signal (123 acquisitions) |
| `cylinder_gamma_packing_test_gamma_distributed_cylinder_list.txt` | the packed cylinders: first line is the file→mm scale, then `x y z r perc T2` per cylinder |
| `cylinder_gamma_packing_test_0.traj.txt` | per-spin trajectories (`N·(T+1)·3` values) |
| `cylinder_gamma_packing_test_simulation_info.txt` | run parameters and warnings |

## Visualize

The cylinders are parallel to **z**, so an x–y view shows the packing as circles with the spin
trajectories weaving through the extra-cylinder space:

```python
import numpy as np, matplotlib.pyplot as plt

pref = "instructions/demos/output/cylinder_gamma_packing_test"
T = 1000                                            # = T in the .conf

# packed cylinders (first line = file->mm scale; columns x y z r perc T2)
scale = float(open(pref + "_gamma_distributed_cylinder_list.txt").readline())
cyl = np.loadtxt(pref + "_gamma_distributed_cylinder_list.txt", skiprows=1)
cx, cy, r = cyl[:, 0]*scale, cyl[:, 1]*scale, cyl[:, 3]*scale   # -> mm (matches trajectories)

fig, ax = plt.subplots(figsize=(6, 6))
for x, y, rad in zip(cx, cy, r):
    ax.add_patch(plt.Circle((x, y), rad, color="0.8", ec="0.5"))

traj = np.loadtxt(pref + "_0.traj.txt").reshape(-1, T+1, 3)
for p in traj[:60]:                                 # first 60 spins, x-y projection
    ax.plot(p[:, 0], p[:, 1], lw=0.4)

ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
ax.autoscale(); plt.tight_layout(); plt.savefig("gamma_cylinders.png", dpi=120)
```

(Requires `numpy` and `matplotlib`.) The signal can be read and plotted vs. b exactly as in the
[free-diffusion tutorial](tutorial_free_diffusion.md#visualize); here it will sit **above** the
free-diffusion curve at high b, because the cylinders restrict the spins.

---
Previous: [free diffusion](tutorial_free_diffusion.md) · Next: [PLY meshes](tutorial_ply_meshes.md)
