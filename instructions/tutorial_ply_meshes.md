# Tutorial: Simulation in PLY models

Restricted diffusion **inside** a closed, triangulated mesh of arbitrary shape (e.g. the
intra-axonal/intra-cellular space). Spins are seeded inside the mesh and reflect off its membrane.

> Units are SI by default (metres, seconds, Tesla). A mesh scale factor (`ply_scale`) is "metres per
> file unit". See the parameter reference in [Getting started](GettingStarted.md).

## The configuration

`docs/conf_file_examples/meshIntraInitialization.conf` (body):

```
N 1000
T 100
duration 0.050
diffusivity 0.6e-10

exp_prefix instructions/demos/output/mesh_initialization_test
scheme_file docs/scheme_files/PGSE_sample_scheme.scheme

write_txt 1
write_bin 0
write_traj_file 1

<obstacle>
ply instructions/meshes/PorusMedia.ply
ply_scale 1e-6
</obstacle>

<voxels>
-1.0e-5 -1.0e-5 -5.0e-5
 1.0e-5  1.0e-5  5.0e-5
</voxels>

<spawning_area>
-5.0e-6 -5.0e-6 -5.0e-6
 5.0e-6  5.0e-6  5.0e-6
</spawning_area>

ini_walkers_pos intra
num_process 1
<END>
```

`PorusMedia.ply` is stored in **micrometres**, so `ply_scale 1e-6` maps its vertices to metres. The
`<voxels>` and `<spawning_area>` are in metres; `ini_walkers_pos intra` keeps only spins that land
**inside** the mesh (seeded uniformly in the spawning area, then filtered).

## Run it

```bash
./MC-DC_Simulator docs/conf_file_examples/meshIntraInitialization.conf
```

(The "highly irregular triangles" warning for this demo mesh is benign.)

## Expected output (under `instructions/demos/output/`)

| file | contents |
|---|---|
| `mesh_initialization_test_DWI.txt` | real part of the DW signal (123 acquisitions) |
| `mesh_initialization_test_0.traj.txt` | per-spin trajectories (`N·(T+1)·3` values) |
| `mesh_initialization_test_simulation_info.txt` | run parameters and warnings |

## Visualize

**The container mesh:** open `instructions/meshes/PorusMedia.ply` in a mesh viewer
(MeshLab, ParaView, or `trimesh`); the repository also ships `src/visIntersected.py` (trimesh-based)
for inspecting a mesh together with a per-face mask.

**The intra trajectories** — confined inside the mesh:

```python
import numpy as np, matplotlib.pyplot as plt

T = 100                                  # = T in the .conf
traj = np.loadtxt("instructions/demos/output/mesh_initialization_test_0.traj.txt").reshape(-1, T+1, 3)
ax = plt.figure().add_subplot(111, projection="3d")
for p in traj[:80]:                      # first 80 spins (positions in mm)
    ax.plot(p[:, 0], p[:, 1], p[:, 2], lw=0.4)
ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")
plt.tight_layout(); plt.savefig("ply_intra_trajectories.png", dpi=120)
```

(Requires `numpy` and `matplotlib`.) The trajectories trace out the interior of the mesh, confirming
the intra-compartment initialization. The signal can be read/plotted vs. b as in the
[free-diffusion tutorial](tutorial_free_diffusion.md#visualize).

---
Previous: [gamma-distributed cylinders](tutorial_gamma_cylinders.md) · Back to [free diffusion](tutorial_free_diffusion.md)
