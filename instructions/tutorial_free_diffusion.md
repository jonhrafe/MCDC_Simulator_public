# Tutorial: Simulation in free diffusion

The simplest simulation: spins diffusing in an unrestricted medium (no obstacles). The
diffusion-weighted signal must follow the free-diffusion law **S(b)/S₀ = exp(−b·D)**, which makes
this a good first run and a sanity check of your build.

> Units are SI by default (metres, seconds, Tesla). See the parameter reference in
> [Getting started](GettingStarted.md).

## The configuration

`docs/conf_file_examples/freeDiffusion.conf` (body):

```
N 1000
T 1000
duration 0.100
diffusivity 0.6e-9

exp_prefix instructions/demos/output/free_diffusion_test
scheme_file docs/scheme_files/PGSE_sample_scheme.scheme

write_txt 1
write_bin 1
write_traj_file 1

<voxel>
0 0 0
0.5e-3 0.5e-3 0.5e-3
</voxel>

num_process 1
<END>
```

1000 spins diffuse for 100 ms (1000 steps) with D = 0.6×10⁻⁹ m²/s, in a 0.5 mm cubic voxel. The
PGSE scheme is the bundled multi-shell protocol (3 b0 + b = 500…2000 s/mm², TE = 45 ms).

## Run it

```bash
./MC-DC_Simulator docs/conf_file_examples/freeDiffusion.conf
```

## Expected output (under `instructions/demos/output/`)

| file | contents |
|---|---|
| `free_diffusion_test_DWI.txt` | real part of the DW signal — one value per acquisition (123) |
| `free_diffusion_test_DWI_img.txt` | imaginary part (≈ 0 for isotropic free diffusion) |
| `free_diffusion_test_simulation_info.txt` | run parameters and any warnings |
| `free_diffusion_test_0.traj.txt` | per-spin trajectories (`N·(T+1)·3` values) |

## Visualize

**1. The signal decay** — compute the b-value of each acquisition and check S/S₀ ≈ exp(−bD):

```python
import numpy as np, matplotlib.pyplot as plt

GAMMA = 2.6751525e8                      # gyromagnetic ratio [rad/(s·T)]
sch = np.loadtxt("docs/scheme_files/PGSE_sample_scheme.scheme", skiprows=1)
G, Delta, delta = sch[:, 3], sch[:, 4], sch[:, 5]
b = (GAMMA * G * delta) ** 2 * (Delta - delta / 3.0)        # s/m^2

S  = np.loadtxt("instructions/demos/output/free_diffusion_test_DWI.txt")
S0 = S[b < 1e6].mean()                   # average of the b0 (G=0) acquisitions
D  = 0.6e-9

order = np.argsort(b)
plt.semilogy(b[order]*1e-6, (S/S0)[order], "o", label="simulated")
plt.semilogy(b[order]*1e-6, np.exp(-b[order]*D), "-", label="exp(-bD)")
plt.xlabel("b  [s/mm$^2$]"); plt.ylabel("S / S$_0$"); plt.legend(); plt.tight_layout()
plt.savefig("free_diffusion_signal.png", dpi=120)
```

The points fall on the `exp(-bD)` line (≈ 0.30 at b = 2000 s/mm² for this D).

**2. The trajectories** — the diffusing cloud of spins:

```python
import numpy as np, matplotlib.pyplot as plt

T = 1000                                 # = T in the .conf
traj = np.loadtxt("instructions/demos/output/free_diffusion_test_0.traj.txt").reshape(-1, T+1, 3)
ax = plt.figure().add_subplot(111, projection="3d")
for p in traj[:40]:                      # first 40 spins
    ax.plot(p[:, 0], p[:, 1], p[:, 2], lw=0.4)
plt.tight_layout(); plt.savefig("free_diffusion_trajectories.png", dpi=120)
```

(Requires `numpy` and `matplotlib`.)

---
Next: [gamma-distributed cylinders](tutorial_gamma_cylinders.md) · [PLY meshes](tutorial_ply_meshes.md)
