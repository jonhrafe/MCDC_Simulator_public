# Frozen gamma-packing substrate lists

Realised tissue-scale packings written by the simulator's gamma generator
(`<exp_prefix>_gamma_distributed_{cylinder,sphere}_list.txt`), kept here as
reusable **native list substrates** and human-inspectable examples of what the
generator produces.

| File | From conf | Realised |
|------|-----------|----------|
| `gamma_cylinder_100um_list.txt` | `../gamma_cylinder_100um.conf` | 340 cylinders, mean r 1.50 µm, ICVF 0.40, voxel ~103 µm |
| `gamma_sphere_100um_list.txt`   | `../gamma_sphere_100um.conf`   | 900 spheres,  mean r 3.00 µm, ICVF 0.40, voxel ~102 µm |

Format (native list): line 1 is the scale factor (`0.001`, i.e. file units are µm
→ ×scale → internal mm); each remaining row is `x y z radius percolation T2` in
file units, including periodic-image replicas. Feed one back as a substrate with

```
<obstacle>
cylinders_list tests/accuracy/gamma_lists/gamma_cylinder_100um_list.txt
</obstacle>
<voxels>
0 0 0
1.031e-4 1.031e-4 1.031e-4      # the auto-sized voxel, in metres (SI)
</voxels>
```

**These are same-build examples, not a bit-exact golden.** The generator is
deterministic for a fixed `seed` + `num_process`, but the exact float stream can
differ across platforms/compilers, so regenerate from the `.conf` on your build
rather than diffing against these byte-for-byte. `tests/validation/test_gamma_packing.py`
regenerates live and checks the packing's *invariants* (ICVF, radius
distribution, non-overlap, auto-voxel, determinism) instead of exact bytes.

Regenerate:

```bash
./MC-DC_Simulator tests/accuracy/gamma_cylinder_100um.conf   # writes *_gamma_distributed_cylinder_list.txt
./MC-DC_Simulator tests/accuracy/gamma_sphere_100um.conf     # writes *_gamma_distributed_sphere_list.txt
```
