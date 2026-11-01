# Dev / CI tools (produce data artifacts)

Utilities that *generate* committed reference data or synthetic inputs. They are
not tests themselves and are not run by CTest; you run them by hand when a
reference needs (re)freezing, then commit the result.

| tool | purpose | needs |
|------|---------|-------|
| `gen_analytic.py` | Freeze dmipy analytic curves (e.g. the van Gelderen restricted-cylinder) to `tests/accuracy/analytic/*.npy` + provenance | `pip install .[analytic]` (dmipy) |
| `build_old_version.py` *(Phase C)* | Build a historical release binary (v1.5, v2.0) from THIS repo's own git history (`git worktree` + compile) and cache it under `tests/crossversion/old_binaries/` | `git`, `g++` |
| `gen_reference.py` *(Phase C)* | Run a version's binary at a fixed seed → freeze `tests/accuracy/reference_xv/` + manifest. Auto-builds v1.5/v2.0 via `build_old_version.py` on first use; v2.2 still needs an external binary | nothing for v1.5/v2.0; external binary for v2.2 |
| `compute_ecvf.py` *(Phase D)* | ECVF from the b0 intra/extra split (weak mesh oracle) | numpy |
| `make_meshes.py` *(Phase D)* | Generate convex / non-convex synthetic PLYs for bbox + classification tests | numpy |

Freeze the default cylinder oracle (run in an env that has dmipy):

    pip install .[analytic]
    python3 tests/tools/gen_analytic.py

This writes `tests/accuracy/analytic/cylinder_r5um_d0p6.npy`; commit it so CI /
the validation tier (numpy-only) can load it without dmipy.
