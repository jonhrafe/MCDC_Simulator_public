#!/usr/bin/env python3
"""
Freeze analytic reference curves for the validation suite (DEV-ONLY).

Runs dmipy ONCE to compute a closed-form diffusion signal for a fixed
(geometry, D, scheme) and writes it as a committed ``.npy`` under
``tests/accuracy/analytic/`` plus a ``.provenance.json``. The test-time suite
then loads the frozen array (numpy only) -- dmipy is NOT needed to RUN the tests.

dmipy is heavy and version-fragile, so it lives in the ``analytic`` extra:
    pip install .[analytic]        # dmipy + numpy
    python3 tests/tools/gen_analytic.py

Two oracles, selected with ``--model`` (both are Gaussian-Phase-Approximation,
small-b, impermeable restricted-diffusion models):

  --model cylinder  (default, for cylinder_restricted.conf):
      C4CylinderGaussianPhaseApproximation (van Gelderen), radius 5 um,
      D 0.6e-9 m^2/s, on PGSE_sample_scheme.scheme  ->  cylinder_r5um_d0p6.npy
  --model sphere    (for sphere_restricted.conf):
      S4SphereGaussianPhaseApproximation (Neuman / Murday-Cotts, the sphere
      analogue), radius 5 um, D 0.6e-9 m^2/s, same scheme  ->  sphere_r5um_d0p6.npy

The sphere radius (5 um) matches tests/accuracy/single_sphere.list and D matches
sphere_restricted.conf, so the frozen curve lines up with what test_sphere_neuman
runs. The cylinder mirror is tests/accuracy/single_cylinder.list / cylinder_restricted.conf.

VALIDITY: the Gaussian Phase Approximation is an INFINITE (cylinder) / closed
(sphere), IMPERMEABLE geometry valid only at small b (<~5000 s/mm^2). The tool
refuses to freeze a scheme whose max b exceeds that (use --force to override
knowingly).

The frozen array is FULL-LENGTH S/S0 (one value per scheme row, in scheme order,
b0 entries ~1.0), normalised as ``signal / mean(signal[b0])`` -- identical to how
``mcdc_testkit.metrics.normalize_by_b0`` normalises the simulator's DWI, so the
test compares like with like.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # tests/
from mcdc_testkit import paths, provenance  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402

import numpy as np  # noqa: E402


def _patch_dmipy():
    """Fix dmipy's ``if deltas == []`` (raises on modern numpy: array vs []).

    Runtime monkeypatch for the current dmipy/numpy combination.
    """
    import inspect
    import textwrap

    import dmipy.core.acquisition_scheme as _acq
    src = textwrap.dedent(inspect.getsource(_acq.DmipyAcquisitionScheme.__init__))
    src = src.replace("if deltas == []:", "if len(deltas) == 0:")
    ns = {}
    exec(src, _acq.__dict__, ns)  # noqa: S102 - deliberate, mirrors the reference script
    _acq.DmipyAcquisitionScheme.__init__ = ns["__init__"]


def cylinder_curve(scheme_path, radius, diffusivity):
    """Return (full-length S/S0 array, dmipy_version, model_name) for a cylinder."""
    _patch_dmipy()
    import dmipy
    from dmipy.core.acquisition_scheme import acquisition_scheme_from_schemefile
    from dmipy.signal_models.cylinder_models import C4CylinderGaussianPhaseApproximation

    acq = acquisition_scheme_from_schemefile(scheme_path)
    model = C4CylinderGaussianPhaseApproximation(
        mu=[0.0, 0.0], lambda_par=diffusivity,
        diffusion_perpendicular=diffusivity, diameter=2.0 * radius)
    signal = np.asarray(model(acq), dtype=np.float64)

    sc = kscheme.parse_scheme(scheme_path)
    if signal.size != sc.n:
        raise AssertionError(f"dmipy produced {signal.size} values != {sc.n} scheme rows")
    ss0 = signal / np.mean(signal[sc.b0_mask])   # full length, b0 ~ 1.0
    return ss0, dmipy.__version__, "C4CylinderGaussianPhaseApproximation"


def sphere_curve(scheme_path, radius, diffusivity):
    """Return (full-length S/S0 array, dmipy_version, model_name) for a sphere.

    The sphere analogue of :func:`cylinder_curve`: dmipy's Neuman / Murday-Cotts
    Gaussian Phase Approximation for restricted diffusion inside an impermeable
    sphere. A sphere is ISOTROPIC, so the signal depends only on (G, Delta,
    delta) -- every gradient direction in a shell shares one value (unlike the
    cylinder, which is angle-dependent). Same monkeypatch and scheme loading /
    normalisation as the cylinder path so the frozen array is directly
    comparable to the simulator's normalised DWI.
    """
    _patch_dmipy()
    import dmipy
    from dmipy.core.acquisition_scheme import acquisition_scheme_from_schemefile
    from dmipy.signal_models.sphere_models import S4SphereGaussianPhaseApproximation

    acq = acquisition_scheme_from_schemefile(scheme_path)
    # Isotropic restricted sphere: one diffusion constant, no orientation (no mu).
    model = S4SphereGaussianPhaseApproximation(
        diameter=2.0 * radius, diffusion_constant=diffusivity)
    signal = np.asarray(model(acq), dtype=np.float64)

    sc = kscheme.parse_scheme(scheme_path)
    if signal.size != sc.n:
        raise AssertionError(f"dmipy produced {signal.size} values != {sc.n} scheme rows")
    ss0 = signal / np.mean(signal[sc.b0_mask])   # full length, b0 ~ 1.0
    return ss0, dmipy.__version__, "S4SphereGaussianPhaseApproximation"


# model key -> (curve function, default output basename, GPA note for the manifest)
MODELS = {
    "cylinder": (cylinder_curve, "cylinder_r5um_d0p6",
                 "dmipy Gaussian Phase Approximation (van Gelderen): infinite "
                 "impermeable cylinder, small-b only."),
    "sphere":   (sphere_curve, "sphere_r5um_d0p6",
                 "dmipy Gaussian Phase Approximation (Neuman / Murday-Cotts): "
                 "impermeable sphere, small-b only. Isotropic."),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", choices=sorted(MODELS), default="cylinder",
                    help="restricted geometry to freeze (default: cylinder)")
    ap.add_argument("--radius", type=float, default=5e-6, help="cylinder/sphere radius [m]")
    ap.add_argument("--diffusivity", type=float, default=0.6e-9, help="D [m^2/s]")
    ap.add_argument("--scheme", default=paths.accuracy("PGSE_sample_scheme.scheme"))
    ap.add_argument("--name", default=None,
                    help="output basename (no extension); default depends on --model")
    ap.add_argument("--out-dir", default=paths.ANALYTIC_DIR)
    ap.add_argument("--force", action="store_true",
                    help="freeze even if max b exceeds the GPA small-b validity limit")
    a = ap.parse_args()

    curve_fn, default_name, gpa_note = MODELS[a.model]
    name = a.name or default_name

    sc = kscheme.parse_scheme(a.scheme)
    max_b_mm = float(np.max(sc.bvals) / 1e6)
    if max_b_mm > 5000.0 and not a.force:
        print(f"REFUSE: max b = {max_b_mm:.0f} s/mm^2 exceeds the GPA small-b validity "
              f"(~5000 s/mm^2). The dmipy {a.model} curve would be unreliable. "
              f"Re-run with --force only if you know what you are doing.", file=sys.stderr)
        return 1

    curve, dmipy_ver, model_name = curve_fn(a.scheme, a.radius, a.diffusivity)

    os.makedirs(a.out_dir, exist_ok=True)
    npy = os.path.join(a.out_dir, name + ".npy")
    np.save(npy, curve)
    manifest = provenance.make_manifest(
        kind="analytic",
        model=model_name,
        dmipy_version=dmipy_ver,
        radius_m=a.radius,
        diffusivity_m2s=a.diffusivity,
        scheme=os.path.relpath(a.scheme, paths.REPO_ROOT),
        scheme_sha256=provenance.sha256_file(a.scheme),
        n_measurements=int(sc.n),
        max_b_s_per_mm2=round(max_b_mm, 1),
        normalization="full-length S/S0 = signal / mean(signal[b0]); b0 entries ~1.0; scheme row order",
        note=gpa_note,
    )
    with open(os.path.join(a.out_dir, name + ".provenance.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {os.path.relpath(npy, paths.REPO_ROOT)} ({curve.size} values) + provenance "
          f"[{a.model}: dmipy {dmipy_ver}, R={a.radius:g} m, D={a.diffusivity:g} m^2/s, "
          f"max b={max_b_mm:.0f} s/mm^2]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
