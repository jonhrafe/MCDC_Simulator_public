#!/usr/bin/env python3
"""
Permeable / two-compartment diffusion in a GAMMA PACKING (experiments 3 & 4).

Same calibrated substrates as ``test_packing_restricted.py``, now with membrane
permeability and/or an extra compartment. There is NO closed form for the DWI of a
polydisperse *permeable* packing, so the oracle is a combination of limits,
conservation and monotonicity anchored on the validated restricted GPA:

  EXP 3 -- permeable, intra-seeded, kappa sweep (sphere + cylinder):
    * kappa->0 recovers the volume-weighted restricted GPA (the strong anchor:
      at kappa=0 nothing leaks, so intra == the impermeable restricted signal ==
      oracles.packing_signal, checked to <2%);
    * the extra fraction at b0 increases monotonically with kappa (exchange out of
      the intra compartment grows with permeability);
    * the total apparent diffusivity increases monotonically with kappa (exchange
      relaxes the restriction toward the faster extra compartment).

  EXP 4 -- two-compartment, seeded uniformly over the voxel (intra+extra ~ ICVF):
    * sphere IMPERMEABLE: the intra signal (via separate_signals) still equals the
      volume-weighted restricted GPA (strong, same oracle as exp 1); the extra
      signal is HINDERED (0 < ADC_extra < d_extra); and the separate-signal
      conservation identity S_tot = f_i S_intra + f_e S_extra holds.
    * cylinder PERMEABLE: conservation + hindered-extra bound (full exchange; the
      restricted anchor no longer applies to the total, only bounds/limits do).

A Karger two-compartment curve (oracles.karger) is drawn on the plot as a REFERENCE
only -- it is a Gaussian-compartment model, so it does not describe the restricted
(non-Gaussian) intra compartment and is never a pass/fail oracle here.

  python3 tests/validation/test_packing_exchange.py --bin ./MC-DC_Simulator --n 20000
"""
import argparse
import os
import re
import shutil
import sys
import tempfile
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mcdc_testkit import io, metrics, oracles, paths, run, substrate  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402

import numpy as np  # noqa: E402

CASES = {
    "sphere":   dict(conf="gamma_sphere_100um.conf",   geometry="sphere",   dims=(0, 1, 2)),
    "cylinder": dict(conf="gamma_cylinder_100um.conf", geometry="cylinder", dims=(0, 1)),
}
D_EXTRA = 2.0e-9            # extra-compartment diffusivity (m^2/s)
KAPPA_SWEEP = [0.0, 1e-5, 1e-4]     # membrane permeability sweep (m/s); 1e-5 is the anchor
KAPPA_ANCHOR = 1e-5

Variant = namedtuple("Variant", "sc D radii geometry intra extra total fi fe")


# --------------------------------------------------------------------------- #
#  Run one (substrate, kappa, seeding) variant                                #
# --------------------------------------------------------------------------- #
def run_variant(binary, name, kappa, seeding, n=20000, gpu=0, timeout=1200):
    """Run a packing with given permeability + seeding; return a Variant of S/S0.

    ``seeding`` is 'intra' (keep the conf's ``ini_walkers_pos intra``) or
    'uniform' (strip it -> volume-uniform over the voxel). intra/extra/total are
    each normalised by their OWN b0 mean (so each is a compartment S/S0 starting
    at 1); fi/fe are the intra/extra volume fractions from the total-b0 split.
    """
    case = CASES[name]
    text = open(paths.accuracy(case["conf"])).read()
    D = float(io.read_conf_value(text, "diffusivity"))
    sc = kscheme.parse_scheme(os.path.join(paths.REPO_ROOT,
                                           io.read_conf_value(text, "scheme_file")))
    if seeding == "uniform":
        text = re.sub(r"(?m)^ini_walkers_pos.*\n", "", text)
    ov = {"N": n, "separate_signals": 1, "d_extra": repr(D_EXTRA)}
    if kappa > 0:
        ov["obstacle_permeability"] = repr(kappa)

    work = tempfile.mkdtemp(prefix=f"mcdc_ex_{name}_")
    try:
        prefix = run.run_case(binary, text, work, name, paths.REPO_ROOT,
                              overrides=ov, gpu=gpu, timeout=timeout)
        if prefix is None:
            return None
        tot = io.load_dwi(prefix + "_DWI.txt", n_walkers=n)
        intra = io.load_dwi(prefix + "_DWI_intra.txt", n_walkers=n)
        extra = io.load_dwi(prefix + "_DWI_extra.txt", n_walkers=n)
        voxel_um = substrate.voxel_side_um(prefix + "_simulation_info.txt")
        centres, radii = substrate.read_list(f"{prefix}_gamma_distributed_{name}_list.txt")
        _, pr = substrate.primary(centres, radii, voxel_um, case["dims"])
        b0 = sc.b0_mask
        fi = float(intra[b0].mean() / tot[b0].mean())
        fe = float(extra[b0].mean() / tot[b0].mean())
        # per-compartment S/S0 (own b0); extra may be ~0 when nothing leaked
        intra_n = intra / max(intra[b0].mean(), 1e-30)
        extra_n = extra / max(extra[b0].mean(), 1e-30)
        total_n = metrics.normalize_by_b0(tot, b0)
        return Variant(sc, D, pr * 1e-6, case["geometry"], intra_n, extra_n, total_n, fi, fe)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def _shell_adc(ss0, sc, b_target=1000.0):
    """Apparent diffusivity from the per-shell mean at the b closest to b_target."""
    nb = ~sc.b0_mask
    means = metrics.per_shell_mean(ss0[nb], sc.bvals[nb])
    b = min(means, key=lambda x: abs(x - b_target))
    return float(-np.log(max(means[b], 1e-12)) / (b * 1e6)), b


# --------------------------------------------------------------------------- #
#  Experiment 3: permeable, kappa sweep                                       #
# --------------------------------------------------------------------------- #
def evaluate_exp3(binary, name, n=20000, gpu=0, tol=0.02):
    """kappa->0 GPA anchor + monotone extra-fraction + monotone ADC. (ok, msgs)."""
    vs = [(k, run_variant(binary, name, k, "intra", n=n, gpu=gpu)) for k in KAPPA_SWEEP]
    if any(v is None for _, v in vs):
        return None, [f"  {name}: GPU-skipped"]
    msgs, ok = [f"[exp3 {name} permeable | intra-seeded | N={n}]"], True

    # anchor: kappa=0 intra == volume-weighted restricted GPA
    v0 = dict(vs)[0.0]
    gpa = oracles.packing_signal(v0.sc, v0.radii, v0.D, v0.geometry)
    nb = ~v0.sc.b0_mask
    a_mc = metrics.per_shell_mean(v0.intra[nb], v0.sc.bvals[nb])
    a_an = metrics.per_shell_mean(gpa[nb], v0.sc.bvals[nb])
    worst = max(abs(a_mc[b] - a_an[b]) for b in a_mc)
    ok = ok and worst <= tol
    msgs.append(f"  [{'OK ' if worst <= tol else 'BAD'}] kappa=0 intra vs GPA: worst |Δ|={worst:.4f} (tol {tol})")

    # monotone extra fraction and total ADC with kappa
    fes = [v.fe for _, v in vs]
    adcs = [_shell_adc(v.total, v.sc)[0] for _, v in vs]
    mono_fe = all(x < y + 1e-6 for x, y in zip(fes, fes[1:]))
    mono_adc = all(x < y + 1e-13 for x, y in zip(adcs, adcs[1:]))
    ok = ok and mono_fe and mono_adc
    msgs.append(f"  [{'OK ' if mono_fe else 'BAD'}] extra fraction ↑ with kappa: "
                + ", ".join(f"{k:g}:{f:.3f}" for (k, _), f in zip(vs, fes)))
    msgs.append(f"  [{'OK ' if mono_adc else 'BAD'}] total ADC ↑ with kappa: "
                + ", ".join(f"{k:g}:{a*1e9:.3f}" for (k, _), a in zip(vs, adcs)) + " (1e-9 m²/s)")
    return ok, msgs


# --------------------------------------------------------------------------- #
#  Experiment 4: two-compartment (uniform seeding)                            #
# --------------------------------------------------------------------------- #
def evaluate_exp4(binary, name, kappa, n=20000, gpu=0, tol=0.02):
    """Conservation + hindered-extra bound (+ intra==GPA when impermeable)."""
    v = run_variant(binary, name, kappa, "uniform", n=n, gpu=gpu)
    if v is None:
        return None, [f"  {name}: GPU-skipped"]
    label = "impermeable" if kappa == 0 else f"kappa={kappa:g}"
    msgs, ok = [f"[exp4 {name} two-compartment | uniform seed | {label} | N={n}]"], True
    nb = ~v.sc.b0_mask

    # conservation: total S/S0 == fi*intra + fe*extra (own-b0 normalised)
    recon = v.fi * v.intra + v.fe * v.extra
    cons = metrics.max_abs_dev(v.total[nb], recon[nb])
    ok = ok and cons <= 5e-3
    msgs.append(f"  [{'OK ' if cons <= 5e-3 else 'BAD'}] conservation S_tot=fi·Si+fe·Se: "
                f"max|Δ|={cons:.4f}  (fi={v.fi:.3f} fe={v.fe:.3f})")

    # extra is hindered: 0 < ADC_extra < d_extra
    adc_e, be = _shell_adc(v.extra, v.sc)
    hindered = 0.0 < adc_e < D_EXTRA
    ok = ok and hindered
    msgs.append(f"  [{'OK ' if hindered else 'BAD'}] extra hindered: ADC_extra={adc_e*1e9:.3f} "
                f"< d_extra={D_EXTRA*1e9:.1f} (1e-9 m²/s, at b={be:.0f})")

    # impermeable sphere: intra still matches the restricted GPA (strong)
    if kappa == 0:
        gpa = oracles.packing_signal(v.sc, v.radii, v.D, v.geometry)
        a_mc = metrics.per_shell_mean(v.intra[nb], v.sc.bvals[nb])
        a_an = metrics.per_shell_mean(gpa[nb], v.sc.bvals[nb])
        worst = max(abs(a_mc[b] - a_an[b]) for b in a_mc)
        ok = ok and worst <= tol
        msgs.append(f"  [{'OK ' if worst <= tol else 'BAD'}] intra vs GPA: worst |Δ|={worst:.4f} (tol {tol})")
    return ok, msgs


# --- pytest entry points (analytic tier; needs pytest) -----------------------
def _binary_or_skip():
    import pytest
    binary = paths.binary_path()
    if not os.path.exists(binary):
        pytest.skip(f"simulator binary not built: {binary}")
    return binary


def _assert(ok, msgs):
    import pytest
    if ok is None:
        pytest.skip("binary built without CUDA / no GPU available")
    assert ok, "\n".join(msgs)


def test_exp3_sphere_permeable():
    """Exp 3: permeable sphere packing -- kappa->0 GPA + monotone exchange."""
    _assert(*evaluate_exp3(_binary_or_skip(), "sphere"))


def test_exp3_cylinder_permeable():
    """Exp 3: permeable cylinder packing -- kappa->0 GPA + monotone exchange."""
    _assert(*evaluate_exp3(_binary_or_skip(), "cylinder"))


def test_exp4_sphere_two_compartment():
    """Exp 4: impermeable sphere intra+extra -- intra==GPA + conservation + hindered extra."""
    _assert(*evaluate_exp4(_binary_or_skip(), "sphere", kappa=0.0))


def test_exp4_cylinder_two_compartment():
    """Exp 4: permeable cylinder intra+extra -- conservation + hindered extra."""
    _assert(*evaluate_exp4(_binary_or_skip(), "cylinder", kappa=KAPPA_ANCHOR))


# --- GPU equivalents + CPU/GPU parity (self-skip without CUDA) ----------------
def test_exp3_sphere_permeable_gpu():
    """Exp 3 on the GPU: permeable sphere packing (self-skips without CUDA)."""
    _assert(*evaluate_exp3(_binary_or_skip(), "sphere", gpu=1))


def test_exp3_cylinder_permeable_gpu():
    """Exp 3 on the GPU: permeable cylinder packing (self-skips without CUDA)."""
    _assert(*evaluate_exp3(_binary_or_skip(), "cylinder", gpu=1))


def test_exp4_sphere_two_compartment_gpu():
    """Exp 4 on the GPU: impermeable sphere intra+extra (self-skips without CUDA)."""
    _assert(*evaluate_exp4(_binary_or_skip(), "sphere", kappa=0.0, gpu=1))


def test_exp4_cylinder_two_compartment_gpu():
    """Exp 4 on the GPU: permeable cylinder intra+extra (self-skips without CUDA)."""
    _assert(*evaluate_exp4(_binary_or_skip(), "cylinder", kappa=KAPPA_ANCHOR, gpu=1))


def _parity(name, kappa, seeding):
    """CPU ≡ GPU on the identical permeable packing: total signal per-shell agree."""
    import pytest
    binary = _binary_or_skip()
    g = run_variant(binary, name, kappa, seeding, gpu=1)
    if g is None:
        pytest.skip("binary built without CUDA / no GPU available")
    c = run_variant(binary, name, kappa, seeding, gpu=0)
    nb = ~c.sc.b0_mask
    worst = metrics.per_shell_max_dev(c.total[nb], g.total[nb], c.sc.bvals[nb])
    assert worst <= 0.03, f"{name} permeable CPU/GPU parity off: worst per-shell |Δ|={worst:.4f}"


def test_exp3_sphere_cpu_gpu_parity():
    """Exp 3: permeable sphere CPU ≡ GPU (same answer)."""
    _parity("sphere", KAPPA_ANCHOR, "intra")


def test_exp3_cylinder_cpu_gpu_parity():
    """Exp 3: permeable cylinder CPU ≡ GPU (same answer)."""
    _parity("cylinder", KAPPA_ANCHOR, "intra")


# --- standalone CLI ----------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--gpu", type=int, default=0, choices=(0, 1))
    ap.add_argument("--tol", type=float, default=0.02)
    a = ap.parse_args()
    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1

    all_ok = True
    blocks = [
        ("exp3", evaluate_exp3, [("sphere",), ("cylinder",)]),
        ("exp4", evaluate_exp4, [("sphere", 0.0), ("cylinder", KAPPA_ANCHOR)]),
    ]
    for _, fn, argsets in blocks:
        for argset in argsets:
            ok, msgs = fn(binary, *argset, n=a.n, gpu=a.gpu, tol=a.tol)
            for m in msgs:
                print(m)
            if ok is None:
                print("  SKIP (no GPU)\n")
                continue
            all_ok = all_ok and ok
            print(f"  => {'PASS' if ok else 'FAIL'}\n")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
