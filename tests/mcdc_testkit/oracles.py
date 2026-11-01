"""Analytical / frozen reference signals (the "oracles" tests compare against).

Two kinds live here:

* **closed-form** oracles computed live in pure numpy (no heavy deps) -- e.g.
  free diffusion ``exp(-bD)``;
* **frozen** oracles loaded from committed ``.npy`` arrays under
  ``tests/accuracy/analytic/`` -- e.g. the dmipy restricted-cylinder (van
  Gelderen) curve, produced once by ``tests/tools/gen_analytic.py`` (which needs
  dmipy) so the test-time path needs only numpy. See that tool and
  ``tests/validation/README.md`` for provenance and validity limits.
"""
import os

import numpy as np

from . import paths
from .scheme import GAMMA, angle_to_axis_deg


def free_exp(bvals, D):
    """Free (Gaussian) diffusion signal S/S0 = exp(-b*D). ``b`` and ``D`` SI."""
    return np.exp(-np.asarray(bvals, dtype=np.float64) * D)


# --------------------------------------------------------------------------- #
#  Restricted-diffusion Gaussian Phase Approximation (pure numpy)             #
#                                                                             #
#  Closed-form PGSE signal for restricted diffusion inside an impermeable     #
#  cylinder (van Gelderen) / sphere (Neuman-Murday-Cotts). This is the SAME   #
#  model dmipy's C4Cylinder.../S4Sphere...GaussianPhaseApproximation compute, #
#  re-implemented here so the polydisperse (packing) oracle needs only numpy  #
#  and is evaluated live from a packing's own radii -- there is no single     #
#  frozen curve to store because it is specific to the exact packing.         #
#                                                                             #
#  VALIDITY (identical to the dmipy path): infinite (cyl) / closed (sph),     #
#  IMPERMEABLE geometry, small b (<~5000 s/mm^2). Validates restricted        #
#  diffusion, NOT permeability or high b.                                     #
# --------------------------------------------------------------------------- #

# Zeros of J1'(x) (cylinder) and of the spherical-Bessel j1'(x) (sphere, the
# reflecting-sphere Neumann eigenvalues). 15 roots -> the tail terms fall off
# ~root^-6, so this is far past convergence at these b-values / diffusion times.
_ROOTS_CYL = np.array([
    1.8411837813, 5.3314427735, 8.5363163663, 11.7060049026, 14.8635886339,
    18.0155278622, 21.1643698591, 24.3113268578, 27.4570505710, 30.6019229552,
    33.7461812201, 36.8899866977, 40.0334439409, 43.1766274197, 46.3195966130])
_ROOTS_SPH = np.array([
    2.0815759779, 5.9403699890, 9.2058418752, 12.4044450219, 15.5792364104,
    18.7426455848, 21.8996964795, 25.0528252115, 28.2033610039, 31.3520917266,
    34.4995149214, 37.6459603231, 40.7916552313, 43.9367614714, 47.0813974122])


def _gpd_perp(Gperp, Delta, delta, D, R, roots, geom_const):
    """Perpendicular-gradient GPA attenuation E = exp(-2 gamma^2 Gperp^2 S).

    All args SI; ``Gperp``/``Delta``/``delta`` are per-measurement arrays (n,),
    ``roots`` is (M,), ``geom_const`` is 1 (cylinder: root^2 R^2 - 1) or 2
    (sphere: root^2 R^2 - 2). Returns E of shape (n,).

        S = sum_m 1/(a_m^2 (a_m^2 R^2 - c)) *
              [ 2 delta / B_m
                - (2 + e^{-B(Delta-delta)} - 2 e^{-B delta} - 2 e^{-B Delta}
                     + e^{-B(Delta+delta)}) / B_m^2 ]
        a_m = root_m / R,   B_m = a_m^2 D
    """
    Gperp = np.asarray(Gperp, float)[:, None]
    Dl = np.asarray(Delta, float)[:, None]
    dl = np.asarray(delta, float)[:, None]
    a2 = (roots / R) ** 2                      # (M,)
    B = a2 * D                                 # (M,) rate 1/s
    pref = 1.0 / (a2 * (a2 * R ** 2 - geom_const))          # (M,)
    bracket = (2.0 * dl / B
               - (2.0 + np.exp(-B * (Dl - dl)) - 2.0 * np.exp(-B * dl)
                  - 2.0 * np.exp(-B * Dl) + np.exp(-B * (Dl + dl))) / B ** 2)
    S = (pref[None, :] * bracket).sum(axis=1)              # (n,)
    return np.exp(-2.0 * GAMMA ** 2 * Gperp[:, 0] ** 2 * S)


def sphere_neuman(scheme, R, D):
    """Restricted-sphere GPA signal S/S0 for radius ``R`` (m), diffusivity ``D``.

    Isotropic: depends only on (G, Delta, delta), so every direction in a shell
    shares one value. Returns a full-length array (one per scheme row; b0 ~ 1).
    """
    return _gpd_perp(scheme.G, scheme.Delta, scheme.delta, D, R, _ROOTS_SPH, 2)


def cylinder_vangelderen(scheme, R, D, axis=(0.0, 0.0, 1.0)):
    """Restricted-cylinder GPA signal S/S0 for radius ``R`` (m), diffusivity ``D``.

    The cylinder axis defaults to +z (the gamma packer lays all cylinders along
    z). The signal factorises into a free PARALLEL part exp(-b cos^2(theta) D)
    and a restricted PERPENDICULAR van Gelderen part evaluated at Gperp =
    G sin(theta), where theta is the gradient-to-axis angle. Full-length (b0 ~ 1).
    """
    theta = np.radians(angle_to_axis_deg(scheme.gdirs, axis))   # 0=|| axis, 90=perp
    sin_t, cos_t = np.sin(theta), np.cos(theta)
    e_perp = _gpd_perp(scheme.G * sin_t, scheme.Delta, scheme.delta, D, R, _ROOTS_CYL, 1)
    e_par = np.exp(-scheme.bvals * cos_t ** 2 * D)
    return e_par * e_perp


def packing_signal(scheme, radii, D, geometry, axis=(0.0, 0.0, 1.0)):
    """Volume-weighted restricted signal of a polydisperse packing (S/S0).

    ``geometry`` is 'sphere' or 'cylinder'. Each obstacle contributes its
    single-obstacle GPA signal weighted by its intra volume: spheres by r^3
    (V=4/3 pi r^3), cylinders by r^2 (V=pi r^2 L, common L). This mirrors the
    Monte-Carlo signal when walkers are seeded UNIFORMLY in the intra volume
    (``ini_walkers_pos intra``): the walker count per obstacle is ~proportional
    to its volume, so MC ~ the volume-weighted average of per-obstacle signals.
    Returns a full-length array (one per scheme row; b0 ~ 1).
    """
    radii = np.asarray(radii, float)
    if geometry == "sphere":
        w = radii ** 3
        sig = np.stack([sphere_neuman(scheme, r, D) for r in radii])       # (K, n)
    elif geometry == "cylinder":
        w = radii ** 2
        sig = np.stack([cylinder_vangelderen(scheme, r, D, axis) for r in radii])
    else:
        raise ValueError(f"geometry must be 'sphere' or 'cylinder', got {geometry!r}")
    return (w[:, None] * sig).sum(axis=0) / w.sum()


def karger(bvals, f_a, D_a, D_b, k, k_ba=None):
    r"""Kärger two-compartment exchange signal S/S0 (all SI: b in s/m^2,
    D in m^2/s, k in 1/s).

    Two exchanging Gaussian compartments ``a`` and ``b`` with equilibrium
    fractions ``f_a`` and ``f_b = 1 - f_a``, intrinsic diffusivities ``D_a`` /
    ``D_b`` and an exchange process with rates ``k_a`` (a->b) and ``k_b``
    (b->a). Detailed balance requires ``f_a k_a = f_b k_b``.

    ``k`` is the TOTAL exchange rate (apparent exchange rate AXR = k_a + k_b);
    the rates are then split as ``k_a = k f_b``, ``k_b = k f_a`` (which satisfies
    detailed balance and gives an exchange time tau_ex = 1/k). Pass ``k_ba`` to
    give the two rates explicitly instead: then ``k`` == k_a (a->b) and
    ``k_ba`` == k_b (b->a).

    The echo attenuation is the bi-exponential (Kärger 1988)::

        S/S0 = c1 exp(-b D1') + c2 exp(-b D2')

    with the apparent diffusivities being the eigenvalues of the 2x2 Kärger
    matrix ``L = [[D_a + k_a/b, -k_b/b], [-k_a/b, D_b + k_b/b]]``::

        D1',D2' = 1/2 [ D_a + D_b + (k_a+k_b)/b ]
                  -/+ 1/2 sqrt{ [D_a - D_b + (k_a-k_b)/b]^2 + 4 k_a k_b / b^2 }

    and weights ``c1 = (D2' - Dbar)/(D2' - D1')``, ``c2 = 1 - c1``, where
    ``Dbar = f_a D_a + f_b D_b`` is the population-weighted mean diffusivity.
    (Both derived by diagonalising L; c1 + c2 = 1 so S(b=0) = 1. See the
    self-test in tests/validation/test_karger_exchange.py, which cross-checks
    this closed form against a direct matrix exponential and the no-/fast-
    exchange limits.)

    REGIME / VALIDITY. This is the Gaussian-phase (narrow-pulse) Kärger model:
    each compartment is assumed *Gaussian* (free/unrestricted over the diffusion
    time) with a well-defined intrinsic D, exchanging in the slow/intermediate
    regime. It does NOT describe restricted (non-Gaussian) intra-compartment
    diffusion (e.g. a small sphere), where the apparent compartment diffusivity
    is time-dependent -- see test_karger_exchange.py for why the MC permeable-
    sphere fixture is validated by geometry/limits, not a raw Kärger fit.
    """
    b = np.asarray(bvals, dtype=np.float64)
    f_b = 1.0 - f_a
    if k_ba is None:
        k_a, k_b = k * f_b, k * f_a          # AXR split obeying detailed balance
    else:
        k_a, k_b = k, k_ba
    Dbar = f_a * D_a + f_b * D_b
    out = np.ones_like(b)                     # S(b=0) = 1 exactly
    nz = b > 0.0
    bb = b[nz]
    half_sum = 0.5 * (D_a + D_b + (k_a + k_b) / bb)
    root = 0.5 * np.sqrt((D_a - D_b + (k_a - k_b) / bb) ** 2 + 4.0 * k_a * k_b / bb ** 2)
    D1 = half_sum - root                      # smaller apparent diffusivity
    D2 = half_sum + root                      # larger apparent diffusivity
    # weights from diagonalisation (guard the degenerate D1==D2 case, e.g. D_a==D_b, k==0)
    denom = np.where(np.abs(D2 - D1) > 0.0, D2 - D1, 1.0)
    c1 = np.where(np.abs(D2 - D1) > 0.0, (D2 - Dbar) / denom, f_a)
    c2 = 1.0 - c1
    out[nz] = c1 * np.exp(-bb * D1) + c2 * np.exp(-bb * D2)
    return out


def load_analytic(name, analytic_dir=None):
    """Load a frozen analytic curve ``<analytic_dir>/<name>.npy`` (float64).

    Raises FileNotFoundError with a pointer to the generator when the curve has
    not been frozen yet (e.g. a fresh checkout, or an env without dmipy).
    """
    analytic_dir = analytic_dir or paths.ANALYTIC_DIR
    path = os.path.join(analytic_dir, name + ".npy")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"frozen analytic curve not found: {path}\n"
            f"Generate it with:  python3 tests/tools/gen_analytic.py  "
            f"(requires dmipy; see tests/validation/README.md).")
    return np.load(path).astype(np.float64)
