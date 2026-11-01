#!/usr/bin/env python3
"""
Extra-Cellular Volume Fraction (ECVF) from a ``separate_signals`` run.

A weak, analytic-free oracle for mesh substrates that have no closed-form signal.
When a config sets ``separate_signals 1`` the simulator splits the DWI by the
compartment each walker was classified into at seeding and writes
``<prefix>_DWI_intra.txt`` and ``<prefix>_DWI_extra.txt`` (one summed value per
scheme measurement). At a b0 measurement (``G == 0``) the phase is zero, so every
walker contributes ``exp(-TE/T2)`` -- with ``T2 = inf`` that is exactly 1, so the
b0 signal is simply the WALKER COUNT in that compartment. With uniform seeding
those counts are proportional to the compartment VOLUMES, hence

    ECVF = mean(extra[b0]) / (mean(intra[b0]) + mean(extra[b0]))         (geometric)

is the extra-cellular volume fraction and ``ICVF = 1 - ECVF``.

T2-WEIGHTING CAVEAT
-------------------
If the intra and extra compartments have DIFFERENT T2 (``t2_intra != t2_extra``),
the b0 signals are ``count * exp(-TE/T2)``-weighted, so the raw ratio above is a
T2-BIASED ECVF, not the geometric one. Given the (known) T2s and TE you can divide
the weighting out to recover the clean geometric fraction:

    count = b0_signal * exp(+TE/T2)          # undo the exp(-TE/T2) weighting

Pass ``--t2-intra`` and ``--t2-extra`` (both, in seconds) to apply this correction;
``--te`` overrides the echo time (default: the TE of the b0 rows read from the
scheme, corrected per-measurement so multi-TE b0 sets are handled). With a single
compartment T2 (or ``T2 = inf``, the usual impermeable-mesh setup) no correction is
needed and the raw ratio already IS the geometric ECVF.

Pure post-processing: numpy only, no simulator, no GPU.

    python3 tests/tools/compute_ecvf.py <exp_prefix> <scheme.scheme>
    python3 tests/tools/compute_ecvf.py run --scheme s.scheme --t2-intra 0.04 --t2-extra 0.08
"""
import argparse
import os
import sys
from collections import namedtuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # tests/
from mcdc_testkit import io, scheme as kscheme  # noqa: E402

import numpy as np  # noqa: E402

# ecvf/icvf: geometric fractions AFTER any T2 correction; ecvf_raw: the uncorrected
# b0 ratio. intra_b0/extra_b0: raw b0 means (walker-count units). n_b0: #b0 rows.
Result = namedtuple("EcvfResult", "ecvf ecvf_raw icvf intra_b0 extra_b0 n_b0 corrected")


def _find_signal(prefix, suffix):
    """Path of ``<prefix>_DWI<suffix>.txt`` (exact first, then a glob fallback)."""
    exact = f"{prefix}_DWI{suffix}.txt"
    if os.path.exists(exact):
        return exact
    import glob
    cands = glob.glob(f"{prefix}*_DWI{suffix}.txt")
    if cands:
        return max(cands, key=os.path.getmtime)
    raise FileNotFoundError(
        f"missing {exact} -- run the config with `separate_signals 1` so the "
        f"simulator writes _DWI_intra/_DWI_extra.")


def load_b0(prefix, scheme_path):
    """Return ``(intra_b0, extra_b0, te_b0)`` arrays over the b0 rows of the scheme.

    ``intra_b0``/``extra_b0`` are the raw per-measurement b0 signals (summed over
    walkers); ``te_b0`` is the echo time of each b0 row (for the T2 correction).
    """
    intra = io.load_signal(_find_signal(prefix, "_intra"))
    extra = io.load_signal(_find_signal(prefix, "_extra"))
    sc = kscheme.parse_scheme(scheme_path)
    if intra.size != sc.n or extra.size != sc.n:
        raise AssertionError(
            f"signal length ({intra.size} intra / {extra.size} extra) != scheme "
            f"rows ({sc.n}); wrong scheme for this run?")
    if not sc.b0_mask.any():
        raise ValueError(f"no b0 (G==0) rows in {scheme_path}; cannot split by b0")
    return intra[sc.b0_mask], extra[sc.b0_mask], sc.TE[sc.b0_mask]


def compute_ecvf(intra_b0, extra_b0, te_b0=None, t2_intra=None, t2_extra=None, te=None):
    """ECVF from the b0 intra/extra split, optionally undoing T2 weighting.

    ``intra_b0``/``extra_b0``: b0 signals (arrays or scalars). Correction is applied
    only when BOTH ``t2_intra`` and ``t2_extra`` are given (the ratio depends on the
    T2 *contrast*); ``te`` (scalar) overrides ``te_b0`` (per-row, from the scheme).
    """
    intra_b0 = np.atleast_1d(np.asarray(intra_b0, dtype=np.float64))
    extra_b0 = np.atleast_1d(np.asarray(extra_b0, dtype=np.float64))
    mi_raw, me_raw = float(intra_b0.mean()), float(extra_b0.mean())
    denom_raw = mi_raw + me_raw
    if denom_raw <= 0.0:
        raise ValueError(f"non-positive total b0 signal ({denom_raw}); nothing seeded?")
    ecvf_raw = me_raw / denom_raw

    corrected = False
    if t2_intra is not None and t2_extra is not None:
        if te is not None:
            te_arr = np.full(intra_b0.shape, float(te), dtype=np.float64)
        elif te_b0 is not None:
            te_arr = np.atleast_1d(np.asarray(te_b0, dtype=np.float64))
        else:
            raise ValueError("T2 correction needs --te or a scheme TE (te_b0)")
        # Undo exp(-TE/T2): recover the geometric (count-proportional) b0 weight.
        mi = float((intra_b0 * np.exp(te_arr / float(t2_intra))).mean())
        me = float((extra_b0 * np.exp(te_arr / float(t2_extra))).mean())
        ecvf, corrected = me / (mi + me), True
    else:
        ecvf = ecvf_raw

    return Result(float(ecvf), float(ecvf_raw), float(1.0 - ecvf),
                  mi_raw, me_raw, int(intra_b0.size), corrected)


def ecvf_from_prefix(prefix, scheme_path, t2_intra=None, t2_extra=None, te=None):
    """Convenience: load the b0 split for ``prefix`` and compute the ECVF."""
    intra_b0, extra_b0, te_b0 = load_b0(prefix, scheme_path)
    return compute_ecvf(intra_b0, extra_b0, te_b0=te_b0,
                        t2_intra=t2_intra, t2_extra=t2_extra, te=te)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="exp_prefix of the run (reads <prefix>_DWI_intra.txt/_extra.txt)")
    ap.add_argument("scheme", help="scheme file (for the b0 mask + per-b0 TE)")
    ap.add_argument("--t2-intra", type=float, default=None, help="intra T2 [s] (enables T2 correction)")
    ap.add_argument("--t2-extra", type=float, default=None, help="extra T2 [s] (enables T2 correction)")
    ap.add_argument("--te", type=float, default=None, help="echo time [s] override (default: scheme b0 TE)")
    a = ap.parse_args()

    if (a.t2_intra is None) != (a.t2_extra is None):
        print("NOTE: T2 correction needs BOTH --t2-intra and --t2-extra; "
              "reporting the raw (T2-weighted) ECVF.", file=sys.stderr)

    try:
        r = ecvf_from_prefix(a.prefix, a.scheme,
                             t2_intra=a.t2_intra, t2_extra=a.t2_extra, te=a.te)
    except (FileNotFoundError, AssertionError, ValueError) as e:
        print(f"FAIL: {e}", file=sys.stderr)
        return 1

    print(f"b0 split ({r.n_b0} b0 rows): intra={r.intra_b0:.4g}  extra={r.extra_b0:.4g}")
    if r.corrected:
        print(f"ECVF (raw, T2-weighted) = {r.ecvf_raw:.4f}")
        print(f"ECVF (T2-corrected)     = {r.ecvf:.4f}   ICVF = {r.icvf:.4f}")
    else:
        print(f"ECVF = {r.ecvf:.4f}   ICVF = {r.icvf:.4f}   (no T2 correction)")

    if not (0.0 < r.ecvf < 1.0):
        print(f"FAIL: ECVF {r.ecvf:.4f} not in (0,1).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
