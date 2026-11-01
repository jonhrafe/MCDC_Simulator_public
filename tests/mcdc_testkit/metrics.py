"""Signal normalisation and comparison metrics.

The DWI the simulator writes is a sum over walkers; the physically meaningful
quantity is S/S0. The pipeline is: ``load_dwi(path, N)`` -> ``normalize_by_b0``
-> per-shell means / NRMSE / ADC.
"""
import numpy as np

_TINY = np.finfo(float).tiny


def normalize_by_b0(dwi, b0_mask):
    """S/S0: divide the (already /N) signal by the mean of its b0 measurements."""
    b0 = np.asarray(dwi, dtype=np.float64)[b0_mask]
    denom = np.mean(b0)
    if not np.isfinite(denom) or denom == 0.0:
        raise ValueError("b0 mean is zero/non-finite; cannot normalise")
    return np.asarray(dwi, dtype=np.float64) / denom


def shells(bvals, decimals=0, unit=1e6):
    """Distinct b-shells (rounded, in ``bvals/unit`` units; default s/mm^2)."""
    return np.unique(np.round(np.asarray(bvals, dtype=np.float64) / unit, decimals))


def per_shell_mean(values, bvals, decimals=0, unit=1e6):
    """Mean of ``values`` grouped by b-shell -> ``{b_rounded: mean}``.

    ``b`` keys are in ``bvals/unit`` units (default s/mm^2), rounded to
    ``decimals`` so identical shells collapse to one key.
    """
    values = np.asarray(values, dtype=np.float64)
    key = np.round(np.asarray(bvals, dtype=np.float64) / unit, decimals)
    return {float(b): float(values[key == b].mean()) for b in np.unique(key)}


def per_shell_max_dev(a, b, bvals, decimals=0, unit=1e6):
    """Largest per-shell mean difference |mean_a - mean_b| over shared shells.

    The standard cross-hardware / cross-version parity metric: group both signals
    by b-shell and take the worst absolute difference of the shell means (robust to
    the per-measurement MC noise that a raw elementwise diff would pick up).
    """
    ma = per_shell_mean(a, bvals, decimals, unit)
    mb = per_shell_mean(b, bvals, decimals, unit)
    return max(abs(ma[k] - mb[k]) for k in ma if k in mb)


def nrmse(run, ref, norm="range"):
    """Normalised RMSE of ``run`` vs ``ref``.

    ``norm='range'`` divides by peak-to-peak of ``ref`` (robust when the signal
    spans ~0..1); ``'mean'`` divides by mean(|ref|); ``'none'`` returns plain RMSE.
    """
    run, ref = np.asarray(run, float), np.asarray(ref, float)
    rmse = float(np.sqrt(np.mean((run - ref) ** 2)))
    if norm == "range":
        denom = float(np.ptp(ref))
    elif norm == "mean":
        denom = float(np.mean(np.abs(ref)))
    elif norm == "none":
        return rmse
    else:
        raise ValueError(f"unknown norm {norm!r}")
    return rmse / max(denom, _TINY)


def max_abs_dev(run, ref):
    return float(np.max(np.abs(np.asarray(run, float) - np.asarray(ref, float))))


def max_rel_dev(run, ref):
    run, ref = np.asarray(run, float), np.asarray(ref, float)
    return float(np.max(np.abs(run - ref) / np.maximum(np.abs(ref), _TINY)))


def within_tol(run, ref, atol, rtol):
    """Golden-style elementwise check ``|run-ref| <= atol + rtol*|ref|``.

    Returns ``(ok, n_bad, max_abs, max_rel)``.
    """
    run, ref = np.asarray(run, float), np.asarray(ref, float)
    abs_dev = np.abs(run - ref)
    ok_mask = abs_dev <= (atol + rtol * np.abs(ref))
    n_bad = int(np.count_nonzero(~ok_mask))
    rel = abs_dev / np.maximum(np.abs(ref), _TINY)
    return (n_bad == 0, n_bad, float(abs_dev.max()), float(rel.max()))


def adc(ss0, bvals):
    """Apparent diffusion coefficient ADC = -ln(S/S0)/b (NaN where S/S0<=0)."""
    ss0 = np.asarray(ss0, dtype=np.float64).copy()
    ss0[ss0 <= 0] = np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.log(ss0) / np.asarray(bvals, dtype=np.float64)


def binomial_sigma(p, n):
    """Std-dev of an empirical probability estimate from ``n`` Bernoulli trials."""
    return float(np.sqrt(max(p * (1.0 - p), 1e-15) / max(n, 1)))
