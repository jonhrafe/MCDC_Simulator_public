"""Config parsing / rewriting and DWI signal loading.

These are the helpers every test script used to re-implement inline
(``read_exp_prefix`` / ``write_conf_with_prefix`` / the DWI ``fromfile``).
Centralised here so there is one correct implementation.
"""
import numpy as np


def read_conf_value(text, key):
    """First value of ``key`` in a ``.conf`` (whitespace-delimited), or None.

    Comment lines (``#``) are skipped. Only matches the token exactly as the
    first field of a line -- it does NOT reproduce the simulator's fuzzy
    Levenshtein matching (tests should use exact keys).
    """
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if parts and parts[0] == key:
            return parts[1] if len(parts) > 1 else ""
    return None


def read_exp_prefix(text):
    """Return the ``exp_prefix`` value, or raise if the config lacks one."""
    prefix = read_conf_value(text, "exp_prefix")
    if prefix is None:
        raise ValueError("no exp_prefix line in config")
    return prefix


def write_conf_override(text, overrides):
    """Return ``text`` with ``overrides`` applied as ``key value`` lines.

    A key already present (as the first token of a non-comment line) is
    *replaced* in place; a key not present is *inserted before* the mandatory
    ``<END>`` tag (the parser stops there), or appended if there is no ``<END>``.
    This one function covers both historical needs: rewriting ``exp_prefix`` and
    injecting a ``gpu`` flag before ``<END>``.
    """
    remaining = {k: v for k, v in overrides.items()}
    out = []
    for line in text.splitlines():
        s = line.strip()
        first = s.split()[0] if (s and not s.startswith("#")) else None
        if first is not None and first in remaining:
            out.append(f"{first} {remaining.pop(first)}")
            continue
        if s == "<END>" and remaining:
            out.extend(f"{k} {v}" for k, v in remaining.items())
            remaining = {}
        out.append(line)
    if remaining:  # no <END> in the file
        out.extend(f"{k} {v}" for k, v in remaining.items())
    return "\n".join(out) + "\n"


def load_signal(path):
    """Load a text signal (one float per line) as a flat float64 array.

    Used for the golden ``_DWI*.txt`` files (bit-level regression).
    """
    return np.loadtxt(path, dtype=np.float64).ravel()


def load_dwi(path, n_walkers=None):
    """Load a DWI signal (``.txt`` or ``.bfloat``) as a flat float64 array.

    The simulator writes the DWI as a SUM over walkers, so pass ``n_walkers`` to
    get the normalised S = signal / N (still needs dividing by the b0 mean to get
    S/S0 -- see :func:`mcdc_testkit.metrics.normalize_by_b0`). ``.bfloat`` is
    little-endian float32, one value per measurement.
    """
    if path.endswith(".bfloat") or path.endswith(".bin"):
        arr = np.fromfile(path, dtype="<f4").astype(np.float64)
    else:
        arr = np.loadtxt(path, dtype=np.float64).ravel()
    if n_walkers is not None:
        arr = arr / float(n_walkers)
    return arr
