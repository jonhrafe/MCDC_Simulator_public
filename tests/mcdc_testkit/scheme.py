"""PGSE (Stejskal-Tanner) scheme parsing and b-value computation.

A STEJSKALTANNER scheme row is 7 whitespace-separated floats:
``gx gy gz  G  Delta  delta  TE`` (gradient direction, gradient amplitude in
T/m, big/small delta and echo time in s -- all SI). b0 rows have ``G == 0``.
"""
from collections import namedtuple

import numpy as np

# Gyromagnetic ratio of the proton, rad/(s*T). Matches the simulator's
# giro = 267.51525e3 rad/(ms*T) (see src/constants.h) and test_gpu_free_diffusion.
GAMMA = 2.6751525e8


def bvalue(G, Delta, delta):
    """PGSE b-value b = (gamma*G*delta)^2 * (Delta - delta/3), SI (s/m^2)."""
    return (GAMMA * G * delta) ** 2 * (Delta - delta / 3.0)


Scheme = namedtuple("Scheme", "gdirs G Delta delta TE bvals b0_mask n")


def parse_scheme(path):
    """Parse a STEJSKALTANNER scheme file into arrays (SI units).

    ``bvals`` is in s/m^2 (divide by 1e6 for s/mm^2). ``b0_mask`` marks the
    ``G == 0`` rows. Header / non-numeric lines are skipped.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 7:
                try:
                    rows.append([float(v) for v in parts[:7]])
                except ValueError:
                    continue  # header (e.g. "VERSION: STEJSKALTANNER")
    if not rows:
        raise ValueError(f"no PGSE rows parsed from {path}")
    a = np.asarray(rows, dtype=np.float64)
    gdirs, G, Delta, delta, TE = a[:, 0:3], a[:, 3], a[:, 4], a[:, 5], a[:, 6]
    bvals = bvalue(G, Delta, delta)
    return Scheme(gdirs, G, Delta, delta, TE, bvals, G == 0.0, len(a))


def angle_to_axis_deg(gdirs, axis=(0.0, 0.0, 1.0)):
    """Angle (deg, 0..90) between each gradient direction and ``axis``.

    Uses ``|cos|`` so antiparallel directions map to the same angle. For a
    z-cylinder plot the conventional x-axis is ``90 - angle_to_axis_deg`` (0 =
    perpendicular to the axis, 90 = along it).
    """
    axis = np.asarray(axis, dtype=np.float64)
    axis = axis / np.linalg.norm(axis)
    cos = np.abs(np.asarray(gdirs, dtype=np.float64) @ axis)
    return np.degrees(np.arccos(np.clip(cos, 0.0, 1.0)))
