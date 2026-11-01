"""MC/DC test kit -- shared helpers for the validation & regression suite.

One import for the plumbing every test needs: locate the binary, run a config,
load and normalise the DWI, compare against an oracle. Kept dependency-light
(only numpy) so the smoke / fast-regression tier needs nothing else.

    from mcdc_testkit import io, run, scheme, metrics, oracles, paths

Add ``tests/`` to ``sys.path`` first (see ``tests/conftest.py`` for pytest, or
the ``__main__`` block in each standalone script).
"""
from . import io, metrics, oracles, paths, provenance, run, scheme, substrate

__all__ = ["io", "run", "scheme", "metrics", "oracles", "provenance", "paths",
           "substrate"]
