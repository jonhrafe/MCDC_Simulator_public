"""Filesystem locations for the MC/DC test suite.

Single source of truth for "where is the repo / the binary / the fixtures", so
no test script has to re-derive them from ``__file__`` gymnastics. The binary
can be overridden with the ``MCDC_BIN`` environment variable (used by CTest,
which passes ``$<TARGET_FILE:MC-DC_Simulator>``, and by anyone testing a build
that does not live at the repo root).
"""
import os

# tests/mcdc_testkit/paths.py -> repo root is three levels up.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TESTS_DIR = os.path.join(REPO_ROOT, "tests")
ACCURACY_DIR = os.path.join(TESTS_DIR, "accuracy")
ANALYTIC_DIR = os.path.join(ACCURACY_DIR, "analytic")
REFERENCE_XV_DIR = os.path.join(ACCURACY_DIR, "reference_xv")
CROSSVERSION_DIR = os.path.join(TESTS_DIR, "crossversion")
CROSSVERSION_CONFS_DIR = os.path.join(CROSSVERSION_DIR, "confs")
# Locally-built old-version binaries (gitignored, machine-local cache) --
# see tests/tools/build_old_version.py.
OLD_BINARIES_DIR = os.path.join(CROSSVERSION_DIR, "old_binaries")


def binary_path(explicit=None):
    """Resolve the simulator binary: explicit arg > $MCDC_BIN > repo-root default."""
    return explicit or os.environ.get("MCDC_BIN") or os.path.join(REPO_ROOT, "MC-DC_Simulator")


def accuracy(*parts):
    """Path under tests/accuracy/ (fixtures, configs, lists, schemes)."""
    return os.path.join(ACCURACY_DIR, *parts)
