"""pytest configuration for the MC/DC validation suite.

Puts ``tests/`` on ``sys.path`` so ``import mcdc_testkit`` works from any test
module, and registers the tier/selection markers so ``pytest -m ...`` and the
CTest label wrappers agree.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # tests/

MARKERS = {
    "gpu": "requires a CUDA build + GPU; self-skips otherwise",
    "slow": "large-N analytical validation (minutes on CPU, seconds on GPU)",
    "crossversion": "compares against frozen previous-version reference signals",
    "heavy": "large-mesh tier (on-demand)",
}


def pytest_configure(config):
    for name, desc in MARKERS.items():
        config.addinivalue_line("markers", f"{name}: {desc}")
