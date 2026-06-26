"""Console entry points that locate and run the bundled compiled binaries."""
import os
import sys
from pathlib import Path


def _binary(name):
    """Return the path to a bundled binary, or fall back to PATH."""
    here = Path(__file__).resolve().parent / "_bin"
    for cand in (here / name, here / f"{name}.exe"):
        if cand.exists():
            return str(cand)
    found = __import__("shutil").which(name)
    if found:
        return found
    sys.exit(f"[mcdc] '{name}' binary not found. Reinstall with `pip install .` "
             f"from the repository root, or build it with CMake (see instructions/compilation.md).")


def _run(name):
    argv = sys.argv[1:]
    try:
        os.execv(_binary(name), [name, *argv])          # replace the process (POSIX)
    except (AttributeError, OSError):
        import subprocess
        sys.exit(subprocess.call([_binary(name), *argv]))


def main():
    """`mcdc <config.conf>` -> MC-DC_Simulator."""
    _run("MC-DC_Simulator")


def datasynth():
    """`mcdc-datasynth ...` -> dataSynth."""
    _run("dataSynth")
