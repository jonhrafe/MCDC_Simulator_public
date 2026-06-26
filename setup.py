"""Build the MC/DC C++ simulator via its CMake project during `pip install`.

All package metadata lives in pyproject.toml. This file only adds the CMake build step:
it configures + builds the existing CMakeLists, then bundles the resulting binaries
(MC-DC_Simulator, dataSynth) inside the `mcdc` package so the console scripts can find them.

Only a C++17 compiler is required on the system (pip provides CMake via build-system.requires).
By default the build is portable; set MCDC_NATIVE=1 to compile with -march=native for this machine.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

HERE = Path(__file__).resolve().parent
BINARIES = ["MC-DC_Simulator", "dataSynth"]


class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        build_dir = Path(self.build_temp).resolve()
        build_dir.mkdir(parents=True, exist_ok=True)

        native = os.environ.get("MCDC_NATIVE", "0") == "1"
        cfg = [
            "cmake", "-S", str(HERE), "-B", str(build_dir),
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DMCDC_NATIVE_ARCH={'ON' if native else 'OFF'}",
        ]
        try:
            subprocess.run(cfg, check=True)
            subprocess.run(["cmake", "--build", str(build_dir), "-j"], check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            sys.exit(f"[mcdc] CMake build failed: {e}\n"
                     f"Ensure a C++17 compiler is installed (see instructions/compilation.md).")

        # The CMakeLists places the binaries at the project root (RUNTIME_OUTPUT_DIRECTORY).
        dst = Path(self.build_lib) / "mcdc" / "_bin"
        dst.mkdir(parents=True, exist_ok=True)
        for name in BINARIES:
            src = HERE / name
            if src.exists():
                shutil.copy2(src, dst / name)
                os.chmod(dst / name, 0o755)
            elif name == "MC-DC_Simulator":
                sys.exit(f"[mcdc] expected binary not produced: {src}")


setup(
    ext_modules=[CMakeExtension("mcdc._build")],
    cmdclass={"build_ext": CMakeBuild},
)
