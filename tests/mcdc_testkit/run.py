"""Run the simulator and detect failure / GPU-skip from its output.

Mirrors the error handling every script grew independently: the simulator
signals config-validation failure via ``assert()``, which is a no-op under
``-DNDEBUG`` (release builds), so an exit-code-only check misses it -- we also
scan stdout for ``[ERROR]`` / ``Assertion`` / ``error:``.
"""
import os
import subprocess
import sys

from . import io

# Substrings the GPU path prints when CUDA is absent / no device -> the test
# should SKIP (not fail), so CPU-only machines still pass the suite.
GPU_SKIP_MARKERS = ("WITHOUT CUDA", "CUDA error", "no CUDA-capable device")


def is_gpu_skip(log):
    return any(m in log for m in GPU_SKIP_MARKERS)


def has_error(log):
    """True if the simulator's stdout reports an error even on a 0 exit code."""
    low = log.lower()
    return ("[ERROR]" in log) or ("Assertion" in log) or (
        "error:" in low and "no error" not in low)


def run_sim(binary, conf_path, workdir, timeout=600):
    """Run ``binary --conf conf_path`` with cwd=workdir; return CompletedProcess.

    cwd matters: relative paths inside the config (``scheme_file``,
    ``cylinders_list`` ...) resolve against ``workdir`` (normally the repo root).
    """
    return subprocess.run(
        [binary, "--conf", conf_path], cwd=workdir,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)


def raise_on_error(res):
    """Raise RuntimeError if the run failed (non-zero exit OR error in output)."""
    log = res.stdout or ""
    if res.returncode != 0 or has_error(log):
        sys.stderr.write(log)
        raise RuntimeError(
            f"simulator failed (exit {res.returncode}; error_in_output={has_error(log)})")


def run_case(binary, conf_text, out_dir, base, workdir, overrides=None,
             gpu=None, timeout=600):
    """Materialise a config, run it, and return the ``exp_prefix`` path.

    ``exp_prefix`` is forced to ``out_dir/base`` (outputs land in the temp dir);
    ``gpu`` (0/1) and any ``overrides`` (dict of ``key -> value``) are applied via
    :func:`io.write_conf_override`. Relative fixture paths resolve against
    ``workdir``. Returns the prefix on success, or ``None`` if the GPU path
    skipped (CUDA absent). Raises on a real failure.
    """
    ov = {"exp_prefix": os.path.join(out_dir, base)}
    if gpu is not None:
        ov["gpu"] = int(gpu)
    if overrides:
        ov.update(overrides)
    conf_path = os.path.join(out_dir, base + ".conf")
    with open(conf_path, "w") as fh:
        fh.write(io.write_conf_override(conf_text, ov))
    res = run_sim(binary, conf_path, workdir, timeout=timeout)
    log = res.stdout or ""
    if gpu and is_gpu_skip(log):
        return None
    raise_on_error(res)
    return os.path.join(out_dir, base)
