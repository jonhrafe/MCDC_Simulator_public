#!/usr/bin/env python3
"""
Freeze a cross-version reference DWI signal (DEV-ONLY, Phase C).

Runs a chosen simulator build (the current one, or an old v1.5 / v2.0
binary) on a given config at a FIXED positive seed and FIXED num_process, then
copies the resulting ``_DWI.txt`` into
``tests/accuracy/reference_xv/<version>/<case>_DWI.txt`` and records a manifest
entry in that directory's ``MANIFEST.json``. The run-time cross-version tier then
compares the current build against these frozen signals with numpy only.

    python3 tests/tools/gen_reference.py --version current \
        --config tests/accuracy/free_diffusion.conf \
        --seed 12345 --num-process 1 --gpu 0

DETERMINISM GUARD -- this tool REFUSES to run without a fixed positive ``--seed``
AND a ``--num-process``: a MC run is reproducible only for both fixed together
(the per-worker RNG stream depends on the thread count). Even so, results are
NOT bit-exact ACROSS versions (different RNGs); the frozen signals are compared
statistically (per-shell means), never for equality -- see
tests/crossversion/README.md.

PER-VERSION CONFIGS -- the old binaries use DIFFERENT ``.conf`` / unit
conventions (v1.5 & v2.0 need ``scale_from_stu 1`` and a different cylinder-list
scale). Do NOT assume the current SI config runs on v1.5: pass the matching
per-version config (tests/crossversion/confs/<version>/). The
binary and its argstyle (``--conf FILE`` vs positional ``FILE``) come from
tests/crossversion/binaries.toml unless overridden with ``--bin`` / ``--argstyle``.

NEVER freezes trajectories -- ``write_traj_file`` is forced to 0 and only the
tiny ``_DWI.txt`` signal (~2 KB) is copied.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # tests/
from mcdc_testkit import io, paths, provenance, run  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402

import build_old_version  # noqa: E402 (sibling module, tests/tools/)
import numpy as np  # noqa: E402

BINARIES_TOML = os.path.join(paths.TESTS_DIR, "crossversion", "binaries.toml")


def load_versions(path=BINARIES_TOML):
    """`version -> {path, argstyle, ...}` from binaries.toml, or {} if unreadable.

    Uses the python 3.11+ stdlib ``tomllib``; on older pythons (or a malformed /
    missing file) returns an empty registry so ``--bin`` still works.
    """
    try:
        import tomllib
    except ModuleNotFoundError:  # py < 3.11
        sys.stderr.write("[warn] tomllib unavailable (python < 3.11); use --bin.\n")
        return {}
    try:
        with open(path, "rb") as fh:
            return tomllib.load(fh).get("versions", {})
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"[warn] could not read {path}: {exc}\n")
        return {}


def parse_version_string(info_path):
    """Extract 'Software Version' from a ``*_simulation_info.txt``, or None.

    The line reads ``[INFO] Software Version: -------2.2.0`` -- strip the label
    and the dashed fill. This string is UNRELIABLE across builds (e.g. the
    current build self-reports 2.1.3.beta); the manifest also pins the git
    commit. Captured only for reference.
    """
    try:
        with open(info_path) as fh:
            for line in fh:
                if "Software Version" in line:
                    tail = line.split("Software Version", 1)[1].lstrip(": ").strip()
                    return tail.lstrip("-").strip() or None
    except OSError:
        return None
    return None


def resolve_binary(version, explicit, versions_cfg):
    """(binary_path, argstyle) from --bin/--argstyle or the registry, or (None, ...)."""
    entry = versions_cfg.get(version, {})
    argstyle = entry.get("argstyle", "conf")
    if explicit:
        return explicit, argstyle
    if version == "current":
        return paths.binary_path(), argstyle  # honours $MCDC_BIN
    p = entry.get("path")
    if not p:
        return None, argstyle
    return (p if os.path.isabs(p) else os.path.join(paths.REPO_ROOT, p)), argstyle


def resolve_scheme(conf_text, workdir):
    """Absolute path to the config's ``scheme_file`` (resolved against workdir)."""
    rel = io.read_conf_value(conf_text, "scheme_file")
    if not rel:
        return None
    return rel if os.path.isabs(rel) else os.path.join(workdir, rel)


def run_binary(binary, conf_path, workdir, argstyle, timeout):
    """Invoke the binary the way this version expects (--conf vs positional)."""
    cmd = [binary, "--conf", conf_path] if argstyle == "conf" else [binary, conf_path]
    return subprocess.run(cmd, cwd=workdir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=timeout)


def _int_or_none(text, key):
    v = io.read_conf_value(text, key)
    if v is None:
        return None
    try:
        return int(float(v))
    except ValueError:
        return v  # keep the raw string if it is not numeric


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", required=True,
                    help="registry key (current / v1.5 / v2.0); "
                         "names the reference_xv/<version>/ subdir")
    ap.add_argument("--config", required=True,
                    help="per-version .conf to run (old versions need their own -- "
                         "see tests/crossversion/confs/)")
    # Determinism guard: intentionally NOT argparse-required so we can emit a
    # helpful refusal (and a clean exit code) rather than an argparse usage error.
    ap.add_argument("--seed", type=int, default=None,
                    help="fixed POSITIVE RNG seed (required; determinism guard)")
    ap.add_argument("--num-process", type=int, default=None,
                    help="fixed worker count (required; determinism guard)")
    ap.add_argument("--gpu", type=int, choices=(0, 1), default=None,
                    help="override the gpu flag (0 = force CPU); omit to leave the config's")
    ap.add_argument("--case", default=None,
                    help="frozen basename (default: the config's exp_prefix basename)")
    ap.add_argument("--bin", default=None,
                    help="explicit binary path (overrides binaries.toml)")
    ap.add_argument("--argstyle", choices=("conf", "positional"), default=None,
                    help="how the binary takes the config (default: from binaries.toml)")
    ap.add_argument("--workdir", default=paths.REPO_ROOT,
                    help="CWD for the run (relative paths in the conf resolve here)")
    ap.add_argument("--out-root", default=paths.REFERENCE_XV_DIR,
                    help="root for reference_xv/<version>/")
    ap.add_argument("--timeout", type=int, default=1800, help="per-run timeout [s]")
    a = ap.parse_args()

    # --- determinism guard -------------------------------------------------
    if a.seed is None or a.seed <= 0:
        print("REFUSE: --seed must be a fixed POSITIVE integer. A MC run is only "
              "reproducible for a fixed seed (default seed -1 is random). "
              "Cross-version references must be deterministic.", file=sys.stderr)
        return 2
    if a.num_process is None or a.num_process <= 0:
        print("REFUSE: --num-process is required (a fixed positive worker count). "
              "The per-worker RNG stream depends on the thread count, so a run is "
              "reproducible only for a fixed seed AND a fixed num_process.", file=sys.stderr)
        return 2

    versions_cfg = load_versions()
    binary, argstyle = resolve_binary(a.version, a.bin, versions_cfg)
    if a.argstyle:
        argstyle = a.argstyle
    if not binary:
        print(f"ERROR: no binary for version {a.version!r}: pass --bin or add it to "
              f"{os.path.relpath(BINARIES_TOML, paths.REPO_ROOT)}.", file=sys.stderr)
        return 1
    # Old versions built from THIS repo's own history (v1.5, v2.0 -- see
    # tests/tools/build_old_version.py) are compiled on first use and cached
    # under tests/crossversion/old_binaries/; later runs just reuse the cache.
    if not a.bin and a.version in build_old_version.KNOWN_COMMITS and not os.path.exists(binary):
        print(f"[gen_reference] {a.version} binary not cached locally; building it from "
              f"this repo's own history ...", file=sys.stderr)
        binary = build_old_version.ensure_binary(a.version)
    if not (os.path.isabs(binary) and os.path.exists(binary)) and shutil.which(binary) is None:
        print(f"ERROR: binary not found for version {a.version!r}: {binary}\n"
              f"       (old builds are machine-specific; see binaries.toml). Skipping.",
              file=sys.stderr)
        return 1
    if not os.path.exists(a.config):
        print(f"ERROR: config not found: {a.config}", file=sys.stderr)
        return 1

    with open(a.config) as fh:
        conf_text = fh.read()
    try:
        case = a.case or os.path.basename(io.read_exp_prefix(conf_text))
    except ValueError:
        case = a.case or os.path.splitext(os.path.basename(a.config))[0]

    import tempfile
    work_out = tempfile.mkdtemp(prefix=f"mcdc_ref_{a.version}_")
    try:
        overrides = {
            "exp_prefix": os.path.join(work_out, case),
            "seed": a.seed,
            "num_process": a.num_process,
            "write_txt": 1,        # we freeze the .txt DWI
            "write_traj_file": 0,  # NEVER freeze trajectories
        }
        if a.gpu is not None:
            overrides["gpu"] = int(a.gpu)
        materialized = io.write_conf_override(conf_text, overrides)
        conf_path = os.path.join(work_out, case + ".conf")
        with open(conf_path, "w") as fh:
            fh.write(materialized)

        print(f"[run] {a.version}: {os.path.basename(binary)} ({argstyle}) on "
              f"{os.path.relpath(a.config, paths.REPO_ROOT) if a.config.startswith(paths.REPO_ROOT) else a.config} "
              f"(seed={a.seed}, num_process={a.num_process}, gpu={a.gpu})")
        res = run_binary(binary, conf_path, a.workdir, argstyle, a.timeout)
        log = res.stdout or ""
        if a.gpu and run.is_gpu_skip(log):
            print("SKIP: the GPU path reports no CUDA device; nothing frozen.", file=sys.stderr)
            return 0
        run.raise_on_error(res)

        prefix = os.path.join(work_out, case)
        txt, bfl = prefix + "_DWI.txt", prefix + "_DWI.bfloat"

        frozen_dir = os.path.join(a.out_root, a.version)
        os.makedirs(frozen_dir, exist_ok=True)
        frozen = os.path.join(frozen_dir, case + "_DWI.txt")
        if os.path.exists(txt):
            shutil.copyfile(txt, frozen)  # verbatim -> signal_sha256 matches the binary
            src_format = "txt"
        elif os.path.exists(bfl):
            # Old builds may only write .bfloat; convert to the committed .txt form.
            np.savetxt(frozen, io.load_dwi(bfl), fmt="%.10g")
            src_format = "bfloat->txt"
        else:
            print(f"ERROR: no _DWI.txt or _DWI.bfloat produced under {prefix}. "
                  f"Does this build honour write_txt? Check the per-version config.",
                  file=sys.stderr)
            return 1

        scheme_abs = resolve_scheme(conf_text, a.workdir)
        scheme_ok = scheme_abs and os.path.exists(scheme_abs)
        n_meas = None
        if scheme_ok:
            try:
                n_meas = int(kscheme.parse_scheme(scheme_abs).n)
            except (ValueError, OSError):
                n_meas = None

        entry = provenance.make_manifest(
            case=case,
            version=a.version,
            binary_path=os.path.abspath(binary) if os.path.isabs(binary) else binary,
            argstyle=argstyle,
            git_commit=provenance.git_commit(binary),
            version_string=parse_version_string(prefix + "_simulation_info.txt"),
            seed=a.seed,
            num_process=a.num_process,
            gpu=(int(a.gpu) if a.gpu is not None else None),
            N=_int_or_none(conf_text, "N"),
            T=_int_or_none(conf_text, "T"),
            source_config=(os.path.relpath(a.config, paths.REPO_ROOT)
                           if os.path.abspath(a.config).startswith(paths.REPO_ROOT)
                           else os.path.abspath(a.config)),
            conf_sha256=provenance.sha256_text(materialized),  # exact config that ran
            scheme=(os.path.relpath(scheme_abs, paths.REPO_ROOT)
                    if scheme_ok and scheme_abs.startswith(paths.REPO_ROOT) else scheme_abs),
            scheme_sha256=(provenance.sha256_file(scheme_abs) if scheme_ok else None),
            n_measurements=n_meas,
            signal=os.path.basename(frozen),
            signal_source_format=src_format,
            signal_sha256=provenance.sha256_file(frozen),
        )

        manifest_path = os.path.join(frozen_dir, "MANIFEST.json")
        manifest = {}
        if os.path.exists(manifest_path):
            try:
                with open(manifest_path) as fh:
                    manifest = json.load(fh)
            except (OSError, ValueError):
                manifest = {}
        manifest[os.path.basename(frozen)] = entry  # append / update this case
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=True)

        rel_frozen = os.path.relpath(frozen, paths.REPO_ROOT)
        print(f"\nfroze {rel_frozen}  ({provenance.sha256_file(frozen)[:12]}..., "
              f"{src_format})")
        print(f"manifest -> {os.path.relpath(manifest_path, paths.REPO_ROOT)}")
        print("\n" + json.dumps(entry, indent=2, sort_keys=True))
        return 0
    finally:
        shutil.rmtree(work_out, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
