#!/usr/bin/env python3
"""
Build a historical MC-DC_Simulator release binary from THIS repo's own git
history (DEV-ONLY, Phase C support tool -- not part of the pip package).

The old cross-version binaries (v1.5, v2.0) do not need an external clone: the
released commits are still in this repo's own history. They were found by
bisecting ``VERSION_ID`` in ``src/constants.h`` (see ``KNOWN_COMMITS`` below and
``tests/crossversion/README.md``). ``compile.sh`` and the vendored Eigen headers
are identical at those commits and at HEAD, so the build recipe is the same
plain ``g++`` invocation used today -- no CMake, no CUDA, no second checkout of
the repo.

Mechanism: ``git worktree add --detach`` checks out the historical commit into
a scratch directory that shares this repo's object store (fast, no full clone,
never touches your current working tree or index), compiles ``src/main.cpp``
against every other ``*.cpp`` in ``src/`` there, copies the resulting binary
into a persistent local cache under ``tests/crossversion/old_binaries/<version>/``
(gitignored -- these are machine-local build artifacts, never committed), and
removes the worktree.

NOTE: this does NOT just run the checked-out ``compile.sh`` verbatim. At
fc3607c (v1.5) that script's own file list omits ``AABBFixedGrid.cpp`` /
``FixedGrid.cpp`` even though the AABB grid feature landed in that same
commit and other translation units call into it (a historical bug in that
commit's build script, not in this tool) -- linking fails with "undefined
reference to AABBFixedGrid::...". Sources are discovered dynamically instead:
every ``*.cpp`` in ``src/`` except ``main.cpp`` (added explicitly) and any file
that defines its own ``int main(`` (standalone tools bundled in the same
directory, e.g. ``dataSynth.cpp``, ``CPUPLYInt.cpp``) -- see
``_discover_sources`` below.

Usage:
    python3 tests/tools/build_old_version.py --version v1.5
    python3 tests/tools/build_old_version.py --version v2.0 --force
    python3 tests/tools/build_old_version.py --version v1.5 --keep-worktree

Importable as a library too -- ``gen_reference.py`` calls ``ensure_binary()``
directly so old binaries are built automatically the first time they are
needed and simply reused (no rebuild) on every run after that.
"""
import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # tests/
from mcdc_testkit import paths  # noqa: E402

# version label -> commit that produced that release, IN THIS REPO'S OWN
# HISTORY. Resolved by bisecting VERSION_ID (src/constants.h) with:
#   git log --all --oneline --follow -- src/constants.h
#   git show <commit>:src/constants.h | grep VERSION_ID
KNOWN_COMMITS = {
    "v1.5": "fc3607cd0326d510a51eb137799b95f3cbe0b5fc",  # VERSION_ID "1.50.000"
    "v2.0": "db833fbaa9462532c3c894fd8f609e8eba31063d",  # VERSION_ID "2.0.3.beta"
    # v2.2 has no matching VERSION_ID anywhere in this repo's history (the
    # string jumps 2.1.0 -> 2.1.1.beta -> 2.1.2.beta -> 2.1.3.beta) -- add an
    # entry here with the exact commit if/when one is identified.
}

WORKTREE_ROOT = os.path.join(paths.REPO_ROOT, "build_old")


def _git(*args, cwd=None, check=False):
    res = subprocess.run(["git", *args], cwd=cwd or paths.REPO_ROOT,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if check and res.returncode != 0:
        raise SystemExit(f"git {' '.join(args)} failed:\n{res.stdout}")
    return res


def binary_path(version):
    return os.path.join(paths.OLD_BINARIES_DIR, version, "MC-DC_Simulator")


def is_built(version):
    return os.path.exists(binary_path(version))


def _remove_worktree(worktree_dir):
    _git("worktree", "remove", "--force", worktree_dir)
    _git("worktree", "prune")


_MAIN_RE = re.compile(r'^\s*int\s+main\s*\(', re.MULTILINE)


def _discover_sources(src_dir):
    """Every *.cpp in src_dir for the MC-DC_Simulator target: all of them
    except main.cpp (added explicitly, first) and any file defining its own
    `int main(` -- standalone tools living in the same directory (dataSynth.cpp,
    CPUPLYInt.cpp, ...), auto-detected rather than a hardcoded exclude list so
    this doesn't need updating per-commit.
    """
    sources = []
    for name in sorted(os.listdir(src_dir)):
        if not name.endswith(".cpp") or name == "main.cpp":
            continue
        with open(os.path.join(src_dir, name), errors="ignore") as fh:
            text = fh.read()
        if _MAIN_RE.search(text):
            continue
        sources.append(name)
    return sources


def _compile(worktree_dir):
    """g++ main.cpp + discovered sources -> <worktree_dir>/MC-DC_Simulator.

    Same flags as the historical compile.sh (-O3 -std=c++17 -march=native
    -flto -pthread); Eigen is vendored under src/Eigen at every commit this
    tool targets, so -I. (from within src/) is the only include path needed.
    """
    src_dir = os.path.join(worktree_dir, "src")
    sources = _discover_sources(src_dir)
    cmd = ["g++", "-O3", "-std=c++17", "-march=native", "-flto",
           "-lpthread", "-pthread", "-w", "-I.", "main.cpp", *sources,
           "-o", "../MC-DC_Simulator"]
    return subprocess.run(cmd, cwd=src_dir, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)


def ensure_binary(version, commit=None, force=False, keep_worktree=False):
    """Return a locally-built old-version binary, building it only if missing.

    Cached at tests/crossversion/old_binaries/<version>/MC-DC_Simulator. Set
    force=True to rebuild even if a cached binary already exists.
    """
    out_bin = binary_path(version)
    if os.path.exists(out_bin) and not force:
        return out_bin

    commit = commit or KNOWN_COMMITS.get(version)
    if not commit:
        raise SystemExit(
            f"no known release commit for version {version!r}; pass --commit "
            f"explicitly (known: {', '.join(sorted(KNOWN_COMMITS))})")

    resolved = _git("rev-parse", "--verify", f"{commit}^{{commit}}")
    if resolved.returncode != 0:
        raise SystemExit(f"commit {commit!r} not found in this repo:\n{resolved.stdout}")
    commit_sha = resolved.stdout.strip()

    worktree_dir = os.path.join(WORKTREE_ROOT, version)
    if os.path.exists(worktree_dir):
        _remove_worktree(worktree_dir)
    os.makedirs(WORKTREE_ROOT, exist_ok=True)

    print(f"[build_old_version] {version}: checking out {commit_sha[:10]} into "
          f"{os.path.relpath(worktree_dir, paths.REPO_ROOT)} ...", file=sys.stderr)
    _git("worktree", "add", "--detach", worktree_dir, commit_sha, check=True)

    try:
        print(f"[build_old_version] {version}: compiling ...", file=sys.stderr)
        build = _compile(worktree_dir)
        if build.returncode != 0:
            sys.stderr.write(build.stdout)
            raise SystemExit(f"build failed for {version} @ {commit_sha[:10]}")

        built = os.path.join(worktree_dir, "MC-DC_Simulator")
        if not os.path.exists(built):
            raise SystemExit(
                f"compile.sh reported success but {built} was not produced")

        os.makedirs(os.path.dirname(out_bin), exist_ok=True)
        shutil.copy2(built, out_bin)
        os.chmod(out_bin, 0o755)

        gxx = subprocess.run(["g++", "--version"], stdout=subprocess.PIPE, text=True)
        stamp = {
            "version": version,
            "commit": commit_sha,
            "built_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "compiler": (gxx.stdout.splitlines() or [""])[0],
        }
        with open(os.path.join(os.path.dirname(out_bin), "BUILD_INFO.json"), "w") as fh:
            json.dump(stamp, fh, indent=2)

        print(f"[build_old_version] {version}: built "
              f"{os.path.relpath(out_bin, paths.REPO_ROOT)} from {commit_sha[:10]}",
              file=sys.stderr)
        return out_bin
    finally:
        if keep_worktree:
            print(f"[build_old_version] kept worktree at "
                  f"{os.path.relpath(worktree_dir, paths.REPO_ROOT)} (--keep-worktree)",
                  file=sys.stderr)
        else:
            _remove_worktree(worktree_dir)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", required=True,
                     help=f"which historical release to build (known: "
                          f"{', '.join(sorted(KNOWN_COMMITS))}; others need --commit)")
    ap.add_argument("--commit", default=None,
                     help="override the commit (default: the known release commit)")
    ap.add_argument("--force", action="store_true",
                     help="rebuild even if a cached binary already exists")
    ap.add_argument("--keep-worktree", action="store_true",
                     help="do not remove the temporary git worktree after building "
                          "(handy for iterating on a build failure)")
    a = ap.parse_args()
    path = ensure_binary(a.version, commit=a.commit, force=a.force,
                          keep_worktree=a.keep_worktree)
    print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
