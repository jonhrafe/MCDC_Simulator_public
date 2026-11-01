"""Provenance capture for frozen reference data.

Every frozen signal (analytic curve, cross-version reference) ships with a
manifest recording exactly how it was produced, so a stale or mis-attributed
baseline is detectable. Builds are pinned by **git commit**, not the
``Software Version`` string (which is unreliable in this repo).
"""
import datetime
import hashlib
import os
import subprocess


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def git_commit(repo):
    """Short git commit of ``repo`` (dir or a file inside it), or None."""
    try:
        cwd = repo if os.path.isdir(repo) else os.path.dirname(repo)
        out = subprocess.run(
            ["git", "-C", cwd, "rev-parse", "--short", "HEAD"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=10)
        return out.stdout.strip() or None
    except Exception:
        return None


def now_iso():
    return datetime.datetime.now().replace(microsecond=0).isoformat()


def make_manifest(**fields):
    """Assemble a manifest dict, stamping ``date`` and ``host`` if absent."""
    m = {"date": now_iso(), "host": os.uname().nodename}
    m.update(fields)
    return m
