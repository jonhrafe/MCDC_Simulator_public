# Frozen cross-version reference signals (previous-version "golden" references)

DWI signals produced by running the **older stable simulator binaries**
(v1.5, v2.0) on the canonical validation configs with a **fixed seed
and fixed `num_process`**, then frozen and committed here. This is the Tier-3
cross-version baseline: it lets the current build be checked against what the
previous stable versions produced.

Populated by [`tests/tools/gen_reference.py`](../../tools/gen_reference.py)
(Phase C). Layout: `reference_xv/<version>/<case>_DWI.txt` plus a
`MANIFEST.json` recording, per file: binary path, **git commit** (version
strings in this repo are unreliable), seed, `num_process`, N, T, scheme/conf/
signal sha256, date, host.

Comparison is **statistical, never bit-exact** — different versions use different
RNGs, so results only agree within a loose per-shell band (MC noise in
quadrature). The **analytic curve is the primary oracle**; these previous-version
signals are a secondary overlay showing historical agreement. See
[`tests/crossversion/README.md`](../../crossversion/README.md).
