#!/usr/bin/env python3
"""
Cross-version DWI comparator: current build vs the previous releases and the analytic oracle.

Given a set of LABELLED ``_DWI`` files (the current build + frozen previous-
version references under tests/accuracy/reference_xv/) plus a scheme and an
OPTIONAL analytic ``.npy`` curve, it computes per-shell mean S/S0, per-shell
ADC and NRMSE-vs-reference tables and writes a **markdown report**. The report
is ALWAYS produced (numpy only). PNG plots are produced ONLY if matplotlib
imports; otherwise they are skipped with a note (never a crash).

    python3 tests/crossversion/compare_versions.py \
        --scheme tests/accuracy/PGSE_sample_scheme.scheme \
        --signal current=tests/accuracy/output/free_diffusion_DWI.txt \
        --signal v1.5=tests/accuracy/reference_xv/v1.5/free_diffusion_DWI.txt \
        --analytic cylinder_r5um_d0p6 \
        --out report.md --plot-dir plots

NORMALISATION -- the simulator writes the DWI as a SUM over walkers; S/S0 is
signal / mean(signal[b0]), so the walker count N CANCELS and signals from runs
with different N are directly comparable (no --n-walkers needed). The analytic
``.npy`` is already full-length S/S0 in scheme order (see gen_analytic.py).

DETERMINISM -- cross-version and CPU<->GPU signals are NEVER bit-exact
(different RNGs); this compares per-shell MEANS within a statistical band, with
the analytic curve as the primary oracle. Missing inputs (a reference not frozen
on this host, matplotlib absent) are reported and skipped, not fatal.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # tests/
from mcdc_testkit import io, metrics, paths, provenance  # noqa: E402
from mcdc_testkit import scheme as kscheme  # noqa: E402

import numpy as np  # noqa: E402

# matplotlib is optional: the markdown report never depends on it.
try:
    import matplotlib
    matplotlib.use("Agg")  # headless: only writes PNGs
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:  # noqa: BLE001 - any import problem -> just skip plots
    HAVE_MPL = False


def parse_labeled(items):
    """``["label=path", ...]`` -> ``[(label, path), ...]`` (preserving order)."""
    out = []
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"--signal must be LABEL=PATH, got {it!r}")
        label, path = it.split("=", 1)
        out.append((label.strip(), os.path.expanduser(path.strip())))
    return out


def resolve_analytic(name_or_path):
    """Resolve an analytic curve given a path or a bare name in analytic/."""
    if not name_or_path:
        return None
    cand = name_or_path if name_or_path.endswith(".npy") else name_or_path + ".npy"
    for p in (cand, os.path.join(paths.ANALYTIC_DIR, os.path.basename(cand))):
        if os.path.exists(p):
            return p
    return None


def load_curves(labeled, analytic_path, sc):
    """Return ``{label: full-length S/S0 (scheme order)}`` + notes for skips.

    Simulator signals are normalised by their own b0 mean; the analytic curve is
    already S/S0. Length mismatches / missing files are recorded and skipped.
    """
    curves, notes = {}, []
    for label, path in labeled:
        if not os.path.exists(path):
            notes.append(f"- `{label}`: **skipped** — file not found (`{path}`).")
            continue
        arr = io.load_dwi(path)
        if arr.size != sc.n:
            notes.append(f"- `{label}`: **skipped** — {arr.size} values != {sc.n} scheme rows.")
            continue
        try:
            curves[label] = metrics.normalize_by_b0(arr, sc.b0_mask)
        except ValueError as exc:
            notes.append(f"- `{label}`: **skipped** — {exc}.")
    if analytic_path:
        ana = np.load(analytic_path).astype(np.float64).ravel()
        if ana.size == sc.n:
            curves["analytic"] = ana
        else:
            notes.append(f"- `analytic`: **skipped** — {ana.size} values != {sc.n} scheme rows.")
    return curves, notes


def fmt_table(headers, rows):
    """Minimal GitHub-markdown table."""
    line = "| " + " | ".join(headers) + " |"
    sep = "| " + " | ".join("---" for _ in headers) + " |"
    body = ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join([line, sep, *body])


def _cell(x, nan="—"):
    return nan if (x is None or not np.isfinite(x)) else f"{x:.4f}"


def build_report(curves, sc, ref_label, title, plot_notes, decimals):
    """Assemble the markdown report string (numpy-only, always available)."""
    labels = list(curves)
    shells = metrics.shells(sc.bvals, decimals=decimals)          # s/mm^2, sorted
    per_mean = {lab: metrics.per_shell_mean(curves[lab], sc.bvals, decimals=decimals)
                for lab in labels}
    nonb0 = ~sc.b0_mask

    md = []
    md.append(f"# {title}\n")
    md.append(f"_Generated {provenance.now_iso()} · host `{os.uname().nodename}` · "
              f"repo `{provenance.git_commit(paths.REPO_ROOT) or 'n/a'}`_\n")
    md.append(f"- **Measurements:** {sc.n} ({int(np.count_nonzero(sc.b0_mask))} b0)")
    md.append(f"- **Shells [s/mm²]:** {', '.join(f'{b:g}' for b in shells)}")
    md.append(f"- **NRMSE reference:** `{ref_label}`")
    md.append(f"- **Plots:** {'written (see below)' if HAVE_MPL else '**skipped** — matplotlib not importable'}\n")

    # 1. per-shell mean S/S0
    md.append("## Per-shell mean S/S0\n")
    headers = ["b [s/mm²]", *[f"`{l}`" for l in labels]]
    rows = [[f"{b:g}", *[_cell(per_mean[l].get(b)) for l in labels]] for b in shells]
    md.append(fmt_table(headers, rows) + "\n")

    # 2. per-shell ADC (skip the b=0 shell; ADC = -ln(S/S0)/b)
    md.append("## Per-shell ADC  [×10⁻⁹ m²/s]\n")
    nz = [b for b in shells if b > 0]
    rows = []
    for b in nz:
        b_si = b * 1e6
        cells = []
        for l in labels:
            m = per_mean[l].get(b)
            adc = metrics.adc(np.array([m]), np.array([b_si]))[0] if m else np.nan
            cells.append(_cell(adc * 1e9))
        rows.append([f"{b:g}", *cells])
    md.append(fmt_table(["b [s/mm²]", *[f"`{l}`" for l in labels]], rows) + "\n")

    # 3. NRMSE + max|Δ| vs the reference (over non-b0 measurements)
    md.append(f"## Agreement vs `{ref_label}`  (non-b0 measurements)\n")
    if ref_label in curves:
        ref = curves[ref_label][nonb0]
        rows = []
        for l in labels:
            if l == ref_label:
                continue
            run = curves[l][nonb0]
            rows.append([f"`{l}`", f"{metrics.nrmse(run, ref):.4f}",
                         f"{metrics.max_abs_dev(run, ref):.4f}"])
        md.append(fmt_table(["version", "NRMSE (range-norm)", "max abs Δ(S/S0)"], rows) + "\n")
    else:
        md.append(f"_Reference `{ref_label}` not among the loaded curves — NRMSE omitted._\n")

    if plot_notes:
        md.append("## Notes\n")
        md.extend(plot_notes)
        md.append("")

    md.append("---")
    md.append("_Cross-version and CPU↔GPU signals are never bit-exact (different RNGs); "
              "agreement is statistical (per-shell means within an MC-noise band). The "
              "analytic curve, when present, is the primary oracle; previous-version "
              "signals are a secondary historical overlay._")
    return "\n".join(md)


def make_plots(curves, sc, plot_dir, decimals):
    """Write mean-S/S0-vs-b (lin+log) and ADC-vs-b PNGs; return relative links."""
    os.makedirs(plot_dir, exist_ok=True)
    shells = metrics.shells(sc.bvals, decimals=decimals)
    per_mean = {lab: metrics.per_shell_mean(curves[lab], sc.bvals, decimals=decimals)
                for lab in curves}

    def series(lab):
        return np.array([per_mean[lab].get(b, np.nan) for b in shells])

    saved = []

    def _emph(lab):  # analytical drawn as a bold black reference line
        return dict(color="black", lw=2.4, marker="D", ms=5, zorder=5) if lab == "analytic" \
            else dict(lw=1.6, marker="o", ms=5, alpha=0.85, zorder=3)

    for yscale in ("linear", "log"):
        fig = plt.figure(figsize=(8, 5))
        for lab in curves:
            y = series(lab).copy()
            if yscale == "log":
                y[y <= 0] = np.nan
            plt.plot(shells, y, label=lab, **_emph(lab))
        plt.yscale(yscale)
        plt.xlabel("b-value [s/mm²]")
        plt.ylabel(f"mean S/S0{' (log)' if yscale == 'log' else ''}")
        plt.title(f"Mean signal per shell ({yscale})")
        plt.grid(True, which="both", ls=":", alpha=0.4)
        plt.legend(fontsize=8)
        name = f"mean_signal_per_shell_{yscale}.png"
        fig.savefig(os.path.join(plot_dir, name), dpi=120, bbox_inches="tight")
        plt.close(fig)
        saved.append(name)

    nz = np.array([b for b in shells if b > 0])
    if nz.size:
        fig = plt.figure(figsize=(8, 5))
        for lab in curves:
            m = np.array([per_mean[lab].get(b, np.nan) for b in nz])
            adc = metrics.adc(m, nz * 1e6) * 1e9
            plt.plot(nz, adc, label=lab, **_emph(lab))
        plt.xlabel("b-value [s/mm²]")
        plt.ylabel("ADC [×10⁻⁹ m²/s]")
        plt.title("ADC per shell")
        plt.grid(True, ls=":", alpha=0.4)
        plt.legend(fontsize=8)
        fig.savefig(os.path.join(plot_dir, "adc_per_shell.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)
        saved.append("adc_per_shell.png")
    return saved


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scheme", required=True, help="STEJSKALTANNER scheme file")
    ap.add_argument("--signal", action="append", metavar="LABEL=PATH", required=True,
                    help="a labelled _DWI file; repeat for each version")
    ap.add_argument("--analytic", default=None,
                    help="analytic S/S0 .npy (path or a bare name in accuracy/analytic/)")
    ap.add_argument("--ref", default=None,
                    help="NRMSE reference label (default: 'analytic' if present, else the first signal)")
    ap.add_argument("--out", default="crossversion_report.md", help="markdown report path")
    ap.add_argument("--plot-dir", default=None,
                    help="write PNGs here (needs matplotlib; omit to skip plots)")
    ap.add_argument("--title", default="Cross-version DWI comparison")
    ap.add_argument("--decimals", type=int, default=0, help="b-shell rounding (s/mm^2)")
    a = ap.parse_args()

    if not os.path.exists(a.scheme):
        print(f"ERROR: scheme not found: {a.scheme}", file=sys.stderr)
        return 1
    sc = kscheme.parse_scheme(a.scheme)

    analytic_req = a.analytic
    analytic_path = resolve_analytic(analytic_req)
    curves, notes = load_curves(parse_labeled(a.signal), analytic_path, sc)
    if analytic_req and not analytic_path:
        notes.append(f"- `analytic`: **skipped** — curve `{analytic_req}` not found "
                     f"(generate with tests/tools/gen_analytic.py).")
    if not curves:
        print("ERROR: no usable signals loaded (all missing or wrong length):", file=sys.stderr)
        for n in notes:
            print("  " + n, file=sys.stderr)
        return 1

    ref_label = a.ref or ("analytic" if "analytic" in curves else next(iter(curves)))

    plot_notes = list(notes)
    if a.plot_dir:
        if HAVE_MPL:
            saved = make_plots(curves, sc, a.plot_dir, a.decimals)
            rel = os.path.relpath(a.plot_dir, os.path.dirname(os.path.abspath(a.out)) or ".")
            plot_notes.append("Plots written: " + ", ".join(
                f"`{os.path.join(rel, s)}`" for s in saved))
        else:
            plot_notes.append("- Plots **skipped**: matplotlib is not importable in this "
                              "environment (`pip install .[viz]`). Tables above are complete.")

    report = build_report(curves, sc, ref_label, a.title, plot_notes, a.decimals)
    with open(a.out, "w") as fh:
        fh.write(report + "\n")

    print(f"wrote {a.out}  ({len(curves)} curve(s): {', '.join(curves)})")
    if a.plot_dir and not HAVE_MPL:
        print("note: matplotlib not available -> plots skipped (report still written).")
    for n in notes:
        print("  " + n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
