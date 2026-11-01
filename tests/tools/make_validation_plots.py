#!/usr/bin/env python3
"""
Generate CPU-vs-GPU-vs-analytic validation plots, organised by the suite's TEST TIERS.

Most validation tests print a PASS/FAIL and a worst-case number; this tool draws the
picture behind that number. For each experiment it runs the simulator on the CPU and
(when CUDA is present) the GPU at a modest N, overlays the two Monte-Carlo signals on
the analytic oracle, and adds a residual panel:

    * MC(CPU) - analytic  and  MC(GPU) - analytic   -> physics error vs the oracle
    * MC(CPU) - MC(GPU)                             -> CPU/GPU hardware parity

Output mirrors the CTest tier taxonomy (see tests/README.md), one sub-folder per
experiment under its tier, in `--out` (default `tests/output/`, git-ignored):

    tests/output/
      tier1_smoke/       free_diffusion
      tier2_fast/        (scaffold: golden is bit-exact, not an overlay)
      tier3_analytic/    restricted_sphere  restricted_cylinder  cross_representation
                         packing_restricted  packing_exchange  mesh_cell  time_dependence
      tier4_crossversion/ free_diffusion  (current vs v1.5 vs v2.0; old binaries build locally)
      tier5_heavy/       (scaffold: large meshes, no analytic oracle)

`gpu` is a cross-cutting MODIFIER, not a tier: every experiment runs on the CPU AND the
GPU, so the CPU/GPU comparison lives in each figure's residual panel. Each experiment is
independent — a failure (or an absent GPU) is logged and the rest still run. A markdown
table is always written; PNGs need matplotlib (the `[test]` extra installs it).

  python3 tests/tools/make_validation_plots.py --bin ./MC-DC_Simulator --out tests/output
  python3 tests/tools/make_validation_plots.py --only packing_restricted --n 20000
"""
import argparse
import os
import re
import shutil
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tests"))
sys.path.insert(0, os.path.join(ROOT, "tests", "validation"))
sys.path.insert(0, HERE)  # tests/tools -- gen_reference / build_old_version (crossversion)

from mcdc_testkit import io, metrics, oracles, paths, run as krun, scheme as kscheme  # noqa: E402
import test_free_exp as free_mod            # noqa: E402
import test_sphere_neuman as sph_mod        # noqa: E402
import test_cylinder_vangelderen as cyl_mod  # noqa: E402
import test_cross_representation as xrep_mod  # noqa: E402
import test_packing_restricted as pack_mod   # noqa: E402
import test_packing_exchange as exch_mod      # noqa: E402
import test_mesh_cell as mesh_mod             # noqa: E402
import test_time_dependence as td_mod         # noqa: E402
import numpy as np                            # noqa: E402


# --------------------------------------------------------------------------- #
#  matplotlib helpers                                                         #
# --------------------------------------------------------------------------- #
def _plt():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except Exception:
        return None


def _shell_means(ss0, sc):
    """(bs, means-dict) of per-shell mean S/S0 over the non-b0 measurements."""
    nb = ~sc.b0_mask
    means = metrics.per_shell_mean(ss0[nb], sc.bvals[nb])
    return sorted(means), means


def _write_table(path, title, bs, cols):
    head = "| b (s/mm²) | " + " | ".join(lab for lab, _ in cols) + " |"
    rule = "|" + "---|" * (len(cols) + 1)
    rows = [f"# {title}\n", head, rule]
    for b in bs:
        cells = " | ".join(f"{d.get(b, float('nan')):.4f}" for _, d in cols)
        rows.append(f"| {b:.0f} | {cells} |")
    with open(path, "w") as fh:
        fh.write("\n".join(rows) + "\n")


def overlay_ss0(out_dir, fname, title, sc, cpu_ss0, gpu_ss0, oracle_ss0,
                oracle_label="analytic"):
    """Shared plot: per-shell mean S/S0 vs b for CPU/GPU markers on the oracle line,
    plus a residual panel (CPU-oracle, GPU-oracle, CPU-GPU). Writes .md + .png.
    Returns a one-line summary string."""
    os.makedirs(out_dir, exist_ok=True)
    bs, m_cpu = _shell_means(cpu_ss0, sc)
    have_gpu = gpu_ss0 is not None
    m_gpu = _shell_means(gpu_ss0, sc)[1] if have_gpu else {}
    m_ora = None
    if oracle_ss0 is not None:
        m_ora = _shell_means(oracle_ss0, sc)[1]

    cols = [("MC CPU", m_cpu)]
    if have_gpu:
        cols.append(("MC GPU", m_gpu))
    if m_ora is not None:
        cols.append((oracle_label, m_ora))
    if have_gpu:
        cols.append(("CPU−GPU", {b: m_cpu[b] - m_gpu[b] for b in bs}))
    _write_table(os.path.join(out_dir, fname + ".md"), title, bs, cols)

    worst_hw = max((abs(m_cpu[b] - m_gpu[b]) for b in bs), default=float("nan")) if have_gpu else float("nan")
    worst_or = max((abs(m_cpu[b] - m_ora[b]) for b in bs), default=float("nan")) if m_ora is not None else float("nan")
    summary = f"worst |CPU−GPU|={worst_hw:.4f}" + (f"  worst |CPU−{oracle_label}|={worst_or:.4f}" if m_ora is not None else "")

    plt = _plt()
    if plt is None:
        return summary + "  (no matplotlib -> md only)"
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.2, 7.0),
                                  gridspec_kw=dict(height_ratios=[3, 1.2]), sharex=True)
    if m_ora is not None:
        ax.plot(bs, [m_ora[b] for b in bs], "-k", lw=2, label=oracle_label)
    ax.plot(bs, [m_cpu[b] for b in bs], "o", ms=7, color="tab:blue", label="MC CPU")
    if have_gpu:
        ax.plot(bs, [m_gpu[b] for b in bs], "x", ms=9, mew=2, color="tab:red", label="MC GPU")
    ax.set_ylabel("S/S0 (per-shell mean)")
    ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend()

    axr.axhline(0, color="k", lw=0.6)
    if m_ora is not None:
        axr.plot(bs, [m_cpu[b] - m_ora[b] for b in bs], "o-", color="tab:blue", label=f"CPU − {oracle_label}")
        if have_gpu:
            axr.plot(bs, [m_gpu[b] - m_ora[b] for b in bs], "x-", color="tab:red", label=f"GPU − {oracle_label}")
    if have_gpu:
        axr.plot(bs, [m_cpu[b] - m_gpu[b] for b in bs], "s--", color="tab:green", label="CPU − GPU")
    axr.set_xlabel("b (s/mm²)")
    axr.set_ylabel("Δ S/S0")
    axr.grid(True, ls=":", alpha=0.4)
    axr.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, fname + ".png"), dpi=120)
    plt.close(fig)
    return summary


# --------------------------------------------------------------------------- #
#  Per-experiment adapters (each returns a status string)                     #
# --------------------------------------------------------------------------- #
def exp_free(binary, n, out_dir):
    r_cpu = xrep_mod.run_conf(binary, "free_diffusion.conf", n=n, gpu=0)
    r_gpu = xrep_mod.run_conf(binary, "free_diffusion.conf", n=n, gpu=1)
    D = float(io.read_conf_value(open(paths.accuracy("free_diffusion.conf")).read(), "diffusivity"))
    ora = oracles.free_exp(r_cpu.sc.bvals, D)
    gpu = None if (r_gpu is None or r_gpu.skipped) else r_gpu.ss0
    return overlay_ss0(out_dir, "free_diffusion", f"Free diffusion vs exp(−bD)  (D={D*1e9:.2f}e-9, N={n})",
                       r_cpu.sc, r_cpu.ss0, gpu, ora, oracle_label="exp(−bD)")


def exp_sphere(binary, n, out_dir):
    r_cpu = sph_mod.run_sphere(binary, n=n, gpu=0)
    r_gpu = sph_mod.run_sphere(binary, n=n, gpu=1)
    ora = oracles.sphere_neuman(r_cpu.sc, sph_mod.RADIUS, r_cpu.D)
    gpu = None if r_gpu.skipped else r_gpu.ss0
    return overlay_ss0(out_dir, "restricted_sphere",
                       f"Restricted sphere (r={sph_mod.RADIUS*1e6:.0f}µm) vs Neuman GPA  (N={n})",
                       r_cpu.sc, r_cpu.ss0, gpu, ora, oracle_label="Neuman GPA")


def exp_cylinder(binary, n, out_dir):
    r_cpu = cyl_mod.run_cylinder(binary, n=n, gpu=0)
    r_gpu = cyl_mod.run_cylinder(binary, n=n, gpu=1)
    ora = oracles.cylinder_vangelderen(r_cpu.sc, cyl_mod.RADIUS, r_cpu.D)
    gpu = None if r_gpu.skipped else r_gpu.ss0
    return overlay_ss0(out_dir, "restricted_cylinder",
                       f"Restricted cylinder (r={cyl_mod.RADIUS*1e6:.0f}µm) vs van Gelderen GPA  (N={n})",
                       r_cpu.sc, r_cpu.ss0, gpu, ora, oracle_label="van Gelderen GPA")


def exp_cross_rep(binary, n, out_dir):
    r_list = xrep_mod.run_conf(binary, xrep_mod.LIST_CONF, n=n, gpu=0)
    r_ply = xrep_mod.run_conf(binary, xrep_mod.PLY_CONF, n=n, gpu=0)
    ora = xrep_mod._vangelderen_curve(r_list.sc)
    # here "CPU" = native list, "GPU"-slot reused to overlay the PLY representation
    return overlay_ss0(out_dir, "cross_representation",
                       f"Cross-representation: native list vs PLY vs analytic  (N={n})",
                       r_list.sc, r_list.ss0, r_ply.ss0, ora, oracle_label="van Gelderen")


def exp_packing(binary, n, out_dir):
    lines = []
    for name in ("sphere", "cylinder"):
        cpu = pack_mod.run_case(binary, name, n=n, gpu=0)
        gpu = pack_mod.run_case(binary, name, n=n, gpu=1)
        gpu_ss0 = None if gpu.skipped else gpu.mc
        s = overlay_ss0(out_dir, f"packing_{name}",
                        f"{name} packing (restricted) vs volume-weighted GPA  ({cpu.n_obs} obs, N={n})",
                        cpu.sc, cpu.mc, gpu_ss0, cpu.ana, oracle_label="GPA (vw)")
        lines.append(f"{name}: {s}")
    return " | ".join(lines)


def exp_exchange(binary, n, out_dir):
    """Sphere packing, total S/S0 vs b at increasing permeability (CPU), + CPU/GPU at high κ."""
    kappas = [0.0, 1e-5, 1e-4]
    variants = {k: exch_mod.run_variant(binary, "sphere", k, "intra", n=n, gpu=0) for k in kappas}
    sc = variants[0.0].sc
    os.makedirs(out_dir, exist_ok=True)
    plt = _plt()
    bs = _shell_means(variants[0.0].total, sc)[0]
    cols = [(f"κ={k:g}", _shell_means(variants[k].total, sc)[1]) for k in kappas]
    cols.append(("f_extra(κ)", {kappas[i]: variants[k].fe for i, k in enumerate(kappas)}))  # informational
    _write_table(os.path.join(out_dir, "packing_exchange.md"),
                 f"Sphere packing exchange — total S/S0 vs b (N={n})", bs,
                 [(lab, d) for lab, d in cols if all(isinstance(x, float) for x in d.values())][:len(kappas)])
    # CPU vs GPU at the highest kappa
    gpu_hi = exch_mod.run_variant(binary, "sphere", kappas[-1], "intra", n=n, gpu=1)
    fe = ", ".join(f"κ={k:g}:{variants[k].fe:.3f}" for k in kappas)
    if plt is None:
        return f"f_extra grows: {fe}  (no matplotlib)"
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.2, 7.0),
                                  gridspec_kw=dict(height_ratios=[3, 1.2]), sharex=True)
    colors = ["tab:blue", "tab:orange", "tab:red"]
    for k, c in zip(kappas, colors):
        m = _shell_means(variants[k].total, sc)[1]
        ax.plot(bs, [m[b] for b in bs], "o-", color=c, label=f"CPU κ={k:g} (f_e={variants[k].fe:.2f})")
    if gpu_hi is not None:
        mg = _shell_means(gpu_hi.total, sc)[1]
        ax.plot(bs, [mg[b] for b in bs], "x", ms=9, mew=2, color="k", label=f"GPU κ={kappas[-1]:g}")
    ax.set_ylabel("total S/S0 (per-shell mean)")
    ax.set_title(f"Sphere packing exchange: permeability raises the signal toward free  (N={n})")
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(fontsize=8)
    m0 = _shell_means(variants[0.0].total, sc)[1]
    axr.axhline(0, color="k", lw=0.6)
    for k, c in zip(kappas[1:], colors[1:]):
        mk = _shell_means(variants[k].total, sc)[1]
        axr.plot(bs, [mk[b] - m0[b] for b in bs], "o-", color=c, label=f"κ={k:g} − κ=0")
    if gpu_hi is not None:
        mg = _shell_means(gpu_hi.total, sc)[1]
        mc = _shell_means(variants[kappas[-1]].total, sc)[1]
        axr.plot(bs, [mc[b] - mg[b] for b in bs], "s--", color="tab:green", label=f"CPU−GPU (κ={kappas[-1]:g})")
    axr.set_xlabel("b (s/mm²)")
    axr.set_ylabel("Δ S/S0")
    axr.grid(True, ls=":", alpha=0.4)
    axr.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "packing_exchange.png"), dpi=120)
    plt.close(fig)
    return f"f_extra grows with κ: {fe}"


def exp_mesh(binary, n, out_dir):
    cpu, sc, _ = mesh_mod.run_mesh(binary, n=n, gpu=0)
    gpu, _, _ = mesh_mod.run_mesh(binary, n=n, gpu=1)
    return overlay_ss0(out_dir, "mesh_cpu_gpu",
                       f"Real cell mesh — CPU vs GPU (no closed form; parity is the oracle)  (N={n})",
                       sc, cpu, gpu, None)


def exp_time_dependence(binary, n, out_dir):
    sc = kscheme.parse_scheme(os.path.join(paths.REPO_ROOT, td_mod.SCHEME))
    imp_cpu = td_mod._run(binary, "native", 0.0, n=n, gpu=0)
    imp_gpu = td_mod._run(binary, "native", 0.0, n=n, gpu=1)
    gpa = oracles.sphere_neuman(sc, td_mod.R, td_mod.D)
    os.makedirs(out_dir, exist_ok=True)
    plt = _plt()

    # (a) restriction: S/S0 vs b, one curve per Delta, CPU markers + GPA line
    if plt is not None:
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        cmap = plt.get_cmap("viridis")
        for i, dl in enumerate(td_mod.DELTAS):
            m = td_mod._by(sc, dl)
            order = np.argsort(sc.bvals[m])
            x = (sc.bvals[m] / 1e6)[order]
            col = cmap(i / max(len(td_mod.DELTAS) - 1, 1))
            ax.plot(x, imp_cpu.ss0[m][order], "o-", color=col, label=f"Δ={dl}ms CPU")
            ax.plot(x, gpa[m][order], "--", color=col, alpha=0.7)
            if imp_gpu is not None and not imp_gpu.skipped:
                ax.plot(x, imp_gpu.ss0[m][order], "x", color=col, ms=8, mew=1.5)
        ax.set_xlabel("b (s/mm²)")
        ax.set_ylabel("S/S0")
        ax.set_title(f"Time-dependence (κ=0): restriction grows with Δ; markers=MC, dashed=Neuman GPA  (N={n})")
        ax.grid(True, ls=":", alpha=0.4)
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "timedep_restriction.png"), dpi=120)
        plt.close(fig)

    # (b) exchange: f_extra vs kappa (CPU and GPU)
    fes_cpu, fes_gpu = [], []
    for k in td_mod.KAPPA:
        fes_cpu.append(td_mod._run(binary, "native", k, n=n, gpu=0).f_extra)
        g = td_mod._run(binary, "native", k, n=n, gpu=1)
        fes_gpu.append(None if (g is None or g.skipped) else g.f_extra)
    if plt is not None:
        fig, ax = plt.subplots(figsize=(7.0, 4.6))
        ax.plot(range(len(td_mod.KAPPA)), fes_cpu, "o-", color="tab:blue", label="CPU")
        if any(v is not None for v in fes_gpu):
            ax.plot(range(len(td_mod.KAPPA)), [v if v is not None else np.nan for v in fes_gpu],
                    "x--", color="tab:red", label="GPU")
        ax.set_xticks(range(len(td_mod.KAPPA)))
        ax.set_xticklabels([f"{k:g}" for k in td_mod.KAPPA])
        ax.set_xlabel("permeability κ (m/s)")
        ax.set_ylabel("extra-compartment fraction f_extra")
        ax.set_title("Exchange: f_extra rises with permeability")
        ax.grid(True, ls=":", alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "timedep_exchange.png"), dpi=120)
        plt.close(fig)
    fe = ", ".join(f"{k:g}:{f:.3f}" for k, f in zip(td_mod.KAPPA, fes_cpu))
    return f"f_extra(κ)=[{fe}]"


# --------------------------------------------------------------------------- #
#  Registry + driver                                                          #
# --------------------------------------------------------------------------- #
def _run_version_ss0(version, conf_relpath, n, explicit_binary=None, seed=12345, timeout=1800):
    """Run one build (current / v1.5 / v2.0) on its per-version conf; return (ss0, sc).

    Reuses the crossversion pipeline (tests/tools/gen_reference.py + build_old_version.py):
    resolve the binary from binaries.toml -- old versions COMPILE from this repo's own git
    history on first use and cache under old_binaries/ -- run it the way that version expects
    (--conf vs positional), then normalise its _DWI.txt by the b0 shell.
    """
    import gen_reference as gr          # tests/tools sibling
    import build_old_version            # tests/tools sibling

    cfg = gr.load_versions()
    binary, argstyle = gr.resolve_binary(version, explicit_binary, cfg)
    if not explicit_binary and version in build_old_version.KNOWN_COMMITS and (
            not binary or not os.path.exists(binary)):
        binary = build_old_version.ensure_binary(version)   # builds + caches on first use
    if not binary:
        raise RuntimeError(f"no binary registered for {version!r}")

    conf_abs = os.path.join(paths.REPO_ROOT, conf_relpath)
    conf_text = open(conf_abs).read()
    sc = kscheme.parse_scheme(os.path.join(paths.REPO_ROOT, io.read_conf_value(conf_text, "scheme_file")))
    N = int(n) if n else int(io.read_conf_value(conf_text, "N"))

    work = tempfile.mkdtemp(prefix=f"mcdc_xv_{version.replace('.', 'p')}_")
    try:
        prefix = os.path.join(work, "xv")
        materialized = io.write_conf_override(conf_text, {
            "exp_prefix": prefix, "N": N, "seed": seed, "num_process": 1,
            "write_txt": 1, "write_bin": 0, "write_traj_file": 0,
        })
        cpath = os.path.join(work, "c.conf")
        with open(cpath, "w") as fh:
            fh.write(materialized)
        res = gr.run_binary(binary, cpath, paths.REPO_ROOT, argstyle, timeout)
        krun.raise_on_error(res)
        dwi = io.load_dwi(prefix + "_DWI.txt", n_walkers=N)
        return metrics.normalize_by_b0(dwi, sc.b0_mask), sc
    finally:
        shutil.rmtree(work, ignore_errors=True)


def exp_crossversion(binary, n, out_dir):
    """free_diffusion: the current build vs v1.5 vs v2.0, all overlaid on exp(-bD)."""
    confs = [
        ("current", "tests/accuracy/free_diffusion.conf", binary),
        ("v1.5", "tests/crossversion/confs/v1.5/free_diffusion.conf", None),
        ("v2.0", "tests/crossversion/confs/v2.0/free_diffusion.conf", None),
    ]
    markers = {"current": ("o", "tab:blue"), "v1.5": ("s", "tab:orange"), "v2.0": ("^", "tab:green")}
    got, sc = {}, None
    for ver, crel, exbin in confs:
        try:
            ss0, s = _run_version_ss0(ver, crel, n, explicit_binary=exbin)
            got[ver] = ss0
            sc = sc or s
        except Exception as exc:  # noqa: BLE001 - a missing/failed version must not kill the plot
            print(f"    [{ver}] skipped: {exc}")
    if "current" not in got:
        raise RuntimeError("current build did not run")

    D = float(io.read_conf_value(open(paths.accuracy("free_diffusion.conf")).read(), "diffusivity"))
    means = {ver: _shell_means(ss0, sc)[1] for ver, ss0 in got.items()}
    bs = sorted(means["current"])
    ora = metrics.per_shell_mean(oracles.free_exp(sc.bvals, D)[~sc.b0_mask], sc.bvals[~sc.b0_mask])

    os.makedirs(out_dir, exist_ok=True)
    cols = [(f"S/S0 {v}", means[v]) for v in means] + [("exp(−bD)", ora)]
    _write_table(os.path.join(out_dir, "free_diffusion.md"),
                 f"Cross-version free diffusion — current vs v1.5 vs v2.0 (N={n})", bs, cols)
    worst = {v: max(abs(means[v][b] - ora[b]) for b in bs) for v in means}
    summary = "versions: " + ", ".join(f"{v}(|Δoracle|≤{worst[v]:.4f})" for v in means)

    plt = _plt()
    if plt is None:
        return summary + "  (no matplotlib -> md only)"
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.2, 7.0),
                                  gridspec_kw=dict(height_ratios=[3, 1.2]), sharex=True)
    ax.plot(bs, [ora[b] for b in bs], "-k", lw=2, label="exp(−bD)")
    for ver in means:
        mk, col = markers.get(ver, ("x", "tab:gray"))
        ax.plot(bs, [means[ver][b] for b in bs], mk, ms=7, color=col, label=f"MC {ver}")
    ax.set_ylabel("S/S0 (per-shell mean)")
    ax.set_title(f"Cross-version: free diffusion, current vs v1.5 vs v2.0  (N={n})")
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend()

    axr.axhline(0, color="k", lw=0.6)
    for ver in means:
        if ver == "current":
            continue
        mk, col = markers.get(ver, ("x", "tab:gray"))
        axr.plot(bs, [means["current"][b] - means[ver][b] for b in bs], mk + "-", color=col,
                 label=f"current − {ver}")
    axr.set_xlabel("b (s/mm²)")
    axr.set_ylabel("Δ S/S0")
    axr.grid(True, ls=":", alpha=0.4)
    axr.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "free_diffusion.png"), dpi=120)
    plt.close(fig)
    return summary


# (tier, experiment-name, adapter). The tier is the suite's CTest rung; the folder is
# tests/output/<tier>/<name>/. `gpu` is not a tier -- every adapter runs CPU + GPU.
EXPERIMENTS = [
    ("tier1_smoke",    "free_diffusion",       exp_free),
    ("tier3_analytic", "restricted_sphere",    exp_sphere),
    ("tier3_analytic", "restricted_cylinder",  exp_cylinder),
    ("tier3_analytic", "cross_representation", exp_cross_rep),
    ("tier3_analytic", "packing_restricted",   exp_packing),
    ("tier3_analytic", "packing_exchange",     exp_exchange),
    ("tier3_analytic", "mesh_cell",            exp_mesh),
    ("tier3_analytic", "time_dependence",      exp_time_dependence),
    ("tier4_crossversion", "free_diffusion",   exp_crossversion),
]

# Tiers with no per-b overlay plot yet -- a README placeholder documents what belongs there.
SCAFFOLD_TIERS = {
    "tier2_fast": "Fast-regression tier (< ~8 min, CPU, deterministic): golden-master signals "
                  "(bit-exact -> a pass/fail, not an overlay), same-seed reproducibility, Powles "
                  "crossing counters, per-obstacle properties, waveform ≡ PGSE. The plottable ones "
                  "(waveform-vs-PGSE overlay, Powles p̂ vs prob) are not wired into this tool yet.",
    "tier5_heavy": "Heavy tier: WM two-mesh and ~50M-triangle meshes, validated by mass/ECVF/memory "
                  "invariants (no analytic oracle). On-demand — scaffold, no plots yet.",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=None, help="simulator binary (default: repo root)")
    ap.add_argument("--n", type=int, default=10000, help="walkers per run (both CPU and GPU)")
    ap.add_argument("--out", default=os.path.join(ROOT, "tests", "output"),
                    help="output folder (git-ignored; default tests/output/)")
    ap.add_argument("--only", default=None,
                    help="run a single experiment by name (e.g. packing_restricted) or a whole tier (e.g. tier3_analytic)")
    a = ap.parse_args()

    binary = paths.binary_path(a.bin)
    if not os.path.exists(binary):
        print(f"FAIL: simulator binary not found: {binary}", file=sys.stderr)
        return 1

    todo = [(tier, name, fn) for tier, name, fn in EXPERIMENTS
            if not a.only or a.only in (name, tier)]
    if not todo:
        names = [n for _, n, _ in EXPERIMENTS] + sorted({t for t, _, _ in EXPERIMENTS})
        print(f"no experiment matches --only {a.only}. choices: {names}")
        return 1

    results = []
    for tier, name, fn in todo:
        out_dir = os.path.join(a.out, tier, name)
        print(f"\n=== {tier}/{name} ===")
        try:
            summary = fn(binary, a.n, out_dir)
            print(f"  {summary}")
            results.append((tier, name, "ok", summary))
        except Exception as exc:  # noqa: BLE001 - one experiment must not kill the rest
            print(f"  ERROR: {exc}")
            traceback.print_exc()
            results.append((tier, name, "error", str(exc)))

    # scaffold-tier placeholders (so the tier folder exists and is self-documenting)
    for tier, blurb in SCAFFOLD_TIERS.items():
        td = os.path.join(a.out, tier)
        os.makedirs(td, exist_ok=True)
        with open(os.path.join(td, "README.md"), "w") as fh:
            fh.write(f"# {tier}\n\n{blurb}\n")

    # index, grouped by tier
    idx = ["# Validation plots — CPU vs GPU vs analytic, by test tier\n",
           f"Binary: `{binary}`  |  N = {a.n} walkers per run\n",
           "`gpu` is a modifier, not a tier: every plot overlays the CPU and GPU runs, so the "
           "CPU/GPU comparison is in each figure's residual panel. Each experiment folder holds "
           "a PNG (top: per-shell S/S0 with the analytic line + CPU/GPU markers; bottom: "
           "residuals incl. CPU−GPU) and a markdown data table.\n"]
    tiers_order = ["tier1_smoke", "tier2_fast", "tier3_analytic", "tier4_crossversion", "tier5_heavy"]
    for tier in tiers_order:
        rows = [r for r in results if r[0] == tier]
        idx.append(f"\n## {tier}")
        if not rows and tier in SCAFFOLD_TIERS:
            idx.append(f"_scaffold_ — {SCAFFOLD_TIERS[tier]}")
            continue
        for _, name, status, summary in rows:
            mark = "✅" if status == "ok" else "⚠️"
            idx.append(f"- {mark} **{name}** — {summary}")
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "index.md"), "w") as fh:
        fh.write("\n".join(idx) + "\n")
    n_ok = sum(1 for r in results if r[2] == "ok")
    print(f"\n{n_ok}/{len(results)} experiments plotted into {a.out}/ (by tier)")
    return 0 if n_ok else 1


if __name__ == "__main__":
    sys.exit(main())
