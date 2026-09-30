#!/usr/bin/env python3
"""Analyse one apo MD run (first N ps only).

Steps
    1. gmx trjconv: make molecules whole and centre on the protein
       (REQUIRED: without it the RMSD is wrong by several angstroms).
    2. Per frame: backbone RMSD, core RMSD, radius of gyration, His NE2-NE2
       distances. Per residue: C-alpha RMSF. Reference = frame 0 (the first
       frame after equilibration).
    3. Write CSV files and a 4-panel PNG.

MDAnalysis cannot read GROMACS 2026 .tpr files, so a .gro file is used as
the topology.

Example:
    python scripts/analyze_md.py DrPPO1
"""
import argparse
import shutil
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import calc_bonds

from common import REPO, get_protein, load_config


def make_whole(rundir, xtc):
    """Write md_whole.xtc / md_whole.gro (PBC made whole, protein centred)."""
    if shutil.which("gmx") is None:
        sys.exit("'gmx' not found; needed for trjconv.")
    # stdin: 1st group = centring, 2nd group = output
    for src, dst in [(xtc, "md_whole.xtc"), ("md.gro", "md_whole.gro")]:
        res = subprocess.run(
            ["gmx", "trjconv", "-s", "md.tpr", "-f", src, "-o", dst,
             "-pbc", "mol", "-center"],
            cwd=rundir, input="Protein\nSystem\n", text=True, capture_output=True)
        if res.returncode != 0:
            sys.exit(f"trjconv failed for {src}:\n{res.stderr[-1500:]}")


def analyse(gro, xtc, his, core_sel, tmax):
    u = mda.Universe(str(gro), str(xtc))
    bb = u.select_atoms("backbone")
    core = u.select_atoms(core_sel) if core_sel else None
    if core is not None and len(core) == 0:
        sys.exit(f"Core selection '{core_sel}' matched no atoms.")
    ca = u.select_atoms("name CA")
    prot = u.select_atoms("protein")
    ne2 = {h: u.select_atoms(f"resid {h} and name NE2") for h in his}
    for h, sel in ne2.items():
        if len(sel) != 1:
            sys.exit(f"Expected one NE2 atom for His{h}, found {len(sel)}. "
                     "Check the residue numbering in config/proteins.yaml.")
    pairs = list(combinations(his, 2))

    u.trajectory[0]
    ref_bb = bb.positions - bb.center_of_mass()
    ref_core = core.positions - core.center_of_mass() if core is not None else None

    t, rb, rc, rg, cas = [], [], [], [], []
    dist = {p: [] for p in pairs}
    for ts in u.trajectory:
        if ts.time > tmax:
            break
        t.append(ts.time)
        _bb_c = bb.center_of_mass()
        R, r = rotation_matrix(bb.positions - _bb_c, ref_bb)
        rb.append(r)
        if core is not None:
            _, r = rotation_matrix(core.positions - core.center_of_mass(), ref_core)
            rc.append(r)
        else:
            rc.append(np.nan)
        cas.append((ca.positions - _bb_c) @ R.T)   # C-alpha in the aligned frame
        rg.append(prot.radius_of_gyration())
        for a, b in pairs:
            dist[(a, b)].append(float(calc_bonds(
                ne2[a].positions[0], ne2[b].positions[0], box=ts.dimensions)))

    t, rb, rc, rg = map(np.array, (t, rb, rc, rg))
    cas = np.array(cas)
    rmsf = np.sqrt(((cas - cas.mean(0)) ** 2).sum(2).mean(0))
    return dict(t=t, rb=rb, rc=rc, rg=rg, rmsf=rmsf, resids=ca.resids,
                dist=dist, pairs=pairs)


def save(name, res, outdir):
    t, rb, rc, rg, rmsf = res["t"], res["rb"], res["rc"], res["rg"], res["rmsf"]
    dist, pairs = res["dist"], res["pairs"]
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"{name}: {len(t)} frames, {t[0]:.0f}-{t[-1]:.0f} ps")
    print(f"  RMSD backbone {rb.mean():.2f} ± {rb.std():.2f} (max {rb.max():.2f}) A")
    print(f"  RMSD core     {np.nanmean(rc):.2f} ± {np.nanstd(rc):.2f} A")
    print(f"  Rg            {rg.mean():.2f} ± {rg.std():.2f} A")
    print(f"  RMSF mean {rmsf.mean():.2f} A, max {rmsf.max():.2f} A "
          f"at resid {res['resids'][rmsf.argmax()]}")

    rows = []
    for a, b in pairs:
        x = np.array(dist[(a, b)])
        rows.append(dict(pair=f"{a}-{b}", start=x[0], mean=x.mean(), sd=x.std(), max=x.max()))
        print(f"  His {a}-{b}: start {x[0]:.1f}  mean {x.mean():.1f} ± {x.std():.1f}  max {x.max():.1f} A")

    pd.DataFrame(rows).round(2).to_csv(outdir / f"{name}_his_distances.csv", index=False)
    pd.DataFrame(dict(time_ps=t, rmsd_bb=rb, rmsd_core=rc, rg=rg)).to_csv(
        outdir / f"{name}_timeseries.csv", index=False)
    pd.DataFrame(dict(resid=res["resids"], rmsf=rmsf)).to_csv(
        outdir / f"{name}_rmsf.csv", index=False)
    pd.DataFrame({f"{a}-{b}": dist[(a, b)] for a, b in pairs},
                 index=pd.Index(t, name="time_ps")).to_csv(outdir / f"{name}_his_timeseries.csv")

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(t / 1000, rb, label="backbone")
    ax[0, 0].plot(t / 1000, rc, label="core")
    ax[0, 0].set(xlabel="Time (ns)", ylabel="RMSD (Å)")
    ax[0, 0].legend()
    ax[0, 1].plot(t / 1000, rg)
    ax[0, 1].set(xlabel="Time (ns)", ylabel="Rg (Å)")
    ax[1, 0].plot(res["resids"], rmsf)
    ax[1, 0].set(xlabel="Residue", ylabel="RMSF (Å)")
    cmap = plt.get_cmap("tab10")
    for i, (a, b) in enumerate(pairs):
        ax[1, 1].plot(t / 1000, dist[(a, b)], color=cmap(i % 10), lw=1.2, label=f"{a}-{b}")
    ax[1, 1].set(xlabel="Time (ns)", ylabel="His NE2-NE2 (Å)")
    ax[1, 1].legend(fontsize=8, ncol=2)
    fig.suptitle(f"{name} (apo, 0-{t[-1] / 1000:.1f} ns, single replicate)")
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_md.png", dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="protein name defined in the config")
    ap.add_argument("--config", help="path to a proteins.yaml")
    ap.add_argument("--rundir", help="run directory (default: runs/<name>)")
    ap.add_argument("--xtc", default="md.xtc",
                    help="trajectory inside the run directory (default: md.xtc)")
    ap.add_argument("--outdir", default=str(REPO / "results"),
                    help="where to write CSV/PNG files (default: results/)")
    ap.add_argument("--tmax", type=float, default=None,
                    help="analyse frames up to this time in ps "
                         "(default: analysis_window_ps from the config)")
    ap.add_argument("--skip-whole", action="store_true",
                    help="reuse existing md_whole.* files")
    args = ap.parse_args()

    cfg = load_config(args.config)
    prot = get_protein(cfg, args.name)
    rundir = Path(args.rundir) if args.rundir else REPO / "runs" / args.name
    tmax = args.tmax if args.tmax is not None else cfg["analysis_window_ps"]
    if not prot.get("core"):
        print("WARNING: no 'core' selection set for this protein; core RMSD will be NaN.")

    if not args.skip_whole:
        make_whole(rundir, args.xtc)
    res = analyse(rundir / "md_whole.gro", rundir / "md_whole.xtc",
                  prot["his"], prot.get("core"), tmax)
    save(args.name, res, Path(args.outdir))


if __name__ == "__main__":
    main()
