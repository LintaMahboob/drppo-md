#!/usr/bin/env python3
"""Prepare a protein model and run the apo MD workflow with GROMACS.

Workflow (identical for every protein):
    trim model, all His -> HSD
    pdb2gmx (CHARMM36 jul2022, TIP3P)
    dodecahedral box (1.2 nm margin), solvate, 0.15 M NaCl
    energy minimisation -> NVT (100 ps) -> NPT (200 ps) -> production MD

Stages whose output file already exists are skipped, so the script can be
re-run after an interruption. If the production run was interrupted, re-run
with --resume so that it continues from md.cpt.

Example:
    python scripts/run_md.py DrPPO1 --gpu
"""
import argparse
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

from common import REPO, get_protein, load_config

FF_NAME = "charmm36-jul2022.ff"
FF_URL = ("http://mackerell.umaryland.edu/download.php?filename="
          "CHARMM_ff_params_files/charmm36-jul2022.ff.tgz")
MDP_FILES = ["ions.mdp", "em.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]


def ensure_forcefield():
    """Return the CHARMM36 force-field directory, downloading it if needed."""
    ff = REPO / "forcefield" / FF_NAME
    if ff.exists():
        return ff
    ff.parent.mkdir(exist_ok=True)
    tgz = ff.parent / "ff.tgz"
    print(f"Downloading {FF_NAME} ...")
    req = urllib.request.Request(FF_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp, open(tgz, "wb") as out:
        shutil.copyfileobj(resp, out)
    with tarfile.open(tgz) as tar:
        tar.extractall(ff.parent)
    tgz.unlink()
    if not ff.exists():
        sys.exit(f"Download finished but {ff} was not found. Download the force "
                 "field manually from the MacKerell lab site and place the "
                 f"{FF_NAME} folder in {ff.parent}.")
    return ff


def trim_model(src, dst, first, last):
    """Keep residues first..last, rename HIS -> HSD, reset B-factors to 1.00."""
    out = []
    with open(src) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                if not first <= int(line[22:26]) <= last:
                    continue
                line = line[:17] + line[17:20].replace("HIS", "HSD") + line[20:]
                line = line[:60] + "  1.00" + line[66:]   # pLDDT -> 1.00
            out.append(line)
    Path(dst).write_text("".join(out))


def set_mdp(text, key, value):
    """Replace the value of an existing mdp option."""
    pattern = re.compile(rf"^{re.escape(key)}\s*=.*$", re.M)
    if not pattern.search(text):
        sys.exit(f"Option '{key}' not found in mdp file")
    return pattern.sub(f"{key:<22} = {value}", text)


class Gmx:
    """Run gmx commands in a working directory and log their output."""

    def __init__(self, workdir):
        self.workdir = Path(workdir)
        self.log = self.workdir / "run.log"

    def __call__(self, args, stdin=None):
        cmd = ["gmx"] + [str(a) for a in args]
        print("  $", " ".join(cmd))
        res = subprocess.run(cmd, cwd=self.workdir, input=stdin,
                             text=True, capture_output=True)
        with open(self.log, "a") as fh:
            fh.write(f"\n### {' '.join(cmd)}\n{res.stdout}\n{res.stderr}")
        if res.returncode != 0:
            print(res.stderr[-1500:])
            sys.exit(f"FAILED: {' '.join(cmd)}  (full output in {self.log})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="protein name defined in the config (e.g. DrPPO1)")
    ap.add_argument("--config", help="path to a proteins.yaml (default: config/proteins.yaml)")
    ap.add_argument("--outdir", default=str(REPO / "runs"),
                    help="parent folder for run directories (default: runs/). "
                         "Use a LOCAL disk; network/cloud-synced folders can "
                         "truncate trajectories.")
    ap.add_argument("--gpu", action="store_true", help="offload nonbonded and PME to the GPU")
    ap.add_argument("--maxh", type=float, default=None,
                    help="stop the production run after this many hours")
    ap.add_argument("--md-nsteps", type=int, default=None,
                    help="override production length (default 5,000,000 = 10 ns)")
    ap.add_argument("--seed", type=int, default=None,
                    help="fixed velocity-generation seed (default: random, as in the paper)")
    ap.add_argument("--resume", action="store_true",
                    help="continue an interrupted production run from md.cpt")
    args = ap.parse_args()

    if shutil.which("gmx") is None:
        sys.exit("'gmx' not found. Install GROMACS (see environment.yml) and retry.")

    cfg = load_config(args.config)
    prot = get_protein(cfg, args.name)
    model = REPO / prot["model"]
    if not model.exists():
        sys.exit(f"Model not found: {model}\nPlace the PDB there or edit config/proteins.yaml.")

    wd = Path(args.outdir) / args.name
    wd.mkdir(parents=True, exist_ok=True)
    gmx = Gmx(wd)

    # --- inputs: force field, mdp files, trimmed model --------------------
    ff = ensure_forcefield()
    if not (wd / FF_NAME).exists():
        shutil.copytree(ff, wd / FF_NAME)
    for name in MDP_FILES:
        text = (REPO / "mdp" / name).read_text()
        if name == "nvt.mdp" and args.seed is not None:
            text = set_mdp(text, "gen_seed", args.seed)
        if name == "md.mdp" and args.md_nsteps is not None:
            text = set_mdp(text, "nsteps", args.md_nsteps)
        (wd / name).write_text(text)
    trim = f"{args.name}_trim.pdb"
    trim_model(model, wd / trim, prot["first"], prot["last"])

    gpu = ["-nb", "gpu", "-pme", "gpu"] if args.gpu else []
    gpu_nb = ["-nb", "gpu"] if args.gpu else []

    # --- stages: (expected output, gmx arguments, stdin) ------------------
    stages = [
        ("conf.gro", ["pdb2gmx", "-f", trim, "-o", "conf.gro", "-p", "topol.top",
                      "-ff", "charmm36-jul2022", "-water", "tip3p", "-ignh"], None),
        ("box.gro", ["editconf", "-f", "conf.gro", "-o", "box.gro", "-c",
                     "-d", "1.2", "-bt", "dodecahedron"], None),
        ("solv.gro", ["solvate", "-cp", "box.gro", "-cs", "spc216.gro",
                      "-o", "solv.gro", "-p", "topol.top"], None),
        # -maxwarn 1: the un-neutralised system triggers a charge warning
        ("ions.tpr", ["grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol.top",
                      "-o", "ions.tpr", "-maxwarn", "1"], None),
        ("ions.gro", ["genion", "-s", "ions.tpr", "-o", "ions.gro", "-p", "topol.top",
                      "-pname", "NA", "-nname", "CL", "-neutral", "-conc", "0.15"], "SOL\n"),
        ("em.tpr", ["grompp", "-f", "em.mdp", "-c", "ions.gro", "-p", "topol.top",
                    "-o", "em.tpr"], None),
        ("em.gro", ["mdrun", "-deffnm", "em"] + gpu_nb, None),
        ("nvt.tpr", ["grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro",
                     "-p", "topol.top", "-o", "nvt.tpr"], None),
        ("nvt.gro", ["mdrun", "-deffnm", "nvt"] + gpu, None),
        ("npt.tpr", ["grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro",
                     "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr"], None),
        ("npt.gro", ["mdrun", "-deffnm", "npt"] + gpu, None),
        ("md.tpr", ["grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt",
                    "-p", "topol.top", "-o", "md.tpr"], None),
    ]
    for output, cmd, stdin in stages:
        if (wd / output).exists():
            print(f"skip {cmd[0]} ({output} exists)")
            continue
        gmx(cmd, stdin)

    # --- production --------------------------------------------------------
    prod = ["mdrun", "-deffnm", "md", "-cpt", "15"] + gpu
    if args.maxh:
        prod += ["-maxh", args.maxh]
    if args.resume and (wd / "md.cpt").exists():
        prod += ["-cpi", "md.cpt"]
    elif (wd / "md.gro").exists():
        print("skip production (md.gro exists; use --resume to continue an interrupted run)")
        return
    gmx(prod)
    print(f"{args.name}: production run finished. Output in {wd}")


if __name__ == "__main__":
    main()
