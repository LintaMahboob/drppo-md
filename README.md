# drppo-md

Reproducible GROMACS workflow for the short apo molecular dynamics (MD) analyses of three yam (*Dioscorea rotundata*) polyphenol oxidase models, **DrPPO1, DrPPO2 and DrPPO7**, used to assess fold stability and the behaviour of the active-site histidine cluster.

Supporting code for: [PAPER TITLE], [AUTHORS], [JOURNAL, YEAR]. Data: [FIGSHARE DOI].

## Scope and limitations

- **Analysis window: 0-2 ns**, 201 frames per protein (one frame every 10 ps), **single replicate**, **apo** (no copper).
- The production input (`mdp/md.mdp`) is set to 10 ns, but only the first 2 ns are analysed (`analysis_window_ps` in `config/proteins.yaml`). In the original runs, the DrPPO1 and DrPPO7 trajectories were incomplete, so a common 0-2 ns window was used for all three. These are not 10 ns results.
- The apo simulations alone do not separate DrPPO1 from the other models.
- Results are illustrative structural-stability checks, not converged sampling.

## Repository layout

```
config/proteins.yaml      per-protein settings (trim range, histidines, core region)
mdp/                      GROMACS parameter files (ions, em, nvt, npt, md)
models/                   put the input PDB models here (not included)
scripts/run_md.py         prepare + run the MD workflow for one protein
scripts/analyze_md.py     trajectory analysis for one protein
scripts/combine_results.py  summary table and combined figures
environment.yml           conda environment
```

## Installation

```bash
conda env create -f environment.yml
conda activate drppo-md
gmx --version
```

The paper used **GROMACS 2026.3**. Other recent versions should work, but results will not be bit-identical.

For GPU runs, install a CUDA build of GROMACS (for example `conda install -c conda-forge "gromacs=*=nompi_cuda*"`, setting `CONDA_OVERRIDE_CUDA` if needed) and pass `--gpu` to `run_md.py`.

## Input models

Place the ColabFold/AlphaFold models in `models/` with the names given in `config/proteins.yaml` (for example `models/DrPPO1_rank_001.pdb`). The script keeps residues `first`-`last`, renames all `HIS` to `HSD`, and resets the B-factor column to 1.00. Residue numbering stays the original model numbering.

## Usage

```bash
# 1. Build the system and run MD (repeat for DrPPO2 and DrPPO7)
python scripts/run_md.py DrPPO1 --gpu

# 2. Analyse the first 2 ns (repeat for each protein)
python scripts/analyze_md.py DrPPO1

# 3. Combine all three proteins
python scripts/combine_results.py
```

Useful options of `run_md.py`:

| Option | Purpose |
|---|---|
| `--gpu` | offload nonbonded and PME to the GPU |
| `--maxh H` | stop the production run after H hours |
| `--md-nsteps N` | change production length (default 5,000,000 steps = 10 ns) |
| `--seed N` | fix the velocity-generation seed (default random, as in the paper) |
| `--resume` | continue an interrupted production run from `md.cpt` |
| `--outdir DIR` | run directory parent (default `runs/`) |

Stages whose output file already exists are skipped, so an interrupted setup can simply be re-run. The CHARMM36 (jul2022) force field is downloaded automatically into `forcefield/`; if the download fails, get it from the MacKerell lab website and place the `charmm36-jul2022.ff` folder there.

## Simulation protocol

GROMACS, CHARMM36 (jul2022), TIP3P; all His as HSD; dodecahedral box with 1.2 nm margin; 0.15 M NaCl (neutralised), no copper. Energy minimisation, then 100 ps NVT and 200 ps NPT with position restraints, then unrestrained production MD: 2 fs time step, V-rescale (300 K, tau 0.1 ps), Parrinello-Rahman (1 bar, tau 2 ps), PME with 1.2 nm cutoff, van der Waals cutoff 1.2 nm with force-switch from 1.0 nm, LINCS on bonds to hydrogen, frames every 10 ps. Full details are in `mdp/`.

## Analysis

`analyze_md.py` first makes molecules whole and centres the protein (`gmx trjconv -pbc mol -center`). **This step is essential:** without it the DrPPO1 RMSD was 17.6 Å instead of 2.64 Å. It then computes, against frame 0 (the first frame after equilibration):

- backbone RMSD and "core" backbone RMSD (core residue selection per protein in the config)
- radius of gyration
- C-alpha RMSF
- His NE2-NE2 distances for all histidine pairs listed in the config

Outputs go to `results/`:

| File | Content |
|---|---|
| `<name>_timeseries.csv` | `time_ps`, `rmsd_bb`, `rmsd_core`, `rg` |
| `<name>_rmsf.csv` | `resid`, `rmsf` |
| `<name>_his_distances.csv` | `pair`, `start` (frame 0), `mean`, `sd`, `max` |
| `<name>_his_timeseries.csv` | per-frame distance of each pair |
| `<name>_md.png` | four-panel overview |
| `combined/` | `summary_table.csv`, `stability.png`, `histidine_heatmap.png`, `histidine_source_data.csv` |

Units: Å and ps.

## Known issues and tips

- **Write trajectories to a local disk.** On Google Colab, write to `/content` and copy to Drive at the end. Writing the trajectory directly to a mounted Drive folder truncated `md.xtc` in one case, and check file sizes after copying.
- MDAnalysis cannot read GROMACS 2026 `.tpr` files (tpx version 138), so the scripts use `.gro` files as topology.
- "Start" distances are taken after equilibration, so they differ from the distances in the AlphaFold models.
- GPU runs are not bit-reproducible, and the default random seed means repeated runs will differ in detail.
- The published DrPPO2 analysis simulated His137 but excluded it from the tables and figures, since it is not part of the five-histidine cluster. The config analyses only the five cluster histidines, which gives identical values for the others.
- Core-RMSD residue ranges for DrPPO2 and DrPPO7 are not set yet in `config/proteins.yaml`. Fill them in, or leave them as `null` to skip the core RMSD.

## Running on Google Colab

```python
!pip -q install condacolab
import condacolab; condacolab.install()          # the runtime restarts
!mamba install -y -c conda-forge mdanalysis pyyaml pandas matplotlib "gromacs=*=nompi_cuda*"
!git clone https://github.com/[USER]/drppo-md && cd drppo-md
```

Then run the commands from the Usage section. Reinstall GROMACS and MDAnalysis after every runtime reset, and use `--maxh` and `--resume` to work around session limits.

## Citation

If you use this workflow, please cite the paper above and the tools it relies on: GROMACS (Abraham et al., 2015, *SoftwareX* 1-2, 19-25), CHARMM36 (Best et al., 2012, *J. Chem. Theory Comput.* 8, 3257-3273), TIP3P (Jorgensen et al., 1983, *J. Chem. Phys.* 79, 926-935) and MDAnalysis (Michaud-Agrawal et al., 2011, *J. Comput. Chem.* 32, 2319-2327).

## License

[MIT / choose a license]
