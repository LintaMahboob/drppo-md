#!/usr/bin/env python3
"""Combine the per-protein results into the summary table and figures.

Reads  results/<name>_{timeseries,rmsf,his_distances}.csv
Writes results/combined/
    summary_table.csv       mean ± SD of RMSD, Rg and mean RMSF
    stability.png           RMSD / Rg / RMSF for all proteins
    histidine_heatmap.png   His NE2-NE2 distances (upper = mean, lower = start)
    histidine_source_data.csv  numbers behind the heatmap

Example:
    python scripts/combine_results.py
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import REPO, load_config


def fmt(series):
    return f"{series.mean():.2f} ± {series.std():.2f}"


def summary_table(names, ts, rm):
    rows = [{"Protein": n,
             "Backbone RMSD (Å)": fmt(ts[n].rmsd_bb),
             "Core RMSD (Å)": fmt(ts[n].rmsd_core),
             "Rg (Å)": fmt(ts[n].rg),
             "Mean RMSF (Å)": f"{rm[n].rmsf.mean():.2f}"} for n in names]
    return pd.DataFrame(rows)


def stability_figure(names, ts, rm, path):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for n in names:
        x = ts[n].time_ps / 1000
        ax[0, 0].plot(x, ts[n].rmsd_core, label=n)
        ax[0, 1].plot(x, ts[n].rmsd_bb, label=n)
        ax[0, 2].plot(x, ts[n].rg, label=n)
    ax[0, 0].set(title="Core backbone RMSD", xlabel="Time (ns)", ylabel="RMSD (Å)")
    ax[0, 0].legend()
    ax[0, 1].set(title="Backbone RMSD", xlabel="Time (ns)", ylabel="RMSD (Å)")
    ax[0, 2].set(title="Radius of gyration", xlabel="Time (ns)", ylabel="Rg (Å)")
    for i, n in enumerate(names):
        ax[1, i].plot(rm[n].resid, rm[n].rmsf)
        ax[1, i].set(title=f"{n} RMSF", xlabel="Residue", ylabel="RMSF (Å)", ylim=(0, 5))
        ax[1, i].text(0.98, 0.95, f"max {rm[n].rmsf.max():.1f} Å",
                      transform=ax[1, i].transAxes, ha="right", va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def heatmap_figure(names, hd, path):
    fig, ax = plt.subplots(1, len(names), figsize=(5.7 * len(names), 5.5))
    ax = np.atleast_1d(ax)
    im = None
    for i, n in enumerate(names):
        df = hd[n]
        his = sorted({int(v) for p in df.pair for v in p.split("-")})
        M = pd.DataFrame(np.nan, index=his, columns=his)
        for _, r in df.iterrows():
            a, b = map(int, r.pair.split("-"))
            M.loc[a, b] = r["mean"]     # upper triangle
            M.loc[b, a] = r["start"]    # lower triangle
        im = ax[i].imshow(M, cmap="viridis_r", vmin=3, vmax=14)
        for r_ in range(len(his)):
            for c_ in range(len(his)):
                v = M.iloc[r_, c_]
                if not np.isnan(v):
                    ax[i].text(c_, r_, f"{v:.1f}", ha="center", va="center", fontsize=8,
                               color="white" if v > 9 else "black")
        ax[i].set_xticks(range(len(his)))
        ax[i].set_xticklabels(his, rotation=90)
        ax[i].set_yticks(range(len(his)))
        ax[i].set_yticklabels(his)
        ax[i].set_title(f"{n}: upper = mean, lower = start")
    fig.colorbar(im, ax=ax, label="NE2-NE2 (Å)", shrink=0.8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", help="path to a proteins.yaml")
    ap.add_argument("--results", default=str(REPO / "results"),
                    help="folder with the per-protein CSV files (default: results/)")
    ap.add_argument("--proteins", nargs="+",
                    help="proteins to include (default: all in the config)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    names = args.proteins or list(cfg["proteins"])
    res = Path(args.results)
    out = res / "combined"
    out.mkdir(parents=True, exist_ok=True)

    ts = {n: pd.read_csv(res / f"{n}_timeseries.csv") for n in names}
    rm = {n: pd.read_csv(res / f"{n}_rmsf.csv") for n in names}
    hd = {n: pd.read_csv(res / f"{n}_his_distances.csv") for n in names}

    summary = summary_table(names, ts, rm)
    print(summary.to_string(index=False))
    summary.to_csv(out / "summary_table.csv", index=False)

    stability_figure(names, ts, rm, out / "stability.png")
    heatmap_figure(names, hd, out / "histidine_heatmap.png")

    frames = []
    for n in names:
        df = hd[n].copy()
        df.insert(0, "protein", n)
        frames.append(df)
    pd.concat(frames).to_csv(out / "histidine_source_data.csv", index=False)
    print(f"Wrote combined outputs to {out}")


if __name__ == "__main__":
    main()
