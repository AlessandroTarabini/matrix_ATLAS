#!/usr/bin/env python3
"""
Plot the three correlation estimates vs N_sig / N_bkg from results_*.json
files produced by run.py.

Usage:
  python3 plotter.py
  python3 plotter.py --indir .
  python3 plotter.py --rho-sig 0.5 --rho-bkg 0.2
"""

import argparse
import glob
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def format_suffix_number(x):
    if float(x).is_integer():
        s = str(int(round(float(x))))
    else:
        s = f"{float(x):g}"
    return s.replace("-", "m").replace(".", "p")


def load_results(indir):
    """Load all results_*.json files and return a list of records."""
    paths = sorted(glob.glob(os.path.join(indir, "results_*.json")))
    records = []
    for path in paths:
        with open(path) as f:
            data = json.load(f)
        cfg = data["config"]
        records.append({
            "path": path,
            "nsig": float(cfg["nsig_tot"]),
            "nbkg": float(cfg["nbkg_tot"]),
            "rho_sig": float(cfg["rho_sig"]),
            "rho_bkg": float(cfg["rho_bkg"]),
            "sb_ratio": float(cfg["nsig_tot"]) / float(cfg["nbkg_tot"]),
            "rho_truth": float(data["rho_truth"]),
            "rho_truth_err": float(data["rho_truth_err"]),
            "rho_bootfit": float(data["bootstrap"]["rho_bootfit"]),
            "rho_bootfit_err": float(data["bootstrap_errs"]["rho_bootfit"]),
            "rho_sideband": float(data["bootstrap"]["rho_sideband"]),
            "rho_sideband_err": float(data["bootstrap_errs"]["rho_sideband"]),
        })
    return records


def filter_records(records, rho_sig=None, rho_bkg=None):
    out = records
    if rho_sig is not None:
        out = [r for r in out if math.isclose(r["rho_sig"], rho_sig, abs_tol=1e-9)]
    if rho_bkg is not None:
        out = [r for r in out if math.isclose(r["rho_bkg"], rho_bkg, abs_tol=1e-9)]
    return out


def group_by_correlations(records):
    """Split into series with the same (rho_sig, rho_bkg)."""
    groups = {}
    for r in records:
        key = (r["rho_sig"], r["rho_bkg"])
        groups.setdefault(key, []).append(r)
    for key in groups:
        groups[key].sort(key=lambda r: r["sb_ratio"])
    return groups


def plot_series(records, outdir):
    """Make one PDF for one (rho_sig, rho_bkg) series."""
    rho_sig = records[0]["rho_sig"]
    rho_bkg = records[0]["rho_bkg"]

    x = np.array([r["sb_ratio"] for r in records])
    truth = np.array([r["rho_truth"] for r in records])
    truth_err = np.array([r["rho_truth_err"] for r in records])
    bootfit = np.array([r["rho_bootfit"] for r in records])
    bootfit_err = np.array([r["rho_bootfit_err"] for r in records])
    sideband = np.array([r["rho_sideband"] for r in records])
    sideband_err = np.array([r["rho_sideband_err"] for r in records])

    fig, ax = plt.subplots(1, 1, figsize=(8.5, 5.0))

    ax.errorbar(x, truth, yerr=truth_err, fmt="o-", color="tab:blue",
                capsize=3, ms=5, lw=1.2, label="toys + fit (truth)")
    ax.errorbar(x, bootfit, yerr=bootfit_err, fmt="s-", color="tab:red",
                capsize=3, ms=5, lw=1.2, label="bootstrap + refit")
    ax.errorbar(x, sideband, yerr=sideband_err, fmt="^-", color="tab:green",
                capsize=3, ms=5, lw=1.2, label="bootstrap + sideband count")

    ax.axhline(rho_bkg, color="tab:green", ls=":", lw=1.2, alpha=0.9,
               label=rf"input $\rho_{{\mathrm{{bkg}}}}$ = {rho_bkg:g}")
    ax.axhline(rho_sig, color="tab:purple", ls=":", lw=1.2, alpha=0.9,
               label=rf"input $\rho_{{\mathrm{{sig}}}}$ = {rho_sig:g}")

    ax.set_xscale("log")
    ax.set_xlabel(r"$N_S^{\mathrm{tot}} / N_B^{\mathrm{tot}}$")
    ax.set_ylabel("correlation")
    ax.set_title(rf"Correlation vs $S/B$  ($\rho_{{\mathrm{{sig}}}}={rho_sig:g}$, "
                 rf"$\rho_{{\mathrm{{bkg}}}}={rho_bkg:g}$)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()

    suffix = (f"_{format_suffix_number(rho_sig)}"
              f"_{format_suffix_number(rho_bkg)}")
    path = os.path.join(outdir, f"correlation_vs_sb{suffix}.pdf")
    fig.savefig(path)
    plt.close(fig)
    return path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--indir", default=os.path.dirname(os.path.abspath(__file__)),
                   help="directory containing results_*.json files")
    p.add_argument("--outdir", default=None,
                   help="output directory for the plot (default: same as --indir)")
    p.add_argument("--rho-sig", type=float, default=None,
                   help="keep only results with this rho_sig")
    p.add_argument("--rho-bkg", type=float, default=None,
                   help="keep only results with this rho_bkg")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = args.outdir or args.indir

    records = load_results(args.indir)
    if not records:
        raise SystemExit(f"no results_*.json found in {args.indir}")

    records = filter_records(records, args.rho_sig, args.rho_bkg)
    if not records:
        raise SystemExit("no results left after filtering on rho_sig / rho_bkg")

    groups = group_by_correlations(records)
    print(f"found {len(records)} result files in {len(groups)} series")
    for (rho_sig, rho_bkg), series in sorted(groups.items()):
        print(f"  rho_sig={rho_sig:g}, rho_bkg={rho_bkg:g}: "
              f"{len(series)} points, "
              f"S/B in [{series[0]['sb_ratio']:.2e}, {series[-1]['sb_ratio']:.2e}]")
        path = plot_series(series, outdir)
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
