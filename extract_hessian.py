#!/usr/bin/env python3
## To plot the correlation matrix from the hessian matrix of the combined fit.
## python3 /afs/cern.ch/work/a/atarabin/matrix_ATLAS/extract_hessian.py --output-pdf /afs/cern.ch/work/a/atarabin/matrix_ATLAS/r_correlation_subset.pdf

import argparse
import sys
import ROOT

# Hardcoded axis order for labels starting with "r_".
# Edit this list directly to change plotting/printing order.
LABEL_ORDER = [
    "r_PTH_0p0_15p0",
    "r_PTH_15p0_30p0",
    "r_PTH_30p0_45p0",
    "r_PTH_45p0_80p0",
    "r_PTH_80p0_120p0",
    "r_PTH_120p0_200p0",
    "r_PTH_200p0_350p0",
    "r_PTH_350p0_10000p0",
    "r_NJ_0p0_1p0",
    "r_NJ_1p0_2p0",
    "r_NJ_2p0_3p0",
    "r_NJ_3p0_100p0",
    "r_PTJ0_0p0_30p0",
    "r_PTJ0_30p0_75p0",
    "r_PTJ0_75p0_120p0",
    "r_PTJ0_120p0_200p0",
    "r_PTJ0_200p0_10000p0",
]


def _axis_labels(axis):
    labels = []
    for i in range(1, axis.GetNbins() + 1):
        label = axis.GetBinLabel(i)
        labels.append(label.strip() if label else "")
    return labels


def _select_indices(labels, prefix):
    return [i for i, label in enumerate(labels, start=1) if label.startswith(prefix)]


def _apply_custom_order(indices, labels, requested):
    if not requested:
        return indices

    label_to_index = {labels[i - 1]: i for i in indices}

    ordered = []
    used = set()
    missing = []

    for label in requested:
        idx = label_to_index.get(label)
        if idx is None:
            missing.append(label)
            continue
        ordered.append(idx)
        used.add(idx)

    # Keep any non-specified labels in their original axis order.
    ordered.extend(i for i in indices if i not in used)

    if missing:
        print(
            "Warning: these labels were not found and were ignored: "
            + ", ".join(missing),
            file=sys.stderr,
        )

    return ordered


def _print_subset(hist, x_indices, y_indices, x_labels, y_labels):
    # Header
    header = ["row\\col"] + [x_labels[ix - 1] for ix in x_indices]
    print("\t".join(header))

    # Matrix rows
    for iy in y_indices:
        row_label = y_labels[iy - 1]
        values = [f"{hist.GetBinContent(ix, iy): .10f}" for ix in x_indices]
        print("\t".join([row_label] + values))


def _build_subset_hist(hist, x_indices, y_indices, x_labels, y_labels):
    # Reverse Y order for plotting so the top row follows the requested order.
    y_plot_indices = list(reversed(y_indices))

    subset = ROOT.TH2F(
        "h_correlation_subset",
        "h_correlation subset",
        len(x_indices),
        0,
        len(x_indices),
        len(y_plot_indices),
        0,
        len(y_plot_indices),
    )

    for x_bin, ix in enumerate(x_indices, start=1):
        subset.GetXaxis().SetBinLabel(x_bin, x_labels[ix - 1])
    for y_bin, iy in enumerate(y_plot_indices, start=1):
        subset.GetYaxis().SetBinLabel(y_bin, y_labels[iy - 1])

    for y_bin, iy in enumerate(y_plot_indices, start=1):
        for x_bin, ix in enumerate(x_indices, start=1):
            subset.SetBinContent(x_bin, y_bin, hist.GetBinContent(ix, iy))

    subset.GetXaxis().LabelsOption("v")
    subset.GetXaxis().SetLabelSize(0.03)
    subset.GetYaxis().SetLabelSize(0.03)
    subset.GetZaxis().SetRangeUser(-1.0, 1.0)
    return subset


def _save_subset_pdf(subset_hist, output_pdf):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    canvas = ROOT.TCanvas("c_corr_subset", "Correlation subset", 1800, 1200)
    canvas.SetLeftMargin(0.26)
    canvas.SetBottomMargin(0.28)
    canvas.SetRightMargin(0.16)
    canvas.SetTopMargin(0.08)
    subset_hist.SetStats(0)
    subset_hist.Draw("COLZ TEXT")
    canvas.SaveAs(output_pdf)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file",
        default="/afs/cern.ch/work/a/atarabin/matrix_ATLAS/datacards/differentials/combined_fit/robustHesseTest.root",
        help="Path to ROOT file.",
    )
    parser.add_argument(
        "--hist",
        default="h_correlation",
        help="TH2 name inside ROOT file.",
    )
    parser.add_argument(
        "--prefix",
        default="r_",
        help="Keep labels starting with this prefix.",
    )
    parser.add_argument(
        "--output-pdf",
        default="",
        help="Optional output PDF path for the selected submatrix.",
    )
    args = parser.parse_args()

    root_file = ROOT.TFile.Open(args.file)
    if not root_file or root_file.IsZombie():
        print(f"Error: cannot open ROOT file '{args.file}'", file=sys.stderr)
        sys.exit(1)

    hist = root_file.Get(args.hist)
    if not hist:
        print(f"Error: histogram '{args.hist}' not found in '{args.file}'", file=sys.stderr)
        sys.exit(1)

    x_labels = _axis_labels(hist.GetXaxis())
    y_labels = _axis_labels(hist.GetYaxis())

    x_indices = _select_indices(x_labels, args.prefix)
    y_indices = _select_indices(y_labels, args.prefix)

    if not x_indices or not y_indices:
        print(f"No labels starting with '{args.prefix}' were found.")
        sys.exit(0)

    x_indices = _apply_custom_order(x_indices, x_labels, LABEL_ORDER)
    y_indices = _apply_custom_order(y_indices, y_labels, LABEL_ORDER)

    _print_subset(hist, x_indices, y_indices, x_labels, y_labels)
    if args.output_pdf:
        subset_hist = _build_subset_hist(hist, x_indices, y_indices, x_labels, y_labels)
        _save_subset_pdf(subset_hist, args.output_pdf)
        print(f"\nSaved PDF: {args.output_pdf}")


if __name__ == "__main__":
    main()
