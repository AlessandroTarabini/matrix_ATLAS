#!/usr/bin/env python3
#python bootstrap.py
#If your mass column is not called mass (e.g. mgg), run: python bootstrap.py --mass-column mgg
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


DEFAULT_INPUT_DIR = (
    "/eos/cms/store/group/phys_higgs/cmshgg/earlyRun3Hgg/"
    "analysis_Ntuples/forPaper_24_11_11"
)

BINS_PTH = [0, 15, 30, 45, 80, 120, 200, 350, 10000]
BINS_NJ = [0, 1, 2, 3, 1000]
BINS_PTJ0 = [0, 30, 75, 120, 200, 10000]

PTH_COLUMN = "pt"
NJ_COLUMN = "NJ"
PTJ0_COLUMN = "PTJ0"


def edge_bin_mask(x: np.ndarray, edges: np.ndarray, b: int) -> np.ndarray:
    """Bin b: [edges[b], edges[b+1]) except last bin [edges[-2], edges[-1]] closed on the right."""
    lo = float(edges[b])
    hi = float(edges[b + 1])
    if b == len(edges) - 2:
        return (x >= lo) & (x <= hi)
    return (x >= lo) & (x < hi)


def write_binned_counts_txt(
    path: Path,
    header_comment: str,
    edges: np.ndarray,
    counts: np.ndarray,
) -> None:
    """counts shape (n_bins, n_replicas)."""
    n_bins, n_rep = counts.shape
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(header_comment.rstrip() + "\n")
        f.write("# edges " + " ".join(str(float(x)) for x in edges) + "\n")
        f.write(f"# n_bins {n_bins}\n")
        f.write("# replica_idx " + " ".join(f"bin{b}" for b in range(n_bins)) + "\n")
        for r in range(n_rep):
            f.write(
                str(r)
                + " "
                + " ".join(str(int(counts[b, r])) for b in range(n_bins))
                + "\n"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read parquet files recursively, keep sidebands, and assign "
            "N Poisson(1) bootstrap weights per event; each event uses one "
            "seed derived from run/lumi/event."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(DEFAULT_INPUT_DIR),
        help="Base directory to scan recursively for parquet files.",
    )
    parser.add_argument(
        "--mass-column",
        default="mass",
        help="Column name containing diphoton mass.",
    )
    parser.add_argument(
        "--low-cut",
        type=float,
        default=115.0,
        help="Lower sideband cut (keep events below this value).",
    )
    parser.add_argument(
        "--high-cut",
        type=float,
        default=130.0,
        help="Upper sideband cut (keep events above this value).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("bootstrap_sideband_counts.txt"),
        help="Output txt file with sideband event count per bootstrap replica.",
    )
    parser.add_argument(
        "--run-column",
        default="run",
        help="Parquet column for run number.",
    )
    parser.add_argument(
        "--lumi-column",
        default="lumi",
        help="Parquet column for luminosity block (LS).",
    )
    parser.add_argument(
        "--event-column",
        default="event",
        help="Parquet column for event number.",
    )
    parser.add_argument(
        "--seed-offset",
        type=int,
        default=0,
        help=(
            "Optional extra entropy XORed into the per-event seed (e.g. to "
            "define independent bootstrap replicates while keeping run/lumi/event)."
        ),
    )
    parser.add_argument(
        "--n-replicas",
        type=int,
        default=100,
        help=(
            "Number of independent Poisson(1) bootstrap weights per event "
            "(used to build sideband event count for each replica)."
        ),
    )
    return parser.parse_args()


def collect_parquet_files(base_dir: Path) -> list[Path]:
    return sorted(base_dir.rglob("*.parquet"))


def seeds_from_run_lumi_event(
    run: np.ndarray,
    lumi: np.ndarray,
    event: np.ndarray,
    seed_offset: int,
) -> np.ndarray:
    """
    Deterministic 32-bit seed per row from (run, lumi, event).
    Does not use Python's hash() (not stable across processes).
    """
    r = np.asarray(run, dtype=np.uint64)
    l = np.asarray(lumi, dtype=np.uint64)
    e = np.asarray(event, dtype=np.uint64)
    # Mix into a single word (FNV-1a style constants)
    x = np.uint64(14695981039346656037)
    x = (x ^ r) * np.uint64(1099511628211)
    x = (x ^ l) * np.uint64(1099511628211)
    x = (x ^ e) * np.uint64(1099511628211)
    x = x ^ np.uint64(seed_offset & 0xFFFFFFFFFFFFFFFF)
    return (x & np.uint64(0xFFFFFFFF)).astype(np.uint64)


def poisson1_n_per_seed(seeds: np.ndarray, n_replicas: int) -> np.ndarray:
    """Draw N Poisson(1) values from each event-specific RNG seed."""
    out = np.empty((len(seeds), n_replicas), dtype=np.int64)
    for i, s in enumerate(seeds):
        out[i, :] = np.random.default_rng(int(s)).poisson(1.0, size=n_replicas)
    return out


def main() -> None:
    args = parse_args()
    if args.n_replicas < 1:
        raise ValueError("--n-replicas must be >= 1")
    parquet_files = collect_parquet_files(args.input_dir)

    if not parquet_files:
        raise FileNotFoundError(f"No parquet files found under {args.input_dir}")

    id_cols = (args.run_column, args.lumi_column, args.event_column)
    edges_pth = np.asarray(BINS_PTH, dtype=np.float64)
    n_bins_pth = len(edges_pth) - 1
    edges_nj = np.asarray(BINS_NJ, dtype=np.float64)
    n_bins_nj = len(edges_nj) - 1
    edges_ptj0 = np.asarray(BINS_PTJ0, dtype=np.float64)
    n_bins_ptj0 = len(edges_ptj0) - 1

    total_events = 0
    total_sideband_events = 0
    replica_counts = np.zeros(args.n_replicas, dtype=np.int64)
    replica_counts_pth = np.zeros((n_bins_pth, args.n_replicas), dtype=np.int64)
    replica_counts_nj = np.zeros((n_bins_nj, args.n_replicas), dtype=np.int64)
    replica_counts_ptj0 = np.zeros((n_bins_ptj0, args.n_replicas), dtype=np.int64)

    iterator = parquet_files
    if tqdm is not None:
        iterator = tqdm(parquet_files, desc="Processing parquet files", unit="file")
    else:
        print(f"Processing {len(parquet_files)} parquet files...")

    for i, path in enumerate(iterator, start=1):
        if tqdm is None:
            print(f"[{i}/{len(parquet_files)}] {path.name}")
        df = pd.read_parquet(path)
        #Check that all columns are there in the parquet file
        if args.mass_column not in df.columns:
            raise KeyError(
                f"Mass column '{args.mass_column}' not found in file: {path}"
            )
        for c in id_cols:
            if c not in df.columns:
                raise KeyError(
                    f"Column '{c}' not found in file: {path}"
                )
        for col, label in (
            (PTH_COLUMN, "PTH (pt)"),
            (NJ_COLUMN, "NJ"),
            (PTJ0_COLUMN, "PTJ0"),
        ):
            if col not in df.columns:
                raise KeyError(
                    f"Column '{col}' ({label}) not found in file: {path}"
                )
        total_events += len(df)
        mask = ((df[args.mass_column] < args.low_cut) | (df[args.mass_column] > args.high_cut)) & (df["lead_mvaID"] > 0.25) & (df["sublead_mvaID"] > 0.25)
        df_sb = df.loc[mask].copy()

        if df_sb.empty:
            continue

        n_sb = len(df_sb)
        r = df_sb[args.run_column].to_numpy()
        lum = df_sb[args.lumi_column].to_numpy()
        ev = df_sb[args.event_column].to_numpy()
        ## seeds is a 1D array that holds the seed for every sideband event
        seeds = seeds_from_run_lumi_event(r, lum, ev, args.seed_offset)
        ## weight_matrix is a 2D array that holds the bootstrap weights for every sideband event and every replica
        ## rows (n_sb) -> onw row per event that passed the sideband cut
        ## columns (args.n_replicas) -> one column per replica
        weight_matrix = poisson1_n_per_seed(seeds, args.n_replicas)
        replica_counts += weight_matrix.sum(axis=0)

        pt = df_sb[PTH_COLUMN].to_numpy(dtype=np.float64)
        finite_pt = np.isfinite(pt)
        for b in range(n_bins_pth):
            m = finite_pt & edge_bin_mask(pt, edges_pth, b)
            if not np.any(m):
                continue
            replica_counts_pth[b] += weight_matrix[m].sum(axis=0)

        nj = df_sb[NJ_COLUMN].to_numpy(dtype=np.float64)
        finite_nj = np.isfinite(nj)
        for b in range(n_bins_nj):
            m = finite_nj & edge_bin_mask(nj, edges_nj, b)
            if not np.any(m):
                continue
            replica_counts_nj[b] += weight_matrix[m].sum(axis=0)

        ptj0 = df_sb[PTJ0_COLUMN].to_numpy(dtype=np.float64)
        finite_ptj0 = np.isfinite(ptj0)
        for b in range(n_bins_ptj0):
            m = finite_ptj0 & edge_bin_mask(ptj0, edges_ptj0, b)
            if not np.any(m):
                continue
            replica_counts_ptj0[b] += weight_matrix[m].sum(axis=0)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        f.write("# inclusive: sum of Poisson(1) weights in sidebands per replica\n")
        f.write("# replica_index sideband_event_count\n")
        for idx, count in enumerate(replica_counts):
            f.write(f"{idx} {int(count)}\n")

    stem = args.output.stem
    parent = args.output.parent
    pth_path = parent / f"{stem}_pth.txt"
    nj_path = parent / f"{stem}_nj.txt"
    ptj0_path = parent / f"{stem}_ptj0.txt"
    write_binned_counts_txt(
        pth_path,
        "# PTH bins (column pt); sum of Poisson weights per bin per replica",
        edges_pth,
        replica_counts_pth,
    )
    write_binned_counts_txt(
        nj_path,
        "# NJ bins (column NJ); sum of Poisson weights per bin per replica",
        edges_nj,
        replica_counts_nj,
    )
    write_binned_counts_txt(
        ptj0_path,
        "# PTJ0 bins (column PTJ0); sum of Poisson weights per bin per replica",
        edges_ptj0,
        replica_counts_ptj0,
    )

    print(f"PTH-binned output     : {pth_path}")
    print(f"NJ-binned output      : {nj_path}")
    print(f"PTJ0-binned output    : {ptj0_path}")


if __name__ == "__main__":
    main()
