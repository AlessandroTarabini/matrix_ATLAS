#!/usr/bin/env python3
"""
Toy study: statistical correlation between the fitted H->gamgam signal yields
in two overlapping kinematic bins (A, B), e.g. a pT(H) bin and a |y(H)| bin.

GOAL
----
In differential H->gamgam measurements the signal yields extracted in bins of
different observables (pT, rapidity, ...) are statistically correlated,
because the SAME events enter more than one observable's bin. The standard
way to estimate these correlations is the bootstrap: reweight the data with
Poisson(1) event weights many times and REDO THE FULL S+B FIT on every
replica. This script tests whether the same correlation can be obtained by
simply COUNTING events in the mass sidebands of each replica (no fits).

MODEL (per event)
-----------------
  * diphoton mass  m in [mlo, mhi] GeV
      - signal:     Gaussian(mu=125, sigma=2)
      - background: falling exponential  exp(-m/tau)
  * bin membership: each event belongs to bin A, bin B, both, or neither.
    The overlap probability is set to  p_AB = rho * sqrt(pA * pB),
    which for a Poisson process makes  corr(N_A, N_B) = rho  exactly
    (see membership_probs for the math). rho can be chosen independently
    for signal and background.

THREE ESTIMATES OF THE CORRELATION
---------------------------------
  1. "toys + fit"           : many INDEPENDENT pseudo-experiments, extended ML
                              S+B fit in each bin, correlate the fitted signal
                              yields. Only possible in a toy study --> ground truth.
  2. "bootstrap + refit"    : ONE dataset, Poisson(1) per-event bootstrap weights,
                              full refit of every replica (the expensive ATLAS way).
  3. "bootstrap + counting" : same bootstrap replicas, but only COUNT events in
                              the mass sidebands of each bin (no fits at all).
                              Correlate the sideband counts across replicas.

Usage examples:
  python3 run.py                      # defaults (S/B ~ H->gg like)
  python3 run.py --nsig 12000         # higher S/B
  python3 run.py --rho-sig 0.8 --rho-bkg 0.2
  python3 run.py --quick              # fast sanity check
"""

import argparse
import json
import math
import os
import sys
import time
from dataclasses import dataclass, asdict

import numpy as np

import ROOT

# no graphics windows, and silence RooFit's very chatty INFO/WARNING output
# (we run thousands of fits; the default printout would flood the terminal)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
@dataclass
class Config:
    # ---- mass model ------------------------------------------------------
    mlo: float = 100.0         # lower edge of the m_gg spectrum [GeV]
    mhi: float = 180.0         # upper edge of the m_gg spectrum [GeV]
    mu: float = 125.0          # signal peak position (Higgs mass) [GeV]
    sigma: float = 2.0         # signal Gaussian resolution [GeV]
    tau: float = 40.0          # background exponential decay constant [GeV]
    win_lo: float = 115.0      # signal window; everything outside
    win_hi: float = 135.0      #   [win_lo, win_hi] counts as "sideband"
    nbins: int = 320           # 0.25 GeV bins -> fine-binned ML fit ~ unbinned

    # ---- expected yields of the INCLUSIVE sample (full mass range) -------
    # The per-bin yields follow from these via the membership probabilities
    # below. nsig_tot / nbkg_tot is the main S/B knob (--nsig / --nbkg).
    nsig_tot: float = 300.0
    nbkg_tot: float = 100000.0

    # ---- bin membership --------------------------------------------------
    # pA / pB  = marginal probability for an event to fall in bin A / bin B.
    # A and B are bins of two DIFFERENT observables, so an event can be in
    # both (that overlap is what creates the statistical correlation) or in
    # neither (it falls in other bins of both observables). pA + pB does NOT
    # have to sum to 1 since we are taking about two different observables.
    # rho      = desired correlation of the event COUNTS in A and B,
    #            settable separately for signal and background.
    pA_sig: float = 0.40
    pB_sig: float = 0.50
    rho_sig: float = 0.50
    pA_bkg: float = 0.40
    pB_bkg: float = 0.50
    rho_bkg: float = 0.20

    # ---- statistics of the study ------------------------------------------
    ntoys: int = 10000          # independent toys for the ground truth
    nboot_fit: int = 10000       # bootstrap replicas that are refitted
    nboot_count: int = 10000   # bootstrap replicas for counting (cheap)
    seed: int = 20260715


def membership_probs(pA, pB, rho):
    """Translate (pA, pB, rho) into category probabilities (AB, A-only, B-only).

    Math behind it:
      The total number of events is Poisson(lambda). Splitting a Poisson
      total multinomially into exclusive categories gives INDEPENDENT Poisson
      counts per category (Poisson thinning theorem). With
          N_A = N_AB + N_Aonly     and     N_B = N_AB + N_Bonly,
      the only shared piece is N_AB, hence
          cov(N_A, N_B)  = Var(N_AB)      = lambda * p_AB
          Var(N_A)       = lambda * pA,   Var(N_B) = lambda * pB
          corr(N_A, N_B) = p_AB / sqrt(pA * pB).
      So choosing  p_AB = rho * sqrt(pA * pB)  makes the count correlation
      exactly rho.

    The check below ensures the four category probabilities
    (AB, A-only, B-only, neither) form a valid distribution; in particular
    p_AB <= min(pA, pB) limits the reachable rho to sqrt(min/max of pA,pB).
    """
    pAB = rho * math.sqrt(pA * pB)
    probs = (pAB, pA - pAB, pB - pAB)
    if min(probs) < 0 or sum(probs) > 1:
        raise ValueError(f"invalid membership probabilities for pA={pA}, pB={pB}, rho={rho}")
    return probs


def fit_start_yields(cfg):
    """Rough starting values for the S+B fit in each kinematic bin.

    Derived from the membership probabilities in Config (pA_sig, pB_sig, ...),
    so changing those knobs automatically updates every fit's starting point.
    """
    return {
        "nsig": {"A": cfg.nsig_tot * cfg.pA_sig, "B": cfg.nsig_tot * cfg.pB_sig},
        "nbkg": {"A": cfg.nbkg_tot * cfg.pA_bkg, "B": cfg.nbkg_tot * cfg.pB_bkg},
    }


def format_suffix_number(x):
    """Convert a number into a filename-safe token.

    Examples:
      10000 -> "10000"
      0.5   -> "0p5"
      -0.2  -> "m0p2"
    """
    if float(x).is_integer():
        s = str(int(round(float(x))))
    else:
        s = f"{float(x):g}"
    return s.replace("-", "m").replace(".", "p")


def output_suffix(cfg):
    """Common suffix for output files based on the main physics knobs."""
    parts = [
        format_suffix_number(cfg.nsig_tot),
        format_suffix_number(cfg.nbkg_tot),
        format_suffix_number(cfg.rho_sig),
        format_suffix_number(cfg.rho_bkg),
    ]
    return "_" + "_".join(parts)


# --------------------------------------------------------------------------
# event generation
# --------------------------------------------------------------------------
def sample_signal_mass(rng, n, cfg):
    """Draw n signal masses from Gaussian(mu, sigma), truncated to [mlo, mhi].

    Truncation by redrawing: values outside the range are simply drawn again
    until all are inside (with mu +- 4 sigma well inside the range this loop
    almost never runs more than once).
    """
    m = rng.normal(cfg.mu, cfg.sigma, size=n)
    bad = (m < cfg.mlo) | (m > cfg.mhi)
    while np.any(bad):
        m[bad] = rng.normal(cfg.mu, cfg.sigma, size=int(bad.sum()))
        bad = (m < cfg.mlo) | (m > cfg.mhi)
    return m


def sample_bkg_mass(rng, n, cfg):
    """Draw n background masses from exp(-m/tau) truncated to [mlo, mhi].

    Uses inverse-CDF sampling: for a truncated exponential the CDF can be
    inverted in closed form, so one uniform random number u in [0,1) maps
    directly to a mass. 'span' is the total probability mass of the
    untruncated exponential inside the window, needed for the normalisation.
    """
    span = 1.0 - math.exp(-(cfg.mhi - cfg.mlo) / cfg.tau)
    u = rng.random(n)
    return cfg.mlo - cfg.tau * np.log(1.0 - u * span)


def sample_component(rng, n_expected, probs3, mass_sampler, cfg):
    """Generate one component (signal OR background) of one pseudo-dataset.

    Steps:
      1. Fluctuate the total number of events: n ~ Poisson(n_expected).
      2. Draw a mass for every event from the component's mass pdf.
      3. Assign every event to one of 4 exclusive membership categories
         (both bins / A only / B only / neither) with the probabilities from
         membership_probs. Technically: one uniform number per event is
         placed in the cumulative-probability intervals via searchsorted.
      4. Convert the category index into two boolean flags inA / inB.

    Note that mass and bin membership are drawn independently: the mass
    shape is the same in bin A and bin B (a simplification; in reality the
    background slope varies a bit from bin to bin).
    """
    n = rng.poisson(n_expected)
    m = mass_sampler(rng, n, cfg)
    u = rng.random(n)
    cat = np.searchsorted(np.cumsum(probs3), u)   # 0=AB 1=A-only 2=B-only 3=none
    inA = (cat == 0) | (cat == 1)
    inB = (cat == 0) | (cat == 2)
    return m, inA, inB


def generate_dataset(rng, cfg, probs_sig, probs_bkg):
    """One full pseudo-dataset = signal component + background component.

    Returns three parallel arrays over ALL events:
      masses, inA flags, inB flags.
    The signal/background label is intentionally NOT returned: just like in
    real data, the analysis downstream must not know which event is which.
    """
    ms, sA, sB = sample_component(rng, cfg.nsig_tot, probs_sig, sample_signal_mass, cfg)
    mb, bA, bB = sample_component(rng, cfg.nbkg_tot, probs_bkg, sample_bkg_mass, cfg)
    return (np.concatenate([ms, mb]),
            np.concatenate([sA, bA]),
            np.concatenate([sB, bB]))


# --------------------------------------------------------------------------
# RooFit: extended ML fit of one mass spectrum
#
# Fit model per kinematic bin (like a H->gg analysis category):
#     n_sig * Gaussian(m; mu, sigma)  +  n_bkg * Exponential(m; c)
# Free parameters: n_sig, n_bkg, background slope c.
# Fixed parameters: mu, sigma (in the real analysis the signal shape comes
# from simulation, so it is not floated here either).
#
# The fit is binned in 0.25 GeV bins, which is statistically equivalent to
# an unbinned fit (the bin width is << sigma) but MUCH faster - important
# because this script performs thousands of fits.
# --------------------------------------------------------------------------
class BinFitter:
    def __init__(self, name, cfg):
        self.cfg = cfg
        # observable and pdfs are built ONCE per bin and reused for every fit
        self.m = ROOT.RooRealVar(f"m_{name}", "m_{#gamma#gamma}", cfg.mlo, cfg.mhi)
        # signal shape: fixed-parameter Gaussian (constants, not floated)
        self.mu = ROOT.RooRealVar(f"mu_{name}", "mu", cfg.mu)
        self.sg = ROOT.RooRealVar(f"sigma_{name}", "sigma", cfg.sigma)
        self.gauss = ROOT.RooGaussian(f"sig_{name}", "sig", self.m, self.mu, self.sg)
        # background shape: exponential exp(c*m) with free slope c (< 0)
        self.c = ROOT.RooRealVar(f"c_{name}", "c", -1.0 / cfg.tau, -0.2, 0.0)
        self.expo = ROOT.RooExponential(f"bkg_{name}", "bkg", self.m, self.c)
        # yields; n_sig may go negative so that downward background
        # fluctuations under the peak do not bias the fit at low S/B
        self.nsig = ROOT.RooRealVar(f"nsig_{name}", "nsig", 100.0, -2.0e4, 1.0e6)
        self.nbkg = ROOT.RooRealVar(f"nbkg_{name}", "nbkg", 1.0e4, 0.0, 1.0e7)
        # RooAddPdf with yield coefficients = extended likelihood model
        self.model = ROOT.RooAddPdf(f"model_{name}", "model",
                                    [self.gauss, self.expo], [self.nsig, self.nbkg])
        # one reusable TH1 as container for the (possibly weighted) spectrum
        self.hist = ROOT.TH1D(f"h_{name}", "", cfg.nbins, cfg.mlo, cfg.mhi)
        self.hist.SetDirectory(ROOT.nullptr)
        self.edges = np.linspace(cfg.mlo, cfg.mhi, cfg.nbins + 1)

    def fit_counts(self, counts, nsig0, nbkg0):
        """Fit a spectrum given as an array of bin contents.

        Returns (fitted n_sig, its HESSE error, fitted slope c, fit status);
        status == 0 means the minimisation converged.
        """
        # 1. copy the numpy bin contents into the reusable TH1
        for i, cval in enumerate(counts):
            self.hist.SetBinContent(i + 1, float(cval))
        # 2. wrap it into a RooDataHist (cheap: just 'nbins' numbers)
        dh = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(self.m), self.hist)
        # 3. reset the floating parameters to sensible starting values -
        #    important because the same RooFit objects are reused for
        #    thousands of fits and would otherwise start from the previous
        #    fit's minimum
        self.nsig.setVal(nsig0)
        self.nbkg.setVal(nbkg0)
        self.c.setVal(-1.0 / self.cfg.tau)
        # 4. run the extended binned ML fit (Minuit2 migrad)
        res = self.model.fitTo(dh, Save=True, PrintLevel=-1, Extended=True,
                               Minimizer=("Minuit2", "migrad"), Strategy=0,
                               Offset=True, Warnings=False, Verbose=False)
        status = res.status()
        out = (self.nsig.getVal(), self.nsig.getError(), self.c.getVal(), status)
        del res, dh
        return out

    def fit_masses(self, masses, nsig0, nbkg0, weights=None):
        """Fit a spectrum given as raw event masses (optionally weighted).

        The optional per-event weights are how bootstrap replicas enter:
        a replica is the SAME events histogrammed with Poisson(1) weights.
        """
        counts, _ = np.histogram(masses, bins=self.edges, weights=weights)
        return self.fit_counts(counts, nsig0, nbkg0)


# --------------------------------------------------------------------------
# counting helpers (the cheap alternative to refitting)
# --------------------------------------------------------------------------
def cell_counts(masses, inA, inB, cfg):
    """Split the dataset into 6 DISJOINT cells and count events in each.

    The cells are (bin membership) x (mass region):
        membership: AB = in both bins, Ao = A only, Bo = B only
        mass:       win = inside [win_lo, win_hi], sb = sidebands
    Events in neither bin are irrelevant and simply dropped.

    Everything the counting methods need later can be assembled from these
    6 numbers, e.g. sideband count of bin A = AB_sb + Ao_sb.
    """
    win = (masses >= cfg.win_lo) & (masses < cfg.win_hi)
    cells = {}
    for tag, mask in (("AB", inA & inB), ("Ao", inA & ~inB), ("Bo", ~inA & inB)):
        cells[f"{tag}_win"] = int(np.count_nonzero(mask & win))
        cells[f"{tag}_sb"] = int(np.count_nonzero(mask & ~win))
    return cells


def sideband_counts_from_weights(weights, masks):
    """Count sideband events in bins A and B for ONE bootstrap replica.

    This is the literal bootstrap-counting procedure:
      1. take one Poisson(1) weight per event,
      2. sum those weights over the sideband events of bin A and of bin B.

    The masks are disjoint sideband categories:
      AB_sb = event is in both bins and in the sidebands
      Ao_sb = event is only in A and in the sidebands
      Bo_sb = event is only in B and in the sidebands
    """
    CA_sb = float(weights[masks["AB_sb"]].sum() + weights[masks["Ao_sb"]].sum())
    CB_sb = float(weights[masks["AB_sb"]].sum() + weights[masks["Bo_sb"]].sum())
    return CA_sb, CB_sb


def corr(x, y):
    """Pearson correlation coefficient of two samples."""
    return float(np.corrcoef(np.asarray(x), np.asarray(y))[0, 1])


def fisher_err(rho, n):
    """Statistical uncertainty of a correlation estimated from n samples
    (from the Fisher z-transformation): sigma_rho ~ (1-rho^2)/sqrt(n-3)."""
    return (1.0 - rho * rho) / math.sqrt(max(n - 3, 1))


# --------------------------------------------------------------------------
# analytic expectation for the sideband-count correlation
# --------------------------------------------------------------------------
def predicted_sideband_corr(cfg, probs_sig, probs_bkg):
    """Closed-form prediction of corr(C_A, C_B) for the raw sideband counts.

    Same formula as in membership_probs (corr = shared / sqrt(tot_A*tot_B)),
    but restricted to sideband events: each expected count is
        N = nsig * P(membership) * f_sb(signal) + nbkg * P(membership) * f_sb(bkg)
    where f_sb is the fraction of the component's mass pdf falling in the
    sidebands (tiny for the signal, ~70% for the background). Because the
    sidebands are background-dominated, the result is essentially rho_bkg.
    """
    # sideband fraction of the (truncated) exponential background pdf
    span = lambda lo, hi: math.exp(-lo / cfg.tau) - math.exp(-hi / cfg.tau)
    f_sb_bkg = 1.0 - span(cfg.win_lo, cfg.win_hi) / span(cfg.mlo, cfg.mhi)
    # sideband fraction of the Gaussian signal pdf (leakage outside +-4 sigma)
    from math import erf, sqrt
    gcdf = lambda x: 0.5 * (1 + erf((x - cfg.mu) / (cfg.sigma * sqrt(2))))
    f_sb_sig = 1.0 - (gcdf(cfg.win_hi) - gcdf(cfg.win_lo))

    def sb_yield(idx_probs):
        i_sig, i_bkg = idx_probs
        return cfg.nsig_tot * i_sig * f_sb_sig + cfg.nbkg_tot * i_bkg * f_sb_bkg

    # expected sideband counts: shared (AB), total in A, total in B
    pAB = sb_yield((probs_sig[0], probs_bkg[0]))
    pA = sb_yield((probs_sig[0] + probs_sig[1], probs_bkg[0] + probs_bkg[1]))
    pB = sb_yield((probs_sig[0] + probs_sig[2], probs_bkg[0] + probs_bkg[2]))
    return pAB / math.sqrt(pA * pB)


# --------------------------------------------------------------------------
# the two big procedures
# --------------------------------------------------------------------------
def run_toys_and_fit(rng, cfg, probs_sig, probs_bkg, fitters, log=print,
                     plot_every=0, outdir=None):
    """GROUND TRUTH: many independent pseudo-experiments.

    For each toy:
      1. generate a completely new dataset from the true model,
      2. fit the mass spectrum of the events in bin A -> S_A,
      3. fit the mass spectrum of the events in bin B -> S_B.
    The correlation of (S_A, S_B) across toys is the TRUE sampling
    correlation of the estimator. This is only possible because we know the
    true model - on real data one would have to use the bootstrap instead.

    If plot_every > 0, the fitted spectra of every plot_every-th toy are
    saved as PDFs in <outdir>/toy_fits/ (useful to eyeball fit quality).
    """
    fitA, fitB = fitters
    # starting values for the fits = expected yields per bin from Config
    y0 = fit_start_yields(cfg)
    nsig0, nbkg0 = y0["nsig"], y0["nbkg"]
    plotdir = None
    if plot_every > 0:
        plotdir = os.path.join(outdir or ".", "toy_fits")
        os.makedirs(plotdir, exist_ok=True)
    suffix = output_suffix(cfg)
    SA, SB, errA, errB = [], [], [], []
    nfail = 0
    t0 = time.time()
    for itoy in range(cfg.ntoys):
        m, inA, inB = generate_dataset(rng, cfg, probs_sig, probs_bkg)
        sA, eA, cA, stA = fitA.fit_masses(m[inA], nsig0["A"], nbkg0["A"])
        sB, eB, cB, stB = fitB.fit_masses(m[inB], nsig0["B"], nbkg0["B"])
        if stA != 0 or stB != 0:      # drop toys where a fit failed
            nfail += 1
            continue
        SA.append(sA); SB.append(sB); errA.append(eA); errB.append(eB)
        if plotdir is not None and (itoy + 1) % plot_every == 0:
            # save this toy's fitted spectra; n_bkg is reconstructed as
            # (total events in the bin) - (fitted n_sig), which is what the
            # extended fit returns up to rounding
            fitres = [(sA, inA.sum() - sA, cA), (sB, inB.sum() - sB, cB)]
            path = plot_example_fit(cfg, (m, inA, inB), fitters, fitres, plotdir,
                                    filename=f"toy_{itoy + 1:04d}{suffix}.pdf",
                                    tag=f"toy {itoy + 1}, ")
            log(f"    saved fit plot: {path}")
        if (itoy + 1) % 200 == 0:
            log(f"    toy {itoy + 1}/{cfg.ntoys}  ({time.time() - t0:.0f}s)")
    return (np.array(SA), np.array(SB), np.array(errA), np.array(errB), nfail)


def run_bootstrap_on_dataset(rng, cfg, dataset, fitters, log=print,
                             plot_every=0, outdir=None):
    """Bootstrap of ONE observed dataset - refit (expensive) vs count (cheap).

    This is what one can actually do on real data: the dataset is fixed,
    and replicas are made by reweighting its events with independent
    Poisson(1) weights (mean 1, so on average the dataset reproduces itself,
    but each replica fluctuates like an independent Poisson sample).

    Two analyses are run on the SAME bootstrap replicas:
      (a) counting (cheap)  - count sideband events only, no fits;
      (b) refit (expensive) - histogram the weighted events of bin A and of
          bin B, and refit both spectra.

    The crucial point is that one Poisson(1) weight vector w is generated per
    replica and reused for BOTH methods. This means the sideband counts and
    the refits see exactly the same fluctuated dataset replica, removing any
    bias from comparing different bootstrap samples.

    If plot_every > 0, the fitted spectra of every plot_every-th accepted
    replica are saved in <outdir>/boot_fits/.
    """
    m, inA, inB = dataset
    fitA, fitB = fitters

    # masks for the literal sideband counting on each replica
    win = (m >= cfg.win_lo) & (m < cfg.win_hi)
    sb_masks = {
        "AB_sb": inA & inB & ~win,
        "Ao_sb": inA & ~inB & ~win,
        "Bo_sb": ~inA & inB & ~win,
    }

    # ---- bootstrap loop: SAME replica for counting and refit ---------------
    mA, mB = m[inA], m[inB]                       # masses per bin (fixed)
    idxA, idxB = np.where(inA)[0], np.where(inB)[0]  # to slice the weights
    # starting values from Config membership probabilities (pA_sig, pB_sig, ...)
    y0 = fit_start_yields(cfg)
    nsig0 = y0["nsig"]
    # for n_bkg, the observed event count in the bin is a better guess than
    # the model expectation, because this is the one dataset we actually have
    nbkg0 = {"A": float(inA.sum()), "B": float(inB.sum())}
    plotdir = None
    if plot_every > 0:
        plotdir = os.path.join(outdir or ".", "boot_fits")
        os.makedirs(plotdir, exist_ok=True)
    suffix = output_suffix(cfg)
    SA, SB, CA, CB = [], [], [], []
    nfail = 0
    t0 = time.time()
    for irep in range(cfg.nboot_fit):
        # one Poisson(1) weight per event of the FULL dataset
        w = rng.poisson(1.0, size=m.size).astype(float)
        # count sideband events in this exact same replica
        countA, countB = sideband_counts_from_weights(w, sb_masks)
        # ... and refit the two mass spectra built from the same weights
        sA, _, slopeA, stA = fitA.fit_masses(mA, nsig0["A"], nbkg0["A"], weights=w[idxA])
        sB, _, slopeB, stB = fitB.fit_masses(mB, nsig0["B"], nbkg0["B"], weights=w[idxB])
        if stA != 0 or stB != 0:
            nfail += 1
            continue
        # keep both methods on the exact same accepted replica subset
        SA.append(sA); SB.append(sB); CA.append(countA); CB.append(countB)
        if plotdir is not None and (irep + 1) % plot_every == 0:
            nA, nB = float(w[idxA].sum()), float(w[idxB].sum())
            fitres = [(sA, nA - sA, slopeA), (sB, nB - sB, slopeB)]
            path = plot_example_fit(cfg, (m, inA, inB), fitters, fitres, plotdir,
                                    filename=f"replica_{irep + 1:05d}{suffix}.pdf",
                                    tag=f"replica {irep + 1}, ",
                                    weights=w)
            log(f"    saved bootstrap fit plot: {path}")
        if (irep + 1) % 100 == 0:
            log(f"    replica {irep + 1}/{cfg.nboot_fit}  ({time.time() - t0:.0f}s)")
    out = {"rho_sideband": corr(CA, CB), "nrep_used": len(SA)}
    out["rho_bootfit"] = corr(SA, SB)
    out["std_SA_bootfit"] = float(np.std(SA))
    out["std_SB_bootfit"] = float(np.std(SB))
    out["nfail"] = nfail
    return out


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
def plot_example_fit(cfg, dataset, fitters, results, outdir,
                     filename="example_fit.pdf", tag="", weights=None):
    """Draw one dataset in each bin with its fitted S+B model overlaid
    (data points, total fit, background component, window edges).

    Used for toy fits, bootstrap refit replicas, etc. If weights is given,
    the displayed histogram uses those per-event bootstrap weights."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m, inA, inB = dataset
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=False)
    grid = np.linspace(cfg.mlo, cfg.mhi, 500)

    for ax, mask, label, res in zip(axes, (inA, inB), ("bin A", "bin B"), results):
        nsig, nbkg, cslope = res
        # coarse binning just for display (the fit itself used cfg.nbins)
        edges = np.linspace(cfg.mlo, cfg.mhi, 56)
        bw = edges[1] - edges[0]
        cnt, _ = np.histogram(m[mask], bins=edges,
                              weights=None if weights is None else weights[mask])
        ctr = 0.5 * (edges[:-1] + edges[1:])
        ax.errorbar(ctr, cnt, yerr=np.sqrt(cnt), fmt="ko", ms=3, lw=1, label="toy data")
        # rebuild the fitted curves analytically from the fitted parameters:
        # normalised Gaussian pdf and normalised (truncated) exponential pdf,
        # scaled by the fitted yields and the display bin width
        gpdf = np.exp(-0.5 * ((grid - cfg.mu) / cfg.sigma) ** 2) / (cfg.sigma * math.sqrt(2 * math.pi))
        enorm = cslope / (math.exp(cslope * cfg.mhi) - math.exp(cslope * cfg.mlo))
        epdf = enorm * np.exp(cslope * grid)
        ax.plot(grid, bw * (nsig * gpdf + nbkg * epdf), "r-", lw=1.5, label="S+B fit")
        ax.plot(grid, bw * nbkg * epdf, "b--", lw=1.2, label="B component")
        for x in (cfg.win_lo, cfg.win_hi):     # mark the window/sideband split
            ax.axvline(x, color="gray", ls=":", lw=1)
        ax.set_xlabel(r"$m_{\gamma\gamma}$ [GeV]")
        ax.set_ylabel(f"events / {bw:.1f} GeV")
        ax.set_title(f"{tag}{label}:  $N_S$={nsig:.0f}, $N_B$={nbkg:.0f}")
        ax.legend()
    fig.tight_layout()
    path = os.path.join(outdir, filename)
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_summary(cfg, toys, boot, boot_errs, rho_truth, rho_truth_err, outdir):
    """Single-panel summary of the correlation estimates."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8.2, 4.8))

    # one point per method, compared to the toy truth band and to the
    # input signal/background count correlations configured in the model
    methods = [
        ("toys + fit\n(truth)", rho_truth, rho_truth_err, "tab:blue"),
        ("bootstrap\n+ refit", None, None, "tab:red"),
        ("bootstrap +\nsideband count", None, None, "tab:green"),
    ]
    keys = [None, "rho_bootfit", "rho_sideband"]
    xs, ys, es, cols, labels = [], [], [], [], []
    for i, ((lab, v, e, col), key) in enumerate(zip(methods, keys)):
        if key is not None:
            v, e = boot[key], boot_errs[key]
        xs.append(i); ys.append(v); es.append(e); cols.append(col); labels.append(lab)
    for x, y, e, col in zip(xs, ys, es, cols):
        ax.errorbar([x], [y], yerr=[e], fmt="o", color=col, capsize=4, ms=7)
    ax.axhline(rho_truth, color="tab:blue", ls="--", lw=1, alpha=0.7)
    ax.axhspan(rho_truth - rho_truth_err, rho_truth + rho_truth_err, color="tab:blue", alpha=0.12)
    ax.axhline(cfg.rho_bkg, color="tab:green", ls=":", lw=1.2, alpha=0.9,
               label=rf"input $\rho_{{bkg}}$ = {cfg.rho_bkg:.3f}")
    ax.axhline(cfg.rho_sig, color="tab:purple", ls=":", lw=1.2, alpha=0.9,
               label=rf"input $\rho_{{sig}}$ = {cfg.rho_sig:.3f}")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(r"correlation")
    ax.set_title(f"$\\rho_{{sig}}$={cfg.rho_sig}, $\\rho_{{bkg}}$={cfg.rho_bkg}, "
                 f"$N_S$={cfg.nsig_tot:.0f}, $N_B$={cfg.nbkg_tot:.0f}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    path = os.path.join(outdir, f"correlation_summary{output_suffix(cfg)}.pdf")
    fig.savefig(path)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nsig", type=float, default=Config.nsig_tot,
                   help="expected inclusive signal yield (controls S/B)")
    p.add_argument("--nbkg", type=float, default=Config.nbkg_tot,
                   help="expected inclusive background yield")
    p.add_argument("--rho-sig", type=float, default=Config.rho_sig,
                   help="count correlation between bins A,B for signal")
    p.add_argument("--rho-bkg", type=float, default=Config.rho_bkg,
                   help="count correlation between bins A,B for background")
    p.add_argument("--ntoys", type=int, default=Config.ntoys)
    p.add_argument("--nboot-fit", type=int, default=Config.nboot_fit)
    p.add_argument("--nboot-count", type=int, default=Config.nboot_count)
    p.add_argument("--seed", type=int, default=Config.seed)
    p.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    p.add_argument("--quick", action="store_true", help="small statistics for testing")
    p.add_argument("--plot-toys", type=int, default=0, metavar="N",
                   help="save the fitted spectra of every N-th toy as a PDF "
                        "in <outdir>/toy_fits/ (0 = off, e.g. 200)")
    p.add_argument("--plot-boot", type=int, default=0, metavar="N",
                   help="save the fitted spectra of every N-th bootstrap "
                        "replica as a PDF in <outdir>/boot_fits/ (0 = off)")
    return p.parse_args()


def main():
    # ---- setup -------------------------------------------------------------
    args = parse_args()
    cfg = Config(nsig_tot=args.nsig, nbkg_tot=args.nbkg,
                 rho_sig=args.rho_sig, rho_bkg=args.rho_bkg,
                 ntoys=args.ntoys,
                 nboot_fit=args.nboot_fit, nboot_count=args.nboot_count,
                 seed=args.seed)
    if args.quick:
        cfg.ntoys, cfg.nboot_fit, cfg.nboot_count = 100, 50, 5000
    # The comparison now uses the SAME bootstrap replicas for refit and
    # sideband counting, so the effective number of counting replicas is the
    # number of refit replicas.
    cfg.nboot_count = cfg.nboot_fit

    rng = np.random.default_rng(cfg.seed)
    # translate the (pA, pB, rho) knobs into category probabilities,
    # separately for signal and background
    probs_sig = membership_probs(cfg.pA_sig, cfg.pB_sig, cfg.rho_sig)
    probs_bkg = membership_probs(cfg.pA_bkg, cfg.pB_bkg, cfg.rho_bkg)

    # one fitter (= one RooFit model) per kinematic bin, reused for all fits
    fitters = (BinFitter("A", cfg), BinFitter("B", cfg))

    # ---- [1] ground truth: independent toys + fit ---------------------------
    print(f"[1/3] ground truth: {cfg.ntoys} independent toys, S+B fit per bin")
    SA, SB, errA, errB, nfail = run_toys_and_fit(rng, cfg, probs_sig, probs_bkg, fitters,
                                                 plot_every=args.plot_toys,
                                                 outdir=args.outdir)
    rho_truth = corr(SA, SB)
    rho_truth_err = fisher_err(rho_truth, len(SA))
    print(f"  rho(fit) = {rho_truth:.4f} +- {rho_truth_err:.4f}   "
          f"({nfail} failed toys)")

    # ---- [2] bootstrap of one observed dataset ------------------------------
    print(f"[2/3] bootstrap of ONE observed dataset "
          f"({cfg.nboot_fit} shared replicas for refit and counting)")
    # this dataset plays the role of "the data we recorded" - everything
    # from here on could be done identically on real data
    dataset = generate_dataset(rng, cfg, probs_sig, probs_bkg)
    boot = run_bootstrap_on_dataset(rng, cfg, dataset, fitters,
                                    plot_every=args.plot_boot,
                                    outdir=args.outdir)
    # replica-statistics uncertainty on each bootstrap correlation estimate
    # (note: this does NOT include the dataset-to-dataset fluctuation of the
    # bootstrap central value - rerun with another --seed to probe that)
    boot_errs = {"rho_bootfit": fisher_err(boot["rho_bootfit"], boot["nrep_used"]),
                 "rho_sideband": fisher_err(boot["rho_sideband"], boot["nrep_used"])}
    print(f"  rho(bootstrap refit) = {boot['rho_bootfit']:.4f}, "
          f"rho(sideband count) = {boot['rho_sideband']:.4f}")

    # ---- [3] summary table, plots, json --------------------------------------
    print("\n[3/3] summary")
    rows = [
        ("independent toys + fit (truth)", rho_truth, rho_truth_err),
        ("bootstrap + refit (expensive)", boot["rho_bootfit"], boot_errs["rho_bootfit"]),
        ("bootstrap + sideband counting", boot["rho_sideband"], boot_errs["rho_sideband"]),
    ]
    w = max(len(r[0]) for r in rows)
    print(f"  {'method':<{w}}   correlation")
    for name, v, e in rows:
        print(f"  {name:<{w}}   {v:.4f} +- {e:.4f}")

    # cross-check of the yield UNCERTAINTY (not correlation): the bootstrap
    # refit spread should match the toy spread; the counting estimator is
    # expected to be wider (it is a less efficient estimator than the fit)
    print(f"\n  yield spreads bin A: toys+fit = {SA.std():.1f}, "
          f"bootstrap refit = {boot['std_SA_bootfit']:.1f}")

    p2 = plot_summary(cfg, (SA, SB), boot, boot_errs, rho_truth, rho_truth_err,
                      args.outdir)
    print(f"\n  plot: {p2}")

    out = {
        "config": asdict(cfg),
        "rho_truth": rho_truth, "rho_truth_err": rho_truth_err,
        "bootstrap": boot,
        "bootstrap_errs": boot_errs,
    }
    jpath = os.path.join(args.outdir, f"results{output_suffix(cfg)}.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  results: {jpath}")


if __name__ == "__main__":
    main()
