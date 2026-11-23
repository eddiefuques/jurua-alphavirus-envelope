#!/usr/bin/env python3
"""Figure 2 - E2 divergence relative to E1 in Juruá alphavirus.

Panel a: pairwise amino-acid p-distances, E2 against E1, for the nine complex taxa.
Panel b: externally studentized residuals for all 67 branches, ranked.
Panel c: the JAV residual across 200 bootstrap replicates.

Builds two versions from the same data. The dissertation version must be placed at 100 %
of the 6.5 in text block; scaling it down reduces the effective type size below its 10 pt
floor.

Usage  : python figures/figure2_divergence.py          (run from the repository root)
Inputs : data/tables/e2_e1_pairwise_distances_9taxa.csv
         data/tables/branch_lengths_35taxa.csv
         data/tables/bootstrap_jav_residual_200reps.csv
Outputs: Figure2_divergence_{paper,dissertation}.{png,pdf}
"""
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress, t as tdist

PANEL_A = "data/tables/e2_e1_pairwise_distances_9taxa.csv"
PANEL_B = "data/tables/branch_lengths_35taxa.csv"
PANEL_C = "data/tables/bootstrap_jav_residual_200reps.csv"
JAVC, OTHC, ALT, HIST = "#B8411F", "#8C9BA8", "#1F6F4A", "#A9C3D9"

D  = pd.read_csv(PANEL_A)
S  = pd.read_csv(PANEL_B).sort_values("t").reset_index(drop=True)
BT = pd.read_csv(PANEL_C)
bg, jv = D[~D.JAV], D[D.JAV]
fit = linregress(bg.E1, bg.E2)                      # background fit, JAV excluded
tcrit = tdist.isf(0.05 / len(S) / 2, len(S) - 3)    # Bonferroni over all branches
lo95, hi95 = np.percentile(BT.resid, [2.5, 97.5])
OBS = 0.210


def style(sizes):
    """Three-size role ladder: base / annotation / tick."""
    base, ann, tick = sizes
    mpl.rcParams.update({
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "legend.fontsize": ann, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.titlelocation": "left", "axes.titleweight": "normal",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 300,
        "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def build(path, figsize, nrows, ncols, sizes, pl_size, lw=1.0, ms=(26, 34), bs=16, legy=0.06):
    base, ann, _ = sizes
    style(sizes)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    ax = axes.ravel()

    # (a) pairwise p-distances - description only, not used for inference
    xs = np.linspace(D.E1.min() * .96, D.E1.max() * 1.04, 50)
    ax[0].plot(xs, fit.slope * xs + fit.intercept, color="0.35", lw=lw, zorder=2)
    ax[0].scatter(bg.E1, bg.E2, s=ms[0], color=OTHC, ec="white", lw=.5, zorder=3,
                  label="within-complex pairs")
    ax[0].scatter(jv.E1, jv.E2, s=ms[1], color=JAVC, ec="white", lw=.6, zorder=4,
                  label="JAV pairs")
    ax[0].set_xlabel("E1 p-distance"); ax[0].set_ylabel("E2 p-distance")
    ax[0].set_title("Pairwise divergence, complex", fontsize=base)
    ax[0].legend(loc="lower right", bbox_to_anchor=(1.0, legy), handletextpad=.4,
                 borderpad=.3, labelspacing=.55)
    ax[0].text(.04, .96, f"$R^2$ = {fit.rvalue ** 2:.2f}", transform=ax[0].transAxes,
               va="top", fontsize=ann, color="0.3")
    ax[0].margins(0.07)

    # (b) leave-one-out studentized residuals, ranked; identity by key, not by
    #     labels on the points - the two significant branches sit at the right
    #     edge where inline labels collide with the title and threshold text
    col = [JAVC if n == "Jurua_alphavirus" else (ALT if n == "Barmah_Forest_virus" else OTHC)
           for n in S.taxa]
    ax[1].axhspan(-tcrit, tcrit, color="0.92", lw=0, zorder=0)
    ax[1].scatter(range(len(S)), S.t, s=bs, c=col, ec="white", lw=.4, zorder=3)
    ax[1].axhline(tcrit, color="0.35", ls="--", lw=lw * .8, zorder=2)
    ax[1].text(0.5, tcrit + 0.12, f"Bonferroni t = {tcrit:.2f}", va="bottom", ha="left",
               fontsize=ann, color="0.35")
    keys = [Line2D([], [], ls="none", marker="o", ms=np.sqrt(bs) * 1.05, mec="white",
                   mew=.4, color=c, label=l)
            for c, l in [(JAVC, "Jurua alphavirus"), (ALT, "Barmah Forest virus"),
                         (OTHC, "other branches")]]
    ax[1].legend(handles=keys, loc="lower right", bbox_to_anchor=(1.0, legy),
                 handletextpad=.4, borderpad=.3, labelspacing=.55)
    ax[1].set_xlabel(f"Branches, ranked (n = {len(S)})")
    ax[1].set_ylabel("Studentized residual")
    ax[1].set_title("Genus-wide outlier scan", fontsize=base); ax[1].margins(0.07)

    # (c) bootstrap distribution. The numbers (0.21, CI, 198/200) live in the
    #     caption; only the IDENTITY of the red line stays in the panel, so a
    #     reader never has to leave the figure to learn what a mark is.
    ax[2].axvspan(lo95, hi95, color=JAVC, alpha=.12, lw=0, zorder=0)
    ax[2].hist(BT.resid, bins=22, color=HIST, ec="white", lw=.5, zorder=2)
    ax[2].axvline(OBS, color=JAVC, lw=lw * 1.4, zorder=4)
    ax[2].axvline(0, color="0.45", lw=lw * .8, ls=":", zorder=3)
    ax[2].annotate("observed", (OBS, ax[2].get_ylim()[1] * 0.97), xytext=(5, 0),
                   textcoords="offset points", ha="left", va="top",
                   fontsize=ann, color=JAVC)
    ax[2].set_xlabel("JAV E2 branch-length residual")
    ax[2].set_ylabel("Bootstrap replicates")
    ax[2].set_title("Bootstrap uncertainty", fontsize=base); ax[2].margins(0.07)

    for a, letter in zip(ax, "abc"):
        a.text(-0.14 if ncols == 1 else -0.22, 1.02, letter, transform=a.transAxes,
               fontweight="bold", fontsize=pl_size, va="bottom", ha="left")
    fig.tight_layout(pad=0.9, w_pad=1.7, h_pad=1.5)
    fig.savefig(path); fig.savefig(path.replace(".png", ".pdf"))

    # geometric check: no visible text box may overlap another, ignoring
    # legend-internal rows and tick labels sitting on their own axis
    r = fig.canvas.get_renderer()
    leg = {t for a in fig.axes if a.get_legend() for t in a.get_legend().get_texts()}
    tl = {t for a in fig.axes for t in a.get_xticklabels() + a.get_yticklabels()}
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    return [(a.get_text()[:24], b.get_text()[:24])
            for i, (a, ba) in enumerate(tx) for b, bb in tx[i + 1:]
            if ba.overlaps(bb) and not ({a, b} <= leg) and not ({a, b} <= tl)]


if __name__ == "__main__":
    print("paper       :", build("Figure2_divergence_paper.png", (7.2, 2.7), 1, 3,
                                 (9, 8, 7.5), 11, bs=16) or "no text overlaps")
    print("dissertation:", build("Figure2_divergence_dissertation.png", (6.5, 8.4), 3, 1,
                                 (12, 11, 10), 14, lw=1.3, ms=(48, 62), bs=30) or "no text overlaps")
