#!/usr/bin/env python3
"""Figure 5 - sequence variability across the JAV E2 ectodomain.

Plots the effective number of amino acids (2^H, H = Shannon entropy) averaged over a
centred nine-residue window, for each of the 371 ectodomain residues. Marks the two
JAV-specific insertions, the four predicted sequons, the hypervariable loop and the
antigenic segment.

Builds two versions from the same data. The dissertation version must be placed at 100 %
of the text block to keep its 10 pt type floor.

Usage  : python figures/figure5_e2_variability.py      (run from the repository root)
Input  : data/tables/e2_variability_windowed.csv
Outputs: Figure5_E2_variability_{paper,dissertation}.{png,pdf}
"""
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

DATA = "data/tables/e2_variability_windowed.csv"
INS = {117: "E116-V117-G118", 206: "A206"}                 # lineage-specific insertions
SEQ = {203: "N203", 223: "N223", 266: "N266", 321: "N321"} # predicted sequons
BANDS = [(119, 127, "#E8C56A", "hypervariable loop 119-127", "#7A6118", 1.26),
         (173, 224, "#A9C9E0", "antigenic segment 173-224", "#265A7A", 1.12)]
OFF = {"E116-V117-G118": (-62, 48), "A206": (-16, 40), "N223": (40, 12),
       "N203": (-22, -24), "N266": (4, -28), "N321": (16, 22)}
INSC, SEQC = "#B8411F", "#1F6F7A"

V = pd.read_csv(DATA)
x = V.JAV_E2_position.values
y = V.window9_effective_residues.values
MEAN = float(np.nanmean(y))


def style(sizes):
    base, ann, tick = sizes
    mpl.rcParams.update({
        "font.size": base, "axes.labelsize": base, "legend.fontsize": ann,
        "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42})


def build(path, W, H, sizes):
    base, ann, tick = sizes
    style(sizes)
    fig, ax = plt.subplots(figsize=(W, H))
    for a, b, c, lab, tc, ty in BANDS:
        ax.axvspan(a, b, color=c, alpha=.38, lw=0, zorder=0)
        ax.annotate(lab, xy=((a + b) / 2, y.max()), xytext=((a + b) / 2, y.max() * ty),
                    ha="center", va="center", fontsize=ann, color=tc,
                    arrowprops=dict(arrowstyle="-", lw=.6, color=tc, shrinkA=1, shrinkB=1))
    ax.axhline(MEAN, color="0.45", ls=":", lw=.9, zorder=1)
    ax.plot(x, y, color="#1D2228", lw=.9, zorder=3)
    for d, col in ((INS, INSC), (SEQ, SEQC)):
        for p, lab in d.items():
            ax.plot(p, y[p - 1], "o", ms=3.4, color=col, mec="white", mew=.5, zorder=5)
            ax.annotate(lab, (p, y[p - 1]), xytext=OFF[lab], textcoords="offset points",
                        fontsize=ann, color=col, ha="center", va="center", zorder=6,
                        arrowprops=dict(arrowstyle="-", lw=.7, color=col, shrinkA=1, shrinkB=2))
    ax.set_xlabel("JAV E2 residue")
    ax.set_ylabel("Effective no. of amino acids")
    ax.set_xlim(-6, 378); ax.set_ylim(0.9, y.max() * 1.34)
    fig.tight_layout(pad=0.6)
    fig.savefig(path + ".png", dpi=600); fig.savefig(path + ".pdf")
    r = fig.canvas.get_renderer()
    tl = {t for t in ax.get_xticklabels() + ax.get_yticklabels()}
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    return [(a.get_text()[:18], b.get_text()[:18])
            for i, (a, ba) in enumerate(tx) for b, bb in tx[i + 1:]
            if ba.overlaps(bb) and not ({a, b} <= tl)]


if __name__ == "__main__":
    print(f"mean = {MEAN:.2f} effective amino acids; range {y.min():.2f}-{y.max():.2f}")
    print("paper       :", build("Figure5_E2_variability_paper", 7.2, 2.8, (9, 8, 7.5)))
    print("dissertation:", build("Figure5_E2_variability_dissertation", 6.5, 3.3, (12, 11, 10)))
