#!/usr/bin/env python3
"""Figure 6 - JAV envelope residues at the receptor-contact surfaces of relatives.

Per-residue minimum heavy-atom distance from JAV E2 to MXRA8 (superposed on CHIKV
E2-MXRA8, 6JO8) and from JAV E1 to ApoER2 LA5 (superposed on SFV E1-ApoER2, 9RBQ).
Contacts are defined at 5 A in the reference complex. Markers distinguish residues
conserved or substituted in JAV; open circles mark predicted sequon positions, and bands
mark the two JAV-specific E2 insertions.

Builds two versions from the same data. The dissertation version must be placed at 100 %
of the text block to keep its 10 pt type floor.

Usage  : python figures/figure6_superposition.py       (run from the repository root)
Input  : data/tables/receptor_contact_distances.csv
Outputs: Figure6_receptor_superposition_{paper,dissertation}.{png,pdf}
"""
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

DATA = "data/tables/receptor_contact_distances.csv"
CUT = 5.0                                   # contact definition, angstrom
CONS, SUBS, FEAT, BAND = "#1F5FA8", "#B8411F", "#1F6F4A", "#9BD0A8"
PANELS = [dict(protein="E2", title="E2 superposed on CHIKV E2\u2013MXRA8 (6JO8)",
               ylab="E2 \u2192 MXRA8\nmin distance (\u00c5)",
               bands=[(116, 118), (206, 206)],
               feats={203: ("N203", -17, 16), 223: ("N223", 17, 16),
                      266: ("N266", 0, 17), 321: ("N321", 0, 17)}),
          dict(protein="E1", title="E1 superposed on SFV E1\u2013ApoER2 LA5 (9RBQ)",
               ylab="E1 \u2192 ApoER2 LA5\nmin distance (\u00c5)",
               bands=[],
               feats={141: ("D141", 0, 17)})]

D = pd.read_csv(DATA)


def style(sizes):
    base, ann, tick = sizes
    mpl.rcParams.update({
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "legend.fontsize": ann, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.titlelocation": "left", "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
        "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42})


def build(path, W, H, sizes, ms=22, pl=12, ncol=4, legy=0.015, bottom=0.11):
    base, ann, tick = sizes
    style(sizes)
    fig, axes = plt.subplots(2, 1, figsize=(W, H))
    for ax, P in zip(axes, PANELS):
        d = D[D.protein == P["protein"]].sort_values("jav_residue")
        x, yv = d.jav_residue.values, d.min_distance_A.values
        ax.axhspan(0, CUT, color="0.90", lw=0, zorder=0)
        for a, b in P["bands"]:
            ax.axvspan(a - 0.8, b + 0.8, color=BAND, alpha=.75, lw=0, zorder=1)
        ax.plot(x, yv, color="#55585C", lw=.8, zorder=3)
        for flag, col, mk in [("yes", CONS, "o"), ("no", SUBS, "X")]:
            s = d[(d.receptor_contact == "yes") & (d.conserved_in_jav == flag)]
            ax.scatter(s.jav_residue, s.min_distance_A, s=ms, c=col, marker=mk,
                       ec="white", lw=.4, zorder=5)
        for pos, (lab, dx, dy) in P["feats"].items():
            row = d[d.jav_residue == pos]
            if row.empty:
                continue
            yy = float(row.min_distance_A.iloc[0])
            ax.scatter([pos], [yy], s=ms * 2.6, facecolors="none", edgecolors=FEAT,
                       lw=1.1, zorder=6)
            ax.annotate(lab, (pos, yy), xytext=(dx, dy), textcoords="offset points",
                        ha="center", va="center", fontsize=ann, color=FEAT, zorder=7)
        ax.set_title(P["title"], fontsize=base, pad=5)
        ax.set_ylabel(P["ylab"])
        ax.set_xlim(-6, x.max() + 6); ax.set_ylim(-4, yv.max() * 1.18)
    axes[1].set_xlabel("JAV residue number")
    for ax, s in zip(axes, "ab"):
        ax.text(-0.11, 1.04, s, transform=ax.transAxes, fontweight="bold",
                fontsize=pl, va="bottom", ha="left")
    handles = [Line2D([], [], ls="none", marker="o", ms=np.sqrt(ms) * .85, mfc=CONS,
                      mec="white", mew=.4),
               Line2D([], [], ls="none", marker="X", ms=np.sqrt(ms) * .95, mfc=SUBS,
                      mec="white", mew=.4),
               Patch(fc=BAND, ec="none", alpha=.75),
               Line2D([], [], ls="none", marker="o", ms=np.sqrt(ms * 2.6) * .9,
                      mfc="none", mec=FEAT, mew=1.1)]
    labels = ["receptor contact, conserved in JAV",
              "receptor contact, substituted in JAV",
              "JAV-specific insertion (E2)",
              "N-glycosylation sequon position"]
    fig.tight_layout(pad=0.7, h_pad=1.6, rect=(0, bottom, 1, 1))
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, legy),
               ncol=ncol, handletextpad=.5, columnspacing=1.5, borderpad=.3,
               labelspacing=.5)
    fig.savefig(path + ".png", dpi=600); fig.savefig(path + ".pdf")
    r = fig.canvas.get_renderer()
    leg = {t for lg in fig.legends for t in lg.get_texts()}
    tl = {t for a in fig.axes for t in a.get_xticklabels() + a.get_yticklabels()}
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    return [(a.get_text()[:22], b.get_text()[:22])
            for i, (a, ba) in enumerate(tx) for b, bb in tx[i + 1:]
            if ba.overlaps(bb) and not ({a, b} <= leg) and not ({a, b} <= tl)]


if __name__ == "__main__":
    for nm in ("E2", "E1"):
        d = D[D.protein == nm]
        c = d[d.receptor_contact == "yes"]
        print(f"{nm}: {len(d)} residues, {len(c)} contacts "
              f"({(c.conserved_in_jav=='yes').sum()} conserved, "
              f"{(c.conserved_in_jav=='no').sum()} substituted)")
    print("paper       :", build("Figure6_receptor_superposition_paper", 7.2, 5.2,
                                 (9, 8, 7.5), ncol=4, bottom=0.08) or "clean")
    print("dissertation:", build("Figure6_receptor_superposition_dissertation", 6.5, 6.0,
                                 (12, 11, 10), ms=34, pl=14, ncol=2, bottom=0.13) or "clean")
