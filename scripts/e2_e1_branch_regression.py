#!/usr/bin/env python3
"""E2-versus-E1 divergence as a partitioned branch-length analysis.

Steps:
  1. extract the E1 and E2 domains from the structural polyproteins with the profile
     HMMs (hmmsearch, best-scoring domain, envelope coordinates)
  2. align each domain set separately (MAFFT L-INS-i) and concatenate
  3. fit an edge-unlinked two-partition model on the fixed topology, giving an
     independent E1 and E2 branch length for every branch. IQ-TREE is called with -sp
     (alias -Q); -p / -spp is edge-proportional and would force E2 = const x E1.
     Per-partition models are set in MODELS and read from parts.nex, since -m is not
     parsed alongside -sp.
  4. regress E2 branch length on E1 branch length and evaluate every branch as an
     externally studentized residual from a leave-one-out fit, Bonferroni-corrected
  5. bootstrap alignment columns within each partition, re-estimating branch lengths on
     the fixed topology, for a confidence interval on the JAV residual

Requires hmmsearch, mafft and iqtree on PATH. Intermediate files are written to work/.

Usage  : python scripts/e2_e1_branch_regression.py     (run from the repository root)
Inputs : data/alignments/sp_dataset.fasta, data/scan/E1.hmm, data/scan/E2.hmm,
         data/trees/sp_tree_35taxa.treefile
Outputs: data/tables/branch_lengths_35taxa.csv         one row per branch
         data/tables/bootstrap_jav_residual_200reps.csv  one row per replicate
"""
import io, random, subprocess, os
import numpy as np, pandas as pd
from Bio import SeqIO, Phylo
from scipy.stats import linregress, t as tdist

TOPO = "data/trees/sp_tree_35taxa.treefile"
SP35 = "data/alignments/sp_dataset.fasta"
WORK = "work"
OUT_BRANCHES = "data/tables/branch_lengths_35taxa.csv"
OUT_BOOT = "data/tables/bootstrap_jav_residual_200reps.csv"
NBOOT, SEED = 200, 1
MODELS = "Q.INSECT+I+G4:E1, LG+I+G4:E2"      # BIC-selected, then fixed for the bootstrap
JAV = "Jurua_alphavirus"

# ---------- 1-2. domains -> per-gene alignments -> concatenate ----------
def extract(hmm, fasta, out, domE="1e-10"):
    subprocess.run(f"hmmsearch --domE {domE} --domtblout {out}.dt -o /dev/null {hmm} {fasta}",
                   shell=True, check=True)
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(fasta, "fasta")}
    best = {}
    for line in open(f"{out}.dt"):
        if line.startswith("#"): continue
        f = line.split(); tname, score, ef, et = f[0], float(f[13]), int(f[19]), int(f[20])
        if tname not in best or score > best[tname][0]: best[tname] = (score, ef, et)
    with open(out, "w") as fh:
        for tname, (_, ef, et) in best.items(): fh.write(f">{tname}\n{seqs[tname][ef-1:et]}\n")
    return len(best)

def write_phy(path, aln, order, i1, i2):
    W = max(len(t) for t in order) + 3           # names can exceed any fixed width
    with open(path, "w") as fh:
        fh.write(f" {len(order)} {len(i1)+len(i2)}\n")
        for t in order:
            s = "".join(aln["E1"][t][j] for j in i1) + "".join(aln["E2"][t][j] for j in i2)
            fh.write(f"{t.ljust(W)}{s}\n")

os.makedirs(WORK, exist_ok=True)
for g in ("E1", "E2"):
    extract(f"data/scan/{g}.hmm", SP35, f"{WORK}/{g}.faa")
    subprocess.run(f"mafft --localpair --maxiterate 1000 --quiet {WORK}/{g}.faa > {WORK}/{g}.aln",
                   shell=True, check=True)
aln = {g: {r.id: str(r.seq) for r in SeqIO.parse(f"{WORK}/{g}.aln", "fasta")} for g in ("E1", "E2")}
order = sorted(set(aln["E1"]) & set(aln["E2"]))
n1, n2 = len(next(iter(aln["E1"].values()))), len(next(iter(aln["E2"].values())))
write_phy(f"{WORK}/concat.phy", aln, order, list(range(n1)), list(range(n2)))
# models go INSIDE the nexus file; -m on the command line is not parsed with -sp/-Q
open(f"{WORK}/parts.nex", "w").write(
    f"#nexus\nbegin sets;\n  charset E1 = 1-{n1};\n  charset E2 = {n1+1}-{n1+n2};\n"
    f"  charpartition mine = {MODELS};\nend;\n")

# ---------- 3. edge-unlinked branch lengths on the fixed topology ----------
def fit_branches(phy, pre):
    subprocess.run(f"iqtree -s {phy} -sp {WORK}/parts.nex -te {TOPO} -pre {pre} -nt 1 -quiet -redo",
                   shell=True, check=True)
    trees = [l.strip() for l in open(f"{pre}.parttrees") if l.strip()]
    out = []
    for nwk in trees:
        tr = Phylo.read(io.StringIO(nwk), "newick")
        out.append({frozenset(x.name for x in c.get_terminals()): c.branch_length
                    for c in tr.find_clades() if c.branch_length is not None})
    return out            # [E1_lengths, E2_lengths], keyed by descendant-tip set

b1, b2 = fit_branches(f"{WORK}/concat.phy", f"{WORK}/main")
keys = sorted(set(b1) & set(b2), key=lambda k: (len(k), sorted(k)))
B = pd.DataFrame([{"taxa": "|".join(sorted(k)) if len(k) == 1 else f"clade_{len(k)}",
                   "terminal": len(k) == 1, "E1": b1[k], "E2": b2[k]} for k in keys])

# ---------- 4. leave-one-out externally studentized residuals ----------
def loo(df):
    rows = []
    for i in df.index:
        bg = df.drop(i)
        f = linregress(bg.E1, bg.E2)
        resid = df.loc[i, "E2"] - (f.slope * df.loc[i, "E1"] + f.intercept)
        s = (bg.E2 - (f.slope * bg.E1 + f.intercept)).std(ddof=2)
        n, xb = len(bg), bg.E1.mean(); Sxx = ((bg.E1 - xb) ** 2).sum()
        se = s * np.sqrt(1 + 1/n + (df.loc[i, "E1"] - xb) ** 2 / Sxx)
        tv = resid / se; p = tdist.sf(abs(tv), n - 2) * 2
        rows.append(dict(resid=resid, t=tv, p_raw=p, p_bonf=min(1, p * len(df)),
                         leverage=1/n + (df.loc[i, "E1"] - xb) ** 2 / Sxx))
    return pd.DataFrame(rows, index=df.index)

B = B.join(loo(B)).sort_values("t", ascending=False)
B.round(6).to_csv(OUT_BRANCHES, index=False)

# ---------- 5. column bootstrap, branch lengths only ----------
def jav_residual(d1, d2):
    kk = [k for k in keys if k in d1 and k in d2]
    x = np.array([d1[k] for k in kk]); y = np.array([d2[k] for k in kk])
    m = np.array([k != frozenset({JAV}) for k in kk])
    f = linregress(x[m], y[m]); r = y - (f.slope * x + f.intercept)
    s = r[m].std(ddof=2); n, xb = m.sum(), x[m].mean(); Sxx = ((x[m] - xb) ** 2).sum()
    i = kk.index(frozenset({JAV}))
    return r[i], r[i] / (s * np.sqrt(1 + 1/n + (x[i] - xb) ** 2 / Sxx)), f.slope

random.seed(SEED); os.makedirs(f"{WORK}/boot", exist_ok=True); reps = []
for b in range(NBOOT):
    write_phy(f"{WORK}/boot/r{b}.phy", aln, order,
              [random.randrange(n1) for _ in range(n1)],
              [random.randrange(n2) for _ in range(n2)])
    try:
        d1, d2 = fit_branches(f"{WORK}/boot/r{b}.phy", f"{WORK}/boot/b{b}")
        reps.append(jav_residual(d1, d2))
    except Exception:
        pass
BT = pd.DataFrame(reps, columns=["resid", "t", "slope"])
BT.round(6).to_csv(OUT_BOOT, index=False)

obs = jav_residual(b1, b2)
lo, hi = np.percentile(BT.resid, [2.5, 97.5])
print(f"branches {len(B)} | background slope {linregress(B.E1[B.taxa != JAV], B.E2[B.taxa != JAV]).slope:.3f}")
print(f"surviving Bonferroni: {list(B.taxa[B.p_bonf < 0.05])}")
print(f"JAV residual {obs[0]:.4f} (t = {obs[1]:.2f}) | bootstrap 95% CI [{lo:.3f}, {hi:.3f}] "
      f"| {int((BT.resid > 0).sum())}/{len(BT)} replicates > 0")
