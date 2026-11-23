#!/usr/bin/env python3
"""Per-site variability across the JAV E2 ectodomain.

Computes Shannon entropy per alignment column over unambiguous, non-gap residues, maps
each column to its JAV E2 residue number, and smooths with a centred nine-residue window.
Columns with fewer than MIN_N scorable residues are left undefined; this includes the
JAV-only insertions, where no other taxon has a residue to compare. Entropy is reported
in bits, as the effective number of amino acids (2^H), and normalised to log2(20).

Usage : python scripts/e2_entropy.py            (run from the repository root)
Input : data/alignments/E2_plus_midd.aln        MAFFT L-INS-i, 14-taxon E2 set
Output: data/tables/e2_variability_windowed.csv
"""
import numpy as np, pandas as pd
from collections import Counter
from scipy.stats import entropy
from Bio import SeqIO

ALN, DROP, JAVK = "data/alignments/E2_plus_midd.aln", ("Whataroa",), "Jurua"
SKIP  = set("-XBZJUO")   # gaps and ambiguity codes, excluded from frequencies
ECTO  = 371              # last ectodomain residue (TM from 372, DeepTMHMM)
W, MIN_N, MIN_W = 9, 2, 5

aln = {r.id: str(r.seq) for r in SeqIO.parse(ALN, "fasta")
       if not any(d in r.id for d in DROP)}
jav = next(k for k in aln if JAVK in k)
keys, L = list(aln), len(aln[jav])

H = np.full(L, np.nan)
for i in range(L):
    aa = [aln[k][i] for k in keys if aln[k][i] not in SKIP]
    if len(aa) >= MIN_N:                       # JAV-only insertions fall below this
        H[i] = entropy(np.array(list(Counter(aa).values()), float), base=2)

jpos, m = {}, 0
for i, ch in enumerate(aln[jav]):
    if ch != "-":
        m += 1
        jpos[i] = m
Hjav = {jpos[i]: H[i] for i in jpos}

ecto, winm = [p for p in sorted(Hjav) if p <= ECTO], {}
for p in ecto:
    seg = [Hjav[q] for q in range(p - W // 2, p + W // 2 + 1) if q in Hjav and q <= ECTO]
    if np.sum(~np.isnan(seg)) >= MIN_W:
        winm[p] = float(np.nanmean(seg))

L20 = np.log2(20)
nn = lambda v, f: None if np.isnan(v) else f(v)
T = pd.DataFrame({
    "JAV_E2_position":            ecto,
    "shannon_entropy_bits":       [nn(Hjav[p], lambda v: round(v, 3))     for p in ecto],
    "effective_residues":         [nn(Hjav[p], lambda v: round(2**v, 2))  for p in ecto],
    "normalized_entropy_log20":   [nn(Hjav[p], lambda v: round(v/L20, 3)) for p in ecto],
    "window9_mean_bits":          [round(winm[p], 3)     if p in winm else None for p in ecto],
    "window9_effective_residues": [round(2**winm[p], 2)  if p in winm else None for p in ecto],
    "window9_normalized_log20":   [round(winm[p]/L20, 3) if p in winm else None for p in ecto],
})
T["window_variability_rank"] = T.window9_mean_bits.rank(ascending=False).round().astype("Int64")
T.to_csv("data/tables/e2_variability_windowed.csv", index=False)

print(f"{len(aln)} taxa x {L} columns | ectodomain mean {np.mean(list(winm.values())):.3f} bits "
      f"| max attainable {np.log2(len(aln)):.3f}")
print("undefined columns (JAV-unique):", sorted(p for p in Hjav if np.isnan(Hjav[p])))
