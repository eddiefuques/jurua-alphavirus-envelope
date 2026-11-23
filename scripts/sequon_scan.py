#!/usr/bin/env python3
"""Conservation of predicted N-glycosylation sequons across the Semliki Forest complex.

Scores the four JAV E2 sequons (N203, N223, N266, N321) and the E1 domain I sequon at
position 141 against every unique E2 and E1 sequence of the eight complex species.

Scoring runs in HMM match-state coordinates rather than raw protein coordinates, because
JAV's two E2 insertions offset its residue numbers from its relatives'. hmmalign assigns
every sequence to a common set of match states, which restores the correspondence.

A2M input convention: uppercase = match state, '-' = deletion, lowercase = insertion
(occupies no match state). A sequence is scored at a site only if all three positions are
unambiguous and ungapped, so denominators vary between sites.

Usage : python scripts/sequon_scan.py           (run from the repository root)
Inputs : data/scan/E2_all.a2m, data/scan/E1_all.a2m
Output : data/tables/sequon_conservation_sfc.csv
"""
import csv

E2_SEQUONS = {"N203": 203, "N223": 223, "N266": 266, "N321": 321}
E1_LABEL, E1_RESIDUE = "N141 (E1)", 141
SPECIES = ["Bebaru_virus", "Chikungunya_virus", "Getah_virus", "Mayaro_virus",
           "Onyong-nyong_virus", "Ross_River_virus", "Semliki_Forest_virus", "Una_virus"]
AMBIGUOUS = set("XBZJU")


def read_a2m(path):
    seqs, key = {}, None
    for line in open(path):
        if line.startswith(">"):
            key = line[1:].strip(); seqs[key] = ""
        elif key is not None:
            seqs[key] += line.strip()
    return seqs


def match_states(seq):
    return [c for c in seq if c.isupper() or c == "-"]


def species_of(header):
    return next((sp for sp in SPECIES if sp in header), None)


def call_site(trimer):
    """True, False, or None if the site cannot be scored."""
    if len(trimer) < 3 or "-" in trimer or (AMBIGUOUS & set(trimer)):
        return None
    return trimer[0] == "N" and trimer[1] != "P" and trimer[2] in "ST"


def residue_to_match_state(jav_seq, residue):
    n, state = 0, -1
    for c in jav_seq:
        if c == "-":
            state += 1; continue
        if c.isupper():
            state += 1
        n += 1
        if n == residue:
            return state
    raise ValueError(f"residue {residue} not present")


def scan(path, jav_key, sites):
    seqs = read_a2m(path)
    states = {lab: residue_to_match_state(seqs[jav_key], pos) for lab, pos in sites.items()}
    tally = {sp: {lab: [0, 0] for lab in sites} for sp in SPECIES}
    jav = match_states(seqs[jav_key])
    jav_calls = {lab: call_site("".join(jav[s:s + 3])) for lab, s in states.items()}
    for header, seq in seqs.items():
        sp = species_of(header)
        if sp is None:
            continue
        chars = match_states(seq)
        for lab, s in states.items():
            call = call_site("".join(chars[s:s + 3]))
            if call is None:
                continue
            tally[sp][lab][1] += 1
            tally[sp][lab][0] += int(call)
    return tally, jav_calls, states, len(seqs) - 1


if __name__ == "__main__":
    t2, jav2, s2, n2 = scan("data/scan/E2_all.a2m", "JAV_E2", E2_SEQUONS)
    t1, jav1, s1, n1 = scan("data/scan/E1_all.a2m", "JAV_E1", {E1_LABEL: E1_RESIDUE})

    labels = list(E2_SEQUONS) + [E1_LABEL]
    rows = []
    for sp in SPECIES:
        row = {"Virus": sp}
        for lab in E2_SEQUONS:
            row[lab] = "{}/{}".format(*t2[sp][lab])
        row[E1_LABEL] = "{}/{}".format(*t1[sp][E1_LABEL])
        rows.append(row)
    jav = {"Virus": "Jurua alphavirus"}
    for lab in labels:
        jav[lab] = "present" if (jav2 | jav1)[lab] else "absent"
    rows.append(jav)

    out = "data/tables/sequon_conservation_sfc.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["Virus"] + labels)
        w.writeheader(); w.writerows(rows)

    print(f"E2 sequences examined: {n2} | E1: {n1}")
    print("match states:", {**s2, **s1})
    width = max(len(r["Virus"]) for r in rows)
    for r in rows:
        print("  " + r["Virus"].ljust(width) + "  " + "  ".join(f"{r[l]:>10}" for l in labels))
    for lab in labels:
        src = t2 if lab in E2_SEQUONS else t1
        w_ = sum(src[sp][lab][0] for sp in SPECIES); n_ = sum(src[sp][lab][1] for sp in SPECIES)
        print(f"  TOTAL {lab:10} {w_}/{n_}  ({100*w_/n_:.1f}%)")
