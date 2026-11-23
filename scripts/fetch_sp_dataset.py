#!/usr/bin/env python3
"""Rebuild the 35-taxon structural polyprotein FASTA from GenBank.

Reads the accessions listed in data/tables/sequence_accessions.csv and fetches each
protein record. The Juruá alphavirus entry is given as a nucleotide range rather than a
protein accession; it is fetched from the genome record and translated from the CDS
coordinates in JAV_CDS.

Usage : python scripts/fetch_sp_dataset.py      (run from the repository root)
Input : data/tables/sequence_accessions.csv
Output: data/alignments/sp_dataset.fasta
"""
import csv, re, time, urllib.parse, urllib.request

TABLE = "data/tables/sequence_accessions.csv"
OUT = "data/alignments/sp_dataset.fasta"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
JAV_CDS = (7187, 10873)          # 1-based inclusive, from the genome annotation


def efetch(db, acc, rettype):
    q = urllib.parse.urlencode(dict(db=db, id=acc, rettype=rettype, retmode="text"))
    for attempt in range(4):
        try:
            with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=60) as fh:
                return fh.read().decode()
        except Exception:
            time.sleep(2 + 3 * attempt)
    raise RuntimeError(f"could not fetch {acc} from {db}")


def seq_of(fasta):
    return "".join(l.strip() for l in fasta.splitlines() if not l.startswith(">"))


CODONS = {}
for i, b1 in enumerate("TCAG"):
    for b2 in "TCAG":
        for b3 in "TCAG":
            CODONS[b1 + b2 + b3] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
                len(CODONS)]


def translate(nt):
    aa = "".join(CODONS.get(nt[i:i + 3].upper(), "X") for i in range(0, len(nt) - 2, 3))
    return aa[:-1] if aa.endswith("*") else aa


rows = list(csv.DictReader(open(TABLE)))
sp = [r for r in rows if "polyprotein" in r["analysis"].lower()]
written = 0
with open(OUT, "w") as out:
    for r in sp:
        name = r["virus"].strip().replace(" ", "_")
        acc = r["accession"].split()[0]
        if re.search(r"translated CDS", r["accession"]):
            lo, hi = JAV_CDS
            aa = translate(seq_of(efetch("nuccore", acc, "fasta"))[lo - 1:hi])
        else:
            aa = seq_of(efetch("protein", acc, "fasta"))
        out.write(f">{name}|{acc}\n{aa}\n")
        written += 1
        time.sleep(0.4)
print(f"wrote {written} structural polyproteins to {OUT}")
