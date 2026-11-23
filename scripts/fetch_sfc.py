#!/usr/bin/env python3
"""Retrieve all Semliki Forest complex protein records from GenBank.

Queries NCBI protein for each of the eight complex species and writes one FASTA per
species. Species already present as a non-empty file are skipped, so the script can be
re-run to resume. Set EMAIL below before running; NCBI asks that E-utilities requests
identify the user.

The result depends on the retrieval date, since the database grows. The query is
"<species>"[Organism], recorded here so the counts reported in the paper can be
regenerated as retrieved.

Usage : python scripts/fetch_sfc.py             (run from the repository root)
Output: data/scan/<species>.faa, one file per species
"""
import urllib.request, urllib.parse, json, time, sys, os
EMAIL = None
SPECIES = ["Semliki Forest virus","Chikungunya virus","O'nyong-nyong virus","Ross River virus",
           "Getah virus","Bebaru virus","Mayaro virus","Una virus"]
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def params(d):
    if EMAIL: d = dict(d, email=EMAIL, tool="alphavirus_envelope_scan")
    return urllib.parse.urlencode(d)
def get(url, tries=6):
    for t in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=180).read().decode("utf-8","replace")
        except Exception as e:
            sys.stderr.write(f"retry {t}: {e}\n"); time.sleep(3+3*t)
    raise RuntimeError("failed: "+url[:120])

OUTDIR = "data/scan"
os.makedirs(OUTDIR, exist_ok=True)
for sp in SPECIES:
    out = f"{OUTDIR}/{sp.replace(' ','_').replace(chr(39),'')}.faa"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print("skip", sp, flush=True); continue
    q = f'"{sp}"[Organism]'
    r = json.loads(get(BASE+"esearch.fcgi?"+params(dict(db="protein", term=q, retmode="json",
                                                        retmax=0, usehistory="y"))))["esearchresult"]
    n, webenv, qk = int(r["count"]), r["webenv"], r["querykey"]
    print(f"{sp}: {n} records", flush=True)
    with open(out,"w") as fh:
        for start in range(0, n, 500):
            fh.write(get(BASE+"efetch.fcgi?"+params(dict(db="protein", WebEnv=webenv, query_key=qk,
                         retstart=start, retmax=500, rettype="fasta", retmode="text"))))
            time.sleep(0.4)
            if start % 2500 == 0: print(f"   {sp} {start}/{n}", flush=True)
    print(f"   done {sp}: {os.path.getsize(out)/1e6:.1f} MB", flush=True)
print("ALL DONE", flush=True)
