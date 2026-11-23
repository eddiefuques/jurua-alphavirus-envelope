# jurua-alphavirus-envelope

Analysis code and data for the comparative genomic and structural characterization of the
envelope glycoproteins of Juruá alphavirus (JAV), a novel alphavirus recovered from
mosquitoes in the western Amazon Basin.

This repository accompanies:

> Fuques E, Massey AL, Levi T, Trubl G, Vega Thurber RL. Genomic and structural
> characterization of a novel Amazonian alphavirus reveals lineage-specific changes in
> both envelope glycoproteins. Manuscript submitted.

The virus genome is GenBank **PV335208**; the underlying reads are BioProject
**PRJNA1223432**. Both were deposited with the original description of the virus and are
not new to this study.

## Layout

```
scripts/     analyses that produce the results
figures/     standalone scripts for the figures drawn directly from derived tables
model/       the predicted JAV E1/E2 heterodimer, with its confidence data
data/
  alignments/   input alignments and sequence sets
  trees/        maximum likelihood trees with support values
  tables/       derived tables consumed and produced by the scripts
  scan/         profile HMMs and extracted domain sets for the conservation scan
```

## Quickstart

Every script is run from the repository root and resolves its inputs relative to it.

```bash
conda env create -f environment.yml && conda activate jav-envelope
python figures/figure2_divergence.py        # writes both figure versions as PNG and PDF
```

Each script states its inputs, outputs and usage in its own header.
`scripts/e2_e1_branch_regression.py` additionally needs `hmmsearch`, `mafft` and `iqtree`
on PATH, and writes its intermediate files to `work/`.

## What the scripts do

| script | what it does |
|---|---|
| `scripts/fetch_sp_dataset.py` | rebuilds the 35-taxon structural polyprotein set from GenBank using the accessions in `data/tables` |
| `scripts/fetch_sfc.py` | retrieves the Semliki Forest complex protein records used for the conservation scans, one FASTA per species |
| `scripts/e2_e1_branch_regression.py` | fits E1 and E2 as edge-unlinked partitions on the 35-taxon topology, scans all 67 branches as externally studentized residuals from a leave-one-out fit, and bootstraps the JAV residual |
| `scripts/e2_entropy.py` | per-column Shannon entropy across the E2 ectodomain, reported as effective number of amino acids over a nine-residue window |
| `scripts/sequon_scan.py` | scores the four JAV E2 sequons and the E1 domain I sequon against every unique complex E2 and E1 sequence, in HMM match-state coordinates |
| `figures/*.py` | one standalone script per scripted figure (Figures 2, 5, 6), each reading only files in `data/` |

## Scope

This repository contains the analysis code written specifically for this study: the steps
that are not fully specified by the Methods and therefore cannot be reconstructed from the
published description alone. Standard steps run with established tools under stated
parameters, such as alignment with MAFFT, tree inference with IQ-TREE, recombination
screening, and alignment rendering with ESPript, are documented in the Methods and are not
reproduced here as code. What is included are the custom analyses whose implementation
details affect the result, together with the inputs they consume and the derived tables
they produce.

`figures/` covers the figures drawn directly from those derived tables: Figure 2, Figure 5
and Figure 6. Figures built by arranging alignment renderings and structural views produced
by the tools named in the Methods are not scripted here.

## Data files and what they correspond to

Files in `data/tables/` are the inputs and outputs of the scripts, named for their
contents. The numbered, captioned versions of the supplementary tables are published with
the paper; they are not duplicated here.

| file | contents | appears in the paper as |
|---|---|---|
| `sequence_accessions.csv` | accessions for all three comparison sets, by analysis | Supplementary Table S1 |
| `e2_e1_pairwise_distances_9taxa.csv` | pairwise amino-acid p-distances, E2 against E1, nine complex taxa | Figure 2a |
| `branch_lengths_35taxa.csv` | edge-unlinked E1 and E2 branch lengths and residual statistics for all 67 branches | Figure 2b, Supplementary Table S4 |
| `bootstrap_jav_residual_200reps.csv` | JAV residual across 200 column-resampled replicates | Figure 2c |
| `e2_variability_windowed.csv` | per-residue effective number of amino acids across the E2 ectodomain | Figure 5 |
| `receptor_contact_distances.csv` | per-residue minimum heavy-atom distance from JAV to each reference receptor | Figure 6 |
| `sequon_conservation_sfc.csv` | sequon counts per species at the four E2 sites and E1 N141 | Table 1 |

The superposition distances in `receptor_contact_distances.csv` were computed from the
model in `model/` against PDB 6JO8 and 9RBQ; that step is not scripted here.

## The structural model

`model/` holds the AlphaFold 3 heterodimer underlying the structural analyses, together
with its confidence data and the job request that produced it, described in
`model/README.md`. These files are AlphaFold Server output and carry their own terms of
use; they are not covered by the licence that applies to `data/`.

## External tools

The alignments, trees and scans in `data/` were produced with MAFFT v7, HMMER 3,
IQ-TREE 2, OpenRDP and 3SEQ, MEGA11, DeepTMHMM, NetSurfP-3.0, AlphaFold 3 and ChimeraX.
Versions and parameters are given in the paper's Methods. The scripts here consume those
outputs rather than re-running the external tools, except where noted.

## Reproducibility note

The conservation scans in `data/scan` are drawn from a GenBank query whose result grows as
new alphavirus sequences are deposited. `scripts/fetch_sfc.py` records that query verbatim,
and the retrieval date is stated in the Methods, so the denominators reported in the paper
can be regenerated exactly as retrieved. Re-running the fetch at a later date will return a
larger set and correspondingly different denominators.

## License

Code is released under the MIT License. Data files are released under CC BY 4.0.
Sequences retrieved from GenBank remain subject to their original terms.
