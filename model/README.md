# Predicted JAV E1/E2 heterodimer

The AlphaFold 3 model underlying the structural analyses in this study.

| file | contents |
|---|---|
| `jav_E1_E2_heterodimer.cif` | coordinates. Chain A is E1 (438 residues), chain B is E2 (426 residues). The B-factor column holds per-atom pLDDT. |
| `jav_E1_E2_summary_confidences.json` | model-level pTM, ipTM, per-chain and chain-pair values |
| `jav_E1_E2_full_data.json` | per-atom pLDDT and the predicted aligned error matrix |
| `jav_E1_E2_job_request.json` | the submitted job: both input sequences and the model seed |
| `AlphaFold_Server_terms_of_use.md` | the terms these files are distributed under |

Generated with the AlphaFold 3 server (Abramson et al. 2024), job
`dimmer_E1_E2_complete_seqs`, seed 2025359455, run 20 February 2026. The server returned
five ranked predictions; this is model 0. Residue numbering matches the paper throughout
(E2 1-426, E1 1-438).

Per-residue pLDDT is the mean over all atoms of a residue, which is how the values in
Supplementary Table S6 were computed. Model-level scores are pTM 0.80 and ipTM 0.81.

The server-generated multiple sequence alignments and structure templates are not included
here. The model can be used as supplied, but regenerating it from the job request would
require those inputs.

## Terms

These files are AlphaFold Server output and are distributed under the AlphaFold Server
terms of use, reproduced in `AlphaFold_Server_terms_of_use.md`. Those terms permit
non-commercial use only and require that any further distribution carry the same terms and
notice. They are not covered by the CC BY 4.0 licence that applies to `data/`.
