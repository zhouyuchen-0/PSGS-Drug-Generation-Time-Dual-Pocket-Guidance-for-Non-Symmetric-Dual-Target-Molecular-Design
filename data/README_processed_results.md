# Processed result CSV files

This folder contains processed CSV files associated with the PSGS-Drug manuscript.

## formal_sets/

- `Prior-guided population set.csv`: formal protein-pocket-prior-guided population set used for population-level directionality and global risk-profile analysis. The uploaded version removes the non-portable RDKit object-string column, if present, and includes `dock_sum` as an alias of `dock_combined` when needed.
- `Contact-guided final set.csv`: formal contact-guided final-generation set used for full-set docking, Top-k enrichment, candidate prioritization, and binding-region analysis.

## controls_and_ablation/

- `result_3fap-7pqv_A0_noprior_noseed.csv`: A0 ablation output; no pocket prior and no seed.
- `result_3fap-7pqv_A1_prior_only.csv`: A1 ablation output; dual-pocket prior without seed.
- `result_3fap-7pqv_A2_prior_plus_seed.raw.csv`: A2 raw ablation output; dual-pocket prior with contact-derived seed fragments. This file is provided as the raw A2 generation output before final deduplication.
- `result_B1_3fap_only_postscreen_7pqv.csv`: B1 single-pocket control; 3FAP-first generation followed by retrospective 7PQV screening.
- `result_B2_7pqv_only_postscreen_3fap.csv`: B2 single-pocket control; 7PQV-first generation followed by retrospective 3FAP screening.

## Notes

These files contain processed molecule-level outputs used for manuscript tables, figures, and Supporting Information analyses. They do not contain local absolute paths, credentials, API keys, or private user information. Docking scores are computational prioritization signals and should not be interpreted as experimental binding measurements.

