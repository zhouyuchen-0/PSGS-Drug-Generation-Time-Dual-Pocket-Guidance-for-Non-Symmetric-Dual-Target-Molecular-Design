# Configuration files

This folder contains public configuration records for the formal PSGS-Drug runs, controlled ablation groups, and single-pocket controls reported in the manuscript.

All local absolute paths have been converted to repository-relative paths for public release. If users run the workflow in a different directory, update the receptor, pocket-prior, output, and log paths accordingly.

## Files

| File | Role |
|---|---|
| `setting_prior.yaml` | Formal prior-guided population set. |
| `setting_prior_contactfrag.yaml` | Formal contact-guided final-generation set. |
| `setting_ablation_A0_noprior_noseed.yaml` | A0 no-prior/no-seed ablation. |
| `setting_ablation_A1_prior_only.yaml` | A1 prior-only ablation. |
| `setting_ablation_A2_prior_plus_seed.yaml` | A2 prior-plus-contact-prefix ablation. |
| `setting_baseline_B1_3fap_only.yaml` | B1 3FAP-first single-pocket control. |
| `setting_baseline_B2_7pqv_only.yaml` | B2 7PQV-first single-pocket control. |

## Path conventions

- Receptor and pocket-prior files are expected under `receptors/`.
- Generated outputs are written under `outputs/`.
- The pocket-prior service is assumed to run locally at `http://127.0.0.1:26974/prior`.

## Notes

These files document the reported experiments and can be used as configuration templates for rerunning the workflow. Some upstream generation components may require local infrastructure or third-party model code, as described in the main repository README and Supporting Information.
