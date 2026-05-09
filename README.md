# PSGS-Drug Supporting Data

This package contains the data, configuration files, and analysis scripts needed to reproduce the main computational results reported in:

**PSGS-Drug: Generation-Time Dual-Pocket Guidance for Non-Symmetric Dual-Target Molecular Design**

PSGS-Drug evaluates generation-time dual-pocket molecular design for a non-symmetric MEK1 ATP-site / mTOR FRB-site task represented by **7PQV** and **3FAP**. The package is organized to support reproduction of the reported formal result sets, controls, ablations, docking summaries, Top-k analyses, medicinal-chemistry/ADMET profiling, PLIP hotspot-overlap analysis, and higher-exhaustiveness repeated-docking checks.

---

## 1. Package contents

Recommended directory structure:

```text
PSGS_Drug_Supporting_Data/
│
├── README.md
├── environment/
│   ├── environment.yml
│   └── requirements.txt
│
├── configs/
│   ├── setting_prior.yaml
│   ├── setting_prior_contactfrag_run_contact_v2_round2_formal.yaml
│   ├── setting_ablation_A0_noprior_noseed.yaml
│   ├── setting_ablation_A1_prior_only.yaml
│   ├── setting_ablation_A2_prior_plus_seed.yaml
│   ├── setting_baseline_B1_3fap_only.yaml
│   ├── setting_baseline_B2_7pqv_only.yaml
│   └── setting_pocket2mol_baseline_3fap_7pqv.yaml
│
├── receptors/
│   ├── 3fap_receptor.pdbqt
│   ├── 7pqv_receptor.pdbqt
│   ├── docking_boxes.yaml
│   └── receptor_preparation_README.md
│
├── formal_sets/
│   ├── Prior-guided population set.csv
│   └── Contact-guided final set.csv
│
├── controls_and_ablation/
│   ├── result_3fap-7pqv_A0_noprior_noseed.csv
│   ├── result_3fap-7pqv_A1_prior_only.csv
│   ├── result_3fap-7pqv_A2_prior_plus_seed.raw.csv
│   ├── A0_A1_global_ablation_summary.csv
│   ├── A1_A2_topk_priority_summary.csv
│   ├── B1_B2_contact_full_summary.csv
│   └── B1_B2_contact_main_matrix_data.csv
│
├── pocket2mol/
│   ├── setting_pocket2mol_baseline_3fap_7pqv.yaml
│   ├── result_pocket2mol_merged_1889.csv
│   ├── result_pocket2mol_dual_docked_1889.csv
│   └── metrics_pocket2mol_1889.csv
│
├── topk/
│   ├── topk_by_docksum_prior_vs_contact.csv
│   ├── topk_by_integrated_priority_prior_vs_contact.csv
│   ├── contact_top100_by_integrated_priority.csv
│   ├── contact_top100_molecules_by_docksum.csv
│   ├── contact_top100_molecules_by_priority.csv
│   └── top20_candidates_for_PLIP_interaction_analysis.csv
│
├── admet_descriptors/
│   ├── Contact-guided final-generation set-top100.csv
│   ├── descriptor_distribution_summary_prior_vs_contact.csv
│   ├── full_quality_summary_prior_vs_contact.csv
│   └── admetlab3_contact_top100_raw.csv
│
├── scaffold_alerts/
│   ├── Contact-guided final set-bm_scaffold_summary.csv
│   ├── Contact-guided final set-ring-metric-summary.txt
│   ├── Contact-guided final set-summary_pains.txt
│   ├── seed_fragment_scaffolds.csv
│   ├── seed_similarity_molecule_level.csv
│   └── seed_similarity_scaffold_summary.csv
│
├── plip/
│   ├── top6_balanced_multi_criteria_candidates.csv
│   ├── top20_candidates_for_PLIP_interaction_analysis.csv
│   ├── plip_summary_top20.csv
│   ├── reference_hotspot_residues.csv
│   └── representative_plip_reports/
│
├── repeated_docking/
│   ├── validation_high_exh_summary_table.csv
│   └── validation_high_exh_per_molecule.csv
│
├── baseline_generation_quality/
│   ├── baseline_reproduction_settings.csv
│   ├── generation_quality_baselines_summary.csv
│   └── generated_molecules/
│       ├── jtvae_1889_generated.csv
│       ├── mars_1889_generated.csv
│       ├── rationalerl_1889_generated.csv
│       ├── reinvent2_1889_generated.csv
│       └── vega_1889_generated.csv
│
└── scripts/
    ├── evaluate_generation_quality.py
    ├── calc_descriptors.py
    ├── prepare_ligands.py
    ├── run_dual_vina_docking.py
    ├── parse_vina_results.py
    ├── compute_dual_hit_rates.py
    ├── compute_topk_metrics.py
    ├── compute_integrated_priority.py
    ├── seed_similarity_scaffold_analysis.py
    ├── plip_hotspot_overlap_analysis.py
    ├── summarize_high_exh_redocking.py
    └── draw_figures.py
```

Some scripts, folders, or file names may differ slightly depending on the local execution environment. When this occurs, use the file-role descriptions below as the authoritative mapping.

---

## 2. Key result files

### 2.1 Formal PSGS-Drug result sets

| File | Description |
|---|---|
| `formal_sets/Prior-guided population set.csv` | Formal protein-pocket-prior-guided population set. This file corresponds to the main prior-guided population set used for population-level directionality and global risk-profile analysis. |
| `formal_sets/Contact-guided final set.csv` | Formal contact-guided final-generation set. This file corresponds to the final PSGS-Drug candidate-generation setting used for full-set docking, Top-k enrichment, candidate prioritization, and representative binding-region analysis. |

### 2.2 Controls and ablations

| File | Description |
|---|---|
| `controls_and_ablation/A0_A1_global_ablation_summary.csv` | Global ablation summary for no-prior/no-seed versus dual-pocket-prior-only generation. |
| `controls_and_ablation/A1_A2_topk_priority_summary.csv` | Top-k integrated-priority comparison between prior-only and prior-plus-contact-prefix ablation settings. |
| `controls_and_ablation/result_3fap-7pqv_A2_prior_plus_seed.raw.csv` | Raw prior-plus-contact-prefix ablation output. This file contains 1889 generated records and 1864 unique SMILES before final deduplication. |
| `controls_and_ablation/B1_B2_contact_full_summary.csv` | Summary for 3FAP-first and 7PQV-first single-pocket post-screening controls. |
| `controls_and_ablation/B1_B2_contact_main_matrix_data.csv` | Matrix data used for the B1/B2 control visualization. |

### 2.3 Pocket2Mol contextual comparator

| File | Description |
|---|---|
| `pocket2mol/setting_pocket2mol_baseline_3fap_7pqv.yaml` | Configuration record for Pocket2Mol contextual comparator generation and downstream evaluation. |
| `pocket2mol/result_pocket2mol_dual_docked_1889.csv` | Pocket2Mol molecules after merging, standardization, deduplication, and dual-target docking. |
| `pocket2mol/metrics_pocket2mol_1889.csv` | Summary metrics for Pocket2Mol comparator under the same downstream dual-target evaluation protocol. |

---

## 3. Configuration files

The `configs/` folder contains run-level configuration files for the formal PSGS-Drug result sets, ablation groups, single-pocket controls, and Pocket2Mol contextual comparator.

| Configuration file | Role |
|---|---|
| `setting_prior.yaml` | Formal prior-guided population set configuration. |
| `setting_prior_contactfrag_run_contact_v2_round2_formal.yaml` | Formal contact-guided final-generation set configuration. |
| `setting_ablation_A0_noprior_noseed.yaml` | A0 ablation: no pocket prior and no seed. |
| `setting_ablation_A1_prior_only.yaml` | A1 ablation: dual-pocket prior without contact-derived prefix. |
| `setting_ablation_A2_prior_plus_seed.yaml` | A2 ablation: dual-pocket prior with contact-derived prefix. |
| `setting_baseline_B1_3fap_only.yaml` | B1 control: 3FAP-only generation followed by retrospective dual-target screening. |
| `setting_baseline_B2_7pqv_only.yaml` | B2 control: 7PQV-only generation followed by retrospective dual-target screening. |
| `setting_pocket2mol_baseline_3fap_7pqv.yaml` | Pocket2Mol contextual comparator configuration. |

The formal prior-guided population set and A1 share the same conceptual module status but were generated under different complete settings. A1 should therefore be interpreted only within the matched A0–A1 ablation comparison.

---

## 4. Receptor and docking files

The `receptors/` folder contains receptor files and docking-box definitions used for the unified downstream docking protocol.

| File | Description |
|---|---|
| `3fap_receptor.pdbqt` | Prepared 3FAP receptor file for AutoDock Vina docking. |
| `7pqv_receptor.pdbqt` | Prepared 7PQV receptor file for AutoDock Vina docking. |
| `docking_boxes.yaml` | Docking box centers and sizes for 3FAP and 7PQV. |
| `receptor_preparation_README.md` | Notes describing receptor preparation, hydrogen addition, conversion to PDBQT, and removed nonreceptor molecules. |

Main docking parameters:

```yaml
software: AutoDock Vina v1.2.7
exhaustiveness: 4
num_modes: 5
energy_range: 5

3FAP:
  center: [-11.655, 24.756, 33.139]
  size: [37.5, 37.5, 45.0]

7PQV:
  center: [-5.223, 67.587, 35.180]
  size: [15.144, 16.921, 10.290]
```

Docking scores are used only as computational prioritization signals and should not be interpreted as experimental evidence of binding affinity or biological activity.

---

## 5. Analysis scripts

The `scripts/` folder contains or documents the scripts used to reproduce the main analyses. If a script is unavailable because it depends on local infrastructure, the corresponding processed output file is provided.

| Script | Purpose |
|---|---|
| `evaluate_generation_quality.py` | Calculates validity, uniqueness, novelty, IntDiv1, and IntDiv2. |
| `calc_descriptors.py` | Calculates QED, Synth/SA, MW, logP, TPSA, HBD, HBA, rotatable bonds, and Lipinski compliance. |
| `prepare_ligands.py` | Converts SMILES to 3D conformers and PDBQT ligand files. |
| `run_dual_vina_docking.py` | Runs 3FAP and 7PQV Vina docking under the unified downstream protocol. |
| `parse_vina_results.py` | Parses docking outputs and records best Vina scores. |
| `compute_dual_hit_rates.py` | Calculates simultaneous dual-hit rates at predefined thresholds. |
| `compute_topk_metrics.py` | Generates Top-k docking, QED, SA, and dual-hit summaries. |
| `compute_integrated_priority.py` | Calculates the integrated priority score and ranked head-candidate summaries. |
| `seed_similarity_scaffold_analysis.py` | Computes seed-fragment similarity and scaffold-overlap diagnostics. |
| `plip_hotspot_overlap_analysis.py` | Extracts PLIP interactions and calculates key-residue and reference-region recovery. |
| `summarize_high_exh_redocking.py` | Summarizes higher-exhaustiveness repeated-docking results. |
| `draw_figures.py` | Generates manuscript figures from processed tables. |

---

## 6. Reproducing main analyses

The package is organized around processed outputs. The following steps reproduce the reported numerical summaries from the provided data files.

### 6.1 Basic generation-quality metrics

Inputs:

```text
baseline_generation_quality/generated_molecules/*.csv
formal_sets/Prior-guided population set.csv
formal_sets/Contact-guided final set.csv
```

Command:

```bash
python scripts/evaluate_generation_quality.py \
  --input_dir baseline_generation_quality/generated_molecules \
  --reference_library data/250k_rndm_zinc_drugs_clean.smi \
  --out baseline_generation_quality/generation_quality_baselines_summary.csv
```

Expected output:

```text
baseline_generation_quality/generation_quality_baselines_summary.csv
```

This reproduces the validity, uniqueness, internal diversity, and novelty values reported for representative molecular-generation models and PSGS-Drug.

### 6.2 Descriptor calculation

Inputs:

```text
formal_sets/Prior-guided population set.csv
formal_sets/Contact-guided final set.csv
```

Command:

```bash
python scripts/calc_descriptors.py \
  --prior "formal_sets/Prior-guided population set.csv" \
  --contact "formal_sets/Contact-guided final set.csv" \
  --out admet_descriptors/descriptor_distribution_summary_prior_vs_contact.csv
```

Expected outputs include descriptor summaries for QED, Synth/SA, MW, logP, TPSA, HBD, HBA, rotatable bonds, and Lipinski compliance.

### 6.3 Dual-target docking and dual-hit analysis

Inputs:

```text
formal_sets/*.csv
receptors/3fap_receptor.pdbqt
receptors/7pqv_receptor.pdbqt
receptors/docking_boxes.yaml
```

Example workflow:

```bash
python scripts/prepare_ligands.py \
  --input "formal_sets/Contact-guided final set.csv" \
  --smiles_col smiles \
  --out_dir docking/ligands_contact

python scripts/run_dual_vina_docking.py \
  --ligand_dir docking/ligands_contact \
  --receptor_3fap receptors/3fap_receptor.pdbqt \
  --receptor_7pqv receptors/7pqv_receptor.pdbqt \
  --box_config receptors/docking_boxes.yaml \
  --out_dir docking/contact_vina_outputs

python scripts/parse_vina_results.py \
  --input_dir docking/contact_vina_outputs \
  --out formal_sets/Contact-guided final set_dual_docked.csv

python scripts/compute_dual_hit_rates.py \
  --input formal_sets/Contact-guided final set_dual_docked.csv \
  --thresholds -8.5 -9.0 -9.5 \
  --out topk/contact_dual_hit_summary.csv
```

These steps reproduce full-set docking and dual-hit summaries for formal and control result sets.

### 6.4 Top-k enrichment analysis

Inputs:

```text
formal_sets/Contact-guided final set.csv
topk/topk_by_docksum_prior_vs_contact.csv
topk/topk_by_integrated_priority_prior_vs_contact.csv
```

Command:

```bash
python scripts/compute_topk_metrics.py \
  --input "formal_sets/Contact-guided final set.csv" \
  --rank_by dock_sum \
  --k 10 20 50 100 \
  --out topk/contact_topk_by_docksum.csv
```

For integrated-priority ranking:

```bash
python scripts/compute_integrated_priority.py \
  --input "formal_sets/Contact-guided final set.csv" \
  --out topk/contact_topk_by_integrated_priority.csv
```

Expected outputs reproduce the Top-k dock_sum, dual-hit, QED, SA, and integrated-priority summaries.

### 6.5 Ablation and target-order control analysis

Inputs:

```text
controls_and_ablation/A0_A1_global_ablation_summary.csv
controls_and_ablation/A1_A2_topk_priority_summary.csv
controls_and_ablation/B1_B2_contact_full_summary.csv
```

These files reproduce:

- A0→A1 global prior-contribution metrics;
- A1→A2 Top-k prefix-contribution metrics;
- B1/B2 target-order sensitivity metrics.

### 6.6 Seed-fragment similarity and scaffold-overlap analysis

Inputs:

```text
scaffold_alerts/seed_fragment_scaffolds.csv
scaffold_alerts/seed_similarity_molecule_level.csv
formal_sets/Contact-guided final set.csv
```

Command:

```bash
python scripts/seed_similarity_scaffold_analysis.py \
  --input "formal_sets/Contact-guided final set.csv" \
  --seed_file scaffold_alerts/seed_fragment_scaffolds.csv \
  --out scaffold_alerts/seed_similarity_scaffold_summary.csv
```

This reproduces the seed-fragment similarity and scaffold-overlap analysis reported in Table S9.

### 6.7 ADMET and structural-alert summaries

Inputs:

```text
admet_descriptors/Contact-guided final-generation set-top100.csv
admet_descriptors/admetlab3_contact_top100_raw.csv
scaffold_alerts/Contact-guided final set-summary_pains.txt
scaffold_alerts/Contact-guided final set-bm_scaffold_summary.csv
scaffold_alerts/Contact-guided final set-ring-metric-summary.txt
```

These files reproduce the medicinal-chemistry, structural-alert, and ADMET summaries reported in Tables S4 and S5 and Figure 6.

### 6.8 PLIP hotspot-overlap analysis

Inputs:

```text
plip/top20_candidates_for_PLIP_interaction_analysis.csv
plip/plip_summary_top20.csv
plip/reference_hotspot_residues.csv
```

Command:

```bash
python scripts/plip_hotspot_overlap_analysis.py \
  --candidates plip/top20_candidates_for_PLIP_interaction_analysis.csv \
  --reference_hotspots plip/reference_hotspot_residues.csv \
  --out plip/plip_summary_top20.csv
```

This reproduces the candidate-level PLIP interaction and hotspot-overlap summary reported in Table S8.

### 6.9 Higher-exhaustiveness repeated docking

Inputs:

```text
repeated_docking/validation_high_exh_summary_table.csv
repeated_docking/validation_high_exh_per_molecule.csv
```

These files reproduce the higher-exhaustiveness redocking sensitivity analysis reported in Table S10.

---

## 7. Software environment

The analyses were performed with Python and standard cheminformatics/docking tools.

Recommended environment:

```yaml
name: psgs_drug_repro
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - pandas
  - numpy
  - scipy
  - scikit-learn
  - rdkit
  - matplotlib
  - openbabel
  - pip
  - pip:
      - useful-rdkit-utils
```

External software:

| Software | Role |
|---|---|
| AutoDock Vina v1.2.7 | Molecular docking |
| RDKit | SMILES parsing, canonicalization, descriptors, scaffolds, PAINS filtering |
| PLIP | Protein-ligand interaction profiling |
| ADMETlab 3.0 | ADMET endpoint prediction |
| useful_rdkit_utils | Ring-system frequency analysis |

If exact software versions differ from those listed above, rerun the processed-output scripts and compare against the supplied summary tables.

---

## 8. Notes on reproducibility

1. Docking scores are stochastic or implementation-sensitive to ligand preparation, protonation, receptor preparation, and search settings. The provided processed docking summaries should be treated as the reference outputs used for the manuscript.
2. The A2 full-set file is provided as a raw generation output and contains 1889 generated records with 1864 unique SMILES before final deduplication. This is explicitly indicated in the Supporting Information.
3. Pocket2Mol was used as a contextual external structure-conditioned comparator and was not intended as an exhaustive benchmark suite for all structure-based generative models.
4. JT-VAE, MARS, RationaleRL, REINVENT2.0, and VeGA were used only for basic generation-quality contextualization and were not used as dual-target structure-conditioned baselines.
5. ADMETlab predictions are computational risk diagnostics and are not experimental safety validation.
6. PLIP-derived residue contacts and hotspot-overlap summaries are used for interpretability and do not establish experimental binding affinity or target engagement.

---

## 9. Citation

If this package is used, please cite the associated manuscript:

```text
PSGS-Drug: Generation-Time Dual-Pocket Guidance for Non-Symmetric Dual-Target Molecular Design.
```

---

## 10. Contact

For questions about the data package or to request additional implementation details, contact the corresponding author listed in the manuscript.
