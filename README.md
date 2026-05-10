# PSGS-Drug

**PSGS-Drug: Generation-Time Dual-Pocket Guidance for Non-Symmetric Dual-Target Molecular Design**

This repository provides the code, configuration files, receptor/pocket-prior inputs, processed result tables, and analysis scripts associated with the PSGS-Drug manuscript. PSGS-Drug is a generation-time dual-pocket molecular design workflow for a non-symmetric 3FAP/7PQV dual-target task. The repository is intended to support reproducibility of the reported computational results, including the formal prior-guided population set, contact-guided final-generation set, ablation studies, single-pocket target-order controls, Top-k analysis, ADMET/medicinal-chemistry summaries, PLIP hotspot-overlap analysis, and higher-exhaustiveness redocking summaries.

Docking scores, PLIP contacts, and ADMET predictions in this repository are computational prioritization signals only. They should not be interpreted as experimental evidence of binding affinity, target engagement, biological activity, pharmacokinetics, or safety.

---

## Repository structure

```text
PSGS-Drug/
├── code/
│   ├── README.md
│   ├── generation/
│   │   ├── prior_client.py
│   │   ├── sbmolgen_A0A1A2.py
│   │   ├── sbmolgen_B1B2.py
│   │   ├── sbmolgen_contactfrag.py
│   │   └── sbmolgen_prior.py
│   └── token_prior/
│       ├── tokenmol_policy.py
│       └── tokenmol_prior_server.py
│
├── configs/
│   ├── README_configs.md
│   ├── setting_prior.yaml
│   ├── setting_prior_contactfrag.yaml
│   ├── setting_ablation_A0_noprior_noseed.yaml
│   ├── setting_ablation_A1_prior_only.yaml
│   ├── setting_ablation_A2_prior_plus_seed.yaml
│   ├── setting_baseline_B1_3fap_only.yaml
│   └── setting_baseline_B2_7pqv_only.yaml
│
├── data/
│   ├── README_processed_results.md
│   ├── formal_sets/
│   │   ├── Prior-guided population set.csv
│   │   └── Contact-guided final set.csv
│   └── controls_and_ablation/
│       ├── result_3fap-7pqv_A0_noprior_noseed.csv
│       ├── result_3fap-7pqv_A1_prior_only.csv
│       ├── result_3fap-7pqv_A2_prior_plus_seed.raw.csv
│       ├── result_B1_3fap_only_postscreen_7pqv.csv
│       └── result_B2_7pqv_only_postscreen_3fap.csv
│
├── environment/
│   ├── README_environment.md
│   ├── environment_psgs_drug.yml
│   └── environment_psgs_token.yml
│
├── receptors/
│   ├── README_receptors.md
│   ├── 3fap.pkl
│   ├── 7pqv.pkl
│   ├── 3fappro.pdbqt
│   ├── 7pqvpro_clean_fix2.pdbqt
│   ├── 3fapconfig.txt
│   └── 7pqvconfig.txt
│
└── scripts/
    ├── README_scripts.md
    └── analysis/
        ├── admet_prior_profile.py
        ├── analyze_a1_a2_topk.py
        ├── analyze_b1_b2_contact.py
        ├── contact_guided_topk_analysis.py
        ├── higher_exhaustiveness_redocking.py
        ├── plip_hotspot_overlap_analysis.py
        ├── prior_contact_summary.py
        ├── run_plip_for_top20_candidates.py
        ├── seed_similarity_scaffold_analysis.py
        └── select_balanced_representatives.py
```

**Before uploading to GitHub:** remove IDE-specific files such as `.idea/`, `__pycache__/`, and temporary output folders.

---

## Main components

### 1. Formal result sets

The two formal result files are in `data/formal_sets/`:

| File | Description |
|---|---|
| `Prior-guided population set.csv` | Formal protein-pocket-prior-guided population set used for population-level directionality analysis. |
| `Contact-guided final set.csv` | Formal contact-guided final-generation set used for full-set docking, Top-k enrichment, candidate prioritization, ADMET/descriptor analysis, and PLIP interpretation. |

### 2. Ablation and target-order controls

The ablation and control outputs are in `data/controls_and_ablation/`:

| File | Description |
|---|---|
| `result_3fap-7pqv_A0_noprior_noseed.csv` | A0 ablation: no pocket prior and no contact-derived seed. |
| `result_3fap-7pqv_A1_prior_only.csv` | A1 ablation: dual-pocket prior without contact-derived seed. |
| `result_3fap-7pqv_A2_prior_plus_seed.raw.csv` | A2 raw ablation output: dual-pocket prior with contact-derived seed fragments. |
| `result_B1_3fap_only_postscreen_7pqv.csv` | B1 target-order control: 3FAP-first generation followed by 7PQV post-screening. |
| `result_B2_7pqv_only_postscreen_3fap.csv` | B2 target-order control: 7PQV-first generation followed by 3FAP post-screening. |

The A2 file is intentionally provided as a raw output. It contains generated records before final deduplication and should be interpreted as an ablation record rather than as the formal final result set.

### 3. Receptor and pocket-prior files

The `receptors/` folder contains:

- prepared 3FAP and 7PQV receptor PDBQT files;
- docking box configuration records;
- pocket representation files (`3fap.pkl` and `7pqv.pkl`) used by the pocket-prior module.

### 4. Generation and pocket-prior code

The `code/` folder contains PSGS-Drug-specific generation wrappers and Token-Mol prior-service wrappers. These files depend on third-party model code and local generation infrastructure from the original SFG-Drug and Token-Mol implementations.

This repository does not redistribute third-party model repositories or pretrained weights. Users should obtain third-party dependencies from their original repositories or publications.

### 5. Analysis scripts

The `scripts/analysis/` folder contains cleaned scripts for reproducing processed numerical summaries from the supplied CSV files. These scripts use command-line arguments and do not require local absolute paths.

---

## Environment setup

This project uses two conda environments because the molecular-generation workflow and the pocket-prior inference module rely on different dependencies.

### PSGS-Drug environment

Used for processed-result analysis, RDKit-based molecule processing, descriptor calculation, docking-result parsing, Top-k analysis, ADMET/PLIP result processing, and table/figure summary generation.

```bash
conda env create -f environment/environment_psgs_drug.yml
conda activate psgs-drug
```

### PSGS-token environment

Used for Token-Mol pocket encoding, target-specific pocket-prior inference, and local prior-server execution.

```bash
conda env create -f environment/environment_psgs_token.yml
conda activate psgs-token
```

See `environment/README_environment.md` for details.

---

## External tools and third-party dependencies

Some parts of the workflow require external tools or third-party code:

| Tool or codebase | Role |
|---|---|
| AutoDock Vina v1.2.7 | Molecular docking |
| RDKit | SMILES parsing, canonicalization, descriptors, scaffolds, and filters |
| OpenBabel | File-format conversion for docking/PLIP workflows |
| PLIP | Protein-ligand interaction profiling |
| ADMETlab 3.0 | ADMET endpoint prediction |
| Token-Mol | Pocket encoding and protein-aware token prior |
| SFG-Drug original code | MCTS-GRU molecular generation backbone |
| Pocket2Mol | Contextual external structure-conditioned comparator |
| JT-VAE, MARS, RationaleRL, REINVENT2.0, VeGA | Basic generation-quality contextualization |

Third-party model implementations and pretrained weights are not redistributed in this repository. Users should follow the installation instructions and licenses of the original projects.

---

## Example analysis commands

### Formal prior-guided versus contact-guided comparison

```bash
python scripts/analysis/prior_contact_summary.py \
  --prior_csv "data/formal_sets/Prior-guided population set.csv" \
  --contact_csv "data/formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/prior_contact"
```

### Contact-guided Top-k and integrated-priority analysis

```bash
python scripts/analysis/contact_guided_topk_analysis.py \
  --input_csv "data/formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/contact_topk"
```

### A1/A2 Top-k ablation analysis

```bash
python scripts/analysis/analyze_a1_a2_topk.py \
  --a1_csv "data/controls_and_ablation/result_3fap-7pqv_A1_prior_only.csv" \
  --a2_csv "data/controls_and_ablation/result_3fap-7pqv_A2_prior_plus_seed.raw.csv" \
  --out_dir "outputs/a1_a2_topk"
```

### B1/B2 target-order sensitivity analysis

```bash
python scripts/analysis/analyze_b1_b2_contact.py \
  --b1_csv "data/controls_and_ablation/result_B1_3fap_only_postscreen_7pqv.csv" \
  --b2_csv "data/controls_and_ablation/result_B2_7pqv_only_postscreen_3fap.csv" \
  --contact_csv "data/formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/b1_b2_contact"
```

### Seed-fragment similarity and scaffold-overlap analysis

```bash
python scripts/analysis/seed_similarity_scaffold_analysis.py \
  --input_csv "data/formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/seed_similarity"
```

### PLIP hotspot-overlap table preparation

```bash
python scripts/analysis/plip_hotspot_overlap_analysis.py \
  --plip_summary_csv "plip/plip_summary_top20.csv" \
  --out_csv "outputs/table_s8_plip_hotspot_overlap.csv"
```

### Higher-exhaustiveness repeated-docking summary

```bash
python scripts/analysis/higher_exhaustiveness_redocking.py \
  --input_csv "repeated_docking/validation_repeated_docking_high_exh_raw.csv" \
  --out_summary "outputs/validation_high_exh_summary_table.csv" \
  --out_per_molecule "outputs/validation_high_exh_per_molecule.csv"
```

---

## Data and software availability

This repository provides:

1. PSGS-Drug-specific generation wrappers and prior-service interfaces;
2. configuration files for formal runs, ablations, and target-order controls;
3. receptor files, pocket representation files, and docking-box records;
4. processed molecule-level result sets used in the manuscript;
5. analysis scripts used to reproduce processed numerical summaries.

The complete upstream generation workflow may require third-party model code, pretrained weights, and local infrastructure. For such components, this repository provides PSGS-Drug-specific configuration files and processed outputs sufficient to verify the reported numerical summaries.

---

## Reproducibility notes

1. Docking scores are sensitive to ligand preparation, receptor preparation, protonation, search settings, and software version. The processed docking summaries supplied here are the reference outputs used in the manuscript.
2. PLIP-derived residue contacts are used for interaction-pattern interpretation and do not establish experimental binding.
3. ADMETlab predictions are computational risk diagnostics and do not establish safety.
4. The baseline generation-quality models are used only for basic molecular-generation quality contextualization and are not dual-target structure-conditioned baselines.
5. Some generation components depend on external model resources. When exact reruns are not possible, processed outputs and configuration files are provided for verification of the manuscript results.

---

## License

Add the appropriate license for your repository before public release. If third-party code is used or referenced, follow the licenses of the original repositories.

---

## Citation

If this repository is used, please cite the associated manuscript:

```text
PSGS-Drug: Generation-Time Dual-Pocket Guidance for Non-Symmetric Dual-Target Molecular Design.
```

---

## Contact

For questions about this repository, contact the corresponding author listed in the manuscript.
