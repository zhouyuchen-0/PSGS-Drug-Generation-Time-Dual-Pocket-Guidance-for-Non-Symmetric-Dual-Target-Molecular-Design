# PSGS-Drug

**Generation-Time Dual-Pocket Guidance for Non-Symmetric Dual-Target Molecular Design**

This repository contains the code, configuration files, receptor files, model components, generated-molecule result sets, and analysis scripts used in the PSGS-Drug study.

PSGS-Drug is a generation-time dual-pocket molecular generation framework. It includes:

- an MCTS-GRU molecular generation module;
- a PSGS pocket-prior service for target-conditioned token-level prior inference;
- equal-weight dual-pocket prior fusion;
- contact-guided prefix initialization;
- downstream analysis scripts for docking, descriptors, ADMET, PLIP, redocking, and ablation results.

All scripts are intended to use project-relative paths. Local path editing is not required if the repository structure is kept unchanged.

## Model assets

The GitHub repository contains the source code, configuration files, receptor files, result data, and analysis scripts. Large model files are provided separately through the GitHub Release associated with this repository.

Please download the release asset:

```text
psgs_model_assets.zip
```

Unzip it into the repository root. After extraction, the repository root should contain:

```text
configs/
data/
environment/
psgs_model/
receptors/
scripts/
README.md
run_prior_server.py
run_generation_demo.py
```


## 1. Main folders

| Folder | Description |
|---|---|
| `configs/` | YAML configuration files for formal generation, contact-guided generation, ablation experiments, and single-pocket controls. |
| `data/` | Input molecular library, lightweight test subset, generated-molecule result sets, and evaluation data. |
| `environment/` | Conda environment files and dependency records. |
| `psgs_model/` | PSGS model-related code and files, including the generation module and the pocket-prior service. |
| `receptors/` | Receptor files, docking configuration files, and pocket feature files for 3FAP and 7PQV. |
| `scripts/` | Analysis scripts for reproducing reported result metrics from supplied CSV files. |

---

## 2. Environment setup

PSGS-Drug uses two conda environments:

- `psgs-prior`: environment for the PSGS pocket-prior service.
- `psgs-core`: environment for PSGS-Drug molecular generation, docking-related processing, and downstream result analysis.

The YAML files are the recommended installation method. The requirement text files are provided only as dependency records or fallback references.

### 2.1 Create conda environments

```bash
conda env create -f environment/environment_psgs_prior.yml
conda env create -f environment/environment_psgs_core.yml
```

### 2.2 Verify installation

```bash
conda activate psgs-prior
python -c "import torch, transformers, fastapi; print('psgs-prior OK')"
```

```bash
conda activate psgs-core
python -c "import rdkit, pandas, numpy; print('psgs-core OK')"
```

---

## 3. Data files

### 3.1 Input molecular library

```text
data/250k_rndm_zinc_drugs_clean.smi
```

This file is the input molecular library used by the generation scripts.

### 3.2 Lightweight test subset

```text
data/test_subset_1000.smi
```

This file contains 1000 SMILES sampled from the input molecular library. It is provided for quick smoke testing of the generation workflow.

### 3.3 Reported result files

```text
data/formal_sets/
data/controls_and_ablation/
```

These folders contain generated-molecule result sets and evaluation data used for reported-result analysis.

### 3.4 Additional structural-transfer receptor assets: CDK6 (5L2I)–ERα (3ERT)

The repository additionally archives the receptor assets used for the 5L2I–3ERT structural-transfer validation on a structurally heterogeneous kinase–nuclear-receptor pair. These files document the target-specific pocket-prior inputs and the fixed AutoDock Vina receptor/grid definitions used for this additional validation setting.

Place the following files in the repository at the stated paths:

```text
receptors/
├── 5l2i.pkl
├── 3ert.pkl
├── 5l2ipro.pdbqt
├── 3ertpro.pdbqt
├── 5l2iconfig.txt
└── 3ertconfig.txt
```

| File | Role in the structural-transfer validation |
|---|---|
| `receptors/5l2i.pkl` | Precomputed pocket-prior feature input for the cyclin-dependent kinase 6 (CDK6) ATP-binding-pocket context represented by Protein Data Bank (PDB) entry 5L2I. |
| `receptors/3ert.pkl` | Precomputed pocket-prior feature input for the estrogen receptor alpha (ERα) ligand-binding-domain context represented by PDB entry 3ERT. |
| `receptors/5l2ipro.pdbqt` | Prepared CDK6 receptor used for fixed-protocol AutoDock Vina docking. |
| `receptors/3ertpro.pdbqt` | Prepared ERα receptor used for fixed-protocol AutoDock Vina docking. |
| `receptors/5l2iconfig.txt` | Fixed AutoDock Vina grid configuration for 5L2I. |
| `receptors/3ertconfig.txt` | Fixed AutoDock Vina grid configuration for 3ERT. |

The 5L2I–3ERT receptor PDBQT files, docking grids, and Vina settings are fixed receptor-side assets for the structural-transfer validation. They are provided to document the additional target contexts and the corresponding fixed docking protocol; they should not be modified after release.

> **Repository scope.** This release provides the receptor-side structural resources used for the 5L2I–3ERT structural-transfer validation, including the precomputed pocket-prior feature inputs, prepared receptor files, and fixed AutoDock Vina grid definitions. Together, these assets document the structural conditioning and fixed receptor/docking protocol used for the additional kinase–nuclear-receptor validation.

---

## 4. Running PSGS-Drug generation

Generation runs that use a pocket prior require the PSGS pocket-prior service to be running.

### 4.1 Start the PSGS pocket-prior service

Open Terminal 1:

```bash
conda activate psgs-prior
python run_prior_server.py
```

The service should start at:

```text
http://127.0.0.1:26974/prior
```

A successful startup indicates that the PSGS pocket-prior service has loaded the prior model, vocabulary file, pretrained model directory, and receptor pocket feature files.

### 4.2 Run a generation configuration

Open Terminal 2:

```bash
conda activate psgs-core
python run_generation_demo.py
```

By default, the runner executes the contact-guided final-generation configuration:

```text
configs/setting_prior_contactfrag.yaml
```

The runner supports multiple generation modes through `--mode`.

```bash
python run_generation_demo.py --mode prior
python run_generation_demo.py --mode contact
python run_generation_demo.py --mode a0
python run_generation_demo.py --mode a1
python run_generation_demo.py --mode a2
python run_generation_demo.py --mode b1
python run_generation_demo.py --mode b2
```

The mapping between mode names and configuration files is:

| Mode | Configuration file | Generation script |
|---|---|---|
| `prior` | `configs/setting_prior.yaml` | `psgs_model/generation/sbmolgen_prior.py` |
| `contact` | `configs/setting_prior_contactfrag.yaml` | `psgs_model/generation/sbmolgen_contactfrag.py` |
| `a0` | `configs/setting_ablation_A0_noprior_noseed.yaml` | `psgs_model/generation/sbmolgen_A0A1A2.py` |
| `a1` | `configs/setting_ablation_A1_prior_only.yaml` | `psgs_model/generation/sbmolgen_A0A1A2.py` |
| `a2` | `configs/setting_ablation_A2_prior_plus_seed.yaml` | `psgs_model/generation/sbmolgen_A0A1A2.py` |
| `b1` | `configs/setting_baseline_B1_3fap_only.yaml` | `psgs_model/generation/sbmolgen_B1B2.py` |
| `b2` | `configs/setting_baseline_B2_7pqv_only.yaml` | `psgs_model/generation/sbmolgen_B1B2.py` |

A custom YAML configuration file can also be supplied:

```bash
python run_generation_demo.py --cfg configs/setting_prior_contactfrag.yaml
```

When `--cfg` is supplied, the runner infers the generation mode from the configuration filename. If needed, a mode can also be specified explicitly:

```bash
python run_generation_demo.py --mode contact --cfg configs/setting_prior_contactfrag.yaml
```

To check resolved paths without launching generation:

```bash
python run_generation_demo.py --mode contact --dry-run
```

Notes:

- `a0` does not require the pocket-prior service.
- `prior`, `contact`, `a1`, `a2`, `b1`, and `b2` require the PSGS pocket-prior service.
- All paths are resolved relative to the project root.

---

## 5. Result analysis from supplied CSV files

The `scripts/` folder contains analysis scripts for reproducing reported result metrics from the supplied CSV files.

**B1/B2 terminology.** `analyze_b1_b2_contact.py` evaluates B1 and B2 as
single-pocket **conditioning-source asymmetry controls** under the shared
retrospective two-context docking workflow. B1 denotes **3FAP-only**
conditioning during generation, and B2 denotes **7PQV-only** conditioning
during generation. These controls do **not** encode or test a sequential
first-target/second-target generation order.

### 5.1 Structural-transfer asset scope

For the 5L2I–3ERT validation, the repository provides the receptor-side pocket-prior features, prepared receptor PDBQT files, and fixed AutoDock Vina grid configurations listed in Section 3.4. These files document the additional CDK6–ERα structural contexts used for the reported structural-transfer analysis.

For interpretation of the 5L2I–3ERT validation, the manuscript and Supporting Information, Table S16, report the aggregate five-seed paired endpoints, 95% confidence intervals, and Holm-adjusted inference. The receptor-side assets provided here identify the exact structural contexts and fixed docking definitions used for that reported validation.

Available scripts include:

```text
scripts/prior_contact_summary.py
scripts/contact_guided_topk_analysis.py
scripts/analyze_a1_a2_topk.py
scripts/analyze_b1_b2_contact.py
scripts/admet_prior_profile.py
scripts/seed_similarity_scaffold_analysis.py
scripts/higher_exhaustiveness_redocking.py
scripts/plip_hotspot_overlap_analysis.py
scripts/run_plip_for_top20_candidates.py
scripts/select_balanced_representatives.py
```

Example commands:

```bash
conda activate psgs-core

python scripts/prior_contact_summary.py
python scripts/contact_guided_topk_analysis.py
python scripts/analyze_a1_a2_topk.py
python scripts/analyze_b1_b2_contact.py --b1_csv data/controls_and_ablation/result_B1_3fap_only_postscreen_7pqv.csv --b2_csv data/controls_and_ablation/result_B2_7pqv_only_postscreen_3fap.csv --contact_csv "data/formal_sets/Contact-guided final set.csv" --out_dir outputs/b1_b2_conditioning_source_asymmetry --topk 50
python scripts/admet_prior_profile.py
python scripts/seed_similarity_scaffold_analysis.py
python scripts/higher_exhaustiveness_redocking.py
python scripts/plip_hotspot_overlap_analysis.py
```

The exact input and output paths are defined inside each script or through project-relative configuration paths.

---

## 6. Configuration files

Main configuration files are located in:

```text
configs/
```

Important files include:

```text
configs/setting_prior.yaml
configs/setting_prior_contactfrag.yaml
configs/setting_ablation_A0_noprior_noseed.yaml
configs/setting_ablation_A1_prior_only.yaml
configs/setting_ablation_A2_prior_plus_seed.yaml
configs/setting_baseline_B1_3fap_only.yaml
configs/setting_baseline_B2_7pqv_only.yaml
```

All paths in the configuration files are project-relative.

---

## 7. Receptor and pocket-prior files

The receptor and pocket-prior files are located in:

```text
receptors/
```

Key files:

```text
receptors/3fap.pkl
receptors/7pqv.pkl
receptors/3fappro.pdbqt
receptors/7pqvpro_clean_fix2.pdbqt
receptors/3fapconfig.txt
receptors/7pqvconfig.txt
receptors/5l2i.pkl
receptors/3ert.pkl
receptors/5l2ipro.pdbqt
receptors/3ertpro.pdbqt
receptors/5l2iconfig.txt
receptors/3ertconfig.txt
```

- `3fap.pkl` and `7pqv.pkl` are pocket feature files used by the PSGS pocket-prior service for the primary 7PQV–3FAP study setting.
- `3fappro.pdbqt` and `7pqvpro_clean_fix2.pdbqt` are receptor files used for primary-study AutoDock Vina docking.
- `3fapconfig.txt` and `7pqvconfig.txt` define the primary-study docking-box settings.
- `5l2i.pkl`, `3ert.pkl`, `5l2ipro.pdbqt`, `3ertpro.pdbqt`, `5l2iconfig.txt`, and `3ertconfig.txt` are the receptor-side pocket-prior and fixed-docking assets for the additional 5L2I–3ERT structural-transfer validation.

---

## 8. PSGS model package

Model-related files are placed under:

```text
psgs_model/
```

### 8.1 Generation module

```text
psgs_model/generation/
```

This module contains the MCTS-GRU molecular generation code, RNN model files, and generation utilities.

Key files:

```text
psgs_model/generation/RNN-model/model.json
psgs_model/generation/RNN-model/model.h5
psgs_model/generation/sbmolgen_prior.py
psgs_model/generation/sbmolgen_contactfrag.py
psgs_model/generation/sbmolgen_A0A1A2.py
psgs_model/generation/sbmolgen_B1B2.py
psgs_model/generation/prior_client.py
```

### 8.2 Pocket-prior module

```text
psgs_model/prior/
```

This module contains the PSGS pocket-prior service and model files.

Key files:

```text
psgs_model/prior/tokenmol_prior_server.py
psgs_model/prior/tokenmol_policy.py
psgs_model/prior/ada_model.py
psgs_model/prior/bert_tokenizer.py
psgs_model/prior/pocket_fine_tuning_rmse.py
psgs_model/prior/pocket_generation.pt
psgs_model/prior/torsion_voc_pocket.csv
psgs_model/prior/Pretrained_model/
```

---

## 9. Reproducibility scope

This repository supports two levels of reproducibility.

### 9.1 Result-level analysis

This mode reproduces reported result metrics from supplied CSV files. It does not require rerunning molecular generation.

```bash
conda activate psgs-core
python scripts/prior_contact_summary.py
python scripts/contact_guided_topk_analysis.py
python scripts/analyze_a1_a2_topk.py
python scripts/analyze_b1_b2_contact.py --b1_csv data/controls_and_ablation/result_B1_3fap_only_postscreen_7pqv.csv --b2_csv data/controls_and_ablation/result_B2_7pqv_only_postscreen_3fap.csv --contact_csv "data/formal_sets/Contact-guided final set.csv" --out_dir outputs/b1_b2_conditioning_source_asymmetry --topk 50
```

### 9.2 Generation workflow execution

This mode starts the PSGS pocket-prior service and runs the PSGS-Drug generation scripts with the provided receptor files, model files, configuration files, and input molecular library.

Terminal 1:

```bash
conda activate psgs-prior
python run_prior_server.py
```

Terminal 2:

```bash
conda activate psgs-core
python run_generation_demo.py --mode contact
```

---

Also ensure that AutoDock Vina and OpenBabel are installed in the `psgs-core` environment.

---

## 10. Notes

The repository uses project-relative paths. No local path editing is required if the repository structure is kept unchanged.

For generation workflow testing, first start the PSGS pocket-prior service in `psgs-prior`, then run the generation runner in `psgs-core`.

---

## 11. Citation

Please cite the associated manuscript when using this repository.

---

## 12. Release checklist for the 5L2I–3ERT structural-transfer receptor assets

Before claiming that the 5L2I–3ERT receptor assets are publicly available, confirm that the repository or its associated GitHub Release contains all of the following:

- [ ] `receptors/5l2i.pkl` and `receptors/3ert.pkl`;
- [ ] `receptors/5l2ipro.pdbqt` and `receptors/3ertpro.pdbqt`;
- [ ] `receptors/5l2iconfig.txt` and `receptors/3ertconfig.txt`;
- [ ] an explicit README statement that these files document the receptor-side structural inputs and fixed docking grids for the 5L2I–3ERT structural-transfer validation.

The primary 7PQV–3FAP assets and the additional 5L2I–3ERT receptor assets should be versioned with the corresponding manuscript revision to avoid ambiguity about the structural and docking protocol used for each reported analysis.
