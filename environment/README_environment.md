# Environment setup

This project uses two conda environments because the PSGS-Drug molecular-generation workflow and the Token-Mol pocket-prior inference module rely on different dependencies.

## 1. PSGS-Drug environment

The `psgs-drug` environment is used for the main PSGS-Drug analysis and generation workflow, including MCTS-GRU molecular generation, RDKit-based molecule processing, docking-result parsing, descriptor calculation, Top-k analysis, scaffold/ring-system analysis, ADMET result processing, PLIP result summarization, and figure/table generation.

Create and activate the environment:

```bash
conda env create -f environment_psgs_drug.yml
conda activate psgs-drug
```

## 2. PSGS-token environment

The `psgs-token` environment is used for Token-Mol pocket encoding, target-specific pocket-prior inference, and generation of the 3FAP/7PQV pocket-prior outputs used by PSGS-Drug.

Create and activate the environment:

```bash
conda env create -f environment_psgs_token.yml
conda activate psgs-token
```

## External tools

Some analyses also require external tools that are not fully managed by these conda files:

- AutoDock Vina v1.2.7 for molecular docking.
- PLIP for protein-ligand interaction profiling.
- ADMETlab 3.0 for ADMET endpoint prediction.
- Third-party model code or checkpoints for Token-Mol, SFG-Drug, Pocket2Mol, and other contextual baseline methods should be obtained from their original repositories or publications.

## Notes

The two environments are separated to avoid dependency conflicts between the MCTS-GRU molecular generation workflow and the Token-Mol pocket-prior module.

If users only want to reproduce processed tables and figures from the supplied CSV files, the `psgs-drug` environment is usually sufficient.

If users want to rerun pocket encoding or target-specific pocket-prior inference, the `psgs-token` environment is required.

These files are intended as practical reproducibility environments rather than strict lockfiles. If exact package versions are required, users should record and share the fully resolved environment from the original execution platform.
