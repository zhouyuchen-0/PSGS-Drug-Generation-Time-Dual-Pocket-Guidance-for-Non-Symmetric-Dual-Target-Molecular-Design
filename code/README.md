# PSGS-Drug generation and pocket-prior code

This folder contains cleaned code associated with the PSGS-Drug molecular generation workflow and the Token-Mol pocket-prior service.

## Scope

These scripts document the implementation used for:

- formal prior-guided population generation;
- formal contact-guided generation from contact-derived seed prefixes;
- A0/A1/A2 controlled ablation runs;
- B1/B2 single-pocket target-order controls;
- Token-Mol pocket-prior inference through a local FastAPI service.

Some modules still depend on the original SFG-Drug utilities and Token-Mol model files. Third-party code and pretrained weights are not redistributed here. Users should obtain external dependencies from the original repositories or publications.

## Directory structure

```text
generation/
├── sbmolgen_prior.py
├── sbmolgen_contactfrag.py
├── sbmolgen_A0A1A2.py
├── sbmolgen_B1B2.py
└── prior_client.py

token_prior/
├── tokenmol_policy.py
└── tokenmol_prior_server.py
```

## Configuration

All generation scripts are intended to be run with a YAML configuration file:

```bash
python generation/sbmolgen_contactfrag.py --cfg configs/setting_prior_contactfrag.yaml
```

If `--cfg` is not provided, scripts look for the `SFGDRUG_CFG_PATH` environment variable and then fallback to repository-relative defaults.

## Pocket-prior server

Start the Token-Mol prior server in the `psgs-token` environment:

```bash
uvicorn token_prior.tokenmol_prior_server:app --host 127.0.0.1 --port 26974
```

The model paths can be configured with environment variables:

```bash
export TOKENMOL_MODEL_CKPT=./Trained_model/pocket_generation.pt
export TOKENMOL_VOCAB_PATH=./data/torsion_version/torsion_voc_pocket.csv
export TOKENMOL_PRETRAIN_DIR=./Pretrained_model
```

## Notes for public release

- Local absolute paths were replaced with repository-relative paths where possible.
- Comments and docstrings were converted to English.
- Legacy target names are not used.
- The scripts require original SFG-Drug utility modules and Token-Mol model resources.
- For manuscript reproducibility, processed outputs and configuration files are provided in the Supporting Data package.
