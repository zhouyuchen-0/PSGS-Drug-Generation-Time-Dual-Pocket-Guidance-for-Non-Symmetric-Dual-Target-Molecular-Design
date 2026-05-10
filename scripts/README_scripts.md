# PSGS-Drug analysis scripts

This folder contains cleaned, upload-ready scripts associated with the PSGS-Drug manuscript.

## Notes

- Hard-coded Windows paths have been replaced by command-line arguments.
- Comments and docstrings are written in English.
- Legacy target names have been removed.
- These scripts are intended to reproduce processed summaries from supplied CSV files.
- Some original generation components may depend on local infrastructure or third-party repositories; in those cases, configuration files and processed outputs are supplied for verification.

## Typical usage

```bash
python analysis/prior_contact_summary.py \
  --prior_csv "formal_sets/Prior-guided population set.csv" \
  --contact_csv "formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/prior_contact"

python analysis/contact_guided_topk_analysis.py \
  --input_csv "formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/contact_topk"

python analysis/analyze_b1_b2_contact.py \
  --b1_csv "controls_and_ablation/result_B1_3fap_only_postscreen_7pqv.csv" \
  --b2_csv "controls_and_ablation/result_B2_7pqv_only_postscreen_3fap.csv" \
  --contact_csv "formal_sets/Contact-guided final set.csv" \
  --out_dir "outputs/b1_b2"

python analysis/plip_hotspot_overlap_analysis.py \
  --plip_summary_csv "plip/plip_summary_top20.csv" \
  --out_csv "outputs/table_s8_plip_hotspot_overlap.csv"
```
