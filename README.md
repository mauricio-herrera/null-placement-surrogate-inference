# Where the null is imposed matters
## Reproducibility package

**Manuscript:** *Where the null is imposed matters: surrogate inference after thresholding and temporal aggregation*

**Author:** Mauricio Herrera-Marín  
Faculty of Engineering, Universidad del Desarrollo, Santiago, Chile  
ORCID: 0000-0002-9604-3077  
Contact: mherrera@udd.cl
Zenodo v1.0.0 DOI: 10.5281/zenodo.21873328

## Scope

This repository contains the fixed computational protocols, synthetic-analysis scripts,
row-level and summary outputs, independent implementation audit, final manuscript
figure data and plotting scripts, and the public-CMIP6 reconstruction workflow used
for the Amazon illustration.

Raw CMIP6 files are not redistributed.

The confirmatory computational program was closed before manuscript finalization:
no additional lambda/kappa grid points or post-result tuning are part of the reported evidence.

## Repository structure

- `scripts/synthetic/` — mechanism audit, confirmatory experiment, sensitivity experiment,
  probability-equalization experiment, and independent implementation audit.
- `scripts/figures/` — final scripts used to generate the four main manuscript figures.
- `scripts/amazon/` — public Pangeo CMIP6 Amazon reconstruction and provenance audit.
- `protocols/` — protocols fixed before the corresponding scientific runs.
- `freeze_hashes/` — SHA256 records of the frozen confirmatory packages.
- `results/stage3/` — 13,500-trajectory confirmatory outputs.
- `results/stage4/` — mechanism-specific eight-unit sensitivity outputs.
- `results/stage5/` — probability-equalization and estimator-diagnostic outputs.
- `results/stage5C/` — independent implementation audit.
- `data/figure_data/` — CSV tables used by the final figure scripts.
- `figures/` — final manuscript figures.
- `supplementary/` — Supplemental Material source.
- `STAGE5C_SELECTION_V060.csv` — frozen independent-audit selection.

## Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Reproducing the final figures

```bash
python scripts/figures/make_figure1.py
python scripts/figures/make_figures234.py
```

The plotting scripts read the CSV files in `data/figure_data/`. If you run them from
another directory, preserve the repository structure or adjust the relative data path
documented in the scripts.

## Re-running the synthetic analyses

These calculations are computationally intensive. Use smoke tests first where supported.

```bash
python scripts/synthetic/02_stage3_confirmatory_v040.py --smoke --jobs 1 --out smoke_stage3
python scripts/synthetic/03_stage4_power_v050.py --smoke --jobs 1 --out smoke_stage4
python scripts/synthetic/04_stage5AB_mitigation_v060.py --smoke --jobs 1 --out smoke_stage5
```

Full runs:

```bash
python scripts/synthetic/02_stage3_confirmatory_v040.py --jobs 4 --out stage3_run
python scripts/synthetic/03_stage4_power_v050.py --jobs 4 --out stage4_run
python scripts/synthetic/04_stage5AB_mitigation_v060.py --jobs 4 --out stage5_run
```

Independent closure audit:

```bash
python scripts/synthetic/05_stage5C_independent_audit_v060.py \
  --stage3-dir results/stage3 \
  --selection STAGE5C_SELECTION_V060.csv \
  --jobs 4 \
  --out stage5C_run
```

## Amazon / CMIP6 reconstruction

Raw CMIP6 data are not included. To reconstruct the climate illustration from the
public Pangeo CMIP6 catalog:

```bash
python scripts/amazon/extract_monthly_amazon_pr_gate_v013.py
```

See `scripts/amazon/README.md` for the model list, experiments, threshold, regional
definition, asset-selection rules, and provenance-audit script.

## Headline archived quantities

Confirmatory resolution experiment:
- 13,500 trajectories.
- Index-resolution rejections: 1,697.
- Native-resolution rejections: 629.
- Index-only discordances: 1,084.
- Native-only discordances: 16.

Seasonal concentration:
- index-only discordance: 4.41% at lambda = 0;
- index-only discordance: 16.22% at lambda = 4.

Probability equalization:
- concentrated profile: 15.85%;
- oracle uniform profile: 3.93%;
- independent 30-year percentile profile: 3.48%.

Mechanism-specific sensitivity:
- strongest prespecified history-feedback alternative:
  30.7% detection rate (95% Wilson CI 23.8%–38.5%).

## Integrity

`FILE_MANIFEST.csv` lists all archived scientific files and sizes, excluding the two integrity files themselves.  
`SHA256SUMS.txt` contains SHA256 hashes for the archived repository contents and for `FILE_MANIFEST.csv`.

After intentionally changing any archived file, rebuild the integrity records with:

```bash
python scripts/build_integrity_files.py
```

The PDF figure scripts suppress creation/modification timestamps in PDF metadata so repeated figure generation is byte-reproducible on the same plotting stack.

## License

Code is released under the MIT License. Source datasets accessed from CMIP6/Pangeo
retain their original licenses and terms.




