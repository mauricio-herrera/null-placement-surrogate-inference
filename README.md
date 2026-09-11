# Null transport under temporal coarse-graining
## Reproducibility package — v2.0.0

This package contains the submission-oriented manuscript reconstruction after replacing the former CMIP6 Amazon illustration with a prospectively frozen SINCA PM2.5 application.

## Archived releases

- **v2.0.0** — null transport, lumpability, mechanism diagnostics, and prospectively frozen SINCA PM2.5 application; DOI `10.5281/zenodo.22713989`.
- **v1.1.0** — SERRA metadata release; DOI `10.5281/zenodo.21910629`.
- **v1.0.0** — original reproducibility release; DOI `10.5281/zenodo.21873328`.

## Core deliverables

- `manuscript_final.pdf` / `manuscript_final.tex`
- `body_final.tex`
- `Supplementary_Material_Final.pdf` / `.tex`
- `references_final.bib` and compiled `references_final.bbl`
- publication figures in `figures/`

## Main scientific changes relative to the previous reconstruction

1. Added the 31-station SINCA PM2.5 application, with cohort and inferential protocol frozen before application-level surrogate outputs were inspected.
2. Under the fixed 50 µg m^-3 event threshold, the primary Ferro-Segers analysis gives 29 index-only, 0 native-only, and 2 joint rejections; the region sign-flip p-value is 0.003906.
3. Probability equalization reduces index-only decisions from 29 to 5 and reduces the index-minus-native null-center gap at all 31 stations.
4. The prespecified SINCA null-width-contraction prediction is reported as not supported (rho=-0.379; region-bootstrap 95% CI [-0.624, 0.136]).
5. The mechanism is therefore framed as full projected-null distributional distortion, not universal variance contraction.
6. The exposure-standardized missingness sensitivity was completed: among 22 technically defined stations, 22 index-level vs 1 native rejection (21 index-only, 1 both).
7. The non-discriminating Amazon-CMIP6 illustration has been removed from the final manuscript and supplement.

## Reproducibility

`reproducibility/SINCA/` contains the frozen cohort, protocol, machine-readable manifest, analysis-ready input, station-level results, probability-equalization results, exposure sensitivity, and executable scripts. The synthetic/within-fiber mechanism diagnostics are under `reproducibility/diagnostics/`.

The archived v2.0.0 reproducibility package is available at https://doi.org/10.5281/zenodo.22713989.

The PDFs compile with the included `build.sh` without LaTeX warnings or overfull/underfull box warnings in the current environment.
