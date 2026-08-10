# Amazon / CMIP6 reconstruction

The manuscript's climate illustration uses monthly CMIP6 precipitation from the public
Pangeo CMIP6 catalog and does not redistribute raw CMIP6 data.

## Primary reconstruction

Run:

```bash
python scripts/amazon/extract_monthly_amazon_pr_gate_v013.py
```

The script:

- queries the public Pangeo CMIP6 catalog;
- uses the eight model families listed in the manuscript;
- selects `r1i1p1f1` when available and applies the documented grid priority;
- uses historical, SSP1-2.6, SSP2-4.5, and SSP5-8.5 Amon/pr;
- converts precipitation to mm day^-1;
- applies the SREX-AMZ mask;
- reproduces the unweighted grid-cell mean used in the frozen application;
- defines a dry month as regional precipitation < 3.3 mm day^-1;
- writes the monthly reconstruction and deterministic selected-asset table.

Internet access is required because the script reads public CMIP6 Zarr stores.

## Provenance audit

`audit_amazon_cmip6_assets_v016.py` is included to enumerate or compare candidate
CMIP6 assets when provenance auditing is needed. It is not a model-selection step
and is not used to tune the reported scientific result.

Raw CMIP6 files are intentionally excluded from this repository.
