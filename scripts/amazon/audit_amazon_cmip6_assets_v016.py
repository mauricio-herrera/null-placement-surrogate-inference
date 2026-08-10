from __future__ import annotations

"""
Amazon CMIP6 asset provenance audit v0.16.

Purpose
-------
After ruling out:
  * SREX AMZ off-by-one region-number bug;
  * Giorgi AMZ;
  * simple weighting/bounding-box alternatives,

test whether the archived Amazon diagnostics were generated from a different
CMIP6 Amon/pr asset (grid/version) than the one selected in the v0.13 rebuild.

This is a provenance audit, not a model-selection exercise.

Modes
-----
1) INVENTORY ONLY (recommended first; no field computation):
   python audit_amazon_cmip6_assets_v016.py \
       --legacy ./cmip6_annual_amazon_master.csv \
       --inventory-only

2) SMOKE all available ACCESS-CM2 assets:
   python audit_amazon_cmip6_assets_v016.py \
       --legacy ./cmip6_annual_amazon_master.csv \
       --model ACCESS-CM2

3) ALL 8 models:
   python audit_amazon_cmip6_assets_v016.py \
       --legacy ./cmip6_annual_amazon_master.csv

Outputs
-------
amazon_cmip6_asset_audit_v016/
  asset_inventory.csv
  asset_inventory_summary.csv
  candidate_asset_metrics.csv
  candidate_asset_by_year.csv
  best_asset_by_model.csv
  failed.csv
"""

import argparse
import gc
import re
import traceback
from pathlib import Path

import gcsfs
import intake
import intake_esm
import numpy as np
import pandas as pd
import regionmask
import xarray as xr


CATALOG = "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"

MODELS = [
    "ACCESS-CM2",
    "CESM2-WACCM",
    "CanESM5",
    "INM-CM4-8",
    "INM-CM5-0",
    "MIROC6",
    "MPI-ESM1-2-HR",
    "MPI-ESM1-2-LR",
]

PREFERRED_MEMBER = "r1i1p1f1"
DRY_THRESHOLD_MMDAY = 3.3


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--legacy", required=True)
    p.add_argument("--out", default="amazon_cmip6_asset_audit_v016")
    p.add_argument(
        "--inventory-only",
        action="store_true",
        help="Only enumerate candidate Pangeo CMIP6 assets; do not open Zarr fields.",
    )
    p.add_argument(
        "--model",
        choices=MODELS,
        default=None,
        help="Restrict field evaluation to one model.",
    )
    return p.parse_args()


def detect_col(df, names):
    for name in names:
        if name in df.columns:
            return name
    raise KeyError(f"None of {names} found. Columns={list(df.columns)}")


def load_legacy(path):
    df = pd.read_csv(path)
    mc = detect_col(df, ["cmip6_model", "model", "source_id"])
    ec = detect_col(df, ["scenario", "experiment", "experiment_id"])
    yc = detect_col(df, ["year"])
    pc = detect_col(df, ["Amazon_pr_mmday", "amazon_pr_mmday"])
    dc = detect_col(df, ["DSL_months", "dsl_months"])

    x = df[[mc, ec, yc, pc, dc]].copy()
    x.columns = [
        "cmip6_model", "experiment", "year", "legacy_pr", "legacy_dsl"
    ]
    x["cmip6_model"] = x["cmip6_model"].astype(str)
    x["experiment"] = x["experiment"].astype(str)
    x["year"] = x["year"].astype(int)

    x = x[
        x["cmip6_model"].isin(MODELS)
        & (x["experiment"] == "historical")
    ].copy()

    # Collapse repeated historical branches.
    x = (
        x.groupby(["cmip6_model", "year"], as_index=False)
        .agg(
            legacy_pr=("legacy_pr", "mean"),
            legacy_dsl=("legacy_dsl", "mean"),
        )
    )
    return x


def extract_version(zstore: str) -> str:
    # Usually .../gn/v20191108/
    m = re.search(r"/(v\d{8,})/?$", str(zstore))
    return m.group(1) if m else ""


def inventory_assets(cat):
    q = cat.search(
        activity_id="CMIP",
        source_id=MODELS,
        experiment_id="historical",
        variable_id="pr",
        table_id="Amon",
        member_id=PREFERRED_MEMBER,
    )
    df = q.df.copy()

    keep = [
        c for c in [
            "activity_id",
            "institution_id",
            "source_id",
            "experiment_id",
            "member_id",
            "table_id",
            "variable_id",
            "grid_label",
            "zstore",
            "dcpp_init_year",
        ] if c in df.columns
    ]
    df = df[keep].copy()

    # zstore uniquely identifies the physical candidate used here.
    df = df.drop_duplicates(subset=["source_id", "zstore"]).copy()
    df["version"] = df["zstore"].map(extract_version)

    return df.sort_values(
        ["source_id", "grid_label", "version", "zstore"]
    ).reset_index(drop=True)


def pr_to_mmday(da):
    units = str(da.attrs.get("units", "")).lower().replace(" ", "")
    if ("kg" in units and "s-1" in units) or "kgm-2s-1" in units:
        return da * 86400.0
    if "mm" in units and ("day" in units or "d-1" in units):
        return da
    print(
        f"[WARN] pr units={da.attrs.get('units')!r}; "
        "assuming kg m-2 s-1.",
        flush=True,
    )
    return da * 86400.0


def srex_amz_mask(ds):
    regions = regionmask.defined_regions.srex
    pos = list(regions.abbrevs).index("AMZ")
    number = int(regions.numbers[pos])
    try:
        mask = regions.mask(ds["lon"], ds["lat"], wrap_lon=None)
    except TypeError:
        mask = regions.mask(ds["lon"], ds["lat"])
    return mask == number


def annual_series(ds):
    pr = pr_to_mmday(ds["pr"])
    mask = srex_amz_mask(ds)

    spatial_dims = [d for d in pr.dims if d != "time"]
    if len(spatial_dims) < 2:
        raise RuntimeError(f"Unexpected pr dimensions: {pr.dims}")

    # Match the recovered intended operation:
    # spatial mean -> monthly regional series -> dry threshold.
    s = pr.where(mask).mean(spatial_dims, skipna=True).compute()

    years = s["time"].dt.year.values.astype(int)
    vals = np.asarray(s.values, dtype=float)

    monthly = pd.DataFrame({
        "year": years,
        "pr": vals,
    })
    monthly["dry"] = (monthly["pr"] < DRY_THRESHOLD_MMDAY).astype(int)

    return (
        monthly.groupby("year", as_index=False)
        .agg(
            candidate_pr=("pr", "mean"),
            candidate_dsl=("dry", "sum"),
            n_months=("dry", "size"),
        )
    )


def metrics(g):
    g = g.dropna().copy()

    return {
        "n": int(len(g)),
        "pr_corr": float(g["candidate_pr"].corr(g["legacy_pr"])),
        "pr_rmse": float(
            np.sqrt(np.mean((g["candidate_pr"] - g["legacy_pr"]) ** 2))
        ),
        "pr_bias": float(
            np.mean(g["candidate_pr"] - g["legacy_pr"])
        ),
        "dsl_corr": float(
            g["candidate_dsl"].corr(g["legacy_dsl"])
        ) if g["candidate_dsl"].std() > 0 and g["legacy_dsl"].std() > 0
        else np.nan,
        "dsl_mae": float(
            np.mean(np.abs(g["candidate_dsl"] - g["legacy_dsl"]))
        ),
        "dsl_exact_fraction": float(
            np.mean(g["candidate_dsl"] == g["legacy_dsl"])
        ),
        "dsl_total_candidate": float(g["candidate_dsl"].sum()),
        "dsl_total_legacy": float(g["legacy_dsl"].sum()),
    }


def main():
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("[1/4] Opening Pangeo catalog...", flush=True)
    cat = intake.open_esm_datastore(CATALOG)

    print("[2/4] Building exact historical Amon/pr asset inventory...", flush=True)
    inv = inventory_assets(cat)
    inv.to_csv(out / "asset_inventory.csv", index=False)

    summary = (
        inv.groupby("source_id", as_index=False)
        .agg(
            n_assets=("zstore", "nunique"),
            n_grids=("grid_label", "nunique"),
            grids=("grid_label", lambda x: ";".join(sorted(set(map(str, x))))),
            versions=("version", lambda x: ";".join(sorted(
                set(v for v in map(str, x) if v)
            ))),
        )
    )
    summary.to_csv(out / "asset_inventory_summary.csv", index=False)

    print("\nAsset inventory:")
    print(summary.to_string(index=False), flush=True)

    if args.inventory_only:
        print("\nINVENTORY ONLY complete.", flush=True)
        print(f"Outputs: {out.resolve()}", flush=True)
        return

    legacy = load_legacy(args.legacy)
    fs = gcsfs.GCSFileSystem(token="anon")

    work = inv.copy()
    if args.model:
        work = work[work["source_id"] == args.model].copy()

    print(
        f"\n[3/4] Evaluating {len(work)} unique Zarr assets "
        "with fixed correct SREX-AMZ...",
        flush=True,
    )

    all_years = []
    metric_rows = []
    failures = []

    for i, row in work.reset_index(drop=True).iterrows():
        model = str(row["source_id"])
        zstore = str(row["zstore"])
        grid = str(row.get("grid_label", ""))
        version = str(row.get("version", ""))
        asset_id = f"{model}|{grid}|{version}|{i}"

        print(
            f"\n[{i+1}/{len(work)}] {model} grid={grid} version={version}",
            flush=True,
        )
        print(f"    {zstore}", flush=True)

        try:
            ds = xr.open_zarr(
                fs.get_mapper(zstore),
                consolidated=True,
                decode_times=True,
                chunks={},
            )
            try:
                annual = annual_series(ds)
            finally:
                ds.close()
                del ds
                gc.collect()

            ref = legacy[legacy["cmip6_model"] == model].copy()
            joined = annual.merge(ref, on="year", how="inner")

            joined["cmip6_model"] = model
            joined["grid_label"] = grid
            joined["version"] = version
            joined["zstore"] = zstore
            joined["asset_id"] = asset_id
            all_years.append(joined)

            mm = metrics(joined)
            metric_rows.append({
                "cmip6_model": model,
                "grid_label": grid,
                "version": version,
                "zstore": zstore,
                "asset_id": asset_id,
                **mm,
            })

            print(
                "    "
                f"rP={mm['pr_corr']:.6f}, "
                f"RMSE={mm['pr_rmse']:.6f}, "
                f"rDSL={mm['dsl_corr'] if np.isfinite(mm['dsl_corr']) else np.nan:.6f}, "
                f"MAE_DSL={mm['dsl_mae']:.6f}, "
                f"exact={mm['dsl_exact_fraction']:.3f}",
                flush=True,
            )

        except Exception as exc:
            print(f"    FAILED: {type(exc).__name__}: {exc}", flush=True)
            traceback.print_exc()
            failures.append({
                "cmip6_model": model,
                "grid_label": grid,
                "version": version,
                "zstore": zstore,
                "error_type": type(exc).__name__,
                "error": str(exc),
            })

    pd.DataFrame(failures).to_csv(out / "failed.csv", index=False)

    if not metric_rows:
        raise RuntimeError("No assets were successfully evaluated.")

    metrics_df = pd.DataFrame(metric_rows)

    # Rank independently within each model. Precedence: annual precipitation
    # checksum first, then dry-month checksum. This is provenance ranking,
    # not a fitted scientific model.
    metrics_df["rank_score"] = (
        metrics_df["pr_rmse"]
        + metrics_df["dsl_mae"]
    )
    metrics_df = metrics_df.sort_values(
        ["cmip6_model", "rank_score", "pr_rmse"]
    ).reset_index(drop=True)
    metrics_df.to_csv(out / "candidate_asset_metrics.csv", index=False)

    if all_years:
        pd.concat(all_years, ignore_index=True).to_csv(
            out / "candidate_asset_by_year.csv",
            index=False,
        )

    best = (
        metrics_df.groupby("cmip6_model", as_index=False)
        .first()
    )
    best.to_csv(out / "best_asset_by_model.csv", index=False)

    print("\n[4/4] Best candidate by model:", flush=True)
    print(
        best[[
            "cmip6_model", "grid_label", "version",
            "pr_corr", "pr_rmse", "dsl_corr",
            "dsl_mae", "dsl_exact_fraction",
        ]].to_string(index=False),
        flush=True,
    )

    print("\nDecision rule:", flush=True)
    print(
        "  If an alternative grid/version approaches exact annual reproduction, "
        "use its zstore provenance for the monthly rebuild.",
        flush=True,
    )
    print(
        "  If every catalog asset remains far from the archived diagnostics, "
        "the February pipeline used a different local/ESGF asset or a different "
        "spatial implementation; do not run the final monthly null yet.",
        flush=True,
    )
    print(f"\nOutputs: {out.resolve()}", flush=True)


if __name__ == "__main__":
    main()
