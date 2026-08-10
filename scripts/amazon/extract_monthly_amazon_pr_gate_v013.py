from __future__ import annotations

"""
Monthly Amazon precipitation extractor v0.13 for the manuscript decision gate.

What it does
------------
1. Opens the public Pangeo CMIP6 catalog.
2. Selects one deterministic Amon/pr asset per model and experiment.
3. Opens each Zarr store anonymously, ONE AT A TIME.
4. Computes the AR6-AMZ basin-mean monthly precipitation.
5. Reconstructs the original annual dry-month diagnostic:
       DSL_months = number of months per year with Amazon mean pr < 3.3 mm/day
6. Saves a checkpoint after every model/experiment, so a failure does not erase
   completed work.

The original recovered pipeline used an unweighted grid-cell mean inside the
SREX AMZ mask. This extractor reproduces that definition for comparability.
"""

import argparse
import gc
import sys
import traceback
from pathlib import Path

import gcsfs
import intake
import intake_esm  # registers the esm_datastore plugin in intake
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

EXPERIMENTS = ["historical", "ssp126", "ssp245", "ssp585"]

PREFERRED_MEMBER = "r1i1p1f1"
GRID_PRIORITY = ["gn", "gr", "gr1", "gr2"]

# Definition recovered from the original Amazon diagnostic notebook:
DRY_THRESHOLD_MMDAY = 3.3  # about 100 mm/month

OUT = Path("monthly_gate_output_v013_srex")
PARTS = OUT / "parts"
OUT.mkdir(parents=True, exist_ok=True)
PARTS.mkdir(parents=True, exist_ok=True)


def pr_to_mmday(da: xr.DataArray) -> xr.DataArray:
    """Convert CMIP6 precipitation flux to mm/day."""
    units = str(da.attrs.get("units", "")).lower().replace(" ", "")

    # Standard CMIP6 pr: kg m-2 s-1 == mm s-1.
    if ("kg" in units and "s-1" in units) or "kgm-2s-1" in units:
        out = da * 86400.0
    elif "mm/day" in units or "mmday-1" in units or "mmd-1" in units:
        out = da
    elif units in {"mm", "mmday"}:
        # Amon/pr should normally not arrive this way, but do not silently
        # multiply if the source already declares millimetres.
        out = da
    else:
        print(
            f"[WARN] Unrecognized pr units {da.attrs.get('units')!r}; "
            "assuming kg m-2 s-1 and multiplying by 86400.",
            flush=True,
        )
        out = da * 86400.0

    out = out.copy()
    out.attrs["units"] = "mm day-1"
    return out


def choose_asset_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pick exactly one deterministic Zarr asset for every model/experiment.

    Preference:
      1. r1i1p1f1 when available;
      2. otherwise lexicographically first member;
      3. preferred grid gn, gr, gr1, gr2;
      4. otherwise lexicographically first grid/zstore.
    """
    rows = []

    for model in MODELS:
        for experiment in EXPERIMENTS:
            sub = df[
                (df["source_id"] == model)
                & (df["experiment_id"] == experiment)
            ].copy()

            if sub.empty:
                print(f"[MISSING] {model:16s} {experiment}", flush=True)
                continue

            members = sorted(sub["member_id"].dropna().astype(str).unique())
            member = (
                PREFERRED_MEMBER
                if PREFERRED_MEMBER in members
                else members[0]
            )
            sub = sub[sub["member_id"].astype(str) == member].copy()

            available_grids = sorted(
                sub["grid_label"].dropna().astype(str).unique()
            )
            grid = next(
                (g for g in GRID_PRIORITY if g in available_grids),
                available_grids[0],
            )
            sub = sub[sub["grid_label"].astype(str) == grid].copy()

            # The Pangeo catalog should normally have one Zarr store here.
            # If more than one remains, deterministic lexical order is used.
            sub = sub.sort_values("zstore")
            row = sub.iloc[0].copy()

            rows.append(
                {
                    "source_id": model,
                    "experiment_id": experiment,
                    "member_id": member,
                    "grid_label": grid,
                    "activity_id": row.get("activity_id", ""),
                    "institution_id": row.get("institution_id", ""),
                    "zstore": row["zstore"],
                }
            )

    out = pd.DataFrame(rows)

    expected = len(MODELS) * len(EXPERIMENTS)
    print(
        f"\n[CATALOG] selected {len(out)}/{expected} model-experiment assets.",
        flush=True,
    )
    return out


def amazon_mask(ds: xr.Dataset) -> xr.DataArray:
    """
    Return a boolean mask for the SREX Amazon (AMZ) region.

    IMPORTANT
    ---------
    The earlier prototype referred to an "AR6 AMZ" region, but the IPCC AR6
    reference-region set does not contain an AMZ abbreviation. The standard
    regionmask set that contains AMZ is SREX.

    We therefore use SREX-AMZ here as an explicit reconstruction candidate and
    later compare the resulting annual diagnostics against the archived annual
    Amazon series. That comparison is mandatory before scientific use.
    """
    if "lon" not in ds.coords or "lat" not in ds.coords:
        raise RuntimeError(
            f"Dataset lacks lon/lat coordinates. Coords={list(ds.coords)}"
        )

    srex = regionmask.defined_regions.srex
    abbrevs = list(srex.abbrevs)

    if "AMZ" not in abbrevs:
        raise RuntimeError(
            "SREX AMZ region unexpectedly missing from installed regionmask."
        )

    idx = abbrevs.index("AMZ")
    region_number = int(srex.numbers[idx])

    # Explicit lon/lat arrays are compatible with both older and newer
    # regionmask APIs. Do not use lon_name=/lat_name= keyword arguments.
    try:
        mask = srex.mask(
            ds["lon"],
            ds["lat"],
            wrap_lon=None,
        )
    except TypeError:
        # Compatibility fallback for regionmask versions without wrap_lon.
        mask = srex.mask(
            ds["lon"],
            ds["lat"],
        )

    out = (mask == region_number)

    count = out.sum()
    if hasattr(count.data, "compute"):
        count = count.compute()
    n_cells = int(count.item())
    if n_cells == 0:
        raise RuntimeError(
            "SREX-AMZ mask selected zero grid cells; longitude handling must be inspected."
        )

    return out


def basin_monthly_pr(ds: xr.Dataset) -> xr.DataArray:
    """
    Reproduce the historical pipeline:
      mask to SREX AMZ -> simple grid-cell mean -> mm/day.
    """
    if "pr" not in ds:
        raise RuntimeError(f"'pr' not found. Variables={list(ds.data_vars)}")

    pr = ds["pr"]

    # Some stores may expose a member dimension even after direct asset opening.
    if "member_id" in pr.dims:
        if pr.sizes["member_id"] != 1:
            raise RuntimeError(
                f"Unexpected member_id dimension size={pr.sizes['member_id']}"
            )
        pr = pr.isel(member_id=0, drop=True)

    mask = amazon_mask(ds)
    pr = pr_to_mmday(pr).where(mask)

    spatial_dims = [d for d in pr.dims if d != "time"]
    if len(spatial_dims) < 2:
        raise RuntimeError(
            f"Could not identify two spatial dimensions from {pr.dims}"
        )

    # IMPORTANT: deliberately unweighted to match the recovered original
    # diagnostic notebook.
    return pr.mean(dim=spatial_dims, skipna=True)


def process_one(
    fs: gcsfs.GCSFileSystem,
    row: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    model = str(row["source_id"])
    experiment = str(row["experiment_id"])
    member = str(row["member_id"])
    grid = str(row["grid_label"])
    zstore = str(row["zstore"])

    mapper = fs.get_mapper(zstore)

    # Direct anonymous open avoids the repeated ADC credential lookup that
    # generated the warnings in the earlier run.
    ds = xr.open_zarr(
        mapper,
        consolidated=True,
        decode_times=True,
        chunks={},
    )

    try:
        series = basin_monthly_pr(ds)

        # Trigger only the small basin-mean 1-D series, not the full field.
        series = series.compute()

        years = series["time"].dt.year.values.astype(int)
        months = series["time"].dt.month.values.astype(int)
        values = np.asarray(series.values, dtype=float)

        monthly = pd.DataFrame(
            {
                "cmip6_model": model,
                "experiment": experiment,
                "member_id": member,
                "grid_label": grid,
                "year": years,
                "month": months,
                "amazon_pr_mmday": values,
                "is_dry_month": (values < DRY_THRESHOLD_MMDAY).astype(int),
                "zstore": zstore,
            }
        )

        # Guard against duplicate year-month cells.
        key = ["cmip6_model", "experiment", "year", "month"]
        if monthly.duplicated(key).any():
            dup = monthly.loc[monthly.duplicated(key, keep=False), key]
            raise RuntimeError(
                "Duplicate monthly cells detected:\n"
                + dup.head(20).to_string(index=False)
            )

        annual = (
            monthly.groupby(
                ["cmip6_model", "experiment", "member_id", "grid_label", "year"],
                as_index=False,
            )
            .agg(
                Amazon_pr_mmday=("amazon_pr_mmday", "mean"),
                DSL_months=("is_dry_month", "sum"),
                n_months=("month", "size"),
            )
        )

        incomplete = annual[annual["n_months"] != 12]
        if not incomplete.empty:
            print(
                f"[WARN] {model} {experiment}: "
                f"{len(incomplete)} incomplete calendar years. "
                "They are retained and flagged by n_months.",
                flush=True,
            )

        return monthly, annual

    finally:
        ds.close()
        del ds
        gc.collect()



def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract monthly CMIP6 Amazon precipitation for the decision gate."
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Process only the first selected model/experiment asset.",
    )
    parser.add_argument(
        "--clean-smoke",
        action="store_true",
        help="Delete only smoke-test checkpoint files before running.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print("[1/5] Opening Pangeo CMIP6 catalog...", flush=True)
    cat = intake.open_esm_datastore(CATALOG)

    print("[2/5] Searching monthly precipitation assets...", flush=True)
    sub = cat.search(
        activity_id=["CMIP", "ScenarioMIP"],
        source_id=MODELS,
        experiment_id=EXPERIMENTS,
        variable_id="pr",
        table_id="Amon",
    )
    catalog_df = sub.df.copy()

    if catalog_df.empty:
        raise RuntimeError("Catalog search returned zero assets.")

    selected = choose_asset_rows(catalog_df)
    selected.to_csv(OUT / "selected_assets.csv", index=False)

    missing_pairs = {
        (m, e) for m in MODELS for e in EXPERIMENTS
    } - set(zip(selected["source_id"], selected["experiment_id"]))

    if args.smoke:
        selected = selected.head(1).copy()
        print(
            "\n[SMOKE] Processing only the first asset. "
            "If it succeeds, rerun without --smoke.",
            flush=True,
        )
    if missing_pairs:
        pd.DataFrame(
            sorted(missing_pairs),
            columns=["cmip6_model", "experiment"],
        ).to_csv(OUT / "missing_assets.csv", index=False)
        print(
            f"[WARN] {len(missing_pairs)} requested model-experiment pairs "
            "are absent from the catalog selection.",
            flush=True,
        )

    print(
        "[3/5] Creating anonymous Google Cloud filesystem "
        "(no credentials required)...",
        flush=True,
    )
    fs = gcsfs.GCSFileSystem(token="anon")

    completed = []
    failed = []

    print(
        "[MASK] Reconstruction candidate: regionmask SREX-AMZ. "
        "Annual values MUST be compared with the archived Amazon master before use.",
        flush=True,
    )
    print("[4/5] Processing one asset at a time...", flush=True)

    for i, row in selected.reset_index(drop=True).iterrows():
        model = str(row["source_id"])
        experiment = str(row["experiment_id"])
        stem = f"{model}_{experiment}"

        monthly_file = PARTS / f"{stem}_monthly.csv"
        annual_file = PARTS / f"{stem}_annual.csv"

        print(
            f"\n[{i+1:02d}/{len(selected):02d}] "
            f"{model:16s} {experiment:10s} "
            f"{row['member_id']} {row['grid_label']}",
            flush=True,
        )

        # Resume safely after interruption.
        if monthly_file.exists() and annual_file.exists():
            print("    checkpoint exists -> SKIP", flush=True)
            completed.append(
                {
                    "cmip6_model": model,
                    "experiment": experiment,
                    "status": "checkpoint",
                }
            )
            continue

        try:
            monthly, annual = process_one(fs, row)

            monthly.to_csv(monthly_file, index=False)
            annual.to_csv(annual_file, index=False)

            print(
                f"    OK: {len(monthly)} months, "
                f"{annual['year'].min()}-{annual['year'].max()}, "
                f"mean DSL={annual['DSL_months'].mean():.2f}",
                flush=True,
            )
            completed.append(
                {
                    "cmip6_model": model,
                    "experiment": experiment,
                    "status": "ok",
                    "n_months": len(monthly),
                    "year_min": int(annual["year"].min()),
                    "year_max": int(annual["year"].max()),
                }
            )

        except Exception as exc:
            print(f"    FAILED: {type(exc).__name__}: {exc}", flush=True)
            traceback.print_exc()
            failed.append(
                {
                    "cmip6_model": model,
                    "experiment": experiment,
                    "member_id": row["member_id"],
                    "grid_label": row["grid_label"],
                    "zstore": row["zstore"],
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )

    pd.DataFrame(completed).to_csv(OUT / "completed.csv", index=False)
    pd.DataFrame(failed).to_csv(OUT / "failed.csv", index=False)

    print("\n[5/5] Combining checkpoints...", flush=True)

    monthly_parts = [
        pd.read_csv(p) for p in sorted(PARTS.glob("*_monthly.csv"))
    ]
    annual_parts = [
        pd.read_csv(p) for p in sorted(PARTS.glob("*_annual.csv"))
    ]

    if not monthly_parts or not annual_parts:
        if args.smoke:
            print(
                "[SMOKE FAILED] No checkpoint was produced. "
                "Inspect the single traceback above; no combine step was attempted.",
                flush=True,
            )
            return
        raise RuntimeError("No successful monthly/annual checkpoints were produced.")

    monthly_all = pd.concat(monthly_parts, ignore_index=True)
    annual_all = pd.concat(annual_parts, ignore_index=True)

    monthly_all = monthly_all.sort_values(
        ["cmip6_model", "experiment", "year", "month"]
    ).reset_index(drop=True)
    annual_all = annual_all.sort_values(
        ["cmip6_model", "experiment", "year"]
    ).reset_index(drop=True)

    monthly_all.to_csv(
        OUT / "amazon_monthly_pr_matched.csv",
        index=False,
    )
    annual_all.to_csv(
        OUT / "amazon_annual_reconstructed_from_monthly.csv",
        index=False,
    )

    # Parquet is optional.
    try:
        monthly_all.to_parquet(
            OUT / "amazon_monthly_pr_matched.parquet",
            index=False,
        )
        annual_all.to_parquet(
            OUT / "amazon_annual_reconstructed_from_monthly.parquet",
            index=False,
        )
    except Exception as exc:
        print(f"[INFO] Parquet not written: {exc}", flush=True)

    coverage = (
        annual_all.groupby(["cmip6_model", "experiment"], as_index=False)
        .agg(
            year_min=("year", "min"),
            year_max=("year", "max"),
            n_years=("year", "size"),
            incomplete_years=("n_months", lambda s: int((s != 12).sum())),
            mean_dsl=("DSL_months", "mean"),
            mean_pr=("Amazon_pr_mmday", "mean"),
        )
    )
    coverage.to_csv(OUT / "coverage_summary.csv", index=False)

    print("\nDONE.", flush=True)
    print(f"Output directory: {OUT.resolve()}", flush=True)
    print(
        f"Successful model-experiment series: "
        f"{coverage[['cmip6_model','experiment']].drop_duplicates().shape[0]}",
        flush=True,
    )
    print(f"Failed assets: {len(failed)}", flush=True)
    print("\nCoverage:", flush=True)
    print(coverage.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
