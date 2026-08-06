"""
Compute dye statistics (mean and standard deviation)
for years 450-500 of all experiments.

Author: Laura Endres
"""

from pathlib import Path
import xarray as xr
from tqdm import tqdm


# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path("/nfs/annie/earpal/database/experiments")

FIRST_YEAR = 450
LAST_YEAR = 500

TIME_CHUNK = 50


# =============================================================================
# Helper functions
# =============================================================================

def open_experiment(exp):
    """
    Open one experiment and rename dye variables to dye00...dye08.
    """

    files = str(BASE_DIR / exp / "time_series" / f"{exp}.dye0?.annual.nc")

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        parallel=True,
        chunks={"t": TIME_CHUNK},
    )

    rename_dict = {
        old: f"dye{i:02d}"
        for i, old in enumerate(ds.data_vars)
    }

    return ds.rename_vars(rename_dict)


def get_time_dim(ds):
    """
    Return the name of the time dimension.
    """

    if "t" in ds.dims:
        return "t"

    if "time" in ds.dims:
        return "time"

    raise ValueError("Could not identify the time dimension.")


def select_analysis_period(ds):
    """
    Select years 450-499 (50 years).

    NOTE:
    Adjust if your time coordinate is 1...500.
    """

    time_dim = get_time_dim(ds)

    if time_dim in ds.coords:
        try:
            return ds.sel({time_dim: slice(FIRST_YEAR, LAST_YEAR)})
        except Exception:
            pass

    return ds.isel({time_dim: slice(FIRST_YEAR, LAST_YEAR)})


def compute_statistics(ds):
    """
    Compute mean and standard deviation.
    """

    ds = select_analysis_period(ds)

    time_dim = get_time_dim(ds)

    # IMPORTANT: skip NaNs safely
    mean = ds.mean(time_dim, skipna=True)

    std = ds.std(time_dim, skipna=True)

    out = xr.Dataset()

    for dye in ds.data_vars:

        out[f"{dye}_mean"] = mean[dye]
        out[f"{dye}_std"] = std[dye]

    return out


def save_statistics(ds_stats, outfile):
    """
    Save compressed NetCDF.
    """

    outfile.parent.mkdir(parents=True, exist_ok=True)

    encoding = {
        var: {
            "zlib": True,
            "complevel": 4,
        }
        for var in ds_stats.data_vars
    }

    ds_stats.to_netcdf(outfile, encoding=encoding)




# =============================================================================
# Main pipeline
# =============================================================================

def run_pipeline(
    experiments,
    outdir,
    keep_in_memory=False,
    save=True,
    overwrite=False,
    exclude=None,
    preflight=True,
):
    """
    Main dye/statistics pipeline with optional preflight check.
    """

    
    
    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)

    outdir = Path(outdir)
    results = {}

    # =========================================================
    # 1. PREFLIGHT CHECK (VERY LIGHTWEIGHT)
    # =========================================================
    if preflight:
        print("\n================ PREFLIGHT CHECK ================\n")

        for state, scenarios in experiments.items():

            print(f"\n--- {state.upper()} ---")

            for sce, info in scenarios.items():

                exp = info["exp"]

                if sce in exclude or exp in exclude:
                    print(f"{sce:15s} -> SKIPPED (excluded)")
                    continue

                try:
                    ds = open_experiment(exp)  # lightweight open
                    tlen = ds.sizes.get("t", None)
                    ndyes = len([v for v in ds.data_vars if v.startswith("dye")])

                    print(f"{sce:15s} | exp={exp:8s} | time={tlen} | dyes={ndyes}")

                    ds.close()

                except Exception as e:
                    print(f"{sce:15s} | exp={exp:8s} | ERROR: {e}")

        print("\n================================================\n")

    # =========================================================
    # 2. MAIN COMPUTATION
    # =========================================================

    for state, scenarios in experiments.items():

        print(f"\n===== {state.upper()} =====")

        results[state] = {}

        for sce, info in tqdm(scenarios.items(), desc=state, unit="experiment"):

            exp = info["exp"]

            # -----------------------------------------------------
            # EXCLUSION
            # -----------------------------------------------------
            if sce in exclude or exp in exclude:
                continue

            outfile = outdir / state / f"{sce}_mean_std.nc"
            outfile.parent.mkdir(parents=True, exist_ok=True)

            # -----------------------------------------------------
            # SKIP IF EXISTS
            # -----------------------------------------------------
            if outfile.exists() and not overwrite:

                if keep_in_memory:
                    results[state][exp] = xr.open_dataset(outfile, chunks="auto")
                else:
                    results[state][exp] = outfile

                continue

            # -----------------------------------------------------
            # COMPUTE
            # -----------------------------------------------------
            ds = open_experiment(exp)

            stats = compute_statistics(ds)

            # -----------------------------------------------------
            # SAVE
            # -----------------------------------------------------
            if save:
                stats.to_netcdf(outfile)

            # -----------------------------------------------------
            # STORE
            # -----------------------------------------------------
            if keep_in_memory:
                results[state][exp] = stats
            else:
                results[state][exp] = outfile

            ds.close()

    print("\n✓ Pipeline finished.")

    return results


### Sea ice lookup
from pathlib import Path
from tqdm import tqdm
import xarray as xr
import pickle


def build_seaice_lookup(
    experiments,
    outdir="myintermediates",
    overwrite=False,
    exclude=None,
):
    """
    Build or load a lookup table containing climatological
    March and September sea-ice concentration fields.

    Returns
    -------
    seaice_lookup : dict

    seaice_lookup[state][exp] = {
        "march": DataArray,
        "sept": DataArray,
    }
    """

    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / "seaice_lookup.pkl"

    # -------------------------------------------------
    # Load existing lookup
    # -------------------------------------------------
    if outfile.exists() and not overwrite:

        print(f"Loading existing lookup: {outfile}")

        with open(outfile, "rb") as f:
            return pickle.load(f)

    # -------------------------------------------------
    # Build lookup
    # -------------------------------------------------
    print("Building sea-ice lookup...")

    seaice_lookup = {}

    for state, scen_dict in tqdm(experiments.items(), desc="States"):

        seaice_lookup[state] = {}

        for scenario, info in tqdm(
            scen_dict.items(),
            leave=False,
            desc=state
        ):

            exp = info["exp"]

            if scenario in exclude or exp in exclude:
                continue

            fn = (
                f"/nfs/annie/earpal/database/experiments/"
                f"{exp}/time_series/{exp}.iceconc.monthly.nc"
            )

            ds = xr.open_dataset(fn)

            clim = ds.groupby("t.month").mean("t")

            seaice_lookup[state][exp] = {
                "march": clim["iceconc_mm_srf"]
                    .sel(month=3)
                    .isel(surface=0)
                    .load(),

                "sept": clim["iceconc_mm_srf"]
                    .sel(month=9)
                    .isel(surface=0)
                    .load(),
            }

            ds.close()

    # -------------------------------------------------
    # Save lookup
    # -------------------------------------------------
    with open(outfile, "wb") as f:
        pickle.dump(seaice_lookup, f)

    print(f"✓ Saved lookup to {outfile}")

    return seaice_lookup



### mld lookup
from pathlib import Path
from tqdm import tqdm
import pickle
import xarray as xr


def build_mld_lookup(
    experiments,
    outdir="myintermediates",
    overwrite=False,
    exclude=None,
):
    """
    Build or load lookup table of winter (DJFM) mean mixed-layer depth.

    Returns
    -------
    mld_lookup[state][exp] = DataArray
    """

    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / "mld_lookup.pkl"

    # -------------------------------------------------
    # Load existing lookup
    # -------------------------------------------------
    if outfile.exists() and not overwrite:

        print(f"Loading {outfile}")

        with open(outfile, "rb") as f:
            mld_lookup = pickle.load(f)

        if exclude:
            for state in list(mld_lookup.keys()):
                for exp in list(mld_lookup[state].keys()):
                    if exp in exclude:
                        del mld_lookup[state][exp]

        return mld_lookup

    # -------------------------------------------------
    # Build lookup
    # -------------------------------------------------

    database = "/nfs/see-fs-01_users/eelse/database"

    mld_lookup = {}

    for state, scen_dict in tqdm(experiments.items(), desc="States"):

        mld_lookup[state] = {}

        for scenario, info in tqdm(
            scen_dict.items(),
            desc=state,
            leave=False,
        ):

            exp = info["exp"]

            if scenario in exclude or exp in exclude:
                continue

            fn = (
                f"{database}/{exp}/time_series/"
                f"{exp}.oceanmixedpf.monthly.nc"
            )

            ds = xr.open_dataset(fn)

            mld = (
                ds.mixLyrDpth_mm_uo
                .isel(unspecified=0, drop=True)[1:]
            )

            # winter climatology (DJFM)
            winter = (
                mld.where(
                    mld.t.dt.month.isin([12, 1, 2, 3]),
                    drop=True
                )
                .mean("t")
                .load()
            )

            mld_lookup[state][exp] = winter

            ds.close()

    with open(outfile, "wb") as f:
        pickle.dump(mld_lookup, f)

    print(f"✓ Saved {outfile}")

    return mld_lookup


# AMOC time series

from pathlib import Path
from tqdm import tqdm
import pickle
import xarray as xr
# pylaeoclim_leeds is only needed to *rebuild* the AMOC lookup from raw HadCM3
# output (util.ButterLowPass, below). The figure scripts ship the precomputed
# amoc_lookup.pkl, so this import is not required for plotting. Import it lazily
# so the module loads without pylaeoclim_leeds installed; the error is only
# raised if a rebuild is actually triggered (cache missing / overwrite=True).
try:
    import pylaeoclim_leeds.util_hadcm3 as util
except ImportError:
    util = None



def build_amoc_lookup(
    experiments,
    outdir="myintermediates",
    overwrite=False,
    exclude=None,
):
    """
    Build or load AMOC time series lookup.

    Returns
    -------
    amoc_lookup[state][exp] = {
        "raw": DataArray,
        "filtered": DataArray
    }

    AMOC is calculated as maximum overturning at 26.5N.
    """

    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / "amoc_lookup.pkl"


    # -------------------------------------------------
    # Load existing lookup
    # -------------------------------------------------
    if outfile.exists() and not overwrite:

        print(f"Loading {outfile}")

        with open(outfile, "rb") as f:
            amoc_lookup = pickle.load(f)

        if exclude:
            for state in list(amoc_lookup.keys()):
                for exp in list(amoc_lookup[state].keys()):
                    if exp in exclude:
                        del amoc_lookup[state][exp]

        return amoc_lookup


    # -------------------------------------------------
    # Build lookup
    # -------------------------------------------------

    database = "/nfs/see-fs-01_users/eelse/database"

    if util is None:
        raise ImportError(
            "Rebuilding the AMOC lookup requires the external 'pylaeoclim_leeds' "
            "package (util.ButterLowPass), which is not installed. The figure "
            "scripts ship the precomputed amoc_lookup.pkl in data/intermediates/ "
            "and do not need a rebuild; ensure that file is present and call with "
            "overwrite=False."
        )

    filt = util.ButterLowPass(
        order=1,
        fc=2*10**-3,
        fs=1,
        mult=2,
    )

    amoc_lookup = {}


    for state, scen_dict in tqdm(experiments.items(), desc="States"):

        amoc_lookup[state] = {}

        for scenario, info in tqdm(
            scen_dict.items(),
            desc=state,
            leave=False,
        ):

            exp = info["exp"]

            if scenario in exclude or exp in exclude:
                continue


            fn = (
                f"{database}/{exp}/time_series/"
                f"{exp}.merid.annual.nc"
            )

            ds = xr.open_dataset(fn)

            ds = ds.drop_duplicates(dim="t")


            # AMOC at 26.5N
            amoc = (
                ds.Merid_Atlantic
                .sel(latitude=26.5, method="nearest")
                .max("depth")
                .sortby("t")
                .load()
            )


            # filtered time series
            amoc_filt = xr.DataArray(
                filt.process(amoc.values),
                coords=amoc.coords,
                dims=amoc.dims,
                name="AMOC_filtered",
            )


            amoc_lookup[state][exp] = {
                "raw": amoc,
                "filtered": amoc_filt,
            }


            ds.close()


    # -------------------------------------------------
    # Save
    # -------------------------------------------------

    with open(outfile, "wb") as f:
        pickle.dump(amoc_lookup, f)


    print(f"✓ Saved {outfile}")

    return amoc_lookup


import xarray as xr

def build_dye_field_lookup(dye_data, depth_index=0):
    """
    Build dye lookup:
    
    dye_field_lookup[mode]["mean"][dye]
    dye_field_lookup[mode]["std"][dye]

    Each dye is kept independent.
    Averaging is only across experiments within a circulation mode.
    """

    dye_field_lookup = {}

    for mode, experiments in dye_data.items():

        dye_field_lookup[mode] = {
            "mean": {},
            "std": {}
        }

        # identify all dyes available across experiments
        dyes = sorted({
            v.split("_")[0]
            for ds in experiments.values()
            for v in ds.data_vars
            if v.endswith("_mean")
        })

        for dye in dyes:

            dye_experiment_fields = []

            for exp, ds in experiments.items():

                varname = f"{dye}_mean"

                # skip if dye not available in this experiment
                if varname not in ds:
                    continue

                # surface field
                field = ds[varname].isel(depth_1=depth_index)

                dye_experiment_fields.append(field)

            # stack experiments for this dye
            stack = xr.concat(
                dye_experiment_fields,
                dim="experiment"
            )

            # average and spread across experiments
            dye_field_lookup[mode]["mean"][dye] = (
                stack.mean(dim="experiment")
            )

            dye_field_lookup[mode]["std"][dye] = (
                stack.std(dim="experiment")
            )

    return dye_field_lookup



