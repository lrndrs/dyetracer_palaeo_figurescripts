from pathlib import Path
import xarray as xr
from tqdm import tqdm

def build_d18O_lookup(df, column, dye_index=None):
    """
    Build dye → d18O lookup dictionary from a pandas DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Your gdg dataframe

    column : str
        Column containing d18O values (or region/scenario column)

    dye_index : list or None
        Optional list of dye names (dye00...dye08).
        If None, assumes 9 dyes in order.
    """

    import numpy as np
    import xarray as xr

    if dye_index is None:
        dye_index = [f"dye{i:02d}" for i in range(9)]

    values = df[column].values

    if len(values) != len(dye_index):
        raise ValueError(
            f"Mismatch: {len(values)} values vs {len(dye_index)} dyes"
        )

    return dict(zip(dye_index, values))



def build_d18O_mean_lookup(
    experiments,
    d18O_results,
    lon_min,lon_max,
    lat_min,lat_max,
    lsm=None,
    scenarios=("17.8k", "18.2k", "19.4k", "20.7k"),
):
    """
    Precompute area-weighted mean d18O per mode & scenario.
    Skips missing experiments safely.
    """

    mean_lookup = {}
    mean_lookup_old = {}

    # cosine latitude weights (precompute once)
    #weights = np.cos(np.deg2rad(lat))
    #weights = xr.DataArray(weights, coords={lat.name: lat}, dims=[lat.name])

    for mode, scen_dict in tqdm(experiments.items(), desc="modes"):

        mean_lookup[mode] = {}
        mean_lookup_old[mode] = {}

        for scenario in scenarios:

            # -----------------------------------------
            # skip missing scenario (IMPORTANT)
            # -----------------------------------------
            if scenario not in scen_dict:
                continue

            exp = scen_dict[scenario]["exp"]

            # -----------------------------------------
            # skip missing computed results
            # -----------------------------------------
            #if exp not in d18O_results[mode]:
            #    continue


            # -----------------------------------------
            # handle both d18O simulated and predicted fields
            # -----------------------------------------

            field = None

            if exp in d18O_results[mode]:
                print("Using experiment format for", exp)
                field = d18O_results[mode][exp]["d18O"].isel(depth_1=0)
            
            else:
                pred_key = f"{scenario}_d18O_mean"
                print("Looking for", pred_key)
            
                if pred_key in d18O_results[mode]:
                    print("Found prediction")
                    field = d18O_results[mode][pred_key]
                else:
                    print("Not found!")
                    continue

            


            #field = d18O_results[mode][exp]["d18O"].isel(depth_1=0)

            # ---
            # Apply land-sea mask
            # ---

            
            if lsm is not None:
                field = xr.where(lsm==0, field,np.nan) 
            

            # -----------------------------------------
            # area-weighted mean
            # -----------------------------------------

            # Detect coordinate names
            lat_name = "latitude" if "latitude" in field.coords else "lat"
            lon_name = "longitude" if "longitude" in field.coords else "lon"

            # -----------------------------------------
            # Handle longitude convention
            # -----------------------------------------
            
            lon_values = field[lon_name]
            
            # data uses 0-360 convention
            if lon_values.max() > 180:
            
                lon_min_sel = lon_min % 360
                lon_max_sel = lon_max % 360
            
            else:
                lon_min_sel = lon_min
                lon_max_sel = lon_max


            if lon_min_sel <= lon_max_sel:

                region = field.sel(
                    {
                        lon_name: slice(lon_min_sel, lon_max_sel),
                        lat_name: slice(lat_min, lat_max),
                    }
                )
            
            else:
                # region crosses 0/360 boundary
                region = xr.concat(
                    [
                        field.sel({lon_name: slice(lon_min_sel, 360),
                                   lat_name: slice(lat_min, lat_max)}),
                        field.sel({lon_name: slice(0, lon_max_sel),
                                   lat_name: slice(lat_min, lat_max)})
                    ],
                    dim=lon_name
                )


            # Apply Area-Weighting
            weights = np.cos(np.deg2rad(region[lat_name]))



            mean_val = region.weighted(weights).mean().item()
            mean_val_old = region.mean().item()

            mean_lookup[mode][scenario] = mean_val
            mean_lookup_old[mode][scenario] = mean_val_old


    return mean_lookup, mean_lookup_old






import numpy as np
import xarray as xr


def compute_d18O(ds, lookup, dyes=None, suffix=""):
    """
    Compute total d18O field from dye tracers.

    Parameters
    ----------
    ds : xarray.Dataset
        Must contain dye variables like:
        dye00, dye01 ... OR dye00_mean etc (via suffix)

    lookup : dict
        dye -> d18O value

    dyes : list or None
        list of dye names (default: inferred from lookup)

    suffix : str
        e.g. "", "_mean", "_std"
    """

    if dyes is None:
        dyes = list(lookup.keys())

    varnames = [f"{d}{suffix}" for d in dyes]

    missing = [v for v in varnames if v not in ds.data_vars]
    if missing:
        raise KeyError(
            f"Missing variables in dataset: {missing}\n"
            f"Available variables: {list(ds.data_vars)}"
        )

    dye_stack = xr.concat(
        [ds[v] for v in varnames],
        dim=xr.IndexVariable("dye", dyes)
    )

    weights = xr.DataArray(
        np.array([lookup[d] for d in dyes]),
        dims="dye",
        coords={"dye": dyes},
    )

    return (dye_stack * weights).sum("dye")




### Main Pipeline



from pathlib import Path
import xarray as xr
from tqdm import tqdm

from mymodules.d18O_computation import compute_d18O


def build_d18O_results(
    experiments,
    result_input,
    d18O_lookup,
    outdir,
    suffix="",
    keep_in_memory=True,
    save=True,
    overwrite=False,
    exclude=(),
):
    """
    Compute d18O fields from an existing results dictionary.

    Parameters
    ----------
    experiments : dict
        EXPERIMENTS configuration dictionary.

    result_input : dict
        Dictionary produced by run_pipeline():
            result_input[state][experiment] -> xarray.Dataset

    d18O_lookup : dict
        Dictionary:
            d18O_lookup[scenario][dye] = isotope value

    outdir : str or Path
        Directory where d18O NetCDF files are stored.

    keep_in_memory : bool
        If True, return xarray Datasets.
        If False, return filenames.

    save : bool
        Save NetCDF files.

    overwrite : bool
        Recompute existing files.

    exclude : iterable
        Experiment IDs to skip.

    Returns
    -------
    results_d18O : dict
        results_d18O[state][experiment]
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results_d18O = {}

    for state, scenarios in tqdm(experiments.items(), desc="Ocean states"):

        results_d18O[state] = {}

        for scenario, info in tqdm(
            scenarios.items(),
            desc=state,
            leave=False,
        ):

            exp = info["exp"]

            if exp in exclude:
                continue

            outfile = outdir / state / f"{scenario}_d18O.nc"

            # -------------------------------
            # Existing file
            # -------------------------------
            if outfile.exists() and not overwrite:

                if keep_in_memory:
                    results_d18O[state][exp] = xr.open_dataset(outfile)
                else:
                    results_d18O[state][exp] = outfile

                continue

            # -------------------------------
            # Compute
            # -------------------------------
            ds = result_input[state][exp]

            d18O = compute_d18O(
                ds,
                d18O_lookup[scenario],suffix=suffix
            )

            ds_out = xr.Dataset({"d18O": d18O})

            # -------------------------------
            # Save
            # -------------------------------
            if save:
                ds_out.to_netcdf(outfile)

            if keep_in_memory:
                results_d18O[state][exp] = ds_out
            else:
                results_d18O[state][exp] = outfile

    return results_d18O