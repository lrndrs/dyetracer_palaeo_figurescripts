#!/usr/bin/env python3
"""
Precompute light atmospheric fields for the AMOC-mode atmosphere figures
(Surface Air Temperature, 850 hPa wind, storm-track intensity).

Reduces heavy HadCM3 time-series (esp. daily MSL) to small 2-D fields on the
atmospheric grid, tagged by AMOC mode (cold / zonal / merid), and writes them
to  data/intermediates/atmo/ .

Run this ON foe-linux (where the raw /nfs/annie/... time-series live):

    conda activate <your dyetracer-like env>   # needs xarray, scipy, netCDF4, (dask optional)
    python precompute_atmo.py --out /path/to/dyetracer_data_bundle/data/intermediates/atmo

The output (~2-3 MB total) is what the plotting script ships with.
"""
import argparse, os
import numpy as np
import xarray as xr
import scipy.signal as signal

BASE = "/nfs/annie/earpal/database/experiments/{exp}/time_series/{exp}.{field}.nc"

# AMOC mode  ->  experiment id  (all 17.8 ka)
MODES = {"cold": "xpraj", "zonal": "xprak", "merid": "xpral"}

# ---- storm-track band-pass filter (2-6 day synoptic band) -------------------
DT = 1.0                      # daily data
FS = 1.0 / DT
NYQ = 0.5 * FS
LOWCUT, HIGHCUT = 1/6, 1/2    # 6-day (low) .. 2-day (high)
LOW = LOWCUT / NYQ
HIGH = min(HIGHCUT / NYQ, 0.99)
ORDER = 3

def _bandpass(data, low=LOW, high=HIGH, order=ORDER):
    if np.all(np.isnan(data)):
        return data
    b, a = signal.butter(order, [low, high], btype="band")
    return signal.filtfilt(b, a, np.asarray(data), axis=0)

def _pick_msl_var(ds):
    """Return the MSL DataArray. Prefer the notebook's `p_dm_msl`; otherwise
    fall back to the single time-varying pressure field in the file."""
    if "p_dm_msl" in ds.data_vars:
        return ds["p_dm_msl"]
    cands = [v for v in ds.data_vars
             if "t" in ds[v].dims and ("msl" in v.lower() or "slp" in v.lower()
                                       or "pmsl" in v.lower())]
    if len(cands) == 1:
        return ds[cands[0]]
    raise KeyError(
        f"Could not identify the MSL variable in {list(ds.data_vars)}. "
        "Expected `p_dm_msl`; edit _pick_msl_var() to name the right field.")

def _timedelta_to_days(td):
    """Gap between two timesteps -> days, for numpy.timedelta64, python
    datetime.timedelta, or cftime timedeltas (non-standard calendars)."""
    if isinstance(td, np.timedelta64):
        return td / np.timedelta64(1, "D")
    if hasattr(td, "total_seconds"):          # datetime.timedelta / cftime
        return td.total_seconds() / 86400.0
    if hasattr(td, "days"):
        return td.days + getattr(td, "seconds", 0) / 86400.0
    return float(td)

def _assert_daily(ds):
    """Storm-track band-pass is only meaningful on sub-monthly timesteps.
    Stop loudly if the file turns out to hold true monthly means."""
    t = ds["t"].values
    if t.size < 4:
        raise ValueError(f"MSL file has only {t.size} timesteps; need daily data.")
    gaps = [_timedelta_to_days(d) for d in np.diff(t)]
    dt_days = float(np.median(gaps))
    if dt_days > 20:
        raise ValueError(
            f"MSL timestep is ~{float(dt_days):.0f} days — this looks like monthly "
            "means, not daily data. The 2-6 day storm-track band-pass requires "
            "daily timesteps and would be meaningless here. Point --base at a file "
            "with daily resolution.")
    return float(dt_days)

def storm_track_variance(ds_slp, months=None):
    """2-6 day band-pass variance of MSL (hPa^2). months=[12,1,2,3] for DJFM."""
    da = _pick_msl_var(ds_slp) / 100.0      # Pa -> hPa
    if months is not None:
        da = da.sel(t=da.t.dt.month.isin(months))
    anom = da - da.mean(dim="t")
    filt = xr.apply_ufunc(
        _bandpass, anom,
        input_core_dims=[["t"]], output_core_dims=[["t"]],
        vectorize=True, dask="parallelized", keep_attrs=True,
    )
    var = filt.var(dim="t")                 # (msl?, lat, lon)
    # drop any singleton vertical/level dim so output is 2-D (lat, lon)
    for d in list(var.dims):
        if d not in ("latitude", "longitude") and var.sizes[d] == 1:
            var = var.isel({d: 0}, drop=True)
    return var

def _concat_modes(per_mode):
    """Stack a dict {mode: DataArray} along a new 'mode' coord."""
    modes = list(per_mode.keys())
    out = xr.concat([per_mode[m] for m in modes], dim="mode")
    out = out.assign_coords(mode=("mode", modes))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/intermediates/atmo",
                    help="output directory for the light atmo bundle")
    ap.add_argument("--base", default=BASE,
                    help="template path with {exp} and {field}")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sat, u850, v850, st_ann, st_djfm, seaice = {}, {}, {}, {}, {}, {}

    for mode, exp in MODES.items():
        print(f"[{mode}] {exp}", flush=True)
        p = lambda f: args.base.format(exp=exp, field=f)

        # --- Surface Air Temperature: time-mean surface field (kept in K) ----
        ds_t = xr.open_dataset(p("tempsurf.monthly"))
        sat[mode] = ds_t["temp_mm_srf"].mean(dim="t").isel(surface=0)

        # --- 850 hPa wind: time-mean u,v at p=0 ------------------------------
        ds_u = xr.open_dataset(p("u850.monthly"))
        ds_v = xr.open_dataset(p("v850.monthly"))
        u850[mode] = ds_u["u_mm_p"].mean(dim="t").isel(p=0)
        v850[mode] = ds_v["v_mm_p"].mean(dim="t").isel(p=0)

        # --- Storm track: band-pass MSL variance, annual + DJFM --------------
        # NB: the daily MSL time series is stored in a file *named* .mslp.monthly
        # (the "monthly" tag is a naming convention; the timesteps are daily).
        ds_slp = xr.open_dataset(p("mslp.monthly"))
        dt = _assert_daily(ds_slp)          # stop if this is truly monthly means
        print(f"       MSL timestep ~{dt:.1f} day(s)", flush=True)
        st_ann[mode]  = storm_track_variance(ds_slp)
        st_djfm[mode] = storm_track_variance(ds_slp, months=[12, 1, 2, 3])

        # --- Sea-ice monthly climatology (for 50% extent overlay) ------------
        ds_ice = xr.open_dataset(p("iceconc.monthly"))
        seaice[mode] = ds_ice["iceconc_mm_srf"].groupby("t.month").mean("t").isel(surface=0)

    # ---- assemble & write -------------------------------------------------
    xr.Dataset({"temp_mm_srf": _concat_modes(sat)}).to_netcdf(
        os.path.join(args.out, "sat_mean.nc"))
    xr.Dataset({"u": _concat_modes(u850), "v": _concat_modes(v850)}).to_netcdf(
        os.path.join(args.out, "wind850_mean.nc"))
    xr.Dataset({"annual": _concat_modes(st_ann), "djfm": _concat_modes(st_djfm)}).to_netcdf(
        os.path.join(args.out, "stormtrack.nc"))
    xr.Dataset({"iceconc": _concat_modes(seaice)}).to_netcdf(
        os.path.join(args.out, "seaice_monthly.nc"))

    print("wrote:", sorted(os.listdir(args.out)), flush=True)

if __name__ == "__main__":
    main()
