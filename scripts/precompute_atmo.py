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

def storm_track_variance(ds_slp, months=None):
    """2-6 day band-pass variance of MSL (hPa^2). months=[12,1,2,3] for DJFM."""
    da = ds_slp.p_dm_msl / 100.0            # Pa -> hPa
    if months is not None:
        da = da.sel(t=da.t.dt.month.isin(months))
    anom = da - da.mean(dim="t")
    filt = xr.apply_ufunc(
        _bandpass, anom,
        input_core_dims=[["t"]], output_core_dims=[["t"]],
        vectorize=True, dask="parallelized", keep_attrs=True,
    )
    return filt.var(dim="t")                # (msl, lat, lon)

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
        ds_slp = xr.open_dataset(p("msl.daily"))
        st_ann[mode]  = storm_track_variance(ds_slp).isel(msl=0)
        st_djfm[mode] = storm_track_variance(ds_slp, months=[12, 1, 2, 3]).isel(msl=0)

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
