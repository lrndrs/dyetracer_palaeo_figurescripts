#!/usr/bin/env python3
"""
Precompute light atmospheric fields for the AMOC-mode atmosphere figures
(Surface Air Temperature, 850 hPa wind, storm-track intensity).

Reduces heavy HadCM3 time-series (esp. daily MSL) to small 2-D fields on the
atmospheric grid, tagged by AMOC mode (cold / zonal / merid), and writes them
to  data/intermediates/atmo/ .

The storm-track step is STREAMED directly from the raw per-month daily-MSL
files (thousands of them, ~500 model-years): each chunk is band-passed in
memory and its variance accumulated, so the full multi-GB daily record is
never concatenated or saved. The only storm-track output is the small 2-D
band-pass-variance climatology (annual + DJFM), one field per AMOC mode.

Run this ON foe-linux (where the raw /nfs/annie/... time-series live):

    conda activate <your dyetracer-like env>   # needs xarray, scipy, netCDF4, (dask optional)
    python precompute_atmo.py --out /path/to/dyetracer_data_bundle/data/intermediates/atmo

By default the storm track reads per-month MSL files from
  --msl-dir  /nfs/annie/zgcq387/work/um/{exp}/pa
  --msl-glob {exp}a#pa*.nc
Override --chunk-files (default 60 = 5 model-years/chunk) to trade RAM vs
edge losses. If you already concatenated a single daily file, pass
  --msl-single /path/to/{exp}.mslp.monthly.nc
to use it instead of streaming.

The output (~2-3 MB total) is what the plotting script ships with.
"""
import argparse, os, warnings
import numpy as np
import xarray as xr
import scipy.signal as signal

# xarray's use_cftime kwarg is deprecated but still works; silence the notice
warnings.filterwarnings("ignore", message=".*use_cftime.*")

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

def _list_msl_files(msl_dir, exp, glob_pat):
    """Sorted (= chronological, thanks to zero-padded step numbers) list of the
    raw per-month daily-MSL files for one experiment."""
    import glob as _glob
    pat = os.path.join(msl_dir.format(exp=exp), glob_pat.format(exp=exp))
    files = sorted(_glob.glob(pat))
    return files

def stormtrack_climatology_streaming(files, chunk_files=60, trim=15,
                                     djfm_months=(12, 1, 2, 3), verbose=True):
    """Storm-track band-pass variance CLIMATOLOGY, streamed from per-month files.

    Reads the raw monthly files in contiguous chunks (default 60 files = 5
    model-years = 1800 continuous days), band-passes each chunk in memory, and
    accumulates the sum-of-squares of the band-passed field for the annual and
    DJFM masks. Never materialises the full daily record -> peak memory is one
    chunk (~100 MB), and the only output is the small 2-D climatology.

    band-pass FIRST, then mask DJFM (the correct order): the continuous record
    is filtered, then the DJFM days are selected -- so no spurious energy is
    injected at season seams. `trim` days are dropped from each chunk edge to
    remove filtfilt transients.

    Returns (annual_var, djfm_var) as 2-D DataArrays (latitude, longitude), hPa^2.
    """
    sum_ann = cnt_ann = sum_djf = cnt_djf = None
    lat = lon = None
    n = len(files)
    checked_daily = False
    for i0 in range(0, n, chunk_files):
        grp = files[i0:i0 + chunk_files]
        ds = xr.open_mfdataset(grp, combine="nested", concat_dim="t",
                               use_cftime=True, data_vars="minimal",
                               coords="minimal", compat="override").sortby("t")
        if not checked_daily:
            _assert_daily(ds); checked_daily = True
        da = _pick_msl_var(ds) / 100.0                     # Pa -> hPa
        for d in list(da.dims):                            # squeeze singleton level
            if d not in ("t", "latitude", "longitude") and da.sizes[d] == 1:
                da = da.isel({d: 0}, drop=True)
        da = da.load()
        months = da["t"].dt.month.values
        vals = da.values                                   # (t, lat, lon)
        anom = vals - np.nanmean(vals, axis=0, keepdims=True)
        filt = _bandpass(anom)                             # continuous band-pass
        if trim > 0 and filt.shape[0] > 2 * trim:
            sl = slice(trim, filt.shape[0] - trim)
        else:
            sl = slice(None)
        f, m = filt[sl], months[sl]
        sq = f * f
        valid = ~np.isnan(sq)
        if sum_ann is None:
            sum_ann = np.nansum(sq, axis=0); cnt_ann = valid.sum(axis=0)
            sum_djf = np.zeros_like(sum_ann); cnt_djf = np.zeros_like(cnt_ann)
            lat, lon = da["latitude"], da["longitude"]
        else:
            sum_ann += np.nansum(sq, axis=0); cnt_ann += valid.sum(axis=0)
        dmask = np.isin(m, djfm_months)
        if dmask.any():
            sqd = sq[dmask]
            sum_djf += np.nansum(sqd, axis=0)
            cnt_djf += (~np.isnan(sqd)).sum(axis=0)
        if verbose:
            print(f"       chunk {i0//chunk_files+1}: files "
                  f"{i0+1}-{min(i0+chunk_files,n)} of {n}", flush=True)
        ds.close()

    var_ann = xr.DataArray(sum_ann / np.where(cnt_ann == 0, np.nan, cnt_ann),
                           coords={"latitude": lat, "longitude": lon},
                           dims=("latitude", "longitude"), name="annual")
    var_djf = xr.DataArray(sum_djf / np.where(cnt_djf == 0, np.nan, cnt_djf),
                           coords={"latitude": lat, "longitude": lon},
                           dims=("latitude", "longitude"), name="djfm")
    return var_ann, var_djf

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
    ap.add_argument("--msl-dir", default="/nfs/annie/zgcq387/work/um/{exp}/pa",
                    help="dir holding the raw per-month daily-MSL files "
                         "({exp} placeholder). Storm track is streamed from "
                         "these -- no concatenated daily file needed.")
    ap.add_argument("--msl-glob", default="{exp}a#pa*.nc",
                    help="glob for the per-month MSL files within --msl-dir")
    ap.add_argument("--chunk-files", type=int, default=60,
                    help="per-month files band-passed per streaming chunk "
                         "(60 = 5 model-years; larger = fewer edge losses, more RAM)")
    ap.add_argument("--trim", type=int, default=15,
                    help="days dropped from each chunk edge (filter transients)")
    ap.add_argument("--msl-single", default=None,
                    help="fallback: a single pre-concatenated daily MSL file "
                         "template ({exp}); use instead of streaming --msl-dir")
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

        # --- Storm track: band-pass MSL variance CLIMATOLOGY ----------------
        # Streamed directly from the raw per-month daily-MSL files: band-pass
        # each chunk in memory and accumulate variance. No giant concatenated
        # daily file is ever built or saved -- the only output is the small 2-D
        # climatology below. (band-pass FIRST, then mask DJFM = correct order.)
        if args.msl_single:
            ds_slp = xr.open_dataset(args.msl_single.format(exp=exp), use_cftime=True)
            dt = _assert_daily(ds_slp)
            print(f"       MSL timestep ~{dt:.1f} day(s) (single file)", flush=True)
            st_ann[mode]  = storm_track_variance(ds_slp)
            st_djfm[mode] = storm_track_variance(ds_slp, months=[12, 1, 2, 3])
        else:
            files = _list_msl_files(args.msl_dir, exp, args.msl_glob)
            if not files:
                raise FileNotFoundError(
                    f"no MSL files at {args.msl_dir.format(exp=exp)}/"
                    f"{args.msl_glob.format(exp=exp)}")
            print(f"       streaming storm track from {len(files)} monthly "
                  f"files ({args.chunk_files}/chunk)", flush=True)
            st_ann[mode], st_djfm[mode] = stormtrack_climatology_streaming(
                files, chunk_files=args.chunk_files, trim=args.trim)

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
