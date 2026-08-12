#!/usr/bin/env bash
# =============================================================================
# concat_msl_daily.sh
#
# Rebuild a continuous DAILY mean-sea-level-pressure time series from the
# per-month UM output files, for the storm-track precompute step.
#
# Each input file (e.g.  xpraka#pa000005300mr+.nc ) holds one 30-day month of
# daily-mean MSL (variable p_dm_msl, 360_day calendar) with its OWN time_origin.
# We decode each file's calendar and rebuild ONE monotonic daily axis -- so we
# use `cdo mergetime` (or an xarray fallback), NOT `ncrcat` (ncrcat keeps only
# the first file's time units and would corrupt the axis).
#
# Usage:
#   ./concat_msl_daily.sh                 # does xpraj and xprak (defaults)
#   ./concat_msl_daily.sh xprak           # one experiment
#   INROOT=/some/where OUTDIR=/out ./concat_msl_daily.sh xpraj xprak
#
# Env overrides:
#   INROOT   root holding <exp>/pa/         (default /nfs/annie/zgcq387/work/um)
#   SUBDIR   stream subfolder under <exp>    (default pa)
#   OUTDIR   where to write <exp>.mslp.monthly.nc (default current dir)
#
# Output:  <OUTDIR>/<exp>.mslp.monthly.nc   (daily timesteps; the "monthly" tag
#          matches the filename convention the precompute script expects)
# =============================================================================
set -euo pipefail

INROOT="${INROOT:-/nfs/annie/zgcq387/work/um}"
SUBDIR="${SUBDIR:-pa}"
OUTDIR="${OUTDIR:-.}"
EXPS=("$@"); [ ${#EXPS[@]} -eq 0 ] && EXPS=(xpraj xprak)

mkdir -p "$OUTDIR"

# pick a merge backend once
BACKEND=""
if command -v cdo >/dev/null 2>&1; then
    BACKEND="cdo"
elif python -c "import xarray, cftime" >/dev/null 2>&1; then
    BACKEND="xarray"
else
    echo "ERROR: need either 'cdo' on PATH or python with xarray+cftime." >&2
    exit 1
fi
echo "merge backend: $BACKEND"

for exp in "${EXPS[@]}"; do
    dir="$INROOT/$exp/$SUBDIR"
    out="$OUTDIR/${exp}.mslp.monthly.nc"
    echo "=============================================================="
    echo "[$exp]  input dir: $dir"

    # gather monthly files in chronological (= zero-padded lexicographic) order
    shopt -s nullglob
    files=( "$dir/${exp}a#pa"*".nc" )
    shopt -u nullglob
    if [ ${#files[@]} -eq 0 ]; then
        echo "  no files matching ${exp}a#pa*.nc in $dir -- skipping" >&2
        continue
    fi
    # sort defensively (glob is already sorted, but be explicit)
    IFS=$'\n' files=( $(printf '%s\n' "${files[@]}" | sort) ); unset IFS
    echo "  found ${#files[@]} monthly files"
    echo "  first: $(basename "${files[0]}")"
    echo "  last : $(basename "${files[${#files[@]}-1]}")"

    # --- merge, rebuilding one consistent daily axis ------------------------
    if [ "$BACKEND" = "cdo" ]; then
        cdo -s -O mergetime "${files[@]}" "$out"
    else
        python - "$out" "${files[@]}" <<'PY'
import sys, xarray as xr
out, files = sys.argv[1], sys.argv[2:]
# decode each file's own 360_day origin -> cftime, then concat along t
ds = xr.open_mfdataset(files, combine="nested", concat_dim="t",
                       use_cftime=True, decode_times=True,
                       data_vars="minimal", coords="minimal", compat="override")
ds = ds.sortby("t")
enc = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
ds.to_netcdf(out, encoding=enc)
ds.close()
PY
    fi

    # --- verify the result is a clean daily axis ----------------------------
    python - "$out" <<'PY'
import sys, numpy as np, xarray as xr
ds = xr.open_dataset(sys.argv[1], use_cftime=True)
t = ds["t"].values
gap_days = np.median([ (t[i+1]-t[i]).total_seconds()/86400.0
                       if hasattr(t[i+1]-t[i], "total_seconds")
                       else (t[i+1]-t[i])/np.timedelta64(1,"D")
                       for i in range(len(t)-1) ])
print(f"  -> wrote {sys.argv[1]}")
print(f"     timesteps : {t.size}")
print(f"     span      : {t[0]}  ..  {t[-1]}")
print(f"     median dt : {gap_days:.2f} day(s)   (should be ~1.0)")
ds.close()
PY
done

echo "=============================================================="
echo "done. point precompute at these via --msl-base '$OUTDIR/{exp}.mslp.monthly.nc'"
