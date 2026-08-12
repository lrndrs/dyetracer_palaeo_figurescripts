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
# STAGED MERGE: with thousands of monthly files, a single `cdo mergetime *`
# opens every file at once and dies with "Too many open files". Instead we
# merge in BATCHES (default 100 files) into intermediates, then merge the
# intermediates, repeating until one file remains -- so only ~BATCH files are
# ever open simultaneously. Works for any file count.
#
# Usage:
#   ./concat_msl_daily.sh                 # does xpraj and xprak (defaults)
#   ./concat_msl_daily.sh xprak           # one experiment
#   INROOT=/some/where OUTDIR=/out ./concat_msl_daily.sh xpraj xprak
#
# Env overrides:
#   INROOT   root holding <exp>/pa/          (default /nfs/annie/zgcq387/work/um)
#   SUBDIR   stream subfolder under <exp>     (default pa)
#   OUTDIR   where to write <exp>.mslp.monthly.nc (default current dir)
#   BATCH    files merged per stage           (default 100)
#   KEEPTMP  set to 1 to keep intermediates    (default: clean up)
#
# Output:  <OUTDIR>/<exp>.mslp.monthly.nc    (daily timesteps; the "monthly" tag
#          matches the filename convention the precompute script expects)
# =============================================================================
set -euo pipefail

INROOT="${INROOT:-/nfs/annie/zgcq387/work/um}"
SUBDIR="${SUBDIR:-pa}"
OUTDIR="${OUTDIR:-.}"
BATCH="${BATCH:-100}"
EXPS=("$@"); [ ${#EXPS[@]} -eq 0 ] && EXPS=(xpraj xprak)

mkdir -p "$OUTDIR"

# best-effort raise of the open-file limit (belt-and-suspenders; batching is
# what actually fixes the problem, but a higher limit never hurts)
ulimit -n 4096 2>/dev/null || true

# pick a merge backend once
if command -v cdo >/dev/null 2>&1; then
    BACKEND="cdo"
elif python -c "import xarray, cftime" >/dev/null 2>&1; then
    BACKEND="xarray"
else
    echo "ERROR: need either 'cdo' on PATH or python with xarray+cftime." >&2
    exit 1
fi
echo "merge backend: $BACKEND   (batch size: $BATCH)"

# merge_group OUT IN1 IN2 ...  ->  concatenate along time into OUT
merge_group() {
    local out="$1"; shift
    if [ "$BACKEND" = "cdo" ]; then
        cdo -s -O mergetime "$@" "$out"
    else
        python - "$out" "$@" <<'PY'
import sys, xarray as xr
out, files = sys.argv[1], sys.argv[2:]
ds = xr.open_mfdataset(files, combine="nested", concat_dim="t",
                       use_cftime=True, data_vars="minimal",
                       coords="minimal", compat="override")
ds = ds.sortby("t")
ds.to_netcdf(out)
ds.close()
PY
    fi
}

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
    IFS=$'\n' files=( $(printf '%s\n' "${files[@]}" | sort) ); unset IFS
    echo "  found ${#files[@]} monthly files"
    echo "  first: $(basename "${files[0]}")"
    echo "  last : $(basename "${files[${#files[@]}-1]}")"

    tmp="$OUTDIR/.mergetmp_${exp}"
    rm -rf "$tmp"; mkdir -p "$tmp"

    # --- staged reduce: merge in batches until one file remains -------------
    level=0
    cur=( "${files[@]}" )
    while [ ${#cur[@]} -gt "$BATCH" ]; do
        next=(); i=0; b=0
        while [ $i -lt ${#cur[@]} ]; do
            grp=( "${cur[@]:$i:$BATCH}" )
            interm="$tmp/L${level}_$(printf '%04d' $b).nc"
            merge_group "$interm" "${grp[@]}"
            next+=( "$interm" )
            i=$(( i + BATCH )); b=$(( b + 1 ))
        done
        echo "  stage $level: ${#cur[@]} files -> ${#next[@]} intermediates"
        # drop the previous level's intermediates (keep original inputs)
        [ $level -gt 0 ] && rm -f "$tmp/L$((level-1))_"*.nc
        cur=( "${next[@]}" ); level=$(( level + 1 ))
    done

    # final merge of the remaining <= BATCH files
    merge_group "$out" "${cur[@]}"
    echo "  final merge: ${#cur[@]} files -> $(basename "$out")"
    [ "${KEEPTMP:-0}" = "1" ] || rm -rf "$tmp"

    # --- verify the result is a clean daily axis ----------------------------
    python - "$out" <<'PY'
import sys, numpy as np, xarray as xr
ds = xr.open_dataset(sys.argv[1], use_cftime=True)
t = ds["t"].values
diffs = [t[i+1]-t[i] for i in range(len(t)-1)]
gap = np.median([d.total_seconds()/86400.0 if hasattr(d,"total_seconds")
                 else d/np.timedelta64(1,"D") for d in diffs])
# duplicate-timestep check (overlapping files would show dt=0 somewhere)
ndup = int(np.sum([ (d.total_seconds() if hasattr(d,"total_seconds")
                     else d/np.timedelta64(1,"s"))==0 for d in diffs]))
print(f"  -> wrote {sys.argv[1]}")
print(f"     timesteps : {t.size}")
print(f"     span      : {t[0]}  ..  {t[-1]}")
print(f"     median dt : {gap:.2f} day(s)   (should be ~1.0)")
print(f"     dup steps : {ndup}   (should be 0)")
ds.close()
PY
done

echo "=============================================================="
echo "done. point precompute at these via --msl-base '$OUTDIR/{exp}.mslp.monthly.nc'"
