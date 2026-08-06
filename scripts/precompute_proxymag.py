#!/usr/bin/env python
"""
Precompute proxy-site dye/d18O contributions (``proxymag``) for Figure 7.

Figure 7 summarises, for each marine/land palaeoclimate proxy site, the
per-dye d18O contribution under four meltwater scenarios. The extraction
(weighted spatial averaging of every dye field at every site) is the slow
part of the figure; this script runs it once and caches the result as a
small nested-dict pickle so that ``fig7_proxy_summary.py`` only has to plot.

Output
------
data/intermediates/dyestuff_modelpaper/proxymag.pkl
    Nested dict: proxymag[scenario][mode][site][dye] -> float (d18O contribution)
    plus proxymag[scenario][mode][site]["total"].
    Scenarios: 17.8k, 19.4k  (simulated, from the mean_std dye fields)
               18.2k, 20.7k  (predicted, from the mean_dye fields)
    Modes: cold, zonal.

Notes
-----
* The two extraction phases (simulated vs predicted) are merged into ONE
  proxymag dict here. In the original notebook the second phase reset
  ``proxymag = {}``, which silently discarded the simulated scenarios; merging
  fixes that.
* Surface dye fields are loaded into memory once (``.load()``) before the
  site loops, so the per-site/per-dye extraction is pure in-memory maths
  rather than repeated lazy reads.

Run from the repository root:
    python scripts/precompute_proxymag.py
"""

import os
import sys
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# --- Anchor to repo root so relative data/ paths resolve ---------------------
_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
os.chdir(_REPO)

from myconfig.DYES import DYE_TABLE
from myconfig.PROXYSITES import PROXYSITES
# Llarga Cave has no uptake mask in the shipped land_uptakemasks.pkl and is
# excluded from Figure 7 (matches the notebook's Fig7 c3 filter).
PROXYSITES = {k: v for k, v in PROXYSITES.items() if k != "Llarga Cave"}
from myconfig.EXPERIMENTS import EXPERIMENTS
from myconfig.EXPERIMENTS_prediction import EXPERIMENTS_prediction
from mymodules import d18O_computation as mf
from mymodules.dyefield_computation import run_pipeline

INTERMED = "data/intermediates/dyestuff_modelpaper"
OUT_PKL = os.path.join(INTERMED, "proxymag.pkl")

MODES = ["cold", "zonal"]
UPTAKE_MODES = {"zonal": "xqeie", "cold": "xqeic"}
SIM_SCENARIOS = ["17.8k", "19.4k"]      # simulated (mean_std dye fields)
PRED_SCENARIOS = ["18.2k", "20.7k"]     # predicted (mean_dye fields)


# =============================================================================
# d18O per-region lookup (scenario -> {dye: d18O weight})
# =============================================================================
def build_lookup():
    gdg = pd.read_pickle(os.path.join(INTERMED, "..", "regionalmeltdischarge_withd18O.pkl"))
    gdg["Newnames"] = DYE_TABLE["region"]
    old_keys = ["17.8ka", "18.2 ka", "19.4 ka", "20.7 ka"]
    new_keys = ["17.8k", "18.2k", "19.4k", "20.7k"]
    lookup = {
        nk: mf.build_d18O_lookup(gdg, column=f"mean (-35.0) region d18O anomaly {ok}")
        for ok, nk in zip(old_keys, new_keys)
    }
    return gdg, lookup


# =============================================================================
# Proxy extraction helpers (identical maths to the notebook)
# =============================================================================
def extract_ocean_proxy(field, lat, lon, size=2):
    """Latitude-weighted mean of ``field`` in a +/- size box around (lat, lon)."""
    if lon < 0:
        lon = lon + 360
    subset = field.sel(
        latitude=slice(lat - size, lat + size),
        longitude=slice(lon - size, lon + size),
    )
    weights = np.cos(np.deg2rad(subset.latitude))
    return float(subset.weighted(weights).mean(skipna=True))


def extract_land_proxy(field, uptake, area):
    """Uptake-weighted land proxy value."""
    uptake = uptake.rename({"lat": "latitude", "lon": "longitude"})
    weights = uptake * area
    weights = weights / weights.sum(skipna=True)
    return float((field * weights).sum(skipna=True))


def cell_area_2d(dye_lat, template):
    """1.25-degree grid-cell area (m^2), broadcast to the model 2-D grid."""
    R = 6.371e6
    lat_rad = np.deg2rad(dye_lat)
    dlat = np.deg2rad(1.25)
    dlon = np.deg2rad(1.25)
    ca = (R ** 2) * np.abs(np.sin(lat_rad + dlat / 2) - np.sin(lat_rad - dlat / 2)) * dlon
    ca = xr.DataArray(ca, coords=[dye_lat], dims=["latitude"])
    return ca * xr.ones_like(template)


# =============================================================================
# One extraction phase (writes into the shared proxymag dict)
# =============================================================================
def run_phase(proxymag, scenarios, field_getter, exp_getter, lookup,
              uptake_masks, area_2d):
    for scen in scenarios:
        proxymag.setdefault(scen, {})
        for mode in MODES:
            proxymag[scen][mode] = {}
            print(f"  {mode} - {scen}")
            for site, info in PROXYSITES.items():
                proxymag[scen][mode][site] = {}
                total = 0.0
                exp = exp_getter(mode, scen)
                for d in range(9):
                    dye = f"dye0{d}"
                    field = field_getter(mode, exp, d)
                    if info["type"] == "Marine Sediment Core":
                        value = extract_ocean_proxy(field, info["lat"], info["lon"])
                    elif info["type"] in ("Speleothem", "Ice Core"):
                        umask = uptake_masks[site][UPTAKE_MODES[mode]]
                        value = extract_land_proxy(field, umask, area_2d)
                    else:
                        raise ValueError(f"Unknown proxy type {info['type']}")
                    contrib = value * lookup[scen][dye]
                    proxymag[scen][mode][site][dye] = contrib
                    total += contrib
                proxymag[scen][mode][site]["total"] = total
    return proxymag


def main():
    gdg, lookup = build_lookup()

    # --- Uptake masks; remap shipped keys to site display names -------------
    with open(os.path.join(INTERMED, "land_uptakemasks.pkl"), "rb") as f:
        raw = pickle.load(f)
    uptake_masks = {
        "La Vallina (NISA)": raw["NISA_LaVallina"],
        "Cave Without a Name": raw["NonameCave"],
        "NGRIP": raw["NGRIP"],
    }

    # --- Simulated dye fields (mean_std) via the pipeline cache -------------
    print("Loading simulated dye fields (run_pipeline cache)...")
    results = run_pipeline(
        experiments=EXPERIMENTS,
        outdir=INTERMED + "/",
        keep_in_memory=True, save=True, overwrite=False,
        exclude=("xqeic", "xqeie", "xpujm"), preflight=False,
    )
    dye_lat = results["cold"]["xpujf"].latitude
    template = results["cold"]["xpujf"]["dye00_mean"].isel(depth_1=0)
    area_2d = cell_area_2d(dye_lat, template)

    # Pre-load simulated surface fields into memory: (mode, exp, d) -> DataArray
    sim_fields = {}
    for mode in MODES:
        for scen in SIM_SCENARIOS:
            exp = EXPERIMENTS[mode][scen]["exp"]
            for d in range(9):
                sim_fields[(mode, exp, d)] = (
                    results[mode][exp][f"dye0{d}_mean"].isel(depth_1=0).load()
                )

    # --- Predicted dye fields (mean_dye) ------------------------------------
    print("Loading predicted dye fields (mean_dye)...")
    def _load_pred(fname):
        ds = xr.open_dataset(os.path.join(INTERMED, fname)).load()
        # Shipped mean_dye files carry lat/lon; rename to match the simulated
        # grid (latitude/longitude) and area_2d used by the extractors.
        rename = {k: v for k, v in {"lat": "latitude", "lon": "longitude"}.items() if k in ds.dims}
        return ds.rename(rename) if rename else ds

    dict_dye = {
        "zonal": _load_pred("mean_dye_zonal.nc"),
        "cold": _load_pred("mean_dye_cold.nc"),
    }
    pred_fields = {}
    for mode in MODES:
        for d in range(9):
            pred_fields[(mode, mode, d)] = dict_dye[mode][f"dye0{d}"]

    proxymag = {}

    print("Phase 1: simulated scenarios", SIM_SCENARIOS)
    run_phase(
        proxymag, SIM_SCENARIOS,
        field_getter=lambda mode, exp, d: sim_fields[(mode, exp, d)],
        exp_getter=lambda mode, scen: EXPERIMENTS[mode][scen]["exp"],
        lookup=lookup, uptake_masks=uptake_masks, area_2d=area_2d,
    )

    print("Phase 2: predicted scenarios", PRED_SCENARIOS)
    run_phase(
        proxymag, PRED_SCENARIOS,
        field_getter=lambda mode, exp, d: pred_fields[(mode, mode, d)],
        exp_getter=lambda mode, scen: mode,   # predicted fields keyed by mode
        lookup=lookup, uptake_masks=uptake_masks, area_2d=area_2d,
    )

    with open(OUT_PKL, "wb") as f:
        pickle.dump(proxymag, f, protocol=pickle.HIGHEST_PROTOCOL)
    size_kb = os.path.getsize(OUT_PKL) / 1024
    n = sum(len(proxymag[s][m]) for s in proxymag for m in proxymag[s])
    print(f"\nWrote {OUT_PKL} ({size_kb:.1f} KB) — "
          f"{len(proxymag)} scenarios x {len(MODES)} modes, {n} site-entries.")


if __name__ == "__main__":
    main()
