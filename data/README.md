# Plotting data — manifest and provenance

This folder holds the light, plotting-ready data for the
`dyetracer_palaeo_figurescripts` figure scripts. **The data payload is not stored in
git** (`.gitignore`: `data/*`) — download it from the Zenodo record and unpack it here.

> **Zenodo DOI:** {{ZENODO_DOI}}

Total ~45 MB, reduced from ~735 MB of raw HadCM3 output — only the fields needed for
plotting are included.

---

## Expected layout after unpacking

```
data/
├── inputs/
│   ├── temev.qrparm.omask.nc            # HadCM3 land–sea mask (144×288)
│   ├── qrparm.waterfix.hadcm3.nc        # ocean waterfix field (field672)
│   ├── discharge_ts.csv                 # meltwater discharge time series (precomputed; Fig. 1)
│   └── ice_sheet_extent/                # LGM ice-sheet extent shapefile (8 files)
│       └── LGM_best_estimate.{shp,shx,dbf,prj,cpg,sbn,sbx,shp.xml}
├── intermediates/
│   ├── seaice_lookup.pkl                # sea-ice cache, keyed cold/merid/zonal (Fig. 1)
│   ├── mld_lookup.pkl                   # mixed-layer-depth cache (Fig. 1)
│   ├── amoc_lookup.pkl                  # AMOC-strength cache (Fig. 1)
│   ├── gdf_regions.pkl                  # 9 source-region polygons (GeoDataFrame)
│   ├── regionalmeltdischarge_withd18O.pkl   # per-region melt + d18O table (GeoDataFrame)
│   └── dyestuff_modelpaper/
│       ├── dye_regions_norm.nc          # normalized dye-region field
│       ├── mean_dye_{cold,merid,zonal}.nc   # dye-mean fields per AMOC mode (9 dyes)
│       ├── land_uptakemasks.pkl         # land proxy-site uptake masks (Fig. 7)
│       ├── proxymag.pkl                 # precomputed proxy-site contributions (Fig. 7)
│       └── {cold,merid,zonal}/{scenario}_mean_std.nc   # surface mean/std fields
└── trajectories/
    └── {NGRIP,NISA_LaVallina,NonameCave}_{xqeic,xqeie}_th00_UTOT_weighted.nc  # Fig. 6
```

---

## Contents by category

### `inputs/` — model input grids (~2.9 MB)
Model grid definitions and a precomputed discharge time series.
- `temev.qrparm.omask.nc` — land–sea mask on the HadCM3 ocean grid (latitude 144 × longitude 288).
- `qrparm.waterfix.hadcm3.nc` — ocean waterfix (`field672`), 144 × 290.
- `discharge_ts.csv` — 261-row meltwater discharge time series (columns `t, elwg, gin, arc, tot`; `t` = −26000…0 yr; `tot` ≈ 0.04–0.40 Sv). Precomputed offline (originally from `mw_protocol.plotting.create_discharge_ts`) so Fig. 1 runs without the external package.
- `ice_sheet_extent/LGM_best_estimate.*` — LGM ice-sheet extent shapefile (913 polygons, EPSG:4326); after Batchelor et al. (2019).

### `intermediates/` — precomputed caches (~36 MB)
Reduced products so the analysis pipeline short-circuits to loading cached results
instead of re-reducing raw model output.
- `seaice_lookup.pkl`, `mld_lookup.pkl`, `amoc_lookup.pkl` — per-mode (cold/merid/zonal) caches used by Fig. 1.
- `gdf_regions.pkl` — 9 dye source-region polygons (GeoDataFrame, 9 × 7).
- `regionalmeltdischarge_withd18O.pkl` — per-region melt and d18O table (GeoDataFrame, 9 × 30).
- `dyestuff_modelpaper/dye_regions_norm.nc` — normalized dye-region field.
- `dyestuff_modelpaper/mean_dye_{cold,merid,zonal}.nc` — dye-mean surface fields (9 dyes, means only) per AMOC mode.
- `dyestuff_modelpaper/land_uptakemasks.pkl` — land proxy-site uptake masks (keys: `NISA_LaVallina`, `NonameCave`, `NGRIP`).
- `dyestuff_modelpaper/proxymag.pkl` — precomputed per-site, per-dye d18O contributions (built by `scripts/precompute_proxymag.py`).
- `dyestuff_modelpaper/{cold,merid,zonal}/{17.8k,18.2k,19.4k,20.7k}_mean_std.nc` — **surface-only** per-scenario dye mean/std fields (18 variables `dye00_mean/std … dye08_mean/std`).

### `trajectories/` — atmospheric back-trajectories (~6 MB)
The six trajectory files actually plotted in Fig. 6 (of 24 total): three proxy locations
(`NGRIP`, `NISA_LaVallina`, `NonameCave`) × two high-resolution experiments
(`xqeic` = cold, `xqeie` = zonal), threshold `th00`, weighted by total moisture uptake
(`UTOT_weighted`).

---

## Runtime-generated files

The figure scripts write additional per-scenario d18O files under
`intermediates/dyestuff_modelpaper/` (top level and `min/`, `max/` subfolders) as
`{scenario}_d18O.nc`. These are **computed at run time** by
`mymodules.d18O_computation.build_d18O_results` from the shipped mean/std and dye-mean
fields — they are not part of the download and can be deleted freely; the scripts
recreate them.

---

## Provenance notes

- **Cached pickles were generated at a slightly newer upstream code revision than the
  original notebooks.** Two keys therefore differ from the notebook code — proxy-site
  names in `land_uptakemasks.pkl` and the discharge-scenario keys — and the shipped
  scripts have been adapted to the shipped data (the data is the ground truth for these
  figures).
- **`intermediates/dyestuff_modelpaper/merid/18.2k_mean_std.nc` is intentionally
  absent.** Its experiment (xpujm) ran only 245 of the required 500 years, so the
  years 450–499 statistics window is entirely missing (all-NaN); no figure uses it.
- **Per-scenario mean/std fields are surface-only** (top depth level) to keep the bundle
  light. This is all the figures require.
- **Proxy-site coverage.** Fig. 7 uses six marine sediment cores, three speleothems
  (La Vallina / NISA, Cave Without a Name, NGRIP as the ice core), and drops "Llarga
  Cave", which has no uptake mask in `land_uptakemasks.pkl`.

## Source / citation

Base climate simulations: HadCM3 (BRIDGE) with GLAC-1D meltwater forcing, following
Rome et al. (2022). The meltwater-routing toolbox (`mw_protocol`) is external:
Olnavy (2022), Zenodo, https://doi.org/10.5281/zenodo.6788389.
