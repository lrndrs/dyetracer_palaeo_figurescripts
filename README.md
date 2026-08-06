# dyetracer_palaeo_figurescripts

Figure-generating scripts and light plotting data for the manuscript

> **Tracing Meltwater from Northern Ice Sheets to Palaeoclimate Archives During the
> Early Last Deglaciation: a conservative tracer approach**
> L. Endres, R. Ivanovic, Y. Rome, J. Tindall, I. Thurnherr, M. Sprenger, H. Stoll
> (submitted to *Paleoceanography and Paleoclimatology*, AGU).

Each script reproduces one main-text figure (Fig. 1–7) from a small, plotting-ready
subset of the HadCM3 dye-tracer model output. The heavy raw model data have been
reduced offline to just the fields needed for plotting (~45 MB total).

---

## Repository layout

```
dyetracer_palaeo_figurescripts/
├── scripts/                     # one runnable .py per figure
│   ├── precompute_proxymag.py       # (run first for Fig. 7) caches proxy-site contributions
│   ├── fig1_maps_amocmode.py        # Fig. 1  AMOC modes + discharge variants
│   ├── fig2_dye_regions.py          # Fig. 2  dye source regions
│   ├── fig3_dye_propagation.py      # Fig. 3  dye propagation / scenario maps
│   ├── fig4_fig5_scenario_maps.py   # Fig. 4 & 5  GIN/IRD scenario maps, ratios, scatter
│   ├── fig6_trajectories.py         # Fig. 6  atmospheric back-trajectories
│   └── fig7_proxy_summary.py        # Fig. 7  proxy-site d18O contribution summary
├── mymodules/                   # analysis helpers (see "Bundled code" below)
│   ├── dyefield_computation.py      # dye-field / lookup pipeline (run_pipeline etc.)
│   ├── d18O_computation.py          # d18O results + area-weighted means
│   ├── grid_utils.py                # grid-edge helper (see "External code" below)
│   └── myfunc.py                    # misc numerical helpers
├── myconfig/                    # experiment tables + site definitions
│   ├── EXPERIMENTS.py                # AMOC-mode experiment IDs (cold/merid/zonal)
│   ├── EXPERIMENTS_prediction.py     # predicted-scenario mapping
│   ├── DYES.py                       # 9-region dye table
│   └── PROXYSITES.py                 # marine-core / speleothem / ice-core coordinates
├── data/                        # plotting data — NOT in git; download from Zenodo
│   └── README.md                    # data manifest + provenance
├── figures/                     # created at run time; holds the output PDFs
├── environment.yml              # conda environment (recommended)
├── requirements.txt             # pip alternative
└── README.md
```

---

## Setup

### 1. Environment

**Conda (recommended** — resolves the compiled geospatial stack cleanly**):**

```bash
conda env create -f environment.yml
conda activate dyetracer
```

**pip alternative:**

```bash
pip install -r requirements.txt
```

Tested with Python 3.11 on conda-forge. `cartopy` and `geopandas` pull in GEOS/PROJ/GDAL;
if pip cannot build them on your platform, use the conda route.

### 2. Data

The plotting data live on Zenodo, not in this git repository (see **Data** below).
Download the data archive and unpack it so that the `data/` folder sits at the repo root:

```bash
# from the repo root, after downloading the Zenodo data archive:
tar -xzf dyetracer_data_bundle.tar.gz      # creates ./data/...
```

The scripts expect data under `data/inputs/`, `data/intermediates/`, and
`data/trajectories/` (relative to the repo root — the scripts `chdir` to the repo root
automatically).

### 3. Cartopy coastlines (first run only)

Fig. 1, 2, 4/5, 6 and 7 draw Natural Earth coastlines via cartopy. On first use cartopy
downloads these once. If your run machine has no internet, pre-fetch them on a connected
machine and point cartopy at the cache:

```bash
export CARTOPY_DATA_DIR=/path/to/cartopy_cache
```

---

## Running

From the repo root, with the `dyetracer` environment active:

```bash
# Figure 7 needs its cache built first:
python scripts/precompute_proxymag.py

# then any figure (each writes PDFs into ./figures/):
python scripts/fig1_maps_amocmode.py
python scripts/fig2_dye_regions.py
python scripts/fig3_dye_propagation.py
python scripts/fig4_fig5_scenario_maps.py
python scripts/fig6_trajectories.py
python scripts/fig7_proxy_summary.py
```

Each script is self-contained: it adds the repo root to `sys.path`, changes to the repo
root, reads from `data/`, and writes its figure(s) to `figures/`. Outputs:

| Script | Main outputs (in `figures/`) |
|---|---|
| `fig1_maps_amocmode.py` | `Fig1_Rev_AMOCmodes.pdf`, `Fig1_Rev_VarianteDischarge.pdf` |
| `fig2_dye_regions.py` | `Fig2_Rev_DyeRegions.pdf`, `Fig2_Rev_Supp_Discharge.pdf` |
| `fig3_dye_propagation.py` | `Fig3_DyePropagation.pdf` |
| `fig4_fig5_scenario_maps.py` | `Fig5_Rev_ScenarioMap.pdf`, `Fig5_Rev_ScenarioMap_predict.pdf`, `Fig5_Rev_Ratio.pdf`, `Fig4_Rev_GIN_IRD.pdf`, `Fig4_Rev_GIN_IRD_predict_combined.pdf`, `Fig4_Rev_Scatter_RegionCompared_predicted.pdf` |
| `fig6_trajectories.py` | `Fig6_Rev_Trajectories.pdf` |
| `fig7_proxy_summary.py` | `Fig7_a_ProxySites_Map.pdf`, `Fig7_b_ProxySites_Bars_var.pdf`, `SFig7_b_ProxySites_Bars_var.pdf` + per-scenario `*_proxy_d18O_contributions.csv` |

---

## Data

Plotting data are archived on Zenodo alongside this code:

> **DOI:** {{ZENODO_DOI}}

The data payload is deliberately excluded from git (`.gitignore`: `data/*`); only
`data/README.md` — the data manifest and provenance — is tracked. Download the data
archive from the Zenodo record and unpack it at the repo root as described in **Setup**.

The bundle is ~45 MB (reduced from ~735 MB of raw HadCM3 output) and contains only the
fields required for plotting: model input grids, precomputed dye-field / sea-ice / MLD /
AMOC lookup caches, dye-mean fields, surface-only per-scenario mean/std fields, region and
proxy-site tables, and the six atmospheric trajectory files actually plotted. See
[`data/README.md`](data/README.md) for the full manifest, field descriptions, and
provenance notes.

---

## Bundled and external code

**Bundled analysis code** (`mymodules/`, `myconfig/`) is provided so the scripts run
standalone. Note that most analysis is *precomputed*: the shipped lookup/mean caches let
the pipeline functions short-circuit to loading cached results, so the raw-data reduction
steps are not re-executed here.

**External code — `mw_protocol`.** The original notebooks depended on the `mw_protocol`
package (meltwater-routing toolbox) published by Rome et al. Only two calls were used
across all figures, and both have been removed as external dependencies:

- `glac1d_toolbox.create_coordinate_edges` → reimplemented locally as
  `mymodules/grid_utils.create_coordinate_edges` (a short grid-edge helper).
- `plotting.create_discharge_ts` → its output is precomputed and shipped as
  `data/inputs/discharge_ts.csv`, so Fig. 1 runs without the package.

If you need the full `mw_protocol` toolbox, cite and obtain it from:

> Olnavy (2022). *Olnavy/rome2022_paleoceanography_oscillations: Reviews round 1 – v1.4*
> [Software]. Zenodo. https://doi.org/10.5281/zenodo.6788389

Base climate simulations follow Rome et al. (2022), HadCM3 (BRIDGE), with GLAC-1D
meltwater forcing.

---

## Notes on provenance / reproducibility

- **Cached pickles were generated at a slightly newer upstream revision than the original
  notebooks.** Two data keys differ from the notebook code (proxy-site names in
  `land_uptakemasks.pkl`, and discharge-scenario keys); the shipped scripts have been
  adapted to the shipped data, which is the ground truth for these figures.
- **`data/intermediates/dyestuff_modelpaper/merid/18.2k_mean_std.nc` is intentionally
  absent.** That experiment (xpujm) had only 245 of the required 500 years, so its
  years 450–499 statistics window is entirely missing; it is not used by any figure.
- The per-scenario mean/std fields are **surface-only** (top depth level) to keep the
  bundle light — this is all the figures require.

---

## Citation

If you use these scripts or data, please cite the manuscript (above) and the Zenodo
archive ({{ZENODO_DOI}}).
