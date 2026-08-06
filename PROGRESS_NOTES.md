# Progress notes — figure-script packaging

_Last updated: session pause (commit `a0d218f`)._

## Status: 6/6 scripts convert; 6/6 run end-to-end against bundled data*

\* Fig7 runs but is **slow** (>13 min) — flagged for a precompute refactor (see Open items).

| Script | Runs offline | Outputs (in `figures/`) |
|--------|--------------|--------------------------|
| `scripts/fig1_maps_amocmode.py`   | ✅ | Fig1_Rev_AMOCmodes.pdf, Fig1_Rev_VarianteDischarge.pdf |
| `scripts/fig2_dye_regions.py`     | ✅ | Fig2_Rev_DyeRegions.pdf, Fig2_Rev_Supp_Discharge.pdf |
| `scripts/fig3_dye_propagation.py` | ✅ | Fig3_DyePropagation.pdf |
| `scripts/fig4_fig5_scenario_maps.py` | ✅ | Fig5_Rev_ScenarioMap.pdf, Fig5_Rev_ScenarioMap_predict.pdf, Fig5_Rev_Ratio.pdf, Fig4_Rev_GIN_IRD.pdf, Fig4_Rev_GIN_IRD_predict_combined.pdf, Fig4_Rev_Scatter_RegionCompared_predicted.pdf |
| `scripts/fig6_trajectories.py`    | ✅ | Fig6_Rev_Trajectories.pdf |
| `scripts/fig7_proxy_summary.py`   | ⚠️ slow | Fig7_a_ProxySites_Map.pdf (+ bars, CSV — run to completion pending) |

## What was done
- **Notebook → `.py` conversion** for all six figures. Each script anchors to repo
  root, rewrites data paths to `data/`, drops interactive display/scratch cells,
  and writes to `figures/`.
- **`mymodules` completed/fixed:**
  - `grid_utils.py` (new) — local `create_coordinate_edges` replacing
    `mw_protocol.glac1d_toolbox` (external; cite Rome et al. 2022, Zenodo 6788389).
  - `d18O_computation.py` — added `area_weighted_mean`; fixed `build_d18O_results`
    to create the `{state}/` sub-directory before writing δ¹⁸O caches.
  - `dyefield_computation.py` — made `pylaeoclim_leeds` import graceful (only needed
    to *rebuild* the AMOC lookup, which the bundled cache avoids).
- **Data**: light bundle under `data/` (gitignored; Zenodo-only). δ¹⁸O `min/max`
  prediction caches are computed on first run of Fig4/5 from bundled `mean_std` data.
- **Env**: `dyetracer` (python 3.11); `dask` added this session (needed by
  `run_pipeline`'s `chunks="auto"`).

## Run instructions (local)
```bash
conda activate dyetracer   # or the env from environment.yml (still to write)
export CARTOPY_DATA_DIR="$PWD/.cartopy_cache"   # cache Natural Earth downloads
python scripts/fig3_dye_propagation.py          # etc.
```

## Open items (resume here)
1. **Fig7 precompute refactor** (user request): split the slow per-site/per-mode/
   per-dye land-proxy extraction into a one-off precompute step that writes a small
   table (e.g. `data/intermediates/proxy_contributions.pkl` / .csv), so the figure
   script only *plots*. Then re-time Fig7.
2. **Step 7 — dependency spec**: write `environment.yml` + `requirements.txt` from
   the `dyetracer` env (include `dask`).
3. **Step 8 — full local test pass** once Fig7 is fast.
4. **Step 9 — README**: dependencies, data provenance (Zenodo DOI placeholder
   `{{ZENODO_DOI}}`), `mw_protocol` external + citation, `CARTOPY_DATA_DIR` note.
5. **Step 10** — final commit/push (this note + env spec + README).
