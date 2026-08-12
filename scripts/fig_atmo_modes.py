#!/usr/bin/env python3
"""
Atmosphere across AMOC modes (17.8 ka): Surface Air Temperature, 850 hPa wind,
and 2-6 day band-pass storm-track intensity (annual + DJFM), for the three
AMOC modes cold (xpraj) / zonal (xprak) / merid (xpral).

Reads the light precomputed bundle written by scripts/precompute_atmo.py:
    data/intermediates/atmo/sat_mean.nc
    data/intermediates/atmo/wind850_mean.nc
    data/intermediates/atmo/stormtrack.nc
    data/intermediates/atmo/seaice_monthly.nc
and the land-sea mask from data/inputs/temev.qrparm.omask.nc .

Outputs (into figures/):
    SFig_Atmo_SAT.pdf
    SFig_Atmo_Wind850.pdf
    SFig_Atmo_StormTrack_annual.pdf
    SFig_Atmo_StormTrack_DJFM.pdf

NOTE: the storm-track panels apply a scaling factor s = 1/10 carried over from
the original analysis notebook (flagged there as a units issue still to be
confirmed). Kept here to reproduce the submitted figure exactly.
"""
import os
import sys
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmcrameri.cm as cmc

# --- repo-anchored paths -----------------------------------------------------
_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "mymodules"))
os.chdir(_REPO)
(_REPO / "figures").mkdir(exist_ok=True)

from grid_utils import create_coordinate_edges  # noqa: E402

ATMO = _REPO / "data" / "intermediates" / "atmo"
INPUTS = _REPO / "data" / "inputs"

MODES = ["cold", "zonal", "merid"]
MODE_LABEL = {"cold": "cold", "zonal": "zonal", "merid": "merid"}
STORM_SCALE = 1.0 / 10.0   # notebook scaling factor (see module docstring)

# --- load light bundle -------------------------------------------------------
sat = xr.open_dataset(ATMO / "sat_mean.nc")["temp_mm_srf"]        # (mode, lat, lon) [K]
wind = xr.open_dataset(ATMO / "wind850_mean.nc")                  # u, v (mode, lat, lon)
storm = xr.open_dataset(ATMO / "stormtrack.nc")                  # annual, djfm (mode, lat, lon)
ice = xr.open_dataset(ATMO / "seaice_monthly.nc")["iceconc"]     # (mode, month, lat, lon)

lsm = xr.open_dataset(INPUTS / "temev.qrparm.omask.nc").lsm

lon_atm = sat.longitude.values
lat_atm = sat.latitude.values
lon_atm_b = create_coordinate_edges(lon_atm)
lat_atm_b = lat_atm


def _sorted_180(da):
    """Return (lon_sorted, da_sorted) with longitude wrapped to [-180,180]."""
    lon = da.longitude.values
    lon180 = np.where(lon > 180, lon - 360, lon)
    order = np.argsort(lon180)
    return lon180[order], da.isel(longitude=order)


# =============================================================================
# 1. Surface Air Temperature (annual mean, deg C) + zonal-minus-cold difference
# =============================================================================
def plot_sat():
    proj = ccrs.Robinson()
    fig, axes = plt.subplots(1, 4, subplot_kw={"projection": proj}, figsize=(20, 5))
    normsat = Normalize(vmin=-20, vmax=40)
    normsat_diff = Normalize(vmin=-5, vmax=10)

    def _panel(ax, field, title, cmap, norm):
        m = ax.pcolormesh(lon_atm_b, lat_atm_b, field,
                          transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, alpha=0.8)
        ax.contour(lsm.longitude, lsm.latitude, lsm,
                   cmap=ListedColormap(["black"]), transform=ccrs.PlateCarree(),
                   zorder=1, alpha=1)
        ax.set_title(title, loc="left", fontsize=10)
        return m

    sat_c = {m: sat.sel(mode=m) - 273.15 for m in MODES}
    m0 = _panel(axes[0], sat_c["cold"], "SAT annual mean cold 17.8ka", cmc.vik, normsat)
    _panel(axes[1], sat_c["zonal"], "SAT annual mean zonal 17.8ka", cmc.vik, normsat)
    _panel(axes[2], sat_c["merid"], "SAT annual mean merid 17.8ka", cmc.vik, normsat)
    md = _panel(axes[3], sat_c["zonal"] - sat_c["cold"],
                "SAT annual zonal-cold diff 17.8ka", cmc.hawaii_r, normsat_diff)

    cb = fig.colorbar(m0, ax=axes, orientation="horizontal", shrink=0.4, extend="both")
    cb.set_label("Annual mean 2m Temperature [\u00b0C]", size="small")
    cb.ax.tick_params(labelsize="small")
    cb2 = fig.colorbar(md, ax=axes, orientation="horizontal", shrink=0.4, extend="both")
    cb2.set_label("Annual 2m Temperature difference [\u00b0C]", size="small")
    cb2.ax.tick_params(labelsize="small")
    fig.savefig("figures/SFig_Atmo_SAT.pdf", bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# 2. 850 hPa wind speed + vectors, and cold-minus-zonal difference
# =============================================================================
def plot_wind850():
    ws = {m: np.sqrt(wind["u"].sel(mode=m) ** 2 + wind["v"].sel(mode=m) ** 2) for m in MODES}
    lon_sorted, _ = _sorted_180(ws["cold"])
    ws_s = {m: _sorted_180(ws[m])[1] for m in MODES}
    u_s = {m: _sorted_180(wind["u"].sel(mode=m))[1] for m in MODES}
    v_s = {m: _sorted_180(wind["v"].sel(mode=m))[1] for m in MODES}
    lat = ws["cold"].latitude.values

    fig, ax = plt.subplots(2, 2, subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(16, 10))
    ax = ax.flatten()
    for a in ax:
        a.set_extent([-100, 20, 0, 88], crs=ccrs.PlateCarree())
        a.add_feature(cfeature.COASTLINE)
        a.add_feature(cfeature.LAND, facecolor="lightgray", alpha=0.5)

    for i, m in enumerate(MODES):
        levels = np.linspace(0, float(ws_s[m].max()), 20)
        cf = ax[i].contourf(lon_sorted, lat, ws_s[m], levels=levels, cmap="coolwarm", extend="both")
        ax[i].quiver(lon_sorted, lat, u_s[m].values, v_s[m].values, scale=200, color="black")
        cbar = fig.colorbar(cf, ax=ax[i], orientation="vertical", pad=0.02, aspect=30)
        cbar.set_label("Wind Speed (m/s)")
        ax[i].set_title(f"{MODE_LABEL[m]} 850 hPa Wind Speed")

    cfd = ax[3].contourf(lon_sorted, lat, ws_s["cold"] - ws_s["zonal"], cmap="coolwarm", extend="both")
    cbar = fig.colorbar(cfd, ax=ax[3], orientation="vertical", pad=0.02, aspect=30)
    cbar.set_label("Diff. Wind Speed (m/s)")
    ax[3].set_title("cold-zonal 850 hPa Wind Speed")

    fig.tight_layout()
    fig.savefig("figures/SFig_Atmo_Wind850.pdf", bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# 3. Storm-track intensity (band-pass MSL variance) with 50% sea-ice contours
# =============================================================================
def _plot_stormtrack(var_name, ice_months, ice_legend, outfile, title_tag):
    from matplotlib.lines import Line2D

    st = storm[var_name]  # (mode, lat, lon)
    lon_sorted, _ = _sorted_180(st.sel(mode="cold"))
    st_s = {m: _sorted_180(st.sel(mode=m))[1] for m in MODES}
    lat = st.sel(mode="cold").latitude.values
    lon = st.sel(mode="cold").longitude.values

    vmax = float(st_s["cold"].max()) * STORM_SCALE
    levels = np.linspace(0, vmax, 20)

    fig, ax = plt.subplots(2, 2, subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(16, 10))
    ax = ax.flatten()
    for a in ax:
        a.set_extent([-100, 20, 0, 88], crs=ccrs.PlateCarree())
        a.add_feature(cfeature.COASTLINE)
        a.add_feature(cfeature.LAND, facecolor="lightgray", alpha=0.5)

    handles = [Line2D([0], [0], color="green", linestyle="-", linewidth=2, label=ice_legend[0])]
    if len(ice_legend) > 1:
        handles.append(Line2D([0], [0], color="red", linestyle="--", linewidth=2, label=ice_legend[1]))

    for i, m in enumerate(MODES):
        cf = ax[i].contourf(lon_sorted, lat, st_s[m] * STORM_SCALE, levels=levels,
                            cmap="coolwarm", extend="both", transform=ccrs.PlateCarree())
        ice_m = ice.sel(mode=m)
        if isinstance(ice_months[0], list):  # DJFM: mean over months
            ice_field = ice_m.sel(month=ice_m.month.isin(ice_months[0])).mean(dim="month")
            ax[i].contour(lon, lat, ice_field, transform=ccrs.PlateCarree(),
                          levels=[0.5], colors="green", linewidths=1.2, linestyles="-")
        else:                                # annual: March (max) + September (min)
            ax[i].contour(lon, lat, ice_m.sel(month=ice_months[0]), transform=ccrs.PlateCarree(),
                          levels=[0.5], colors="green", linewidths=1.2, linestyles="-")
            ax[i].contour(lon, lat, ice_m.sel(month=ice_months[1]), transform=ccrs.PlateCarree(),
                          levels=[0.5], colors="xkcd:red", linewidths=1.2, linestyles="--")
        cbar = fig.colorbar(cf, ax=ax[i], orientation="vertical", pad=0.02, aspect=30)
        cbar.set_label("Storm Track Intensity (hPa\u00b2)")
        ax[i].legend(handles=handles, loc="lower center", title="50% sea ice extent")
        ax[i].set_title(f"{MODE_LABEL[m]}: {title_tag} Storm Track intensity")

    cfd = ax[3].contourf(lon_sorted, lat, (st_s["cold"] - st_s["zonal"]) * STORM_SCALE,
                        cmap="coolwarm", extend="both", transform=ccrs.PlateCarree())
    cbar = fig.colorbar(cfd, ax=ax[3], orientation="vertical", pad=0.02, aspect=30)
    cbar.set_label("Diff in Storm Track Intensity (hPa\u00b2)")
    ax[3].set_title(f"{title_tag} storm-track intensity difference (cold-zonal)")

    fig.tight_layout()
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def plot_stormtrack_annual():
    _plot_stormtrack("annual", [3, 9], ["March (max)", "September (min)"],
                     "figures/SFig_Atmo_StormTrack_annual.pdf", "Mean annual")


def plot_stormtrack_djfm():
    _plot_stormtrack("djfm", [[12, 1, 2, 3]], ["Avg Winter"],
                     "figures/SFig_Atmo_StormTrack_DJFM.pdf", "DJFM")


if __name__ == "__main__":
    plot_sat()
    plot_wind850()
    plot_stormtrack_annual()
    plot_stormtrack_djfm()
    print("done: SFig_Atmo_SAT.pdf, SFig_Atmo_Wind850.pdf, "
          "SFig_Atmo_StormTrack_annual.pdf, SFig_Atmo_StormTrack_DJFM.pdf")
