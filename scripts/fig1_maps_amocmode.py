"""Figure 1 - AMOC-mode maps (sea ice, MLD, overturning) and meltwater discharge.

Reproduces Figure 1 of:
  "Tracing Meltwater from Northern Ice Sheets to Palaeoclimate Archives During
   the Early Last Deglaciation: a conservative tracer approach"

Run from anywhere; the script anchors itself to the repository root, adds it to
sys.path (so ``mymodules`` / ``myconfig`` import), and reads the light plotting
data under ``data/`` (see data/README.md for how to obtain it from Zenodo).
Figures are written to ``figures/``.

Map figures use Cartopy, which downloads Natural Earth coastline data on first
run. Set the environment variable ``CARTOPY_DATA_DIR`` to a writable directory
to cache it (see README).
"""
import os, sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
os.chdir(_REPO)
os.makedirs("figures", exist_ok=True)

# --- cell 0 ---
import os

# --- cell 1 ---
import xarray as xr
import cartopy.feature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import cmcrameri.cm as cmc
from matplotlib.colors import BoundaryNorm, TwoSlopeNorm, LinearSegmentedColormap, ListedColormap, Normalize
import matplotlib as mpl

import regionmask
import geopandas as gpd
import pandas as pd
#Yvans
from mymodules import grid_utils as tb

#Mine
import mymodules.myfunc as mf 


colores = {'0':'tab:blue', '1':'tab:orange', '2':'tab:green', '3':'tab:red', '4':'tab:purple', '5':'tab:brown', '6':'tab:pink',
         '7':'tab:olive', '8':'tab:cyan'}

# --- cell 2 ---
from mymodules.dyefield_computation import build_seaice_lookup
from myconfig.EXPERIMENTS import EXPERIMENTS

seaice_lookup = build_seaice_lookup(
    EXPERIMENTS,
    outdir="data/intermediates",
    overwrite=False,
    exclude=('xqeic','xqeie','xpujm'), #simulations to be excluded here from the total experiments
)

# --- cell 3 ---
import xarray as xr


def compute_seaice_mode_stats(seaice_lookup):
    """
    Compute mean and standard deviation of March and September
    sea-ice concentration for each circulation mode.

    Returns
    -------
    seaice_stats : dict

    seaice_stats[mode]["march"]["mean"]
    seaice_stats[mode]["march"]["std"]
    seaice_stats[mode]["sept"]["mean"]
    seaice_stats[mode]["sept"]["std"]
    """

    seaice_stats = {}

    for mode, exp_dict in seaice_lookup.items():

        seaice_stats[mode] = {}

        for season in ["march", "sept"]:

            fields = [
                exp_dict[exp][season]
                for exp in exp_dict
            ]

            stack = xr.concat(fields, dim="experiment")

            seaice_stats[mode][season] = {
                "mean": stack.mean("experiment"),
                "std": stack.std("experiment"),
                "n": stack.sizes["experiment"],
            }

    return seaice_stats



seaice_stats = compute_seaice_mode_stats(seaice_lookup)

# --- cell 4 ---
#Load Land Sea Mask
with mf.Timer('Land sea mask loader'):
    data_folder = "data/inputs"
    ds_lsm = xr.open_dataset(f"{data_folder}/temev.qrparm.omask.nc")
    lsm = ds_lsm.lsm


# --- cell 5 ---
import mymodules.dyefield_computation as dc

dir(dc)

# --- cell 6 ---
from mymodules.dyefield_computation import build_mld_lookup
from myconfig.EXPERIMENTS import EXPERIMENTS

mld_lookup = build_mld_lookup(
    EXPERIMENTS,
    outdir="data/intermediates",
    overwrite=False,
    exclude=('xqeic','xqeie','xpujm'), #simulations to be excluded here from the total experiments
)

# --- cell 7 ---
# mld stats

import xarray as xr


def compute_mld_mode_stats(mld_lookup):
    """
    Compute mean and standard deviation
    of winter MLD for each circulation mode.

    Returns
    -------
    mld_stats[mode]["mean"]
    mld_stats[mode]["std"]
    """

    mld_stats = {}

    for mode, exp_dict in mld_lookup.items():

        fields = [
            exp_dict[exp]
            for exp in exp_dict
        ]

        stack = xr.concat(
            fields,
            dim="experiment"
        )

        mld_stats[mode] = {
            "mean": stack.mean("experiment"),
            "std": stack.std("experiment"),
            "n": stack.sizes["experiment"],
        }

    return mld_stats


mld_stats = compute_mld_mode_stats(mld_lookup)



# --- cell 8 ---
import mymodules.dyefield_computation as dc

#dir(dc)

from mymodules.dyefield_computation import build_amoc_lookup
from myconfig.EXPERIMENTS import EXPERIMENTS

amoc_lookup = build_amoc_lookup(
    EXPERIMENTS,
    outdir="data/intermediates",
    overwrite=False,
    exclude=('xqeic','xqeie','xpujm'), #simulations to be excluded here from the total experiments
)

# --- cell 9 ---
# AMOC stats

def compute_amoc_mode_stats(amoc_lookup):

    """
    Compute ensemble mean and std AMOC time series.

    Returns
    -------
    amoc_stats[mode]["raw"]["mean"]
    amoc_stats[mode]["raw"]["std"]

    amoc_stats[mode]["filtered"]["mean"]
    amoc_stats[mode]["filtered"]["std"]
    """

    amoc_stats = {}


    for mode, exp_dict in amoc_lookup.items():

        amoc_stats[mode] = {}


        for key in ["raw", "filtered"]:

            fields = [
                exp_dict[exp][key]
                for exp in exp_dict
            ]


            stack = xr.concat(
                fields,
                dim="experiment"
            )


            amoc_stats[mode][key] = {
                "mean": stack.mean("experiment"),
                "std": stack.std("experiment"),
                "n": stack.sizes["experiment"],
            }

            
        ensemble_means = [
        exp_dict[exp]["raw"].mean().item()
        for exp in exp_dict
            ]
    
        amoc_stats[mode]["ensemble summary"] = {
        "mean": np.mean(ensemble_means),
        "std": np.std(ensemble_means),
            }


    #mean_amoc = amoc_stats[mode]["summary"]["mean"]
    #std_amoc  = amoc_stats[mode]["summary"]["std"]

    return amoc_stats


amoc_stats = compute_amoc_mode_stats(amoc_lookup)



# --- cell 11 ---
## Plotting

# --- cell 13 ---
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D


# --------------------------------------------------
# SETTINGS
# --------------------------------------------------

modes = ["merid", "zonal", "cold"]

scenarios =  {'17.8k','18.2k','19.4k','20.7k'}

scenario_styles = {
    "17.8k": "-",
    "18.2k": "--",
    "19.4k": "-.",
    "20.7k": ":",
}

mode_colors = {
    "merid": "tab:blue",
    "zonal": "tab:orange",
    "cold": "tab:green",
}


projection = ccrs.NorthPolarStereo(central_longitude=-30)

# choose threshold for stippling
mld_std_threshold = 50.0   # metres

# grid for stippling
stipple_stride = 2

cmap=cmc.lapaz_r

# --------------------------------------------------
# FIGURE LAYOUT
# --------------------------------------------------

fig = plt.figure(figsize=(12,6))

gs = gridspec.GridSpec(
    nrows=2,
    ncols=3,
    height_ratios=[1,1],
    wspace=0.1,
    hspace=0.2
)


# map axes
ax_map = [
    fig.add_subplot(gs[0,i], projection=projection)
    for i in range(3)
]

# time series axes
ax_ts = [
    fig.add_subplot(gs[1,i])
    for i in range(3)
]


# --------------------------------------------------
# TOP ROW: MEAN MLD MAP + STD STIPPLING
# --------------------------------------------------

norm = Normalize(vmin=0,vmax=640)

im = None

for i,mode in enumerate(modes):

    ax = ax_map[i]

    ax.set_extent(
        [-85,5,35,90],
        crs=ccrs.PlateCarree()
    )

    ax.coastlines(zorder=10)

    # land
    #ax.add_feature(
    #    cfeature.LAND,
    #    facecolor="lightgrey",
    #    zorder=1
    #)


    # -------------------------
    # Mean MLD
    # -------------------------

    pcm = ax.pcolormesh(
        mld_stats[mode]["mean"].longitude,
        mld_stats[mode]["mean"].latitude,
        xr.where(
            mld_stats[mode]["mean"] == 0,
            np.nan,
            mld_stats[mode]["mean"],
        ),
        cmap=cmap,
        norm=norm,
        shading="nearest",
        transform=ccrs.PlateCarree(),
        zorder=1,
    )

    # -------------------------
    # MLD stippling
    # -------------------------
    lon = mld_stats[mode]["mean"].longitude,
    lat = mld_stats[mode]["mean"].latitude

    std = mld_stats[mode]["std"]
    
    mask = std > mld_std_threshold

    xx, yy = np.meshgrid(
        lon,
        lat,
    )

    ax.scatter(
        xx[::stipple_stride, ::stipple_stride][mask.values[::stipple_stride, ::stipple_stride]],
        yy[::stipple_stride, ::stipple_stride][mask.values[::stipple_stride, ::stipple_stride]],
        s=2,
        color="k",
        marker=".",
        transform=ccrs.PlateCarree(),
        zorder=4,
    )

    # -------------------------
    # March sea ice
    # -------------------------

    march_mean = seaice_stats[mode]["march"]["mean"]
    march_std = seaice_stats[mode]["march"]["std"]

    lon_si = seaice_stats[mode]["sept"]["mean"].longitude
    lat_si = seaice_stats[mode]["sept"]["mean"].latitude

    ax.contour(
        lon_si,
        lat_si,
        march_mean,
        levels=[0.5],
        colors="tab:red",
        linestyles="--",
        linewidths=2,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )

    ax.contour(
        lon_si,
        lat_si,
        march_mean + march_std,
        levels=[0.5],
        colors="tab:red",
        linewidths=1,
        linestyles="--",
        alpha=0.6,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )

    ax.contour(
        lon_si,
        lat_si,
        march_mean - march_std,
        levels=[0.5],
        colors="tab:red",
        linewidths=1,
        linestyles="--",
        alpha=0.6,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )

    # -------------------------
    # September sea ice
    # -------------------------

    sept_mean = seaice_stats[mode]["sept"]["mean"]
    sept_std = seaice_stats[mode]["sept"]["std"]

    


    ax.contour(
        lon_si,
        lat_si,
        sept_mean,
        levels=[0.5],
        colors="tab:red",
        linewidths=2,
        linestyles="-",
        transform=ccrs.PlateCarree(),
        zorder=5,
    )

    ax.contour(
        lon_si,
        lat_si,
        sept_mean + sept_std,
        levels=[0.5],
        colors="tab:red",
        linewidths=1,
        linestyles="--",
        alpha=0.6,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )

    ax.contour(
        lon_si,
        lat_si,
        sept_mean - sept_std,
        levels=[0.5],
        colors="tab:red",
        linewidths=1,
        linestyles="--",
        alpha=0.6,
        transform=ccrs.PlateCarree(),
        zorder=5,
    )

    ax.set_title(mode.capitalize())


# colorbar
cax = fig.add_axes([0.92,0.55,0.015,0.3])

cb = fig.colorbar(
    pcm,
    cax=cax
)

cb.set_label(
    "Mean winter MLD (m)"
)



# --------------------------------------------------
# BOTTOM ROW: AMOC TIME SERIES
# --------------------------------------------------

for i,mode in enumerate(modes):

    ax = ax_ts[i]

    # --------------------------------------------------
    # LOOP OVER MODES
    # --------------------------------------------------

#for ax, mode in zip(axes, modes):

    color = mode_colors[mode]

    # -----------------------------------------
    # Ensemble mean/std over ALL experiments
    # -----------------------------------------

    #mean_amoc = float(
    #    amoc_stats[mode]["raw"]["mean"].mean()
    #)

    #std_amoc = float(
    #    amoc_stats[mode]["raw"]["mean"].std()
    #)

    mean_amoc = amoc_stats[mode]["ensemble summary"]['mean']
    std_amoc = amoc_stats[mode]["ensemble summary"]['std']


    # mean line
    ax.axhline(
        mean_amoc,
        color="k",
        lw=2,
        zorder=1,
    )

    # ±1 std
    ax.axhspan(
        mean_amoc - std_amoc,
        mean_amoc + std_amoc,
        color="0.8",
        alpha=0.4,
        zorder=0,
    )

    # -----------------------------------------
    # Plot each scenario
    # -----------------------------------------

    for scenario in scenarios:

        if scenario not in EXPERIMENTS[mode]:
            continue

        exp = EXPERIMENTS[mode][scenario]["exp"]

        if exp not in amoc_lookup[mode]:
            continue

        raw = amoc_lookup[mode][exp]["raw"]
        filt = amoc_lookup[mode][exp]["filtered"]

        x = np.arange(raw.size)

        # raw (faint)
        ax.plot(
            x,
            raw,
            color=color,
            alpha=0.25,
            lw=1,
        )

        # filtered
        ax.plot(
            x,
            filt,
            color=color,
            linestyle=scenario_styles[scenario],
            lw=2,
            label=scenario,
            alpha=0.8
        )

    # -----------------------------------------
    # Formatting
    # -----------------------------------------

    #ax.set_title(mode.capitalize())

    ax.set_ylim(5, 20)

    ax.grid(
        alpha=0.3,
        linestyle=":"
    )

    ax.set_xlabel("Simulation year")


ax_ts[0].set_ylabel("AMOC (Sv)")

for ax in ax_ts[1:]:
    ax.set_yticklabels([])



# --------------------------------------------------
# LEGENDS
# --------------------------------------------------

# --------------------------------------------------
# LEGENDS
# --------------------------------------------------

# Sea ice edge legend (top panel, left)
seaice_handles = [
    Line2D(
        [],
        [],
        color="tab:red",
        lw=2,
        ls="--",
        label="March"
    ),
    Line2D(
        [],
        [],
        color="tab:red",
        lw=2,
        ls="-",
        label="Sept"
    ),
    #Line2D(
    #    [],
    #    [],
    #    color="tab:red",
    #    lw=1,
    #    ls="--",
    #    alpha=0.6,
    #    label="±1 std"
    #),
]

fig.legend(
    handles=seaice_handles,
    loc="upper left",
    bbox_to_anchor=(0.125, 0.88),
    frameon=True,
    title="Sea ice edge",
)


# Scenario legend (bottom row, right)
scenario_handles = [
    Line2D(
        [],
        [],
        color="black",
        linestyle=style,
        lw=2,
        label=scen
    )
    for scen, style in scenario_styles.items()
]

fig.legend(
    handles=scenario_handles,
    loc="upper left",
    bbox_to_anchor=(0.815, 0.46),
    frameon=True,
    title="Scenario",
)


#plt.tight_layout(rect=[0,0,1,1])

#plt.tight_layout(rect=[0,0.05,0.9,1])

#plt.show()

# Save the Figure
plt.savefig(
    "figures/Fig1_Rev_AMOCmodes.pdf",
    format="pdf",
    bbox_inches="tight",
    dpi=300,
)



# --- cell 14 ---
# Plot A: the discharge timeseries (precomputed via mw_protocol.plotting.create_discharge_ts,
# stored as data/inputs/discharge_ts.csv so the plot runs without the external package).
df = pd.read_csv("data/inputs/discharge_ts.csv")
t = df["t"].values


# --- cell 15 ---
# (dataframe already loaded from discharge_ts.csv above)


# --- cell 16 ---
# Plot
fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(df["t"]/1000, df["tot"], color="deepskyblue", linewidth=2)

# Interpolation times (the ones you want to mark)
highlight_times = np.array([-17800, -18200, -19400, -20700])

# Interpolate values at those times
interpolated_values = np.interp(highlight_times, df["t"], df["tot"])

# 🔹 Vertical snapshot lines with different styles
line_styles = ["-", "--", "-.", ":"]

for t, ls in zip(highlight_times, line_styles):
    ax.axvline(x=t/1000, color="black", linestyle=ls, linewidth=2)

# Reversed x-axis ⏪
ax.set_xlim(-17, -22)
ax.set_ylim(0.05, 0.3)

# Remove top and right borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_axisbelow(True)

# Labels — explicit ordered list so each label matches its line style
# (highlight_times [-17800,-18200,-19400,-20700] -> line_styles ["-","--","-.",":"]).
# NB: `scenarios` is a set and iterates in arbitrary order, which scrambled the legend.
labels = ["17.8k", "18.2k", "19.4k", "20.7k"]
for t, ls, lab in zip(highlight_times, line_styles, labels):
    ax.axvline(x=t/1000, color="black", linestyle=ls, linewidth=2, label=lab)

ax.legend(frameon=False)

ax.set_xlabel("Time (kyr)")
ax.set_ylabel("Total melt discharge (Sv)")

# Save the Figure
plt.savefig(
    "figures/Fig1_Rev_VarianteDischarge.pdf",
    format="pdf",
    bbox_inches="tight",
    dpi=300,
)

