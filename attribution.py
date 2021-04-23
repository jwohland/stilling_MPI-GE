import matplotlib.cm as cm
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pickle
import xesmf as xe
from scipy.stats import spearmanr, pearsonr
from utils import *

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy
    import cartopy.crs as ccrs


def quantile(da, q):
    tmp = da.values
    tmp = tmp[np.isfinite(tmp)]
    return np.quantile(tmp, q)


def count_quadrant(x, y):
    """
    count the number of occurences of
        x > 0 and
    :param x:
    :param y:
    :return:
    """


def conditional_prob(x, y, x_thr):
    joint_prob = y[(x > x_thr) & (y < 0)].size
    single_prob = y[x > x_thr].size
    if single_prob == 0:
        cond_prob = 0
    else:
        cond_prob = joint_prob / single_prob
    return cond_prob


# load data
slice_dic = {
    "historical": {"end": slice("1990", "2000"), "equiv": slice("1880", "1890")},
    "rcp85": {"end": slice("2090", "2100"), "equiv": slice("1970", "1980")},
    "rcp45": {"end": slice("2090", "2100"), "equiv": slice("1910", "1920")},
    "rcp26": {"end": slice("2090", "2100"), "equiv": slice("1890", "1900")},
}

ref = ann_mean(
    xr.open_dataset(
        "../data/historical/ensmean/sfcWind_Lmon_MPI-ESM_historical_ensmean_185001-200512.nc"
    )
)
ref = sel_time(ref, slice("1850", "1859")).mean(dim="time")

try:
    with open("../output/attribution_maps.pickle", "rb") as handle:
        ds_dict = pickle.load(handle)
except FileNotFoundError:
    ds_dict = {}
    ds_stylized = ann_mean(xr.open_dataset(glob.glob("../data/1pCO2/ensmean/*.nc")[0]))
    LUH_ref = open_LUH_period("../data/LUHa.v1/", 1850, 1860).mean(dim="time")
    for experiment in ["historical", "rcp26", "rcp45", "rcp85"]:
        ds_dict[experiment] = {}
        ds_tmp = (
            ann_mean(
                xr.open_dataset(glob.glob("../data/" + experiment + "/ensmean/*.nc")[0])
            )
            .sel({"time": slice_dic[experiment]["end"]})
            .mean(dim="time")
        )
        ds_dict[experiment]["Full Diff"] = (
            ds_tmp - ref
        )  # full (dynamical + land use) difference
        ds_dict[experiment]["Dyn. Diff"] = (
            ds_stylized.sel({"time": slice_dic[experiment]["equiv"]}).mean(dim="time")
            - ref
        )  # only CO2 forcing
        ds_dict[experiment]["Full - Dyn."] = (
            ds_dict[experiment]["Full Diff"] - ds_dict[experiment]["Dyn. Diff"]
        )
        path = "../data/LUHa.v1/"
        if experiment == "rcp26":
            path += "IMAGE_rcp26/"
        elif experiment == "rcp45":
            path += "MiniCAM_rcp45/"
        elif experiment == "rcp85":
            path += "MESSAGE_rcp85/"
        ds_dict[experiment]["LUH Diff"] = (
            open_LUH_period(
                path,
                int(slice_dic[experiment]["end"].start),
                int(slice_dic[experiment]["end"].stop),
            ).mean(dim="time")
            - LUH_ref
        )
        with open("../output/attribution_maps.pickle", "wb") as handle:
            pickle.dump(ds_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# prep Figure
f, ax = plt.subplots(
    ncols=4,
    nrows=4,
    figsize=(9, 6),
    sharex=True,
    sharey=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)
cbar1_ax = f.add_axes([0.05, 0.11, 0.6, 0.02])
cbar2_ax = f.add_axes([0.8, 0.11, 0.175, 0.02])

for row, experiment in enumerate(ds_dict.keys()):
    plt.text(
        -0.1,
        0.5,
        experiment,
        rotation=90,
        fontdict={"fontsize": 12},
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax[row, 0].transAxes,
    )
    for col, field in enumerate(ds_dict[experiment].keys()):
        if row + col == 0:
            ds_dict[experiment]["Full Diff"]["sfcWind"].plot(
                ax=ax[row, col],
                vmin=-1.5,
                vmax=1.5,
                add_colorbar=True,
                cbar_ax=cbar1_ax,
                extend="both",
                cmap=cm.get_cmap("RdBu_r"),
                cbar_kwargs={
                    "orientation": "horizontal",
                    "label": "Difference in decadal mean wind speed [m/s]",
                },
            )
        else:
            if field != "LUH Diff":
                ds_dict[experiment][field]["sfcWind"].plot(
                    ax=ax[row, col],
                    vmin=-1.5,
                    vmax=1.5,
                    add_colorbar=False,
                    extend="both",
                    cmap=cm.get_cmap("RdBu_r"),
                )
            else:
                ds_dict[experiment][field]["gothr+gsecd"].plot(
                    ax=ax[row, 3],
                    vmin=-1,
                    vmax=1,
                    add_colorbar=True,
                    cbar_ax=cbar2_ax,
                    extend="both",
                    cmap=cm.get_cmap("RdBu"),
                    cbar_kwargs={
                        "orientation": "horizontal",
                        "label": "Difference in primary \n plus secondary land",
                    },
                )

for tmp_ax in ax.flatten():
    tmp_ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
    tmp_ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
ax[0, 0].set_title(r"Full change $\Delta s$")
ax[0, 1].set_title(r"Dynamical change $\Delta_{dyn} s$")
ax[0, 2].set_title(r"Residual change $\Delta_{res} s$")
ax[0, 3].set_title(r"Land use change")

add_letters(ax, x=-0.05, y=0.98)
plt.subplots_adjust(0.03, 0.14, 0.99, 0.96, 0.1, 0.1)
plt.savefig("../plots/attribution_maps.jpeg", dpi=400)

# compare onshore pdfs
landmask = xr.open_dataarray("../data/runoff/landmask.nc")
colors = ["Orange", "Olive"]
f, ax = plt.subplots(nrows=4, sharex=True, figsize=(4, 8))
label_dict = {
    "Full - Dyn.": r"Residual change $\Delta_{res} s$",
    "Dyn. Diff": r"Dynamical change $\Delta_{dyn} s$",
}
for row, experiment in enumerate(ds_dict.keys()):
    total_change = 0
    for i, var in enumerate(["Full - Dyn.", "Dyn. Diff"]):
        tmp_da = ds_dict[experiment][var]["sfcWind"]
        tmp_da = tmp_da.where(np.isfinite(landmask)).values
        ax[row].hist(
            tmp_da[np.isfinite(tmp_da)],
            label=label_dict[var],
            density=True,
            bins=100,
            alpha=0.7,
            color=colors[i],
        )
        ax[row].axvline(tmp_da[np.isfinite(tmp_da)].mean(), ls="--", color=colors[i])
        total_change += tmp_da[np.isfinite(tmp_da)].mean()
    ax[row].axvline(
        total_change, ls="--", color="black", label=r"Full change $\Delta s$"
    )
    ax[row].set_ylabel(experiment + " PDF")
    ax[row].set_xlim(xmin=-0.8, xmax=0.8)
ax[0].legend(bbox_to_anchor=(1.03, 1.6))
ax[3].set_xlabel("Wind speed change [m/s]")
add_letters(ax)
plt.subplots_adjust(0.12, 0.06, 0.95, 0.89)
plt.savefig("../plots/contribution_histograms.jpeg", dpi=400)


"""
scatter plots
"""


lu_variable = "gothr+gsecd"
for wind_type in ["abs", "rel"]:
    f, ax = plt.subplots(nrows=4, sharex=True, sharey=True, figsize=(4, 8))
    for row, experiment in enumerate(ds_dict.keys()):
        ds_wind = ds_dict[experiment]["Full - Dyn."]["sfcWind"].copy()
        if wind_type == "rel":
            ds_wind /= ref["sfcWind"]

        ds_lu = ds_dict[experiment]["LUH Diff"][lu_variable].copy()
        # land use lons cover -180 ... 180, tranform to 0 ... 360
        ds_lu = ds_lu.assign_coords(lon=((ds_lu.lon + 360) % 360)).sortby("lon")
        # aggregate land use data to similar resolution (2x2 degrees compared to 1.875x1.875)
        ds_lu_agg = (
            ds_lu.coarsen(lat=4, boundary="trim")
            .mean()
            .coarsen(lon=4, boundary="trim")
            .mean()
        )

        regridder = xe.Regridder(ds_lu_agg, ds_wind, "bilinear", periodic=True)
        ds_lu_int = regridder(ds_lu_agg)
        y = ds_lu_int.where(np.abs(ds_lu_int) > 0.01)
        x = ds_wind.where(np.abs(ds_lu_int) > 0.01)

        # drop nans
        x, y = x.values, y.values
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        ax[row].scatter(
            x,
            y,
            s=5,
            alpha=0.3,
            label="R = "
            + str(np.round(pearsonr(x, y)[0], 2))
            + "\n"
            + r"$\rho$ = "
            + str(np.round(spearmanr(x, y)[0], 2)),
        )
        ax[row].set_ylabel("Land use change")
        ax[row].axhline(0, ls="--", color="Orange", alpha=0.7)
        ax[row].axvline(0, ls="--", color="Orange", alpha=0.7)

        # add occurence percentage of all four quadrants
        q1 = int(np.round(y[(x < 0) & (y > 0)].size / y.size * 100))
        q2 = int(np.round(y[(x > 0) & (y > 0)].size / y.size * 100))
        q3 = int(np.round(y[(x > 0) & (y < 0)].size / y.size * 100))
        q4 = int(np.round(y[(x < 0) & (y < 0)].size / y.size * 100))

        ax[row].text(
            -1, 0.2, str(q1) + "%", color="Darkorange", ha="center", va="center"
        )
        ax[row].text(
            2.5, 0.2, str(q2) + "%", color="Darkorange", ha="center", va="center"
        )
        ax[row].text(
            2.5, -0.2, str(q3) + "%", color="Darkorange", ha="center", va="center"
        )
        ax[row].text(
            -1, -0.2, str(q4) + "%", color="Darkorange", ha="center", va="center"
        )
        ax[row].legend(loc="upper right", markerscale=0)

        # add conditional probability
        ax2 = ax[row].twinx()
        ys = [
            conditional_prob(x, y, x_thr=x_thr) * 100 for x_thr in np.arange(0, 2, 0.1)
        ]
        ax2.plot(np.arange(0, 2, 0.1), ys, color="Purple")
        ax2.set_ylabel("cond. prob. luc <0 [%]", color="Purple")
        ax2.set_ylim(ymin=95, ymax=107)
        ax2.set_yticks([96, 98, 100])
        [t.set_color("Purple") for t in ax2.yaxis.get_ticklines()]
        [t.set_color("Purple") for t in ax2.yaxis.get_ticklabels()]
    # ax[3].set_xlim(xmin=-0.5, xmax=2)
    if wind_type == "abs":
        ax[3].set_xlabel("Residual wind speed change [m/s]")
    else:
        ax[3].set_xlabel("Relative residual wind speed change [1]")
    add_letters(ax)
    plt.tight_layout()
    plt.savefig(
        "../plots/scatter/scatter_" + lu_variable + "_" + wind_type + ".jpeg", dpi=300,
    )
    plt.close("all")
