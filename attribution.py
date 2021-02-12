# plot global map of wind speed trends over 1900 - 2010

import xarray as xr
import glob
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import warnings
import numpy as np
import pickle

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy
    import cartopy.crs as ccrs


def ann_mean(ds):
    return ds.resample({"time": "1Y"}).mean()


def sel_time(ds, tslice):
    return ds.sel({"time": tslice})


def open_datasets_ensmean(filelist):
    ds = [ann_mean(sel_time(xr.open_dataset(x, use_cftime=True))) for x in filelist]
    ds = xr.concat(
        ds,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in filelist]),
    )
    ds = ds.mean(dim="ensemble_member")
    return ds


def plot_field(data, ax=None, title=None, **kwargs):
    """
    plots maps
    :param data: data xarray that shall be plotted
    :param ax: axes
    :param title: figure title
    :return:
    """
    ax = ax or plt.axes(projection=ccrs.PlateCarree())
    data.plot(ax=ax, **kwargs)
    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
    ax.set_title(title)
    return ax


def open_LUH(filename, year=None, name=None):
    da = (
        xr.open_rasterio(filename)
        .rename({"x": "lon", "y": "lat"})
        .drop("band")
        .squeeze()
    )
    if year:
        da = da.assign_coords({"time": pd.Timestamp(str(year) + "-01-01")})
    if name:
        da.name = name
    return da


def open_LUH_period(path, ts, te):
    """
    open time varying primary and secondary vegeation maps and combine them in a single xarray Dataset
    :param path: path to data
    :param ts: start year of respective experiment, e.g. 1850 for historical
    :param te: end year for respective experiment, e.g. 2000 for historical
    :return:
    """
    da_gothr = xr.concat(
        [
            open_LUH(path + "updated_states/gothr." + str(year) + ".txt", year, "gothr")
            for year in range(ts, te)
        ],
        dim="time",
    )
    da_gsecd = xr.concat(
        [
            open_LUH(path + "updated_states/gsecd." + str(year) + ".txt", year, "gsecd")
            for year in range(ts, te)
        ],
        dim="time",
    )
    ds = xr.merge([da_gothr, da_gsecd])
    ds["gothr+gsecd"] = ds["gothr"] + ds["gsecd"]
    return ds


# load data
slice_dic = {
    "historical": {"end": slice("1990", "2000"), "equiv": slice("1880", "1890")},
    "rcp85": {"end": slice("2090", "2100"), "equiv": slice("1970", "1980")},
    "rcp45": {"end": slice("2090", "2100"), "equiv": slice("1910", "1920")},
    "rcp26": {"end": slice("2090", "2100"), "equiv": slice("1890", "1900")},
}

ref = ann_mean(xr.open_dataset(glob.glob("../data/historical/ensmean/*.nc")[0]))
ref = sel_time(ref, slice("1850", "1859")).mean(dim="time")
ds_stylized = ann_mean(xr.open_dataset(glob.glob("../data/1pCO2/ensmean/*.nc")[0]))
LUH_ref = open_LUH_period("../data/LUHa.v1/", 1850, 1860).mean(dim="time")

try:
    with open("../output/attribution_maps.pickle", "rb") as handle:
        ds_dict = pickle.load(handle)
except FileNotFoundError:
    ds_dict = {}
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
    figsize=(12, 8),
    sharex=True,
    sharey=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)
cbar1_ax = f.add_axes([0.05, 0.08, 0.6, 0.02])
cbar2_ax = f.add_axes([0.8, 0.08, 0.175, 0.02])

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
ax[0, 2].set_title(r"Full minus dynamical $\Delta s - \Delta_{dyn} s$")
ax[0, 3].set_title(r"Land use change")


plt.subplots_adjust(0.04, 0.1, 0.99, 0.96)
plt.savefig("../plots/attribution_maps.jpeg", dpi=300)

# compare onshore pdfs
landmask = xr.open_dataarray("../data/runoff/landmask.nc")
colors = ["Orange", "Olive"]
f, ax = plt.subplots(nrows=4, sharex=True, figsize=(4, 10))
for row, experiment in enumerate(ds_dict.keys()):
    for i, var in enumerate(["Full - Dyn.", "Dyn. Diff"]):
        tmp_da = ds_dict[experiment][var]["sfcWind"]
        tmp_da = tmp_da.where(np.isfinite(landmask)).values
        ax[row].hist(
            tmp_da[np.isfinite(tmp_da)],
            label=var,
            density=True,
            bins=100,
            alpha=0.7,
            color=colors[i],
        )
        ax[row].axvline(tmp_da[np.isfinite(tmp_da)].mean(), ls="--", color=colors[i])
    ax[row].set_ylabel(experiment + " PDF")
    ax[row].set_xlim(xmin=-1, xmax=1)
ax[0].legend()
ax[3].set_xlabel("Wind speed change [m/s]")
plt.tight_layout()
plt.savefig("../plots/contribution_histograms.jpeg", dpi=300)


"""
The following is under construction!!!
"""
from scipy.stats import spearmanr, pearsonr


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


def conditional_prob(x,y,x_thr):
    joint_prob = y[(x>x_thr) & (y < 0)].size
    single_prob = y[x>x_thr].size
    if single_prob == 0:
        cond_prob = 0
    else:
        cond_prob = joint_prob/single_prob
    return cond_prob

lu_variable = "gothr+gsecd"
# todo maybe do masking differently and split by geographic area?
for wind_type in ["abs", "rel"]:
    for analysis_type in [
        "all",
        "hotspots",
    ]:  # hotspots only considers locations with large wind speed changes
        f, ax = plt.subplots(nrows=4, sharex=True, sharey=True, figsize=(4, 10))
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
            ds_lu_int = ds_lu.interp(lat=ds_wind.lat, lon=ds_wind.lon, method="linear")
            ds_lu_int[
                :, 0
            ] = 0  # todo improve this ugly fix. Maybe using xEMSF? Without specifying fill_value=0, nans are added at lon=0 because interpolation doesn't understand periodicity
            # exclude grid boxes where land use change is exactly zero (i.e., oceans and non-used lands)
            y = ds_lu_int.where(np.abs(ds_lu_int) > 0.01)
            x = ds_wind.where(np.abs(ds_lu_int) > 0.01)

            if analysis_type == "hotspots":
                wind_threshold = quantile(ds_wind, 0.95)
                x = x.where(ds_wind > wind_threshold)
                y = y.where(ds_wind > wind_threshold)

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
            ys = [ conditional_prob(x,y,x_thr=x_thr) for x_thr in np.arange(-1,2, .1)]
            ax2.plot(np.arange(-1,2, .1), ys, color="Purple")
        # ax[3].set_xlim(xmin=-0.5, xmax=2)
        ax[3].set_xlabel("Wind speed change [m/s]")
        plt.tight_layout()
        plt.savefig(
            "../plots/scatter/scatter_"
            + analysis_type
            + "_"
            + lu_variable
            + "_"
            + wind_type
            + ".jpeg",
            dpi=300,
        )
        plt.close("all")
