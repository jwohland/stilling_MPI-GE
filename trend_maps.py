# plot global map of wind speed trends over 1900 - 2010

from scipy.stats import linregress
import xarray as xr
import glob
import pandas as pd
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import warnings

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


path = "../data/historical/"

for t_slice in [
    slice("1900", "2005"),
    slice("1850", "1975"),
    slice("1945", "1975"),
    slice("1975", "2005"),
]:
    ds = ann_mean(xr.open_dataset(glob.glob("../data/historical/ensmean/*.nc")[0]))
    ds = sel_time(ds, t_slice)
    # ds = open_datasets_ensmean(sorted(glob.glob(path + "sfcWind*.nc")))

    nt, nlat, nlon = ds["sfcWind"].values.shape
    trends = np.zeros((nlat, nlon)) * np.nan
    for i_lat in range(nlat):
        print(i_lat)
        for i_lon in range(nlon):
            trend = linregress(
                range(ds["sfcWind"].shape[0]), ds["sfcWind"].values[:, i_lat, i_lon]
            )
            if trend.pvalue < 0.05:
                trends[i_lat, i_lon] = trend.slope

    trends = xr.DataArray(
        trends * 10, dims=["lat", "lon"], coords=[ds.coords["lat"], ds.coords["lon"]]
    )
    plt.close("all")
    ax = plot_field(
        trends,
        title="95% significant trends "
        + t_slice.start
        + " - "
        + t_slice.stop
        + " [m/s/decade]",
        vmin=-0.25,
        vmax=0.25,
        levels=11,
        cmap=cm.get_cmap("RdBu_r"),
    )
    plt.tight_layout()
    plt.savefig(
        "../plots/historical_trends_"
        + str(t_slice.start)
        + "_"
        + str(t_slice.stop)
        + ".jpeg",
        dpi=300,
    )


# 10y mean maps
ref = ann_mean(xr.open_dataset(glob.glob("../data/historical/ensmean/*.nc")[0]))
ref = sel_time(ref, slice("1850", "1859")).mean(dim="time")
for t_start in np.arange(1860, 2005, 10):
    ds = ann_mean(xr.open_dataset(glob.glob("../data/historical/ensmean/*.nc")[0]))
    t_slice = slice(
        str(t_start), str(t_start + 9)
    )  # slice includes first and last years, so this is a 10y period
    ds = sel_time(ds, t_slice).mean(dim="time")
    plt.close("all")
    ax = plot_field(
        ds["sfcWind"] - ref["sfcWind"],
        title="Diff mean surface wind speed "
        + t_slice.start
        + " - "
        + t_slice.stop
        + " [m/s]",
        cmap=cm.get_cmap("RdBu_r"),
        vmax=2,
        vmin=-2,
    )
    plt.tight_layout()
    plt.savefig(
        "../plots/map_windspeeds"
        + str(t_slice.start)
        + "_"
        + str(t_slice.stop)
        + ".jpeg",
        dpi=300,
    )


# global mean timeseries
f, ax = plt.subplots(ncols=2, sharey=True, figsize=(12,6))
ds_list = [
    ann_mean(xr.open_dataset(x, use_cftime=True))
    for x in sorted(glob.glob("../data/pi-control/*.nc"))
]
# use_cftime needed after 2200. Otherwise SerializationWarning is raised
ds_picontrol = xr.concat(ds_list, dim="time")
ds_picontrol = ds_picontrol.mean(dim=["lat", "lon"])
ax[0].hist(
    ds_picontrol["sfcWind"].values,
    bins=200,
    orientation="horizontal",
    density=True,
    label="pi-control",
    color="black",
    alpha=0.5,
)
ax[0].legend()
ax[0].set_xlabel("PDF")

for experiment in ["historical", "rcp26", "rcp45", "rcp85"]:
    ds = ann_mean(
        xr.open_dataset(glob.glob("../data/" + experiment + "/ensmean/*.nc")[0])
    )
    ds = ds.mean(dim=["lat", "lon"])
    ds["sfcWind"].plot(ax=ax[1], label=experiment)
ax[0].set_ylabel("Global mean 10m wind speed [m/s]")
ax[1].legend()
ax[1].set_ylabel("")
plt.tight_layout()
plt.savefig(
    "../plots/global_windspeeds.jpeg", dpi=300,
)
