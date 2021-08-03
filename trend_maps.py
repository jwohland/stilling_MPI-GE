# plot global map of wind speed trends over 1900 - 2010
import warnings

from scipy.stats import linregress
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy
    import cartopy.crs as ccrs

from utils import (
    ensemble_mean_wind_speed,
    reference_ensemble_mean_wind_speed,
    open_picontrol,
    mean_sliced_annual_mean,
    open_landmask,
    add_letters,
)


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


def plot_historical_trends(path_to_data, path_to_plots):
    for t_slice in [
        slice("1900", "2005"),
        slice("1850", "1975"),
        slice("1945", "1975"),
        slice("1975", "2005"),
    ]:
        ds = mean_sliced_annual_mean(
            ensemble_mean_wind_speed(path_to_data, t_slice, mean_time=False)
        )

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
            trends * 10,
            dims=["lat", "lon"],
            coords=[ds.coords["lat"], ds.coords["lon"]],
        )
        plt.close("all")
        plot_title = (
            f"95% significant trends {t_slice.start} - {t_slice.stop} [m/s/decade]"
        )
        ax = plot_field(
            trends,
            title=plot_title,
            vmin=-0.25,
            vmax=0.25,
            levels=11,
            cmap=cm.get_cmap("RdBu_r"),
        )
        plt.tight_layout()
        fig_path = f"{path_to_plots}/maps/historical/historical_trends_{t_slice.start}_{t_slice.stop}.jpeg"
        plt.savefig(fig_path, dpi=300)


def plot_windspeed_maps(path_to_data, path_to_plots):
    # 10y mean maps
    ref = reference_ensemble_mean_wind_speed(path_to_data)
    experiment_start_time = {
        "historical": np.arange(1860, 2005, 10),
        "1pCO2": np.arange(1860, 2005, 10),
        "rcp26": np.arange(2010, 2100, 10),
        "rcp45": np.arange(2010, 2100, 10),
        "rcp85": np.arange(2010, 2100, 10),
    }
    for experiment, start_times in experiment_start_time.items():
        for t_start in start_times:
            # slice includes first and last years, so this is a 10y period
            t_slice = slice(str(t_start), str(t_start + 9))
            ds = mean_sliced_annual_mean(ensemble_mean_wind_speed(experiment), t_slice)
            plt.close("all")
            plot_title = f"{experiment}\nDiff mean surface wind speed {t_slice.start} - {t_slice.stop} [m/s]"
            ax = plot_field(
                ds["sfcWind"] - ref["sfcWind"],
                title=plot_title,
                cmap=cm.get_cmap("RdBu_r"),
                vmax=2,
                vmin=-2,
            )
            plt.tight_layout()
            fig_path = f"{path_to_plots}/maps/{experiment}/map_windspeeds/{t_slice.start}_{t_slice.stop}.jpeg"
            plt.savefig(fig_path, dpi=300)


def plot_global_windspeeds(path_to_data, path_to_plots):
    # global mean timeseries

    # prep data
    landmask = open_landmask(path_to_data)
    ds_picontrol = open_picontrol(path_to_data, spatial_averaging=False)

    onshore_mean = ds_picontrol.where(np.isfinite(landmask)).mean(dim=["lat", "lon"])

    # plotting
    f, ax = plt.subplots(
        ncols=2, sharey=True, figsize=(8, 4), gridspec_kw={"width_ratios": [1, 2]}
    )
    ax[0].hist(
        onshore_mean["sfcWind"].values,
        bins=200,
        orientation="horizontal",
        density=True,
        label="pi-control",
        color="black",
        alpha=0.5,
    )
    ax[0].legend(loc="lower left")
    ax[0].set_xlabel("PDF")

    for experiment in ["historical", "rcp26", "rcp45", "rcp85"]:
        ds = ensemble_mean_wind_speed(experiment)
        ds = ds.where(np.isfinite(landmask)).mean(dim=["lat", "lon"])
        ds["sfcWind"].plot(ax=ax[1], label=experiment)
    ax[0].set_ylabel("Global mean onshore 10m wind speed [m/s]")
    ax[1].legend(loc="lower left")
    ax[1].set_ylabel("")
    ax[1].set_xlabel("")

    # add max & min lines for pi-control
    for i in range(2):
        ax[i].axhline(
            onshore_mean["sfcWind"].min().values, ls="--", color="black", alpha=0.5
        )
        ax[i].axhline(
            onshore_mean["sfcWind"].max().values, ls="--", color="black", alpha=0.5
        )
    ax[0].text(
        x=20,
        y=onshore_mean["sfcWind"].min().values + 0.01,
        s="Minimum",
        ha="center",
        va="center",
    )
    ax[0].text(
        x=20,
        y=onshore_mean["sfcWind"].max().values + 0.01,
        s="Maximum",
        ha="center",
        va="center",
    )
    add_letters(ax)
    plt.tight_layout()
    plt.savefig(f"{path_to_plots}/global_windspeeds.jpeg", dpi=300)
