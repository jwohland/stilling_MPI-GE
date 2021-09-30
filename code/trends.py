import warnings
import glob
import os

from scipy.stats import linregress, norm
import xarray as xr
import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy.crs as ccrs
    import cartopy.mpl.geoaxes

from utils import (
    ensemble_mean_wind_speed,
    annual_mean,
    add_letters,
    open_picontrol,
    selbox,
    open_datasets,
    LONMAX,
    LONMIN,
    LATMIN,
    LATMAX,
)

P_THRESHOLD = 5


def selHadISD(ds, path_to_data):
    """
    averages over all station locations used in Zeng et al. (2019) from the HadISD dataset
    :param ds:
    :return:
    """
    # load HadISD station list and limit to region of interest
    station_list = pd.read_excel(
        f"{path_to_data}/HadISD/Zeng_SIData2_HadISDv202.xlsx",
        usecols=["lons", "lats"], sheet_name="stations"
    )
    station_list = station_list.where(
        (station_list.lons < LONMAX)
        & (station_list.lons > LONMIN)
        & (station_list.lats > LATMIN)
        & (station_list.lats < LATMAX)
    ).dropna()
    # interpolate input data to HadISD stations and return the station mean
    ds_stations = []
    for index, row in station_list.iterrows():
        ds_stations.append(ds.interp(lat=row["lats"], lon=row["lons"], method="linear"))
    return xr.concat(ds_stations, dim="station_number").mean(dim="station_number")


def slope_if_significant(y, p_threshold=P_THRESHOLD, trend_length=20):
    """
    Calculates the slope by fitting a linear trend of length trend_length to the timeseries y.
    If the slope is not statistically significant at the p-values provided
    as p_threhold, the function returns nan.
    :param y:
    :param p_threshold:
    :param trend_length:
    :return:
    """
    p_threshold = p_threshold / 100  # from percentage to fraction
    res = linregress(np.arange(trend_length) / 10, y)
    if res[3] < p_threshold:
        return res[0]
    else:
        return np.nan


def calc_frac_partoftrend(y):
    """
    computes the number of timesteps that are part of a 20 timestep trend
    :param y:
    :return:
    """
    y = y.copy()
    y[np.isfinite(y)] = 1  # all values are 1
    y[np.isnan(y)] = 0
    for i in range(y.size):
        if y[i] == 1:
            for j in range(1, 20):
                try:
                    # if next timestep doesn't feature new trend increase weight
                    if y[i + j] == 0:
                        y[i] += 1
                    else:
                        break
                except IndexError:
                    # add remaining years to 20y at the end of the timeseries
                    y[i] += 20 - j
                    break
    return np.round(y.sum() / (y.size + 19) * 100)


def test_calc_frac_partoftrend():
    # test frac_partoftrend
    test_array = (
        np.zeros(81) * np.nan
    )  # 81 year slope timeseries corresponds to 100y input data
    test_array[3] = 3
    assert calc_frac_partoftrend(test_array) == 20.0 / 100 * 100
    test_array[-1] = 2
    assert calc_frac_partoftrend(test_array) == 40.0 / 100 * 100
    test_array[4] = 1
    assert calc_frac_partoftrend(test_array) == 41.0 / 100 * 100
    test_array[:] = 2
    assert calc_frac_partoftrend(test_array) == 100.0
    print("Test of function `calc_frac_partoftrend` completed succesfully")


def plot_histo(
    slopes,
    ax,
    experiment,
    full_output=False,
    bins=50,
    trend_length=20,
    p_threshold=P_THRESHOLD,
):
    """
    Plots histogram of wind speed trends that are significant at a given p value threshold
    along with the observation-based estimates taken from earlier studies.
    :param slopes:
    :param ax:
    :param experiment:
    :param full_output:
    :param bins:
    :param trend_length:
    :param p_threshold:
    :return:
    """
    n, bins, patches = ax.hist(
        slopes[np.isfinite(slopes)],
        bins=bins,
        density=True,
        color="darkorange",
        alpha=0.7,
    )
    ax.set_xlim(xmin=-0.2, xmax=0.2)
    textdic = {
        "horizontalalignment": "center",
        "verticalalignment": "center",
        "rotation": 90,
        "fontsize": 8,
    }
    ax.axvline(x=-0.09, color="purple", ls="--")  # p.1, 2nd column, 1st paragraph
    ax.text(-0.083, n.max() * 3.0 / 4, "Vautard et al. [2010] 1979 - 2008", textdic)
    ax.axvline(x=-0.1, color="purple", ls="--")  # from SI Fig. 4e
    ax.text(-0.107, n.max() * 3.0 / 4, "Zeng et al. [2019] 1978 - 2003", textdic)
    ax.axvline(x=0.11, color="purple", ls="--")  # from SI Fig. 4e
    ax.text(0.103, n.max() * 3.0 / 4, "Zeng et al. [2019] 2004 - 2017", textdic)
    frac_partoftrend = calc_frac_partoftrend(slopes)
    xlabel = f"Significant wind speed trends at {100 - p_threshold}% level [m/s/decade]"
    ax.set_xlabel(xlabel, fontsize=12)
    plot_title = f"{experiment}: {int(frac_partoftrend)}% of years belong to a {trend_length}y trend period"
    ax.set_title(plot_title)
    if full_output:
        return n, bins, frac_partoftrend
    else:
        return bins


def plot_full_timeseries_with_trend_marks(path_to_data, path_to_plots):
    """
    Plots annual and 20y running-mean wind speed timeseries during pre-industrial control.
    Markers or red and blue color denote onsets of significant 20y trend periods.
    A map of the considered domain is added.
    :param path_to_data:
    :param path_to_plots:
    :return:
    """
    # plot full timeseries and mark trends
    ds_picontrol = open_picontrol(path_to_data)

    slopes = np.asarray(
        [
            slope_if_significant(
                ds_picontrol["sfcWind"][x : x + 20], p_threshold=P_THRESHOLD
            )
            for x in range(1980)
        ]
    )
    slopes_ts = xr.DataArray(
        slopes, dims="time", coords={"time": ds_picontrol["sfcWind"].time[:1980]}
    )

    # plot slopes and mark trend onsets
    f, ax = plt.subplots(figsize=(12, 5))
    ds_picontrol["sfcWind"].plot(ax=ax, alpha=0.5, label="Annual mean")
    ds_picontrol["sfcWind"].rolling(time=20, center=True).mean().dropna(
        dim="time"
    ).plot(ax=ax, color="black", label="20y mean")
    ds_picontrol["sfcWind"].where(slopes_ts > 0).plot.line(
        marker="o", linewidth=0, color="red", alpha=0.7, label="onset upward trend"
    )
    ds_picontrol["sfcWind"].where(slopes_ts < 0).plot.line(
        marker="o", linewidth=0, color="green", alpha=0.7, label="onset downward trend"
    )
    # add inset with map of focus region
    axins = inset_axes(
        ax,
        width="10%",
        height="20%",
        loc="upper left",
        axes_class=cartopy.mpl.geoaxes.GeoAxes,
        axes_kwargs=dict(map_projection=ccrs.PlateCarree()),
    )
    axins.set_extent((LONMIN - 1, LONMAX + 1, LATMIN - 1.5, LATMAX + 0.5))
    axins.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
    axins.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.15)

    axins.add_patch(
        mpatches.Rectangle(
            xy=[LONMIN, LATMIN],
            width=LONMAX - LONMIN,
            height=LATMAX - LATMIN,
            facecolor="blue",
            alpha=0.2,
            transform=ccrs.PlateCarree(),
        )
    )
    axins.outline_patch.set_visible(False)
    ax.legend(loc="upper right", ncol=2)
    ax.set_xlabel("Year of pi-control simulation", fontsize=12)
    ax.set_ylabel("European mean wind speed [m/s]", fontsize=12)
    ax.set_title("")
    ax.set_ylim(ymax=5.42)
    ax.set_xlim(xmin=ds_picontrol.time[0].values, xmax=ds_picontrol.time[-1].values)
    plt.tight_layout()
    plt.savefig(f"{path_to_plots}/timeseries_picontrol_Europe.jpeg", dpi=300)
    plt.close("all")


def plot_trend_histograms(path_to_data, path_to_plots):
    """
    Plots trend histograms for different combinations of
        - trend lengths (15, 20, 25 years)
        - aggregation types (box average or interpolated to HadISD station locations)
        - significance levels (p values of 0.05, 0.1, 0.15, 1)
    A Gaussian is fitted to those plots where no significance screening is applied (i.e. p=1)
    :param path_to_data:
    :param path_to_plots:
    :return:
    """
    ds_picontrol = open_picontrol(path_to_data)
    ds_list_HadISD = [
        selHadISD(annual_mean(xr.open_dataset(x, use_cftime=True)), path_to_data)
        for x in sorted(glob.glob(f"{path_to_data}/pi-control/*.nc"))
    ]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
    ds_picontrol_HadISD = xr.concat(ds_list_HadISD, dim="time")

    # PI-CONTROL plot trend histograms for different p-values
    for trend_length in [15, 20, 25]:
        # HadISD is sensitivity test with data averaged to European HadISD stations
        for agg_type in ["HadISD", "box"]:
            for p_threshold in [5, 10, 15, 100]:
                if agg_type == "box":
                    ds_tmp = ds_picontrol.copy()
                else:
                    ds_tmp = ds_picontrol_HadISD.copy()
                slopes = np.asarray(
                    [
                        slope_if_significant(
                            ds_tmp["sfcWind"][x : x + trend_length],
                            p_threshold=p_threshold,
                            trend_length=trend_length,
                        )
                        for x in range(ds_tmp.time.size - trend_length)
                    ]
                )
                f, ax = plt.subplots()
                bins = plot_histo(
                    slopes,
                    ax,
                    "Pi-control",
                    trend_length=trend_length,
                    p_threshold=p_threshold,
                )
                # fit Gaussian to histogram without significance screening
                if p_threshold == 100:
                    mu, std = norm.fit(slopes)
                    ax.plot(bins, norm.pdf(bins, mu, std), color="red")

                ax.set_ylabel("MPI-GE PDF", fontsize=12)
                add_letters(ax)
                plt.tight_layout()
                if agg_type == "box":
                    fig_path = f"{path_to_plots}/picontrol_wind_trends_Europe_{p_threshold}_{trend_length}y.jpeg"
                else:
                    fig_path = f"{path_to_plots}/picontrol_HadISD_wind_trends_Europe_{p_threshold}_{trend_length}y.jpeg"
                plt.savefig(fig_path, dpi=300)
                plt.close("all")


def plot_pi_control_cmip6_trend_histograms(path_to_data, path_to_plots):
    # PI-CONTROL trend histograms for CMIP6 ensemble
    filelist = glob.glob(f"{path_to_data}/CMIP6_annual/*.nc")
    models = np.unique([x.split("/")[-1].split("_")[2] for x in filelist])

    CMIP6_histos = {}
    CMIP6_bins = np.arange(-0.2, 0.2, 0.005)
    for i, model in enumerate(models):
        print(str(int(i / len(models) * 100)) + "%")
        ds_list = [
            selbox(xr.open_dataset(x, use_cftime=True))
            for x in sorted(glob.glob(f"{path_to_data}/CMIP6_annual/*{model}*.nc"))
        ]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
        ds_CMIP6 = xr.concat(ds_list, dim="time")
        slopes = np.asarray(
            [
                slope_if_significant(
                    ds_CMIP6["sfcWind"][x : x + 20], p_threshold=P_THRESHOLD
                )
                for x in range(ds_CMIP6.time.size - 20)
            ]
        )
        f, ax = plt.subplots()
        CMIP6_histos[model] = plot_histo(
            slopes,
            ax,
            "Pi-control " + model,
            full_output=True,
            bins=CMIP6_bins,
            p_threshold=P_THRESHOLD,
        )
        ax.set_ylabel("PDF")
        plt.tight_layout()
        os.makedirs(f"{path_to_plots}/CMIP6", exist_ok=True)
        fig_path = f"{path_to_plots}/CMIP6/{model}_picontrol_wind_trends_Europe_{P_THRESHOLD}.jpeg"
        plt.savefig(fig_path, dpi=300)
        plt.close("all")

    # ensemble mean histo
    del CMIP6_histos["EC-Earth3-CC"]  # has no data
    del CMIP6_histos["AWI-ESM-1-1-LR"]  # only has negative trends of -0.07 m/s/dec
    df_CMIP6 = pd.DataFrame(CMIP6_histos, index=["n", "bins", "fracoftrends"])
    n_mean, frac_mean = df_CMIP6.loc["n"].mean(), df_CMIP6.loc["fracoftrends"].mean()
    f, ax = plt.subplots()
    ax.bar(CMIP6_bins[1:], n_mean, width=0.005, color="Darkorange", alpha=0.7)
    textdic = {
        "horizontalalignment": "center",
        "verticalalignment": "center",
        "rotation": 90,
        "fontsize": 8,
    }
    ax.axvline(x=-0.09, color="purple", ls="--")  # p.1, 2nd column, 1st paragraph
    ax.text(
        -0.083, n_mean.max() * 3.0 / 4, "Vautard et al. [2010] 1979 - 2008", textdic
    )
    ax.axvline(x=-0.1, color="purple", ls="--")  # from SI Fig. 4e
    ax.text(-0.107, n_mean.max() * 3.0 / 4, "Zeng et al. [2019] 1978 - 2003", textdic)
    ax.axvline(x=0.11, color="purple", ls="--")  # from SI Fig. 4e
    ax.text(0.103, n_mean.max() * 3.0 / 4, "Zeng et al. [2019] 2004 - 2017", textdic)
    ax.set_ylabel("CMIP6 ensemble mean PDF", fontsize=12)
    xlabel = f"Significant wind speed trends at {100 - P_THRESHOLD}% level [m/s/decade]"
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_title(f"Pi-control: {int(frac_mean)}% of years belong to a 20y trend period")
    ax.set_xlim(xmin=-0.2, xmax=0.2)
    add_letters(ax, letter_offset=1)
    plt.tight_layout()
    os.makedirs(f"{path_to_plots}/CMIP6", exist_ok=True)
    fig_path = (
        f"{path_to_plots}/CMIP6/Ensmean_picontrol_wind_trends_Europe_{P_THRESHOLD}.jpeg"
    )
    plt.savefig(fig_path, dpi=300)
    plt.close("all")


def plot_experiment_trend_histograms(path_to_data, path_to_plots):
    # trend histograms in other periods
    for letter_index, experiment in enumerate(
        ["historical", "rcp26", "rcp45", "rcp85"]
    ):
        print(experiment)

        windfiles = sorted(glob.glob(f"{path_to_data}/{experiment}/sfcWind*.nc"))
        ds = open_datasets(windfiles)

        ds_ensmean = annual_mean(
            selbox(ensemble_mean_wind_speed(path_to_data, experiment))
        )

        # Get internal variability as wind speed minus ensemble mean wind speed
        ds_internal = ds - ds_ensmean
        for p_threshold in [5, 100]:
            # calculate trend slopes in individual ens members
            slopes = []
            for ens_member in ds_internal.ensemble_member:
                da_internal = ds_internal["sfcWind"].sel(
                    {"ensemble_member": ens_member}
                )
                slopes.extend(
                    [
                        slope_if_significant(
                            da_internal[x : x + 20], p_threshold=p_threshold
                        )
                        for x in range(da_internal.size - 20)
                    ]
                )
            slopes = np.asarray(slopes)

            # calculate trend slopes in ensemble mean
            slopes_ensmean = np.asarray(
                [
                    slope_if_significant(
                        ds_ensmean["sfcWind"][x : x + 20], p_threshold=p_threshold
                    )
                    for x in range(ds_ensmean.time.size - 20)
                ]
            )

            # plotting
            f, ax = plt.subplots()
            ax2 = ax.twinx()
            bins = plot_histo(slopes, ax, experiment, p_threshold=p_threshold)
            ax2.hist(
                slopes_ensmean[np.isfinite(slopes_ensmean)],
                bins=50,
                density=True,
                color="darkgreen",
                alpha=0.7,
                label="ensemble mean",
            )
            ax.set_ylabel("PDF ensemble members", color="darkorange", fontsize=12)
            ax2.set_ylabel("PDF ensemble mean", color="darkgreen", fontsize=12)

            # fit Gaussian to histogram without significance screening
            if p_threshold == 100:
                mu, std = norm.fit(slopes)
                ax.plot(bins, norm.pdf(bins, mu, std), color="red")
            add_letters(ax, letter_offset=letter_index)
            plt.tight_layout()
            fig_path = f"{path_to_plots}/{experiment}_wind_trends_Europe_{p_threshold}_all.jpeg"
            plt.savefig(fig_path, dpi=300)
            plt.close("all")