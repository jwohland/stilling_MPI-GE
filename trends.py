from scipy.stats import linregress, norm, pearsonr, anderson, kstest
import numpy as np
from utils import *


def selHadISD(ds):
    """
    averages over all station locations used in Zeng et al. (2019) from the HadISD dataset
    :param ds:
    :return:
    """
    # load HadISD station list and limit to region of interest
    station_list = pd.read_excel(
        "../data/HadISD/Zeng_SIData2_HadISDv202.xlsx", usecols=["lons", "lats"]
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


def slope_if_significant(y, p_threshold=0.05, trend_length=20):
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
                    if (
                        y[i + j] == 0
                    ):  # if next timestep doesn't feature new trend increase weight
                        y[i] += 1
                    else:
                        break
                except IndexError:
                    y[i] += (
                        20 - j
                    )  # add remaing years to 20y at the end of the timeseries
                    break
    return np.round(y.sum() / (y.size + 19) * 100)


def plot_histo(slopes, ax, experiment, full_output=False, bins=50, trend_length=20):
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
    ax.set_xlabel(
        "Significant wind speed trends at "
        + str(100 - p_threshold)
        + "% level [m/s/decade]"
    )
    ax.set_title(
        experiment
        + ": "
        + str(int(frac_partoftrend))
        + "% of years belong to a "
        + str(trend_length)
        + "y trend period"
    )
    if full_output:
        return n, bins, frac_partoftrend
    else:
        return bins


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
print("Test completed succesfully")

# plot full timeseries and mark trends
ds_picontrol = open_picontrol()


slopes = np.asarray(
    [
        slope_if_significant(ds_picontrol["sfcWind"][x : x + 20], p_threshold=5 / 100.0)
        for x in range(1980)
    ]
)
slopes_ts = xr.DataArray(
    slopes, dims="time", coords={"time": ds_picontrol["sfcWind"].time[:1980]}
)

f, ax = plt.subplots(figsize=(12, 5))
ds_picontrol["sfcWind"].plot(ax=ax, alpha=0.5, label="Annual mean")
ds_picontrol["sfcWind"].rolling(time=20, center=True).mean().dropna(dim="time").plot(
    ax=ax, color="black", label="20y mean"
)
ds_picontrol["sfcWind"].where(slopes_ts > 0).plot.line(
    marker="o", linewidth=0, color="red", alpha=0.7, label="onset upward trend"
)
ds_picontrol["sfcWind"].where(slopes_ts < 0).plot.line(
    marker="o", linewidth=0, color="green", alpha=0.7, label="onset downward trend"
)
ax.legend(loc="upper right", ncol=4)
ax.set_xlabel("Year of pi-control simulation", fontsize=12)
ax.set_ylabel("European mean wind speed [m/s]", fontsize=12)
ax.set_title("")
ax.set_ylim(ymax=5.4)
ax.set_xlim(xmin=ds_picontrol.time[0].values, xmax=ds_picontrol.time[-1].values)
plt.tight_layout()
plt.savefig("../plots/timeseries_picontrol_Europe.jpeg", dpi=300)


"""
Plot trend histograms
"""
ds_list_HadISD = [
    selHadISD(ann_mean(xr.open_dataset(x, use_cftime=True)))
    for x in sorted(glob.glob("../data/pi-control/*.nc"))
]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
ds_picontrol_HadISD = xr.concat(ds_list_HadISD, dim="time")

# PI-CONTROL plot trend histograms for different p-values
for trend_length in [15, 20, 25]:
    for agg_type in [
        "HadISD",
        "box",
    ]:  # HadISD is sensitivity test with data averaged to European HadISD stations
        for p_threshold in [5, 100]:
            if agg_type == "box":
                ds_tmp = ds_picontrol.copy()
            else:
                ds_tmp = ds_picontrol_HadISD.copy()
            slopes = np.asarray(
                [
                    slope_if_significant(
                        ds_tmp["sfcWind"][x : x + trend_length],
                        p_threshold=p_threshold / 100.0,
                        trend_length=trend_length
                    )
                    for x in range(ds_tmp.time.size-trend_length)
                ]
            )
            f, ax = plt.subplots()
            bins = plot_histo(slopes, ax, "Pi-control", trend_length=trend_length)
            # fit Gaussian to histogram without significance screening
            if p_threshold == 100:
                mu, std = norm.fit(slopes)
                ax.plot(bins, norm.pdf(bins, mu, std), color="red")

            ax.set_ylabel("PDF")
            plt.tight_layout()
            if agg_type == "box":
                plt.savefig(
                    "../plots/picontrol_wind_trends_Europe_"
                    + str(p_threshold)
                    + "_"
                    + str(trend_length)
                    + "y.jpeg",
                    dpi=300,
                )
            else:
                plt.savefig(
                    "../plots/picontrol_HadISD_wind_trends_Europe_"
                    + str(p_threshold)
                    + "_"
                    + str(trend_length)
                    + "y.jpeg",
                    dpi=300,
                )

# PI-CONTROL trend histograms for CMIP6 ensemble
filelist = glob.glob("../data/CMIP6_annual/*.nc")
models = np.unique([x.split("/")[-1].split("_")[2] for x in filelist])

p_threshold = 5
CMIP6_histos = {}
CMIP6_bins = np.arange(-0.2, 0.2, 0.005)
for i, model in enumerate(models):
    print(str(int(i / len(models) * 100)) + "%")
    ds_list = [
        selbox(xr.open_dataset(x, use_cftime=True))
        for x in sorted(glob.glob("../data/CMIP6_annual/*" + model + "*.nc"))
    ]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
    ds_CMIP6 = xr.concat(ds_list, dim="time")
    slopes = np.asarray(
        [
            slope_if_significant(
                ds_CMIP6["sfcWind"][x : x + 20], p_threshold=p_threshold / 100.0
            )
            for x in range(ds_CMIP6.time.size - 20)
        ]
    )
    f, ax = plt.subplots()
    CMIP6_histos[model] = plot_histo(
        slopes, ax, "Pi-control " + model, full_output=True, bins=CMIP6_bins
    )
    ax.set_ylabel("PDF")
    plt.tight_layout()
    plt.savefig(
        "../plots/CMIP6/"
        + model
        + "_picontrol_wind_trends_Europe_"
        + str(p_threshold)
        + ".jpeg",
        dpi=300,
    )
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
ax.text(-0.083, n_mean.max() * 3.0 / 4, "Vautard et al. [2010] 1979 - 2008", textdic)
ax.axvline(x=-0.1, color="purple", ls="--")  # from SI Fig. 4e
ax.text(-0.107, n_mean.max() * 3.0 / 4, "Zeng et al. [2019] 1978 - 2003", textdic)
ax.axvline(x=0.11, color="purple", ls="--")  # from SI Fig. 4e
ax.text(0.103, n_mean.max() * 3.0 / 4, "Zeng et al. [2019] 2004 - 2017", textdic)
ax.set_ylabel("CMIP6 ensemble mean PDF")
ax.set_xlabel(
    "Significant wind speed trends at "
    + str(100 - p_threshold)
    + "% level [m/s/decade]"
)
ax.set_title(
    "Pi-control: "
    + str(int(frac_mean))
    + "% of years belong to a 20y trend period"
)
ax.set_xlim(xmin=-0.2, xmax=0.2)
plt.tight_layout()
plt.savefig(
    "../plots/CMIP6/Ensmean_picontrol_wind_trends_Europe_" + str(p_threshold) + ".jpeg",
    dpi=300,
)

# trend histograms in other periods
for experiment in ["historical", "rcp26", "rcp45", "rcp85"]:
    print(experiment)
    path = "../data/" + experiment + "/"

    windfiles = sorted(glob.glob(path + "sfcWind*.nc"))
    ds = open_datasets(windfiles)
    ds_ensmean = ann_mean(selbox(xr.open_dataset(glob.glob(path + "ensmean/*.nc")[0])))

    ds_internal = ds - ds_ensmean
    for p_threshold in [5, 100]:
        # calculate trend slopes in individual ens members
        slopes = []
        for ens_member in ds_internal.ensemble_member:
            da_internal = ds_internal["sfcWind"].sel({"ensemble_member": ens_member})
            slopes.extend(
                [
                    slope_if_significant(
                        da_internal[x : x + 20], p_threshold=p_threshold / 100.0
                    )
                    for x in range(da_internal.size - 20)
                ]
            )
        slopes = np.asarray(slopes)
        frac_trends = np.round(
            slopes[np.isfinite(slopes)].size / (da_internal.size - 20) * 100
        )

        # calculate trend slopes in ensemble mean
        slopes_ensmean = np.asarray(
            [
                slope_if_significant(
                    ds_ensmean["sfcWind"][x : x + 20], p_threshold=p_threshold / 100.0
                )
                for x in range(ds_ensmean.time.size - 20)
            ]
        )

        # plotting
        f, ax = plt.subplots()
        ax2 = ax.twinx()
        bins = plot_histo(slopes, ax, experiment)
        ax2.hist(
            slopes_ensmean[np.isfinite(slopes_ensmean)],
            bins=50,
            density=True,
            color="darkgreen",
            alpha=0.7,
            label="ensemble mean",
        )
        ax.set_ylabel("PDF ensemble members", color="darkorange")
        ax2.set_ylabel("PDF ensemble mean", color="darkgreen")

        # fit Gaussian to histogram without significance screening
        if p_threshold == 100:
            mu, std = norm.fit(slopes)
            ax.plot(bins, norm.pdf(bins, mu, std), color="red")
        plt.tight_layout()
        plt.savefig(
            "../plots/"
            + str(experiment)
            + "_wind_trends_Europe_"
            + str(p_threshold)
            + "_all.jpeg",
            dpi=300,
        )
