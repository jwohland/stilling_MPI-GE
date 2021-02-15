import xarray as xr
from scipy.stats import linregress, norm, pearsonr, anderson, kstest
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import periodogram
from pyts.decomposition import SingularSpectrumAnalysis as SSA
import pandas as pd


def selbox(ds):
    # lats, lons = slice(45, 55), slice(10, 20)  # "Germany"
    lats, lons = slice(37.5, 60), slice(-10, 25)
    return ds.sel({"lat": lats, "lon": lons}).mean(dim=["lat", "lon"])


def ann_mean(ds):
    return ds.resample({"time": "1Y"}).mean()


def slope_if_significant(y, p_threshold=0.05):
    res = linregress(np.arange(20) / 10, y)
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


def open_datasets(filelist):
    ds = [ann_mean(selbox(xr.open_dataset(x, use_cftime=True))) for x in filelist]
    ds = xr.concat(
        ds,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in filelist]),
    )
    return ds


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
path = "../data/pi-control/"
ds_list = [
    ann_mean(selbox(xr.open_dataset(x, use_cftime=True)))
    for x in sorted(glob.glob(path + "*.nc"))
]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
ds_picontrol = xr.concat(ds_list, dim="time")


slopes = np.asarray(
    [
        slope_if_significant(ds_picontrol["sfcWind"][x : x + 20], p_threshold=5 / 100.0)
        for x in range(1980)
    ]
)
slopes_ts = xr.DataArray(
    slopes, dims="time", coords={"time": ds_picontrol["sfcWind"].time[:1980]}
)

f, ax = plt.subplots(figsize=(26, 8))
ds_picontrol["sfcWind"].plot(ax=ax, alpha=0.5)
ds_picontrol["sfcWind"].rolling(time=10, center=True).mean().dropna(dim="time").plot(
    ax=ax, color="grey"
)
ds_picontrol["sfcWind"].rolling(time=20, center=True).mean().dropna(dim="time").plot(
    ax=ax, color="black"
)

ds_picontrol["sfcWind"].where(slopes_ts > 0).plot.line(
    marker="o", linewidth=0, color="red", alpha=0.8, label="onset upward"
)
ds_picontrol["sfcWind"].where(slopes_ts < 0).plot.line(
    marker="o", linewidth=0, color="green", alpha=0.8, label="onset_downward"
)
plt.legend()
plt.ylabel("Wind speed [m/s]")
plt.title("MPI-GE PI control,Europe ")
plt.tight_layout()
plt.savefig("../plots/box_timeseries_picontrol_Europe_slopes.pdf")

# autocorrelation
f, ax = plt.subplots()
lags, corrs = np.arange(1, 50), []
for lag in lags:
    vals = ds_picontrol["sfcWind"].values
    corrs.append(pearsonr(vals[lag:], vals[:-lag])[0])
ax.plot(lags, corrs, marker="o", linewidth=0)
ax.set_xlabel("Lag [years]")
ax.set_ylabel("Autocorrelation (Lagged Pearson correlation coefficient)")
plt.suptitle("Wind autocorrelation Europe pi-control MPI-GE")
plt.savefig("../plots/sfcWind_autocorrelation.pdf")

# periodogramm
freq, Pxx = periodogram(vals)
f, ax = plt.subplots()
ax.scatter(freq, Pxx)
ax.set_xlabel("Frequency [1/y]")
ax.set_ylabel("Power spectral density")
ax.set_xlim(xmin=-0.001, xmax=0.1)
plt.suptitle("Wind periodogram Europe pi-control MPI-GE")
plt.savefig("../plots/periodogram.pdf")


# histogram
f, ax = plt.subplots()
ax.hist(ds_picontrol["sfcWind"].values, bins=100, density=True)
ares = anderson(ds_picontrol["sfcWind"].values)
if (ares.critical_values > ares.statistic).all():
    # H_0 (data come from normal distribution) can not be rejected at 1% level
    title_string = "Anderson says Gaussian (1% significance level) \n "
else:
    title_string = "Anderson says not Gaussian (1% significance level) \n "
ksres = kstest(ds_picontrol["sfcWind"].values, cdf=norm.cdf)
title_string += (
    "KS reports p-value of "
    + str(ksres.pvalue)
    + " and D="
    + str(np.round(ksres.statistic, 1))
)
ax.set_title(title_string)
ax.set_xlabel("sfcWind [m/s]")
ax.set_ylabel("PDF")
plt.savefig("../plots/values_histogram.pdf")


# PI-CONTROL plot trend histograms for different p-values
for p_threshold in [1, 5, 10, 100]:
    slopes = np.asarray(
        [
            slope_if_significant(
                ds_picontrol["sfcWind"][x : x + 20], p_threshold=p_threshold / 100.0
            )
            for x in range(1980)
        ]
    )
    frac_trends = np.round(slopes[np.isfinite(slopes)].size / 1980.0 * 100)
    frac_partoftrend = calc_frac_partoftrend(slopes)

    f, ax = plt.subplots()
    ax.axvline(x=-0.09, color="purple", label="Vautard et al. [2010]")
    ax.axvline(x=-0.1, color="orange", label="Zheng et al. [2019] 1978 - 2003")
    ax.axvline(
        x=0.1, color="orange", label="Zheng et al. [2019] 2004 - 2017"
    )  # from SI Fig. 4e visually derived

    n, bins, patches = ax.hist(slopes[np.isfinite(slopes)], bins=50, density=True)
    ax.set_ylabel("PDF")
    ax.set_xlabel(
        "Significant wind speed trends at "
        + str(100 - p_threshold)
        + "% level [m/s/decade]"
    )
    ax.set_title(
        str(frac_trends)
        + "% of 20y periods feature trends in MPI-GE pi-control \n "
        + str(frac_partoftrend)
        + "% of years belong to a 20y trend period"
    )
    plt.legend()

    # fit Gaussian to histogram without significance screening
    if p_threshold == 100:
        mu, std = norm.fit(slopes)
        ax.plot(bins, norm.pdf(bins, mu, std), color="red")

    plt.savefig("../plots/picontrol_wind_trends_Europe_" + str(p_threshold) + ".pdf")


# trend histograms in other periods  # todo unify with pi-control above
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
        frac_partoftrend = calc_frac_partoftrend(slopes)
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
        n, bins, patches = ax.hist(
            slopes[np.isfinite(slopes)],
            bins=50,
            density=True,
            color="darkorange",
            alpha=0.7,
            label="ensemble members"
        )
        ax.axvline(x=-0.09, color="purple", ls="--")
        textdic = {
            "horizontalalignment": "center",
            "verticalalignment": "center",
            "rotation": 90,
            "fontsize": 8,
        }
        ax.text(-0.083, n.max() *3./4, "Vautard et al. [2010]", textdic)
        ax.axvline(x=-0.1, color="purple", ls="--")
        ax.text(-0.107, n.max() *3./4, "Zheng et al. [2019] 1978 - 2003", textdic)
        ax.axvline(x=0.1, color="purple", ls="--")  # from SI Fig. 4e visually derived
        ax.text(0.093, n.max() *3./4, "Zheng et al. [2019] 2004 - 2017", textdic)
        ax2.hist(
            slopes_ensmean[np.isfinite(slopes_ensmean)],
            bins=50,
            density=True,
            color="darkgreen",
            alpha=0.7,
            label="ensemble mean"
        )
        ax.set_ylabel("PDF ensemble members", color="darkorange")
        ax2.set_ylabel("PDF ensemble mean", color="darkgreen")
        ax.set_xlabel(
            "Significant wind speed trends at "
            + str(100 - p_threshold)
            + "% level [m/s/decade]"
        )
        ax.set_title(
            experiment
            + ": "
            + str(frac_partoftrend)
            + "% of years belong to a 20y trend period"
        )
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
