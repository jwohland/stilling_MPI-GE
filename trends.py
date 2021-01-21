import xarray as xr
from scipy.stats import linregress
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, pearsonr, anderson, kstest
from scipy.signal import periodogram
from pyts.decomposition import SingularSpectrumAnalysis as SSA

# compute trends over 20y locally and globally (and over land)


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
lags, corrs = np.arange(1,50), []
for lag in lags:
    vals = ds_picontrol["sfcWind"].values
    corrs.append(pearsonr(vals[lag:], vals[:-lag])[0])
ax.plot(lags, corrs, marker="o", linewidth=0)
ax.set_xlabel("Lag [years]")
ax.set_ylabel("Autocorrelation (Lagged Pearson correlation coefficient)")
plt.suptitle("Wind autocorrelation Europe pi-control MPI-GE")
plt.savefig("../plots/sfcWind_autocorrelation.pdf")

# periodogramm
freq, Pxx =  periodogram(vals)
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
title_string += "KS reports p-value of " +str(ksres.pvalue) + " and D=" + str(np.round(ksres.statistic, 1))
ax.set_title(title_string)
ax.set_xlabel("sfcWind [m/s]")
ax.set_ylabel("PDF")
plt.savefig("../plots/values_histogram.pdf")

# SSA
#groups = [np.arange(i, i + 5) for i in range(0, 11, 5)]  # no idea what that does
ssa = SSA(window_size=20)
X_ssa = ssa.fit_transform(ds_picontrol["sfcWind"].values.reshape(1,-1))  # reshaping needed because ssa expects 2D imput
f, ax = plt.subplots()
for i in range(X_ssa.shape[0]):
    ax.plot(X_ssa[i, :], label="SSA " + str(i))
plt.legend()
plt.suptitle("SSA analysis of sfcWind Europe, window_size=20")
plt.savefig("../plots/SSA.pdf")

# plot trend histograms for different p-values
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

    plt.savefig("../plots/historical_wind_trends_Europe_" + str(p_threshold) + ".pdf")
