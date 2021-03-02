from utils import *
import numpy as np
from scipy.signal import periodogram
from scipy.stats import norm, pearsonr, anderson, kstest

path = "../data/pi-control/"
ds_list = [
    ann_mean(selbox(xr.open_dataset(x, use_cftime=True)))
    for x in sorted(glob.glob(path + "*.nc"))
]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
ds_picontrol = xr.concat(ds_list, dim="time")


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
