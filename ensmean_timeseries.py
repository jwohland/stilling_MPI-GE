import xarray as xr
import pandas as pd
import glob
import matplotlib.pyplot as plt


def selbox(ds):
    # lats, lons = slice(45, 55), slice(10, 20)  # "Germany"
    lats, lons = slice(37.5, 60), slice(-10, 25)
    return ds.sel({"lat": lats, "lon": lons}).mean(dim=["lat", "lon"])


def ann_mean(ds):
    return ds.resample({"time": "1Y"}).mean()


def open_datasets(filelist):
    ds = [ann_mean(selbox(xr.open_dataset(x, use_cftime=True))) for x in filelist]
    ds = xr.concat(
        ds,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in filelist]),
    )
    return ds


f, ax = plt.subplots(ncols=5, nrows=2, figsize=(25, 12), sharey="row")

# plot pi-control
path = "../data/pi-control/"
ds_list = [
    ann_mean(selbox(xr.open_dataset(x, use_cftime=True)))
    for x in sorted(glob.glob(path + "*.nc"))
]
# use_cftime needed after 2200. Otherwise SerializationWarning is raised
ds_picontrol = xr.concat(ds_list, dim="time")

ds_picontrol["sfcWind"].rolling(time=20, center=True).mean().dropna(dim="time").plot(
    ax=ax[0, 0]
)
for i in range(5):
    ax[0, i].axhline(
        y=ds_picontrol["sfcWind"]
        .rolling(time=20, center=True)
        .mean()
        .dropna(dim="time")
        .max(),
        ls="--",
        color="blue",
        alpha=0.6,
    )

# plot others
for col_plot, experiment in enumerate(["historical", "rcp26", "rcp45", "rcp85"]):
    path = "../data/" + experiment + "/"

    windfiles = sorted(glob.glob(path + "sfcWind*.nc"))
    solarfiles = sorted(glob.glob(path + "rsds*.nc"))

    # There are 12 duplicate timesteps in a number of ensemble members that are excluded from this analysis
    if experiment == "rcp45":
        windfiles.remove(
            "../data/rcp45/sfcWind_Lmon_MPI-ESM_rcp45_r007i2005p3_200601-209912.nc"
        )
        windfiles.remove(
            "../data/rcp45/sfcWind_Lmon_MPI-ESM_rcp45_r021i2005p3_200601-209912.nc"
        )
        solarfiles.remove(
            "../data/rcp45/rsds_Lmon_MPI-ESM_rcp45_r007i2005p3_200601-209912.nc"
        )
        solarfiles.remove(
            "../data/rcp45/rsds_Lmon_MPI-ESM_rcp45_r021i2005p3_200601-209912.nc"
        )
    if experiment == "rcp26":
        windfiles.remove(
            "../data/rcp26/sfcWind_Lmon_MPI-ESM_rcp26_r055i2005p3_200601-209912.nc"
        )
        solarfiles.remove(
            "../data/rcp26/rsds_Lmon_MPI-ESM_rcp26_r055i2005p3_200601-209912.nc"
        )

    # open wind and solar
    ds_wind = open_datasets(windfiles)
    print("wind opened")
    ds_solar = open_datasets(solarfiles)
    print("solar opened")

    ax_wind, ax_solar = ax[0, col_plot + 1], ax[1, col_plot + 1]
    # plot 20y running means
    for ens_member in ds_wind.ensemble_member:
        tmp_wind = ds_wind["sfcWind"].sel({"ensemble_member": ens_member})
        tmp_wind.rolling(time=20, center=True).mean().dropna(dim="time").plot(
            ax=ax_wind, alpha=0.3
        )
    ds_wind["sfcWind"].rolling(time=20, center=True).mean().dropna(dim="time").mean(
        dim="ensemble_member"
    ).plot(ax=ax_wind, alpha=0.75, color="black")
    print("wind plotted")

    for ens_member in ds_solar.ensemble_member:
        tmp_solar = ds_solar["rsds"].sel({"ensemble_member": ens_member})
        tmp_solar.rolling(time=20, center=True).mean().dropna(dim="time").plot(
            ax=ax_solar, alpha=0.3
        )
    ds_solar["rsds"].rolling(time=20, center=True).mean().dropna(dim="time").mean(
        dim="ensemble_member"
    ).plot(ax=ax_solar, alpha=0.75, color="black")
    print("solar plotted")

    ax_wind.set_title(experiment)
ax[0, 0].set_ylabel("10m wind speed [m/s]", fontsize=15)
ax[1, 0].set_ylabel("rsds [W/m**2]", fontsize=15)
ax[1, 0].text(0.3, 143, "No data", fontsize=25)
for i in range(1,5):
    ax[0, i].set_xlim(ax[1,i].get_xlim())
# add titles and save
plt.suptitle("Mean over Europe, 20y running means", fontsize=25)
plt.subplots_adjust(0.05, 0.05, 0.95, 0.9)
plt.savefig("../plots/box_timeseries_Europe.pdf")


"""
    # open wind
    ds_wind = [ann_mean(selbox(xr.open_dataset(x))) for x in windfiles]
    ds_wind = xr.concat(
        ds_wind,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in windfiles]),
    )
    print("wind opened")
    
    # open solar
    ds_solar = [ann_mean(selbox(xr.open_dataset(x))) for x in solarfiles]
    ds_solar = xr.concat(
        ds_solar,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in solarfiles]),
    )
    print("solar opened")
    
    
    for ens_member in ds_wind.ensemble_member:
        tmp_wind = ds_wind["sfcWind"].sel({"ensemble_member": ens_member})
        tmp_wind.plot(ax=ax[0, 1], alpha=0.5)
        tmp_wind.rolling(time=20, center=True).mean().dropna(dim="time").plot(
            ax=ax[0, 2], alpha=0.5
        )
    ds_wind["sfcWind"].mean(dim="ensemble_member").plot(
        ax=ax[0, 1], alpha=0.75, color="black"
    )
    ds_wind["sfcWind"].rolling(time=20, center=True).mean().dropna(dim="time").mean(
        dim="ensemble_member"
    ).plot(ax=ax[0, 2], alpha=0.75, color="black")
    print("wind plotted")
    
    
    for ens_member in ds_solar.ensemble_member:
        tmp_solar = ds_solar["rsds"].sel({"ensemble_member": ens_member})
        tmp_solar.plot(ax=ax[1, 1], alpha=0.5)
        tmp_solar.rolling(time=20, center=True).mean().dropna(dim="time").plot(
            ax=ax[1, 2], alpha=0.5
        )
    ds_solar["rsds"].mean(dim="ensemble_member").plot(
        ax=ax[1, 1], alpha=0.75, color="black"
    )
    ds_solar["rsds"].rolling(time=20, center=True).mean().dropna(dim="time").mean(
        dim="ensemble_member"
    ).plot(ax=ax[1, 2], alpha=0.75, color="black")
    print("solar plotted")
"""


"""
plot 20y means individually
"""

f, ax = plt.subplots(ncols=2, nrows=100, figsize=(10, 200), sharex=True)

for i, ens_member in enumerate(ds_wind.ensemble_member):
    tmp_wind = ds_wind["sfcWind"].sel({"ensemble_member": ens_member})
    tmp_wind.rolling(time=20, center=True).mean().dropna(dim="time").plot(
        ax=ax[i, 0], alpha=0.5
    )
    tmp_solar = ds_solar["rsds"].sel({"ensemble_member": ens_member})
    tmp_solar.rolling(time=20, center=True).mean().dropna(dim="time").plot(
        ax=ax[i, 1], alpha=0.5
    )
plt.tight_layout()
plt.savefig("../plots/box_timeseries_individual.pdf", dpi=300)
