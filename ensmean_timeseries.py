from utils import *


f, ax = plt.subplots(ncols=4, figsize=(10, 5), sharey=True)
# plot others
for col_plot, experiment in enumerate(["historical", "rcp26", "rcp45", "rcp85"]):
    path = "../data/" + experiment + "/"
    # open wind
    ds_wind = open_datasets(sorted(glob.glob(path + "sfcWind*.nc")))
    print("wind opened")

    # plot 20y running means
    for ens_member in ds_wind.ensemble_member:
        tmp_wind = ds_wind["sfcWind"].sel({"ensemble_member": ens_member})
        tmp_wind.rolling(time=20, center=True).mean().dropna(dim="time").plot(
            ax=ax[col_plot], alpha=0.3
        )
    ds_wind["sfcWind"].rolling(time=20, center=True).mean().dropna(dim="time").mean(
        dim="ensemble_member"
    ).plot(ax=ax[col_plot], alpha=0.75, color="black")
    print("wind plotted")

    ax[col_plot].set_title(experiment)
ax[0].set_ylabel("10m wind speed [m/s]", fontsize=15)
plt.subplots_adjust(0.05, 0.05, 0.95, 0.9)
plt.savefig("../plots/box_timeseries_Europe.jpeg", dpi=300)
