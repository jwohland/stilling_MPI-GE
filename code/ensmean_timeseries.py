import glob

import matplotlib.pyplot as plt

from utils import open_picontrol, open_datasets


def plot_ensemble_members_timeseries(path_to_data, path_to_plots):
    """
    Plots 20y running mean timeseries (ensemble mean and ensemble members) of surface wind speeds for the
    historical period and the 3 representative concentration pathways in the future.
    For comparison with the pre-industrial period, a dashed line representing the maximum during the pre-industrial
    period is added.
    :param path_to_data: Path to directory that contains scenarios (historical, rcps) as subdirectories
    :param path_to_plots: Path to plot directory
    :return:
    """
    ds_picontrol = open_picontrol(path_to_data)
    max_pi = float(
        ds_picontrol["sfcWind"].rolling(time=20).mean().dropna(dim="time").max()
    )

    f, ax = plt.subplots(ncols=4, figsize=(12, 4), sharey=True)
    # plot others
    for col_plot, experiment in enumerate(["historical", "rcp26", "rcp45", "rcp85"]):

        # open wind
        ds_wind = open_datasets(
            sorted(glob.glob(f"{path_to_data}/{experiment}/sfcWind*.nc"))
        )
        print("{experiment} wind opened")

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
        ax[col_plot].axhline(max_pi, ls="--", label="PI control maximum")
        ax[col_plot].set_title(experiment)
    ax[0].set_ylabel("10m wind speed [m/s]", fontsize=15)
    ax[3].legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{path_to_plots}/box_timeseries_Europe.jpeg", dpi=300)
