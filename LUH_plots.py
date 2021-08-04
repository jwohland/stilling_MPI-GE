import os
import warnings

from matplotlib import cm
import matplotlib.pyplot as plt
from utils import open_LUH_period, add_letters, sel_time

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy.crs as ccrs
    import cartopy


def open_LUH(path_to_data):
    # open time varying fields
    return open_LUH_period(f"{path_to_data}/LUHa.v1/", 1850, 2000)


def LUH_reference_period(all_LUH):
    return sel_time(all_LUH, slice("1850", "1860"))


def plot_absolute_land_use_change(path_to_data, path_to_plots):
    """
    plot absolute values
    """
    f, ax = plt.subplots(
        ncols=4,
        nrows=4,
        figsize=(12, 10),
        sharex=True,
        sharey=True,
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    ds = open_LUH(path_to_data)
    for j, start_year in enumerate([1850, 1890, 1940, 1990]):
        t_slice = slice(str(start_year), str(start_year + 10))
        for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
            ds_sliced_mean = sel_time(ds[var], t_slice)
            ds_sliced_mean.plot(
                ax=ax[j, i], vmin=0.01, vmax=1.01, levels=11, extend="neither"
            )
            ax[j, i].set_title(str(start_year) + " to " + str(start_year + 10))
            ax[j, i].add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
            ax[j, i].add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)

    plt.tight_layout()
    os.makedirs(f"{path_to_plots}/LUH1", exist_ok=True)
    plt.savefig(f"{path_to_plots}/LUH1/LUH_absolute.jpeg", dpi=300)


def plot_historical_LUH_change(path_to_data, path_to_plots):

    """
    plot ref and changes to ref
    """
    f, ax = plt.subplots(
        ncols=4,
        nrows=4,
        figsize=(12, 10),
        sharex=True,
        sharey=True,
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    ds = open_LUH(path_to_data)
    ref = LUH_reference_period(ds)
    for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
        ref[var].plot(ax=ax[0, i], vmin=0.01, vmax=1.01, levels=11, extend="neither")
        ax[0, i].set_title("1850 to 1860")

    for j, start_year in enumerate([1890, 1940, 1990]):
        t_slice = slice(str(start_year), str(start_year + 10))
        diff = sel_time(ds, t_slice) - ref
        for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
            diff[var].plot(
                ax=ax[j + 1, i],
                vmin=-1,
                vmax=1,
                # levels=11,
                extend="neither",
                cmap=cm.get_cmap("RdBu_r"),
            )
            ax[j + 1, i].set_title(f"Diff {start_year} to {start_year + 10}")

    for tmp_ax in ax.flatten():
        tmp_ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
        tmp_ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
    plt.tight_layout()
    os.makedirs(f"{path_to_plots}/LUH1", exist_ok=True)
    plt.savefig(f"{path_to_plots}/LUH1/LUH_change.jpeg", dpi=300)


def plot_future_LUH_change(path_to_data, path_to_plots):
    """
    plot ref and changes to ref in the future
    """
    ds = open_LUH(path_to_data)
    ref = LUH_reference_period(ds)
    # 2nd reference scenario end of 20th century
    ref_end = sel_time(ds, slice("1990", "2000"))

    for i_ref, reference_data in enumerate([ref, ref_end]):
        f, ax = plt.subplots(
            ncols=3,
            figsize=(8, 2.3),
            sharex=True,
            sharey=True,
            subplot_kw={"projection": ccrs.PlateCarree()},
        )
        cbar_ax = f.add_axes([0.2, 0.19, 0.6, 0.04])

        if i_ref == 0:
            vmax = 1
        else:
            vmax = 0.5
        for j, experiment in enumerate(
            ["IMAGE_rcp26", "MiniCAM_rcp45", "MESSAGE_rcp85"]
        ):
            print(experiment)
            ds_future = open_LUH_period(
                f"{path_to_data}/LUHa.v1/{experiment}", 2090, 2100
            )
            diff = sel_time(ds_future, slice("2090", "2100")) - reference_data

            var = "gothr+gsecd"
            if j != 0:
                diff[var].plot(
                    ax=ax[j],
                    vmin=-vmax,
                    vmax=vmax,
                    add_colorbar=False,
                    extend="both",
                    cmap=cm.get_cmap("RdBu"),
                )
            else:
                diff[var].plot(
                    ax=ax[j],
                    vmin=-vmax,
                    vmax=vmax,
                    add_colorbar=True,
                    cbar_ax=cbar_ax,
                    extend="both",
                    cmap=cm.get_cmap("RdBu"),
                    cbar_kwargs={
                        "orientation": "horizontal",
                        "label": "Difference in primary plus secondary land",
                    },
                )
            ax[j].set_title(experiment)

        for tmp_ax in ax.flatten():
            tmp_ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
            tmp_ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
        add_letters(ax, x=-0.03, y=1.04)
        plt.subplots_adjust(0.015, 0.12, 0.99, 0.99, 0.08, 0.08)
        os.makedirs(f"{path_to_plots}/LUH1", exist_ok=True)
        if i_ref == 0:
            plt.savefig(f"{path_to_plots}/LUH1/LUH_change_future.jpeg", dpi=300)
        else:
            plt.savefig(f"{path_to_plots}/LUH1/LUH_change_future_ref2000.jpeg", dpi=300)
