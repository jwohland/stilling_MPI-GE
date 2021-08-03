import warnings
import pickle
import glob

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xesmf as xe
from scipy.stats import spearmanr, pearsonr

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy
    import cartopy.crs as ccrs

from utils import (
    add_letters,
    mean_sliced_annual_mean,
    open_LUH_period,
    reference_ensemble_mean_wind_speed,
    open_landmask,
)


def quantile(da, q):  # TODO: is this function even needed?
    """Get quantile `q` of finite data in data array `da`"""
    tmp = da.values
    tmp = tmp[np.isfinite(tmp)]
    return np.quantile(tmp, q)


def conditional_prob(x, y, x_thr):
    """
    Calculate conditional probability that `y` is less than zero,
    given that `x` is greater than a given threshold `x_thr`
    """
    joint_prob = y[(x > x_thr) & (y < 0)].size
    single_prob = y[x > x_thr].size
    if single_prob == 0:
        cond_prob = 0
    else:
        cond_prob = joint_prob / single_prob
    return cond_prob


def experiment_wind_speed_components_and_luh(path_to_data, path_to_cache=None):
    path_to_LUH1 = f"{path_to_data}/LUHa.v1"
    path_to_experiments = f"{path_to_data}/{{experiment}}/ensmean/*.nc"
    path_to_1p_CO2 = f"{path_to_data}/1pCO2/ensmean/*.nc"
    experiment_params = {
        "historical": {
            "end": slice("1990", "2000"),
            "equiv": slice("1880", "1890"),
            "LUH1_path": f"{path_to_LUH1}",
        },
        "rcp85": {
            "end": slice("2090", "2100"),
            "equiv": slice("1970", "1980"),
            "LUH1_path": f"{path_to_LUH1}/MESSAGE_rcp85",
        },
        "rcp45": {
            "end": slice("2090", "2100"),
            "equiv": slice("1910", "1920"),
            "LUH1_path": f"{path_to_LUH1}/MiniCAM_rcp45",
        },
        "rcp26": {
            "end": slice("2090", "2100"),
            "equiv": slice("1890", "1900"),
            "LUH1_path": f"{path_to_LUH1}/IMAGE_rcp26",
        },
    }
    if path_to_cache is not None and "attribution_maps.pickle" in os.listdir(
        path_to_cache
    ):
        with open(f"{path_to_cache}/attribution_maps.pickle", "rb") as handle:
            ds_dict = pickle.load(handle)
    else:
        # TODO: better name than 'ref'
        ref = reference_ensemble_mean_wind_speed(path_to_data)
        ds_dict = {}
        CO2_ref = xr.open_dataset(glob.glob(path_to_1p_CO2)[0])
        LUH_ref = open_LUH_period(path_to_LUH1, 1850, 1860).mean(dim="time")
        for experiment in ["historical", "rcp26", "rcp45", "rcp85"]:
            ds_dict[experiment] = {}
            experiment_ensemble_mean = mean_sliced_annual_mean(
                xr.open_dataset(
                    glob.glob(path_to_experiments.format(experiment=experiment))[0]
                ),
                experiment_params[experiment]["end"],
            )
            # full (dynamical + land use) difference
            ds_dict[experiment]["Full Diff"] = experiment_ensemble_mean - ref

            # only CO2 forcing
            ds_dict[experiment]["Dyn. Diff"] = (
                mean_sliced_annual_mean(CO2_ref, experiment_params[experiment]["equiv"])
                - ref
            )
            # TODO: add comment on what is happening here
            ds_dict[experiment]["Full - Dyn."] = (
                ds_dict[experiment]["Full Diff"] - ds_dict[experiment]["Dyn. Diff"]
            )
            # TODO: add comment on what is happening here
            ds_dict[experiment]["LUH Diff"] = (
                open_LUH_period(
                    experiment_params[experiment]["LUH1_path"],
                    int(experiment_params[experiment]["end"].start),
                    int(experiment_params[experiment]["end"].stop),
                ).mean(dim="time")
                - LUH_ref
            )
            # TODO: save these as xarrays (which could even be packaged on zenodo)?
            # E.g. a dataset for wind speed info, with different components as dataarrays and experiments as a coordinate,
            # and then another dataset for LUH1 data
            if path_to_cache is not None:
                with open(f"{path_to_cache}/attribution_maps.pickle", "wb") as handle:
                    pickle.dump(ds_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return ds_dict


def plot_attribution_maps(wind_speed_components_and_luh_dict, path_to_plots):
    # prep Figure
    f, ax = plt.subplots(
        ncols=4,
        nrows=4,
        figsize=(9, 6),
        sharex=True,
        sharey=True,
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    cbar1_ax = f.add_axes([0.05, 0.11, 0.6, 0.02])
    cbar2_ax = f.add_axes([0.8, 0.11, 0.175, 0.02])

    for row, experiment in enumerate(wind_speed_components_and_luh_dict.keys()):
        plt.text(
            -0.1,
            0.5,
            experiment,
            rotation=90,
            fontdict={"fontsize": 12},
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax[row, 0].transAxes,
        )
        for col, field in enumerate(
            wind_speed_components_and_luh_dict[experiment].keys()
        ):
            if row + col == 0:
                wind_speed_components_and_luh_dict[experiment]["Full Diff"][
                    "sfcWind"
                ].plot(
                    ax=ax[row, col],
                    vmin=-1.5,
                    vmax=1.5,
                    add_colorbar=True,
                    cbar_ax=cbar1_ax,
                    extend="both",
                    cmap=cm.get_cmap("RdBu_r"),
                    cbar_kwargs={
                        "orientation": "horizontal",
                        "label": "Difference in decadal mean wind speed [m/s]",
                    },
                )
            else:
                if field != "LUH Diff":
                    wind_speed_components_and_luh_dict[experiment][field][
                        "sfcWind"
                    ].plot(
                        ax=ax[row, col],
                        vmin=-1.5,
                        vmax=1.5,
                        add_colorbar=False,
                        extend="both",
                        cmap=cm.get_cmap("RdBu_r"),
                    )
                else:
                    wind_speed_components_and_luh_dict[experiment][field][
                        "gothr+gsecd"
                    ].plot(
                        ax=ax[row, 3],
                        vmin=-1,
                        vmax=1,
                        add_colorbar=True,
                        cbar_ax=cbar2_ax,
                        extend="both",
                        cmap=cm.get_cmap("RdBu"),
                        cbar_kwargs={
                            "orientation": "horizontal",
                            "label": "Difference in primary \n plus secondary land",
                        },
                    )

    for tmp_ax in ax.flatten():
        tmp_ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
        tmp_ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
    ax[0, 0].set_title(r"Full change $\Delta s$")
    ax[0, 1].set_title(r"Dynamical change $\Delta_{dyn} s$")
    ax[0, 2].set_title(r"Residual change $\Delta_{res} s$")
    ax[0, 3].set_title(r"Land use change")

    add_letters(ax, x=-0.05, y=0.98)
    plt.subplots_adjust(0.03, 0.14, 0.99, 0.96, 0.1, 0.1)
    plt.savefig(f"{path_to_plots}/attribution_maps.jpeg", dpi=400)
    plt.close("all")


def plot_onshore_contribution_histograms(
    wind_speed_components_and_luh_dict, path_to_data, path_to_plots
):
    # compare onshore pdfs
    landmask = open_landmask(path_to_data)
    colors = ["Orange", "Olive"]
    f, ax = plt.subplots(nrows=4, sharex=True, figsize=(4, 8))
    label_dict = {
        "Full - Dyn.": r"Residual change $\Delta_{res} s$",
        "Dyn. Diff": r"Dynamical change $\Delta_{dyn} s$",
    }
    for row, experiment in enumerate(wind_speed_components_and_luh_dict.keys()):
        total_change = 0
        for i, var in enumerate(["Full - Dyn.", "Dyn. Diff"]):
            tmp_da = wind_speed_components_and_luh_dict[experiment][var]["sfcWind"]
            tmp_da = tmp_da.where(np.isfinite(landmask)).values
            ax[row].hist(
                tmp_da[np.isfinite(tmp_da)],
                label=label_dict[var],
                density=True,
                bins=100,
                alpha=0.7,
                color=colors[i],
            )
            ax[row].axvline(
                tmp_da[np.isfinite(tmp_da)].mean(), ls="--", color=colors[i]
            )
            total_change += tmp_da[np.isfinite(tmp_da)].mean()
        ax[row].axvline(
            total_change, ls="--", color="black", label=r"Full change $\Delta s$"
        )
        ax[row].set_ylabel(experiment + " PDF")
        ax[row].set_xlim(xmin=-0.8, xmax=0.8)
    ax[0].legend(bbox_to_anchor=(1.03, 1.6))
    ax[3].set_xlabel("Wind speed change [m/s]")
    add_letters(ax)
    plt.subplots_adjust(0.12, 0.06, 0.95, 0.89)
    plt.savefig(f"{path_to_plots}/contribution_histograms.jpeg", dpi=400)
    plt.close("all")


def plot_luh_vs_wind_speed_scatter(
    wind_speed_components_and_luh_dict, path_to_data, path_to_plots
):
    """
    scatter plots
    """

    # TODO: better name than 'ref'
    ref = reference_ensemble_mean_wind_speed(path_to_data)
    luh1_variable = "gothr+gsecd"
    for wind_type in ["abs", "rel"]:
        f, ax = plt.subplots(nrows=4, sharex=True, sharey=True, figsize=(4, 8))
        for row, experiment in enumerate(wind_speed_components_and_luh_dict.keys()):
            ds_wind = wind_speed_components_and_luh_dict[experiment]["Full - Dyn."][
                "sfcWind"
            ].copy()
            if wind_type == "rel":
                ds_wind /= ref["sfcWind"]

            ds_lu = wind_speed_components_and_luh_dict[experiment]["LUH Diff"][
                luh1_variable
            ].copy()
            # land use lons cover -180 ... 180, tranform to 0 ... 360
            ds_lu = ds_lu.assign_coords(lon=((ds_lu.lon + 360) % 360)).sortby("lon")
            # aggregate land use data to similar resolution (2x2 degrees compared to 1.875x1.875)
            ds_lu_agg = (
                ds_lu.coarsen(lat=4, boundary="trim")
                .mean()
                .coarsen(lon=4, boundary="trim")
                .mean()
            )

            regridder = xe.Regridder(ds_lu_agg, ds_wind, "bilinear", periodic=True)
            ds_lu_int = regridder(ds_lu_agg)
            y = ds_lu_int.where(np.abs(ds_lu_int) > 0.01)
            x = ds_wind.where(np.abs(ds_lu_int) > 0.01)

            # drop nans
            x, y = x.values, y.values
            x = x[np.isfinite(x)]
            y = y[np.isfinite(y)]
            label = (
                f"R = {np.round(pearsonr(x, y)[0], 2)}\n"
                fr"$\rho$ = {np.round(spearmanr(x, y)[0], 2)}"
            )
            ax[row].scatter(x, y, s=5, alpha=0.3, label=label)
            ax[row].set_ylabel("Land use change")
            ax[row].axhline(0, ls="--", color="Orange", alpha=0.7)
            ax[row].axvline(0, ls="--", color="Orange", alpha=0.7)

            # add occurence percentage of all four quadrants
            q1 = int(np.round(y[(x < 0) & (y > 0)].size / y.size * 100))
            q2 = int(np.round(y[(x > 0) & (y > 0)].size / y.size * 100))
            q3 = int(np.round(y[(x > 0) & (y < 0)].size / y.size * 100))
            q4 = int(np.round(y[(x < 0) & (y < 0)].size / y.size * 100))

            ax[row].text(
                -1, 0.2, f"{q1}%", color="Darkorange", ha="center", va="center"
            )
            ax[row].text(
                2.5, 0.2, f"{q2}%", color="Darkorange", ha="center", va="center"
            )
            ax[row].text(
                2.5, -0.2, f"{q3}%", color="Darkorange", ha="center", va="center"
            )
            ax[row].text(
                -1, -0.2, f"{q4}%", color="Darkorange", ha="center", va="center"
            )
            ax[row].legend(loc="upper right", markerscale=0)
            ax[row].set_ylim(ymax=0.99)

            # add conditional probability
            ax2 = ax[row].twinx()
            ys = [
                conditional_prob(x, y, x_thr=x_thr) * 100
                for x_thr in np.arange(0, 2, 0.1)
            ]
            ax2.plot(np.arange(0, 2, 0.1), ys, color="Purple")
            ax2.set_ylabel("cond. prob. luc <0 [%]", color="Purple")
            ax2.set_ylim(ymin=95, ymax=107)
            ax2.set_yticks([96, 98, 100])
            [t.set_color("Purple") for t in ax2.yaxis.get_ticklines()]
            [t.set_color("Purple") for t in ax2.yaxis.get_ticklabels()]
        # ax[3].set_xlim(xmin=-0.5, xmax=2)
        if wind_type == "abs":
            ax[3].set_xlabel(r"Residual wind speed change $\Delta_{res}s$ [m/s]")
        else:
            ax[3].set_xlabel("Relative residual wind speed change [1]")
        add_letters(ax)
        plt.tight_layout()
        fig_path = f"{path_to_plots}/scatter/scatter_{luh1_variable}_{wind_type}.jpeg"
        plt.savefig(fig_path, dpi=300)
