import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=ImportWarning)
    import cartopy.crs as ccrs
    import cartopy


def forest_perarea(ds):
    """
    Fraction of grid box covered by forest.
    Calculated following the suggestions in the FAQ (see below), then multiplied with cell_area
    because relative evolution appears more relevant here.

    From LUH FAQ (https://luh.umd.edu/faq.shtml, Jan 27th 2021):
    "You can compute an estimate of the forest area in the Land-Use Harmonization products.
    For the original LUH products, first you will need to download the forest/non-forest map
    (fnf_map.txt) and the grid-cell area map (cellarea_halfdeg.txt) in addition to the LUH data
    packages.
    The forest area in an individual grid-cell is given by:
        (gothr + gsecd)*fnf*cellarea
    (where gothr and gsecd are the fractions of the grid-cell occupied by primary and secondary land respectively).
    :param ds
    :return:
    """
    return (ds["gothr"] + ds["gsecd"]) * ds["fnf"]


def open_LUH(filename, year=None, name=None):
    da = (
        xr.open_rasterio(filename)
        .rename({"x": "lat", "y": "lon"})
        .drop("band")
        .squeeze()
    )
    if year:
        da = da.assign_coords({"time": pd.Timestamp(str(year) + "-01-01")})
    if name:
        da.name = name
    return da


def open_LUH_period(path, ts, te):
    """
    open time varying primary and secondary vegeation maps and combine them in a single xarray Dataset
    :param path: path to data
    :param ts: start year of respective experiment, e.g. 1850 for historical
    :param te: end year for respective experiment, e.g. 2000 for historical
    :return:
    """
    da_gothr = xr.concat(
        [
            open_LUH(path + "updated_states/gothr." + str(year) + ".txt", year, "gothr")
            for year in range(ts, te)
        ],
        dim="time",
    )
    da_gsecd = xr.concat(
        [
            open_LUH(path + "updated_states/gsecd." + str(year) + ".txt", year, "gsecd")
            for year in range(ts, te)
        ],
        dim="time",
    )
    ds = xr.merge([da_fnf, da_gothr, da_gsecd])
    ds["forest"] = forest_perarea(ds)
    ds["gothr+gsecd"] = ds["gothr"] + ds["gsecd"]
    return ds


# open constant forest/non-forest map
da_fnf = open_LUH("../data/LUHa.v1/fnf_map.txt", name="fnf")
# open time varying fields
ds = open_LUH_period("../data/LUHa.v1/", 1850, 2000)

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
for j, start_year in enumerate([1850, 1890, 1940, 1990]):
    t_slice = slice(str(start_year), str(start_year + 10))
    for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
        ds[var].sel({"time": t_slice}).mean(dim="time").plot(
            ax=ax[j, i], vmin=0.01, vmax=1.01, levels=11, extend="neither"
        )
        ax[j, i].set_title(str(start_year) + " to " + str(start_year + 10))
        ax[j, i].add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
        ax[j, i].add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)

plt.tight_layout()
plt.savefig("../plots/LUH1/LUH_absolute.jpeg", dpi=300)


"""
plot ref and changes to ref
"""
ref = ds.sel({"time": slice("1850", "1860")}).mean(dim="time")
f, ax = plt.subplots(
    ncols=4,
    nrows=4,
    figsize=(12, 10),
    sharex=True,
    sharey=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)

for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
    ref[var].plot(ax=ax[0, i], vmin=0.01, vmax=1.01, levels=11, extend="neither")
    ax[0, i].set_title(str(start_year) + " to " + str(start_year + 10))

for j, start_year in enumerate([1890, 1940, 1990]):
    t_slice = slice(str(start_year), str(start_year + 10))
    diff = ds.sel({"time": t_slice}).mean(dim="time") - ref
    for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
        diff[var].plot(
            ax=ax[j + 1, i],
            vmin=-1,
            vmax=1,
            # levels=11,
            extend="neither",
            cmap=cm.get_cmap("RdBu_r"),
        )
        ax[j + 1, i].set_title(
            "Diff " + str(start_year) + " to " + str(start_year + 10)
        )

for tmp_ax in ax.flatten():
    tmp_ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
    tmp_ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
plt.tight_layout()
plt.savefig("../plots/LUH1/LUH_change.jpeg", dpi=300)


"""
plot ref and changes to ref in the future
"""

f, ax = plt.subplots(
    ncols=4,
    nrows=3,
    figsize=(12, 8),
    sharex=True,
    sharey=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)
cbar_ax = f.add_axes([0.2, 0.08, 0.6, 0.02])

ref_end = ds.sel({"time": slice("1990", "2000")}).mean(dim="time")  # 2nd reference scenario end of 20th century

for i_ref, reference_data in enumerate([ref, ref_end]):
    for j, experiment in enumerate(["IMAGE_rcp26", "MiniCAM_rcp45", "MESSAGE_rcp85"]):
        print(experiment)
        ds_future = open_LUH_period("../data/LUHa.v1/" + experiment + "/", 2090, 2100)
        diff = ds_future.sel({"time": slice("2090", "2100")}).mean(dim="time") - reference_data
        for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
            if i + j != 0:
                diff[var].plot(
                    ax=ax[j, i],
                    vmin=-1,
                    vmax=1,
                    add_colorbar=False,
                    extend="neither",
                    cmap=cm.get_cmap("RdBu_r"),
                )
            else:
                diff[var].plot(
                    ax=ax[j, i],
                    vmin=-1,
                    vmax=1,
                    add_colorbar=True,
                    cbar_ax=cbar_ax,
                    extend="neither",
                    cmap=cm.get_cmap("RdBu_r"),
                    cbar_kwargs={"orientation": "horizontal", "label": ""},
                )
            ax[j, i].set_title(var)
        plt.text(
            0,
            1.4,
            experiment,
            fontdict={"fontsize": 15},
            horizontalalignment="left",
            verticalalignment="center",
            transform=ax[j, 0].transAxes,
        )

        # ax[j, 0].set_yaxis(experiment)

    for tmp_ax in ax.flatten():
        tmp_ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), lw=0.2)
        tmp_ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), lw=0.2)
    plt.subplots_adjust(0.01, 0.1, 0.99, 0.96)
    if i_ref == 0:
        plt.savefig("../plots/LUH1/LUH_change_future.jpeg", dpi=300)
    else:
        plt.savefig("../plots/LUH1/LUH_change_future_ref2000.jpeg", dpi=300)
