"""
Compute forested area from LUH1 following the the LUH FAQ (https://luh.umd.edu/faq.shtml, Jan 27th 2021):

You can compute an estimate of the forest area in the Land-Use Harmonization products. 
For the original LUH products, first you will need to download the forest/non-forest map 
(fnf_map.txt) and the grid-cell area map (cellarea_halfdeg.txt) in addition to the LUH data 
packages. 
The forest area in an individual grid-cell is given by: 
    (gothr + gsecd)*fnf*cellarea 
(where gothr and gsecd are the fractions of the grid-cell occupied by primary and secondary land respectively). 
To compute the global forest area you would sum this quantity over all grid-cells globally.
"""

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm


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


# open constant forest/non-forest map
da_fnf = open_LUH("fnf_map.txt", name="fnf")
# open time varying primary and secondary vegetation maps
da_gothr = xr.concat(
    [
        open_LUH("updated_states/gothr." + str(year) + ".txt", year, "gothr")
        for year in range(1850, 2000)
    ],
    dim="time",
)
da_gsecd = xr.concat(
    [
        open_LUH("updated_states/gsecd." + str(year) + ".txt", year, "gsecd")
        for year in range(1850, 2000)
    ],
    dim="time",
)
ds = xr.merge([da_fnf, da_gothr, da_gsecd])
ds["forest"] = forest_perarea(ds)
ds["gothr+gsecd"] = ds["gothr"] + ds["gsecd"]


"""
plot absolute values
"""
f, ax = plt.subplots(ncols=4, nrows=4, figsize=(12, 10))
for j, start_year in enumerate([1850, 1890, 1940, 1990]):
    t_slice = slice(str(start_year), str(start_year + 10))
    for i, var in enumerate(["gothr", "gsecd", "gothr+gsecd", "forest"]):
        ds[var].sel({"time": t_slice}).mean(dim="time").plot(
            ax=ax[j, i], vmin=0.01, vmax=1.01, levels=11, extend="neither"
        )
        ax[j, i].set_title(str(start_year) + " to " + str(start_year + 10))


plt.tight_layout()
plt.savefig("LUH_absolute.jpeg", dpi=300)


"""
plot ref and changes to ref
"""
ref = ds.sel({"time": slice("1850", "1860")}).mean(dim="time")
f, ax = plt.subplots(ncols=4, nrows=4, figsize=(12, 10))

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
            levels=11,
            extend="neither",
            cmap=cm.get_cmap("RdBu_r"),
        )
        ax[j + 1, i].set_title(
            "Diff " + str(start_year) + " to " + str(start_year + 10)
        )

plt.tight_layout()
plt.savefig("LUH_change.jpeg", dpi=300)
