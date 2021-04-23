import xarray as xr
import pandas as pd
import glob
import matplotlib.pyplot as plt

"""
European bounding box
"""

LATMIN, LATMAX = 37.5, 60
LONMIN, LONMAX = -10, 25

"""
functions related to loading and slicing the data
"""


def selbox(ds):
    lats, lons = slice(LATMIN, LATMAX), slice(LONMIN, LONMAX)
    return ds.sel({"lat": lats, "lon": lons}).mean(dim=["lat", "lon"])


def ann_mean(ds):
    return ds.resample({"time": "1Y"}).mean()


def sel_time(ds, tslice):
    return ds.sel({"time": tslice})


# MPI-GE
def open_datasets(filelist):
    ds = [ann_mean(selbox(xr.open_dataset(x, use_cftime=True))) for x in filelist]
    ds = xr.concat(
        ds,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in filelist]),
    )
    return ds


def open_picontrol():
    path = "../data/pi-control/"
    ds_list = [
        ann_mean(selbox(xr.open_dataset(x, use_cftime=True)))
        for x in sorted(glob.glob(path + "*.nc"))
    ]  # use_cftime needed after 2200. Otherwise SerializationWarning is raised
    return xr.concat(ds_list, dim="time")

# LUH
def open_LUH(filename, year=None, name=None):
    da = (
        xr.open_rasterio(filename)
        .rename({"x": "lon", "y": "lat"})
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
    ds = xr.merge([da_gothr, da_gsecd])
    ds["forest"] = forest_perarea(ds)
    ds["gothr+gsecd"] = ds["gothr"] + ds["gsecd"]
    return ds


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
    # open constant forest/non-forest map
    da_fnf = open_LUH("../data/LUHa.v1/fnf_map.txt", name="fnf")
    return (ds["gothr"] + ds["gsecd"]) * da_fnf



def add_letters(ax, x=-0.08, y=1.02, fs=10, letter_offset=0):
    """
    adds bold letters a,b,c,... to the upper left corner of subplots
    :param ax: axis
    :param x: x location of text
    :param y: ylocation of text
    :param fs: fontsize
    :return:
    """
    import string
    letters = list(string.ascii_lowercase)
    try:
        ax.flat
        for il, tmp_ax in enumerate(ax.flat):
            tmp_ax.text(
                x,
                y,
                letters[il + letter_offset],
                weight="bold",
                horizontalalignment="center",
                verticalalignment="center",
                transform=tmp_ax.transAxes,
                fontsize=fs,
            )
    except AttributeError:
        ax.text(
            x,
            y,
            letters[letter_offset],
            weight="bold",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=fs,
        )