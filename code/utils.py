import glob
import string
import os

import xarray as xr
import pandas as pd

"""
European bounding box
"""

LATMIN, LATMAX = 37.5, 60
LONMIN, LONMAX = -10, 25

"""
functions related to loading and slicing the data
"""


def selbox(ds):
    """Get data averaged over all gridboxes in the European region"""
    lats, lons = slice(LATMIN, LATMAX), slice(LONMIN, LONMAX)
    return ds.sel({"lat": lats, "lon": lons}).mean(dim=["lat", "lon"])


def annual_mean(ds):
    """Get annually average data"""
    return ds.resample({"time": "1Y"}).mean()


def sel_time(ds, tslice, mean_time=True):
    """
    Get a time slice `tslice` from xarray data object `ds`.
    If `mean_time` is True, result will be averaged over the time dimension after slicing
    """
    ds_sliced = ds.sel({"time": tslice})
    if mean_time is True:
        return ds_sliced.mean(dim="time")
    else:
        return ds_sliced


def mean_sliced_annual_mean(ds, tslice, mean_time=True):
    return sel_time(annual_mean(ds), tslice, mean_time)


# MPI-GE
def open_datasets(filelist):
    """Open and concatenate data from all MPI-GE ensemble members"""
    ds = [annual_mean(selbox(xr.open_dataset(x, use_cftime=True))) for x in filelist]
    ds = xr.concat(
        ds,
        dim=pd.Index(name="ensemble_member", data=[x.split("_")[-2] for x in filelist]),
    )
    return ds


def ensemble_mean_wind_speed(path_to_data, experiment):
    """
    Get reference ensemble mean wind speed for a given experiment
    :param experiment: string name for the simulated experiment, one of [historical, rcp26, rcp45, rcp85]
    """
    ensmean_path = glob.glob(f"{path_to_data}/{experiment}/ensmean/*.nc")
    assert len(ensmean_path) == 1

    return xr.open_dataset(ensmean_path[0])


def reference_ensemble_mean_wind_speed(path_to_data):
    """Get reference ensemble mean wind speed for the period 1850-1859 (incl.)"""
    return mean_sliced_annual_mean(
        ensemble_mean_wind_speed(path_to_data, "historical"), slice("1850", "1859")
    )


def open_picontrol(path_to_data, spatial_averaging=True):
    """Open and concatenate all pre-industrial control simulation years"""

    def _open_ds(ds_path):
        # use_cftime needed after 2200. Otherwise SerializationWarning is raised
        if spatial_averaging:
            return selbox(xr.open_dataset(ds_path, use_cftime=True))
        else:
            return xr.open_dataset(ds_path, use_cftime=True)

    ds_list = [
        annual_mean(_open_ds(x))
        for x in sorted(glob.glob(f"{path_to_data}/pi-control/*.nc"))
    ]
    return xr.concat(ds_list, dim="time")


def open_LUH(filename, year=None, name=None):
    """Open LUH1 (land-use change) raster and add data on year and name (if given)"""
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


def open_LUH_period(path_to_LUH1, ts, te, path_to_experiment="."):
    """
    open time varying primary and secondary vegeation maps and combine them in a single xarray Dataset
    :param path: path to data
    :param ts: start year of respective experiment, e.g. 1850 for historical
    :param te: end year for respective experiment, e.g. 2000 for historical
    :return:
    """
    path_to_data = f"{path_to_LUH1}/{path_to_experiment}/updated_states/{{indicator}}.{{year}}.txt"

    ds = xr.merge([
        xr.concat([
            open_LUH(path_to_data.format(indicator=indicator, year=year), year, indicator)
            for year in range(ts, te)
        ], dim="time")
        for indicator in ["gothr", "gsecd"]
    ])

    ds["forest"] = forest_perarea(ds, path_to_LUH1)
    ds["gothr+gsecd"] = ds["gothr"] + ds["gsecd"]
    return ds


def forest_perarea(ds, path_to_LUH1):
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
    da_fnf = open_LUH(f"{path_to_LUH1}/fnf_map.txt", name="fnf")
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
