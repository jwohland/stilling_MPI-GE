"""
preprocessing
- for 3D fluxes (FDS, FDSC), keep the surface level only (1000 hPa)
- for 3D winds (U, V), keep lowest levels only (>960 hPa, i.e.altitudes < 475m for default values, 2 levels)
- cut out Europe (30 to 75 lons -15 to 35 with 1 degree buffer at each end)
- only keep wind speeds S = (U**2 + V**2)**(1/2) instead of components
"""

import xarray as xr
import glob

base_path = "/cluster/work/apatt/wojan/CESM2/"  # todo update


def select_levels(ds):
    levs = ds["U"].lev
    levs = levs[levs > 960]  # wind up to 960 hPa level
    ilev = 1000  # flux only lowest level
    return ds.sel({"lev": levs, "ilev": ilev})


def select_Europe(ds):
    lat_min, lat_max = 29, 76
    lon_min, lon_max = -16 + 360, 36
    lats = slice(lat_min, lat_max)
    # lon handling done manually as 360 degree symmetry not captured
    lons = list(ds.lon[ds.lon > lon_min].values)
    lons.extend(list(ds.lon[ds.lon < lon_max].values))
    return ds.sel({"lat": lats, "lon": lons})


def compute_windspeed(ds):
    ds["S"] = (ds["U"] ** 2 + ds["V"] ** 2) ** (1.0 / 2)
    ds["S"].attrs = {
        "units": "m/s",
        "long_name": "wind speed",
    }
    return ds


def preprocess(ds):
    ds = select_levels(ds)
    ds = select_Europe(ds)
    ds = compute_windspeed(ds)
    ds.compute()
    return ds.drop(["U", "V"])


# h7 is average, h8 is instantaneous
for stream in ["h7", "h8"]:
    filename = glob.glob(base_path + "raw/*" + stream + "*")[0].split("/")[-1]
    ds = xr.open_dataset(base_path + "raw/" + filename)
    ds = preprocess(ds)
    ds.to_netcdf(base_path + "preprocessed/" + filename)
