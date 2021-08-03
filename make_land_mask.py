import xarray as xr
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import numpy as np


def in_greenland(da):
    poly_greenland = Polygon(
        [
            (-46.2, 56.2),
            (-75.7, 78.1),
            (49.5, 84.2),
            (-0.8, 84.2),
            (-1.0, 76.0),
            (-46.2, 56.2),
        ]
    )  # polygon containing greenland, notation is lon, lat
    lon_shifted = da.lon
    if lon_shifted > 180:
        lon_shifted -= 360
    if poly_greenland.contains(Point(lon_shifted, da.lat)):
        return True


def create_landmask(path_to_data):
    runoff = xr.open_dataset(
        f"{path_to_data}/runoff/mrro_Lmon_MPI-ESM_historical_r001i1850p3_timestep1.nc"
    )
    landmask = runoff["mrro"].drop("time").squeeze()

    # step1: exclude Antarctica
    landmask = landmask.where(landmask.lat > -60)

    # step2: exclude Greenland
    for i_lat in range(96):
        for i_lon in range(192):
            # if location is not already masked for another reason
            #    --> mask it not if in Greenland
            if np.isfinite(landmask[i_lat, i_lon]):
                if in_greenland(landmask[i_lat, i_lon]):
                    landmask[i_lat, i_lon] = np.nan
    landmask /= landmask

    landmask.name = "landmask"
    landmask.attrs = {
        "long_name": "Land mask based on non-zero runoff locations excluding ice sheets"
    }

    landmask.to_netcdf(f"{path_to_data}/runoff/landmask.nc")
