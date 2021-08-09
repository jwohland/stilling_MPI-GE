"""
rsds and sfcWind are not consistently available either as Amon or Lmon fields.

Check here whether sfcWinds during historical period identical in Amon and Lmon

"""

import xarray as xr
import glob

ds_Lmon = xr.open_dataset(glob.glob("../data/historical/ensmean/*.nc")[0])
ds_Amon = xr.open_dataset(glob.glob("../data/historical/Amon/ensmean/*.nc")[0])

# Compare values because time steps offset by about half a month
diff = ds_Amon["sfcWind"].values - ds_Lmon["sfcWind"].values
assert (diff == 0).all()
