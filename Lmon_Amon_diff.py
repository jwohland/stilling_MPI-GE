"""
rsds and sfcWind are not consistently available either as Amon or Lmon fields.

Check here whether sfcWinds during historical period identical in Amon and Lmon

"""

from utils import ensemble_mean_wind_speed


def test_historical_ensemble_mean_wind_speeds(path_to_data):
    ds_Lmon = ensemble_mean_wind_speed(path_to_data, "historical")
    ds_Amon = ensemble_mean_wind_speed(path_to_data, "historical/Amon")

    # Compare values because time steps offset by about half a month
    diff = ds_Amon["sfcWind"].values - ds_Lmon["sfcWind"].values
    assert (diff == 0).all()
