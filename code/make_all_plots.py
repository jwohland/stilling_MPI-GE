import argparse
import os

from attribution import (
    experiment_wind_speed_components_and_luh,
    plot_attribution_maps,
    plot_onshore_contribution_histograms,
    plot_luh_vs_wind_speed_scatter,
)
from LUH_plots import plot_future_LUH_change
from trend_maps import plot_global_windspeeds
from trends import (
    plot_full_timeseries_with_trend_marks,
    test_calc_frac_partoftrend,
    plot_trend_histograms,
    plot_pi_control_cmip6_trend_histograms,
    plot_experiment_trend_histograms,
)
from ensmean_timeseries import plot_ensemble_members_timeseries


def make_abs_path(my_path):
    this_dir = os.path.dirname(__file__)
    if not os.path.isabs(my_path):
        return os.path.normpath(os.path.join(this_dir, my_path))
    else:
        return my_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data_path",
        default="../data",
        help="Path to data directory. If not an absolute path, will be assumed as relative to the directory in which `make_all_plots.py` is found.",
        type=make_abs_path,
    )
    parser.add_argument(
        "plots_path",
        default="../plots",
        help="Path to plot output directory. If not an absolute path, will be assumed as relative to the directory in which `make_all_plots.py` is found.",
        type=make_abs_path,
    )
    parser.add_argument(
        "--cache_path",
        default=None,
        help="Path to intermediate data caching directory. If not an absolute path, will be assumed as relative to the directory in which `make_all_plots.py` is found.",
        type=make_abs_path,
    )

    args = parser.parse_args()
    path_to_data = args.data_path
    path_to_plots = args.plots_path
    path_to_cache = args.cache_path

    # This data dictionary is used for Figs. 1-3
    data_dict = experiment_wind_speed_components_and_luh(path_to_data, path_to_cache)

    print("Producing Fig. 1 | attribution_maps.jpeg")
    plot_attribution_maps(data_dict, path_to_plots)

    print("Producing Fig. 2 | contribution_histograms.jpeg")
    plot_onshore_contribution_histograms(data_dict, path_to_data, path_to_plots)

    print("Producing Fig. 3 | scatter_gothr+gsecd_abs.jpeg")
    plot_luh_vs_wind_speed_scatter(data_dict, path_to_data, path_to_plots)

    print("Producing Fig. 4 | LUH_change_future_ref2000.jpeg")
    plot_future_LUH_change(path_to_data, path_to_plots)

    print("Producing Fig. 5 | global_windspeeds.jpeg ()")
    plot_global_windspeeds(path_to_data, path_to_plots)

    print("Producing Fig. 6 | timeseries_picontrol_Europe.jpeg")
    plot_full_timeseries_with_trend_marks(path_to_data, path_to_plots)

    # Test to ensure the calculation of timesteps in a 20y trend is working as expected (used in next three plots)
    test_calc_frac_partoftrend()

    print("Producing Fig. 7a & Appendix Figs. A2, A3a-b, A4a-b | picontrol_wind_trends_Europe_5_20y.jpeg & picontrol_{HadISD|None}_wind_trends_Europe_{p_threshold}_{trend_length}.jpeg")
    plot_trend_histograms(path_to_data, path_to_plots)

    print("Producing Fig. 7b & Appendix Figs. 5-7 | CMIP6/Ensmean_picontrol_wind_trends_Europe_5.jpeg & CMIP6/{Model_name}_picontrol_wind_trends_Europe_5.jpeg")
    plot_pi_control_cmip6_trend_histograms(path_to_data, path_to_plots)

    print("Producing Figs. 8a-d | {experiment}_wind_trends_Europe_5_all.jpeg")
    plot_experiment_trend_histograms(path_to_data, path_to_plots)

    print("Producing Appendix Fig. A1| box_timeseries_Europe.jpeg")
    plot_ensemble_members_timeseries(path_to_data, path_to_plots)
