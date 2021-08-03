# MPI-GE wind stilling

Code underlying analysis performed in `Wohland, J., Folini, D., Pickering, B., 2021. Wind speed stilling and its recovery due to internal climate variability. Earth System Dynamics Discussions 1–27`. [Ref: DOI 10.5194/esd-2021-29](https://doi.org/10.5194/esd-2021-29).
## Figure overview

| Figure | Filename | Creating python script |
|---|---|---|
Fig. 1 | attribution_maps.jpeg | attribution.py:plot_attribution_maps
Fig. 2 | contribution_histograms.jpeg | attribution.py:plot_onshore_contribution_histograms
Fig. 3 | scatter_gothr+gsecd_abs.jpeg | attribution.py:plot_luh_vs_wind_speed_scatter
Fig. 4 | LUH1/LUH_change_future_ref2000.jpeg | LUH_plots.py:plot_future_LUH_change
Fig. 5 | global_windspeeds.jpeg | trend_maps.py:plot_global_windspeeds
Fig. 6 | timeseries_picontrol_Europe.jpeg | trends.py:plot_full_timeseries_with_trend_marks
Fig. 7a | picontrol_wind_trends_Europe_5_20y.jpeg | trends.py:plot_trend_histograms
Fig. 7b | CMIP6/Ensmean_picontrol_wind_trends_Europe_5.jpeg | trends.py:plot_pi_control_cmip6_trend_histograms
Fig. 8a | historical_wind_trends_Europe_5_all.jpeg | trends.py:plot_experiment_trend_histograms
Fig. 8b | rcp26_wind_trends_Europe_5_all.jpeg | trends.py:plot_experiment_trend_histograms
Fig. 8c | rcp45_wind_trends_Europe_5_all.jpeg | trends.py:plot_experiment_trend_histograms
Fig. 8d | rcp85_wind_trends_Europe_5_all.jpeg | trends.py:plot_experiment_trend_histograms
Appendix Fig. A1| box_timeseries_Europe.jpeg | ensmean_timeseries.py:plot_ensemble_members_timeseries
Appendix Fig. A2| picontrol_HadISD_wind_trends_Europe_5_20y.jpeg | trends.py:plot_trend_histograms
Appendix Fig. A3a | picontrol_wind_trends_Europe_5_15y.jpeg | trends.py:plot_trend_histograms
Appendix Fig. A3b | picontrol_wind_trends_Europe_5_25y.jpeg | trends.py:plot_trend_histograms
Appendix Fig. A4a | picontrol_wind_trends_Europe_10_20y.jpeg | trends.py:plot_trend_histograms
Appendix Fig. A4b | picontrol_wind_trends_Europe_15_20y.jpeg | trends.py:plot_trend_histograms
Appendix Figs. 5-7 | CMIP6/{Model_name}_picontrol_wind_trends_Europe_5.jpeg | trends.py:plot_pi_control_cmip6_trend_histograms


All figures can be created at once by running the function `make_all_plots.py` from the command line:

```bash
python make_all_plots.py [data_path] [plots_path] [--cache_path]
```

Where directories for input data and output plots are user-defined (default `data` and `plots`).
`--cache_path` can be also defined to store some intermediate processed data.

## Other Files

make_land_mask.py computes a land mask on the MPI grid using runoff data and excluding Antartica and Greenland. 

Lmon_Amon_diff.py verify that wind speed data provided in the modeling realm Atmosphere (Amon) and Land (Lmon) contain the same values. 

## Input data

Input data is from the Max Planck Institute for Meteorology (MPI-M) Grand Ensemble (MPI-GE) that can be downloaded at https://esgf-data.dkrz.de/search/mpi-ge/ and from the Land Use Harmonization (LUH) Project that is accessable via https://luh.umd.edu/data.shtml#LUH1_Data.

Maher, N., Milinski, S., Suarez‐Gutierrez, L., Botzet, M., Dobrynin, M., Kornblueh, L., Kröger, J., Takano, Y., Ghosh, R., Hedemann, C., Li, C., Li, H., Manzini, E., Notz, D., Putrasahan, D., Boysen, L., Claussen, M., Ilyina, T., Olonscheck, D., Raddatz, T., Stevens, B., Marotzke, J., 2019. The Max Planck Institute Grand Ensemble: Enabling the Exploration of Climate System Variability. J. Adv. Model. Earth Syst. 11, 2050–2069. https://doi.org/10.1029/2019MS001639

Hurtt, G.C., Chini, L.P., Frolking, S., Betts, R.A., Feddema, J., Fischer, G., Fisk, J.P., Hibbard, K., Houghton, R.A., Janetos, A., Jones, C.D., Kindermann, G., Kinoshita, T., Klein Goldewijk, K., Riahi, K., Shevliakova, E., Smith, S., Stehfest, E., Thomson, A., Thornton, P., van Vuuren, D.P., Wang, Y.P., 2011. Harmonization of land-use scenarios for the period 1500–2100: 600 years of global gridded annual land-use transitions, wood harvest, and resulting secondary lands. Climatic Change 109, 117–161. https://doi.org/10.1007/s10584-011-0153-2
