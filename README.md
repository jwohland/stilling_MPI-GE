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


All figures can be created at once by running the function `make_all_plots.py` from the command line, assuming you are in the Anaconda environment given in this repository:

```bash
python make_all_plots.py [data_path] [plots_path] [--cache_path]
```

Where directories for input data and output plots are user-defined (default `data` and `plots`).
`--cache_path` can be also defined to store some intermediate processed data.

## Other Files

`requirements.yaml` can be used to create an Anaconda environment with all the necessary packages to run the scripts in this repository.

`make_land_mask.py` computes a land mask on the MPI grid using runoff data and excluding Antartica and Greenland.

`Lmon_Amon_diff.py` verify that wind speed data provided in the modeling realm Atmosphere (Amon) and Land (Lmon) contain the same values.

## Input data

Input data is from four main sources:

1. [The Max Planck Institute for Meteorology (MPI-M) Grand Ensemble (MPI-GE)](https://esgf-data.dkrz.de/search/mpi-ge/). [Ref: DOI 10.1029/2019MS001639](https://doi.org/10.1029/2019MS001639)
2. [The Land Use Harmonization (LUH) Project](https://luh.umd.edu/data.shtml#LUH1_Data). [Ref: DOI 10.1007/s10584-011-0153-2](https://doi.org/10.1007/s10584-011-0153-2)
3. [The Coupled Model Intercomparison Project Phase 6 (CMIP6)](https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip6). [Ref: DOI 10.5194/gmd-9-1937-2016](https://doi.org/10.5194/gmd-9-1937-2016)
4. Zeng et al. (2019). [Ref: DOI 10.1038/s41558-019-0622-6](https://doi.org/10.1038/s41558-019-0622-6)

Instructions on accessing the exact datasets required from these sources are given in the README files of the respective subdirectories of the `data` directory. All directories must be filled with data to reproduce this study's analysis.
