# MPI-GE wind stilling

Code underlying analysis performed in

Wohland, J., Folini, D., Pickering, B., 2021. Wind speed stilling and its recovery due to internal climate variability (preprint). Under review for Earth System Dynamics. https://doi.org/10.5194/esd-2021-29

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

`requirements.yaml` can be used to create an Anaconda environment with all the necessary packages to run the scripts in this repository (including the command line tool `cdo`).
## Input data
Executing scripts from this repository requires (a) downloading of the required input data and (b) preprocessing the data, and (c) updating the file paths.

### a) Downloads
To download the input, please follow the structure and recommendations given in the `data` subdirectories.
The MPI-GE data can be downloaded using the provided `wget` scripts after you have registered at ESGF (see https://esgf-data.dkrz.de/user/add/?next=http://esgf-data.dkrz.de/projects/esgf-dkrz/ to create a new profile).
The other downloads are explained in `README.md` files in the respective subfolders.

Input data is from 4 different sources:
1. [The Max Planck Institute for Meteorology (MPI-M) Grand Ensemble (MPI-GE)](https://esgf-data.dkrz.de/search/mpi-ge/). [Ref: DOI 10.1029/2019MS001639](https://doi.org/10.1029/2019MS001639)
2. [The Land Use Harmonization (LUH) Project](https://luh.umd.edu/data.shtml#LUH1_Data). [Ref: DOI 10.1007/s10584-011-0153-2](https://doi.org/10.1007/s10584-011-0153-2)
3. [The Coupled Model Intercomparison Project Phase 6 (CMIP6)](https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip6). [Ref: DOI 10.5194/gmd-9-1937-2016](https://doi.org/10.5194/gmd-9-1937-2016)
4. Zeng et al. (2019). [Ref: DOI 10.1038/s41558-019-0622-6](https://doi.org/10.1038/s41558-019-0622-6)

### b) Preprocessing
To speed up Python calculations, the MPI-GE ensemble means are pre-calculated using CDO (see `code/preprocess/preprocess_MPI-GE.sh`).
They are stored in a subfolder `ensmean` of the respective experiment (e.g., `rcp26/ensmean` contains the ensemble mean of the rcp26 scenario).

Similarly, CDO is used to preprocess the CMIP6 data (see `code/preprocess/preprocess_CMIP6.sh`).

The preprocessed landmask based on MPI-GE is available in the data directory (`data/runoff/landmask.nc`), having been produced from raw data using `code/preprocess/make_land_mask.py`.

### References
Maher, N., Milinski, S., Suarez‐Gutierrez, L., Botzet, M., Dobrynin, M., Kornblueh, L., Kröger, J., Takano, Y., Ghosh, R., Hedemann, C., Li, C., Li, H., Manzini, E., Notz, D., Putrasahan, D., Boysen, L., Claussen, M., Ilyina, T., Olonscheck, D., Raddatz, T., Stevens, B., Marotzke, J., 2019. The Max Planck Institute Grand Ensemble: Enabling the Exploration of Climate System Variability. J. Adv. Model. Earth Syst. 11, 2050–2069. https://doi.org/10.1029/2019MS001639

Hurtt, G.C., Chini, L.P., Frolking, S., Betts, R.A., Feddema, J., Fischer, G., Fisk, J.P., Hibbard, K., Houghton, R.A., Janetos, A., Jones, C.D., Kindermann, G., Kinoshita, T., Klein Goldewijk, K., Riahi, K., Shevliakova, E., Smith, S., Stehfest, E., Thomson, A., Thornton, P., van Vuuren, D.P., Wang, Y.P., 2011. Harmonization of land-use scenarios for the period 1500–2100: 600 years of global gridded annual land-use transitions, wood harvest, and resulting secondary lands. Climatic Change 109, 117–161. https://doi.org/10.1007/s10584-011-0153-2

Eyring, V., Bony, S., Meehl, G.A., Senior, C.A., Stevens, B., Stouffer, R.J., Taylor, K.E., 2016. Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. Geosci. Model Dev. 9, 1937–1958. https://doi.org/10.5194/gmd-9-1937-2016
