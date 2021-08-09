# MPI-GE wind stilling

Code underlying analysis performed in 

Wohland, J., Folini, D., Pickering, B., 2021. Wind speed stilling and its recovery due to internal climate variability (preprint). Under review for Earth System Dynamics. https://doi.org/10.5194/esd-2021-29

## Figure overview

| Figure | Filename | Creating python script |
|---|---|---|
Fig 1 | attribution_maps.jpeg | attribution.py
Fig. 2 | contribution_histograms.jpeg | attribution.py
Fig. 3 | scatter_gothr+gsecd_abs.jpeg | attribution.py
Fig. 4 | LUH_change_future_ref2000.jpeg | LUH_plots.py
Fig. 5 | global_windspeeds.jpeg | trend_maps.py
Fig. 6 | timeseries_picontrol_Europe.jpeg | trends.py
Fig. 7a | picontrol_wind_trends_Europe_5_20y.jpeg | trends.py
Fig. 7b | Ensmean_picontrol_wind_trends_Europe_5.jpeg | trends.py
Fig. 8a | historical_wind_trends_Europe_5_all.jpeg | trends.py
Fig. 8b | rcp26_wind_trends_Europe_5_all.jpeg  | trends.py
Fig. 8c | rcp45_wind_trends_Europe_5_all.jpeg | trends.py
Fig. 8d | rcp85_wind_trends_Europe_5_all.jpeg | trends.py
Appendix Fig. A1| box_timeseries_Europe.jpeg | ensmean_timeseries.py
Appendix Fig. A2| picontrol_HadISD_wind_trends_Europe_5_20y.jpeg | trends.py
Appendix Fig. A3a | picontrol_wind_trends_Europe_5_15y.jpeg | trends.py
Appendix Fig. A3b | picontrol_wind_trends_Europe_5_25y.jpeg | trends.py
Appendix Fig. A4a | picontrol_wind_trends_Europe_10_20y.jpeg | trends.py
Appendix Fig. A4b | picontrol_wind_trends_Europe_15_20y.jpeg | trends.py
Appendix Figs. 5-7 | {Model_name_picontrol_wind_trends_Europe_5.jpeg} | trends.py
## Other Files

make_land_mask.py computes a land mask on the MPI grid using runoff data and excluding Antartica and Greenland. 

Lmon_Amon_diff.py verify that wind speed data provided in the modeling realm Atmosphere (Amon) and Land (Lmon) contain the same values. 

## Input data
Executing scripts from this repository requires (a) downloading of the required input data and (b) preprocessing the data, and (c) updating the file paths. 

### a) Downloads
To download the input, please follow the structure and recommendations given in the `data` subdirectories. 
The MPI-GE data can be downloaded using the provided `wget` scripts after you have registered at ESGF (see https://esgf-data.dkrz.de/user/add/?next=http://esgf-data.dkrz.de/projects/esgf-dkrz/ to create a new profile).
The other downloads are explained in `readme.txt` files in the respective subfolders. 

Input data is from 3 different sources:
- the Max Planck Institute for Meteorology (MPI-M) Grand Ensemble (MPI-GE) that can be downloaded at https://esgf-data.dkrz.de/search/mpi-ge/ 
- the Land Use Harmonization (LUH) Project that is accessable via https://luh.umd.edu/data.shtml#LUH1_Data.
- the Climate Model Intercomparison Project Phase 6 and can be downloaded from: https://esgf-data.dkrz.de/search/cmip6-dkrz/

The MPI-GE based landmask is available in `data/runoff`.

###b) Preprocessing
To speed up Python calculations, the MPI-GE ensemble means are pre-calculated using CDO (see `code/preprocess/preprocess_MPI-GE`). 
They are stored in a subfolder `ensmean` of the respective experiment (e.g., `rcp26/ensmean` contains the ensemble mean of the rcp26 scenario).

The CMIP6 preprocessing is performed in `code/preprocess/preprocess.CMIP6`.



### References
Maher, N., Milinski, S., Suarez‐Gutierrez, L., Botzet, M., Dobrynin, M., Kornblueh, L., Kröger, J., Takano, Y., Ghosh, R., Hedemann, C., Li, C., Li, H., Manzini, E., Notz, D., Putrasahan, D., Boysen, L., Claussen, M., Ilyina, T., Olonscheck, D., Raddatz, T., Stevens, B., Marotzke, J., 2019. The Max Planck Institute Grand Ensemble: Enabling the Exploration of Climate System Variability. J. Adv. Model. Earth Syst. 11, 2050–2069. https://doi.org/10.1029/2019MS001639

Hurtt, G.C., Chini, L.P., Frolking, S., Betts, R.A., Feddema, J., Fischer, G., Fisk, J.P., Hibbard, K., Houghton, R.A., Janetos, A., Jones, C.D., Kindermann, G., Kinoshita, T., Klein Goldewijk, K., Riahi, K., Shevliakova, E., Smith, S., Stehfest, E., Thomson, A., Thornton, P., van Vuuren, D.P., Wang, Y.P., 2011. Harmonization of land-use scenarios for the period 1500–2100: 600 years of global gridded annual land-use transitions, wood harvest, and resulting secondary lands. Climatic Change 109, 117–161. https://doi.org/10.1007/s10584-011-0153-2

Eyring, V., Bony, S., Meehl, G.A., Senior, C.A., Stevens, B., Stouffer, R.J., Taylor, K.E., 2016. Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. Geosci. Model Dev. 9, 1937–1958. https://doi.org/10.5194/gmd-9-1937-2016
