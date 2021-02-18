# MPI-GE wind stilling

Code underlying analysis performed in TITLE NAME.


## Figure overview

| Figure | Filename | Creating python script |
|---|---|---|
Fig 1 | attribution_maps.jpeg | attribution.py
Fig. 2a | contribution_histograms.jpeg | attribution.py
Fig. 2b | scatter_gothr+gsecd_abs.jpeg | attribution.py
Fig. 3 | LUH_change_future_ref2000.jpeg | LUH_plots.py
Fig. 4 | global_windspeeds.jpeg | trend_maps.py
Fig. 5 | timeseries_picontrol_Europe.jpeg | trends.py
Fig. 6 | box_timeseries_Europe.pdf | ensmean_timeseries.py
Appendix Fig. A1|  LUH_change.jpeg | LUH_plots.py

## Other Files

make_land_mask.py computes a land mask on the MPI grid using runoff data and excluding Antartica and Greenland. 

Lmon_Amon_diff.py verify that wind speed data provided in the modeling realm Atmosphere (Amon) and Land (Lmon) contain the same values. 

## Input data

Input data is from the Max Planck Institute for Meteorology (MPI-M) Grand Ensemble (MPI-GE) that can be downloaded at https://esgf-data.dkrz.de/search/mpi-ge/ and from the Land Use Harmonization (LUH) Project that is accessable via https://luh.umd.edu/data.shtml#LUH1_Data.

Maher, N., Milinski, S., Suarez‐Gutierrez, L., Botzet, M., Dobrynin, M., Kornblueh, L., Kröger, J., Takano, Y., Ghosh, R., Hedemann, C., Li, C., Li, H., Manzini, E., Notz, D., Putrasahan, D., Boysen, L., Claussen, M., Ilyina, T., Olonscheck, D., Raddatz, T., Stevens, B., Marotzke, J., 2019. The Max Planck Institute Grand Ensemble: Enabling the Exploration of Climate System Variability. J. Adv. Model. Earth Syst. 11, 2050–2069. https://doi.org/10.1029/2019MS001639

Hurtt, G.C., Chini, L.P., Frolking, S., Betts, R.A., Feddema, J., Fischer, G., Fisk, J.P., Hibbard, K., Houghton, R.A., Janetos, A., Jones, C.D., Kindermann, G., Kinoshita, T., Klein Goldewijk, K., Riahi, K., Shevliakova, E., Smith, S., Stehfest, E., Thomson, A., Thornton, P., van Vuuren, D.P., Wang, Y.P., 2011. Harmonization of land-use scenarios for the period 1500–2100: 600 years of global gridded annual land-use transitions, wood harvest, and resulting secondary lands. Climatic Change 109, 117–161. https://doi.org/10.1007/s10584-011-0153-2
