# Data access instructions
This directory contains land use data from the Land Use Harmonization (LUH) Project used in CMIP5 experiments ("historic" data downloaded in January 2021, "future" data in February 2021).
It is accessible via https://luh.umd.edu/data.shtml#LUH1_Data.
To populate this directory with historic data, download [the compressed historic data LUHa.v1 file](http://gsweb1vh2.umd.edu/luh_data/LUHa.v1/LUHa.v1.tgz) and expand it into this directory.
To populate the future data subdirectories, the corresponding files for
[IMAGE](http://gsweb1vh2.umd.edu/luh_data/LUHa.v1_future.v1.1/IMAGE/LUHa.v1_image.v1.1.tgz),
[MiniCAM](http://gsweb1vh2.umd.edu/luh_data/LUHa.v1_future.v1/MiniCAM/LUHa.v1_minicam.v1.tgz), and
[MESSAGE](http://gsweb1vh2.umd.edu/luh_data/LUHa.v1_future.v1/MESSAGE/LUHa.v1_message.v1.tgz)
can be downloaded from the same web portal and data stored in subdirectories corresponding to the model names.

To save diskspace, only update states of `gothr` and `gsecd` are kept for all experiments (historic and future); i.e. the folder `lu` can be deleted.
All years prior to 1850 are also deleted.

# Reference
Hurtt, G.C., Chini, L.P., Frolking, S., Betts, R.A., Feddema, J., Fischer, G., Fisk, J.P., Hibbard, K., Houghton, R.A., Janetos, A., Jones, C.D., Kindermann, G., Kinoshita, T., Klein Goldewijk, K., Riahi, K., Shevliakova, E., Smith, S., Stehfest, E., Thomson, A., Thornton, P., van Vuuren, D.P., Wang, Y.P., 2011. Harmonization of land-use scenarios for the period 1500–2100: 600 years of global gridded annual land-use transitions, wood harvest, and resulting secondary lands. Climatic Change 109, 117–161. https://doi.org/10.1007/s10584-011-0153-2