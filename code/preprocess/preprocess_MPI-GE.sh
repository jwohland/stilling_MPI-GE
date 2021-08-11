# First argument to the
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)
cd "$script_dir/../../data" || stop

for scenario in 1pCO2 historical rcp26 rcp45 rcp85
do
  echo ${scenario}
  cd ${scenario}
  one_name=(sfcWind*.nc)  # a file name
  realm=(`echo $one_name | tr '_' ' '`)  # split file names by "_" to identify realm (i.e., Lmon or Amon)
  cdo ensmean sfcWind*.nc ensmean/sfcWind_${realm[1]}_MPI-ESM_1pctCO2_ensmean_185001-199912_test.nc
  cd ..
done