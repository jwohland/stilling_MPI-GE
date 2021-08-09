module load cdo

cd /cluster/work/apatt/wojan/MPI-GE/data/CMIP6 || stop

ANNUALPATH=/cluster/work/apatt/wojan/MPI-GE/data/CMIP6_annual/

for model in *
do
  echo ${model}
  cd ${model}
  ensemble_members=(*)
  ensemble_member=${ensemble_members[0]}  # choose first ensemble member
  cd $ensemble_member
  for file in */*.nc
  do
    outname="${ANNUALPATH}${file#*/}"
    cdo yearmonmean ${file} ${outname}
  done
  cd ../..
done