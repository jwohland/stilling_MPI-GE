# Data access instructions
CMIP6 pi-control data. Modeling realm `Amon`. Variable `sfcWind`.

Data was copied from the ETH IAC data pool (/net/atmos/data/cmip6/piControl/Amon/sfcWind) on Feb 18th 2021.

CDO throws an error when computing the annual mean
        (Warning (cdfScanVarAttributes) : NetCDF: Variable not found - areacella)
which can be safely ignored according to the CDO developers
        (https://code.mpimet.mpg.de/boards/2/topics/3472)

Only annual mean data is required in this analysis, populated in the directory `CMIP6_annual` using the script `preprocess_CMIP6_annual.sh`. Once `CMIP6_annual` has been populated, this directory is no longer required. To point to `CMIP6` data stored elsewhere, you can update the `preprocess_CMIP6_annual.sh` accordingly.
# Reference
Eyring, V., Bony, S., Meehl, G. A., Senior, C. A., Stevens, B., Stouffer, R. J., and Taylor, K. E.: Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization, Geosci. Model Dev., 9, 1937-1958, doi:10.5194/gmd-9-1937-2016, 2016.