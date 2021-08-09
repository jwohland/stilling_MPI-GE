CMIP6 pi-control data. Modeling realm Amon. Variable sfcWind. 

Copied from the ETH IAC data pool (/net/atmos/data/cmip6/piControl/Amon/sfcWind) on Feb 18th 2021.
From outside ETH, data can be accessed, for example, via the ESGF.

After downloading, aggregate data to annual means and only keep the first ensemble member (using `code/preprocess_CMIP6.sh`.

CDO throws an error when computing the annual mean 
	(Warning (cdfScanVarAttributes) : NetCDF: Variable not found - areacella) 
which can be safely ignored according to the CDO developers
	(https://code.mpimet.mpg.de/boards/2/topics/3472)
