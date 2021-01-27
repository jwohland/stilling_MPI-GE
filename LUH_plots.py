"""
Compute forested area from LUH1 following the the LUH FAQ (https://luh.umd.edu/faq.shtml, Jan 27th 2021):

You can compute an estimate of the forest area in the Land-Use Harmonization products. 
For the original LUH products, first you will need to download the forest/non-forest map 
(fnf_map.txt) and the grid-cell area map (cellarea_halfdeg.txt) in addition to the LUH data 
packages. 
The forest area in an individual grid-cell is given by: 
    (gothr + gsecd)*fnf*cellarea 
(where gothr and gsecd are the fractions of the grid-cell occupied by primary and secondary land respectively). 
To compute the global forest area you would sum this quantity over all grid-cells globally.
"""

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def forest_perarea(gothr, gsecd, fnf):
    """
    Fraction of grid box covered by forest.
    Calculated following the suggestions in the FAQ (see above), then multiplied with cell_area
    because relative evolution appears more relevant here.
    :param gothr:
    :param gsecd:
    :param fnf:
    :return:
    """
    return (gothr + gsecd) * fnf


da_fnf = xr.open_rasterio("fnf_map.txt")

# data_cellarea = pd.read_csv(
#    "cellarea_halfdeg.txt", delimiter=" "
# )  # txt file is missing the header information. Need to define grid manually
#
# data_cellarea = np.genfromtxt("cellarea_halfdeg.txt", delimiter=" ")
# da_cellarea = da_fnf.copy()
# da_cellarea.values[0, :, :] = data_cellarea

years = [str(x) for x in range(1850, 2000, 30)]
for year in years:
    da_gothr = xr.open_rasterio("updated_states/gothr." + year + ".txt")
    da_gsecd = xr.open_rasterio("updated_states/gsecd." + year + ".txt")

    f, ax = plt.subplots(ncols=4, figsize=(12, 2))
    titles = ["da_fnf", "da_gothr", "da_gsecd"]
    for i, da in enumerate([da_fnf, da_gothr, da_gsecd]):
        da.plot(ax=ax[i], vmin=0.01, vmax=1.01, levels=11, extend="neither")
        ax[i].set_title(titles[i])

    da_forest = forest_perarea(da_gothr, da_gsecd, da_fnf)
    da_forest.plot(ax=ax[3], vmin=0.01, vmax=1.01, levels=11, extend="neither")
    ax[3].set_title("forest")
    plt.tight_layout()
    plt.savefig("LUH_" + year + ".jpeg", dpi=300)
