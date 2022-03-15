# Code for discrete water temperature analyses

This folder of the repository contains the code for two publications on water temperature dynamics in the upper San Francisco Estuary. Some utility functions are found up a directory in the main folder. 

Bashevkin, S. M., and B. Mahardja. 2022. Seasonally variable relationships between surface water temperature and inflow in the upper San Francisco Estuary. Limnology and Oceanography n/a. doi:[10.1002/lno.12027](https://doi.org/10.1002/lno.12027)

Bashevkin, S. M., B. Mahardja, and L. R. Brown. Warming in the upper San Francisco Estuary: Patterns of water temperature change from five decades of data. Limnology and Oceanography n/a. doi:[10.1002/lno.12057](https://doi.org/10.1002/lno.12057)

Both papers use the following code and data files:
1. Data processing: [Temperature data processing.R](<Temperature data processing.R>)
2. Utility functions: [Utility_functions.R](/Utility_functions.R)
3. Map plotting: [Climate change data map.R](<Climate change data map.R>)
4. Delta regions: Shapefile found in the "Delta subregions" folder, this is also identical to [deltamapr::R_EDSM_Subregions_Mahardja](https://github.com/InteragencyEcologicalProgram/deltamapr)

## Temperature-inflow relationships

The code for this paper can be found in the [Drought analysis.R](<Drought analysis.R>) file.

## Warming in the San Francisco Estuary

The code for this paper can be found in the [Climate change modeling.R](<Climate change modeling.R>) and [Data simulation.R](<Data simulation.R>) files.
