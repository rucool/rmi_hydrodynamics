## Glider_Data:

glider2composite.ipynb - Glider deployments were retrieved from the Rutgers [ERDDAP](http://slocum-data.marine.rutgers.edu//erddap) server from 2010 to present and subset so that only the portions of the track inside of the nearshore NJ Wind Lease Areas (WLAs) OCS-532, 498, 499, and 549 were analyzed. The deployments were split into three periods during Cold Pool set up (May – June), peak (July – August) and breakdown (September – October). Mixed Layer Depth and Potential Energy Anomaly are calculated to create composite profiles to initialize the carpenter mixing models (steady-state and unsteady). This notebook is available for download only due to the large amount of data being pulled from ERDDAP. 

Metadata for the glider data is provided by ERDDAP.

## Carpenter_Mixing_2016:

This repository is a python implementation of both the steady-state and unsteady models presented in Carpenter et al., (2016) as used in Cassin et al., (2024). These notebook are available to run using colab.


## References:

Cassin, R., Miles, T. N., & Pareja-Roman, L. F. (2024, February). Investigating Potential Impacts of Offshore Wind Turbine Foundations on the Mid-Atlantic Bight Cold Pool. In 2024 Ocean Sciences Meeting. AGU.

Carpenter, J. R., Merckelbach, L., Callies, U., Clark, S., Gaslikova, L., & Baschek, B. (2016). Potential impacts of offshore wind farms on North Sea stratification. PloS one, 11(8), e0160830.
