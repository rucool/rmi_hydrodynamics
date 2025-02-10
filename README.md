# RMI_HYDRODYNAMICS

## Glider_Data:

glider2composite.ipynb - Glider deployments were retrieved from the Rutgers [ERDDAP](http://slocum-data.marine.rutgers.edu//erddap) server from 2010 to present and subset so that only the portions of the track inside of the nearshore NJ Wind Lease Areas (WLAs) OCS-532, 498, 499, and 549 were analyzed. The deployments were split into three periods during Cold Pool set up (May – June), peak (July – August) and breakdown (September – October). Mixed Layer Depth and Potential Energy Anomaly are calculated to create composite profiles to initialize the carpenter mixing models (steady-state and unsteady). 

This notebook is available for download only due to the large amount of data being pulled from ERDDAP. 

Metadata for the glider data is provided by ERDDAP.

## Carpenter_Mixing_2016:

This repository is a python implementation of both the steady-state (steady_state_model.ipynb) and unsteady (unsteady_model.ipynb) models presented in Carpenter et al., (2016) as used in Cassin et al., (2024).  The steady-state and unsteady models are used to approximate how much time the stratification would take to mix, (Tmix) at varying current speeds ranging from 0.01 to 0.8m/s at 0.01 incremenents to represent tidal and storm magnitudes. Experiments are performed at high and low drag coefficients, CD ,  of 0.35 and 1, which account for uncertainty from the monopile roughness, armor layer, among other physical uncertainties in drag. The models utilize the composite profiles generated from glider data as initial conditions.

These notebooks are available to run using colab.

## DOPPIO_DAC

DOPPIO_DAC.ipynb - Extracts depth averaged water velocities across 13 years (2010 – 2023) of numerical model simualtions from the Regional Ocean Modeling System (ROMS) [DOPPIO](https://tds.marine.rutgers.edu/thredds/catalog/catalog.html) at the centroid location of OCS-A 0499. The DOPPIO model is a regional domain at 7km spatial resolution, with 41 vertical layers, designed to represent the complex circulation and hydrodynamics of the MAB coastal areas. The water velocity distribution is compared to the Tmix results from the Carpenter mixing model to highlight the differences in timescales.

This notebook is available to run using colab.

## Work-flow

1. glider2compsite.ipynb: Generate composite potential density profiles from in-situ glider data. This notebook will export a carp_pea_init_profs.csv of profiles for set up peak and breakdown periods.

2. steady_state_model.ipynb and unsteady_model.ipynb:  Initialize both models using carp_pea_init_profs.csv to calculate mixing timescales. These notebooks will output a single csv for all three periods (steady_state_results_H25.csv and unsteady_results_H25.csv).

3. DOPPIO_DAC.ipynb: Use steady_state_results_H25.csv and unsteady_results_H25.csv and the known centroid location of OCS A-0499 to extract depth averaged velocities and generate a histogram comparing the velocity to the mixing timescales.

## References:

Cassin, R., Miles, T. N., & Pareja-Roman, L. F. (2024, February). Investigating Potential Impacts of Offshore Wind Turbine Foundations on the Mid-Atlantic Bight Cold Pool. In 2024 Ocean Sciences Meeting. AGU.

Carpenter, J. R., Merckelbach, L., Callies, U., Clark, S., Gaslikova, L., & Baschek, B. (2016). Potential impacts of offshore wind farms on North Sea stratification. PloS one, 11(8), e0160830.
