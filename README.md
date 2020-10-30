# MPA
Magnetic Plasma Analyzer Instrument concept. 
This repository contains the simulation codes used to design this solar wind ion analyzer.

## SIMION codes
This folder contains the codes to set up the SIMION ion-optics simulations. We find a geometry file (.gem) that define the instrument's permanent magnets and multiple lua codes (.fly2) that generates the particular fluxes going through the instrument. A record file (.rec) defines what data is to be recorded and its format.

## Matlab codes
This folder only contains the Matlab/Python codes to map the relative errors made on plasma moments with physical properties of MPA. Run map_generation_xx.m first, then go to the output_date folder and run the related (xx) Python fitting code and then the mapping_xx.m to visualise the color maps.

## Simion-Matlab
This folder contains the Matlab/Python codes to process results from SIMION simulations. The Python code uses the output of one of the Matlab codes and prints the best value after least-square fitting. Some outputs of SIMION are already present in this folder so that solar_winds_maxwellian.m and solar_winds_kappa.m can be executed straight away. The python files give the estimation of plasma moments from the output of flight_kappa.m, flight_maxwellian.m or proton_alpha.m. 
