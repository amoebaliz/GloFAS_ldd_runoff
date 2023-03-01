# GloFAS_ldd_runoff
This repository contains a python script for writing coastal runoff fields for a regional mom6 model using freswhater discharge (FWD) in the Global Flood Awareness System (GloFAS) based on the local drain direction (ldd). The bash script allows for submitting multiple jobs to expedite processing time (assumes multiple processors are available). 

The python script was adapted from code written by Andrew C. Ross @ GFDL


## Special Software Requirements

xarray

xesmf


## Input Requirements:
GloFAS reanalysis output (netCDF files - recommend spatial subset for region of interest): 
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical?tab=overview

GloFAS ldd map (netCDF file):
https://confluence.ecmwf.int/display/CEMS/Auxiliary+Data?preview=%2F242067380%2F304227214%2Fldd_glofas.nc

MOM6 ocean_hgrid (netcdf file)

MOM6 ocean_mask (netcdf file)

