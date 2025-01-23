#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:40:17 2024

creates anomalies by subtracting the mean pressure data at every grid cell for 20CRv3 files.

@author: ccorbella@giub.local
"""

import xarray as xr

#%% load dataset
ds_20CR = xr.open_dataset('/scratch2/ccorbella/files/20CRv3_ensembles/PRMSL.Europe.allmems/PRMSL.Europe.allmems.allyears.nc')

# Extract PRMSL data
prmsl = ds_20CR['PRMSL']

# Calculate the mean over time for each lat, lon, and ensemble member
prmsl_mean = prmsl.mean(dim='time', skipna=True)  # Mean over time (ignores NaNs)

# Calculate anomalies by subtracting the mean
prmsl_anomalies = prmsl - prmsl_mean

# Create a new dataset for anomalies
anomalies_ds = prmsl_anomalies.to_dataset(name='PRMSL')

# Save anomalies to a new NetCDF file
output_path = '/scratch2/ccorbella/files/20CRv3_ensembles/PRMSL.Europe.allmems/PRMSL.Europe.allmems.anomalies.nc'
anomalies_ds.to_netcdf(output_path)
