#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:40:17 2024

creates anomalies fitting the first 4 harmonics at every grid cell for 20CRv3 files.

@author: ccorbella@giub.local
"""

import numpy as np
import xarray as xr

#%% load dataset
ds_20CR = xr.open_dataset('/scratch2/ccorbella/files/20CRv3_ensembles/TMP2m.Europe.allmems/TMP2m.Europe.allmems.allyears.nc')
lats = ds_20CR['lat'].values
lons = ds_20CR['lon'].values
t2m  = ds_20CR['TMP2m'].values
time = ds_20CR['time'].values
print(t2m.shape)

#%%
doy = ds_20CR['time'].dt.dayofyear.values

freq = 2 * np.pi / 365.25
s = np.column_stack([ # shape: 5844 (time) x 9 (harmonics)
    np.ones(len(doy)),  # Intercept
    np.sin(freq * doy), np.cos(freq * doy),  # Annual
    np.sin(2 * freq * doy), np.cos(2 * freq * doy),  # Semi-annual
    np.sin(3 * freq * doy), np.cos(3 * freq * doy),  # Tertiary
    np.sin(4 * freq * doy), np.cos(4 * freq * doy)   # Quaternary
])

# Create a new dataset for anomalies
anomalies_ds = xr.Dataset(
    {
        "TMP2m": (
            ("ensemble_member", "time", "lat", "lon"),
            np.full(
                (ds_20CR.dims['ensemble_member'], ds_20CR.dims['time'], ds_20CR.dims['lat'], ds_20CR.dims['lon']),
                np.nan,
                dtype=np.float32
            )
        )
    },
    coords={
        "ensemble_member": ds_20CR["ensemble_member"],
        "time": ds_20CR["time"],
        "lat": ds_20CR["lat"],
        "lon": ds_20CR["lon"]
    }
)

anomalies = xr.full_like(ds_20CR['TMP2m'], fill_value=np.nan)

#%%
from sklearn.linear_model import LinearRegression

for ensmem in range(ds_20CR.dims['ensemble_member']): # loop over ensmems
    # ensmem=1
    t2m_mem = ds_20CR['TMP2m'].isel(ensemble_member=ensmem)
    
    t2m_flat = t2m_mem.values.reshape((t2m_mem.shape[0], -1))  # Flatten lat/lon into one dimension
    seasonal_fit = np.zeros_like(t2m_flat)  # Placeholder for seasonal fits
    
    # Fit seasonal cycle for each grid cell
    for cell_idx in range(t2m_flat.shape[1]):  # Loop over grid cells
        y = t2m_flat[:, cell_idx]
        lr = LinearRegression(fit_intercept=False)  # Match lstsq behavior
        lr.fit(s, y)
        seasonal_fit[:, cell_idx] = lr.predict(s)
    # coeffs = np.linalg.lstsq(s, t2m_flat, rcond=None)[0]  # Fit harmonics, shape: 9 x 8094 = 9 harmonics x 117 lon x 41 lat
    # seasonal_fit = s @ coeffs  # Predict seasonal cycle: 5844 x 8094
    seasonal_fit = seasonal_fit.reshape(t2m_mem.shape)  # Reshape to original grid: (5844, 71, 114)
    
    # Subtract seasonal cycle to get anomalies
    anomalies_ds["TMP2m"][ensmem, :, :, :] = t2m_mem - seasonal_fit

anomalies_ds.to_netcdf('/scratch2/ccorbella/files/20CRv3_ensembles/TMP2m.Europe.allmems/TMP2m.Europe.allmems.anomalies.nc')

# inspect
# anomalies_ds.isel(ensemble_member=1, lat=1, lon=2)['TMP2m'].plot()

