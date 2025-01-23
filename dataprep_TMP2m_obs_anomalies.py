#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 12:22:27 2024

@author: ccorbella@giub.local
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

#%% Load CSV files generated from R
csv_filepath = '/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/'
df_all  = pd.read_csv(f'{csv_filepath}station_data.csv')
metadata = pd.read_csv(f'{csv_filepath}station_metadata.csv', index_col=0)

# station data might have repeated columns depending on the code in R
for col in df_all.columns:
    if col + ".y" in df_all.columns:
        if df_all[col].equals(df_all[col + ".y"]):      # Check if all values are equal in both columns
            df_all = df_all.drop(columns=[col + ".y"])  # Drop the ".y" column if they are equal

# one time operation 
if (df_all.columns[0]!='date'):
    df_all = df_all.drop(df_all.columns[[0]],axis=1) # drop stupid first column
    df_all = df_all.drop(columns=['IMPROVE_StPetersburg'])

df_all['date'] = pd.to_datetime(df_all['date'])

# crop the data
df_all = df_all[df_all['date'] <= '1821-12-31']
df_all = df_all[df_all['date'] >= '1806-01-01']
df_all.index -= df_all.index[0] # reset indices

#%% fit harmonics
from sklearn.linear_model import LinearRegression

doy = df_all['date'].dt.dayofyear
freq = 2 * np.pi / 365.25

yobs_plain = df_all.iloc[:,1:].to_numpy()
yobs_anomaly = yobs_plain.copy()
seasonal_cycle = []

for i in range(yobs_plain.shape[1]): # loop through stations
    obs_col = yobs_plain[:, i]
    valid_mask = ~np.isnan(obs_col)
    
    if valid_mask.sum() == 0:  # No valid data for this station
        yobs_anomaly[:, i] = np.nan
        seasonal_cycle.append(np.full_like(obs_col, np.nan))  # Fill with NaNs for consistency
        print("station ", i, " has too little data, no fitting")
        continue

    valid_years = round(valid_mask.sum() / 365.25)  # Estimate years of data

    if valid_years < 1:  # Too little data; exclude
        yobs_anomaly[:, i] = np.nan
        sc = np.full_like(obs_col, np.nan)
        print("station ", i, " has too little data, no fitting")
        
    elif valid_years < 2:  # Use only 1st harmonic for stations with <2 years of data
        s_simple = np.column_stack([
            np.ones(doy.shape),
            np.sin(freq * doy), np.cos(freq * doy)
        ])
        lm_obs = LinearRegression().fit(s_simple[valid_mask], obs_col[valid_mask])
        sc_full = lm_obs.predict(s_simple)
        sc = np.full_like(obs_col, np.nan)  # Initialize as NaNs
        sc[valid_mask] = sc_full[valid_mask]  # Constrain to valid dates
        yobs_anomaly[:, i] = obs_col - sc
        
        print("only first harmonic fit for station", i, lm_obs.coef_)
        
    else:  # Use full harmonics for stations with sufficient data
        s = np.column_stack([
            np.ones(doy.shape),
            np.sin(freq * doy), np.cos(freq * doy),
            np.sin(2 * freq * doy), np.cos(2 * freq * doy),
            np.sin(3 * freq * doy), np.cos(3 * freq * doy),
            np.sin(4 * freq * doy), np.cos(4 * freq * doy)
        ])
        lm_obs = LinearRegression().fit(s[valid_mask], obs_col[valid_mask])
        sc_full = lm_obs.predict(s)
        sc = np.full_like(obs_col, np.nan)  # Initialize as NaNs
        sc[valid_mask] = sc_full[valid_mask]  # Constrain to valid dates
        yobs_anomaly[:, i] = obs_col - sc
        print("four harmonics for station", i, lm_obs.coef_)
    seasonal_cycle.append(sc)


        
#%% make plots
station_index = 27  # Example station
for i in range(yobs_plain.shape[1]):
    plt.figure(figsize=(12, 6))
    plt.plot(df_all['date'], yobs_plain[:, i], label="Observations")
    plt.plot(df_all['date'], seasonal_cycle[i], label="Seasonal Fit")
    plt.plot(df_all['date'], yobs_anomaly[:, i], label="Anomalies", alpha=.8)
    plt.axhline(0, linestyle='dashed', color='k', alpha=.5)
    plt.legend(ncol=3)
    plt.title(f"{df_all.columns[i+1]} Station - Anomaly Adjustment (1806-1821)")
    plt.xlabel("Date")
    plt.ylabel(r"Temperature [$^{\circ}$C]")
    plt.savefig(f'/scratch2/ccorbella/files/1807_USBstick/station_anomalies_plots/{df_all.columns[i+1]}_station_{i}_anomalyplot.png')
    plt.show()

#%% save series
np.savetxt('/scratch2/ccorbella/code/KF_assimilation/obs_anomalies_1806-1821.txt', yobs_anomaly)

