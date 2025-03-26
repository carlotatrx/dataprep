#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 12:22:27 2024
Deseasonalization of temperature observation timeseries 

@author: ccorbella@giub.local
"""

#%%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

#%% Load CSV files generated from R
df = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs.csv')

# temporarily drop Yli.x
df = df.drop(columns={'Ylitornio.x'})

df['Date'] = pd.to_datetime(df[['Year', 'Month', 'Day']])
df = df[df['Date'] <= '1849-12-31']
df = df[df['Date'] >= '1806-01-01']
df.index -= df.index[0] # reset indices

doy = df['Date'].dt.dayofyear
freq = 2 * np.pi / 365.25

yobs_plain = df.iloc[:,3:-1].to_numpy()
yobs_anomaly = yobs_plain.copy()
seasonal_cycle = []

#%% fit harmonics
from sklearn.linear_model import LinearRegression

doy = df['Date'].dt.dayofyear
freq = 2 * np.pi / 365.25

yobs_plain = df.iloc[:,3:-1].to_numpy()
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
for i in range(yobs_plain.shape[1]):
    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], yobs_plain[:, i], label="Observations")
    plt.plot(df['Date'], seasonal_cycle[i], label="Seasonal Fit")
    plt.plot(df['Date'], yobs_anomaly[:, i], label="Anomalies", alpha=.8)
    plt.axhline(0, linestyle='dashed', color='k', alpha=.5)
    plt.legend(ncol=3)
    plt.title(f"{df.columns[i+3]} Station - Anomaly Adjustment (1806-1821)")
    plt.xlabel("Date")
    plt.ylabel(r"Temperature [$^{\circ}$C]")
    plt.savefig(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/{df.columns[i+3]}_station_{i}_anomalyplot.png')
    plt.show()

#%% save series
yobs_anomaly_df = pd.DataFrame(yobs_anomaly, index=df['Date'], columns=df.columns[3:-1])
yobs_anomaly_df.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs_anomalies.csv')

