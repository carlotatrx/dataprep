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
plt.rcParams.update({
    'font.size': 12,               # base font size
    'axes.titlesize': 12,          # figure title
    'axes.labelsize': 12,          # x and y labels
    'xtick.labelsize': 12,         # x tick labels
    'ytick.labelsize': 12,         # y tick labels
    'legend.fontsize': 12,         # legend text
    'axes.xmargin': 0.01,            # x axis margin
    'figure.titlesize': 12         # overall figure title (suptitle)
})

#%% Load CSV files generated from R
df = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs.csv')

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

yobs_anomaly_df = pd.DataFrame(yobs_anomaly, index=df['Date'], columns=df.columns[3:-1])
yobs_anomaly_df.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs_anomalies.csv')

# #%% rename columns for plots
# # modify names for temperature:
# dict_rename = {
#     'AG01_Aarau_Zschogge': 'Aarau',
#     '00033902': 'Kherson',
#     '00033345': 'Kyiv',
#     '00033902_SUBs': 'Kherson_SUBs',
#     'BE01_Bern_Studer': 'Bern',
#     'Burg_Zitenice': 'Zitenice',
#     'Camuffo_Bologna': 'Bologna_Camuffo',
#     'GR04_Marschlins': 'Schloss Marschlins',
#     'JU01_Delemont': 'Delemont',
#     'SH01_Schaffhausen': 'Schaffhausen',
#     'VD01_Vevey': 'Vevey',
#     'DigiHom_Geneva': 'Geneva',
#     'KIT_Karlsruhe': 'Karlsruhe',
#     'Dom_Valencia_SUBs': 'Valencia_SUBs',
#     'KNMI-42_Zwanenburg': 'Zwanenburg',
#     'KNMI-44_Haarlem': 'Haarlem',
#     'KNMI-44_Haarlem_SUBs': 'Haarlem_SUBs',
#     'KNMI-45_Delfft2': 'Delfft',
#     'Europe_Rovereto_1': 'Rovereto',
#     'GCOS_Zurich_Feer': 'Zurich',
#     'IMPROVE_Cadiz': 'Cadiz',
#     'IMPROVE_Stockholm': 'Stockholm',
#     'IMPROVE_StPetersburg': 'StPetersburg_IMPROVE',
#     'Europe_StPetersburg': 'StPetersburg_Europe',
#     'IMPROVE_Uppsala': 'Uppsala',
#     'CRU_Paris': 'Paris',
#     'Dom_Valencia': 'Valencia',
#     'IMPROVE_Milan': 'Milan',
#     'HadCET': 'CET',
#     'Brug_Zitenice': 'Zitenice',
#     'HOH': 'Hohenpeisenberg',
#     'TOR': 'Torino',
#     'GVE': 'Geneva',
#     '169': 'Bologna_ECAD'
# }

# df  = df.rename(columns=dict_rename)

#%% make plots
for i in range(yobs_plain.shape[1]):
    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], yobs_plain[:, i], label="Observations")
    plt.plot(df['Date'], seasonal_cycle[i], label="Seasonal Fit")
    plt.plot(df['Date'], yobs_anomaly[:, i], label="Anomalies", alpha=.8)
    plt.axhline(0, linestyle='dashed', color='k', alpha=.5)
    plt.legend(ncol=3)
    plt.title(f"{df.columns[i+3]} Temperature Series")
    plt.xlabel("Date")
    plt.ylabel(r"Temperature [$^{\circ}$C]")
    plt.savefig(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/ta_{df.columns[i+3]}_anomalyplot.png')
    # plt.show()


# %%
