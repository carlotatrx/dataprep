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


#%% Drop studid 29th February of inexisting years
df_obs = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs.csv')
# drop stupid 29th of february 1807 and 1811
df_obs = df_obs.drop(df_obs[(df_obs['Year'] == 1807) & (df_obs['Month'] == 2) & (df_obs['Day'] == 29)].index)
df_obs = df_obs.drop(df_obs[(df_obs['Year'] == 1811) & (df_obs['Month'] == 2) & (df_obs['Day'] == 29)].index)

df_obs['Date'] = pd.to_datetime(df_obs[['Year', 'Month', 'Day']])
df_obs = df_obs.drop(columns=['Year', 'Month', 'Day'])
df_obs = df_obs[['Date'] + [col for col in df_obs.columns if col != 'Date']] # reorganize columns
df_obs.index -= df_obs.index[0] # reset indices

df_model = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_20cr.csv')
df_model = df_model.rename(columns={'Unnamed: 0':'Date'})
df_model['Date'] = pd.to_datetime(df_model['Date'])
df_model.index -= df_model.index[0] # reset indices
df_model.iloc[:,1:] = df_model.iloc[:,1:] - 273.15 # convert to Celsius

df_obs['Date'] = df_obs['Date'].dt.floor('D') # round down to the nearest day
df_model['Date'] = df_model['Date'].dt.floor('D') # round down to the nearest day

# keep only the common dates
df_obs2 = df_obs[df_obs['Date'].isin(df_model['Date'])]
print("original obs set", df_obs.shape,
      "\ncommon days set", df_obs2.shape,
      "\n20CR set", df_model.shape)

duplicate_dates = df_obs2['Date'].duplicated(keep=False)
df_obs2_duplicates = df_obs2[duplicate_dates]
print(f"Number of duplicated dates: {df_obs2_duplicates.shape[0]}")

unique_dates = df_obs2['Date'].nunique()
print(f"Number of unique dates: {unique_dates}")

# Set date index
try:
    df_obs2.set_index('Date', inplace=True)
    df_model.set_index('Date', inplace=True)
except:
    print("Date column already set as index")

#%% try linear regression
import statsmodels.api as sm

# Day-of-year and frequency
doy = df_obs2.index.dayofyear.to_numpy()
freq = 2 * np.pi / 365.25

# Container for corrected observations
df_corrected = pd.DataFrame(index=df_obs2.index, columns=df_obs2.columns, dtype=float)

for station in df_obs2.columns:
    obs = df_obs2[station]
    model = df_model[station]

    valid_obs = obs.notna()
    valid_model = model.notna()

    if valid_obs.sum() < 365:
        print(f"{station}: <1 year of obs data — skipped")
        continue

    valid_years = round(valid_obs.sum() / 365.25)

    # Select harmonic terms
    if valid_years < 2:
        harmonics = np.column_stack([
            np.ones_like(doy),
            np.sin(freq * doy),
            np.cos(freq * doy)
        ])
        print(f"{station}: using 1st harmonic")
    else:
        harmonics = np.column_stack([
            np.ones_like(doy),
            np.sin(freq * doy), np.cos(freq * doy),
            np.sin(2 * freq * doy), np.cos(2 * freq * doy),
            np.sin(3 * freq * doy), np.cos(3 * freq * doy),
            np.sin(4 * freq * doy), np.cos(4 * freq * doy)
        ])
        print(f"{station}: using 4 harmonics")

    # Fit seasonal cycle to obs
    X_obs = harmonics[valid_obs]
    y_obs = obs[valid_obs].to_numpy()
    beta_obs = sm.OLS(y_obs, X_obs).fit().params
    seasonal_obs = harmonics @ beta_obs

    # Fit seasonal cycle to model (if sufficient data)
    if valid_model.sum() >= 365:
        X_model = harmonics[valid_model]
        y_model = model[valid_model].to_numpy()
        beta_model = sm.OLS(y_model, X_model).fit().params
        seasonal_model = harmonics @ beta_model
    else:
        seasonal_model = seasonal_obs.copy()
        print(f"{station}: insufficient model data, using obs seasonal cycle instead")

    # Apply bias correction
    corrected = obs.to_numpy() + (seasonal_model - seasonal_obs)

    df_corrected[station] = corrected

df_corrected_K = df_corrected.iloc[:, 1:] + 273.15 # convert to Kelvin

# save corrected obs
df_corrected_K.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs_anomalies.csv')

#%% make plots
# Plot anomalies
stations = df_obs2.columns

for station in stations:
    plt.figure(figsize=(12, 4))
    plt.plot(df_obs2.index, df_obs2[station], label='Observation', color='dimgray', alpha=.5, linewidth=1)
    plt.plot(df_model.index, df_model[station], label='20CRv3 Ensemble Mean', color='mediumvioletred', alpha=.5, linewidth=1)
    plt.plot(df_corrected.index, df_corrected[station], label='Corrected', color='orangered', alpha=.7, linewidth=1)
    plt.title(f"Temperature at {station}, Observed vs 20CRv3")
    plt.ylabel(r'Temperature [$^{\circ}$C]')
    plt.legend()
    plt.tight_layout()

    outfile = f"/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/ta_{station}_anomalyplot.png"
    plt.savefig(outfile)
    plt.close()
    print(f"📈 Saved: {station}")
