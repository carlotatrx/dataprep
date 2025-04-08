"""
Process TSV files containing pressure data, extracts metadata and data, 
combines them into a single DataFrame, and pivots the data to have station names as columns.
Resulting dataframe looks like this:
Station      Paris     Valencia     Zitenice
Date                                        
1807-01-01  1037.8  1008.024602   995.514928
1807-01-02  1042.0  1007.521988  1000.281213
1807-01-03  1033.3  1004.523411   992.202052
...
1807-12-30  1033.3  1004.523411   992.202052
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.size': 12,               # base font size
    'axes.titlesize': 12,          # figure title
    'axes.labelsize': 12,          # x and y labels
    'xtick.labelsize': 12,         # x tick labels
    'ytick.labelsize': 12,         # y tick labels
    'legend.fontsize': 12,         # legend text
    'axes.xmargin': 0.01,          # x axis margin
    'figure.titlesize': 12         # overall figure title (suptitle)
})

#%%#####################################################################################################################################
################ PART 2: calculate anomalies at each station #########################################################################
######################################################################################################################################

## Need to have the updated obs file!! created w 04_dataprep_tsv_obs.R 

df_obs = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs.csv')
# drop stupid 29th of february 1818
df_obs = df_obs.drop(df_obs[(df_obs['Year'] == 1818) & (df_obs['Month'] == 2) & (df_obs['Day'] == 29)].index)

df_obs['Date'] = pd.to_datetime(df_obs[['Year', 'Month', 'Day']])
df_obs = df_obs.drop(columns=['Year', 'Month', 'Day'])
df_obs = df_obs[['Date'] + [col for col in df_obs.columns if col != 'Date']]
df_obs.index -= df_obs.index[0] # reset indices

df_model = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_20cr.csv')
df_model = df_model.rename(columns={'Unnamed: 0':'Date'})
df_model['Date'] = pd.to_datetime(df_model['Date'])
df_model.index -= df_model.index[0] # reset indices
df_model.iloc[:, 1:] = df_model.iloc[:, 1:]/100 # convert to hPa

# Merge the two dataframes on the 'Date' column keeping only the common dates
df_obs['Date'] = df_obs['Date'].dt.floor('D') # round down to the nearest day
df_model['Date'] = df_model['Date'].dt.floor('D') # round down to the nearest day

# keep only the common dates
df_obs2 = df_obs[df_obs['Date'].isin(df_model['Date'])]


#%%
yobs_plain = df.iloc[:,3:-1].to_numpy()
# yobs_anomaly = yobs_plain.copy()
yobs_mean = np.nanmean(yobs_plain, axis=0)

# Calculate anomalies
yobs_anomaly = yobs_plain - yobs_mean

yobs_anomaly_df = pd.DataFrame(yobs_anomaly, index=df['Date'], columns=df.columns[3:-1])
yobs_anomaly_df.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs_anomalies.csv')

#%% make plots
# Plot anomalies
for i in range(yobs_anomaly.shape[1]): # loop through stations
    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], yobs_anomaly[:, i], label="Anomalies", alpha=.8)
    plt.axhline(0, linestyle='dashed', color='k', alpha=.5)
    plt.title(f"{df.columns[i+3]} Pressure Series Anomalies")
    plt.ylabel("Pressure [hPa]")
    plt.tight_layout()
    plt.savefig(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/p_{df.columns[i+3]}_anomalyplot.png')
    # plt.show()

# %%
