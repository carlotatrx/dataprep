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
#%%
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

#%%

# Set date index
try:
    df_obs2.set_index('Date', inplace=True)
    df_model.set_index('Date', inplace=True)
except:
    print("Date column already set as index")

# Calculate model - obs difference
delta = df_model - df_obs2

# Compute mean model–obs difference for each station
delta_mean = delta.mean(axis=0, skipna=True)

# Apply correction: for each time point, add the mean delta
df_corrected = df_obs2.add(delta_mean)
df_corrected.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs_anomalies.csv')

#%% make plots
# Plot anomalies
stations = df_obs2.columns

for station in stations:
    plt.figure(figsize=(12, 4))
    plt.plot(df_obs2.index, df_obs2[station], label='Observation', color='dimgray', alpha=.5, linewidth=1)
    plt.plot(df_model.index, df_model[station], label='20CRv3 Ensemble Mean', color='mediumvioletred', alpha=.5, linewidth=1)
    plt.plot(df_corrected.index, df_corrected[station], label='Corrected', color='royalblue', alpha=.7, linewidth=1)
    plt.title(f"Pressure at {station}, Observed vs 20CRv3")
    plt.ylabel('pressure (hPa)')
    plt.legend()
    plt.tight_layout()

    outfile = f"/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/p_{station}_anomalyplot.png"
    plt.savefig(outfile)
    plt.close()
    print(f"📈 Saved: {station}")
