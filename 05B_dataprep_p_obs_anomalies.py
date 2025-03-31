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

#%%#####################################################################################################################################
################ PART 2: calculate anomalies at each station #########################################################################
######################################################################################################################################

df = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs.csv')
# drop stupid 29th of february 1818
df = df.drop(df[(df['Year'] == 1818) & (df['Month'] == 2) & (df['Day'] == 29)].index)

df['Date'] = pd.to_datetime(df[['Year', 'Month', 'Day']])

anomals = df.copy()
obsmean = anomals.iloc[:, 3:-1].mean()

# Calculate anomalies
anomals.iloc[:, 3:-1] = anomals.iloc[:, 3:-1] - obsmean[2:]
anomals.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs_anomalies.csv', index=False)


#%% make plots
for i in range(3, anomals.shape[1]-1):
    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], anomals.iloc[:, i], label="Anomalies", alpha=.8)
    plt.axhline(0, linestyle='dashed', color='k', alpha=.5)
    plt.title(f"{anomals.columns[i]} Pressure Series Anomalies")
    plt.ylabel("Pressure [hPa]")
    plt.savefig(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/{df.columns[i]}_station_{i}_anomalyplot.png')
    # plt.show()

# %%
