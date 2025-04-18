'''
Gdansk
Lat	54.35
Lon	18.53
Alt 8
lat_idx, lon_idx=(36, 68)
'''

#%%
import xarray as xr
import numpy as np
import pandas as pd
from scipy.spatial import KDTree

station_lon, station_lat = 18.53, 54.35


#%% code to find lat lon index
var_short,var_long='ta','TMP2m'
ds = xr.open_dataset(f'/home/ccorbella/scratch2_symboliclink/files/20CRv3_ensembles/{var_long}.Europe.allmems/{var_long}.Europe.allmems.allyears.nc')
ds_lats = ds['lat'].values
ds_lons = ds['lon'].values
grid_lons, grid_lats = np.meshgrid(ds_lons, ds_lats)
grid_points = np.array([grid_lons.flatten(), grid_lats.flatten()]).T
tree = KDTree(grid_points)
_,station_index=tree.query(np.array([station_lon,station_lat]))
lat_idx = station_index // len(ds_lons)
lon_idx = station_index % len(ds_lons)

# extract ensemble mean at location
ens_mean = ds[var_long].isel(lat=lat_idx, lon=lon_idx).mean(dim='ensemble_member').values

#%% code to debias obs from Gdansk
time = ds['time'].values
df_model = pd.DataFrame({'Gdansk': ens_mean-273.15, # convert to Celsius
                         'Date': pd.to_datetime(time)})

df_obs = pd.read_csv('/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/validation/RAJMUND_Gdansk_18070101-18071231_ta.tsv',
                        skiprows=12, sep='\t')
df_obs['Date'] = pd.to_datetime(df_obs[['Year','Month','Day']])

df_obs['Date'] = df_obs['Date'].dt.floor('D') # round down to the nearest day
df_model['Date'] = df_model['Date'].dt.floor('D') # round down to the nearest day

# keep only the common dates
df_obs2 = df_obs[df_obs['Date'].isin(df_model['Date'])]
df_model = df_model[df_model['Date'].isin(df_obs2['Date'])]
print("original obs set", df_obs.shape,
    "\ncommon days set", df_obs2.shape,
    "\n20CR set", df_model.shape)   

# Set date index
try:
    df_obs2.set_index('Date', inplace=True)
    df_model.set_index('Date', inplace=True)
except:
    print("Date column already set as index")

df_obs2 = df_obs2.loc[:, ['Value']]
df_obs2.columns = ['Gdansk']
#df_obs2 = df_obs[['Value']].rename(columns={'Value':'Gdansk'})

# save values of 20cr at gdansk, in Kelvin!
df_model_K = df_model+273.15

df_model_K.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_20cr_at_gdansk.csv')

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

df_corrected_K = df_corrected + 273.15 # convert to Kelvin

# save corrected obs
df_corrected_K.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs_anomalies_gdansk.csv')

# %%
import matplotlib.pyplot as plt
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

# %%
