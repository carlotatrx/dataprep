import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from scipy.stats import pearsonr
from scipy.spatial import KDTree

# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def make_date_col(df):
    df['Date'] = pd.to_datetime(df[['Year', 'Month', 'Day']])

def weighted_tod(hour, center=14, std=4):
    return np.exp(-0.5 * ((hour - center) / std)**2)

# -----------------------------------------------------------------------------
# Core functions
# -----------------------------------------------------------------------------

def create_df_merged(lat_station, lon_station, ds_20CR, filename):
    """
    Load station data, check if subdaily, and merge with 20CRv3 data.
    Returns merged DataFrame with columns: Date, model, obs.
    """
    # Load SEF file
    df = pd.read_csv(filename, skiprows=12, sep='\t')

    # Create date and optionally time columns
    df["Date"] = pd.to_datetime(df[["Year", "Month", "Day"]])
    df["Hour"] = pd.to_numeric(df["Hour"], errors="coerce")

    # Detect if subdaily: more than one row per day
    subdaily = df["Date"].duplicated(keep=False).any()

    if subdaily:
        print("Subdaily data detected. Applying weighted daily mean.")

        # Filter out invalid hours
        df = df[df["Hour"].notna()]

        # Compute adaptive weights
        df["weight"] = df["Hour"].apply(lambda hr: weighted_tod(hr))

        # Normalize weights per day
        df["norm_weight"] = df.groupby("Date")["weight"].transform(lambda w: w / w.sum())

        # Compute weighted daily means
        df_obs = df.groupby("Date").apply(
            lambda g: np.average(g["Value"], weights=g["norm_weight"]),
            include_groups=False # silence stupid warning
        ).reset_index(name='obs')

    else:
        df_obs = df[['Date', 'Value']]
        print(f"Daily data detected.")
        
    # Ensure unique dates
    #    df = df.drop_duplicates(subset="Date")

    # Nearest grid point
    ds_lats = ds_20CR['lat'].values
    ds_lons = ds_20CR['lon'].values
    grid_lons, grid_lats = np.meshgrid(ds_lons, ds_lats)
    grid_points = np.array([grid_lats.flatten(), grid_lons.flatten()]).T
    tree = KDTree(grid_points)
    _, station_index = tree.query(np.array([lat_station, lon_station]).T)
    lat_index = station_index // len(ds_lons)
    lon_index = station_index % len(ds_lons)

    print(f"matched to 20CR grid: lat={ds_lats[lat_index]}, lon={ds_lons[lon_index]}")

    # Extract and average model data
    ds_ta = ds_20CR['TMP2m'][:, :, lat_index, lon_index].mean(dim='ensemble_member')
    df_20CR = ds_ta.to_dataframe().reset_index()
    df_20CR = df_20CR.rename(columns={'time': 'Date', 'TMP2m': 'model'})
    df_20CR = df_20CR.drop(columns=['lat', 'lon'])
    
    # make datecols equal
    df_20CR["Date"] = pd.to_datetime(df_20CR["Date"]).dt.normalize()
    df_obs['Date'] = pd.to_datetime(df_obs['Date']).dt.normalize()
    
    # Merge on date
    df_merged = pd.merge(df_20CR, df_obs, on="Date", how="inner")
    df_merged['model'] -= 273.15  # Convert from Kelvin to Celsius
    return df_merged

def compute_lag_corr(df, max_lag=14):
    """
    Computes Pearson correlation for various time lags.
    """
    lags = np.arange(-max_lag, max_lag + 1)
    corrs = []

    for lag in lags:
        shifted = df.copy()
        shifted['lagged_model'] = shifted['model'].shift(lag)  # Drop NaNs introduced by shifting
        shifted = shifted.dropna(subset=['obs', 'lagged_model'])
        # ensure enough data points for correlation
        corr = pearsonr(shifted['obs'], shifted['lagged_model'])[0] if len(shifted) > 10 else np.nan
        corrs.append(corr)

    return lags, corrs


def rolling_lag_corr(df, window_size=365*3, step_size=90, max_lag=20):
    """
    Computes Pearson correlation for rolling time windows at different lags.
    """
    if len(df) < window_size:
        print("Warning: dataset shorter than window. Reducing window size.")
        window_size = len(df)

    lags = np.arange(-max_lag, max_lag + 1)
    results = []

    for start in range(0, len(df) - window_size, step_size):
        sub_df = df.iloc[start:start + window_size]
        date_center = sub_df.iloc[window_size // 2]['Date']
        best_corr = -np.inf
        best_lag = 0

        for lag in lags:
            shifted = sub_df.copy()
            shifted['Lagged_model'] = shifted['model'].shift(lag)
            valid = shifted.dropna(subset=['obs', 'Lagged_model'])

            if len(valid) > 10:
                corr = pearsonr(valid['obs'], valid['Lagged_model'])[0]
                if corr > best_corr:
                    best_corr = corr
                    best_lag = lag

        results.append({'Date': date_center, 'Best_Lag': best_lag, 'Best_Corr': best_corr})

    return pd.DataFrame(results)


def plot_lag_correlation(lags, corrs, station_name,
                         file_savefig='/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/Ukraine_calendar_shifts'):
    plt.figure(figsize=(8, 5))
    ax = plt.gca()
    plt.plot(lags, corrs, marker="o", linestyle="-")
    plt.axvline(0, color="black", linestyle="--", label="No Shift")
    max_idx = np.nanargmax(corrs)
    plt.plot(lags[max_idx], corrs[max_idx], marker="o", color="red")
    plt.xlabel("Lag (Days)")
    plt.ylabel("Pearson Correlation")
    plt.text(0.05, 0.92, f'max: {corrs[max_idx]:.2f}', transform=ax.transAxes, fontsize=14)
    plt.title(f"Correlation {station_name} observations and 20CRv3, TMP2m")
    plt.legend()
    plt.grid()
    plt.savefig(f'{file_savefig}/{station_name}_lag_correlation.png', dpi=300)
    plt.close()


def plot_rolling_lag(df_result, station_name, expected_shift=None,
                     file_savefig='/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/Ukraine_calendar_shifts'):
    plt.figure(figsize=(12, 6))
    plt.plot(df_result['Date'], df_result['Best_Lag'], marker='o', linestyle='-', label='Best Lag')
    if expected_shift is not None:
        plt.axhline(expected_shift, color="r", linestyle="--", label=f"Expected Shift ({expected_shift} Days)")
    plt.xlabel("Time")
    plt.ylabel("Best Lag (Days)")
    plt.title(f"Changes in Lag Over Time, {station_name}")
    plt.legend()
    plt.grid()
    plt.savefig(f'{file_savefig}/{station_name}_rolling_lag.png', dpi=300)
    plt.close()

ds_20CR = xr.open_dataset('/scratch2/ccorbella/files/20CRv3_ensembles/TMP2m.Europe.allmems/TMP2m.Europe.allmems.allyears.nc')
lat_Kamyanets, lon_Kamyanets = 48.693330, 26.608610  # for Kamyanets
lat_Kharkiv, lon_Kharkiv     = 49.926667, 36.278889  # for Kharkiv
lat_Kherson, lon_Kherson     = 46.73833, 32.70833    # for Kherson
lat_Kyiv, lon_Kyiv           = 50.39222, 30.53639    # for Kyív
lat_Dnipro, lon_Dnipro       = 48.36000, 35.08500    # for Dnipro
lat_Odesa, lon_Odesa         = 46.440833, 30.770278  # for Odesa
lat_Poltava, lon_Poltava     = 49.609444, 34.544722  # for Poltava
lat_Lugansk, lon_Lugansk     = 48.565556, 39.2275    # for Lugansk
# -----------------------------------------------------------------------------

########## DNIPRO ######################################################################
df_merged_Dnipro = create_df_merged(lat_station=lat_Dnipro, lon_station=lon_Dnipro,
                                    ds_20CR=ds_20CR,
                                    filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Dnipro_ta_subdaily.tsv')
                             
print('analyzing Dnipro station...')
lags_Dnipro, corrs_Dnipro = compute_lag_corr(df_merged_Dnipro, max_lag=16)
plot_lag_correlation(lags_Dnipro, corrs_Dnipro, station_name='Dnipro')

# window size of 6 month and advance of 2 month at a time for this shorter series
rolling_Dnipro = rolling_lag_corr(df_merged_Dnipro, window_size=180, step_size=60, max_lag=20)
plot_rolling_lag(rolling_Dnipro, station_name='Dnipro', expected_shift=-13)

########## KAMYANETS ######################################################################
df_merged_Kamyanets = create_df_merged(lat_station=lat_Kamyanets, lon_station=lon_Kamyanets,
                                       ds_20CR=ds_20CR,
                                       filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kamyanets_ta_subdaily.tsv')
print('analyzing Kamyanets station...')
lags_Kamyanets, corrs_Kamyanets = compute_lag_corr(df_merged_Kamyanets)
plot_lag_correlation(lags_Kamyanets, corrs_Kamyanets, station_name='Kamyanets')
# Compute rolling lag correlation for Kamyanets
rolling_Kamyanets = rolling_lag_corr(df_merged_Kamyanets)
plot_rolling_lag(rolling_Kamyanets, station_name='Kamyanets', expected_shift=-12)

########## KHARKIV ######################################################################
df_merged_Kharkiv = create_df_merged(lat_station=lat_Kharkiv, lon_station=lon_Kharkiv,
                                     ds_20CR=ds_20CR,
                                     filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kharkiv_ta_subdaily.tsv')
print('analyzing Kharkiv station...')
lags_Kharkiv, corrs_Kharkiv = compute_lag_corr(df_merged_Kharkiv)
plot_lag_correlation(lags_Kharkiv, corrs_Kharkiv, station_name='Kharkiv')
# Compute rolling lag correlation for Kharkiv
rolling_Kharkiv = rolling_lag_corr(df_merged_Kharkiv)
plot_rolling_lag(rolling_Kharkiv, station_name='Kharkiv', expected_shift=-12)

########## KHERSON ######################################################################
df_merged_Kherson = create_df_merged(lat_station=lat_Kherson, lon_station=lon_Kherson,
                                     ds_20CR=ds_20CR,
                                     filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kherson_ta_subdaily.tsv')
print('analyzing Kherson station...')
lags_Kherson, corrs_Kherson = compute_lag_corr(df_merged_Kherson)
plot_lag_correlation(lags_Kherson, corrs_Kherson, station_name='Kherson')
# Compute rolling lag correlation for Kherson
rolling_Kherson = rolling_lag_corr(df_merged_Kherson)
plot_rolling_lag(rolling_Kherson, station_name='Kherson', expected_shift=-12)

######## KYÍV #############################################################################
df_merged_Kyiv = create_df_merged(lat_station=lat_Kyiv, lon_station=lon_Kyiv,
                                  ds_20CR=ds_20CR,
                                  filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kyiv_ta_subdaily.tsv')
# Compute correlation for lags from -14 to +14 days for Kyív
print('analyzing Kyív station...')
lags_Kyiv, corrs_Kyiv = compute_lag_corr(df_merged_Kyiv)
plot_lag_correlation(lags_Kyiv, corrs_Kyiv, station_name='Kyív')
# Compute rolling lag correlation for Kyív
rolling_Kyiv = rolling_lag_corr(df_merged_Kyiv)
plot_rolling_lag(rolling_Kyiv, station_name='Kyív', expected_shift=0)

######## LUGANSK #############################################################################
df_merged_Lugansk = create_df_merged(lat_station=lat_Lugansk, lon_station=lon_Lugansk,
                                   ds_20CR=ds_20CR,
                                   filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Lugansk_ta_subdaily.tsv')
# Compute correlation for lags from -14 to +14 days for Kyív
print('analyzing Lugansk station...')
lags_Lugansk, corrs_Lugansk = compute_lag_corr(df_merged_Lugansk)
plot_lag_correlation(lags_Lugansk, corrs_Lugansk, station_name='Lugansk')
# Compute rolling lag correlation for Kyív
rolling_Lugansk = rolling_lag_corr(df_merged_Lugansk)
plot_rolling_lag(rolling_Lugansk, station_name='Lugansk', expected_shift=0)

######## ODESA #############################################################################
df_merged_Odesa = create_df_merged(lat_station=lat_Odesa, lon_station=lon_Odesa,
                                   ds_20CR=ds_20CR,
                                   filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Odesa_ta_subdaily.tsv')
# Compute correlation for lags from -14 to +14 days for Kyív
print('analyzing Odesa station...')
lags_Odesa, corrs_Odesa = compute_lag_corr(df_merged_Odesa)
plot_lag_correlation(lags_Odesa, corrs_Odesa, station_name='Odesa')
# Compute rolling lag correlation for Kyív
rolling_Odesa = rolling_lag_corr(df_merged_Odesa)
plot_rolling_lag(rolling_Odesa, station_name='Odesa', expected_shift=0)


######## POLTAVA #############################################################################
df_merged_Poltava = create_df_merged(lat_station=lat_Poltava, lon_station=lon_Poltava,
                                  ds_20CR=ds_20CR,
                                  filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Poltava_ta_subdaily.tsv')
# Compute correlation for lags from -14 to +14 days for Kyív
print('analyzing Poltava station...')
lags_Poltava, corrs_Poltava = compute_lag_corr(df_merged_Poltava)
plot_lag_correlation(lags_Poltava, corrs_Poltava, station_name='Poltava')
# Compute rolling lag correlation for Kyív
rolling_Poltava = rolling_lag_corr(df_merged_Poltava)
plot_rolling_lag(rolling_Poltava, station_name='Poltava', expected_shift=0)
