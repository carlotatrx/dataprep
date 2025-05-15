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

# -----------------------------------------------------------------------------
# Core functions
# -----------------------------------------------------------------------------

def create_df_merged(lat_station, lon_station, ds_20CR, filename):
    """
    Loads station data, extracts the corresponding 20CR grid cell,
    and merges the two into a single DataFrame.
    """
    df = pd.read_csv(filename, skiprows=12, sep='\t')
    make_date_col(df)
    df_short = df[['Date', 'Value']].copy()
    df_short['Date'] = pd.to_datetime(df_short['Date']).dt.date

    assert df_short['Date'].nunique() == df_short.shape[0], "repeated dates"

    # Build KDTree to find nearest grid point
    ds_lats = ds_20CR['lat'].values
    ds_lons = ds_20CR['lon'].values
    grid_lons, grid_lats = np.meshgrid(ds_lons, ds_lats)
    grid_points = np.array([grid_lats.flatten(), grid_lons.flatten()]).T
    tree = KDTree(grid_points)
    _, station_index = tree.query(np.array([lat_station, lon_station]))

    lat_index = station_index // len(ds_lons)
    lon_index = station_index % len(ds_lons)

    ds_ta = ds_20CR['TMP2m'][:, :, lat_index, lon_index].mean(dim='ensemble_member')
    df_20CR = ds_ta.to_dataframe().reset_index()
    df_20CR = df_20CR.rename(columns={'time': 'Date', 'TMP2m': 'Value'})
    df_20CR = df_20CR.drop(columns=['lat', 'lon'])
    df_20CR['Date'] = pd.to_datetime(df_20CR['Date']).dt.date

    df_merged = pd.merge(df_20CR, df_short, on='Date', how='inner')
    df_merged = df_merged.rename(columns={'Value_x': 'model', 'Value_y': 'obs'})
    df_merged['model'] -= 273.15  # Convert from Kelvin to Celsius
    return df_merged


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


        # Apply weights (only keep times that exist in weights)
        default_weights_ToD = {10: 0.3, 14:0.5, 22:0.2}

        # Compute weighted daily means
        df["weight"] = df["Hour"].map(default_weights_ToD)
        df = df.groupby("Date").apply(
            lambda group: np.average(group["Value"], weights=group["weight"]),
            include_groups=False
        ).reset_index(name="obs")
    else:
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
    df_20CR["Date"] = pd.to_datetime(df_20CR["Date"]).dt.date

    # Merge on date
    df["Date"] = df["Date"].dt.date
    df_merged = pd.merge(df_20CR, df, on="Date", how="inner")
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


def rolling_lag_corr(df, window_size=365*3, step_size=90, max_lag=14):
    """
    Computes Pearson correlation for rolling time windows at different lags.
    """
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
                         file_savefig='/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image'):
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
                     file_savefig='/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image'):
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
lat_Kherson, lon_Kherson = 46.73833, 32.70833 # for Kherson
lat_Kyiv, lon_Kyiv       = 50.39222, 30.53639 # for Kyív
lat_Dnipro, lon_Dnipro   = 48.36000, 35.08500 # for Dnipro

df_merged_Kyiv = create_df_merged(lat_station=lat_Kyiv, lon_station=lon_Kyiv,
                                  ds_20CR=ds_20CR,
                                  filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kyiv_ta_daily.tsv')

df_merged_Kherson = create_df_merged(lat_station=lat_Kherson, lon_station=lon_Kherson,
                                     ds_20CR=ds_20CR,
                                     filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kherson_ta_daily_TonS.tsv')

df_merged_Dnipro = create_df_merged(lat_station=lat_Dnipro, lon_station=lon_Dnipro,
                                     ds_20CR=ds_20CR,
                                     filename='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Dnipro_ta_subdaily.tsv')
                             
######## KYÍV ##########
# Compute correlation for lags from -14 to +14 days for Kyív
lags_Kyiv, corrs_Kyiv = compute_lag_corr(df_merged_Kyiv)
plot_lag_correlation(lags_Kyiv, corrs_Kyiv, station_name='Kyív')
# Compute rolling lag correlation for Kyív
rolling_Kyiv = rolling_lag_corr(df_merged_Kyiv)
plot_rolling_lag(rolling_Kyiv, station_name='Kyív', expected_shift=0)

########## KHERSON ##########
lags_Kherson, corrs_Kherson = compute_lag_corr(df_merged_Kherson)
plot_lag_correlation(lags_Kherson, corrs_Kherson, station_name='Kherson')
# Compute rolling lag correlation for Kherson
rolling_Kherson = rolling_lag_corr(df_merged_Kherson)
plot_rolling_lag(rolling_Kherson, station_name='Kherson', expected_shift=-12)

########## DNIPRO ##########
lags_Dnipro, corrs_Dnipro = compute_lag_corr(df_merged_Dnipro, max_lag=16)
plot_lag_correlation(lags_Dnipro, corrs_Dnipro, station_name='Dnipro')
rolling_Dnipro = rolling_lag_corr(df_merged_Dnipro)
plot_rolling_lag(rolling_Dnipro, station_name='Dnipro', expected_shift=-13)