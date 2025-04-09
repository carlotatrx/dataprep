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
import xarray as xr
from scipy.spatial import KDTree
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#%%
def find_closest_gridpoints(var_short: str, var_long: str, plot: bool = False):
    """
    Find the closest grid points in a 20CR dataset for each observation station.

    Parameters:
    - var_short: Short variable code (e.g., 'ta' or 'p')
    - var_long: Long variable name as used in the NetCDF files (e.g., 'TMP2m' or 'PRMSL')
    - dataset_path: Path to the NetCDF file
    - station_csv_path: Path to CSV with station metadata (columns: 'Name', 'lat', 'lon')
    - output_path: Where to save the resulting station-to-gridpoint index file
    """
    ds = xr.open_dataset(f'/home/ccorbella/scratch2_symboliclink/files/20CRv3_ensembles/{var_long}.Europe.allmems/{var_long}.Europe.allmems.allyears.nc')
    ds_lats = ds['lat'].values
    ds_lons = ds['lon'].values

    station_meta = pd.read_csv(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/{var_short}_latlon.csv')
    station_lats = station_meta['lat'].astype(float).values
    station_lons = station_meta['lon'].astype(float).values

    # Build grid of lat-lon points
    grid_lons, grid_lats = np.meshgrid(ds_lons, ds_lats)
    grid_points = np.array([grid_lons.flatten(), grid_lats.flatten()]).T
    tree = KDTree(grid_points)

    lat_index, lon_index, station_name = [], [], []

    for i in range(station_meta.shape[0]):
        distances, station_index = tree.query(np.array([station_lons[i], station_lats[i]]))
        lat_idx = station_index // len(ds_lons)
        lon_idx = station_index % len(ds_lons)
        lat_index.append(lat_idx)
        lon_index.append(lon_idx)
        station_name.append(station_meta['Name'][i])

        print(f"{station_meta['Name'][i]}: Station lat/lon = ({station_lats[i]:.4f}, {station_lons[i]:.4f}) -> "
              f"Grid lat/lon = ({ds_lats[lat_idx]:.4f}, {ds_lons[lon_idx]:.4f})")

    stations_df = pd.DataFrame({'station': station_name, 'lat': lat_index, 'lon': lon_index})
    stations_df = stations_df.sort_values(by='station').reset_index(drop=True)
    stations_df.to_csv(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/{var_short}_latlon_index.csv', index=False)
    print(f"\n✅ Saved gridpoint indices of {var_short}.\n")



    if plot:
        # Reconstruct matched grid lat/lon values
        station_lats_20CR = [ds_lats[i] for i in lat_index]
        station_lons_20CR = [ds_lons[i] for i in lon_index]

        # Make plot
        plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([-30, 50, 30, 75], crs=ccrs.PlateCarree())

        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, edgecolor='black')

        # Grid lines
        for lon in ds_lons:
            plt.plot([lon] * len(ds_lats), ds_lats, color='lightgray', linewidth=0.5, transform=ccrs.PlateCarree())
        for lat in ds_lats:
            plt.plot(ds_lons, [lat] * len(ds_lons), color='lightgray', linewidth=0.5, transform=ccrs.PlateCarree())

        gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
        gl.top_labels = False
        gl.right_labels = False

        # Plot stations and matched grid cells
        plt.scatter(station_lons, station_lats, color='red', label=f'{var_short.upper()} Stations', transform=ccrs.PlateCarree())
        plt.scatter(station_lons_20CR, station_lats_20CR, color='blue', marker='P', label='20CRv3 closest grid', transform=ccrs.PlateCarree())
        plt.title(f'{var_short.upper()} Validation Stations and Their Closest 20CRv3 Points')
        plt.legend()

        plot_path = f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/obs_locationplot.png'
        plt.savefig(plot_path)
        plt.show()
        print(f"\n 🗺️ Map saved to: {plot_path}\n")

def extract_20CR_atlocs(var_short: str, var_long: str) -> pd.DataFrame:
    """
    Extracts the ensemble mean timeseries at each station's grid location.

    Parameters:
    - index_csv: Path to CSV containing station, lat index, lon index
    - dataset_path: Path to NetCDF dataset
    - var_name: Variable name in NetCDF (e.g., 'TMP2m', 'PRMSL')

    Returns:
    - DataFrame with dates as index and station names as columns
    """
    index_df = pd.read_csv(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/{var_short}_latlon_index.csv')
    ds_20CR = xr.open_dataset(f'/home/ccorbella/scratch2_symboliclink/files/20CRv3_ensembles/{var_long}.Europe.allmems/{var_long}.Europe.allmems.allyears.nc')
    
    time = ds_20CR['time'].values
    # station_names = index_df['station'].tolist()
    
    data_dict = {}
    
    for _, row in index_df.iterrows():
        name = row['station']
        lat_idx = int(row['lat'])
        lon_idx = int(row['lon'])
        
        # Extract the ensemble data and take mean over members
        ens_mean = ds_20CR[var_long].isel(lat=lat_idx, lon=lon_idx).mean(dim='ensemble_member').values
        data_dict[name] = ens_mean
        print(f"✓ Extracted ensemble mean for {name} at lat_idx={lat_idx}, lon_idx={lon_idx}")

    df_model = pd.DataFrame(data_dict, index=pd.to_datetime(time))
    return df_model

def make_date_col(df):
    df['Date'] = pd.to_datetime(df[['Year', 'Month', 'Day']], errors='coerce')
    return df

def plot_station_timeseries(df_obs: pd.DataFrame, df_model: pd.DataFrame, var_short: str):
    """
    Plot observed vs model ensemble mean time series for each station.

    Parameters:
    - df_obs: DataFrame with observed values (index = time, columns = station)
    - df_model: DataFrame with model ensemble mean (same format)
    - var_short: Variable short name, e.g. 'ta' or 'p' (used in file names and plot titles)
    """
    # Ensure the index is datetime
    make_date_col(df_obs)
    df_obs['Date'] = pd.to_datetime(df_obs['Date'], errors='coerce')
    df_obs.set_index('Date', inplace=True)

    # drop Year, Month, Day columns from df_obs
    df_obs = df_obs.drop(columns=['Year', 'Month', 'Day'])

    # fix date in model
    df_model['Date'] = pd.to_datetime(df_model['Unnamed: 0'])
    df_model.set_index('Date', inplace=True)
    df_model = df_model.drop(columns=['Unnamed: 0'])
    print('model',df_model.columns)
    print('obs',df_obs.columns)
    stations = df_obs.columns.intersection(df_model.columns)

    # adjustment
    if var_short == 'ta':
        df_model = df_model - 273.15
    elif var_short == 'p':
        df_model = df_model / 100
    
    for station in stations:
        plt.figure(figsize=(12, 4))
        plt.plot(df_obs.index, df_obs[station], label='Observation', color='dimgray', linewidth=1)
        plt.plot(df_model.index, df_model[station], label='20CRv3 Ensemble Mean', color='mediumvioletred', linewidth=1)
        title_var = 'Temperature' if var_short == 'ta' else 'Pressure'
        plt.title(f"{title_var} at {station}, Observed vs 20CRv3")
        ylabel = 'temperature (°C)' if var_short == 'ta' else 'pressure (hPa)'
        plt.ylabel(ylabel)
        plt.legend()
        plt.tight_layout()

        outfile = f"/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/{var_short}_{station}.png"
        plt.savefig(outfile)
        plt.close()
        print(f"📈 Saved: {station}")
    print(df_model.shape, df_obs.shape)

if __name__ == '__main__':
    find_closest_gridpoints(
        var_short='ta',
        var_long='TMP2m',
        plot=True
    )
    find_closest_gridpoints(
        var_short='p',
        var_long='PRMSL'
    )
    df_ta_20cr = extract_20CR_atlocs(var_short='ta', var_long='TMP2m')
    df_p_20cr = extract_20CR_atlocs(var_short='p', var_long='PRMSL')
    df_ta_20cr.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_20cr.csv')
    df_p_20cr.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_20cr.csv')
    print("✅ 20CR data extracted and saved.")

    plot_station_timeseries(
        df_obs=pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs.csv'),
        df_model=pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_20cr.csv'),
        var_short='ta'
    )

    plot_station_timeseries(
        df_obs=pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs.csv'),
        df_model=pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_20cr.csv'),
        var_short='p'
    )

# #%%#####################################################################################################################################
# ################ PART 3: find cells in 20CR where stations are located ###############################################################
# ######################################################################################################################################

# # open downloaded dataset from 20CRv3
# ds_20CR = xr.open_dataset('/home/ccorbella/scratch2_symboliclink/files/20CRv3_ensembles/TMP2m.Europe.allmems/TMP2m.Europe.allmems.allyears.nc')
# ds_lats = ds_20CR['lat'].values
# ds_lons = ds_20CR['lon'].values

# station_meta = pd.read_csv(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_latlon.csv')
# station_lats = np.array(station_meta['lat']).astype(float)
# station_lons = np.array(station_meta['lon']).astype(float)

# # Extract grid points from 20CR dataset
# grid_lons, grid_lats = np.meshgrid(ds_lons, ds_lats)

# # Flatten the grid
# grid_points = np.array([grid_lons.flatten(), grid_lats.flatten()]).T

# # Build a KDTree for quick nearest-neighbor lookup
# tree = KDTree(grid_points)

# lat_index, lon_index, station_name = [], [], []

# for i in range(station_meta.shape[0]-1, -1, -1): # reverse order
#     # Find the nearest grid points for each station
#     distances, station_index = tree.query(np.array([station_lons[i], station_lats[i]]).T)

#     # Convert the flattened index back to (lat, lon) indices
#     lat_index.append(station_index // len(ds_lons))
#     lon_index.append(station_index % len(ds_lons))
#     station_name.append(station_meta['Name'][i])
#     print(f"{station_meta['Name'][i]}: Station lat/lon = ({station_lats[i]}, {station_lons[i]}), "
#           f"Grid lat/lon = ({ds_lats[lat_index[-1]]}, {ds_lons[lon_index[-1]]})")


# stations_df = pd.DataFrame({'station' :station_name,'lat':lat_index, 'lon':lon_index})
# stations_df = stations_df.sort_values(by='station').reset_index(drop=True)
# stations_df.to_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_latlon_index.txt')


# ######################################################################################################################################
# ################ PART 4: visualize stations ###############################################################
# ######################################################################################################################################

# # find corresponding cell in 20CR
# station_lats_20CR, station_lons_20CR = [], []
# for i in range(len(lat_index)):
#     station_lats_20CR.append(ds_lats[lat_index[i]])
#     station_lons_20CR.append(ds_lons[lon_index[i]])

# plt.figure(figsize=(10, 8))
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent([-30, 50, 30, 75], crs=ccrs.PlateCarree())  # Focus on Europe

# # Add European contours
# ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
# ax.add_feature(cfeature.BORDERS, linewidth=0.5)
# ax.add_feature(cfeature.LAND, edgecolor='black')

# # Plot grid lines based on ds_lons and ds_lats
# for lon in ds_lons:
#     plt.plot([lon] * len(ds_lats), ds_lats, color='lightgray', linewidth=0.5, transform=ccrs.PlateCarree())
# for lat in ds_lats:
#     plt.plot(ds_lons, [lat] * len(ds_lons), color='lightgray', linewidth=0.5, transform=ccrs.PlateCarree())

# # Add grid lines
# gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
# gl.xlabels_top = False
# gl.ylabels_right = False

# # Add the station points
# plt.scatter(station_lons, station_lats, color='red', label=f'{VAR} Stations', transform=ccrs.PlateCarree())
# plt.scatter(station_lons_20CR, station_lats_20CR, color='blue', marker='P', label='20CRv3 closest grid')
# plt.title(f'{VAR} Validation Stations and Their Closest 20CRv3 Points')
# plt.legend()

# # Show the plot
# plt.savefig(f'/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/station_plot_validation_{VAR}.png')
# plt.show()    

# # %%
