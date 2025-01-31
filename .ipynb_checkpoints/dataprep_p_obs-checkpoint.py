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
import os
import xarray as xr
from scipy.spatial import KDTree # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

######################################################################################################################################
################ PART 1: extract station metadata from TSV files #####################################################################
######################################################################################################################################

# Directory containing the files
filepath = '/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/'

# Initialize list to store extracted data
data = []
combined_data = pd.DataFrame()

# Iterate through files in the directory
for filename in os.listdir(filepath):
    if '_ta' in filename and filename.endswith('.tsv'):
        try:
            with open(f'{filepath}{filename}', 'r') as file:
                lines = file.readlines()
                
                # Extract metadata
                metadata = {}
                for line in lines:
                    if line.startswith("SEF"):
                        continue
                    if line.startswith("Year"):
                        break
                    if '\t' in line:
                        parts = line.strip().split('\t', 1)
                        if len(parts) == 2:
                            key, value = parts
                            metadata[key] = value

                # Extract station name
                station_name = metadata.get('Name', 'Unknown')

                # Read data portion into a DataFrame
                data_start_idx = next(i for i, line in enumerate(lines) if line.startswith("Year"))
                df = pd.read_csv(f'{filepath}{filename}', sep='\t', skiprows=data_start_idx, na_values=['NA'])
                
                # Find smallest and largest dates
                df['Date'] = pd.to_datetime(df[['Year', 'Month', 'Day']].dropna(how='any'), errors='coerce')
                min_date = df['Date'].min()
                max_date = df['Date'].max()
                
                # Append extracted data
                data.append({
                    'Name': metadata.get('Name', ''),
                    'Lat': metadata.get('Lat', ''),
                    'Lon': metadata.get('Lon', ''),
                    'Alt': metadata.get('Alt', ''),
                    'Vbl': metadata.get('Vbl', ''),
                    'Stat': metadata.get('Stat', ''),
                    'Units': metadata.get('Units', ''),
                    'Min Date': min_date.strftime('%Y-%m-%d') if pd.notnull(min_date) else '',
                    'Max Date': max_date.strftime('%Y-%m-%d') if pd.notnull(max_date) else ''
                })
                
                # Select relevant columns
                df['Station'] = station_name
                df = df[['Date', 'Value', 'Station']]

                # Append to combined data
                combined_data = pd.concat([combined_data, df], ignore_index=True)

        except Exception as e:
            print(f"Error processing file: {e}")

# save metadata to a CSV file
df_output = pd.DataFrame(data)
df_output.to_csv(f'{filepath}t2m_obsvalidation_metadata.csv', index=False)

# Pivot the combined data to make stations as columns
pivot_table = combined_data.pivot(index='Date', columns='Station', values='Value')
pivot_table.to_csv(f'{filepath}t2m_obsvalidation_data.csv')

######################################################################################################################################
################ PART 2: calculate anomalies at each station #########################################################################
######################################################################################################################################

df = pd.read_csv(f'{directory}p_obs_data.csv', index_col='Date', parse_dates=True)
df = df.drop('NaT')
df =  df.loc['1806-01-01':'1821-12-31']

anomals = df - df.mean()
anomals.to_csv('p_obs_anomalies.csv')

######################################################################################################################################
################ PART 3: find cells in 20CR where stations are located ###############################################################
######################################################################################################################################

# open downloaded dataset from 20CRv3
ds_20CR = xr.open_dataset('/home/ccorbella/scratch2_symboliclink/files/20CRv3_ensembles/PRMSL.Europe.allmems/PRMSL.Europe.allmems.allyears.nc')
ds_lats = ds_20CR['lat'].values
ds_lons = ds_20CR['lon'].values

station_meta = pd.read_csv('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/p_obs_metadata.csv')
station_lats = np.array(station_meta['Lat']).astype(float)
station_lons = np.array(station_meta['Lon']).astype(float)

# Extract grid points from 20CR dataset
grid_lons, grid_lats = np.meshgrid(ds_lons, ds_lats)

# Flatten the grid
grid_points = np.array([grid_lats.flatten(), grid_lons.flatten()]).T

# Build a KDTree for quick nearest-neighbor lookup
tree = KDTree(grid_points)

lat_index, lon_index, station_name = [], [], []

for i in range(station_meta.shape[0]-1, -1, -1): # reverse order
    # Find the nearest grid points for each station
    distances, station_index = tree.query(np.array([station_lats[i], station_lons[i]]).T)

    # Convert the flattened index back to (lat, lon) indices
    lat_index.append(station_index // len(ds_lons))
    lon_index.append(station_index % len(ds_lons))
    station_name.append(station_meta['Name'][i])

np.savetxt('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/lat_index_p.txt', lat_index)
np.savetxt('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/lon_index_p.txt', lon_index)

######################################################################################################################################
################ PART 4: visualize stations ###############################################################
######################################################################################################################################

# find corresponding cell in 20CR
station_lats_20CR, station_lons_20CR = [], []
for i in range(len(lat_index)):
    station_lats_20CR.append(ds_lats[lat_index[i]])
    station_lons_20CR.append(ds_lons[lon_index[i]])

plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-30, 50, 30, 75], crs=ccrs.PlateCarree())  # Focus on Europe

# Add European contours
ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.add_feature(cfeature.LAND, edgecolor='black')

# Plot grid lines based on ds_lons and ds_lats
for lon in ds_lons:
    plt.plot([lon] * len(ds_lats), ds_lats, color='lightgray', linewidth=0.5, transform=ccrs.PlateCarree())
for lat in ds_lats:
    plt.plot(ds_lons, [lat] * len(ds_lons), color='lightgray', linewidth=0.5, transform=ccrs.PlateCarree())

# Add grid lines
gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
gl.xlabels_top = False
gl.ylabels_right = False

# Add the station points
plt.scatter(station_lons, station_lats, color='red', label='Pressure Stations', zorder=5, transform=ccrs.PlateCarree())
plt.scatter(station_lons_20CR, station_lats_20CR, color='blue', label='20CRv3 closest grid')
plt.title('Pressure Stations and Their Closest 20CRv3 Points')
plt.legend()

# Show the plot
plt.savefig('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/station_plot_p.png')
plt.show()    
