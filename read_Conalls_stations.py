import pandas as pd
import glob
import os

file_path = "/scratch3/PALAEO-RA/daily_data/tmp/Files_for_Carlota_20260216/NonScratch_HistObs_Part01_CR_daily.xlsx"


def extract_pre_1900_single_file(filepath):
    print(f"Reading {filepath}...")
    
    # Read the excel file without headers so we can navigate by raw index
    df = pd.read_excel(filepath, sheet_name="Sheet1")

    # drop last two columns (they are empty)
    df = df.iloc[:, :-2]
    
    # 1. Extract the Metadata Block
    # Metadata labels are in Column 6 (index 5)
    # Metadata values start from Column 7 (index 6)
    # Metadata is in the first 21 rows (indices 0 to 20)
    meta_labels = df.iloc[0:20, 5].tolist()
    meta_values_df = df.iloc[0:20, 6:].copy()
    meta_values_df.index = meta_labels # Set the row names to 'StnNumberCR', etc.
    
    # 2. Extract the Data Block
    # Data starts at line 23 (index 22)
    data_df = df.iloc[21:].copy()
    
    data_df.rename(
        columns=dict(zip(data_df.columns[:5], ["DayIndex", "Year", "Month", "Day", "DoY"])),
        inplace=True
    )
    data_df["Year"] = pd.to_numeric(data_df["Year"], errors="coerce")
    data_df["Month"] = pd.to_numeric(data_df["Month"], errors="coerce")
    data_df["Day"] = pd.to_numeric(data_df["Day"], errors="coerce")

    data_df = data_df.dropna(subset=["Year"])
    pre_1900_df = data_df[data_df["Year"] < 1900]
    
    results = []
    
    # 3. Check each station's data column
    # Iterate over the data columns (Column 7 onwards -> index 6 onwards)
    for col_idx in range(6, df.shape[1]):
        
        # Grab this column's pre-1900 data and drop the 'NaN's
        station_data = pd.to_numeric(pre_1900_df.iloc[:, col_idx], errors='coerce')
        valid_data = station_data.dropna()
        
        # If there is valid data left, record it!
        if not valid_data.empty:
            first_idx = valid_data.index[0]
            last_idx = valid_data.index[-1]
            
            # Format the dates as YYYYMMDD
            start_date = f"{int(pre_1900_df.loc[first_idx, "Year"])}{int(pre_1900_df.loc[first_idx, "Month"]):02d}{int(pre_1900_df.loc[first_idx, "Day"]):02d}"
            end_date  = f"{int(pre_1900_df.loc[last_idx, "Year"])}{int(pre_1900_df.loc[last_idx, "Month"]):02d}{int(pre_1900_df.loc[last_idx, "Day"]):02d}"
            
            # Pull the corresponding metadata for this specific column
            station_info = meta_values_df.iloc[:,col_idx].to_dict()
            station_info['start_date'] = start_date
            station_info['end_date'] = end_date
            
            results.append(station_info)

    # 4. Build the final output table
    final_df = pd.DataFrame(results)
    
    # Reorder the columns to put your new date fields at the front
    if not final_df.empty:
        front_cols = ['start_date', 'end_date', 'StationName']
        remaining_cols = [c for c in final_df.columns if c not in front_cols]
        final_df = final_df[front_cols + remaining_cols]
        
    return final_df

# --- Execution ---

# list all excel files in the folder
folder_path = "/scratch3/PALAEO-RA/daily_data/tmp/Files_for_Carlota_20260216/"
excel_files = glob.glob(os.path.join(folder_path, "*.xlsx"))

# delete the files that end in "QCd" from the list
excel_files = [f for f in excel_files if not f.endswith("QCd.xlsx")]

# iterate over the files and extract pre-1900 stations
for file in excel_files:
    print(f"Processing file: {file}")

    # exract the number before _CR_daily.xlsx
    file_number = os.path.basename(file).split('_')[2].replace('Part', '')
    
    result_table = extract_pre_1900_single_file(file_path)
    # now delete the first part of the name, e.g.  'NonScratch_HistObs_Part09' and keep the number only

    print("Pre-1900 Stations Found:")
    print(result_table.head())

    # Optional: Save it so you can inspect the output
    result_table.to_csv(f"/scratch3/PALAEO-RA/daily_data/tmp/Files_for_Carlota_20260216/conalls_stations_file_{file_number}.csv", index=False)
