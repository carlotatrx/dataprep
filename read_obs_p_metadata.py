import os
import pandas as pd

filepath = '/scratch2/ccorbella/files/1807_USBstick/'

data = []

# Iterate through files in the directory
for filename in os.listdir(filepath):
    if '_p_' in filename and filename.endswith('.tsv'):
        try:
            with open(f"{filepath}{filename}", 'r') as file:
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

                # Read data portion into a DataFrame
                data_start_idx = next(i for i, line in enumerate(lines) if line.startswith("Year"))
                df = pd.read_csv(filename, sep='\t', skiprows=data_start_idx, na_values=['NA'])
                
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

        except Exception as e:
            print(f"Error processing file: {e}")

# Create DataFrame and display results
df_output = pd.DataFrame(data)
df_output.to_csv('p_obs_metadata.csv', index=False)