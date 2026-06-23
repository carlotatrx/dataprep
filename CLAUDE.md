# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Data preparation pipeline for the **PALAEO-RA** project: processing ~100 historical European weather station records (~1750–1875) into **SEF (Structured Environmental Format)** files for use in a Kalman-filter climate reanalysis. Scripts are in R, Python, and Bash; Jupyter notebooks are used for exploratory analysis and plotting.

## Key external data paths

- `/scratch3/PALAEO-RA/daily_data/final/<StationName>/` — output SEF `.tsv` files (one per station/variable)
- `/scratch3/PALAEO-RA/daily_data/original/<StationName>/` — raw input data
- `/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/` — WeaR and other pre-processed series
- `/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/` — intermediate preprocessed files
- `/scratch2/ccorbella/` — same filesystem as `/mnt/climcal2/scratch2/ccorbella/`

## Running scripts

No build system. Scripts are run individually:

```bash
Rscript 02_dataprep_timeseries2tsv.R
Rscript 04_dataprep_tsv_obs.R
python3 05A_dataprep_ta_obs_anomalies.py
python3 merge_wrdv2.py /scratch3/PALAEO-RA/daily_data/final/StationName
python3 merge_wrdv2.py /path/to/dir --dry-run   # preview without writing
jupyter notebook plots4paper.ipynb
```

QC pipeline (run in order):
```bash
python3 qc/01_qc_find_files.py
python3 qc/02_qc_check.py
Rscript qc/03_qc_create_files.R
python3 qc/04_qc_list_empty_files.py
Rscript qc/05_qc_replace_files.R
```

## Main pipeline

The numbered scripts define the processing order:

| Script | Language | Purpose |
|--------|----------|---------|
| `01_dataprep_stations.ipynb` | Python/nb | Station inventory and metadata setup |
| `02_dataprep_timeseries2tsv.R` | R | Convert raw station files (XLS, CSV, etc.) to SEF `.tsv` |
| `03_daymean_ERA5_adjustment.R` | R | Daily mean computation with ERA5 adjustment |
| `04_dataprep_tsv_obs.R` | R | Reads all SEF files, extracts lat/lon/variable metadata |
| `05A_dataprep_ta_obs_anomalies.py` | Python | Temperature deseasonalization vs. 20CR model |
| `05B_dataprep_p_obs_anomalies.py` | Python | Pressure anomaly computation |
| `05_find_20CR_indices.py` | Python | Find nearest 20CR grid cell for each station |
| `06_dates_span_availability_plot.ipynb` | Python/nb | Data availability visualization |

Individual station scripts in `indiv_series/` (one `.R` file per station) handle station-specific raw-to-SEF conversions and are called or run stand-alone.

## SEF file format

SEF 1.0.0 is a tab-separated format with a 12-line header followed by data rows:

```
SEF    1.0.0
ID     <station_code>
Name   <station_name>
Lat    <decimal degrees>
Lon    <decimal degrees>
Alt    <metres>
Source <citation>
Link   <url>
Vbl    <variable: ta|p|rr|dd|rh|Tx|Tn|w>
Stat   <point|mean|...>
Units  <C|hPa|mm|degree|%|m/s>
Meta   <QC and other metadata flags>
Year   Month  Day  Hour  Minute  Period  Value  Meta
...
```

QC'd files are identified by `QC software=dataresqc` in the `Meta` header line and end with `_qc.tsv`.

## Shared R utilities (`helpfun.R`)

Always `source("helpfun.R")` (or its absolute path) in R scripts. Key functions:

- `write_sef_f(Data, outpath, outfile, variable, cod, ...)` — write a data frame to SEF format
- `write_flags_f(infile, qcfile, outpath, ...)` — apply QC flags from a `dataresqc` output file onto an existing SEF
- `read_meta_nonofficial(file, parameter)` — read the 12-line SEF header into a named vector
- `outfile.name(name, var, df, subdaily, obs_name)` — construct standardized output filenames
- `dd_normalize(x)` / `dd2deg(x)` — normalize wind direction strings and convert to degrees
- `units(var)` — canonical units string for a variable code

## Variable codes

| Code | Variable | Units |
|------|----------|-------|
| `ta` | air temperature | C |
| `p` | pressure | hPa |
| `rr` | precipitation | mm |
| `dd` | wind direction | degree |
| `rh` | relative humidity | % |
| `Tx` / `Tn` | max / min temperature | C |
| `w` | wind speed | m/s |

## Merge script (`merge_wrdv2.py`)

Merges consecutive WRDv2-2 segments per station+variable into a single SEF file. Keeps the first file's data for overlap periods. Moves originals to `/scratch3/PALAEO-RA/daily_data/tmp/merge_WRDv2-2/` and logs to `changes_log.txt`. Always test with `--dry-run` first.

## Calendar and unit conversion notes

- Ukrainian data: possible Julian/Gregorian calendar ambiguity (12-day offset in 19th century); see `dataprep_ukraine.ipynb`
- Russian semi-lines (R.s.l.) to hPa: 1000 hPa = 590.60 R.s.l.
- Réaumur to Celsius: `C = Réaumur × 1.25`
- Paris inches to hPa: `P_hPa = (inches + lines/12) × 27.07 mmHg/inch × 1.3322 hPa/mmHg`

## Known data issues (tracked in `readme.md`)

- **CBT**: duplicate date entries on Feb 8 in 1804, 1808, 1812 — currently dropped, should be fixed in source
- **Ylitornio**: pressure unit uncertain (assumed Paris inches); two conflicting versions from different sources
- **València**: large fraction (~23%) of pressure values removed as outliers (>5σ from mean)
- Items marked `🔴` in `readme.md` are open questions requiring follow-up with data providers
