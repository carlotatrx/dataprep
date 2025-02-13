## Description of the stations

- _Barcelona_: From Stefan, he said original data and person of contact is Mariano Barriendos from Universitat de Barcelona
- _vor1763_: folders with various stations from IMPROVE, whose records start before 1763 but may continue to 19th century
- _supplementarydata_: we don't know where they are from
- _Ylitornio, Valencia_: from Stefan

## Changes applied

"Original" files (i.e. as I got them) can be found in `/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig`.

All processed stations and extra files can be found in `/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/`.

The preprocessing and controls or whatever is done with the `daily_means.ipynb` script.

### València

All changes are applied to `Valencia_concatenated.csv`, the original files are left untouched.

Some values seem to be wrongly typed, e.g. rows 2192-2194, 3288-3290 of 1804-1813.csv. I changed the month to 1 instead of 2. I changed the column "Date", which I have created myself on Excel, but not the original "Month" column. The dates which are changed are those that are duplicated, see file `Valencia_orig_obsperday_gt3.csv`. Usually the morning measurement, at 09:00h, was duplicated, so I erased it. At 1860 the dataset is quite meh so I just did until 1859. It is a weight average with $0.3*morning + 0.5*noon + 0.2*evening$

Only 36 days have less than 3 measurements. This is 0.23% of the entire dataset so I decided not to compare with ERA-5 climatology. The days with fewer measurements are stored in the file `Valencia_orig_obsperday_lt3.csv` and the weighted average is $0.5*obs_A+0.5*obs_B$.

### Ylitornio

The following dates are duplicated:

| index | Year | Month | Day | Air pressure | Temperature | Date       |
| ----- | ---- | ----- | --- | ------------ | ----------- | ---------- |
| 66    | 1800 | 12    | 6   | 25.60        | -1.0        | 1800-12-06 |
| 3147  | 1800 | 12    | 6   | 25.48        | -4.0        | 1800-12-06 |
| 3512  | 1810 | 12    | 6   | 25.40        | -10.0       | 1810-12-06 |
| 6799  | 1810 | 12    | 6   | 25.48        | -10.0       | 1810-12-06 |
| 8960  | 1826 | 11    | 5   | 25.55        | -6.7        | 1826-11-05 |
| 9325  | 1826 | 11    | 5   | NaN          | -2.0        | 1826-11-05 |

This is most likely a typo in the writing of the year. I correct this directly in the code.

### Cádiz

We have 29244-13=29231 days with values in Peter's data from IMPROVE. Now let's count how many we have from the raw. We have no duplicates in dataset, and we have 29231 non-NaNs. Yay! it's the same. It looks like `CSF-TP801-819.txt` is just a small subset of `CSF-TP786-879.txt`. No changes made, simply stored the files of temperature and pressure separately.
