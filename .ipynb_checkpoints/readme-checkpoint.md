This README refers to the script `dataprep_stations.ipynb` and to the QC.

# SEF QC Workflow

How the qualityf control files in this directory fit together and the order to run them:

1. `qc_find_files.py`: Normalises filenames based on whether QC metadata in the SEF header. → Run this first to ensure header ↔ filename consistency.
2. `qc_check.py`: Checks headers of all SEFs, and produces three lists:
 - files needing QC
 - files already QC’d
 - files with unit or variable name problems
    → The important output for the next step is files_no_qc.txt.
3. `qc_files.R`: Runs dataresqc QC on all files listed in files_no_qc.txt, writes QC results, and applies the flags.
4. `qc_empty_files_list.py` After QC, lists which QC output files are empty vs non-empty.

**TO-DO**: run the qcs on step 3 and do a step 5 so that the files with emtpy qc count as qc'd.

## Description of the stations

* _Barcelona_: From Stefan, he said original data and person of contact is Mariano Barriendos from Universitat de Barcelona
* _vor1763_: folders with various stations from IMPROVE, whose records start before 1763 but may continue to 19th century
* _supplementarydata_: we don't know where they are from

In the following it's described the changes applied to original files. The output files, saved in `station_timeseries_preprocessed` directory, all have the same format, `stationname_var.csv`, where `var` can be `PRMSL` for pressure or `TMP2m` for temperature. This doesn't mean necessarily that pressures have been reduced to mean sea level or that temperature is at 2m, but I'm using these abbreviations for consistency.

"Original" files (i.e. as I got them) can be found in `/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/`.

All processed stations and extra files can be found in `/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/`.

The preprocessing and controls or whatever is done with the `daily_means.ipynb` script. If for some station one of the variables is found in two pre-processed files (i.e. in USB_stick folder and in station_timeseries_preprocessed folder), it means that the data are the same (unless otherwise specified) and any can therefore be used.

Data in the subdirectory `station_timeseries_preprocessed/not_to_use/` are there because we already have an alternative (from WeaR, KNMI or whatever) and we choose to use that one, namely found in the `1807_USBstick` folder. Similarly, the directory `1807_USBstick/not_to_use/` contains timeseries where a newer or better timeseries exists, namely in the timeseries shared by Stefan or Yuri.

Descriptions and problems in OneDrive excel `stations_1806-1850`, including lat, lon, altitude, source of series from `station_timeseries_orig`.

### Bologna

File `TG_SOUID100862.txt` from folder `ECA_non-blended_custom-1`, where metadata can also be found. These are specified as:

```
17-21 TG   : mean temperature in 0.1 °C

23-27 Q_TG : Quality code for TG (0='valid'; 1='suspect'; 9='missing')
```

All cells with data have a quailty code of 0 (valid). I converted the `TG` column to °C units by multiplying by 0.1.

`🔴 To ask for: `Stefan's Bologna starts in 1814, it's from ECA&D (`1807raw_andmore/ECA_non-blendd_custom-1/TG_SOUID100862`) but online website doesn't exist anymore. Peter's Bologna is from IMROVE and has only 1807. Where are the rest from Peter/Noemi?

Noemi has some Bologna data in `Noemi_scratch3/laki1783/data/01_stations/bologna/Camuffo_Camuffo_Bologna_17150101-18151231_ta.tsv` which is very different from Stefan's.

In the history of 1796-1812: First observations at the Meteorological and Astronomical Observatory, later discarded as deemed of a too low quality

### Cádiz

We have 29244-13=29231 days with values in Peter's data from IMPROVE. Now let's count how many we have from the raw. We have no duplicates in dataset, and we have 29231 non-NaNs. Yay! it's the same. It looks like `CSF-TP801-819.txt` is just a small subset of `CSF-TP786-879.txt`, which we can check with `np.sum(cadiz2['TMP2m']!=cadiz1[:3323]['TMP2m'])`. As they are the same, we are only using the longer dataset. No changes made otherwise, simply stored the files of temperature and pressure separately. I'll consequently use the IMPROVE file from Cadiz TMP2m as it should be the same.

### CBT

`🔴` there are duplicate values I need to check for in the original file. For now, I am just dropping them in `04_dataprep_tsv_obs.R`. Fix them in the original file and remove the dropping of the second value, which is probably a typo of the year or the month improperly written.

```
⚠️  Duplicates found in: Central Belgium Temperature 
# A tibble: 6 × 6
# Groups:   Year, Month, Day [3]
  row_num  Year Month   Day Value Meta 
    <int><int><int><int><dbl><chr>
13642  1804     2     8 -2""   
23643  1804     2     8  2.27 ""   
35103  1808     2     8  2.4  ""   
45104  1808     2     8  1.7  ""   
56564  1812     2     8  1.63 ""   
66565  1812     2     8  6.37 ""   
```

Replace this 

```
      df_val <- x %>%
        group_by(Year, Month, Day) %>%
        mutate(n = n()) %>%
        filter(row_number() == 1) %>%
        ungroup() %>%
        select(Year, Month, Day, Value)
  
      dropped_vals <- x %>%
        group_by(Year, Month, Day) %>%
        filter(n() > 1 & row_number() > 1) %>%
        select(Year, Month, Day, Value)
  
      if (nrow(dropped_vals) > 0) {
        cat("\n🗑️  Dropping repeated values in:", colname, "\n")
        print(dropped_vals)
      }
```

with just `df_val <- x[, c("Year", "Month", "Day", "Value")]`.

### CET

Extracted from the [website](https://www.metoffice.gov.uk/hadobs/hadcet/), data from there are the same as 1807 selected data from Peter.

### London

From the same link where I extracted Paris, I got London and to take the mean pressure I simply took the arithmetic average of the observations of the day (sometimes 2, sometimes 3). But! The values from `London_p_daily.tsv` are very similar, but not exactly the same, as simply taking the arithmetic daily means from `London_11_17870101-18221231_mslp.tsv`. This version from Stefan is newer, so I should use the newer version.

### Milan

We already had temperature from `PALAEO-RA_IMPROVE_Milan_17630101-18621231_ta`, which is the same as the third column of the raw file, and now I add pressure.

### Padova

Ja teníem pressure i temperature fins a 1809, ara s'ha d'allargar i per no tenir 2 series de la mateixa station les concatenejo i en faig una de sola, que la utilitzaré per reemplaçar la de la carpeta USB stick Peter.

### Paris

The data from `Paris_p_daily.tsv` in Peter's USB is almost the same but not exactly to the data from `Paris_4_17850101-1872061_mslp.tsv`, downloaded from  (LDL)[[https://figshare.com/articles/dataset/Sub-daily_sea-level_pressure_series_for_London_GB_Paris_FR_and_De_Bilt_NL_/24242302/3?file=45614718]()]. This link has only pressure data. We use the LDL data, not Peter's data.

`🔴 To ask for:` Moreover, the data in the column _T°C_ on the tab _Tjour_ in the file `paris_Daily_Updated_meteo_2024_127_33_supp.xlsx` (from Stefan) is different than the data from `Paris_ta_noon.tsv` from Peter's stick. If anything, it seems to be closer to the _TX_ column, even if not exactly the same. I'll use Stefan's file. The _T°C_ is $(TX+TN)/2$ with 0.5 precision. **WHAT ARE THE COORDINATES?**

### Prague

Same as Bologna. Only the following has a suspect code:

| 30641 | 1858-11-23 | -150.0 | 1 |
| ----- | ---------- | ------ | - |

Also converted to hPa through multiplication by 0.1. Then my data are the same as Peter's folder. I can then use Peter's folder.

### Stockholm

ATTENTION!! There's another file with explanations, `stockholm_SLP_1756_2012_hPa_hom.txt`, from `1807raw_andmore/database-web-stockholm-historical-weather-observations-1.0/database-web-stockholm-historical-weather-observations-1.0/air_pressure/`.

### Uppsala

Same as Padova, without need for concatenating files.

### Ukraine

to prep for this data, file is `dataprep_ukraine.ipynb`.  Yuri says:

> The Ukrainian data should also be in the folder `/scratch3/PALAEO-RA/DataRescue/Incoming/Ukraine`. The reason we didn't use them is that it is not clear which calendar is used. Russia adopted the Gregorian calendar very late (I think in 1917), and Ukraine has always been a bit of a hybrid between Russia and Western Europe… Therefore some stations used Julian calendar, others Gregorian, or a mixture of the two. There are 12 days of offset between Julian and Gregorian calendar  in the 19$^{th}$ century. It is certainly possible to figure out which calendar is used where, but it implies some work. If could use for example correlation, but the main  difficult is that a station  may have started using Julian and then switched to Gregorian at some unknown time (in that case you should see a gap of 12 days). I don't know if the opposite also happens, less likely but not impossible (then you would see duplicated dates).

Ukrainian data conversion of Russian semi lines (R.s.l.) to hPa:
 1 R.s.l. = 1.27 mm. That is, 1,000 hPa = 1,000 mbar = 750.06 mmHg = 590.60 R.s.l.

Celsius = Réaumur * 1.25

Station IDs were taken from the supplement of teh paper _Ukrainian early (pre-1850) historical weather observations_, last page.

**Dnipro**

> the measurements were performed at specified times three times per day. Air temperature was recorded in Reaumur degrees. Atmospheric pressure was firstly (1833–1838) recorded in inches, but later (1839–1842, 1850) R.s.l. were used as units for this.

No days with less than 3 observations. Daily means is weighted average.

**Kherson**

I'm only looking until 1825 because afterwards there are different measurements in various places.

> Two values of air temperature during 1808-February 1825 were recorded every time of measurement, and both of them were digitized. These two values of temperature were marked in the paper source as ‘on S’ and ‘on N’, probably meaning ‘on South’ and ‘on North’ walls. The differences between these two temperatures are very small. Such temperature time series seem to be something like ‘parallel’ measurements. However, no detailed information was found regarding these slightly different temperature values.

Look at `dataprep_ukraine.ipynb` for explanations.

1809-02-14, for the three measurements, has values of 90 in (3000 hPa), so we discard this day. Since we can't trust the pressure, we discard temperature values as well. 1809-04-03, we assume a typo in the input of the pressure: it should be 29.82 instead of 39.82.

Delete also all days where pressure is < 950hPa. According to Duncan, these are very low pressure values and he assumes that cyclones would already be quite weakened by the time they reach the Ukraine.

**Kyív**

Out of 1128 days with measurements, 66 days (a 1.83%) had less than 3 observations. We simply discard these days.

> two periods (1820–1822 and 1826–1836) are missing. During the whole period, observations were conducted three times per day, though the exact time of the day is reported only since 1837. In earlier records, only the fraction of the day (morning, midday and evening) was specified. During the period of 1812–1825, only temperature records were reported (in Reaumur degree).

### Torino

If I use $t_{daily}=(t_{min}+t_{max})/2$ it's a bit different from the WeaR data `WeaR_TOR_17530101-18621231_ta_hom.tsv`, because of the homogenization component. Noemi suggests I use the original, not the homogenized. Stefan says I should use the homogenized. So I use the homogenized, `/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/WeaR_TOR_17530101-18631231_ta_hom.tsv`.

### València

All changes are applied to `Valencia_concatenated.csv`, the original files are left untouched.

The data are described in Domínguez Castro, Fernando, et al. "Early Spanish meteorological records (1780-1850)." (2014). Each daily observation includes five variables: air temperature (Reaumur scale), atmospheric pressure (in inches and lines), humidity (with a hygrometer made by the observer and expressed in 24 degrees of humidity), etc. Bar(p) are the inches, Paris inch + Bar(l) are lines (1 line=1/12 of an inch).  We assume that the air temperature given is Reaumur scale and the pressure, in Paris inches.

$$
P_{hPa}=\mathrm{(Bar(p) + Bar(l)/12) * 27.07 \frac{mmHg}{inch}*1.3322\frac{hPa}{mmHg}}
$$

Some values seem to be wrongly typed, e.g. rows 2192-2194, 3288-3290 of 1804-1813.csv. I changed the month to 1 instead of 2. 1806-03-25, 1830-09-07 and 1841-10-18 are also wrongly typed.

I changed the column "Date", which I have created myself on Excel, but not the original "Month" column. The dates which are changed are those that are duplicated, see file `Valencia_orig_obsperday_gt3.csv`. Usually the morning measurement, at 09:00h, was duplicated, so I erased it. At 1860 the dataset is quite meh so I just did until 1859.·Only 36 days have less than 3 measurements. This is 0.20% of the entire dataset. The days with fewer measurements are stored in the file `Valencia_orig_obsperday_lt3.csv` and they are removed from the dataset.

Some València temperatures are a bit suspicious, 14°C difference between day and night in the summer (*e.g.* 1830-08-23). I remove pressure values that are 5$$ away from the mean. This removes 13749 (23%) of the values in the period 1806-1850!!

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

The pressure data are, presumably, in inches of mercury. I converted to hPa using $1013.25 * ( P_{inHg} / 29.92)$. It might not be the proper conversion (I couldn't find what the units were, but since I'll assimilate anomalies, it won't matter too much anyway).

`🔴 To ask for: ` `Ylitornio_ta_daily.tsv` from Peter is very different from Stefan. 

`🔴` Changing now pressure as if it were Paris inches, 1 Paris inch = 27.07mmHg. 1mmHg = 1.3332239hPa. But we'll need to compare w 20CR to see what an appropriate comparison would be.

`🔴` WHAT ARE THE COORDS OF STEFAN'S YLI?? This is how different they are:

![Yli_fig](https://github.com/carlotatrx/KF_assimilation/blob/main/dataprep/image/Yli_diff.png)

The only additional info we have is that in Peter's data it says "Stat: mean", and Stefan's data comes from the paper: 'Reference: Helama, S., Holopainen, J., Timonen, M., Ogurtsov, M. G., Lindholm, M., Meriläinen, J., Eronen, M. (2004). Comparison of living-tree and subfossil ringwidths with summer temperatures from 18th, 19th and 20th centuries in Northern Finland. Dendrochronologia 21/3, 147 - 154.'
They say they homogenize daily temperatures. "Homogenization was performed by method of Moberg and Bergström (1997). This approach accounts for changing measurement time of daily observations and is expected to produce much more homogeneous time series of monthly temperatures than a simple average of observations (Holopainen, Vesajoki 2001)."

I study how different they are in `dataprep_stations.ipynb`.

### Zwanenburg

Did arithmetic daily means of pressure, same issue as London and Paris above, as same source. How were the temperatures from KNMI extracted?
