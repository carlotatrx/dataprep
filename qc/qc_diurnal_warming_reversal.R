## Diurnal-warming-reversal QC check for sub-daily air temperature (ta)
##
## dataresqc's wmo_time_consistency test only compares an observation to the
## single PREVIOUS one, with a tolerance that scales with the hour gap, and
## skips any pair more than 12h apart. It also has no notion of expected
## diurnal shape, so it cannot catch a same-day reading that rises when it
## should be falling.
##
## This check looks specifically for that signature: a midday reading
## (10:00-13:59) followed, later the same day, by an evening reading
## (>=15:00) that is warmer by more than half of dataresqc's own ta
## tolerance for that time gap. Physically, temperature should be flat or
## falling from midday to evening in a typical single-peaked diurnal cycle,
## so a rise of this size is a useful review signal (though not necessarily
## an error - e.g. a genuine multi-day warm/foehn spell can also trigger it,
## which is why candidates are reviewed by a human before being applied).
##
## An earlier version of this script compared each reading symmetrically to
## both its previous AND next neighbour ("symmetric_spike"); that produced
## ~1500 false positives on this file alone because it could not distinguish
## a real error from ordinary day/night temperature swing. This directional,
## time-window-restricted version replaced it.
##
## Two modes:
##  - "report": read-only, writes a candidate list for human review
##  - "apply" : splices qc=diurnal_warming_reversal into the Meta field of
##              the rows listed in the (reviewed) candidate file, in place

library(dataresqc)

MODE <- "apply"  # "report" or "apply"

inpath <- "/scratch3/PALAEO-RA/daily_data/SEF/Turin/"
files <- list.files(inpath, pattern="_ta_subdaily_qc\\.tsv$", full.names=TRUE)

outpath <- "/scratch3/PALAEO-RA/DataRescue/Scripts/"
candidate_file <- paste0(outpath, "qc_diurnal_warming_reversal_candidates.txt")

MIDDAY_HOUR_MIN <- 10
MIDDAY_HOUR_MAX <- 14   # exclusive upper bound: [10,14)
EVENING_HOUR_MIN <- 15  # inclusive lower bound: [15,24)
HALF <- 0.5             # evening rise must exceed HALF * dataresqc's own ta tolerance

## dataresqc's own ta tolerance table (from wmo_time_consistency), keyed by
## hour gap dt, capped at 12h for lookup (dt > 11 uses the 11-12h bin: 25 deg)
tol_ta <- function(dt) {
  dt <- pmin(dt, 12)
  tol <- rep(NA_real_, length(dt))
  tol[dt <= 1] <- 4
  tol[dt > 1  & dt <= 2]  <- 7
  tol[dt > 2  & dt <= 3]  <- 9
  tol[dt > 3  & dt <= 4]  <- 11
  tol[dt > 4  & dt <= 5]  <- 13
  tol[dt > 5  & dt <= 6]  <- 15
  tol[dt > 6  & dt <= 7]  <- 17
  tol[dt > 7  & dt <= 8]  <- 18
  tol[dt > 8  & dt <= 9]  <- 20
  tol[dt > 9  & dt <= 10] <- 22
  tol[dt > 10 & dt <= 11] <- 23
  tol[dt > 11]            <- 25
  tol
}

if (MODE == "report") {

  out <- data.frame(file=character(), id=character(),
                     Year=integer(), Month=integer(), Day=integer(),
                     midday_Hour=integer(), midday_Minute=integer(), midday_Value=numeric(),
                     evening_Hour=integer(), evening_Minute=integer(), evening_Value=numeric(),
                     dt=numeric(), diff=numeric(),
                     stringsAsFactors=FALSE)

  for (f in files) {
    x <- read_sef(f, all=TRUE)
    st <- read_meta(f)["id"]
    x <- x[!is.na(x$Value), ]
    x <- x[order(x$Year, x$Month, x$Day, x$Hour, x$Minute), ]
    x$daykey <- paste(x$Year, x$Month, x$Day)

    for (dk in unique(x$daykey)) {
      d <- x[x$daykey == dk, ]
      if (nrow(d) < 2) next
      midday  <- d[d$Hour >= MIDDAY_HOUR_MIN  & d$Hour < MIDDAY_HOUR_MAX, ]
      evening <- d[d$Hour >= EVENING_HOUR_MIN, ]
      if (nrow(midday) == 0 || nrow(evening) == 0) next

      for (i in seq_len(nrow(midday))) {
        for (j in seq_len(nrow(evening))) {
          t1 <- ISOdatetime(midday$Year[i],  midday$Month[i],  midday$Day[i],
                             midday$Hour[i],  midday$Minute[i],  0, tz="UTC")
          t2 <- ISOdatetime(evening$Year[j], evening$Month[j], evening$Day[j],
                             evening$Hour[j], evening$Minute[j], 0, tz="UTC")
          dt <- as.numeric(difftime(t2, t1, units="hours"))
          if (dt <= 0) next

          diff <- evening$Value[j] - midday$Value[i]
          thresh <- HALF * tol_ta(dt)

          if (diff > thresh) {
            out <- rbind(out, data.frame(
              file=basename(f), id=st,
              Year=midday$Year[i], Month=midday$Month[i], Day=midday$Day[i],
              midday_Hour=midday$Hour[i], midday_Minute=midday$Minute[i], midday_Value=midday$Value[i],
              evening_Hour=evening$Hour[j], evening_Minute=evening$Minute[j], evening_Value=evening$Value[j],
              dt=round(dt,2), diff=round(diff,2),
              stringsAsFactors=FALSE))
          }
        }
      }
    }
  }

  cat(nrow(out), "candidate diurnal-warming-reversal rows\n")
  write.table(out, file=candidate_file, row.names=FALSE, quote=FALSE, sep="\t")

} else if (MODE == "apply") {

  cand <- read.table(candidate_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

  for (f in files) {
    fname <- basename(f)
    sub <- cand[cand$file == fname, ]
    if (nrow(sub) == 0) next

    lines <- readLines(f, warn=FALSE)
    hdr_idx <- grep("^Year\tMonth", lines)[1]

    for (i in seq_len(nrow(sub))) {
      for (j in (hdr_idx+1):length(lines)) {
        parts <- strsplit(lines[j], "\t", fixed=TRUE)[[1]]
        if (as.integer(parts[1]) == sub$Year[i]  &&
            as.integer(parts[2]) == sub$Month[i] &&
            as.integer(parts[3]) == sub$Day[i]   &&
            as.integer(parts[4]) == sub$evening_Hour[i]  &&
            as.integer(parts[5]) == sub$evening_Minute[i]) {

          meta <- parts[8]
          if (grepl("qc=", meta, fixed=TRUE)) {
            meta <- sub("(qc=[^|]*)", "\\1;diurnal_warming_reversal", meta)
          } else {
            meta <- paste0(meta, " | qc=diurnal_warming_reversal")
          }
          parts[8] <- meta
          lines[j] <- paste(parts, collapse="\t")
          break
        }
      }
    }

    writeLines(lines, f, sep="\r\n")
  }

  cat("Applied", nrow(cand), "diurnal_warming_reversal flags across", length(unique(cand$file)), "file(s)\n")
}
