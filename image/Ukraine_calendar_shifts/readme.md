bestlag.txt files have date, Best\_lag, Best\_corr as ccolumns.
If a station has more than one best_lag file, it's because it was tested for different window sizes and step sizes.
Subtracting two consecutive days in the Date column tells you the step size, e.g.

1833-05-31 -13 0.741
1833-06-30 -13 0.312
1833-07-30 -14 0.413

means the step size is 30 days.


**Dnipro**: transition occurs on 1838-12-15
