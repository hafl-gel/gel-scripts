This repository contains a collection of R scripts from the BFH/HAFL group GEL.

The easiest way to source individual scripts is with the aid of the package `devtools`.
E.g. to source the script to evaluate turbulence statistics from the GILL Windmaster sonic anemometer run:

```r
# script <- 'sonic-turbulence.r'
# raw_file <- paste0('https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/', script)
# devtools::source_url(raw_file)
devtools::source_url('https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/sonic-turbulence.r')
```

