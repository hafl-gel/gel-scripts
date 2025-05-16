This repository contains a collection of R scripts from the BFH/HAFL group GEL,
which is slowly transformed into a package called 'gel' ([-> get latest release](https://github.com/hafl-gel/gel-scripts/releases/latest)).

At the moment, the package 'gel' contains:
- functions to process raw data from our instruments
    - sonic anemometers
    - HT-8700
    - LI-7500
- functions to process turbulence and flux measurements

All functions which are not yet included in the package are still available in the 'remaining-scripts' directory.

The easiest way to source individual scripts is with the aid of the package `devtools`.
E.g. to source the script to handle raw data downloaded from the MeteoSwiss IDAweb, run:

```r
devtools::source_url('https://rawgithubusercontent.com/hafl-gel/gel-scripts/main/remaining-scripts/idaweb.r')
```
