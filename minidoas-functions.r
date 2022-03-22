### ******************************************************************************
### ******************************************************************************
### ******************* functions for miniDOAS data evaluation *******************
### ******************************************************************************
### ******************************************************************************

programVersion <- "v5.2 (26.12.2020)"

### libraries etc
### ******************************************************************************
# library(data.table)
library(lubridate)
library(robustbase)
require(MASS)
# options(shiny.trace=FALSE)


if (dir.exists('~/repos/5_GitHub/gel-scripts')) {
    # local
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-filters.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-fitting.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-plotting.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-spectrometer.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-others.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-tools.r")
} else {
    # remote
    devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-filters.r?raw=TRUE")
    devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-fitting.r?raw=TRUE")
    devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-plotting.r?raw=TRUE")
    devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-spectrometer.r?raw=TRUE")
    devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-others.r?raw=TRUE")
    devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-tools.r?raw=TRUE")
}

cat("\n*** Skript", programVersion, "***\n")
