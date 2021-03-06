### ******************************************************************************
### ******************************************************************************
### ******************* functions for miniDOAS data evaluation *******************
### ******************************************************************************
### ******************************************************************************

programVersion <- "v5.2 (26.12.2020)"

### libraries etc
### ******************************************************************************
library(lubridate)
library(robustbase)
require(MASS)
require(data.table)
require(ibts)
# options(shiny.trace=FALSE)


if (dir.exists('~/repos/5_GitHub/gel-scripts')) {
    # local to gel
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-filters.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-fitting.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-plotting.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-spectrometer.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-others.r")
    source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-tools.r")
} else {
    # remote
    devtools::source_url("https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/minidoas/miniDOAS-functions-filters.r")
    devtools::source_url("https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/minidoas/miniDOAS-functions-fitting.r")
    devtools::source_url("https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/minidoas/miniDOAS-functions-plotting.r")
    devtools::source_url("https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/minidoas/miniDOAS-functions-spectrometer.r")
    devtools::source_url("https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/minidoas/miniDOAS-functions-others.r")
    devtools::source_url("https://rawgithubusercontent.com/ChHaeni/gel-scripts/main/minidoas/miniDOAS-functions-tools.r")
}

cat("\n*** Skript", programVersion, "***\n")
