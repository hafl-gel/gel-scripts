### ******************************************************************************
### ******************************************************************************
### ******************* functions for miniDOAS data evaluation *******************
### *************************** author: Joerg Sintermann *************************
### this works with the Avaspec and QE65-pro from Jan.2014 on, given converted raw-files by python script
### ******************************************************************************
### ******************************************************************************

programVersion <- "v5.2 (26.12.2020)"

### libraries etc
### ******************************************************************************
# library(data.table)
library(lubridate)
library(IDPmisc)
library(robustbase)
# options(shiny.trace=FALSE)


source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS_functions_filters.R")
source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS_functions_fitting.R")
source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS_functions_plotting.R")
source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS_functions_spectrometer.R")
source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS_functions_others.R")
source("~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS_functions_tools.R")
# devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-filters.r?raw=TRUE")
# devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-fitting.r?raw=TRUE")
# devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-plotting.r?raw=TRUE")
# devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-spectrometer.r?raw=TRUE")
# devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-others.r?raw=TRUE")
# devtools::source_url("https://github.com/ChHaeni/gel-scripts/blob/main/minidoas/miniDOAS-functions-tools.r?raw=TRUE")

cat("\n*** Skript", programVersion, "***\n")
