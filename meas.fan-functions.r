
# NOTES:
# file paths:
# 'messventilator-630mm-1/frequi-1_63_2023_11_14.csv'

# TODO:

## 0. header ----------------------------------------

if (Sys.info()['sysname'] == 'Linux') {
    # Christoph
    path_calibration <- '~/repos/5_GitHub/gel-scripts/meas.fan'
} else {
    # Steffu
    path_calibration <- '~/git/gel-scripts/meas.fan'
}

## 1. functions ----------------------------------------

## 2. testing ----------------------------------------

path_lfe <- '~/LFE/02_Daten/9-Messventilator'
fan_id <- '1'
path_prep <- file.path(path_lfe, 'messventilator-%s0mm-%s')

# get diameter
fan_dia <- switch(as.character(fan_id)
    , '1' = 63
    , '2' = 
    , '3' = 
    , '4' = 
    , '5' = 82
    , 92
)

# create path
path_fan <- sprintf(path_prep, fan_dia, fan_id)

## 3. Calibration data ----------------------------------------

##  • fan calibration ====================

if (FALSE) {
    library(tabulizer)
    path_manual <- '~/LFE/02_Daten/9-Messventilator/AQC-HZ-G-B-AL00006.pdf'
    all_tabs <- extract_tables(file = path_manual)
    # str(all_tabs)
    tabs_mat <- cbind(all_tabs[[11]][, 8], all_tabs[[12]][, 4:5])
    fan_cfs <- apply(tabs_mat, 2, function(x) {
        tab <- apply(do.call(rbind, strsplit(
                sub(',', '.', x[-(1:4)])
                    , split = ' ')), 2, as.numeric)
        x11()
        plot(tab[, 2], tab[, 1])
        m <- lm(tab[, 1] ~ tab[, 2])
        abline(m)
        cfs <- coef(m)
        names(cfs)[2] <- 'Hz'
        attr(cfs, 'limits') <- range(tab[, 2])
        cfs
    }, simplify = FALSE)
    fan_cfs <- setNames(fan_cfs, sub('0', '', tabs_mat[1, ]))
    saveRDS(fan_cfs, file.path(path_calibration, 'fan-calibration.rds'))
} else {
    fan_cfs <- readRDS(file.path(path_calibration, 'fan-calibration.rds'))
}
    

##  • 630 mm ====================

