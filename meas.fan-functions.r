
# NOTES:
# file paths:
# 'messventilator-630mm-1/frequi-1_63_2023_11_14.csv'
# all file times in UTC
# from/to => time_zone = '' => CET or CEST == local time

# TODO:

## 0. header ----------------------------------------

library(data.table)
library(ibts)

switch(Sys.info()['user']
    , christoph = 
    , hac5 = {
        options(gel.scripts.path = '~/repos/5_GitHub/gel-scripts')
    }
)

## 1. functions ----------------------------------------

#
read_frequi <- function(path_data, rds_fan_calib,
    from = NULL, to = NULL, fan_id = NULL, time_zone = '', 
    as_ibts = TRUE) {
    # check path_data
    if (!file.exists(path_data)) stop('path_data: ', path_data, ' is not accessible')
    if (missing(rds_fan_calib)) {
        rds_fan_calib <- file.path(getOption('gel.scripts.path', ''), 'meas.fan/fan-calibration.rds')
    }
    if (!file.exists(rds_fan_calib)) stop('rds_fan_calib: ', rds_fan_calib, ' is not accessible')
    # get fan_id
    if (grepl('mm-\\d$', path_data)) {
        fan_id <- sub('.*(\\d)$', '\\1', path_data)
        # don't know how to do it in a better way
        path_data <- dirname(path_data)
    }
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
    path_prep <- file.path(path_data, 'messventilator-%s0mm-%s')
    path_fan <- sprintf(path_prep, fan_dia, fan_id)
    # get coefficients
    fan_cfs <- readRDS(rds_fan_calib)[[as.character(fan_dia)]]
    # get files
    files_fan <- dir(path_fan, full.names = TRUE)
    # get dates
    dates <- as.integer(sub('.*(\\d{4})_(\\d{2})_(\\d{2})[.]csv$', '\\1\\2\\3', files_fan))
    ### get folders in timerange
    # from
    if (is.null(from)) {
        from_int <- dates[1]
        from <- 0
    } else {
        from <- parse_date_time3(from, tz = time_zone)
        # dates as integers (easier to select range)
        from_int <- as.integer(format(from, '%Y%m%d', tz = 'UTC'))
    }
    # to
    if (is.null(to)) {
        to_int <- dates[length(dates)]
        to <- Inf
    } else {
        to <- parse_date_time3(to, tz = time_zone)
        # dates as integers (easier to select range)
        to_int <- as.integer(format(to, '%Y%m%d', tz = 'UTC'))
    }
    # index
    in_range <- from_int <= dates & to_int >= dates
    if (!any(in_range)) {
        cat('no data available in specified time range!\n')
        return(invisible())
    }
    files_fan <- files_fan[in_range]
    # read data
    fan_data <- rbindlist(lapply(files_fan, fread))
    setnames(fan_data, c('Time', 'Hz', 'V'))
    # fix time
    max_interval_secs <- 40
    fan_data[, c('st', 'et', 'Time') := {
        et <- as.POSIXct(round(Time))
        diff_et <- diff(et)
        diff_et[diff_et > max_interval_secs] <- 30
        st <- et - c(30, diff_et)
        list(st, et, NULL)
    }]
    # select from/to range
    fan_data <- fan_data[st >= from & et <= to]
    # add speed (m/s) & air flow (m3/s)
    fan_data[, c('speed_m_s', 'flow_m3_s', 'extrapol') := {
        speed <- Hz * fan_cfs[2] + fan_cfs[1]
        speed[Hz == 0] <- 0
        list(
            speed,
            speed * pi * (fan_dia / 200) ^ 2,
            Hz < attr(fan_cfs, 'limits')[1] | Hz > attr(fan_cfs, 'limits')[2]
        )
    }]
    # column order
    setcolorder(fan_data, c('st', 'et', 'flow_m3_s', 'extrapol', 'Hz', 'speed_m_s', 'V'))
    if (as_ibts) {
        fan_data <- as.ibts(fan_data)
    }
    # return
    fan_data
}


## 2. Calibration data ----------------------------------------

##  • fan calibration ====================

if (FALSE) {
    if (Sys.info()['sysname'] == 'Linux') {
        # Christoph
        path_calibration <- '~/repos/5_GitHub/gel-scripts/meas.fan'
    } else {
        # Steffu
        path_calibration <- '~/git/gel-scripts/meas.fan'
    }
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
}
    
## 3. testing ----------------------------------------

# path_lfe <- '~/LFE/02_Daten/9-Messventilator'
# id <- '1'
# file_calib <- '~/repos/5_GitHub/gel-scripts/meas.fan/fan-calibration.rds'

# xx <- read_frequi(path_lfe, '14.11.2023 12:00', '16.11.2023 20:00', file_calib, id)

