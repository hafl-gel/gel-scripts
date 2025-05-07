
# NOTES:
# file paths:
# 'messventilator-630mm-1/frequi-1_63_2023_11_14.csv'
# all file times in UTC
# from/to => time_zone = '' => CET or CEST == local time

# TODO:

## 0. header ----------------------------------------

library(data.table)
library(ibts)
library(sodium)

## 1. functions ----------------------------------------

#
read_frequi <- function(path_data, from = NULL, to = NULL, time_zone = 'UTC', 
    V_crit = 3, max_interval_secs = 90, as_ibts = TRUE, fan_dia = NULL) {
    # check path_data
    if (!file.exists(path_data)) stop('path_data: ', path_data, ' is not accessible')
    # get files
    files_fan <- dir(path_data, full.names = TRUE)
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
    # get files within time range
    files_fan <- files_fan[in_range]
    # read data
    fan_data <- rbindlist(lapply(files_fan, fread, 
                fill = TRUE, blank.lines.skip = TRUE))
    # new path structure, check unique diameter
    if (is.null(fan_dia)) {
        suppressWarnings(
            dias <- as.integer(sub('.*_(\\d{3})mm-.*', '\\1', files_fan)) / 10
        )
        if (all(is.na(dias))) {
            # old file structure
            dias <- as.integer(sub('.*-(\\d{3})mm-.*', '\\1', files_fan)) / 10
        }
        if (length(fan_dia <- unique(dias)) != 1) {
            stop('fan diameter is not unique for provided time range (', paste(fan_dia, collapse = ' and '), ')')
        }
    } 
    if (fan_data[, is.character(V1)]) {
        fan_data[, V1 := fast_strptime(V1, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    }
    # sort by date/time since sorting is not guaranteed from filenames
    fan_data <- fan_data[order(V1)]
    
    # check columns 2 & 3
    if (fan_data[, !is.numeric(V2)]) {
        # filter bad entries
        fan_data <- fan_data[grep('^\\d+([.]\\d+)?$', V2)]
        # convert to numeric
        fan_data[, V2 := as.numeric(V2)]
    }
    if (fan_data[, !is.numeric(V3)]) {
        # filter bad entries
        fan_data <- fan_data[grep('^\\d+([.]\\d+)?$', V3)]
        # convert to numeric
        fan_data[, V3 := as.numeric(V3)]
    }
    # omit any NA values
    fan_data <- na.omit(fan_data)
    setnames(fan_data, c('Time', 'Hz', 'V'))
    # fix time
    fan_data[, c('st', 'et', 'Time') := {
        et <- as.POSIXct(round(Time))
        diff_et <- as.numeric(diff(et), units = 'secs')
        diff_et[diff_et > max_interval_secs] <- 60
        st <- et - c(60, diff_et)
        list(st, et, NULL)
    }]
    # select from/to range
    fan_data <- fan_data[st >= from & et <= to]
    # get coefficients
    fan_coefficients <- attr(read_frequi, 'fan_coefficients')[[as.character(fan_dia)]]
    # add speed (m/s) & air flow (m3/s)
    fan_data[, c('speed_m_s', 'flow_m3_s', 'extrapol') := {
        speed <- Hz * fan_coefficients[2] + fan_coefficients[1]
        speed[Hz == 0] <- 0
        list(
            speed,
            speed * pi * (fan_dia / 200) ^ 2,
            Hz > 0 & (Hz < attr(fan_coefficients, 'limits')[1] | Hz > attr(fan_coefficients, 'limits')[2])
        )
    }]
    # set flow values below critical V to NA
    fan_data[V <= V_crit, flow_m3_s := NA_real_]
    # column order
    setcolorder(fan_data, c('st', 'et', 'flow_m3_s', 'extrapol', 'Hz', 'speed_m_s', 'V'))
    if (as_ibts) {
        fan_data <- as.ibts(fan_data)
    }
    # return & attach fan dia
    structure(fan_data, 'fan_dia' = fan_dia)
}


## 2. Calibration data ----------------------------------------

dat2string <- function(x) {
    key <- keygen(as.raw(rep(1, 32)))
    nonce <- as.raw(rep(1, 24))
    xs <- qs2::qd_serialize(x)
    xr <- data_encrypt(xs, key, nonce = nonce)
    paste0(xr, collapse = '')
}
string2dat <- function(x) {
    key <- keygen(as.raw(rep(1, 32)))
    nonce <- as.raw(rep(1, 24))
    xs <- strsplit(x, split = '(?<=.)(?!(..)*.$)', perl = TRUE)
    xr <- as.raw(paste0('0x', unlist(xs)))
    qs2::qd_deserialize(data_decrypt(xr, key, nonce = nonce))
}
splitN <- function(x, n) {
    unlist(strsplit(x, 
            split = paste0('(?<=.{', n, '})(?!(.{', n, '})*.$)'), perl = TRUE))
}
add_to_eof <- function(x, attr_name, header, f_name, ec_file) {
    con <- file(ec_file, open = 'a')
    xs <- dat2string(x)
    on.exit({close(con)})
    writeLines(paste0('# ', header), con)
    writeLines(paste0(
            'attr(', f_name, ', "', attr_name, '") <- string2dat(paste0(c(\n"',
            paste(splitN(xs, 70), collapse = '",\n"'), '"), collapse = ""))\n'), con)
}

##  • fan calibration ====================

if (FALSE) {
    # remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
    # library(tabulizer)
    library(tabulapdf)
    path_manual <- '~/LFE/02_Daten/9-Messventilator/AQC-HZ-G-B-AL00006.pdf'
    all_tabs <- extract_tables(file = path_manual)
    # str(all_tabs)
    tabs_mat <- cbind(all_tabs[[11]][, 8], all_tabs[[12]][, 4:5])
    fan_coefficients <- apply(tabs_mat, 2, function(x) {
        tab <- apply(do.call(rbind, strsplit(
                sub(',', '.', x[-(1:4)])
                    , split = ' ')), 2, as.numeric)
        x11()
        plot(tab[, 2], tab[, 1])
        m <- lm(tab[, 1] ~ tab[, 2])
        abline(m)
        cfs <- coef(m)
        names(cfs)[2] <- 'Hz'
        attr(cfs, 'limits') <- range(tab[, 2], na.rm = TRUE)
        cfs
    }, simplify = FALSE)
    fan_coefficients <- setNames(fan_coefficients, sub('0', '', colnames(tabs_mat)))
    # path_calibration <- '~/repos/5_GitHub/gel-scripts/meas.fan'
    # saveRDS(fan_coefficients, file.path(path_calibration, 'fan-calibration.rds'))
    # readRDS(file.path(path_calibration, 'fan-calibration.rds'))
    # add coefficients to function
    # add values to script
    add_to_eof(fan_coefficients, 'fan_coefficients', '~~~ fan coefficients ~~~', 
        'read_frequi', '~/repos/5_GitHub/gel-scripts/meas.fan-functions.r')
}
    
## 3. Testing ----------------------------------------

# path_lfe <- '~/LFE/02_Daten/9-Messventilator'
# id <- '1'
# file_calib <- '~/repos/5_GitHub/gel-scripts/meas.fan/fan-calibration.rds'

# xx <- read_frequi(path_lfe, '14.11.2023 12:00', '16.11.2023 20:00', file_calib, id)

## 4. Calibration Results ----------------------------------------

# ~~~ fan coefficients ~~~
attr(read_frequi, "fan_coefficients") <- string2dat(paste0(c(
"1c15782abf437a97048e915491c9e4cff9f65c9cabef3907dd7f3e396a07fd6f882ef0",
"2b49c0c9ec66aaffb1553c091642ea808ddfd186351b7c57b47dc770f5552241652de2",
"008ec3bfab4d77ec3ba6a31493f8db91fb807b6bbde7e1e337316d90f0643beb1dc5a6",
"3055a7490353625ab49c45a9eecf9536a7ff4c21c99920c456e83d39a66d867f53f945",
"01fcc303d5d21d5a138e45057e7654e8c800592b5528a52e389dc76497e6c997f83442",
"d90aad76776493e6"), collapse = ""))

