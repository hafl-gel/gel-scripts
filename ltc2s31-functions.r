
library(data.table)
library(ibts)


# read s31 data
read_s31 <- function(path_data, sensor, from, to, tz, as_ibts = FALSE) {
    # process from/to
    if (missing(from)) {
        from_psx <- from_int <- 1
    } else {
        # convert to POSIXct
        from_psx <- parse_date_time3(from, tz = tz)
        # convert to integer dates
        from_int <- as.integer(format(from_psx, '%Y%m%d'))
    }
    if (missing(to)) {
        to_psx <- to_int <- Inf
    } else {
        # convert to POSIXct
        to_psx <- parse_date_time3(to, tz = tz)
        # convert to integer dates
        to_int <- as.integer(format(to_psx, '%Y%m%d'))
    }
    # check path
    if (grepl('dragino-s31-\\d{2}/?$', path_data)) {
        path_files <- path_data
    } else {
        # get id
        id <- as.integer(sub('.*-(\\d{1,2})', '\\1', sensor))
        # check id
        if (is.na(id)) {
            stop('sensor name "', sensor, '" not recognized!')
        }
        # get path to files
        path_files <- sprintf('%s/dragino-s31-%02i', path_data, id)
    }
    # get files
    files <- dir(path_files)
    # get dates as integer
    dates <- gsub('.*-|_|\\.csv', '', files) |> 
        as.integer()
    # read files in time range
    dat <- files[dates >= from_int & dates <= to_int] |>
        file.path(path_files, . = _) |>
        lapply(FUN = fread) |>
        rbindlist()
    # fix times (s31 is always UTC) and round to seconds
    dat[, et := fast_strptime(time, format = '%Y-%m-%dT%H:%M:%OSZ+00:00', lt = FALSE) |>
        round() |> as.POSIXct()]
    # fix apparently sometimes occuring wrong time order
    setorder(dat, et)
    # recorded are 1-minute averages
    dat[, st := c(et[1] - 60, et[-.N])]
    # remove time column
    dat[, time := NULL]
    # subset by from/to
    dat <- dat[st >= from_psx & et <= to_psx]
    # st/et as first columns
    setcolorder(dat, c('st', 'et'))
    if (as_ibts) {
        dat <- as.ibts(dat)
    }
    dat
}

# read ltc2 data
read_ltc2 <- function(path_data, sensor, from, to, tz, as_ibts = FALSE) {
    # process from/to
    if (missing(from)) {
        from_psx <- from_int <- 1
    } else {
        # convert to POSIXct
        from_psx <- parse_date_time3(from, tz = tz)
        # convert to integer dates
        from_int <- as.integer(format(from_psx, '%Y%m%d'))
    }
    if (missing(to)) {
        to_psx <- to_int <- Inf
    } else {
        # convert to POSIXct
        to_psx <- parse_date_time3(to, tz = tz)
        # convert to integer dates
        to_int <- as.integer(format(to_psx, '%Y%m%d'))
    }
    # check path
    if (grepl('dragino-ltc2-\\d{2}/?$', path_data)) {
        path_files <- path_data
    } else {
        # get id
        id <- as.integer(sub('.*-(\\d{1,2})', '\\1', sensor))
        # check id
        if (is.na(id)) {
            stop('sensor name "', sensor, '" not recognized!')
        }
        # get path to files
        path_files <- sprintf('%s/dragino-ltc2-%02i', path_data, id)
    }
    # get files
    files <- dir(path_files)
    # get dates as integer
    dates <- gsub('.*-|_|\\.csv', '', files) |> 
        as.integer()
    # read files in time range
    dat <- files[dates >= from_int & dates <= to_int] |>
        file.path(path_files, . = _) |>
        lapply(FUN = fread) |>
        rbindlist()
    # fix times (ltc2 is always UTC) and round to seconds
    dat[, et := fast_strptime(time, format = '%Y-%m-%dT%H:%M:%OSZ+00:00', lt = FALSE) |>
        round() |> as.POSIXct()]
    # fix apparently sometimes occuring wrong time order
    setorder(dat, et)
    # recorded are 1-minute averages
    dat[, st := c(et[1] - 60, et[-.N])]
    # remove time column
    dat[, time := NULL]
    # subset by from/to
    dat <- dat[st >= from_psx & et <= to_psx]
    # st/et as first columns
    setcolorder(dat, c('st', 'et'))
    if (as_ibts) {
        dat <- as.ibts(dat)
    }
    dat
}
