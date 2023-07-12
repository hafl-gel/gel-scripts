
## 0. functions ----------------------------------------

# read zipped data directory
read_zip <- function(zip_dir, .from = NULL, .to = NULL, use_jq = FALSE) {
    list_all <- unzip(zip_dir, list = TRUE)
    files <- list_all[['Name']][list_all[['Length']] > 0]
    if (!is.null(.from) || !is.null(.to)) {
        times <- fast_strptime(sub('.*(\\d{14})[.]json$', '\\1', files),
            format = '%Y%m%d%H%M%S', lt = FALSE, tz = time_zone)
    }
    ind_files <- switch(paste(is.null(.from), is.null(.to))
        , 'FALSE FALSE' = times >= .from & times <= .to
        , 'FALSE TRUE' = times >= .from
        , 'TRUE FALSE' = times >= .to
        , 'TRUE TRUE' = seq_along(files)
    )
    if (use_jq) {
        out <- as.data.table(
            fromJSON(
                system(
                    paste0('unzip -p ', zip_dir, ' | jq -s "."'), 
                    intern = TRUE
                ), 
                flatten = TRUE
            )
        )
        out[ind_files, ]
    } else {
        rbindlist(lapply(files[ind_files], function(f) {
            on.exit(close(con))
            con <- unz(zip_dir, f)
            out <- fromJSON(readLines(con), flatten = TRUE)
            as.data.frame(out)
        }), fill = TRUE, use.names = TRUE)
    }
}

# read unzipped data directory
read_raw <- function(raw_dir, .from = NULL, .to = NULL, use_jq = FALSE) {
    files <- dir(raw_dir, full.names = TRUE)
    if (!is.null(.from) || !is.null(.to)) {
        times <- fast_strptime(sub('.*(\\d{14})[.]json$', '\\1', files),
            format = '%Y%m%d%H%M%S', lt = FALSE, tz = time_zone)
    }
    ind_files <- switch(paste(is.null(.from), is.null(.to))
        , 'FALSE FALSE' = times >= .from & times <= .to
        , 'FALSE TRUE' = times >= .from
        , 'TRUE FALSE' = times >= .to
        , 'TRUE TRUE' = seq_along(files)
    )
    if (use_jq) {
        out <- as.data.table(
            fromJSON(
                system(
                    paste0('jq -s "." ', paste(files, collapse = ' ')),
                    intern = TRUE
                ), 
                flatten = TRUE
            )
        )
        out[ind_files, ]
    } else {
        rbindlist(lapply(files[ind_files], function(f) {
            out <- read_json(f)
            as.data.frame(out[lengths(out) > 0])
        }), fill = TRUE, use.names = TRUE)
    }
}

# check availability of jq
jq_available <- !inherits(system('jq --help', intern = TRUE), 'try-error')

# json format
# 1 file per minute
# data structure of folders: year/mont/day
# file name contains date_time of recording

library(ibts)
library(jsonlite)
library(data.table)


time_zone <- 'UTC'

# path
base_path <- '~/repos/5_GitHub/gel-scripts/.test-data/xnode'
gas <- 'nh3'
sensor <- 'sensor1'

path_data <- file.path(base_path, gas, sensor)

# select time range -> get relevant files
# from <- NULL
# to <- NULL
from <- '06.07.2023 22:00'
to <- '11.07.2023 07:00'

# check available dates
all_dirs <- dir(path_data, pattern = '^\\d{4}_\\d{2}_\\d{2}([.]zip)?$')
if (length(all_dirs) == 0) {
    # TODO: cat no files in path
    return(invisible())
}
dates <- as.integer(gsub('(_|[.]zip)', '', all_dirs))

### get folders in timerange
# from
if (is.null(from)) {
    from_int <- dates[1]
} else {
    from <- parse_date_time3(from, tz = time_zone)
    from_int <- as.integer(format(from, '%Y%m%d'))
}
# to
if (is.null(to)) {
    to_int <- dates[length(dates)]
} else {
    to <- parse_date_time3(to, tz = time_zone)
    to_int <- as.integer(format(to, '%Y%m%d'))
}
# index
in_range <- from_int <= dates & to_int >= dates
if (!any(in_range)) {
    # TODO: cat no data in specified time range
    return(invisible())
}
dirs <- file.path(path_data, all_dirs[in_range])
n_dir <- length(dirs)
xnode_list <- vector('list', n_dir)

# check zip archives @ start & end
read_fun <- ifelse(grepl('[.]zip$', dirs), list(read_zip), list(read_raw))

# read first directory
xnode_list[[1]] <- read_fun[[1]](dirs[1], from, to, jq_available)

# read in-between (if more than 2 directories)
if (n_dir > 2) {
    # loop over directories
    for (i in seq.int(2, n_dir - 1, 1)) {
        xnode_list[[i]] <- read_fun[[i]](dirs[i],
            use_jq = jq_available)
    }
}

# read last directory (if more than 1 directory)
if (n_dir > 1) {
    xnode_list[[n_dir]] <- read_fun[[n_dir]](dirs[n_dir], .to = to, use_jq = jq_available)
}

# bind list
xnode_data <- rbindlist(xnode_list, fill = TRUE, use.names = TRUE)

# rename
setnames(xnode_data, sub('^.*[.]', '', names(xnode_data)))
setnames(xnode_data, c('value', 'name'), c('ppm', 'gas'))

# parse time
xnode_data[, et := fast_strptime(time, format = '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE,
    tz = time_zone)]
# sort by et
setorder(xnode_data, et)
# get st
xnode_data[, st := {
    dtime <- pmin(diff(et), 60)
    et - c(60, dtime)
}]

# add mg/m3?
xnode_data[, mgm3 := {
    # n_G / n_L * 10^6 = ppm
    # n_G = m_G / M_G
    # n_L = p_L * V_L / R / T_L (<- pV=nRT)
    # m_G / V_L = n_G / n_L * p_L / (R * T_L) * M_G -> M_G in mg/mol...
    pressure * ppm * 1e-1 / (8.314 * (273.15 + temperature)) * 
        c('nh3' = 17, 'co2' = 44)[tolower(gas)]
}]

# simplify?
if (simplify) {
    atts <- xnode_data[!is.na(gas), I(list(
        device_id = unique(device_id),
        dev_eui = unique(dev_eui),
        join_eui = unique(join_eui),
        dev_addr = unique(dev_addr),
        application_id = unique(application_id),
        gas = unique(gas),
        unit = unique(unit)
        ))]
    xnode_data <- xnode_data[, .(st, et, humidity, pressure, supplyVoltage, 
        temperature, ppm, mgm3)]
    setattr(xnode_data, 'xnode', atts)
}

# convert to ibts?

