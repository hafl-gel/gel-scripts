
# json format
# 1 file per minute
# data structure of folders: year/mont/day
# file name contains date_time of recording

library(ibts)

tzone <- 'UTC'

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
files <- dir(path_data, pattern = '^\\d{4}_\\d{2}_\\d{2}([.]zip)?$')
if (length(files) == 0) {
    # TODO: cat no files in path
    return(invisible())
}
dates <- as.integer(gsub('(_|[.]zip)', '', files))

### get folders in timerange
# from
if (is.null(from)) {
    from_int <- dates[1]
} else {
    from <- parse_date_time3(from, tz = tzone)
    from_int <- as.integer(format(from, '%Y%m%d'))
}
# to
if (is.null(to)) {
    to_int <- dates[length(dates)]
} else {
    to <- parse_date_time3(to, tz = tzone)
    to_int <- as.integer(format(to, '%Y%m%d'))
}
# index
in_range <- from_int <= dates & to_int >= dates
if (!any(in_range)) {
    # TODO: cat no data in specified time range
    return(invisible())
}

# read files

# reshape data

# convert to ibts?

