
# functions

read_zip <- function(zip_dir, .from = NULL, .to = NULL, use_jq = FALSE) {
    list_all <- unzip(zip_dir, list = TRUE)
    files <- list_all[['Name']][list_all[['Length']] > 0]
    times <- fast_strptime(sub('.*(\\d{14})[.]json$', '\\1', files),
        format = '%Y%m%d%H%M%S', lt = FALSE, tz = time_zone)
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
                    paste0("unzip -p ", zip_dir, ' | jq -s "."'), 
                    intern = TRUE
                ), 
                flatten = TRUE
            )
        )
        out[ind_files, ]
    } else {
        # rbindlist(lapply(files[ind_files], function(f) {
        #     con <- unz(zip_dir, f)
        #     out <- parse_json(readLines(con))
        #     close(con)
        #     as.data.frame(out[lengths(out) > 0])
        # }), fill = TRUE)
        rbindlist(lapply(files[ind_files], function(f) {
            on.exit(close(con))
            con <- unz(zip_dir, f)
            out <- fromJSON(readLines(con), flatten = TRUE)
            as.data.frame(out)
        }), fill = TRUE)
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
dirs <- all_dirs[in_range]
n_dir <- length(dirs)
xnode_data <- vector('list', n_dir)

# check zip archives @ start & end
is_zip <- grepl('[.]zip$', dirs[in_range])

# read first directory
first_dir <- file.path(path_data, dirs[1])
if (is_zip[1]) {
    xnode_data[[1]] <- read_zip(first_dir, from, to, jq_available)
} else {
    browser()
    # TODO:
    first_files <- dir(first_dir, pattern = '^xnode', full.names = TRUE)
    # first_entries <- rbindlist(lapply(
}

# read in-between
if (n_dir > 2) {
    # loop over directories
    for (i in seq.int(2, n_dir - 1, 1)) {
        if (is_zip[i]) {
            xnode_data[[i]] <- read_zip(file.path(path_data, dirs[i]), 
                use_jq = jq_available)
        } else {
            # TODO:
            browser()
        }
    }
}

# check last directory
if (n_dir > 1) {
    if (is_zip[n_dir]) {
    } else {
    }
}

# check files in timerange
files_list <- lapply(unzip(

# reshape data

# rename
setnames(first_entries, sub('[.]', '_', sub('^[^.]*[.]', '', names(first_entries))))

# convert to ibts?

