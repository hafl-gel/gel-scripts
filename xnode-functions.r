
# NOTES:
# json format
# 1 file per minute
# data structure of folders: year/mont/day
# file name contains date_time of recording

# TODO:
# add docu
# add function to read data from several days (read_xnode?)
#   -> check if zipped or not or both
# check if "date" is provided in to -> add 24h

## 0. header ----------------------------------------

require(ibts)
require(jsonlite)
require(data.table)

# check availability of jq
jq_available <- !inherits(try(system('type jq', intern = TRUE)), 'try-error')



## 1. functions ----------------------------------------

# select files in given time range
valid_times <- function(file_names, time_from = NULL, time_to = NULL) {
    # convert possible POSIXlt/POSIXct to integer (done in main function)
    if (length(time_from) || length(time_to)) {
        # careful - (zip) files are not sorted!!!
        times <- sort(
            as.integer(
                fast_strptime(sub('.*(\\d{14})[.]json$', '\\1', file_names),
                    format = '%Y%m%d%H%M%S', lt = FALSE, tz = 'UTC')
            ), index.return = TRUE
        )
    }
    # return index of valid files (times == end of 1-minute interval)
    switch(
        paste(length(time_from), length(time_to))
        , '1 1' = times$ix[times$x >= (time_from + 60) & times$x <= time_to]
        , '1 0' = times$ix[times$x >= (time_from + 60)]
        , '0 1' = times$ix[times$x <= time_to]
        , '0 0' = seq_along(file_names)
    )
}

# read zipped data directory
read_zip <- function(zip_dir, .from = NULL, .to = NULL, use_jq = FALSE) {
    # get file names
    list_all <- unzip(zip_dir, list = TRUE)
    files <- list_all[['Name']][list_all[['Length']] > 0]
    # select valid files
    ind_files <- valid_times(files, .from, .to)
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
    # get file names
    files <- dir(raw_dir, full.names = TRUE)
    # select valid files
    ind_files <- valid_times(files, .from, .to)
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

# read xnode data
read_xnode <- function(path, from = NULL, to = NULL, time_zone = 'Etc/GMT-1',
    simplify = TRUE, as_ibts = TRUE) {
    # check if path is pointing to top directory or data directories or single files
    dir_regex <- '^\\d{4}_\\d{2}_\\d{2}([.]zip)?$'
    nozip_regex <- sub('\\(.*', '$', dir_regex)
    if (length(path) == 1L && file.info(path)$isdir && !grepl(dir_regex, path)) {
        # check available dates
        all_dirs <- dir(path_data, pattern = dir_regex)
        if (length(all_dirs) == 0) {
            # TODO: cat no files in path
            return(invisible())
        }
    } else if (all(grepl(dir_regex, path))) {
        # directories & zips
        all_dirs <- path
    } else if (all(!file.info(path)$isdir)) {
        stop('providing individual files is not yet supported!')
    } else {
        stop('argument "path" is not valid')
    }
    # dates as integers (easier to select range)
    dates <- as.integer(gsub('(_|[.]zip)', '', all_dirs))
    # fix if there are both, zip & unzipped folder
    dups <- duplicated(dates)
    if (any(dups)) {
        remove_me <- NULL
        for (i in which(dups)) {
            # get entries
            dup_index <- which(dates == dates[i])
            remove_me <- c(remove_me, dup_index)
            dup_entries <- all_dirs[dup_index]
            # replace with zip or ignore any other
            if (any(is_zip <- grepl('[.]zip$', dup_entries))) {
                # add first zip entry
                all_dirs <- c(all_dirs, dup_entries[which(is_zip)[1]])
            } else {
                # add entry without extension
                all_dirs <- c(all_dirs, 
                    dup_entries[grep(nozip_regex, dup_entries)])
            }
        }
        # remove previous duplicates
        all_dirs <- all_dirs[-remove_me]
        # rebuild dates again
        dates <- as.integer(gsub('(_|[.]zip)', '', all_dirs))
    }
    ### get folders in timerange
    # from
    if (is.null(from)) {
        from_int <- dates[1]
    } else {
        from <- parse_date_time3(from, tz = time_zone)
        # dates as integers (easier to select range)
        from_int <- as.integer(format(from, '%Y%m%d'))
    }
    # to
    if (is.null(to)) {
        to_int <- dates[length(dates)]
    } else {
        to <- parse_date_time3(to, tz = time_zone)
        # dates as integers (easier to select range)
        to_int <- as.integer(format(to, '%Y%m%d'))
    }
    # convert to integer
    from <- as.integer(from)
    to <- as.integer(to)
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
        xnode_data <- xnode_data[, .(st, et, mgm3, ppm, temperature, pressure,
            humidity, supplyVoltage)]
        setattr(xnode_data, 'xnode', atts)
    }
    # convert to ibts?
    if (as_ibts) {
        xnode_data <- as.ibts(xnode_data)
        if (simplify) attr(xnode_data, 'xnode') <- atts
    }
    # return
    xnode_data
}


## 2. test with data ----------------------------------------

# # path
# base_path <- '~/LFE/02_Daten/5-draeger'
# gas <- 'nh3'
# sensor <- 'sensor-03'
# path_data <- file.path(base_path, gas, sensor)

# # select time range -> get relevant files
# # from <- NULL
# # to <- NULL
# from <- '12.11.2023'
# to <- '14.11.2023 20:00'

# xx <- read_xnode(path_data, from, to)
# yy <- read_xnode(path_data, from, to, as_ibts = TRUE)

