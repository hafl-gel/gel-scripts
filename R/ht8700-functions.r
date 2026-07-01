
#' Read HT8700 Raw Data
#'
#' This function reads raw data from HT8700 files. It supports multiple file formats
#' and handles different structures of the data.
#'
#' @param file_path A string specifying the path to the HT8700 file.
#' @return A data.table containing the processed HT8700 data or NULL if the file is empty.
#' @export
#'
# main function to read HT8700 raw data
read_ht8700 <- function(file_path, from = NULL, to = NULL, tz = 'UTC') {
    # parse from/to
    from <- parse_date_time3(from, tz = tz)
    to <- parse_date_time3(to, tz = tz)
    # check if directory is provided
    if (file.info(file_path)$isdir) {
        files <- dir(file_path)
        if (length(files) == 0L) {
            stop('directory is empty!')
        }
        # check dates in filenames
        suppressWarnings(
            file_dates <- as.integer(
                sub('.*(20[2-9]\\d)_?(\\d{2})_?(\\d{2}).*', '\\1\\2\\3', files)
            )
        )
        # check from & to
        read_me <- rep(TRUE, length(file_dates))
        if (length(from) > 0) {
            read_me <- read_me & file_dates >= as.integer(
                        format(floor_date(from, unit = 'days'), '%Y%m%d')
                        )
        }
        if (length(to) > 0) {
            read_me <- read_me & file_dates <= as.integer(
                        format(floor_date(to, unit = 'days'), '%Y%m%d')
                        )
        }
        # if no dates => drop
        read_me <- which(read_me)
        # loop over files
        if (any(read_me)) {
            out <- lapply(file.path(file_path, files[read_me]), read_ht8700,
                from = from, to = to, tz = tz)
            # return sorted
            rbindlist(out)[order(Time)]
        } else {
            stop('No valid files available within provided time range!')
        }
    } else {
        # be verbose
        cat("File:", path.expand(file_path), "- ")
        # get file name
        bn <- basename(file_path)
        # check file name
        if (grepl('[.](qdata|qd)$', bn)) {
            return(alloc.col(qs2::qd_read(file_path)))
        } else if (grepl('[.]qs$', bn)) {
            if (!requireNamespace('qs')) {
                stop('data is provided as *.qs file -> install qs library',
                    ' running "install.packages("qs")"')
            }
            return(alloc.col(qs::qread(file_path)))
        } else if (grepl('[.]rds$', bn)) {
            return(readRDS(file_path))
        } else if (is_old_structure <- grepl('^ht8700_', bn)) {
            # read with old function
            return(read_ht8700_old(file_path, tz = 'Etc/GMT-1'))
        } else if (!grepl('^(py_)?fnf_0\\d_ht8700', bn)) {
            # wrong file name
            cat('skip (filename not valid)\n')
            warning('data filename not valid -> skipping file "', bn, '"')
            return(NULL)
        }
        ### read File
        if (grepl('[.]gz$', bn)) {
            # gzip-ped data
            raw <- ht8700_read_cpp_gzip(normalizePath(file_path))
        } else {
            # uncompressed data
            raw <- ht8700_read_cpp(normalizePath(file_path))
        }
        # convert to data.table
        out <- as.data.table(raw)
        # check empty
        if (nrow(out) == 0) {
            cat('no valid data!\n')
            return(NULL)
        }
        # fix time column
        out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE, 
            tz = 'UTC')]
        # remove NA lines that come from conversion
        out <- na.omit(out, cols = 'Time')
        # check if empty again
        if (out[, .N == 0]) {
            cat('file empty\n')
            return(NULL)
        }
        # remove V1
        out[, time_string := NULL]
        # filter by from/to
        if (length(from) > 0) {
            out <- out[Time >= from]
        }
        if (length(to) > 0) {
            out <- out[Time <= to]
        }
        # place Time column first
        setcolorder(out, 'Time')
        cat('done\n')
        # return
        out
    }
}

# helper function to decode alarm codes
decode_alarm <- function(lower, upper) {
    lo <- strtoi(lower, 16L)
    up <- strtoi(upper, 16L)
    if (length(lo) != length(up)) {
        stop('arguments "lower" and "upper" must be of equal lengths!')
    }
    nms <- paste(lower, upper, sep = '-')
    out <- decal(lo, up)
    names(out) <- nms
    rownames(attr(out, 'mat')) <- nms
    out
}

# x: raw ht8700 data
# add: FALSE returns alarm codes only, TRUE returns all ht data
# simple: TRUE only provide 1 column with codes, FALSE 1 column per code containing 
#   TRUE/FALSE values
get_alarms <- function(x, add = FALSE, simple = !add) {
    out <- copy(x)
    add_alarms(out, simple = simple)
    if (!add) {
        if (simple) {
            out <- out[, alarm_codes]
        } else {
            out <- out[, alarm_codes:ac_32]
        }
    }
    out
}

# same as `get_alarms()` but modifies data.table (ht data) in-place
add_alarms <- function(x, simple = FALSE) {
    x[, alarm_codes := '']
    if (simple) {
        x[, alarm_codes := {
            alarms <- decode_alarm(.BY[['alarm_lower_bit']], .BY[['alarm_upper_bit']])
            sapply(alarms, paste, collapse = ',')
        }, by = .(alarm_lower_bit, alarm_upper_bit)]
    } else {
        x[, paste0('ac_', 1:32) := FALSE]
        x[, c('alarm_codes', paste0('ac_', 1:32)) := {
            alarms <- decode_alarm(.BY[['alarm_lower_bit']], .BY[['alarm_upper_bit']])
            c(
                sapply(alarms, paste, collapse = ','),
                as.list(attr(alarms, 'mat'))
            )
        }, by = .(alarm_lower_bit, alarm_upper_bit)]
    }
    invisible(x)
}

# old main function to read HT8700 raw data
read_ht8700_old <- function(FilePath, tz = "Etc/GMT-1"){
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	# get file name
	bn <- basename(FilePath)
    # check file name
    if (is_old_structure <- grepl('^ht8700_', bn)) {
        # old data structure on sonic boxes
        data_pattern <- '^\\d{2}[0-9.:]+,([ A-Z0-9.+-]+,){18}[[ 0-9.+-]+$'
        # possibly gzipped files -> need R.utils package
        if (!any(grepl('R.utils', installed.packages()[, 'Package']))) {
            stop('package "R.utils" must be installed to process gz files!')
        }
        # set data time zone
        tz_data <- 'Etc/GMT-1'
    } else if (grepl('^(py_)?fnf_0\\d_ht8700', bn)) {
        # new data structure on loggerbox
        data_pattern <- '^\\d{4}-\\d{2}-\\d{2}T\\d{2}[0-9.:]+Z,([ A-Z0-9.+-]+,){18}[[ 0-9.+-]+$'
        # set data time zone
        tz_data <- 'UTC'
    } else {
        # wrong file name
        stop('data filename not valid')
    }
	### read File
    raw <- readLines(FilePath, warn = FALSE)
    # filter out erroneous multibyte strings
    raw <- raw[grepl(data_pattern, raw, useBytes = TRUE)]
    # fix one single line
    if (length(raw) == 1) {
        raw <- c(raw, '')
    }
    # read from string
    out <- fread(text = raw, blank.lines.skip = TRUE,
        header = FALSE, na.strings = '999.99', showProgress = FALSE,
    )
    # check empty
    if (nrow(out) == 0) {
        cat('no valid data!\n')
        return(NULL)
    }
    # fix time column
    if (is_old_structure) {
        # get date
        Date <- gsub("^ht8700_.*_([0-9]{8})_.*", "\\1", bn)
        # fix time
        out[, Time := fast_strptime(paste(Date, V1), lt = FALSE, 
            format = "%Y%m%d %H:%M:%OS", tz = tz_data)]
    } else if (out[, is.character(V1)]) {
        # V1 was not converted to POSIXct by fread() call
        out[, Time := fast_strptime(V1, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE, tz = tz_data)]
    } else {
        # V1 should be POSIXct from fread() call
        out[, Time := V1]
    }
    # convert column types
    char_cols <- paste0('V', c(2, 17, 18, 20))
    num_cols <- paste0('V', 3:16)
    suppressWarnings(out[, (char_cols) := lapply(.SD, as.character), .SDcols = char_cols])
    suppressWarnings(out[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols])
    suppressWarnings(out[, V19 := as.integer(V19)])
    # lower bits
    suppressWarnings(out[, V17 := as.character(V17)])
    # upper bits
    suppressWarnings(out[, V18 := as.character(V18)])
    # remove NA lines that come from conversion
    out <- na.omit(out)
    # check if empty again
    if (out[, .N == 0]) {
        cat('file empty\n')
        return(NULL)
    }
    # remove V1
    out[, V1 := NULL]
    # fix column names
    setnames(out,
        c(
            'sn', # column 2
            'nh3_ppb', # column 3
            'nh3_ugm3', # column 4
            'rh_int', # column 5
            'temp_int', # column 6
            'temp_amb', # column 7
            'press_amb', # column 8
            'oss', # column 9
            'peak_pos', # column 10
            'temp_leaser_chip', # column 11
            'temp_leaser_housing', # column 12
            'temp_mct', # column 13
            'temp_mct_housing', # column 14
            'laser_current', # column 15
            'ref_road_2f', # column 16
            'alarm_lower_bit', # column 17
            'alarm_upper_bit', # column 18
            'cleaning_flag', # column 19
            'notused', # column 20
            'Time' # column 21
        )
    )
    # place Time column first
    setcolorder(out, 'Time')
    cat('done\n')
    # return
    out
}

