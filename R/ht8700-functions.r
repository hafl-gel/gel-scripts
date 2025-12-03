
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

