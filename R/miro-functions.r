
#' Read MIRO Raw Data
#'
#' This function reads raw data from MIRO files. It supports multiple file formats
#' and handles different structures of the data.
#'
#' @param file_path A string specifying the path to the MIRO file.
#' @return A data.table containing the processed MIRO data or NULL if the file is empty.
#' @export
#'
# main function to read MIRO raw data
read_miro <- function(file_path, from = NULL, to = NULL, tz = 'UTC') {
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
        # if no dates => read
        read_me <- read_me | is.na(read_me)
        # loop over files
        if (any(read_me)) {
            out <- lapply(file.path(file_path, files[read_me]), read_miro,
                from = from, to = to, tz = tz)
            # return sorted
            rbindlist(out)[order(Time)]
        } else {
            stop('No files available within provided time range!')
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
        } else if (grepl('^\\d{4}-\\d{2}-\\d{2} MGA SN\\d{2,}\\.txt$', bn)) {
            stop('todo! (read local miro file => fix Times in file & date in filename!)...')
        } else if (!grepl('^(py_)?fnf_0\\d_miro', bn)) {
            # wrong file name
            cat('skip (filename not valid)\n')
            warning('data filename not valid -> skipping file "', bn, '"')
            return(NULL)
        }
        ### read File
        if (grepl('[.]gz$', bn)) {
            # gzip-ped data
            raw <- miro_read_loggerbox_cpp_gzip(normalizePath(file_path))
        } else {
            # uncompressed data
            raw <- miro_read_loggerbox_cpp(normalizePath(file_path))
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
        # get miro vs. loggerbox time offset
        out[, time_miro := fast_strptime(time_miro, '%d.%m.%Y %H:%M:%OS', lt = FALSE)]
        out[, miro_pc_offset_secs := as.numeric(time_miro - Time, units = 'secs')]
        # remove V1
        out[, c('time_string', 'time_miro') := NULL]
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


# add sonic data to miro data
add_sonic <- function(x, path) {
    # get from/to from ht-data
    tr <- x[, range(Time)]
    ft <- as.integer(format(tr, format = '%Y%m%d'))
    # get sonic file paths
    files <- dir(path, pattern = 'sonic')
    int_files <- as.integer(sub('.*_(\\d{8})_\\d{6}([.]gz)?$', '\\1', files))
    ind_files <- int_files >= ft[1] & int_files <= ft[2]
    # read sonic data
    dat <- rbindlist(lapply(file.path(path, files[ind_files]), read_windmaster_ascii))
    # merge data
    merge_data(x, dat)
}

