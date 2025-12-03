
# R wrapper, main function
read_licor <- function(file_path, from = NULL, to = NULL, tz = 'UTC') {
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
            out <- lapply(file.path(file_path, files[read_me]), read_licor,
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
        } else if (!grepl('^(py_)?fnf_0\\d_licor_', bn)) {
            # wrong file name
            cat('skip (filename not valid)\n')
            warning('data filename not valid -> skipping file "', bn, '"')
            return(NULL)
        }
        if (grepl('[.]gz$', bn)) {
            raw_list <- licor_read_cpp_gzip(normalizePath(file_path))
        } else {
            raw_list <- licor_read_cpp(normalizePath(file_path))
        }
        out <- data.table::as.data.table(raw_list)
        # convert time
        out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
        # remove time_string
        out[, time_string := NULL]
        # filter by from/to
        if (length(from) > 0) {
            out <- out[Time >= from]
        }
        if (length(to) > 0) {
            out <- out[Time <= to]
        }
        # place Time first
        setcolorder(out, 'Time')
        cat('done\n')
        # remove invalid Time entries
        na.omit(out, cols = 'Time')
    }
}
