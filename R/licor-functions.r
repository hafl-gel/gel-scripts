
# R wrapper, main function
read_licor <- function(FilePath) {
	# get file name
	bn <- basename(FilePath)
    # check file name
    if (grepl('[.]qdata$', bn)) {
        return(qs2::qd_read(FilePath))
    } else if (grepl('[.]qs$', bn)) {
        if (!requireNamespace(qs)) {
            stop('data is provided as *.qs file -> install qs library',
                ' running "install.packages("qs")"')
        }
        return(qs::qread(FilePath))
    } else if (grepl('[.]rds$', bn)) {
        return(readRDS(FilePath))
    }
    if (grepl('[.]gz$', bn)) {
        raw_list <- licor_read_cpp_gzip(normalizePath(FilePath))
    } else {
        raw_list <- licor_read_cpp(normalizePath(FilePath))
    }
    out <- data.table::as.data.table(raw_list)
    # convert time
    out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    # remove time_string
    out[, time_string := NULL]
    # place Time first
    setcolorder(out, 'Time')
    # remove NA entries
    na.omit(out)
}
