
# R wrapper, main function
read_licor <- function(file_path) {
	# get file name
	bn <- basename(file_path)
    # check file name
    if (grepl('[.]qdata$', bn)) {
        return(alloc.col(qs2::qd_read(file_path)))
    } else if (grepl('[.]qs$', bn)) {
        if (!requireNamespace(qs)) {
            stop('data is provided as *.qs file -> install qs library',
                ' running "install.packages("qs")"')
        }
        return(alloc.col(qs::qread(file_path)))
    } else if (grepl('[.]rds$', bn)) {
        return(readRDS(file_path))
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
    # place Time first
    setcolorder(out, 'Time')
    # remove NA entries
    na.omit(out)
}
