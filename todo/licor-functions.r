
library(data.table)
library(lubridate)
library(Rcpp)

## get objects before sourcing script
if (.get_diff <- !exists('.licor_functions')) {
    .licor_functions <- ls(all.names = TRUE)
}

## C++ helper function for gzipped licor data
if (Sys.info()['sysname'] == 'Windows') {
    old_pkg_libs <- Sys.getenv('PKG_LIBS')
    # help the linker to include zlib through the environment variable PKG_LIBS
    Sys.setenv(PKG_LIBS = "-lz")       # Fix linking to zlib
}

if (Sys.info()['sysname'] == 'Windows') {
    # reset env var
    Sys.setenv(PKG_LIBS = old_pkg_libs)
}

# R wrapper, main function
read_licor <- function(FilePath) {
	# get file name
	bn <- basename(FilePath)
    # check file name
    if (grepl('[.]qdata$', bn)) {
        if (!require(qs2)) {
            stop('data is provided as *.qdata file -> install qs2 library',
                ' running "install.packages("qs2")"')
        }
        return(qs2::qd_read(FilePath))
    } else if (grepl('[.]qs$', bn)) {
        if (!require(qs)) {
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
    out <- as.data.table(raw_list)
    # convert time
    out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    # remove time_string
    out[, time_string := NULL]
    # place Time first
    setcolorder(out, 'Time')
    # remove NA entries
    na.omit(out)
}

## get objects after sourcing script
if (.get_diff) {
    .licor_functions <- setdiff(ls(all.names = TRUE), .licor_functions)
}
