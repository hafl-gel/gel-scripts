

library(data.table)
library(ibts)

select_files <- function(ht_path, sonic_path, from, to) {
    # read ht folder
    if (!file.exists(ht_path)) {
        stop('ht_path "', ht_path, '" is not accessible!')
    } else if (length(ht_path) == 1 && file.info(ht_path)$isdir) {
        ht_files <- dir(ht_path, pattern = '^ht8700_sonic-.')
    } else {
        ht_files <- basename(ht_path)
        ht_path <- dirname(ht_path)
    }
    if (length(ht_files) == 0) {
        stop('No HT8700 data in directory "', ht_path, '"')
    }
    # read sonic folder
    if (!file.exists(sonic_path)) {
        stop('sonic_path "', sonic_path, '" is not accessible!')
    } else if (length(sonic_path) == 1 && file.info(sonic_path)$isdir) {
        sonic_files <- dir(sonic_path, pattern = '^data_sonic-.')
    } else {
        sonic_files <- basename(sonic_path)
        sonic_path <- dirname(sonic_path)
    }
    if (length(sonic_files) == 0) {
        stop('No sonic data in directory "', sonic_path, '"')
    }
    # get dates
    sonic_dates <- as.integer(sub('.*_(\\d{8})_.*', '\\1', sonic_files))
    # check from
    if (is.null(from) || (is.character(from) && from[1] == 'first')) {
        dfrom <- min(sonic_dates)
    } else {
        # fix to/from to Etc/GMT-1 (data time zone)
        from <- parse_date_time3(from, tz = 'Etc/GMT-1')
        # create date int
        dfrom <- as.integer(unique(format(from, '%Y%m%d')))
    }
    # check to
    if (is.null(to) || (is.character(to) && to[1] == 'last')) {
        dto <- max(sonic_dates)
    } else {
        # fix to to Etc/GMT-1 (data time zone)
        to <- parse_date_time3(to, tz = 'Etc/GMT-1')
        # create date int
        dto <- as.integer(unique(format(to, '%Y%m%d')))
    }
    # get subset index based on time range
    if (length(dfrom) != length(dto)) {
        stop('arguments "from" and "to" must have equal lengths!')
    } else if (length(dfrom) > 1) {
        takeme <- logical(length(sonic_dates))
        takeme[
            unlist(mapply(\(xf, xt) which(sonic_dates >= xf & sonic_dates <= xt),
                xf = dfrom, xt = dto, SIMPLIFY = FALSE))
            ] <- TRUE
    } else {
        takeme <- sonic_dates >= dfrom & sonic_dates <= dto
    }
    # no data within time range
    if (!any(takeme)) {
        verbose <- dynGet('verbose', 1)
        if (verbose > 0) cat('No data for specified date range...\n')
        return(list())
    }
    # get sonic label
    sonic <- unique(sub('^data_(sonic-.)_.*', '\\1', sonic_files[takeme]))
    if (length(sonic) > 1) {
        stop('Folder "', path_sonic, '" contains files from different sonic computers',
            ' for the given time range!')
    }
    # find matches
    #     -> check sonic
    ht_sonic <- grep(sonic, ht_files, fixed = TRUE, value = TRUE)
    ht_dates <- as.integer(sub('.*_(\\d{8})_.*', '\\1', ht_sonic))
    #     -> check date
    dates_ok <- ht_dates %in% sonic_dates[takeme]
    unique_dates <- as.character(unique(ht_dates[dates_ok]))
    # select sonic files
    sonic_ok <- lapply(unique_dates, \(pat) file.path(
            sonic_path,
            grep(pat, x = sonic_files, fixed = TRUE, value = TRUE)
            ))
    names(sonic_ok) <- unique_dates
    # select ht files
    ht_ok <- lapply(unique_dates, \(pat) file.path(
            ht_path,
            grep(pat, x = ht_sonic, fixed = TRUE, value = TRUE)
            ))
    names(ht_ok) <- unique_dates
    # return list of files
    list(
        dates = unique_dates,
        sonic = sonic_ok,
        ht8700 = ht_ok
    )
}

ht_merged_daily <- function(path_ht = NULL, path_sonic = NULL, path_merged = NULL,
    from_date = NULL, to_date = NULL, save_file = FALSE, sonic_basis = TRUE,
    tzone_out = 'Etc/GMT-1', verbose = 1) {
    # get files & dates
    files <- select_files(path_ht, path_sonic, from = from_date, to = to_date)
    # check merged file
    merged_files <- character(0)
    if (!is.null(path_merged)) {
        merged_files <- dir(path_merged, pattern = '^ht_merged_')
    }
    merged_dates <- sub('^ht_merged_.*Hz_sonic-._(\\d{8})_.*', '\\1', merged_files)
    out <- lapply(files[['dates']], \(u_date) {
        if (verbose > 0) cat('merging data from', sub('(\\d{4})(\\d{2})(\\d{2})', '\\3.\\2.\\1', u_date), '')
        u_sonic_files <- files[['sonic']][[u_date]]
        u_ht_files <- files[['ht8700']][[u_date]]
        # gen hash all sonic and ht files
        u_sha_full <- sha1(list(
            # names
            basename(u_sonic_files),
            basename(u_ht_files),
            # sizes
            file.size(u_sonic_files),
            file.size(u_ht_files)
        ))
        u_sha <- substr(u_sha_full, 1, 14)
        m_file <- m_sha <- ''
        # check if existing file hash matches
        if (u_date %in% merged_dates) {
            m_file <- merged_files[merged_dates %in% u_date]
            # -> hash in filename?
            m_sha <- sub('.*_([^_]*)$', '\\1', m_file)
        }
        # identical?
        if (u_sha == m_sha) {
            if (verbose > 0) cat('- merged file is up-to-date...\n')
            m_data <- fread(file.path(path_merged, m_file))
        } else {
            # create new merged file
            if (verbose > 1) {
                cat(':\n')
                vfoo <- eval
            } else {
                if (verbose > 0) {
                    if (m_sha != '') {
                        cat('- updating merged file - ')
                    } else if(!is.null(path_merged)) {
                        cat('- creating merged file - ')
                    } else {
                        cat('- ')
                    }
                }
                vfoo <- capture.output
            }
            vfoo({
                # read sonic
                sonic_data <- rbindlist(lapply(u_sonic_files, read_windmaster_ascii))
                # read ht
                ht_data <- rbindlist(lapply(u_ht_files, read_ht8700))
            })
            # merge data
            if (sonic_basis) {
                m_data <- merge_data(sonic_data, ht_data)
            } else {
                m_data <- merge_data(ht_data, sonic_data)
            }
            Hz <- m_data[, round(1 / median(diff(as.numeric(Time, units = 'secs'))), -1)]
            # save to disk
            if (save_file && !is.null(path_merged)) {
                fwrite(m_data[, 
                    .(Time, nh3_ppb, nh3_ugm3, temp_amb, press_amb, oss, u, v, w, T, sonic)
                    ], file = paste0(path_merged, '/ht_merged_', Hz, 'Hz_', sonic, '_',
                        u_date, '_', u_sha)
                )
            }
            if (verbose == 1) cat('done\n')
        }
        m_data
    })
    names(out) <- files[['dates']]
    out <- rbindlist(out)
    # fix tzone
    if (nrow(out) > 0) {
        out[, Time := with_tz(Time, tzone_out)]
    }
    out[]
}


# read HT8700 data
read_ht8700 <- function(FilePath, tz = "Etc/GMT-1"){
	### get Date
	bn <- basename(FilePath)
	if(!grepl("^ht8700_", bn)){
        stop('data filename not valid')
	}
    if (!any(grepl('R.utils', installed.packages()[, 'Package']))) {
        stop('package "R.utils" must be installed to process gz files!')
    }
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	Date <- gsub("^ht8700_.*_([0-9]{8})_.*", "\\1", bn)
	### read File
    raw <- readLines(FilePath, warn = FALSE)
    # filter out erroneous multibyte strings
    raw <- raw[grepl('^\\d{2}[0-9.:]+,([ A-Z0-9.+-]+,){18}[[ 0-9.+-]+$', raw, useBytes = TRUE)]
    if (length(raw) == 1) {
        raw <- c(raw, '')
    }
    # read from string
    out <- fread(text = raw, blank.lines.skip = TRUE,
        header = FALSE, na.strings = '999.99', showProgress = FALSE)
    # check empty
    if (nrow(out) == 0) {
        cat('no valid data!\n')
        return(NULL)
    }
    # convert column types
    char_cols <- paste0('V', c(1, 2, 17, 18, 20))
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
    # fix time etc.
    out[, Time := fast_strptime(paste(Date, V1), lt = FALSE, format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")]
    # fix column names
    setnames(out,
        c(
            'tstring', # column 1
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
            'Ref_Road_2F', # column 16
            'alarm_lower_bit', # column 17
            'alarm_upper_bit', # column 18
            'cleaning_flag', # column 19
            'notused', # column 20
            'Time' # column 21
        )
    )
    cat('done\n')
    # return
    out
}

# merge sonic & ht8700
library(Rcpp)

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List match_times(NumericVector time1, NumericVector time2, double deltat)
{
    int len1 = time1.size();
    int len2 = time2.size();
    // iterator for time2
    int j = 0;
    // time diffs
    double dj;
    double dj_1;
    // loop over time1 only
	IntegerVector index1 = seq_len(len1);
    IntegerVector index2(len1);
    // loop over time1
    for (int i = 0; i < len1; i++) {
        // advance time2 index until time2[j] > time1[i] or j == len2 - 1
        while (time2[j] < time1[i] && j < len2 - 1) {
            j++;
        }
        // check j
        dj = time2[j] - time1[i];
        if (j == 0) {
            if (dj <= deltat) {
                index2[i] = j + 1;
            }
        } else {
            // check j and j-1
            dj_1 = time1[i] - time2[j - 1];
            if (dj_1 <= dj && std::abs(dj_1) <= deltat) {
                index2[i] = j;
            } else if (std::abs(dj) <= deltat) {
                index2[i] = j + 1;
            }
        }
        // check if index2[i] is 0
        if (index2[i] == 0) {
            index1[i] = 0;
        }
    }
    return List::create(index1, index2);
}
')

# merge sonic & ht8700 data based on time
#   -> output contains the same times as 'basis'
#   -> values from 'draw' will be repeated or dropped to match 'basis' times
merge_data <- function(basis, draw) {
    if (nrow(basis) == 0 || nrow(draw) == 0) {
        return(
            data.table(
                Time = character(0), 
                nh3_ppb = character(0),
                nh3_ugm3 = character(0), 
                temp_amb = character(0), 
                press_amb = character(0), 
                oss = character(0), 
                u = character(0), 
                v = character(0), 
                w = character(0), 
                T = character(0), 
                sonic = character(0)
            )
        )
    }
    # times
    t_basis <- basis[, as.numeric(Time)]
    t_draw <- draw[, as.numeric(Time)]
    t0 <- t_basis[1]
    # ~ 1/Hz
    d_t <- median(diff(t_basis))
    # get matching indices
    indices <- match_times(t_basis - t0, t_draw - t0, d_t)
    # bind together
    out <- cbind(basis[indices[[1]]], draw[indices[[2]]][, Time := NULL])
    # add Hz
    out[, Hz := round(1 / d_t, -1)]
    out
}

# add sonic data to ht8700 data
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

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List decal(IntegerVector x, IntegerVector y)
{
    int len1 = x.size();
    IntegerVector index2(len1);
    List out = List(len1);
    LogicalMatrix mat(len1, 32);
    // loop over x
    for (int i = 0; i < len1; i++) {
        // create empty vector
        IntegerVector l;
        if (y[i] != 0) {
            // loop
            int a = y[i];
            for (int j = 15; j >= 0; j--) {
                int p = std::pow(2, j);
                if (a >= p) {
                    // add value of j + 1 + 16
                    l.push_front(j + 17);
                    // update mat
                    mat(i, j + 16) = TRUE;
                    // update a
                    a = a % p;
                }
            }
        }
        if (x[i] == 0) {
            out[i] = 0;
        } else {
            // loop
            int a = x[i];
            for (int j = 15; j >= 0; j--) {
                int p = std::pow(2, j);
                if (a >= p) {
                    // add value of j + 1
                    l.push_front(j + 1);
                    // update mat
                    mat(i, j) = TRUE;
                    // update a
                    a = a % p;
                }
            }
        }
        // assign vector to list entry
        if (l.size() == 0) {
            l = 0;
        }
        out[i] = l;
    }
    // add matrix attribute
    out.attr("mat") = mat;
    return out;
}
')


# NOTE:
# decode_alarm() could be done in R starting from:
# decode_alarm <- function(lower, upper) {
#     which(rev(unlist(strsplit(paste0(
#         # upper
#         sprintf('%016i', as.integer(R.utils::intToBin(strtoi(upper, 16L)))),
#         # lower
#         sprintf('%016i', as.integer(R.utils::intToBin(strtoi(lower, 16L))))
#         ), split = ''))) == '1')
# }
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
# simple: TRUE only provide 1 column with codes, FALSE 1 column per code containing TRUE/FALSE values
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


# issue with git due to large files => store only latest version locally,
# other versions are stored on a remote file system
# save file locally with correct name
# save with hash as name on "remote" drive
# add 2 functions to read (incl. check) and write
# add sha1 as attribute to object
# library to serialize R data to files
library(qs)
# library to hash R object
library(digest)

# main function to write "hashed" files
write_hashed <- function(x, file_path, remote_path = getOption('remote_path'), 
    save_local = TRUE, ...) {
    # get file name
    file_name <- basename(file_path)
    # add .hashed to filename
    file_name <- sub('([.]hashed)?$', '.hashed', file_name)
    # get local path
    local_path <- dirname(file_path)
    # check directories
    local_path <- check_path(local_path)
    if (is.null(remote_path)) {
        stop('argument "remote_path" is missing!')
    }
    remote_path <- check_path(remote_path)
    # add hash directory to local
    dir.create(file.path(local_path, '.hash'), showWarnings = FALSE)
    # add file hash
    add_hash(x)
    # save file on remote
    cat('writing remote file...\n')
    qs::qsave(x, file.path(remote_path, get_hash(x)), preset = 'fast')
    # save hash locally
    cat('updating hash file...\n')
    writeLines(get_hash(x), file.path(local_path, '.hash', file_name))
    # save file locally
    if (save_local) {
        cat('writing local file...\n')
        write_local(x, file.path(local_path, file_name), ...)
    }
    invisible(TRUE)
}
# main function to read "hashed" files
read_hashed <- function(file_path, remote_path = getOption('remote_path'), 
    update_local = TRUE) {
    # add .hashed to filename
    file_path <- sub('([.]hashed)?$', '.hashed', file_path)
    # get local path
    local_path <- dirname(file_path)
    # check local directory
    local_path <- check_path(local_path)
    # get file name
    file_name <- basename(file_path)
    # check if .hash folder exists
    if (dir.exists(file.path(local_path, '.hash'))) {
        cat('checking hash...\n')
        # get repo hash
        hash_repo <- readLines(file.path(local_path, '.hash', file_name), n = 1L)
        # get hash from local file: NULL -> file is missing
        hash_file <- read_local(file_path, hash_only = TRUE)
        if (is.null(hash_file)) {
            cat('-> local file does not exist\n')
        } else if (hash_repo == hash_file) {
            cat('-> local file is up-to-date\n')
        } else {
            cat('-> local file needs to be updated\n')
        }
        # check file status
        read_remote <- !update_local || is.null(hash_file) || hash_repo != hash_file
    } else {
        if (update_local) {
            warning('.hash folder is missing. Cannot check local file status.')
        } else {
            stop('.hash folder is missing. Cannot read remote file.')
        }
        read_remote <- FALSE
    }
    # update local file
    if (read_remote) {
        if (is.null(remote_path)) {
            stop('argument "remote_path" is missing!')
        }
        remote_path <- check_path(remote_path)
        # read with qread from remote
        cat('reading from remote...\n')
        out <- alloc.col(qread(file.path(remote_path, hash_repo)))
        # save to local system
        if (update_local) {
            cat('updating local file from remote...\n')
            write_local(out, file_path)
        }
        # return remote
        out
    } else {
        # return local
        cat('reading local file...\n')
        read_local(file_path)
    }
}

## helper functions
# sha1 from list of first time, last time, .N, names
hash <- function(x) {
    x[, digest::sha1(
        list(
            # start
            Time[1],
            # end
            Time[.N],
            # .N
            .N,
            # names
            names(x)
        )
    )]
}
add_hash <- function(x) setattr(x, 'hash', hash(x))
get_hash <- function(x) attr(x, 'hash')
check_path <- function(path) {
    path <- normalizePath(path, mustWork = FALSE)
    if (!dir.exists(path)) {
        stop("Can't access directory \"", path, '"')
    }
    path
}
write_local <- function(dat, path, ...) {
    dat_ser <- qs::qserialize(dat, ...)
    con <- file(path, open = 'wb')
    on.exit(close(con))
    # write hash
    writeBin(get_hash(dat), con)
    # write number of bytes
    writeBin(length(dat_ser), con)
    # write serialized dat
    writeBin(dat_ser, con)
}
read_local <- function(path, hash_only = FALSE) {
    # check if file exists
    if (!file.exists(path)) return(NULL)
    # read local file
    con <- file(path, open = 'rb')
    on.exit(close(con))
    # read hash
    hash_file <- readBin(con, 'character')
    if (hash_only) return(hash_file)
    # read number of bytes
    n_bytes <- readBin(con, 'integer')
    # read serialized data
    out_ser <- readBin(con, 'raw', n = n_bytes)
    alloc.col(qs::qdeserialize(out_ser))
}
