

library(data.table)
library(ibts)


# read windmaster data
read_windmaster_ascii <- function(FilePath, tz = "Etc/GMT-1"){
	### get Date
	bn <- basename(FilePath)
	if(!grepl("^data_", bn)){
        stop('sonic data not valid')
	}
    if (grepl('[.]gz$', bn)) {
        if (!any(grepl('R.utils', installed.packages()[, 'Package']))) {
            stop('package "R.utils" must be installed to process gz files!')
        }
    }
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	Date <- gsub("^data_.*_([0-9]{8})_.*", "\\1", bn)
	### read File
    out <- try(
        fread(FilePath, encoding = 'UTF-8', header = FALSE, fill = TRUE, 
            blank.lines.skip = TRUE, na.strings = '999.99', select = 1:9,
            showProgress = FALSE), 
        silent = TRUE)
    # check if file is empty
	if(inherits(out, 'try-error')){
        if (any(grepl('File is empty', out))) {
            cat('empty\n')
        } else {
            cat('error reading file\n')
            cat(out)
        }
		return(NULL)
	}
    # select valid V9 only
    out <- out[grepl('^[\x01-\x1A]', V9)]
    # check which columns to convert columns if necessary
    vnums <- paste0('V', c(3, 4, 5, 7))
    is.char <- out[, sapply(.SD, is.character), .SDcols = vnums]
    # convert to numeric
    if (any(is.char)) {
        op <- getOption('show.error.messages')
        options(show.error.messages = FALSE)
        on.exit(options(show.error.messages = op))
        check <- try(out[, vnums[is.char] := {
            lapply(.SD, as.numeric)
        }, .SDcols = vnums[is.char]])
        if (inherits(check, 'try-error')) {
            # remove multibyte string
            ind <- out[, Reduce('&', lapply(.SD, grepl, pattern = '^[-+]?[0-9]+[.][0-9]+$')), .SDcols = vnums[is.char]]
            out <- out[ind, ]
            # try to convert again
            out[, vnums[is.char] := {
                lapply(.SD, as.numeric)
            }, .SDcols = vnums[is.char]]
        }
        options(show.error.messages = op)
        on.exit()
    }
    # remove NA lines that come from conversion
    out <- na.omit(out)
    if (out[, .N == 0]) {
        cat('no valid data\n')
        return(NULL)
    }
    # be verbose and print sonic names:
    sonic_label <- out[, sub('[\x01-\x1A]', '', unique(V2))]
    sonic_label <- sonic_label[!(sonic_label == '')]
    if (length(sonic_label) == 0) stop('sonic label not available!')
    if (length(sonic_label) > 1) stop('more than one unique sonic label!')
    cat(paste0("data recorded by sonic-", tolower(sonic_label), "\n"))
    sonic_file <- sub("data_(.*)_[0-9]{8}_[0-9]{6}([.]gz)?$", "\\1", bn)
    if(sonic_label != toupper(sub("sonic-", "", sonic_file))){
        warning(paste0("Sonic label '", sonic_label, "', and hostname '", sonic_file, "' don't match!"), call. = FALSE)
    }
    # check units
    if(out[, V6[1]] != "M"){
        stop("Units of recorded data not compatible with evaluation script! Column 6 should contain 'M' for m/s!")
    }
    # fix time etc.
    out[, c(
        # remove columns
        'V1', 'V2', 'V6', 'V8', 'V9',
        # add columns
        'sonic', 'Time', 'Hz',
        # replace °C by K
        'V7'
        ) := {
        # set times correctly
        st.dec <- fast_strptime(paste(Date, V1), lt = FALSE, format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")
        # get Hz (faster than as.factor(sub('[.].*', '', V1)) !
        Hz <- round(median(tabulate(trunc(as.numeric(st.dec) - as.numeric(st.dec[1])))), -1)
        if (
            hour(st.dec[.N]) == 0 && 
            (hour(st.dec[1]) != 0 || .N > (3 * Hz))
        ) {
            ## fix last 3 * Hz entries, where hour == 0
            sub.st <- st.dec[.N - seq_len(3 * Hz) + 1]
            hr <- hour(sub.st)
            st.dec[.N - seq_len(3 * Hz) + 1] <- fifelse(hr == 0, sub.st + 24 * 3600, sub.st)
        }
        # return list
        list(
            # remove columns
            NULL, NULL, NULL, NULL, NULL,
            # sonic
            sonic_label,
            # Time
            st.dec,
            # Hz
            Hz,
            # T
            V7 + 273.15
            )
    }]
    ### set Output names and order
    setnames(out, c("u", "v", "w", "T", "sonic", "Time","Hz"))
    setcolorder(out,c("Time","Hz","u","v","w","T", "sonic"))
    # check if GillBug affected sonic
    setattr(out, 'GillBug', sonic_label %in% c('C', 'D'))
    # return
    out
}

# read HT8700 data
read_ht8700 <- function(FilePath, tz = "Etc/GMT-1"){
	### get Date
	bn <- basename(FilePath)
	if(!grepl("^ht8700_", bn)){
        stop('data filename not valid')
	}
    # check if zcat and awk exist
    has_awk <- length(system('type awk || echo ""', intern = TRUE)) == 1
    has_zcat <- length(system('type zcat || echo ""', intern = TRUE)) == 1
    if (!has_awk || !has_zcat) {
        warning('consider installing',
            switch(paste(has_awk, has_zcat),
                'FALSE FALSE' = 'awk and zcat',
                'FALSE TRUE' = 'awk',
                'TRUE FALSE' = 'zcat'
            ), 'to faster process raw data files!\n', call. = FALSE
        )
    }
    # check if zipped
    if (is_zipped <- grepl('[.]gz$', bn)) {
        # check if zcat and awk exist
        if (!has_awk || !has_zcat) {
            # if not switch to read lines & process afterwards
            #   i.e. filter for 19 commas and fread again (slower?)
            if (!any(grepl('R.utils', installed.packages()[, 'Package']))) {
                stop('package "R.utils" must be installed to process gz files!')
            }
        }
    }
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	Date <- gsub("^ht8700_.*_([0-9]{8})_.*", "\\1", bn)
	### read File
    if (has_awk && (!is_zipped || (is_zipped && has_zcat))) {
        # awk solution
        if (is_zipped) {
            cat_cmd <- 'zcat'
        } else {
            cat_cmd <- 'cat'
        }
        # awk command incl. filtering for
        # 20 columns, field 20 must be 0, field 19 either 0 or 1,
        # field 2 should be 4 digits
        # add more checks
        awk_cmd <- 'NF == 20 && $2 ~ /^\\s*[0-9]{4}\\s*$/ && $19 ~ /^[01]$/ && $20 == 0'
        # read file
        suppressWarnings(out <- try(fread(
            cmd = paste(cat_cmd, FilePath, '| awk -F, \'', awk_cmd, '\' -'),
            header = FALSE, encoding = 'UTF-8', showProgress = FALSE
        ), silent = TRUE))
        # check if file is empty (is this necessary with fread/zcat/awk?)
        if(inherits(out, 'try-error')){
            if (any(grepl('File is empty', out))) {
                cat('empty\n')
            } else {
                cat('error reading file\n')
                cat(out)
            }
            return(NULL)
        }
        if (out[, .N] == 0) {
            cat('file contains no valid data!\n')
            return(NULL)
        }
    } else {
        # no awk
        suppressWarnings(out <- try(
            fread(FilePath, encoding = 'UTF-8', header = FALSE, fill = TRUE, 
                sep = '\n', 
                blank.lines.skip = TRUE, showProgress = FALSE), 
            silent = TRUE))
        # check if file is empty
        if(inherits(out, 'try-error')){
            if (any(grepl('File is empty', out))) {
                cat('empty\n')
            } else {
                cat('error reading file\n')
                cat(out)
            }
            return(NULL)
        }
        # filter for exactly 19 commas and check for id in second column
        out <- out[grepl('^[^,]*,\\s*\\d{4}\\s*,([^,]*,){17}[^,]*$', V1)]
        # fread again
        out[, paste0('V', 1:20) := fread(text = V1)]
        # suppressWarnings(out <- out[grepl('\\d{4}', V2)])
        # check for 0 in V20 and 0/1 in V19
        suppressWarnings(out <- out[as.integer(V20) == 0L & as.integer(V19) %in% c(0L, 1L)])
        if (out[, .N] == 0) {
            cat('file contains no valid data!\n')
            return(NULL)
        }
    }
    # convert column types
    char_cols <- paste0('V', c(1, 2, 17, 18, 20))
    num_cols <- paste0('V', 3:16)
    suppressWarnings(out[, (char_cols) := lapply(.SD, as.character), .SDcols = char_cols])
    suppressWarnings(out[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols])
    suppressWarnings(out[, V19 := as.integer(V19)])
    # check hexmode ok lower bit
    suppressWarnings(out <- out[!is.na(strtoi(V17, 16L))])
    # check hexmode ok upper bit
    suppressWarnings(out <- out[!is.na(strtoi(V18, 16L))])
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

# helper function to decode alarm codes
decode_alarm <- function(lower, upper) {
    lo <- as.integer(as.hexmode(lower))
    up <- as.integer(as.hexmode(upper))
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
