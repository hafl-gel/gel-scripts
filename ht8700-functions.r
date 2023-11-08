

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
        stop('sonic data not valid')
	}
    if (grepl('[.]gz$', bn)) {
        if (!any(grepl('R.utils', installed.packages()[, 'Package']))) {
            stop('package "R.utils" must be installed to process gz files!')
        }
    }
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	Date <- gsub("^ht8700_.*_([0-9]{8})_.*", "\\1", bn)
	### read File
    suppressWarnings(out <- try(
        fread(FilePath, encoding = 'UTF-8', header = FALSE, fill = TRUE, 
            sep = ',', 
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
    # check column number
    if (ncol(out) > 20) {
        # remove extra columns
        suppressWarnings(out <- out[, V1:V20])
        # check for id in second column
        suppressWarnings(out <- out[grepl('\\d{4}', V2)])
        # check for 0 in V20
        suppressWarnings(out <- out[as.integer(V20) == 0L])
        # check for 0 or 1 in V19
        suppressWarnings(out <- out[as.integer(V19) %in% c(0L, 1L)])
        if (out[, .N] == 0) {
            cat('file contains no valid data!\n')
            return(NULL)
        }
    } else if (out[, .N == 0]) {
        cat('file empty\n')
        return(NULL)
    } else if (ncol(out) < 20) {
        cat('file contains no valid data!\n')
        return(NULL)
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
    out <- na.omit(out, cols = 1:20)
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

