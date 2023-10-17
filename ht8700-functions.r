

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

# read_windmaster_ascii("~/repos/3_Scripts/5_shellSonic/test_data/data_sonicb_20210120_170910")

# TODO: read HT8700 data
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
    out <- try(
        fread(FilePath, encoding = 'UTF-8', header = FALSE, fill = TRUE, 
            blank.lines.skip = TRUE, na.strings = '999.99', 
            # select = 1:9,
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
    # out <- out[grepl('^[\x01-\x1A]', V9)]
    # check which columns to convert columns if necessary
    # vnums <- paste0('V', c(3, 4, 5, 7))
    # is.char <- out[, sapply(.SD, is.character), .SDcols = vnums]
    # convert to numeric
    # if (any(is.char)) {
    #     op <- getOption('show.error.messages')
    #     options(show.error.messages = FALSE)
    #     on.exit(options(show.error.messages = op))
    #     check <- try(out[, vnums[is.char] := {
    #         lapply(.SD, as.numeric)
    #     }, .SDcols = vnums[is.char]])
    #     if (inherits(check, 'try-error')) {
    #         # remove multibyte string
    #         ind <- out[, Reduce('&', lapply(.SD, grepl, pattern = '^[-+]?[0-9]+[.][0-9]+$')), .SDcols = vnums[is.char]]
    #         out <- out[ind, ]
    #         # try to convert again
    #         out[, vnums[is.char] := {
    #             lapply(.SD, as.numeric)
    #         }, .SDcols = vnums[is.char]]
    #     }
    #     options(show.error.messages = op)
    #     on.exit()
    # }
    # remove NA lines that come from conversion
    out <- na.omit(out)
    # be verbose and print sonic names:
    # sonic_label <- out[, sub('[\x01-\x1A]', '', unique(V2))]
    # sonic_label <- sonic_label[!(sonic_label == '')]
    # if (length(sonic_label) == 0) stop('sonic label not available!')
    # if (length(sonic_label) > 1) stop('more than one unique sonic label!')
    # cat(paste0("data recorded by sonic-", tolower(sonic_label), "\n"))
    # sonic_file <- sub("data_(.*)_[0-9]{8}_[0-9]{6}([.]gz)?$", "\\1", bn)
    # if(sonic_label != toupper(sub("sonic-", "", sonic_file))){
    #     warning(paste0("Sonic label '", sonic_label, "', and hostname '", sonic_file, "' don't match!"), call. = FALSE)
    # }
    # check units
    # if(out[, V6[1]] != "M"){
    #     stop("Units of recorded data not compatible with evaluation script! Column 6 should contain 'M' for m/s!")
    # }
    # fix time etc.
    out[, Time := fast_strptime(paste(Date, V1), lt = FALSE, format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")]
    # out[, c(
    #     # remove columns
    #     'V1', 'V2', 'V6', 'V8', 'V9',
    #     # add columns
    #     'sonic', 'Time', 'Hz',
    #     # replace °C by K
    #     'V7'
    #     ) := {
    #     # set times correctly
    #     st.dec <- fast_strptime(paste(Date, V1), lt = FALSE, format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")
    #     # get Hz (faster than as.factor(sub('[.].*', '', V1)) !
    #     Hz <- round(median(tabulate(trunc(as.numeric(st.dec) - as.numeric(st.dec[1])))), -1)
    #     if (
    #         hour(st.dec[.N]) == 0 && 
    #         (hour(st.dec[1]) != 0 || .N > (3 * Hz))
    #     ) {
    #         ## fix last 3 * Hz entries, where hour == 0
    #         sub.st <- st.dec[.N - seq_len(3 * Hz) + 1]
    #         hr <- hour(sub.st)
    #         st.dec[.N - seq_len(3 * Hz) + 1] <- fifelse(hr == 0, sub.st + 24 * 3600, sub.st)
    #     }
    #     # return list
    #     list(
    #         # remove columns
    #         NULL, NULL, NULL, NULL, NULL,
    #         # sonic
    #         sonic_label,
    #         # Time
    #         st.dec,
    #         # Hz
    #         Hz,
    #         # T
    #         V7 + 273.15
    #         )
    # }]
    ### set Output names and order
    # setnames(out, c("u", "v", "w", "T", "sonic", "Time","Hz"))
    # setcolorder(out,c("Time","Hz","u","v","w","T", "sonic"))
    # return
    out
}

# file_path <- '~/repos/3_Scripts/5_shellSonic/test_data/ht8700_sonic-a_20231016_155256.gz'
# file_path <- '~/repos/3_Scripts/5_shellSonic/test_data/ht8700_sonic-a_20231017_000000.gz'
x <- read_ht8700(file_path)
# sonic_path <- '~/repos/3_Scripts/5_shellSonic/test_data/data_sonic-a_20231017_000000.gz'
y <- read_windmaster_ascii(sonic_path)

par(mfrow = c(3, 2))
x[, plot(Time, V3, ylim = c(300, 400))]
x[, plot(Time, V4)]
x[, plot(Time, V6)]
x[, plot(Time, V7)]
y[, plot(Time, T - 273.15, ylim = c(14, 17) + 0.5)]

