
#' @title Read Sonic Data
#' @description Main function to read sonic data from various file formats.
#' @param file_path A string specifying the path to the file to be read.
#' @return A data.table containing the processed sonic data or NULL if no valid data is found.
#' @export
## main function to read sonic data
read_sonic <- function(file_path, from = NULL, to = NULL, tz = 'UTC',
    gillbug_method = c('gill', 'nakai2012', 'none')[1]) {
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
            out <- lapply(file.path(file_path, files[read_me]), read_sonic,
                from = from, to = to, tz = tz, gillbug_method = gillbug_method)
            # return sorted
            rbindlist(out)[order(Time)]
        } else {
            stop('No files available within provided time range!')
        }
    } else {
        bn <- basename(file_path)
        # check if provided as qs or rds
        if (grepl('[.]qdata$', bn)) {
            out <- alloc.col(qs2::qd_read(file_path))
        } else if (grepl('[.]qs$', bn)) {
            if (!requireNamespace('qs')) {
                stop('data is provided as *.qs file -> install qs library',
                    ' running "install.packages("qs")"')
            }
            out <- alloc.col(qs::qread(file_path))
        } else if (grepl('[.]rds$', bn)) {
            out <- readRDS(file_path)
        } else if (grepl("^(py_)?fnf_01_", bn)) {
            out <- read_hs_ascii(file_path)
        } else {
            out <- read_windmaster_ascii(file_path)
        }
        # correct for Gill bug
        if (isTRUE(attr(out, 'GillBug'))) {
            switch(match.arg(gillbug_method)
                , gill = {
                    cat("-- Correcting data for Gill software bug. Correction is done as proposed by Gill instruments...\n")
                    out[, w := w * ifelse(w < 0, 1.289, 1.166)]
                    setattr(out, 'GillBug', 'corrected-gill')
                }
                , nakai2012 = {
                    cat("-- Correcting data for Gill software bug. Correction is done as described in Nakai 2012...\n")
                    out[, c('u', 'v', 'w') := nakai_correction_2012(u, v, w)]
                    setattr(out, 'GillBug', 'corrected-nakai2012')
                }
                , none = {
                    warning("-- Data is intentionally NOT corrected for Gill software bug!\n")
                }
                , stop('Gill software bug correction method not valid')
            )
        }
        if (is.null(out)) {
            # not valid
            return(NULL)
        }
        # filter by from/to
        if (length(from) > 0) {
            out <- out[Time >= from]
        }
        if (length(to) > 0) {
            out <- out[Time <= to]
        }
        out
    }
}

#' @title Read New Windmaster ASCII Data
#' @description Reads and processes data from the new Windmaster ASCII format.
#' @param file_path A string specifying the path to the file to be read.
#' @return A data.table containing the processed Windmaster data or NULL if no valid data is found.
# new windmaster format
read_windmaster_ascii <- function(file_path){
	### get Date
	bn <- basename(file_path)
	if(!grepl("^data_|^py_fnf_0", bn)){
		# run old script
		return(read_windmaster_old_ascii(file_path))
	}
    # be verbose
    cat("File:", path.expand(file_path), "- ")
	Date <- gsub("^.*_(20[0-9]{2})_?(\\d{2})_?(\\d{2}).*", "\\1\\2\\3", bn)
	### read File
    raw <- readLines(file_path, warn = FALSE)
    # filter out erroneous multibyte strings
    raw <- raw[grepl('^[^,]+,[^,]+,([0-9.+-]+,){3}M,([0-9.+-]+,){2}[^,]+$', 
        raw, useBytes = TRUE)]
    out <- fread(text = raw, header = FALSE, na.strings = '999.99', 
        showProgress = FALSE, blank.lines.skip = TRUE)
    # check if file is empty
	if(nrow(out) == 0){
        cat('no valid data\n')
		return(NULL)
	}
    # check which columns to convert columns if necessary
    vnums <- paste0('V', c(3, 4, 5, 7))
    is.char <- out[, sapply(.SD, is.character), .SDcols = vnums]
    # convert to numeric
    out[, vnums[is.char] := {
        lapply(.SD, as.numeric)
    }, .SDcols = vnums[is.char]]
    # remove NA lines that come from conversion
    out <- na.omit(out)
    if (out[, .N == 0]) {
        cat('no valid data\n')
        return(NULL)
    }
    # be verbose and print sonic names:
    sonic_label <- out[, sub('^[^A-Z]*', '', unique(V2))]
    sonic_label <- unique(sonic_label[!(sonic_label == '')])
    if (length(sonic_label) == 0) stop('sonic label not available!')
    if (length(sonic_label) > 1) stop('more than one unique sonic label!')
    cat(paste0("data recorded by sonic-", tolower(sonic_label), "\n"))
    if (!(loggerbox <- grepl('^py_fnf_0', bn))) {
        sonic_file <- sub("data_(.*)_[0-9]{8}_[0-9]{6}([.]gz)?$", "\\1", bn)
        if(sonic_label != toupper(sub("sonic-", "", sonic_file))){
            warning(paste0("Sonic label '", sonic_label, "', and hostname '", sonic_file, 
                    "' don't match!"), call. = FALSE)
        }
    }
    # check units
    if(out[, V6[1]] != "M"){
        stop("Units of recorded data not compatible with evaluation script! Column 6 ",
            "should contain 'M' for m/s!")
    }
    # fix time etc.
    out[, c(
        # remove columns
        'V1', 'V2', 'V6', 'V8', 'V9',
        # add columns
        'sonic', 'Time', 'Hz',
        # replace degree Celsius by K
        'V7'
        ) := {
        # set times correctly
        if (loggerbox) {
            st.dec <- V1
        } else {
            st.dec <- fast_strptime(paste(Date, V1), lt = FALSE, 
                format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")
        }
        # get Hz (faster than as.factor(sub('[.].*', '', V1)) !
        Hz <- round(median(tabulate(trunc(as.numeric(st.dec) - 
                        as.numeric(st.dec[1])))), -1)
        if (
            hour(st.dec[.N]) == 0 && 
            (hour(st.dec[1]) != 0 || .N > (3 * Hz))
        ) {
            ## fix last 3 * Hz entries, where hour == 0
            sub.st <- st.dec[.N - seq_len(3 * Hz) + 1]
            hr <- hour(sub.st)
            st.dec[.N - seq_len(3 * Hz) + 1] <- fifelse(hr == 0, sub.st + 24 * 3600, 
                sub.st)
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

#' @title Read Old Windmaster ASCII Data
#' @description Reads and processes data from the old Windmaster ASCII format.
#' @param file_path A string specifying the path to the file to be read.
#' @return A data.table containing the processed Windmaster data or NULL if no valid data is found.
# old windmaster format
read_windmaster_old_ascii <- function(file_path) {
	### get Date
	bn <- basename(file_path)
	Date <- gsub("^..._([0-9]{6})_.*", "\\1", bn)
    ### check if rg is available
    use_rg <- try(system('rg -V', intern = TRUE), silent = TRUE)
    # use_rg <- length(system('command -v rg', intern = TRUE)) > 0
	### read File
	# browser()
    if (inherits(use_rg, 'try-error')) {
        suppressWarnings(out <- fread(cmd = paste0("grep -v -e ',,' -e '[A-Za-z]' '", 
                    path.expand(file_path), "'"), fill = TRUE, blank.lines.skip = TRUE))
    } else {
        suppressWarnings(out <- fread(cmd = paste0("rg -v -e ',,' -e '[A-Za-z]' '", 
                    path.expand(file_path), "'"), fill = TRUE, blank.lines.skip = TRUE))
    }
	if (nrow(out) == 0) {
		cat("File empty:", path.expand(file_path), "\n")
		return(NULL)
	}
	# remove first (empty) column
    if ('V1' %in% names(out)) out[, V1 := NULL]
    # remove NAs
  	out <- na.omit(out)
    ## get Hz
    Hz <- out[, .N, by = V2][, round(median(N), -1)]
    out[, Hz := Hz]
	## set times
	out[, st.dec := fast_strptime(paste0(Date, V2), lt = FALSE, 
        format = "%y%m%d%H:%M:%S", tz = "Etc/GMT-1") + V3][, c("V2", "V3") := NULL]
    ## fix end of day
    if (out[, 
        hour(st.dec[.N]) == 0 && 
        (hour(st.dec[1]) != 0 || .N > (3 * Hz[1]))
        ]) {
        ## fix last 3 * Hz entries, where hour == 0
        out[.N - seq_len(3 * Hz[1]) + 1, st.dec := {
            hr <- hour(st.dec)
            fifelse(hr == 0, st.dec + 24 * 3600, st.dec)
        }]
    }
	### set Output names and order
	setnames(out, c("u", "v", "w", "T", "Hz", "Time"))
	setcolorder(out, c("Time", "Hz", "u", "v", "w", "T"))
	### remove 999.99 entries
	out <- out[!(u %in% 999.99 | v %in% 999.99 | w %in% 999.99 | T %in% 999.99), ]
	### change units from degree Celsius to K
	out[, T := T + 273.15]
    setattr(out, "GillBug", TRUE)
	out
}

#' @title Read HS Sonic ASCII Data
#' @description Reads and processes data from the HS Sonic ASCII format.
#' @param file_path A string specifying the path to the file to be read.
#' @return A data.table containing the processed HS Sonic data or NULL if no valid data is found.
# HS Sonic Agroscope Wauwilermoos
read_hs_ascii <- function(file_path) {
    # be verbose
    cat("File:", path.expand(file_path), "- ")
    if (grepl('[.]gz$', basename(file_path))) {
        # gzipped data
        raw <- hs_read_cpp_gzip(normalizePath(file_path, mustWork = FALSE))
    } else {
        # uncompressed
        raw <- hs_read_cpp(normalizePath(file_path, mustWork = FALSE))
    }
    if (length(raw) == 0) {
        stop('File path: "', file_path, '" is not accessible!')
    }
    out <- data.table::as.data.table(raw)
    out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    out <- na.omit(out)
    cat(paste0("data recorded by HS sonic\n"))
    out <- out[, {
        .(Time, Hz = round(1 / median(as.numeric(diff(Time)))), u = u_string, v = v_string, 
            w = w_string, T = t_string, sonic = 'HS')
    }]
    setattr(out, "GillBug", FALSE)
    out
}
