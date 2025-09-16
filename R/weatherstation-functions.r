
## 1. ws700 functions ----------------------------------------

ws_read <- function(file_name, Ex = c('E0', 'E2', 'E4')) {
    # get cols
    cols <- unique(unlist(list(
            E0 = paste0('V', c(3, 6, 7, 10:12)),
            E2 = c('V7', 'V11'),
            E4 = c('V3', 'V7')
            )[Ex]))
    # read only valid starting rows
    suppressWarnings(raw_lines <- readLines(file_name))
    reg_Ex <- if (length(Ex) > 1) {
        paste0('(',
            paste(Ex, collapse = '|'),
            ')')
    } else {
        Ex
    }
    suppressWarnings(ind <- grep(paste0(
        '^(\\d{4}-\\d{2}-\\d{2}T)?\\d{2}(:\\d{2}){2}([.]\\d{3}Z,|; )', 
        reg_Ex, '(;[^;]+){9}(;[^;]+)?;$'), raw_lines))
    text_in <- if (length(ind) == 1) {
        c(raw_lines[ind], ' ')
    } else {
        raw_lines[ind]
    }
    # check old/new format
    if (is_new <- grepl('py_fnf', file_name)) {
        # fix comma
        text_in <- sub(',', ';', text_in)
        # read data
        fread(text = text_in, sep = ';', fill = TRUE, key = 'V2')[(Ex), 
            c(list(V1, V2), mget(cols))]
    } else {
        # read data
        fread(text = text_in, sep = ';', fill = TRUE, key = 'V2')[(Ex),
            c(list(paste(sub('.*_', '', file_name), V1), V2), mget(cols))]
    }
}

ws_e0 <- function(dat, tz_ws700 = 'Etc/GMT-1') {
    # read E0 (Ta: Tair, Hr: RH, Pa: Pair, Ra: precip sum, Rt: precip type, Ri: precip intens.)
    # units: Tair: °C, RH: %, Pair: hPa, precip sum: mm, precip intens: mm/h
    e0 <- dat[('E0'), {
        if (is.POSIXct(V1)) {
            Time <- with_tz(V1, tz_ws700)
        } else {
            Time <- fast_strptime(V1, format = '%Y%m%d %H:%M:%S', 
                tz = tz_ws700, lt = FALSE)
        }
        c(
            list(
                as.POSIXct(trunc(Time + 30, 'mins')),
                Time
            ), 
            lapply(.SD, function(x) {
                as.numeric(sub(pattern = '^[a-zA-Z]{2}([+-][0-9.]*)[A-Z]$', replacement = '\\1', x))
            })
        )
    }, .SDcols = V3:V12]
    setnames(e0, c('Time', 'recorded.e0', 'T_air', 'RH', 'P_air', 'Prec_sum', 'Prec_type', 'Prec_int'))
    na.omit(e0, cols = 'Time')
}

ws_e2 <- function(dat, tz_ws700 = 'Etc/GMT-1') {
    # read E2 (Sv: U vct, Dv: WD vct)
    # units: U: m/s, WD: °N
    # V7, V11
    e2 <- dat[('E2'), {
        if (is.POSIXct(V1)) {
            Time <- with_tz(V1, tz_ws700)
        } else {
            Time <- fast_strptime(V1, format = '%Y%m%d %H:%M:%S', 
                tz = tz_ws700, lt = FALSE)
        }
        c(
            list(
                as.POSIXct(trunc(Time + 30, 'mins')),
                Time
            ), 
            lapply(.SD, function(x) {
                as.numeric(sub(pattern = '^[a-zA-Z]{2}([+-][0-9.]*)[A-Z]$', replacement = '\\1', x))
            })
        )
    }, .SDcols = c('V7', 'V11')]
    setnames(e2, c('Time', 'recorded.e2', 'U', 'WD'))
    na.omit(e2, cols = 'Time')
}

ws_1minute <- function(file_names, tz_ws700 = 'Etc/GMT-1') {
    if (grepl('^py_fnf', basename(file_names[1]))) {
        ind <- order(as.integer(sub('.*(\\d{4})_(\\d{2})_(\\d{2}).gz$', 
                    '\\1\\2\\3', file_names)))
    } else {
        ind <- order(as.integer(sub('.*(\\d{8})$', '\\1', file_names)))
    }
    rbindlist(
        lapply(file_names[ind], function(file_name) {
            dat <- ws_read(file_name, Ex = c('E0', 'E2'))
            merge(ws_e0(dat, tz_ws700 = tz_ws700), 
                ws_e2(dat, tz_ws700 = tz_ws700), all = TRUE, by = 'Time')
        })
    )
}

ws_e4 <- function(dat, tz_ws700 = 'Etc/GMT-1') {
    # read E4 (Ca: Compass, Gg: avg gRAD)
    # units: Compass: °N, gRAD: W/m2
    # V3, V7
    e4 <- dat[('E4'), {
        if (is.POSIXct(V1)) {
            Time <- with_tz(V1, tz_ws700)
        } else {
            Time <- fast_strptime(V1, format = '%Y%m%d %H:%M:%S', 
                tz = tz_ws700, lt = FALSE)
        }
        c(
            list(
                as.POSIXct(trunc(Time + 30, 'mins')),
                Time
            ), 
            lapply(.SD, function(x) {
                as.numeric(sub(pattern = '^[a-zA-Z]{2}([+-][0-9.]*)[A-Z]$', replacement = '\\1', x))
            })
        )
    }, .SDcols = c('V3', 'V7')]
    setnames(e4, c('Time', 'recorded.e4', 'Compass', 'gRAD'))
    na.omit(e4, cols = 'Time')
}

ws_10minute <- function(file_names, tz_ws700 = 'Etc/GMT-1') {
    if (grepl('^py_fnf', basename(file_names[1]))) {
        ind <- order(as.integer(sub('.*(\\d{4})_(\\d{2})_(\\d{2}).gz$', 
                    '\\1\\2\\3', file_names)))
    } else {
        ind <- order(as.integer(sub('.*(\\d{8})$', '\\1', file_names)))
    }
    rbindlist(lapply(file_names[ind], function(file_name) {
        dat <- ws_read(file_name, Ex = 'E4')
        ws_e4(dat, tz_ws700 = tz_ws700)
            }))
}

# add function read data from / to
read_ws700 <- function(folder, from = NULL, to = NULL, 
    pooled = c('data.table', 'ibts', '1mins')[1],
    ws_label = NULL, tz_ws700 = 'Etc/GMT-1', 
    keep_recorded = FALSE) {
    # check if folder IS folder or files
    is_folder <- dir.exists(folder)
    is_file <- file.exists(folder) & !is_folder
    # check
    if (any(is_invalid <- !is_folder & !is_file)) {
        stop('path ', paste(folder[is_invalid], collapse = ', '), ' cannot be accessed.')
    }
    # add files
    files <- folder[is_file]
    # get files in folder
    if (any(is_folder)) {
        # fix ws_label
        if (!is.null(ws_label)) {
            ws_label <- toupper(ws_label)
            if (!(ws_label %in% c('A', 'B'))) {
                stop('argument "ws_label" must be either "A" or "B"')
            }
        }
        files <- c(files, dir(folder[is_folder], pattern = paste0('ws700-', ws_label), full.names = TRUE))
    }
    # remove duplicates
    files <- unique(files)
    # check file names
    bnames <- basename(files)
    has_new_format <- any(grepl('^py_fnf_0\\d_ws700', bnames))
    if (!has_new_format && any(is_invalid <- !grepl('^ws700-[AB]_\\d{8}$', bnames))) {
        stop('file name ', paste(files[is_invalid], collapse = ', '), ' is not a valid ws700 file name')
    }
    if (!has_new_format) {
        # check ws label consistency
        ulbls <- unique(sub('ws700-(.)_.*', '\\1', bnames))
        if (length(ulbls) != 1) {
            stop('Found both ws700 labels ("A" and "B")!',
                '\n-> Please specify argument "ws_label" to select either')
        }
        # get datetimes & sort
        if (tz_ws700 != 'Etc/GMT-1') stop('Time zone other than "Etc/GMT-1" is not allowed (yet)')
        datetimes <- fast_strptime(sub('ws700-._(\\d{8})$', '\\1', bnames), 
            format = '%Y%m%d', lt = FALSE, tz = tz_ws700)
    } else {
        # loggerbox has UTC
        tz_ws700 <- 'UTC'
        datetimes <- fast_strptime(sub('.*_ws700_(\\d{4})_(\\d{2})_(\\d{2}).gz$', 
                '\\1\\2\\3', bnames), format = '%Y%m%d', lt = FALSE, 
            tz = tz_ws700)
    }
    ind_sort <- order(datetimes)
    datetimes <- datetimes[ind_sort]
    files <- files[ind_sort]
    n_files <- length(files)
    # if from/to missing -> first/last entry
    if (!is.null(from)) {
        from <- parse_date_time3(from, tz = tz_ws700)
        if (!is.POSIXct(from)) {
            stop('Cannot parse argument "from"')
        }
        from_index <- which(datetimes >= trunc(from, units = 'days'))[1]
        if (is.na(from_index)) {
            latest_time <- ws_1minute(files[n_files], tz_ws700 = tz_ws700)[, Time[.N]]
            warning('No data available! (No data after "', from, '" - latest entry at "', latest_time, '")', call. = FALSE)
            return(NULL)
        }
    } else {
        from_index <- 1
    }
    if (!is.null(to)) {
        to <- parse_date_time3(to, tz = tz_ws700)
        if (!is.POSIXct(to)) {
            stop('Cannot parse argument "to"')
        }
        to_index <- which(datetimes <= trunc(to, units = 'days'))
        to_index <- to_index[length(to_index)]
        if (is.na(to_index)) {
            first_time <- ws_1minute(files[1], tz_ws700 = tz_ws700)[, Time[1]]
            warning('No data available! (No data before "', to, '" - first entry at "', first_time, '")', call. = FALSE)
            return(NULL)
        }
    } else {
        to_index <- n_files
    }
    # get data
    ws1 <- ws_1minute(files[from_index:to_index], tz_ws700 = tz_ws700)
    ws10 <- ws_10minute(files[from_index:to_index], tz_ws700 = tz_ws700)
    # fix duplicated times
    ws1 <- unique(ws1, by = 'Time')
    ws10 <- unique(ws10, by = 'Time')
    # strip recorded columns
    if (!keep_recorded) {
        ws1[, c('recorded.e0', 'recorded.e2') := NULL]
        ws10[, recorded.e4 := NULL]
    }
    if (is.null(from)) {
        from <- ws1[, Time[1] - 60]
    }
    if (is.null(to)) {
        to <- ws1[, Time[.N]]
    }
    # add st/et
    setnames(ws1, 'Time', 'et')
    setnames(ws10, 'Time', 'et')
    ws1[, st := et - c(60, pmin(90, as.numeric(diff(et), units = 'secs')))]
    ws10[, st := et - c(600, pmin(900, as.numeric(diff(et), units = 'secs')))]
    # subset to times
    ws1 <- ws1[st >= from & et <= to]
    ws10 <- ws10[st >= from & et <= to]
    setcolorder(ws1, c('st', 'et'))
    setcolorder(ws10, c('st', 'et'))
    # return with switch on argument pooled
    switch(
        pmatch(pooled[1], c('data.table', 'ibts'), nomatch = 3L)
        # data.table -> convert to list with 2 data.table objects
        , list('1mins' = ws1, '10mins' = ws10)
        # ibts -> convert to list with 2 ibts objects
        , list('1mins' = as.ibts(ws1), '10mins' = as.ibts(ws10))
        # pool to one ibts
        , {
            pooled <- parse_time_diff(pooled)
            ws1 <- pool(as.ibts(ws1), pooled, st.to = from, et.to = to)
            ws10 <- pool(as.ibts(ws10), pooled, st.to = from, et.to = to)
            cbind(ws1, ws10)
        }
    )
}
read_ws700A <- function(...) {
    read_ws700(ws_label = 'A', ...)
}
read_ws700B <- function(...) {
    read_ws700(ws_label = 'B', ...)
}

