
library(data.table)
library(ibts)


## 1. ws700 functions ----------------------------------------

ws_read <- function(file_name, Ex = c('E0', 'E2', 'E4')) {
    # get cols
    cols <- unique(unlist(list(
            E0 = paste0('V', c(2, 3, 6, 7, 10:12)),
            E2 = c('V7', 'V11'),
            E4 = c('V3', 'V7')
            )[Ex]))
    # read only valid starting rows
    fread(file = file_name, sep = ';', fill = TRUE, key = 'V2')[(Ex),
            c(list(paste(sub('.*_', '', file_name), V1), V2), mget(cols))]
}

ws_e0 <- function(dat, tz_ws700 = 'Etc/GMT-1') {
    # read E0 (Ta: Tair, Hr: RH, Pa: Pair, Ra: precip sum, Rt: precip type, Ri: precip intens.)
    e0 <- dat[('E0'), {
        Time <- fast_strptime(V1, format = '%Y%m%d %H:%M:%S', tz = tz_ws700, lt = FALSE)
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
    # V7, V11
    e2 <- dat[('E2'), {
        Time <- fast_strptime(V1, format = '%Y%m%d %H:%M:%S', tz = tz_ws700, lt = FALSE)
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

ws_1minute <- function(file_names) {
    ind <- order(as.integer(sub('.*(\\d{8})$', '\\1', file_names)))
    rbindlist(
        lapply(file_names[ind], function(file_name) {
            dat <- ws_read(file_name, Ex = c('E0', 'E2'))
            merge(ws_e0(dat), ws_e2(dat), all = TRUE, by = 'Time')
        })
    )
}

ws_e4 <- function(dat, tz_ws700 = 'Etc/GMT-1') {
    # read E4 (Ca: Compass, Gg: avg gRAD)
    # V3, V7
    e4 <- dat[('E4'), {
        Time <- fast_strptime(V1, format = '%Y%m%d %H:%M:%S', tz = tz_ws700, lt = FALSE)
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

ws_10minute <- function(file_names) {
    ind <- order(as.integer(sub('.*(\\d{8})$', '\\1', file_names)))
    rbindlist(lapply(file_names[ind], function(file_name) {
        dat <- ws_read(file_name, Ex = 'E4')
        ws_e4(dat)
            }))
}

# add function read data from / to
get_ws700 <- function(folder, from = NULL, to = NULL, 
    pooled = c('data.table', 'ibts', '1mins')[1],
    ws_label = NULL, tz_ws700 = 'Etc/GMT-1', 
    keep_recorded = FALSE) {
    # check if folder IS folder or files
    is_folder <- dir.exists(folder)
    is_file <- file.exists(folder) & !is_folder
    # check
    if (any(is_invalid <- !is_folder & !is_file)) {
        stop('path ', paste(folder[is_invalid], collapse = ', '), ' is not valid')
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
    if (any(is_invalid <- !grepl('^ws700-[AB]_\\d{8}$', bnames))) {
        stop('file name ', paste(files[is_invalid], collapse = ', '), ' is not a valid ws700 file name')
    }
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
            latest_time <- ws_1minute(files[n_files])[, Time[.N]]
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
            first_time <- ws_1minute(files[1])[, Time[1]]
            warning('No data available! (No data before "', to, '" - first entry at "', first_time, '")', call. = FALSE)
            return(NULL)
        }
    } else {
        to_index <- n_files
    }
    # get data
    ws1 <- ws_1minute(files[from_index:to_index])
    ws10 <- ws_10minute(files[from_index:to_index])
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
get_ws700A <- function(...) {
    get_ws700(ws_label = 'A', ...)
}
get_ws700B <- function(...) {
    get_ws700(ws_label = 'B', ...)
}

# get_ws700A('test')
# get_ws700A('test', pooled = 'ibts')
# get_ws700A('test', pooled = '1mins')
# get_ws700(c('test/ws700-A_20220813', 'test'), pooled = '10mins')
# get_ws700(c('test/ws700-A_20220813', 'test/ws700-A_20220812'), pooled = '10mins')
# get_ws700(c('test/ws700-A_20220813', 'test/ws700-A_20220812'), from = '15.08.2022 12:00')
# get_ws700(c('test/ws700-A_20220813', 'test/ws700-A_20220812'), to = '10.08.2022 12:00')
# get_ws700(c('test/ws700-A_20220813', 'test/ws700-A_20220812'), from = '10.08.2022 12:00', to = '20.08.2022', pooled = '1mins')
# get_ws700(c('test/ws700-A_20220813', 'test/ws700-A_20220812'), from = '12.08.2022 12:00', to = '12.08.2022 14:00', pooled = '1mins')


## 2. old format ----------------------------------------

# Beispiel ist auskommentiert

# 	PfadDaten <- "~/LFE/01_Projekte/03_FerEVS/Daten"		
# 	if(!dir.exists(PfadDaten)){
# 		# Christoph's Laptop
# 		PfadDaten <- "Y:/01_Projekte/03_FerEVS/Daten"
# 	}
# 	if(!dir.exists(PfadDaten)){
# 		# Marcel's Laptop
# 		PfadDaten <- "//bfh.ch/data/LFE/HAFL/MST-Hofduenger/01_Projekte/03_FerEVS/Daten/"
# 	}


# File1 <- "Wetterstation/OTT_WS1/0000450354_WS700_1_20181127132351.OML"
# File2 <- "Wetterstation/OTT_WS2/0000450354_WS700_2_20181127131340.OML"



gatherSingleOML <- function(x){
	sInt <- sapply(x,function(y)y$channelInfo$storageInterval)
	ids <- sapply(x,function(y)y$channelInfo$channelId)
	name <- sapply(x,function(y)y$channelInfo$name)
	unit <- sapply(x,function(y)y$channelInfo$unit)
	ns <- sapply(x,function(y)nrow(y$data))
	rInt <- rank(sInt, ties.method = "min")
	uRank <- sort(unique(rInt))
	setNames(lapply(uRank,function(i,oml){
		ind <- which(rInt == i)
		if (!all(ns[ind] - ns[ind[1]] == 0)) {
            # browser()
            # stop("data sets are of different lengths (abfangen, falls irgendwann auftreten sollte)")
            mind <- which.max(ns[ind])
            mTime <- oml[[ind[mind]]]$data[, Time]
            oml[ind] <- lapply(oml[ind], function(z) {
                out <- z
                out$data <- merge(data.table(Time = mTime), z$data, all = TRUE)
                out
                })
        }
		Time <- oml[[ind[1]]]$data[,Time]
		list(
			storageInterval = sInt[ind[1]]
			,names = setNames(name[ind], paste0("WS_",ids[ind]))
			,units = setNames(unit[ind], paste0("WS_",ids[ind]))
			,data = do.call(data.table,c(list(Time=Time),setNames(
				lapply(oml[ind],function(z)z$data[,Value]),paste0("WS_",ids[ind]))
				))
			)
	},oml = x), unlist(lapply(sort(sInt)[uRank],
		function(x)trimws(format(make_difftime(x))))))
}


readSingleOML <- function(File){
	# browser()
	## read all lines
	allLines <- suppressWarnings(readLines(File,encoding="UTF-8"))

	# get station info
	SD_i <- grep("^<StationData ",allLines)
	WSname <- sub(".*name=\\\"(.*)\\\" timezone=.*","\\1",allLines[SD_i])
	WStz <- sub(".*timezone=\\\"(.*)\\\" *>","\\1",allLines[SD_i])
	ntz <- unlist(lapply(strsplit(WStz,split = "[:]"),as.numeric))
	stopifnot(ntz[2] == 0)
	Rtz <- sprintf("Etc/GMT%+i",-ntz[1])

	# get channels
	Ch_i <- grep("<ChannelData",allLines)
	V_starts <- grep("<Values>",allLines) + 1
	V_ends <- grep("</Values>",allLines) - 1

	lapply(seq_along(Ch_i),function(i){
		Chinfo <- allLines[Ch_i[i]]
		channelId <- sub(".*channelId=\\\"(.*)\\\" name=.*","\\1",Chinfo)
		name <- sub(".*name=\\\"(.*)\\\" unit=.*","\\1",Chinfo)
		unit <- sub(".*unit=\\\"(.*)\\\" samplingInterval=.*","\\1",Chinfo)
		samplingInterval <- sub(".*samplingInterval=\\\"(.*)\\\" storageInterval=.*","\\1",Chinfo)
		storageInterval <- sub(".*storageInterval=\\\"(.*)\\\" configuredSamplingInterval=.*","\\1",Chinfo)
		# data
		txt <- allLines[V_starts[i]:V_ends[i]]
		Times <- fast_strptime(
			gsub(".*<VT t=\\\"([0-9]{4}[-][0-9]{2}[-][0-9]{2}T[0-9]{2}[:][0-9]{2}[:][0-9]{2})\\\".*","\\1",txt), 
			format = "%Y-%m-%dT%H:%M:%S", lt = FALSE, tz = Rtz)
		Values <- as.numeric(gsub(".*\\\">(.*)</VT>","\\1",txt))
		# check for 999.9 values
		Values[Values == 999.9] <- NA
        isnotna <- !is.na(Times)
		Errorcode <- rep("00",length(Values))
		Errorind <- grep("errorcode",txt)
		Errorcode[Errorind] <- gsub(".*errorcode=\\\"(.*)\\\">.*","\\1",txt[Errorind])
		Values[Errorind] <- NA
		#
		list(
			channelInfo = list(
				channelId = channelId,
				name = name,
				unit = unit,
				samplingInterval = as.numeric(samplingInterval),
				storageInterval = as.numeric(storageInterval)
				),
			data = data.table(
				Time = Times[isnotna],
				Value = Values[isnotna],
				Errorcode = Errorcode[isnotna]
				)
			)
	})
}

readOML <- function(Files, as_ibts = FALSE){
	if(file.info(Files[1])$isdir){
		Files <- list.files(Files[1], full.names=TRUE, pattern = "[.]OML$")
	}
	# read dates from Files
	Dates <- gsub(".*_([0-9]*).OML$","\\1",Files)
	# read all data by sorted date
	for(i in order(Dates)){
		sOML <- readSingleOML(Files[i])
		gOML <- gatherSingleOML(sOML)
		# checks noch einbauen! (Momentan gilt Annahme, dass keine Änderungen an den WS Einstellungen gemacht werden)
		if(i == 1){
			Out <- gOML
			l <- length(gOML)
		} else {
			for(j in seq_len(l)){
				# Zeitcheck!
				Out[[j]]$data <- rbind(Out[[j]]$data,gOML[[j]]$data[Time > TmCheck[[j]],])
			}
		}
		# Zeitcheck!
		TmCheck <- lapply(gOML,function(x){
			x$data[,max(Time)]
		})
	}
	if(as_ibts){
		# Annahme: Zeiten entsprechen Intervall-Ende
		for(j in seq_len(l)){
			# browser()
			# stop("Sind die Zeiten Anfang oder Ende des Intervals?")
			# interval length
			il <- Out[[j]]$storageInterval
			# check act values
			cClasses <- c("avg","num")[as.numeric(grepl("(act)",Out[[j]]$names,fixed = TRUE)) + 1]
			# check sum
			cClasses[grepl("Precip. ",Out[[j]]$names,fixed = TRUE)] <- "sum"
			Out[[j]]$data <- Out[[j]]$data[,{
				as.ibts(.SD, st = Time - il, et = "Time", colClasses = cClasses)
			}]
		}
	}
	
	Out
}


# data1 <- readOML(file.path(PfadDaten,File1))
# data2 <- readOML(file.path(PfadDaten,File2))

# names(data1)

# ##### Temp
# data2[[1]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[1]][,lines(time,value,col="blue")]
# ##### Dew
# data2[[2]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[2]][,lines(time,value,col="blue")]

# ##### RH
# data2[[3]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[3]][,lines(time,value,col="blue")]
# ##### abs hum
# data2[[4]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[4]][,lines(time,value,col="blue")]

# ##### Press
# data2[[5]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[5]][,lines(time,value,col="blue")]
# ##### air dens
# data2[[6]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[6]][,lines(time,value,col="blue")]

# ##### wind speed
# data2[[7]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[7]][,lines(time,value,col="blue")]
# data2[[8]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[8]][,lines(time,value,col="blue")]
# data2[[9]][,plot(time,value,type="l",col="red",ylim=range(value,na.rm=TRUE))]
# data1[[9]][,lines(time,value,col="blue")]

# ### compass
# data1[[12]]
# data2[[12]]
