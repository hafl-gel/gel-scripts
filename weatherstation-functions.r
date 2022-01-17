
library(data.table)
library(lubridate)

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
