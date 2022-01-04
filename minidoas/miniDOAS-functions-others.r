
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; prepare time range:
prepTimeRange <- function(tr,tzDOAS=""){
  
  if(grepl("^GMT[+]",tzDOAS)){
    warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("GMT+","Etc/GMT-",tzDOAS,fixed=TRUE))    
  } else if(grepl("^UTC[+]",tzDOAS)){
    warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("UTC+","Etc/GMT-",tzDOAS,fixed=TRUE))     
  } else if(grepl("^GMT[-]",tzDOAS)){
    warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("GMT-","Etc/GMT+",tzDOAS,fixed=TRUE))     
  } else if(grepl("^UTC[-]",tzDOAS)){
    warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("UTC-","Etc/GMT+",tzDOAS,fixed=TRUE))       
  }

  if(inherits(tr,"character")){
    evalperiod <- rep(tr,2)[1:2]
    timerange <- parse_date_time(evalperiod, c("d.m.Y","d-m-Y","Y-m-d"), tz=tzDOAS, quiet=TRUE)
    if(!is.na(timerange[2])) timerange[2] <- timerange[2] + 24*3600
    timerange2 <- parse_date_time(evalperiod, c("d.m.Y H:M:S","d.m.Y H:M","d-m-Y H:M:S","d-m-Y H:M","Y-m-d H:M:S","Y-m-d H:M"), tz=tzDOAS, quiet=TRUE)
    timerange[is.na(timerange)] <- timerange2[is.na(timerange)]
  } else if(inherits(tr,"POSIXlt")){
    timerange <- as.POSIXct(tr)
  } else if(!inherits(tr,"POSIXct")){
    stop("argument \"tr\" must be either of class character or POSIXt")
  } else {
    timerange <- tr
  }

  return(timerange)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; split time range:
splitTimeRange <- function(tr,n_days=10,length.out=NULL){
	if(!is.null(length.out)){
		stopifnot(is.integer(length.out)|length.out>1)
		n_days <- ceiling(as.numeric(diff(tr)/(length.out - 1),"days"))
	}
	n_days <- as.difftime(n_days,unit="days")
	tr2 <- trunc(tr[1],"days")
	if(tr[1]-tr2>as.difftime(12,unit="hours")){
		tr2 <- tr2 + as.difftime(1,unit="days")
	}
	timerange <- seq(tr2,tr[2],n_days)
	timerange[1] <- tr[1]
	if(timerange[length(timerange)] != tr[2]){
		tzone <- tz(timerange)
		last <- tr[2] - timerange[length(timerange)]
		if(last<n_days/2){
			timerange[length(timerange)] <- tr[2]
		} else {
			timerange <- with_tz(c(timerange,tr[2]),tzone)
		}
	}
  return(timerange)
}


### describe spectrometer's linearity (see Ocean Optics calibration sheets)
### ****************************************************************************** 
linearity.func <- function(x, cfs){
  y <- cfs[1] 
  for (i in 1:(length(cfs)-1)) {
    y <- y + cfs[i+1] * x^i  
  }
  return(y)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; read DOAS raw data:
# DOASinfo <- DOAS.info 
# dataDir <- rawdata.dir
# rawdataOnly=FALSE
# skip.check.daily=skip.check.daily
# force.write.daily=force.write.daily
# timerange=DOASinfo$timerange
# timezone=""
# ncores=ncores
write_daily <- function(files, obj, path) {
    # open connection
    con <- gzfile(paste0(path, '.bin'), 'wb')
    on.exit({
        close(con)
        unlink(paste0(path, '.bin'))
    })
    # write # files
    writeBin(length(files), con, endian = 'big')
    # write files
    writeBin(files, con, endian = 'big')
    # serialize object
    ser <- serialize(obj, NULL, TRUE)
    # write length of ser
    writeBin(length(ser), con, endian = 'big')
    # write list
    writeBin(ser, con, endian = 'big')
    # close connection
    close(con)
    on.exit()
}

read_daily <- function(path, files = FALSE) {
    # open connection
    con <- gzfile(path, 'rb')
    on.exit(close(con))
    # read # files
    N1 <- readBin(con, 'int', 1L, endian = 'big')
    # read files
    out <- readBin(con, 'character', N1, endian = 'big')
    if (!files) {
        # read length of serialized data
        N2 <- readBin(con, 'int', 1L, endian = 'big')
        # read serialized data and unserialize
        out <- unserialize(readBin(con, 'raw', N2, endian = 'big'))
    }
    # close connection
    close(con)
    on.exit()
    # return result
    out
}

process_dailyfiles <- function(folders, path_data, doas_info, RawData) {
    # loop over folders
    lapply(folders, function(i) {
        cat(paste0("..processing daily file ", i, ".bin\n"))
        # get text files
        rawdat_all <- list.files(file.path(path_data, i), pattern = ".txt", full.names = TRUE)
        # check invalid files
        rawdatSize <- file.size(rawdat_all)
        rawdat <- rawdat_all[rawdatSize > 4000]
        if(sum(rawdatSize <= 4000)){
            cat(paste0(sum(rawdatSize <= 4000), " invalid files (< 4000 bytes)!\n"))
            file.copy(rawdat_all[rawdatSize <= 4000], sub(".txt$", ".invalid", rawdat_all[rawdatSize <= 4000]), overwrite = TRUE)
        }
        # read raw data
        if(length(rawdat)){
            if(doas_info$DOASmodel == "S1"&&(as.numeric(i)-20160101)>=0&&(as.numeric(i)-20170101)<0){
                # S1 data between 2016 and 2017
                dat <- data.frame(
                    wavelength=c(seq.int(doas_info$rawdata.structure$"Header Lines"),doas_info$Spectrometer$wavelength[10:1021]),
                    lapply(rawdat, function(x) read.table(x, header=FALSE, nrows=doas_info$Spectrometer$"Pixel Number"+doas_info$rawdata.structure$"Header Lines"-9,sep="\n",as.is=TRUE,quote="",comment.char=""))
                )
                input <- data.frame(wavelength=doas_info$Spectrometer$wavelength[1:9],matrix(NA,nrow=9,ncol=ncol(dat)-1))
                names(input) <- names(dat)
                dat <- rbind(dat[1:8,],
                    input,
                    dat[9:1020,],make.row.names=FALSE)
            } else {
                # define first row
                firstRow <- c(seq.int(doas_info$rawdata.structure$"Header Lines"), doas_info$Spectrometer$wavelength)
                # read data
                fls <- lapply(rawdat, function(x) read.table(x, header=FALSE, nrows=doas_info$Spectrometer$"Pixel Number"+doas_info$rawdata.structure$"Header Lines",sep="\n",as.is=TRUE,quote="",comment.char=""))
                # temporary solution to missing Klima entries:
                if (nrow(fls[[1]]) == (length(firstRow) - 1)) {
                    fls <- lapply(fls, function(x) rbind(x[1:5,, drop = FALSE], 
                            "Board- / Ambient-T (degC) / ambient-RH (perc)= NA, NA, NA,", 
                            x[6:1051,, drop = FALSE]))
                } else if (any(check_fls <- sapply(fls,nrow) != length(firstRow))) {
                    stop("check files:\n", paste(rawdat[check_fls], collapse = "\n"))
                }
                # cbind
                dat <- data.frame(
                    wavelength=firstRow,
                    fls
                )
            }
            # prepare header
            header <- list(
                RevPos = sub("^.*[=]", "", dat[doas_info$rawdata.structure$"Revolver Position", -1]),
                ShuPos = sub("^.*[=]", "", dat[doas_info$rawdata.structure$"Shutter Position", -1]),
                ststr = sub("^.*_", "", dat[doas_info$rawdata.structure$"Acquisition start time", -1]),
                etstr = sub("^.*_", "", dat[doas_info$rawdata.structure$"Acquisition stop time", -1]),
                TECTemp = sub("^.*[=]", "", dat[doas_info$rawdata.structure$"TEC-Temp (degC)", -1]),
                Klima = sub("^.*[=]", "", dat[doas_info$rawdata.structure$"Board- / Ambient-T (degC) / ambient-RH (perc)", -1]),
                Expos = as.numeric(sub("^.*[=]", "", dat[doas_info$rawdata.structure$"Exposure Time (ms)", -1])),
                AccNum = as.numeric(sub("^.*[=]", "", dat[doas_info$rawdata.structure$"Number of accumulations", -1])),
                st = parse_date_time(sub("^.*_", "", dat[doas_info$rawdata.structure$"Acquisition start time", -1]), doas_info$rawdata.structure$"Time Format", tz = tz(doas_info$timerange)),
                et = parse_date_time(sub("^.*_", "", dat[doas_info$rawdata.structure$"Acquisition stop time", -1]), doas_info$rawdata.structure$"Time Format", tz = tz(doas_info$timerange))
                )
            # prepare raw data as list
            out <- list(intensity = apply(dat[-seq.int(doas_info$rawdata.structure$'Header Lines'), ], 2, as.numeric), header = header)
            # write daily file (files, data, path)
            write_daily(basename(rawdat), out, file.path(path_data, i))
            # return list
            out
        } else {
            # no data
            cat("Folder", i, "is empty...\n")
            # return null
            NULL
        }
    })
}

check_dailyfiles <- function(rawdat_folder, skip.check.daily, force.write.daily, 
    lf_raw_dir, path_data, doas_info, rd_list) {
    # target daily files
    target_daily <- paste0(rawdat_folder, ".bin")
    # check existing daily files
    if (force.write.daily) {
        # "no daily files exist"
        daily_exist <- rep(FALSE, length(target_daily))
    } else {
        # get existing daily files
        daily_exist <- target_daily %in% lf_raw_dir
        # check daily files
        if(!skip.check.daily){
            cat("checking daily files...\n")
            for (i in which(daily_exist)) {
                # read 'old' files
                files_daily <- try(read_daily(file.path(path_data, target_daily[i]), files = TRUE), silent = TRUE)
                # read present files
                files_now <- list.files(file.path(path_data, rawdat_folder[i]), pattern = '[.]txt')
                # check for changes
                if (identical(files_daily, files_now)){
                    cat(paste0("..daily file ", target_daily[i], " exists and is ok...\n"))
                    # read data from daily file
                    rd_list[[rawdat_folder[i]]] <- read_daily(file.path(path_data, target_daily[i]))
                } else {
                    # delete daily file (and re-create it later)
                    unlink(file.path(path_data, target_daily[i]))
                    daily_exist[i] <- FALSE
                }
            }   
        } else {
            # read existing daily files
            for (i in which(daily_exist)) {
                cat("reading daily file", target_daily[i], "...\n")
                rd_list[[rawdat_folder[i]]] <- read_daily(file.path(path_data, target_daily[i]))
            }     
        }
    }
    # process daily files
    rd_list[!daily_exist] <- process_dailyfiles(rawdat_folder[!daily_exist], path_data, doas_info, rd_list)
    # return
    rd_list[!is.null(rd_list)]
}

readDOASdata <- function(DOASinfo, dataDir, rawdataOnly = FALSE, skip.check.daily = FALSE, 
    force.write.daily = FALSE, timerange = DOASinfo$timerange, timezone = "", ncores = 1) {

    if(parl <- (ncores > 1 & !rawdataOnly)) {
        require(snowfall)
        if(!sfIsRunning()){     
            on.exit(sfStop())
            sfInit(TRUE, ncores)
        }
        cl <- sfGetCluster()
        sfLibrary(lubridate)
    }

    if(inherits(DOASinfo, "character")){
        timerange <- prepTimeRange(timerange, timezone)
        DOASinfo <- getDOASinfo(DOASinfo, timerange)
    }

    # correct folder for last day 
    tr.to <- trunc(timerange[2], "days")
    if(as.numeric(tr.to - timerange[2]) == 0) tr.to <- tr.to - as.difftime(1,units="days")

    # target folders
    rawdatfolder <- format(seq(from=trunc(timerange[1],"days"), to=tr.to, 
            by=as.difftime(1, units="days")), "%Y%m%d") 

    # get existing files/folders
    lf.raw.dir <- list.files(dataDir, pattern = '[0-9]{8}([.]bin)?$')
    if(!length(lf.raw.dir)){
        if(!dir.exists(dirname(dataDir))){
            stop("Folder \"",dirname(dataDir),"\" does not exist!")
        } else if(!dir.exists(dataDir)){
            stop("Folder \"",dataDir,"\" does not exist!")
        } else {
            stop("Folder \"",dataDir,"\" does not contain any doas data!")
        }
    } else if(!any(rawdatfolder %in% lf.raw.dir)){
        stop("No data available for specified timerange!")
    }

    # prepare raw data output
    RawData <- vector(length(rawdatfolder), mode = "list")
    names(RawData) <- rawdatfolder

    # if rawdataOnly else daily files
    if (rawdataOnly) {

        # newly introduced raw_data
        raw_data <- RawData

        # loop over folders
        for (i in rawdatfolder) {
            cat(paste0("..reading raw data from folder ",i,"\n"))
            # list files
            rawdat_all <- list.files(file.path(dataDir, i), pattern = ".txt", full.names = TRUE)
            # read files
            if(length(rawdat_all)){
                if(DOASinfo$DOASmodel=="S1" && DOASinfo$timerange[1]>strptime("20160101","%Y%m%d") && DOASinfo$timerange[1]<strptime("20170101","%Y%m%d")){
                    rawdat.et <- parse_date_time(gsub("^.*/(.*).txt","\\1",rawdat_all),c("%Y%m%d %H%M%OS", "%y%m%d%H%M%S"), tz=tz(DOASinfo$timerange))
                } else {
                    rawdat.et <- parse_date_time(gsub("^.*_(.*).txt","\\1",rawdat_all),c("%Y%m%d %H%M%OS", "%y%m%d%H%M%S"), tz=tz(DOASinfo$timerange))
                }
                rawdat.st <- rawdat.et - median(diff(rawdat.et))
                takeme <- (rawdat.et >= timerange[1]) & (rawdat.st <= timerange[2])
                rawdat_all <- rawdat_all[takeme]
                rawdatSize <- file.size(rawdat_all)
                rawdat <- rawdat_all[rawdatSize > 4000]
                if(sum(rawdatSize <= 4000)){
                    cat(paste0(sum(rawdatSize <= 4000)," invalid files (< 4000 bytes)!\n"))
                    file.copy(rawdat_all[rawdatSize <= 4000],gsub(".txt$",".invalid",rawdat_all[rawdatSize <= 4000]))
                }
                if(length(rawdat)){
                    if(DOASinfo$DOASmodel == "S1"&&(as.numeric(i)-20160101)>=0&&(as.numeric(i)-20170101)<0){
                        dat <- data.frame(
                            wavelength=c(seq.int(DOASinfo$rawdata.structure$"Header Lines"),DOASinfo$Spectrometer$wavelength[10:1021]),
                            lapply(rawdat, function(x) read.table(x, header=FALSE, nrows=DOASinfo$Spectrometer$"Pixel Number"+DOASinfo$rawdata.structure$"Header Lines"-9,sep="\n",as.is=TRUE,quote="",comment.char=""))
                        )
                        input <- data.frame(wavelength=DOASinfo$Spectrometer$wavelength[1:9],matrix(NA,nrow=9,ncol=ncol(dat)-1))
                        names(input) <- names(dat)
                        dat <- rbind(dat[1:8,],
                            input,
                            dat[9:1020,],make.row.names=FALSE)
                    } else {
                        # read files
                        firstRow <- c(seq.int(DOASinfo$rawdata.structure$"Header Lines"),DOASinfo$Spectrometer$wavelength)
                        fls <- lapply(rawdat, function(x) read.table(x, header=FALSE, nrows=DOASinfo$Spectrometer$"Pixel Number"+DOASinfo$rawdata.structure$"Header Lines",sep="\n",as.is=TRUE,quote="",comment.char=""))
                        # temporary solution to missing Klima entries:
                        if(nrow(fls[[1]]) == (length(firstRow) - 1)){
                            fls <- lapply(fls, function(x) rbind(x[1:5,, drop = FALSE], 
                                    "Board- / Ambient-T (degC) / ambient-RH (perc)= NA, NA, NA,", 
                                    x[6:1051,, drop = FALSE]))
                        } else if(any(check_fls <- sapply(fls,nrow) != length(firstRow))){
                            stop("check files:\n",paste(rawdat[check_fls],collapse="\n"))
                        }
                        dat <- data.frame(
                            wavelength=c(seq.int(DOASinfo$rawdata.structure$"Header Lines"),DOASinfo$Spectrometer$wavelength),
                            fls)
                    }
                    # prepare header
                    header <- list(
                        RevPos = sub("^.*[=]", "", dat[DOASinfo$rawdata.structure$"Revolver Position", -1]),
                        ShuPos = sub("^.*[=]", "", dat[DOASinfo$rawdata.structure$"Shutter Position", -1]),
                        ststr = sub("^.*_", "", dat[DOASinfo$rawdata.structure$"Acquisition start time", -1]),
                        etstr = sub("^.*_", "", dat[DOASinfo$rawdata.structure$"Acquisition stop time", -1]),
                        TECTemp = sub("^.*[=]", "", dat[DOASinfo$rawdata.structure$"TEC-Temp (degC)", -1]),
                        Klima = sub("^.*[=]", "", dat[DOASinfo$rawdata.structure$"Board- / Ambient-T (degC) / ambient-RH (perc)", -1]),
                        Expos = as.numeric(sub("^.*[=]", "", dat[DOASinfo$rawdata.structure$"Exposure Time (ms)", -1])),
                        AccNum = as.numeric(sub("^.*[=]", "", dat[DOASinfo$rawdata.structure$"Number of accumulations", -1])),
                        st = parse_date_time(sub("^.*_", "", dat[DOASinfo$rawdata.structure$"Acquisition start time", -1]), DOASinfo$rawdata.structure$"Time Format", tz = tz(DOASinfo$timerange)),
                        et = parse_date_time(sub("^.*_", "", dat[DOASinfo$rawdata.structure$"Acquisition stop time", -1]), DOASinfo$rawdata.structure$"Time Format", tz = tz(DOASinfo$timerange))
                    )

                    # write raw data as list
                    raw_data[[i]] <- list(intensity = apply(dat[-seq.int(DOASinfo$rawdata.structure$'Header Lines'), ], 2, as.numeric), header = header)
                } else {
                    raw_data <- raw_data[-match(i, names(raw_data))]
                } 
            } else {
                raw_data <- raw_data[-match(i, names(raw_data))]
                cat("No data in folder", i, "\n")
            }
        }
    } else if (parl) {

        # export functions
        clusterExport(cl, list('process_dailyfiles', 'write_daily', 'read_daily'))

        # call check daily in parallel
        cat("checking daily files in parallel...\n")
        raw_data <- clusterApplyLB(cl, seq_along(RawData), function(i, check_fu, rdf, rd, ...) {
            unlist(
                check_fu(rdf[[i]], rd_list = rd[[i]], ...),
                recursive = FALSE)
                            }, 
            check_fu = check_dailyfiles, rdf = rawdatfolder, rd = RawData, skip.check.daily, force.write.daily,
            lf.raw.dir, dataDir, DOASinfo)

    } else {

        raw_data <- check_dailyfiles(rawdatfolder, skip.check.daily, force.write.daily, 
            lf.raw.dir, dataDir, DOASinfo, RawData)

    }

    # remove empty folders
    raw_data <- raw_data[sapply(raw_data, function(x) !is.null(x))]

    if (!length(raw_data)) {
        stop("No data available for specified time range!")
    } else if (length(unique(sapply(raw_data, nrow))) != 1) {
        stop("Daily Files have differing number of rows! Setting 'force.write.daily = TRUE' might solve the problem.")
    }

    # read header
    header <- do.call(rbind, lapply(raw_data, function(x) data.frame(x[[2]])))

    # get eval period subset
    takeme <- (header$et >= timerange[1]) & (header$st <= timerange[2])

    # bind raw data
    RawData <- data.frame(lapply(raw_data, function(x) x[[1]][, -1]))[, takeme, drop = FALSE]

    # get timerange:
    if (!length(RawData)) stop("No data available for specified time range!\n")

    # subset header
    header <- header[takeme, ]

    # stupid S1
    if(DOASinfo$DOASmodel=="S1" && DOASinfo$timerange[1] < as.POSIXct("2016-01-01", tz=tz(DOASinfo$timerange[1]))){
        RawData <- RawData[1:1021,,drop=FALSE]
        DOASinfo$Spectrometer$"Pixel Number" <- 1021
        # wavelength
        DOASinfo$Spectrometer$pixel <- seq.int(1021)
        DOASinfo$Spectrometer$wavelength <- DOASinfo$Spectrometer$"Calibration Coefficients"[1] + DOASinfo$Spectrometer$"Calibration Coefficients"[2]*DOASinfo$Spectrometer$pixel + DOASinfo$Spectrometer$"Calibration Coefficients"[3]*DOASinfo$Spectrometer$pixel^2 + DOASinfo$Spectrometer$"Calibration Coefficients"[4]*DOASinfo$Spectrometer$pixel^3
        warning("S1: old data structure. Cut high end data!")
    }

    # write header
    Header <- data.frame(
        DOASmodel = DOASinfo$DOASmodel,
        Spectrometer = DOASinfo$Spectrometer$"Spectrometer Name",
        header[, 1:8],
        AccTime = round(as.numeric(difftime(header$et, header$st, units = "secs")), 3),
        st = header$st,
        et = header$et,
        stringsAsFactors = FALSE
    )

    # devide total counts by # of accumulations
    RawData <- sweep(RawData, 2, Header[, "AccNum"], "/")

    return(list(RawData=RawData,Header=Header,DOASinfo=DOASinfo))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; average Spectra:
avgSpec <- function(rawdat,type=c("raw","cal","ref","dark"),tracer=c("ambient","NH3","NO","SO2"),saveToFile=FALSE
  ,pathLength=NA,NH3ambient_ug=NA,NOambient_ug=NA,SO2ambient_ug=NA,cuvetteLength=ifelse(tracer[1]=="ambient",NA,0.075),cuvetteConc_mg=NA,Dirname=getwd(),verbose=TRUE){

  # Calc Sum/n:
  # Specs <- apply(rawdat[[1]][-(1:13),],1,as.numeric)/rawdat[[2]][,"AccNum"]
  Specs <- rawdat[[1]]
  # average specs:
  # SpecAvg <- colMeans(Specs,na.rm=TRUE)
  SpecAvg <- rowMeans(Specs,na.rm=TRUE)

  type <- type[1]
  tracer <- tracer[1]

  cuvetteGas <- switch(type,
    cal = {
      if(tracer=="ambient")stop("ambient & cal specified!")
      if(is.na(cuvetteConc_mg))stop("Calibration concentration must be specified!")
      tracer
      },
    ref = ifelse(tracer=="ambient","none","N2")
      ,
    NA
    )
  
  Info <- list(
    type=type,
    tracer=tracer,
    pathLength=pathLength,
    NH3ambient_ug=NH3ambient_ug,
    NOambient_ug=NOambient_ug,
    SO2ambient_ug=SO2ambient_ug,
    cuvetteLength=cuvetteLength,
    cuvetteConc_mg=cuvetteConc_mg,
    cuvetteGas=cuvetteGas,
    RevolverPosition=paste(unique(rawdat[[2]][,"RevPos"]),sep=","),
    ShutterPosition=paste(unique(rawdat[[2]][,"ShuPos"]),sep=","),
    st=rawdat[[2]][,"st"],
    et=rawdat[[2]][,"et"],
    NumberOfSpecs=length(rawdat[[2]][,"et"]),
    AccNum=list(range=range(rawdat[[2]][,"AccNum"]),mean=mean(rawdat[[2]][,"AccNum"])),
    Expos=list(range=range(rawdat[[2]][,"Expos"]),mean=mean(rawdat[[2]][,"Expos"])),
    TEC=list(range=range(as.numeric(rawdat[[2]][,"TECTemp"])),mean=mean(as.numeric(rawdat[[2]][,"TECTemp"]))),
    Ambient=list(range=range(as.numeric(gsub(",.*$","",rawdat[[2]][,"Klima"]))),mean=mean(as.numeric(gsub(",.*$","",rawdat[[2]][,"Klima"])))),
    Imax=list(range=range(apply(Specs,2,max)),mean=max(SpecAvg))
    )   

  insertme <- paste(tracer,switch(type,raw="measured",cal="calibration",ref="reference",dark="dark"),"spectra",sep=" ")
  txt <- c(
    paste0(insertme," file")
    ,paste0("calculation performed on: ",format(Now <- Sys.time()))
    ,paste0("miniDOAS: ",rawdat$DOASinfo$DOASmodel)
    ,paste0("spectrometer: ",rawdat$DOASinfo$Spectrometer$"Spectrometer Name")
    ,paste0("revolver position: ",paste(unique(rawdat[[2]][,"RevPos"]),sep=","))
    ,paste0("shutter position: ",paste(unique(rawdat[[2]][,"ShuPos"]),sep=","))
    ,"~~~ cuvette measurement ~~~"
    ,paste0("cuvette installed: ",ifelse((cuvetteGas=="none"|is.na(cuvetteGas)),"no","yes"))
    ,paste0("cuvette gas: ",cuvetteGas)
    ,paste0("cuvette concentration (mg/m3): ",cuvetteConc_mg)
    ,paste0("cuvette length: ",cuvetteLength)
    ,"~~~ path measurement ~~~"
    ,paste0("path length: ",pathLength)
    ,paste0("ambient NH3 concentration (ug/m3): ",NH3ambient_ug)
    ,paste0("ambient SO2 concentration (ug/m3): ",SO2ambient_ug)
    ,paste0("ambient NO concentration (ug/m3): ",NOambient_ug)
    ,"~~~ "
    ,paste0("averaging period: ",format(rawdat[[2]][1,"st"])," to ",format(rev(rawdat[[2]][,"et"])[1]))
    ,paste0("averaged ",length(rawdat[[2]][,"et"])," spectra")
    ,paste0("number of accumulations: ",paste(range(rawdat[[2]][,"AccNum"]),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(rawdat[[2]][,"AccNum"])),")")
    ,paste0("exposure time: ",paste(range(rawdat[[2]][,"Expos"]),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(rawdat[[2]][,"Expos"])),")")
    ,paste0("TEC temperature (deg C): ",paste(range(as.numeric(rawdat[[2]][,"TECTemp"])),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(as.numeric(rawdat[[2]][,"TECTemp"]))),")")
    ,paste0("ambient temperature (deg C): ",paste(range(as.numeric(gsub(",.*$","",rawdat[[2]][,"Klima"]))),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(as.numeric(gsub(",.*$","",rawdat[[2]][,"Klima"])))),")")
    ,paste0("Imax: ",sprintf("%1.1f to %1.1f",min(apply(Specs,2,max)),max(apply(Specs,2,max)))," (avg: ",sprintf("%1.1f",max(SpecAvg)),")")
    ,"-----")

  if(saveToFile){
    filename <- paste0("miniDOAS_",rawdat$DOASinfo$DOASmodel,"_",paste(tracer,type,"spec",sep="_"),"_",paste(format(c(rawdat[[2]][1,"st"],rev(rawdat[[2]][,"et"])[1]),"%y%m%d%H%M%S"),collapse="-"),".txt")
    write.table(txt,paste(Dirname,filename,sep="/"),quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
    write.table(data.frame(rawdat$DOASinfo$Spectrometer$wavelength,SpecAvg),paste(Dirname,filename,sep="/"),quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
    cat("Averaged spectra saved to",paste(Dirname,filename,sep="/"),"\n")
  }

  txt[1] <- paste0("Averaging ",insertme,":\n") 
  txt <- paste(c(txt,"\n"),collapse="\n")

  if(verbose){
    cat(txt)
  }

  return(list(SpecAvg=SpecAvg,Info=Info,Specs=Specs,DOASinfo=rawdat$DOASinfo,txt=txt))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; get dark/filter/fit windows and filter strength:
getWindows <- function(DOASinfo,filter.type="BmHarris",timerange=Sys.time(),straylight.window=NULL,filter.window=NULL,fit.window=NULL,filter.strength=NULL,tau.shift=NULL,double=TRUE
  ,special.control=list(b=3.5,Scale=function(r) median(abs(r))/0.6745,delta=0,filter.strength_multiplication=1.25,filter.strength_loess=0.2,fam="gaussian",maxit=c(10,0))){
  
  if(inherits(DOASinfo,"character"))DOASinfo <- getDOASinfo(DOASinfo,timerange)

  check.fit.win <- FALSE
  if (DOASinfo$Spectrometer$Serial == 'AvantesS1') {
    if(any(filter.type == c("loess","special"))) {
      if(is.null(filter.window))filter.window <- c(202.2,226) 
      if(is.null(filter.strength))filter.strength <- switch(filter.type, loess=0.22, special=0.17)
      if(is.null(fit.window))fit.window <- c(204.6,223.6)
    } else if(any(filter.type == c("Rect","Hann","Hamming","Blackman","BmNuttall"
      ,"Sin","Gauss","Kaiser","DolphChebyshev","BmHarris","Tukey","Poisson","Exp","ExpHamming"))){
      check.fit.win <- TRUE
      if(is.null(filter.window))filter.window <- c(203, 226)                                                                                                                            ### wavelength window of evaluation: c(min,max) [nm] (from this, in case of moving average lowpass filter applications, the edge of the filter is thrown away in addition)
      if(is.null(filter.strength)){
        filter.strength <- 41
        } else if(!as.logical(filter.strength%%2)){
          filter.strength <- round(filter.strength)
          if(!as.logical(filter.strength%%2)){
            filter.strength <- filter.strength + 1
          }
          warning("filter.strength must be an uneven number and has been rounded to the nearest uneven number.")
        }  
      if(is.null(fit.window))fit.window <- c(204.1,224.85)
    } else if(filter.type == "FlatTop"){
      check.fit.win <- TRUE
      if(is.null(filter.window))filter.window <- c(203, 226)                                                                                                                            ### wavelength window of evaluation: c(min,max) [nm] (from this, in case of moving average lowpass filter applications, the edge of the filter is thrown away in addition)
      if(is.null(filter.strength)){
        filter.strength <- 71
        } else if(!as.logical(filter.strength%%2)){
          filter.strength <- round(filter.strength)
          if(!as.logical(filter.strength%%2)){
            filter.strength <- filter.strength + 1
          }
          warning("filter.strength must be an uneven number and has been rounded to the nearest uneven number.")
        } 
      if(is.null(fit.window))fit.window <- c(204.1,224.85)      
    } else {
      stop("filter.type = \"",filter.type,"\" is not a valid filter type")
    }

    if(is.null(straylight.window))straylight.window <- c(250, 260)                                                                                                        ### wavelength (nm) window at which the dark-drift of spectrometer is monitored (at ccd pixels > ca. 230 nm that are covered)                                                                                                                     ### performs the multiple linear fit with weights according to the absolut contribution of all absorption bands over the wavelength window

  } else {
    if (any(filter.type == c("loess","special"))) {
      # if(is.null(filter.window))filter.window <- c(203,229.7) 
      if(is.null(filter.window))filter.window <- c(202.2,230.6) 
      if(is.null(filter.strength))filter.strength <- switch(filter.type, loess=0.22, special=0.17)
      if(is.null(fit.window))fit.window <- c(204.3,228.6)
    } else if(any(filter.type == c("Rect","Hann","Hamming","Blackman","BmNuttall"
      ,"Sin","Gauss","Kaiser","DolphChebyshev","BmHarris","Tukey","Poisson","Exp","ExpHamming"))){
      check.fit.win <- TRUE
      if(is.null(filter.window))filter.window <- c(202.2,230.6)                                                                                                                            ### wavelength window of evaluation: c(min,max) [nm] (from this, in case of moving average lowpass filter applications, the edge of the filter is thrown away in addition)
      if(is.null(filter.strength)){
        filter.strength <- 25
        } else if(!as.logical(filter.strength%%2)){
          filter.strength <- round(filter.strength)
          if(!as.logical(filter.strength%%2)){
            filter.strength <- filter.strength + 1
          }
          warning("filter.strength must be an uneven number and has been rounded to the nearest uneven number.")
        }  
      if(is.null(fit.window))fit.window <- c(204.5,228.4)
    } else if(filter.type == "FlatTop"){
      check.fit.win <- TRUE
      if(is.null(filter.window))filter.window <- c(202.2,230.6)                                                                                                                            ### wavelength window of evaluation: c(min,max) [nm] (from this, in case of moving average lowpass filter applications, the edge of the filter is thrown away in addition)
      if(is.null(filter.strength)){
        filter.strength <- 51
        } else if(!as.logical(filter.strength%%2)){
          filter.strength <- round(filter.strength)
          if(!as.logical(filter.strength%%2)){
            filter.strength <- filter.strength + 1
          }
          warning("filter.strength must be an uneven number and has been rounded to the nearest uneven number.")
        } 
      if(is.null(fit.window))fit.window <- c(204.5,228.4)     
    } else {
      stop("filter.type = \"",filter.type,"\" is not a valid filter type")
    }
    if(is.null(straylight.window))straylight.window <- c(150,160) # c(141.5059, 145.3158) alte Version   ### with the current ccd wavelength calibration (pixel 50:90) these correspond to a range at which the straylight induced offset might be monitored 
  }
  if (is.null(tau.shift)||is.na(tau.shift)) tau.shift <- 0                                                      ### if tau consideration is switched off => tau=0

  ### wavelength/pixel indices
  ### ******************************************************************************
  pixel1 <- which(DOASinfo$Spectrometer$wavelength >= filter.window[1] & DOASinfo$Spectrometer$wavelength <= filter.window[2])

  pixel2 <- max(
    which(DOASinfo$Spectrometer$wavelength >= fit.window[1])[1]
    ,pixel1[1]-tau.shift
    ):min(
    rev(which(DOASinfo$Spectrometer$wavelength <= fit.window[2]))[1]
    ,rev(pixel1)[1]-tau.shift
    )
  pixel3 <- which(DOASinfo$Spectrometer$wavelength >= straylight.window[1] & DOASinfo$Spectrometer$wavelength <= straylight.window[2])
  if(check.fit.win){
    fs <- if(double) filter.strength else (filter.strength - 1)/2
    if(abs(pixel1[1]-pixel2[1])<fs){
      pixel2Lo <- pixel1[1]+fs
      fit.window[1] <- DOASinfo$Spectrometer$wavelength[pixel2Lo]
      warning(sprintf("Too broad filtering: lower limit of argument \"fit.window\" is increased from %1.4f to %1.4f nm",DOASinfo$Spectrometer$wavelength[pixel2[1]],fit.window[1]))
    } else {
      pixel2Lo <- pixel2[1]
    }
    if(abs(rev(pixel1)[1]-rev(pixel2)[1])<fs){
      pixel2Hi <- rev(pixel1)[1]-fs
      fit.window[2] <- DOASinfo$Spectrometer$wavelength[pixel2Hi]
      warning(sprintf("Too broad filtering: upper limit of argument \"fit.window\" is decreased from %1.4f to %1.4f nm",DOASinfo$Spectrometer$wavelength[rev(pixel2)[1]],fit.window[2]))
    } else {
      pixel2Hi <- rev(pixel2)[1]
    }
    pixel2 <- seq(pixel2Lo,pixel2Hi)
  }
  pixel4 <- match(pixel2,pixel1)
  if(anyNA(pixel4)){
    pixel4 <- pixel4[!is.na(pixel4)]
    pixel2 <- pixel1[pixel4]
    fit.window <- c(DOASinfo$Spectrometer$wavelength[pixel2[1]],DOASinfo$Spectrometer$wavelength[rev(pixel2)[1]])
    warning(sprintf("\"fit.window\" does not lie within the boundaries of \"filter.window\"!\n\targument \"fit.window\" is adjusted to c(%1.4f, %1.4f)",fit.window[1],fit.window[2]))
  }

  special <- list(b=3.5,Scale=function(r) median(abs(r))/0.6745,delta=0,filter.strength_multiplication=1.25,filter.strength_loess=0.2,fam="gaussian",maxit=c(10,0))
  special[names(special.control)] <- special.control

  return(list(
    filter.window=filter.window,
    fit.window=fit.window,
    straylight.window=straylight.window,
    filter.type=filter.type,
    filter.strength=filter.strength,
    double=double,
    tau.shift=tau.shift,
    pixel1=pixel1,
    pixel2=pixel2,
    pixel3=pixel3,
    pixel4=pixel4,
    special.control=special
    ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; get calibration and reference spectra:
getSpec <- function(spec,DOASmodel=NULL,lite=FALSE){
  
  SpecName <- as.character(substitute(spec))

  tracer <- switch(SpecName,
    "ref.spec"="ambient",
    "ref.dark.spec"="ambient",
    "dark.spec"="ambient",
    "NH3.cal.spec"="NH3",
    "SO2.cal.spec"="SO2",
    "NO.cal.spec"="NO",
    "N2.NH3.cal.spec"="NH3",
    "N2.SO2.cal.spec"="SO2",
    "N2.NO.cal.spec"="NO",
    "N2.dark.cal.spec"="N2"
  )

  type <- switch(SpecName,
    "ref.spec"="ref",
    "ref.dark.spec"="dark",
    "dark.spec"="dark",
    "NH3.cal.spec"="cal",
    "SO2.cal.spec"="cal",
    "NO.cal.spec"="cal",
    "N2.NH3.cal.spec"="ref",
    "N2.SO2.cal.spec"="ref",
    "N2.NO.cal.spec"="ref",
    "N2.dark.cal.spec"="dark"
  ) 

  if(is.list(spec)){
    if(all(names(spec)%in%c("dat.spec","ambient","cuvette","timerange"))){
      dat.spec <- spec$dat.spec
      ambient <- spec$ambient
      cuvette <- spec$cuvette
      timerange <- spec$timerange
      lite <- TRUE    
    } else {
      if(!is.null(spec$type)){
        stopifnot(spec$type==type)
      } else {
        spec$type <- type
      }

      if(type=="cal"){
        stopifnot(
          !is.null(spec$cuvetteConc_mg),
          spec$cuvetteConc_mg>0,
          is.null(spec$cuvetteLength)||spec$cuvetteLength>0
          )
      }
      
      if(!is.null(spec$tracer)){
        if(SpecName == "N2.dark.cal.spec") spec$tracer <- "N2"
        stopifnot(spec$tracer==tracer)
      } else {
        spec$tracer <- tracer
      }

      if(is.null(spec$tz)) spec$tz <- "Etc/GMT-1"
      spec$rawdat <- try(readDOASdata(DOASmodel,timerange=spec$timerange,dataDir=spec$dir,rawdataOnly=TRUE,timezone=spec$tz))
      
      if(inherits(spec$rawdat,"try-error")){
        stop(paste0(SpecName,": ",unlist(strsplit(spec$rawdat," : "))[2]))
      }

      AvgSpec <- do.call(avgSpec,spec[names(spec) %in% names(formals(avgSpec))])

      if(lite){
        dat.spec <- AvgSpec$SpecAvg
        ambient <- AvgSpec$Info[3:6]
        cuvette <- AvgSpec$Info[7:8]
        timerange <- c(min(spec$rawdat$Header[,"st"]),max(spec$rawdat$Header[,"et"])) 
      }
    }
  } else {

    SpecName_old <- switch(as.character(SpecName),
        "ref.spec"="_refLampSpec_",
        "ref.dark.spec"="_refdarkLampSpec_",
        "dark.spec"="_darkSpec_",
        "NH3.cal.spec"="_NH3calSpec_",
        "SO2.cal.spec"="_SO2calSpec_",
        "NO.cal.spec"="_NOcalSpec_",
        "N2.NH3.cal.spec"="_N2calSpec_",
        "N2.SO2.cal.spec"="_N2calSpec_",
        "N2.NO.cal.spec"="_N2calSpec_",
        "N2.dark.cal.spec"="_darkN2calSpec_"
      )
    SpecName_new <- paste(tracer,type,"spec",sep="_")

    dir.spec <- list.files(path=spec, pattern="^miniDOAS_", full.names=TRUE)

    if(any(ind <- grepl(SpecName_old,dir.spec,fixed=TRUE))){
      info.spec <- read.table(dir.spec[which(ind)], sep="\t", stringsAsFactors=FALSE, fill=TRUE, colClasses="character",nrows=10)[,1]
      dat.spec <- read.table(dir.spec[which(ind)], sep="\t", stringsAsFactors=FALSE, fill=TRUE, skip=10,colClasses="numeric")[,1]
      ambient <- list(pathLength=NA,NH3ambient_ug=NA,SO2ambient_ug=NA,NOambient_ug=NA)
      if(any(SpecName==c("NH3.cal.spec","SO2.cal.spec","NO.cal.spec"))){
        konz.spec <- as.numeric(unlist(strsplit(gsub("^.*): ","",readLines(dir.spec[which(ind)],n=3)[3]),",")))
        cuvette <- list(cuvetteLength=0.075,cuvetteConc_mg=konz.spec[3]*konz.spec[1]/(8.3144598*konz.spec[2])*switch(tracer,NH3=17,SO2=64,NO=30)/10,cuvetteGas=tracer)
      } else {
        cuvette <- list(cuvetteLength=ifelse(tracer=="ambient",NA,0.075),cuvetteConc_mg=NA,cuvetteGas=ifelse(type=="ref",ifelse(tracer=="ambient","none","N2"),NA))
      }
    } else {
      ind <- grepl(SpecName_new,dir.spec,fixed=TRUE)
      if (all(!ind)) {
        if (type %in% 'dark') {
            ind[grep('dark', dir.spec)[1]] <- TRUE
        }
        if (all(!ind)) {
            stop('Calibration file for ', tracer, '/', type, ' is missing')
        }
      }
      info.spec <- read.table(dir.spec[which(ind)], sep="\t", stringsAsFactors=FALSE, fill=TRUE, colClasses="character",nrows=25)[,1]
      dat.spec <- read.table(dir.spec[which(ind)], sep=",", stringsAsFactors=FALSE, fill=TRUE, skip=25,colClasses="numeric")[,2]
      suppressWarnings(ambient <- list(
        pathLength=as.numeric(gsub("path length: ","",info.spec[13],fixed=TRUE))
        ,NH3ambient_ug=as.numeric(gsub("ambient NH3 concentration (ug/m3): ","",info.spec[14],fixed=TRUE))
        ,SO2ambient_ug=as.numeric(gsub("ambient SO2 concentration (ug/m3): ","",info.spec[15],fixed=TRUE))
        ,NOambient_ug=as.numeric(gsub("ambient NO concentration (ug/m3): ","",info.spec[16],fixed=TRUE))
        ))
       suppressWarnings(cuvette <- list(
        cuvetteLength=as.numeric(gsub("cuvette length: ","",info.spec[11],fixed=TRUE))
        ,cuvetteConc_mg=as.numeric(gsub("cuvette concentration (mg/m3): ","",info.spec[10],fixed=TRUE))
        ,cuvetteGas=gsub("cuvette gas: ","",info.spec[9],fixed=TRUE)
        ))  
    }
    lite <- TRUE
    timerange <- NA
    cat("getSpec: timerange weitergeben!\n")
  }

  if(DOASmodel=="S1"&&lite&&(is.na(timerange[1]) || timerange[1] < parse_date_time("20170101","Ymd",tz=tz(timerange[1])))){
      dat.spec <- dat.spec[seq.int(1021)]
  }

  if(lite){
    return(list(dat.spec=dat.spec,ambient=ambient,cuvette=cuvette,timerange=timerange))  
  } else {
    return(AvgSpec)
  }
}

getSpecSet <- function(
  spec.dir="",  
  ref.spec=NULL,
  ref.dark.spec=NULL,
  dark.spec=NULL,
  NH3.cal.spec=NULL,
  SO2.cal.spec=NULL,
  NO.cal.spec=NULL,
  N2.NH3.cal.spec=N2.cal.spec,
  N2.SO2.cal.spec=N2.cal.spec,
  N2.NO.cal.spec=N2.cal.spec,
  N2.dark.cal.spec=NULL,
  N2.cal.spec=NULL,
  DOAS.model=NULL){

  # if(is.null(ref.spec))ref.spec <- spec.dir
  if(is.null(ref.dark.spec))ref.dark.spec <- spec.dir
  if(is.null(dark.spec))dark.spec <- spec.dir
  if(is.null(NH3.cal.spec))NH3.cal.spec <- spec.dir
  if(is.null(SO2.cal.spec))SO2.cal.spec <- spec.dir
  if(is.null(NO.cal.spec))NO.cal.spec <- spec.dir
  if(is.null(N2.dark.cal.spec))N2.dark.cal.spec <- spec.dir
  if(is.null(N2.NH3.cal.spec))N2.NH3.cal.spec <- spec.dir
  if(is.null(N2.SO2.cal.spec))N2.SO2.cal.spec <- spec.dir
  if(is.null(N2.NO.cal.spec))N2.NO.cal.spec <- spec.dir

  ### read (average) reference, calibration & noise spectra
  ### ******************************************************************************
  return(
    list(
      DOAS.model = DOAS.model,
      ### lamp reference spectrum
      dat.ref = if (!is.null(ref.spec)) {
          getSpec(ref.spec,DOASmodel=DOAS.model,lite=TRUE)
      },
      ### lamp reference dark spectrum
      dat.ref.dark = getSpec(ref.dark.spec,DOASmodel=DOAS.model,lite=TRUE),
      ### actual dark spectrum
      dat.dark = getSpec(dark.spec,DOASmodel=DOAS.model,lite=TRUE),
      ### NH3 calibration spectrum
      dat.NH3 = getSpec(NH3.cal.spec,DOASmodel=DOAS.model,lite=TRUE),
      dat.N2.NH3 = getSpec(N2.NH3.cal.spec,DOASmodel=DOAS.model,lite=TRUE),
      ### SO2 calibration spectrum
      dat.SO2 = getSpec(SO2.cal.spec,DOASmodel=DOAS.model,lite=TRUE),
      dat.N2.SO2 = getSpec(N2.SO2.cal.spec,DOASmodel=DOAS.model,lite=TRUE),
      ### NO calibration spectrum
      dat.NO = getSpec(NO.cal.spec,DOASmodel=DOAS.model,lite=TRUE),
      dat.N2.NO = getSpec(N2.NO.cal.spec,DOASmodel=DOAS.model,lite=TRUE),
      ### cuvette calibration dark spectrum
      dat.N2.dark = getSpec(N2.dark.cal.spec,DOASmodel=DOAS.model,lite=TRUE)
    )
  )
}

correctSpectra <- function(specSet,rawData=NULL,correct.dark=TRUE,correct.linearity=TRUE,correct.straylight=c("avg","linear","none"),straylight.pix=NULL){
  
  if(is.null(straylight.pix)){
    straylight.pix <- getWindows(specSet$DOAS.model)$pixel3
  }
  ### dark-corrected reference spectra
  ### ******************************************************************************
  if(correct.dark){
    out <- list(
        I.N2.NH3 = specSet$dat.N2.NH3[[1]] - specSet$dat.N2.dark[[1]],
        I.N2.SO2 = specSet$dat.N2.SO2[[1]] - specSet$dat.N2.dark[[1]],
        I.N2.NO = specSet$dat.N2.NO[[1]] - specSet$dat.N2.dark[[1]],
        I.NH3 = specSet$dat.NH3[[1]] - specSet$dat.N2.dark[[1]],
        I.SO2 = specSet$dat.SO2[[1]] - specSet$dat.N2.dark[[1]],
        I.NO = specSet$dat.NO[[1]] - specSet$dat.N2.dark[[1]],
        I.ref = specSet$dat.ref[[1]] - specSet$dat.ref.dark[[1]],
        I.meas = rawData[[1]] - specSet$dat.dark[[1]]
      )
  } else {
    out <- list(
        I.N2.NH3 = specSet$dat.N2.NH3[[1]],
        I.N2.SO2 = specSet$dat.N2.SO2[[1]],
        I.N2.NO = specSet$dat.N2.NO[[1]],
        I.NH3 = specSet$dat.NH3[[1]],
        I.SO2 = specSet$dat.SO2[[1]],
        I.NO = specSet$dat.NO[[1]],
        I.ref = specSet$dat.ref[[1]],
        I.meas = rawData[[1]]
      )    
  }

  ### correct raw ccd data linearity
  ### ******************************************************************************
  if(correct.linearity){
    if(is.null(rawData)){
      lin.coef <- getDOASinfo(specSet$DOAS.model)$Spectrometer$"Linearity Coefficients"
    } else {
      lin.coef <- rawData$DOASinfo$Spectrometer$"Linearity Coefficients"
    }
    out$I.NH3 <- out$I.NH3 / linearity.func(out$I.NH3, lin.coef)
    out$I.SO2 <- out$I.SO2 / linearity.func(out$I.SO2, lin.coef)
    out$I.NO <- out$I.NO / linearity.func(out$I.NO, lin.coef)
    out$I.N2.NH3 <- out$I.N2.NH3 / linearity.func(out$I.N2.NH3, lin.coef)
    out$I.N2.SO2 <- out$I.N2.SO2 / linearity.func(out$I.N2.SO2, lin.coef)
    out$I.N2.NO <- out$I.N2.NO / linearity.func(out$I.N2.NO, lin.coef)
    out$I.ref <- out$I.ref / linearity.func(out$I.ref, lin.coef)
    out$I.meas <- out$I.meas / linearity.func(out$I.meas, lin.coef)
  }

  ### straylight-corrected reference spectra
  ### ******************************************************************************
  switch(pmatch(correct.straylight[1],c("none","avg","linear"),nomatch=4),
    {
      NULL
    },
    {
      out$I.N2.NH3 <- out$I.N2.NH3 - mean(out$I.N2.NH3[straylight.pix])
      out$I.N2.SO2 <- out$I.N2.SO2 - mean(out$I.N2.SO2[straylight.pix])
      out$I.N2.NO <- out$I.N2.NO - mean(out$I.N2.NO[straylight.pix])
      out$I.NH3 <- out$I.NH3 - mean(out$I.NH3[straylight.pix])
      out$I.SO2 <- out$I.SO2 - mean(out$I.SO2[straylight.pix])
      out$I.NO <- out$I.NO - mean(out$I.NO[straylight.pix])
      out$I.ref <- out$I.ref - mean(out$I.ref[straylight.pix])      
      if(!is.null(rawData))out$I.meas <- sweep(out$I.meas,2,colMeans(out$I.meas[straylight.pix,,drop=FALSE]),"-")    
    },
    {
      require(MASS)
      out$I.N2.NH3 <- out$I.N2.NH3 - predict(rlm(out$I.N2.NH3[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.N2.NH3))))
      out$I.N2.SO2 <- out$I.N2.SO2 - predict(rlm(out$I.N2.SO2[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.N2.SO2))))
      out$I.N2.NO <- out$I.N2.NO - predict(rlm(out$I.N2.NO[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.N2.NO))))
      out$I.NH3 <- out$I.NH3 - predict(rlm(out$I.NH3[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.NH3))))
      out$I.SO2 <- out$I.SO2 - predict(rlm(out$I.SO2[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.SO2))))
      out$I.NO <- out$I.NO - predict(rlm(out$I.NO[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.NO))))
      out$I.ref <- out$I.ref - predict(rlm(out$I.ref[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(out$I.ref))))      
      if(!is.null(rawData))out$I.meas <-apply(out$I.meas,2,function(x){x - predict(rlm(x[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=seq.int(length(x))))})         
    },
    stop("argument correct.straylight must be either \"none\", \"avg\" or \"linear\"")
    )

  if(!length(out$I.meas)) out$I.meas <- NULL 

  out$dat.ref <- specSet$dat.ref

  return(out)
}

diffSpecs <- function(Iset,use.ref=TRUE){

  suppressWarnings(
    out <- list(
      NH3.diffspec = log(Iset$I.NH3/Iset$I.N2.NH3),
      SO2.diffspec = log(Iset$I.SO2/Iset$I.N2.SO2),
      NO.diffspec = log(Iset$I.NO/Iset$I.N2.NO)
      )
    )

  if(!is.null(Iset$I.meas)){
    if(use.ref){
      suppressWarnings(out$diffspec <- log(Iset$I.meas/Iset$I.ref))
      out$diffspec[is.na(out$diffspec)] <- log(.Machine$double.eps)
      out$dat.ref <- Iset$dat.ref
    } else {
      out$diffspec <- log(Iset$I.meas)
    }
  }

  return(out)
}

getCalCurves <- function(diffspec,DOAS.win,calrefspec,warn=TRUE,...){

	# NH3:
	if(anyNA(c(calrefspec$dat.NH3$ambient$pathLength,calrefspec$dat.N2.NH3$ambient$pathLength))){
		if(warn)warning("NH3 calibration spectrum: unknown path lengths might be different!")
		NH3_total <- calrefspec$dat.NH3$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NH3$cuvette$cuvetteLength
	} else {
		# if(calrefspec$dat.NH3$ambient$pathLength != calrefspec$dat.N2.NH3$ambient$pathLength & anyNA(c(calrefspec$dat.NH3$ambient$NH3ambient_ug,calrefspec$dat.N2.NH3$ambient$NH3ambient_ug))){
    if(anyNA(c(calrefspec$dat.NH3$ambient$NH3ambient_ug,calrefspec$dat.N2.NH3$ambient$NH3ambient_ug))){
			if(warn)warning("NH3 calibration spectrum: unknown ambient concentrations!")
			NH3_total <- calrefspec$dat.NH3$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NH3$cuvette$cuvetteLength
		} else {
			NH3_total <- calrefspec$dat.NH3$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NH3$cuvette$cuvetteLength + 
				calrefspec$dat.NH3$ambient$pathLength*calrefspec$dat.NH3$ambient$NH3ambient_ug -
				calrefspec$dat.N2.NH3$ambient$pathLength*calrefspec$dat.N2.NH3$ambient$NH3ambient_ug
		}
	}
	NH3.dc <- highpass.filter(diffspec$NH3.diffspec,DOAS.win,...)/NH3_total


	# SO2:
	if(anyNA(c(calrefspec$dat.SO2$ambient$pathLength,calrefspec$dat.N2.SO2$ambient$pathLength))){
		if(warn)warning("SO2 calibration spectrum: unknown path lengths might be different!")
		SO2_total <- calrefspec$dat.SO2$cuvette$cuvetteConc_mg*1000*calrefspec$dat.SO2$cuvette$cuvetteLength
	} else {
		# if(calrefspec$dat.SO2$ambient$pathLength != calrefspec$dat.N2.SO2$ambient$pathLength & anyNA(c(calrefspec$dat.SO2$ambient$SO2ambient_ug,calrefspec$dat.N2.SO2$ambient$SO2ambient_ug))){
    if(anyNA(c(calrefspec$dat.SO2$ambient$SO2ambient_ug,calrefspec$dat.N2.SO2$ambient$SO2ambient_ug))){
			if(warn)warning("SO2 calibration spectrum: unknown ambient concentrations!")
			SO2_total <- calrefspec$dat.SO2$cuvette$cuvetteConc_mg*1000*calrefspec$dat.SO2$cuvette$cuvetteLength
		} else {
			SO2_total <- calrefspec$dat.SO2$cuvette$cuvetteConc_mg*1000*calrefspec$dat.SO2$cuvette$cuvetteLength + 
				calrefspec$dat.SO2$ambient$pathLength*calrefspec$dat.SO2$ambient$SO2ambient_ug -
				calrefspec$dat.N2.SO2$ambient$pathLength*calrefspec$dat.N2.SO2$ambient$SO2ambient_ug
		}
	}
	SO2.dc <- highpass.filter(diffspec$SO2.diffspec,DOAS.win,...)/SO2_total
	
	# NO:
	if(anyNA(c(calrefspec$dat.NO$ambient$pathLength,calrefspec$dat.N2.NO$ambient$pathLength))){
		if(warn)warning("NO calibration spectrum: unknown path lengths might be different!")
		NO_total <- calrefspec$dat.NO$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NO$cuvette$cuvetteLength
	} else {
		# if(calrefspec$dat.NO$ambient$pathLength != calrefspec$dat.N2.NO$ambient$pathLength & anyNA(c(calrefspec$dat.NO$ambient$NOambient_ug,calrefspec$dat.N2.NO$ambient$NOambient_ug))){
    if(anyNA(c(calrefspec$dat.NO$ambient$NOambient_ug,calrefspec$dat.N2.NO$ambient$NOambient_ug))){
			if(warn)warning("NO calibration spectrum: unknown ambient concentrations!")
			NO_total <- calrefspec$dat.NO$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NO$cuvette$cuvetteLength
		} else {
			NO_total <- calrefspec$dat.NO$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NO$cuvette$cuvetteLength + 
				calrefspec$dat.NO$ambient$pathLength*calrefspec$dat.NO$ambient$NOambient_ug -
				calrefspec$dat.N2.NO$ambient$pathLength*calrefspec$dat.N2.NO$ambient$NOambient_ug
		}
	}
	NO.dc <- highpass.filter(diffspec$NO.diffspec,DOAS.win,...)/NO_total

	return(list(
		NH3.dc = NH3.dc
		,SO2.dc = SO2.dc
		,NO.dc = NO.dc
		,Xreg = cbind(NH3=NH3.dc,SO2=SO2.dc,NO=NO.dc)
    ,NH3_total = NH3_total
    ,SO2_total = SO2_total
    ,NO_total = NO_total
		))
}

fitConc <- function(meas.doascurve, DOAS.win, path.length, Cal.dc, fit.type="ARIMA", robust=FALSE, 
	ARIMA.order=rbind(expand.grid(0:2,0,0,KEEP.OUT.ATTRS = FALSE),expand.grid(0,0,1:2,KEEP.OUT.ATTRS = FALSE)), dyn.fixed.pattern=NULL, fit.weights=NULL){
	ord <- as.character(as.numeric(NROW(ARIMA.order) > 1) + 1)[fit.type %in% "ARIMA"]
	rob <- ".rob"[robust]
	fitcurve <- get(paste0("fit.curves.",fit.type,ord,rob),mode="function")
	if(is.null(dyn.fixed.pattern))dyn.fixed.pattern <- meas.doascurve[DOAS.win$pixel4]*0
	return(fitcurve(meas.doascurve,DOAS.win$pixel4,dyn.fixed.pattern,Cal.dc$Xreg[DOAS.win$pixel4,],fit.weights,DOAS.win$tau.shift,ARIMA.order,path.length))
}


inspectEvaluation <- function(rawdat,CalRefSpecs, path.length, index = 1,
  filter.type = "BmHarris", fit.type = c("OLS","ARIMA"),doubleAVG=TRUE, robust = TRUE, straylight.window = NULL, filter.window = NULL, fit.window = NULL, 
  filter.strength = NULL, tau.shift = NULL, special.control = list(b = 3.5, 
      Scale = function(r) median(abs(r))/0.6745, delta = 0, 
      filter.strength_multiplication = 1.25, filter.strength_loess = 0.2, 
      fam = "gaussian"),correct.dark = TRUE, correct.linearity = TRUE, 
  correct.straylight = c("avg", "linear", "none"), use.ref = TRUE,
  Edinburgh_correction = TRUE){
  require(shiny)
  
  # # cal/ref specs:
  # CalRefSpecs <- getSpecSet(spec.dir, ref.spec, ref.dark.spec, 
  #     dark.spec, NH3.cal.spec, SO2.cal.spec, 
  #     NO.cal.spec, N2.NH3.cal.spec, N2.SO2.cal.spec, 
  #     N2.NO.cal.spec, N2.dark.cal.spec, 
  #     DOAS.model=rawdat$DOASinfo$DOASmodel
  #     )
  # require(data.table)
  require(IDPmisc)
  require(lubridate)
  require(robustbase)

  # windows first time:
  DOASwindows <- getWindows(rawdat$DOASinfo,filter.type = filter.type, straylight.window = straylight.window, 
    filter.window = filter.window, fit.window = fit.window, filter.strength = filter.strength, double = doubleAVG,
    tau.shift = tau.shift, special.control = special.control)

  # correct cal/ref specs:
  # CalRefCorr <- correctSpectra(CalRefSpecs)
  SpecCorr <- correctSpectra(CalRefSpecs,rawdat, correct.dark = correct.dark, correct.linearity = correct.linearity, 
    correct.straylight = correct.straylight, straylight.pix=DOASwindows$pixel3)

  # diffspec:
  # DiffSpec <- diffSpecs(SpecCorr)
  DiffSpec <- diffSpecs(SpecCorr,use.ref=use.ref)

  # cal doascurves:
  Cal.dc <- getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE)

  wavelength <- rawdat$DOASinfo$Spectrometer$wavelength[DOASwindows$pixel1]

  index.max <- ncol(rawdat$RawData)
  
  frmls <- switch(DOASwindows$filter.type,
    "Rect" = formals(winRect),
    "Hann" = formals(winHann),
    "Hamming" = formals(winHamming),
    "Blackman" = formals(winBlackman),
    "BmNuttall" = formals(winBmNuttall),
    "FlatTop" = formals(winFlatTop),
    "Sin" = formals(winSin),
    "Gauss" = formals(winGauss),
    "Kaiser" = formals(winKaiser),
    "DolphChebyshev" = formals(winDolphChebyshev),
    "BmHarris" = formals(winBmHarris),
    "Tukey" = formals(winTukey),
    "Poisson" = formals(winPoisson),
    "Exp" = formals(winExp),
    "ExpHamming" = formals(winExpHamming),
    NULL
  )
  n <- DOASwindows$filter.strength
  p1 <- eval(frmls$p1)
  p2 <- eval(frmls$p2)
  p3 <- eval(frmls$p3) 
  p4 <- eval(frmls$p4)

  # runApp:
  runApp(
      list(
      # ui:
        ui={
        fluidPage(

          # Application title
          titlePanel(paste0("spectra analysis - miniDOAS ",rawdat$DOASinfo$DOASmodel)),

          sidebarLayout(position="right",

            sidebarPanel(width=2,
              numericInput("index", 
                label = "index:",
                step = 1L,
                min = 1L, 
                max = index.max, 
                value = index
              )             
              ,uiOutput("filter_strength")
              # ,numericInput("filter.strength", 
              #   label = "filter.strength:",
              #   step = 0.01,
              #   value = DOASwindows$filter.strength
              # )
              ,selectInput("filter.type", 
                label = "filter.type:",
                choices = c("Hamming","Rect","Hann","Blackman","BmNuttall","Sin"
                  ,"Gauss","FlatTop","Kaiser","DolphChebyshev","BmHarris","Tukey","Poisson","Exp","ExpHamming","special","loess"),
                selected = DOASwindows$filter.type
              )
              ,selectInput("double", 
                label = "doubleAVG:",
                choices = c(TRUE,FALSE),
                selected = DOASwindows$double
              )
              ,sliderInput("filter.window", 
                label = "filter.window:", 
                min = 190, 
                max = 240,
                step = 0.1, 
                value = DOASwindows$filter.window
              )
              ,uiOutput("fit_window")
              # ,sliderInput("fit.window", 
              #   label = "fit.window:", 
              #   min = filter_window[1], 
              #   max = filter_window[2],
       #          step = 0.1, 
              #   value = DOASwindows$fit.window
              # )
              ,selectInput("fit.type", 
                label = "fit.type:",
                choices = c("OLS","ARIMA"),
                selected = fit.type
              )
              ,checkboxInput("robust", 
                label = "robust:",
                value = robust
              )
              ,numericInput("tau.shift", 
                label = "tau.shift:", 
                value = DOASwindows$tau.shift,
                step = 1L,
                min = -20L,
                max = 20L
              )
              # ,textInput("Scale", 
              #   label = "Scale:", 
              #   value = paste0(deparse(DOASwindows$special.control$Scale,width.cutoff=500),collapse=" ")
              # )
              ,uiOutput("filterArgs")
              ,checkboxInput("correct.dark", 
                label = "correct.dark:",
                value = correct.dark
              )
              ,checkboxInput("correct.linearity", 
                label = "correct.linearity:",
                value = correct.linearity
              )
              ,selectInput("correct.straylight", 
                label = "correct.straylight:",
                choices = c("avg", "linear", "none"),
                selected = correct.straylight[1]
              )
              ,checkboxInput("use.ref", 
                label = "use.ref:",
                value = use.ref
              )
              ,downloadButton('downloadPlot', 'Save Plot to png')                                                           
  
            ),

            mainPanel(
              plotOutput("specs",height="800px")
            )
          )
            )
        },
      # server
      server={
        function(input, output, session) {

          # # observe DOASwindows:
          # observe({
          #   DOASwindows <- DOASwindows_reactive()
          #   # SpecCorr <- SpecCorr_reactive()
          #   # DiffSpec <- DiffSpec_reactive()
          #   # Cal.dc <- Cal.dc_reactive()
          # })


          # windows:
          DOASwindows_reactive <- reactive({
            if(is.null(input$fit.window)){
              fit_window <- DOASwindows$fit.window
            } else {
              fit_window <- input$fit.window
            }
            if(is.null(input$filter.strength)){
              filter.strength <- DOASwindows$filter.strength
            } else {
              filter.strength <- input$filter.strength
            }
            if(is.null(input$filter.strength_multiplication)){
              filter.strength_multiplication <- DOASwindows$special.control$filter.strength_multiplication
            } else {
              filter.strength_multiplication <- input$filter.strength_multiplication
            }
            if(is.null(input$filter.strength_loess)){
              filter.strength_loess <- DOASwindows$special.control$filter.strength_loess
            } else {
              filter.strength_loess <- input$filter.strength_loess
            }
            if(!(input$filter.type %in% c("special","loess")) & filter.strength < 3){
              if(rawdat$DOASinfo$DOASmodel=="S1" & rawdat$DOASinfo$timerange < parse_date_time("20170101","Ymd")){
                filter_strength <- 41
              } else {
                filter_strength <- 25
              }
            } else if(input$filter.type %in% c("special","loess") & filter.strength > 1){
              filter_strength <- 0.17
            } else {
              filter_strength <- filter.strength
            }
            getWindows(rawdat$DOASinfo, filter.type = input$filter.type, straylight.window = straylight.window, 
              filter.window = input$filter.window, fit.window = fit_window, double = input$double,
              filter.strength = filter_strength, tau.shift = input$tau.shift, special.control = list(b = special.control$b, 
                Scale = special.control$Scale, 
                filter.strength_multiplication = filter.strength_multiplication, filter.strength_loess = filter.strength_loess, 
                fam = special.control$fam))
          })

          # # correct cal/ref specs:
          # SpecCorr_reactive <- reactive({
          #   correctSpectra(CalRefSpecs,rawdat,correct.dark = input$correct.dark, correct.linearity = input$correct.linearity, 
      #           correct.straylight = input$correct.straylight, straylight.pix=DOASwindows$pixel3)
          # })

          # # diffspec:
          # DiffSpec_reactive <- reactive({
          #   diffSpecs(SpecCorr,use.ref=input$use.ref)
          # })

          # # cal doascurves:
          # Cal.dc_reactive <- reactive({
          #   getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE)
          # })
          output$filterArgs <- renderUI({
            DOASwindows <- DOASwindows_reactive()
            if(is.null(input$filter.strength_multiplication)){
              filter.strength_multiplication <- DOASwindows$special.control$filter.strength_multiplication
            } else {
              filter.strength_multiplication <- input$filter.strength_multiplication
            }
            if(is.null(input$filter.strength_loess)){
              filter.strength_loess <- DOASwindows$special.control$filter.strength_loess
            } else {
              filter.strength_loess <- input$filter.strength_loess
            }
            if(DOASwindows$filter.type %in% c("special","loess")){
              tagList(
                  numericInput("filter.strength_multiplication", 
                  label = "filter.strength_multiplication:",
                  step = 0.05,
                  value = filter.strength_multiplication
                )
                ,numericInput("filter.strength_loess", 
                  label = "filter.strength_loess:",
                  step = 0.01,
                  value = filter.strength_loess
                )
              )
            } else {
              if(is.null(input$filter.strength)){
                filter.strength <- DOASwindows$filter.strength
              } else {
                filter.strength <- input$filter.strength
              }
              frmls <- switch(DOASwindows$filter.type,
                "Rect" = formals(winRect),
                "Hann" = formals(winHann),
                "Hamming" = formals(winHamming),
                "Blackman" = formals(winBlackman),
                "BmNuttall" = formals(winBmNuttall),
                "FlatTop" = formals(winFlatTop),
                "Sin" = formals(winSin),
                "Gauss" = formals(winGauss),
                "Kaiser" = formals(winKaiser),
                "DolphChebyshev" = formals(winDolphChebyshev),
                "BmHarris" = formals(winBmHarris),
                "Tukey" = formals(winTukey),
                "Poisson" = formals(winPoisson),
                "Exp" = formals(winExp),
                "ExpHamming" = formals(winExpHamming),
                NULL
              )
              n <- filter.strength
              p1 <- eval(frmls$p1)
              p2 <- eval(frmls$p2)
              p3 <- eval(frmls$p3)
              p4 <- eval(frmls$p4)

              tagList(
                if(!is.null(p1))numericInput("p1", 
                  label = "p1:",
                  step = 0.1,
                  value = p1
                )
                ,if(!is.null(p2))numericInput("p2", 
                  label = "p2:",
                  step = 0.1,
                  value = p2
                )
                ,if(!is.null(p3))numericInput("p3", 
                  label = "p3:",
                  step = 0.1,
                  value = p3
                )
                ,if(!is.null(p4))numericInput("p4", 
                  label = "p4:",
                  step = 0.1,
                  value = p4
                )
              )
            }            
          })

          output$filter_strength <- renderUI({
            DOASwindows <- DOASwindows_reactive()

            if(DOASwindows$filter.type %in% c("special","loess")){
              numericInput("filter.strength", 
                label = "filter.strength:",
                step = 0.01,
                value = DOASwindows$filter.strength
              )
            } else {
              numericInput("filter.strength", 
                label = "filter.strength:",
                step = 2,
                value = DOASwindows$filter.strength
              )
            }
          })

          output$fit_window <- renderUI({
            DOASwindows <- DOASwindows_reactive()
            sliderInput("fit.window", 
              label = "fit.window:", 
              min = DOASwindows$filter.window[1], 
              max = DOASwindows$filter.window[2],
              step = 0.1, 
              value = DOASwindows$fit.window
            )
          })


          plot_reactive <- reactive({
            # DOASwindows <- getWindows(rawdat$DOASinfo, filter.type = input$filter.type, straylight.window = straylight.window, 
            #   filter.window = filter_window, fit.window = input$fit.window, filter.strength = filter_strength, 
            #   tau.shift = input$tau.shift, special.control = list(b = special.control$b, Scale = special.control$Scale, 
            #     filter.strength_multiplication = input$filter.strength_multiplication, filter.strength_loess = input$filter.strength_loess, 
            #     fam = special.control$fam)
            # )
            i <- input$index
            RawDat <- rawdat
            RawDat[[1]] <- rawdat[[1]][,i,drop=FALSE]
            DOASwindows <- DOASwindows_reactive()
            SpecCorr <- correctSpectra(CalRefSpecs,RawDat,correct.dark = input$correct.dark, correct.linearity = input$correct.linearity, 
                correct.straylight = input$correct.straylight, straylight.pix=DOASwindows$pixel3)
            DiffSpec <- diffSpecs(SpecCorr,use.ref=input$use.ref)
            Cal.dc <- getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE,input$p1,input$p2,input$p3,input$p4)
            wavelength <- rawdat$DOASinfo$Spectrometer$wavelength[DOASwindows$pixel1]

            
            # plot(1,main="diffspec")
            # plot(1,main="doascurve")
            # plot(1,main="darkspec")
            # reactive:
            xlim <- range(wavelength)
            # main <- rawdat$Header[i,]

            # I.meas + I.ref
            ylim <- range(SpecCorr$I.meas[DOASwindows$pixel1,],SpecCorr$I.ref[DOASwindows$pixel1],na.rm=TRUE)
            if(!all(is.finite(ylim))){
              ylim <- c(0.01,2)
            } else if(any(ylimBelow <- ylim<0)){
              ylim[ylimBelow] <- c(0.01,2)[ylimBelow]
            } 
            # windows(width=10,height=7)
            plot(1,1,xlim=xlim,ylim=ylim,log="y",type="n",ylab="counts",xlab="",main="spectra")
            if(input$use.ref)lines(wavelength,SpecCorr$I.ref[DOASwindows$pixel1],col="darkgrey",lwd=2)
            lines(wavelength,SpecCorr$I.meas[DOASwindows$pixel1,],col="black",lwd=2)
            legend("bottomright",c("meas.","ref."),lwd=2,col=c("black","darkgrey"),bty="n")

            # log(I.meas/I.ref)
            ylim <- range(DiffSpec$diffspec[DOASwindows$pixel1,],na.rm=TRUE)
            if(!all(is.finite(ylim))) ylim <- c(0,1)
            meas.dc <- highpass.filter(DiffSpec$diffspec,DOASwindows,input$p1,input$p2,input$p3,input$p4)
            # isna <- is.na(meas.dc)
            # windows(width=10,height=7)
            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.meas/I.ref)",xlab="",main="diffspec")
            lines(wavelength,DiffSpec$diffspec[DOASwindows$pixel1,] - meas.dc,lwd=2,col="orange")
            lines(wavelength,DiffSpec$diffspec[DOASwindows$pixel1,],col="black")
            fit <- fitConc(meas.dc, DOASwindows, path.length, Cal.dc, fit.type=input$fit.type, robust=input$robust)
            
            cfs <- fit[[5]]
            fit.SO2 <- cfs[2]*path.length*Cal.dc$SO2.dc
            fit.NO <- cfs[3]*path.length*Cal.dc$NO.dc
            fit.NH3 <- cfs[1]*path.length*Cal.dc$NH3.dc

            # fit.SO2[isna] <- NA
            # fit.NO[isna] <- NA
            # fit.NH3[isna] <- NA

            if(input$tau.shift>0){
              meas.dc <- c(meas.dc[-seq.int(input$tau.shift)],rep(NA,abs(input$tau.shift)))
            } else if(input$tau.shift<0){
              tau.shift <- abs(input$tau.shift)
              meas.dc <- c(rep(NA,tau.shift),meas.dc[-(length(meas.dc) - seq.int(tau.shift) + 1)])
            }

            # doascurve
            ylim <- range(meas.dc,na.rm=TRUE)
            if(!all(is.finite(ylim))) ylim <- c(0,1)
            # windows(width=10,height=7)
            plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)},main="doascurve")
            lines(wavelength,fit.SO2,col="#00bb00aa")
            lines(wavelength,fit.SO2+fit.NO,col="#ff0000dd")
            lines(wavelength,fit.SO2+fit.NO+fit.NH3,col="blue")
            lines(wavelength[DOASwindows$pixel4],fit.SO2[DOASwindows$pixel4],lwd=2,col="#00bb00aa")
            lines(wavelength[DOASwindows$pixel4],fit.SO2[DOASwindows$pixel4]+fit.NO[DOASwindows$pixel4],lwd=2,col="#ff0000dd")
            lines(wavelength[DOASwindows$pixel4],fit.SO2[DOASwindows$pixel4]+fit.NO[DOASwindows$pixel4]+fit.NH3[DOASwindows$pixel4],lwd=2,col="blue")
            abline(v=wavelength[range(DOASwindows$pixel4)],lty=3)
            legend("topright",c("SO2","SO2 + NO","SO2 + NO + NH3"),lty=1,col=c("#00bb00aa","#ff0000dd","blue"),cex=0.7)

            # residual
            # windows(width=10,height=7)
            plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="residuals [-]",xlab="",panel.first={grid();abline(h=0)},col="darkgrey",main="residuals")
            lines(wavelength,meas.dc-fit.SO2-fit.NO-fit.NH3) 
            lines(wavelength[DOASwindows$pixel4],meas.dc[DOASwindows$pixel4]-fit.SO2[DOASwindows$pixel4]-fit.NO[DOASwindows$pixel4]-fit.NH3[DOASwindows$pixel4],lwd=2) 
            abline(v=wavelength[range(DOASwindows$pixel4)],lty=3)

            # cal diffspec
            ylim <- range(DiffSpec$NH3.diffspec[DOASwindows$pixel1],na.rm=TRUE)
            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.NH3/I.N2)",xlab="",panel.first={grid();abline(h=0)})
            # #
            # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel1] - Cal.dc$Xreg[,3]*CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NO$cuvette$cuvetteLength,lwd=2,col="red")
            # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel1])
            # #
            # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel1] - Cal.dc$Xreg[,2]*CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.SO2$cuvette$cuvetteLength,lwd=2,col="green")
            # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel1])
            # #
            lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel1] - Cal.dc$Xreg[,1]*CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NH3$cuvette$cuvetteLength,lwd=2,col="blue")
            lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel1])

            # cal doascurve
            ylim <- range(Cal.dc$Xreg,na.rm=TRUE)
            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="cal. DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)})
            lines(wavelength,Cal.dc$Xreg[,3],lwd=2,col="red")
            lines(wavelength,Cal.dc$Xreg[,2],lwd=2,col="green")
            lines(wavelength,Cal.dc$Xreg[,1],lwd=2,col="blue")

            ## test stats:
            nh3wts <- abs(Cal.dc$NH3.dc[DOASwindows$pixel4])/sum(abs(Cal.dc$NH3.dc[DOASwindows$pixel4]))
            so2wts <- abs(Cal.dc$SO2.dc[DOASwindows$pixel4])/sum(abs(Cal.dc$SO2.dc[DOASwindows$pixel4]))
            nowts <- abs(Cal.dc$NO.dc[DOASwindows$pixel4])/sum(abs(Cal.dc$NO.dc[DOASwindows$pixel4]))
            gof.NH3 <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*nh3wts)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nh3wts)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nh3wts),2)
            gof.SO2 <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*so2wts)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*so2wts)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*so2wts),2)
            gof.NO <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*nowts)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nowts)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nowts),2)
            gof.perc <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*(nh3wts+so2wts+nowts)/3)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*(nh3wts+so2wts+nowts)/3)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*(nh3wts+so2wts+nowts)/3),2)
      
            SAE <- sum(abs(fit[[8]]))
            SAT <- sum(abs(fit[[7]]+fit[[8]]))
            SAF <- sum(abs(fit[[7]]))

            SRE <- sum(fit[[8]])
            STAT1 <- round(SRE/SAE,2)
            STAT2 <- round(SAE/SAT,2)
            STAT3 <- round(SAE/SAF,2)

            OLS <- if(input$fit.type=="ARIMA"){
              fit[[4]][grep("0/0/0",fit[[4]][,1]),2]
            } else {
              NA
            }

            if (Edinburgh_correction) {
                cfs <- cfs * 1.16
            }

            # titel:
            # title(paste0(format(rawdat$Header[i,"st"]),sprintf("  index: %i/%i  -  NH3: %1.1f +/- %1.1f, SO2: %1.1f +/- %1.1f, NO: %1.1f +/- %1.1f",i,index.max,cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3])),outer=TRUE,line=-1)          
            # title(paste0(format(rawdat$Header[i,"st"]),sprintf("  index: %i/%i  -  NH3: %1.1f +/- %1.1f, SO2: %1.1f +/- %1.1f, NO: %1.1f +/- %1.1f  SRE/SAE: %1.2f, SAE/SAF: %1.2f, SAE/SAT: %1.2f, NH3-NH3_000: %1.2f, (NH3-NH3_000)/NH3: %1.2f",i,index.max,cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3],STAT1,STAT2,STAT3,cfs[1]-OLS,(cfs[1]-OLS)/cfs[1])),outer=TRUE,line=-1)          

            msg1 <- sprintf("index %i/%i:  %s  --  NH3: %1.1f +/- %1.1f  --  SO2: %1.1f +/- %1.1f  --  NO: %1.1f +/- %1.1f",i,index.max,format(rawdat$Header[i,"st"]),cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3])
            msg2 <- sprintf("NH3-NH3_000: %1.2f, (NH3-NH3_000)/NH3: %1.2f -- pseudo GOF: -total: %1.2f -NH3: %1.2f -SO2: %1.2f -NO: %1.2f -- SRE/SAE: %1.2f, SAE/SAF: %1.2f, SAE/SAT: %1.2f",cfs[1]-OLS,(cfs[1]-OLS)/cfs[1],gof.perc,gof.NH3,gof.SO2,gof.NO,STAT1,STAT2,STAT3)
            # mtext(as.expression(substitute(italic(msg), list(msg=msg))), line=-1.25, outer=TRUE, cex=0.5)
            # mtext(substitute(italic(msg), list(msg=msg1)), line=1, outer=TRUE, cex=0.75)
            # mtext(substitute(italic(msg), list(msg=msg2)), line=-0.5, outer=TRUE, cex=0.75)
            mtext(msg1, line=0, outer=TRUE)#, cex=0.75)
            mtext(msg2, line=-1.5, outer=TRUE, cex=0.75)

          })

          output$specs <- renderPlot({
            par(mfrow=c(3,2),oma=c(0,0,1.5,0))
            plot_reactive()
          })
          output$downloadPlot <- downloadHandler(
            filename=function(){paste0("shinyInspect_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png")},
            content = function(file){
              png(file,width=1200,height=1200)
              on.exit(dev.off())

            par(mfrow=c(3,2),oma=c(0,0,1.5,0))
            # plot_reactive()
           # DOASwindows <- getWindows(rawdat$DOASinfo, filter.type = input$filter.type, straylight.window = straylight.window, 
            #   filter.window = filter_window, fit.window = input$fit.window, filter.strength = filter_strength, 
            #   tau.shift = input$tau.shift, special.control = list(b = special.control$b, Scale = special.control$Scale, 
            #     filter.strength_multiplication = input$filter.strength_multiplication, filter.strength_loess = input$filter.strength_loess, 
            #     fam = special.control$fam)
            # )
            i <- input$index
            RawDat <- rawdat
            RawDat[[1]] <- rawdat[[1]][,i,drop=FALSE]
            DOASwindows <- DOASwindows_reactive()
            SpecCorr <- correctSpectra(CalRefSpecs,RawDat,correct.dark = input$correct.dark, correct.linearity = input$correct.linearity, 
                correct.straylight = input$correct.straylight, straylight.pix=DOASwindows$pixel3)
            DiffSpec <- diffSpecs(SpecCorr,use.ref=input$use.ref)
            Cal.dc <- getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE,input$p1,input$p2,input$p3,input$p4)
            wavelength <- rawdat$DOASinfo$Spectrometer$wavelength[DOASwindows$pixel1]

            
            # plot(1,main="diffspec")
            # plot(1,main="doascurve")
            # plot(1,main="darkspec")
            # reactive:
            xlim <- range(wavelength)
            # main <- rawdat$Header[i,]

            # I.meas + I.ref
            ylim <- range(SpecCorr$I.meas[DOASwindows$pixel1,],SpecCorr$I.ref[DOASwindows$pixel1],na.rm=TRUE)
            if(!all(is.finite(ylim))){
              ylim <- c(0.01,2)
            } else if(any(ylimBelow <- ylim<0)){
              ylim[ylimBelow] <- c(0.01,2)[ylimBelow]
            } 
            # windows(width=10,height=7)
            plot(1,1,xlim=xlim,ylim=ylim,log="y",type="n",ylab="counts",xlab="",main="spectra")
            if(input$use.ref)lines(wavelength,SpecCorr$I.ref[DOASwindows$pixel1],col="darkgrey",lwd=2)
            lines(wavelength,SpecCorr$I.meas[DOASwindows$pixel1,],col="black",lwd=2)
            legend("bottomright",c("meas.","ref."),lwd=2,col=c("black","darkgrey"),bty="n")

            # log(I.meas/I.ref)
            ylim <- range(DiffSpec$diffspec[DOASwindows$pixel1,],na.rm=TRUE)
            if(!all(is.finite(ylim))) ylim <- c(0,1)
            meas.dc <- highpass.filter(DiffSpec$diffspec,DOASwindows,input$p1,input$p2,input$p3,input$p4)
            # isna <- is.na(meas.dc)
            # windows(width=10,height=7)
            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.meas/I.ref)",xlab="",main="diffspec")
            lines(wavelength,DiffSpec$diffspec[DOASwindows$pixel1,] - meas.dc,lwd=2,col="orange")
            lines(wavelength,DiffSpec$diffspec[DOASwindows$pixel1,],col="black")
            fit <- fitConc(meas.dc, DOASwindows, path.length, Cal.dc, fit.type=input$fit.type, robust=input$robust)
            
            cfs <- fit[[5]]
            fit.SO2 <- cfs[2]*path.length*Cal.dc$SO2.dc
            fit.NO <- cfs[3]*path.length*Cal.dc$NO.dc
            fit.NH3 <- cfs[1]*path.length*Cal.dc$NH3.dc

            # fit.SO2[isna] <- NA
            # fit.NO[isna] <- NA
            # fit.NH3[isna] <- NA

            if(input$tau.shift>0){
              meas.dc <- c(meas.dc[-seq.int(input$tau.shift)],rep(NA,abs(input$tau.shift)))
            } else if(input$tau.shift<0){
              tau.shift <- abs(input$tau.shift)
              meas.dc <- c(rep(NA,tau.shift),meas.dc[-(length(meas.dc) - seq.int(tau.shift) + 1)])
            }

            # doascurve
            ylim <- range(meas.dc,na.rm=TRUE)
            if(!all(is.finite(ylim))) ylim <- c(0,1)
            # windows(width=10,height=7)
            plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)},main="doascurve")
            lines(wavelength,fit.SO2,col="#00bb00aa")
            lines(wavelength,fit.SO2+fit.NO,col="#ff0000dd")
            lines(wavelength,fit.SO2+fit.NO+fit.NH3,col="blue")
            lines(wavelength[DOASwindows$pixel4],fit.SO2[DOASwindows$pixel4],lwd=2,col="#00bb00aa")
            lines(wavelength[DOASwindows$pixel4],fit.SO2[DOASwindows$pixel4]+fit.NO[DOASwindows$pixel4],lwd=2,col="#ff0000dd")
            lines(wavelength[DOASwindows$pixel4],fit.SO2[DOASwindows$pixel4]+fit.NO[DOASwindows$pixel4]+fit.NH3[DOASwindows$pixel4],lwd=2,col="blue")
            abline(v=wavelength[range(DOASwindows$pixel4)],lty=3)
            legend("topright",c("SO2","SO2 + NO","SO2 + NO + NH3"),lty=1,col=c("#00bb00aa","#ff0000dd","blue"),cex=0.7)

            # residual
            # windows(width=10,height=7)
            plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="residuals [-]",xlab="",panel.first={grid();abline(h=0)},col="darkgrey",main="residuals")
            lines(wavelength,meas.dc-fit.SO2-fit.NO-fit.NH3) 
            lines(wavelength[DOASwindows$pixel4],meas.dc[DOASwindows$pixel4]-fit.SO2[DOASwindows$pixel4]-fit.NO[DOASwindows$pixel4]-fit.NH3[DOASwindows$pixel4],lwd=2) 
            abline(v=wavelength[range(DOASwindows$pixel4)],lty=3)

            # cal diffspec
            ylim <- range(DiffSpec$NH3.diffspec[DOASwindows$pixel1],na.rm=TRUE)
            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.NH3/I.N2)",xlab="",panel.first={grid();abline(h=0)})
            # #
            # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel1] - Cal.dc$Xreg[,3]*CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NO$cuvette$cuvetteLength,lwd=2,col="red")
            # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel1])
            # #
            # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel1] - Cal.dc$Xreg[,2]*CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.SO2$cuvette$cuvetteLength,lwd=2,col="green")
            # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel1])
            # #
            lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel1] - Cal.dc$Xreg[,1]*CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NH3$cuvette$cuvetteLength,lwd=2,col="blue")
            lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel1])

            # cal doascurve
            ylim <- range(Cal.dc$Xreg,na.rm=TRUE)
            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="cal. DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)})
            lines(wavelength,Cal.dc$Xreg[,3],lwd=2,col="red")
            lines(wavelength,Cal.dc$Xreg[,2],lwd=2,col="green")
            lines(wavelength,Cal.dc$Xreg[,1],lwd=2,col="blue")

            ## test stats:
            nh3wts <- abs(Cal.dc$NH3.dc[DOASwindows$pixel4])/sum(abs(Cal.dc$NH3.dc[DOASwindows$pixel4]))
            so2wts <- abs(Cal.dc$SO2.dc[DOASwindows$pixel4])/sum(abs(Cal.dc$SO2.dc[DOASwindows$pixel4]))
            nowts <- abs(Cal.dc$NO.dc[DOASwindows$pixel4])/sum(abs(Cal.dc$NO.dc[DOASwindows$pixel4]))
            gof.NH3 <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*nh3wts)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nh3wts)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nh3wts),2)
            gof.SO2 <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*so2wts)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*so2wts)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*so2wts),2)
            gof.NO <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*nowts)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nowts)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*nowts),2)
            gof.perc <- round((sum(na.rm=TRUE,(fit[[7]]-1)/(meas.dc[DOASwindows$pixel4]-1)*(nh3wts+so2wts+nowts)/3)*(sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*(nh3wts+so2wts+nowts)/3)-1)+1)/sum(na.rm=TRUE,meas.dc[DOASwindows$pixel4]*(nh3wts+so2wts+nowts)/3),2)
      
            SAE <- sum(abs(fit[[8]]))
            SAT <- sum(abs(fit[[7]]+fit[[8]]))
            SAF <- sum(abs(fit[[7]]))

            SRE <- sum(fit[[8]])
            STAT1 <- round(SRE/SAE,2)
            STAT2 <- round(SAE/SAT,2)
            STAT3 <- round(SAE/SAF,2)

            OLS <- if(input$fit.type=="ARIMA"){
              fit[[4]][grep("0/0/0",fit[[4]][,1]),2]
            } else {
              NA
            }

            # titel:
            # title(paste0(format(rawdat$Header[i,"st"]),sprintf("  index: %i/%i  -  NH3: %1.1f +/- %1.1f, SO2: %1.1f +/- %1.1f, NO: %1.1f +/- %1.1f",i,index.max,cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3])),outer=TRUE,line=-1)          
            # title(paste0(format(rawdat$Header[i,"st"]),sprintf("  index: %i/%i  -  NH3: %1.1f +/- %1.1f, SO2: %1.1f +/- %1.1f, NO: %1.1f +/- %1.1f  SRE/SAE: %1.2f, SAE/SAF: %1.2f, SAE/SAT: %1.2f, NH3-NH3_000: %1.2f, (NH3-NH3_000)/NH3: %1.2f",i,index.max,cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3],STAT1,STAT2,STAT3,cfs[1]-OLS,(cfs[1]-OLS)/cfs[1])),outer=TRUE,line=-1)          

            msg1 <- sprintf("index %i/%i:  %s  --  NH3: %1.1f +/- %1.1f  --  SO2: %1.1f +/- %1.1f  --  NO: %1.1f +/- %1.1f",i,index.max,format(rawdat$Header[i,"st"]),cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3])
            msg2 <- sprintf("NH3-NH3_000: %1.2f, (NH3-NH3_000)/NH3: %1.2f -- pseudo GOF: -total: %1.2f -NH3: %1.2f -SO2: %1.2f -NO: %1.2f -- SRE/SAE: %1.2f, SAE/SAF: %1.2f, SAE/SAT: %1.2f",cfs[1]-OLS,(cfs[1]-OLS)/cfs[1],gof.perc,gof.NH3,gof.SO2,gof.NO,STAT1,STAT2,STAT3)
            # mtext(as.expression(substitute(italic(msg), list(msg=msg))), line=-1.25, outer=TRUE, cex=0.5)
            # mtext(substitute(italic(msg), list(msg=msg1)), line=1, outer=TRUE, cex=0.75)
            # mtext(substitute(italic(msg), list(msg=msg2)), line=-0.5, outer=TRUE, cex=0.75)
            mtext(msg1, line=0, outer=TRUE)#, cex=0.75)
            mtext(msg2, line=-1.5, outer=TRUE, cex=0.75)
 
            },
            contentType="image/png"
            )
        }
      } 
    )
  )
}


fitparallel <- function(index,DiffSpec,DOAS.win,meas.doascurve,dyn.fixed.pattern,Cal.dc,arima.Order,path.length,
  isna,best.tau,delta.AICc.zero,best.order,aicctab,coeffs,se,fitted.doascurve,residual.best){

  for(i in index){

    if(median(DiffSpec$diffspec[,i],na.rm=TRUE) > -5){
      # highpass filtering:
      meas.doascurve[,i] <- highpass.filter2(DiffSpec$diffspec[DOAS.win$pixel1,i],DOAS.win)
      
      ### fit calibration curves to measured DOAS curve (see Stutz & Platt, 1996, Applied Optics 30), determine best fit after shifting over given tau range, no fixed pattern considered
      ### ******************************************************************************#
      fitted <- fitcurve(meas.doascurve[,i], DOAS.win$pixel4, dyn.fixed.pattern, Cal.dc$Xreg[DOAS.win$pixel4,], NULL, DOAS.win$tau.shift, arima.Order, path.length)

      if(is.null(fitted)){
        isna[i] <- TRUE
      } else {
        best.tau[i] <- fitted[[1]]
        delta.AICc.zero[i] <- fitted[[2]]
        best.order[i] <- fitted[[3]]
        aicctab[[i]] <- fitted[[4]]
        coeffs[,i] <- fitted[[5]]
        se[,i] <- fitted[[6]]
        fitted.doascurve[,i] <- fitted[[7]]
        residual.best[,i] <- fitted[[8]]
      }
    } else {
      isna[i] <- TRUE
    }
  }

  out <- list(meas.doascurve,isna,best.tau,delta.AICc.zero,best.order,aicctab,coeffs,se,fitted.doascurve,residual.best)
  return(out)
}





if(FALSE){  
  DOAS.model=DOASmodel
  evalperiod=evalPeriod
  tz.DOAS=tzDOAS
  path.length=pathLength[DOASmodel]
  tz.Output=tzOutput
  filter.type="BmHarris"
  use.arima=useArima
  use.robust=useRobust
  filter.window=NULL 
  filter.strength=NULL
  double.AVG=TRUE
  fit.window=NULL
  straylight.window=NULL
  skip.check.daily=FALSE
  corr.fixed.pattern=FALSE
  fixed.pattern.length=60*24-1
  correct.linearity=TRUE
  correct.straylight=c("avg","linear","none")
  correct.dark=TRUE
  use.ref = TRUE
  plot.results=FALSE
  avg.period=1
  plot.cal=FALSE
  tau.shift=tauShift[DOASmodel]
  save.dir=saveDir
  rawdata.dir=rawdataDir[DOASmodel]
  reference.dir=referenceDir[DOASmodel]
  save.results=FALSE
  arima.Order=NULL
  arima.Order.FP=NULL
  special.Args = list(b=3.5
    ,Scale=function(r) median(abs(r))/0.6745
    ,delta=0
    ,filter.strength_multiplication=1.25 
    ,filter.strength_loess=0.2
    ,fam="gaussian")
  ref.spec=NULL
  ref.dark.spec=NULL
  dark.spec=NULL
  NH3.cal.spec=NULL
  SO2.cal.spec=NULL
  NO.cal.spec=NULL
  N2.cal.spec=NULL
  N2.NH3.cal.spec=N2.cal.spec
  N2.SO2.cal.spec=N2.cal.spec
  N2.NO.cal.spec=N2.cal.spec
  N2.dark.cal.spec=NULL
  force.write.daily=FALSE
  lite=TRUE
  return.AICcTab=FALSE
  return.specData=FALSE
  ncores=1
  add.name=""
}

evalOffline <- function(
	DOAS.model=NULL,
	evalperiod=NULL,
	tz.DOAS='Etc/GMT-1',
	path.length=NULL,
	tz.Output='Etc/GMT-1',
	filter.type="BmHarris",
	use.arima=FALSE,
	use.robust=TRUE,
	filter.window=NULL, 
	filter.strength=NULL,
    double.AVG=TRUE,
	fit.window=NULL,
	straylight.window=NULL,
	skip.check.daily=FALSE,
	corr.fixed.pattern=FALSE,
	fixed.pattern.length=60*24-1,
	correct.linearity=TRUE,
	correct.straylight=c("avg","linear","none"),
	correct.dark=TRUE,
	use.ref = TRUE,
	plot.results=FALSE,
	avg.period=1,
	plot.cal=FALSE,
	tau.shift=NULL,
	save.dir=NULL,
	rawdata.dir=NULL,
    reference.dir=NULL,
	save.results=FALSE,
	arima.Order=NULL,
	arima.Order.FP=NULL,
	special.Args = list(b=3.5
        ,Scale=function(r) median(abs(r))/0.6745
        ,delta=0
        ,filter.strength_multiplication=1.25 
        ,filter.strength_loess=0.2
        ,fam="gaussian"),
	ref.spec=NULL,
	ref.dark.spec=NULL,
	dark.spec=NULL,
	NH3.cal.spec=NULL,
	SO2.cal.spec=NULL,
	NO.cal.spec=NULL,
	N2.NH3.cal.spec=N2.cal.spec,
	N2.SO2.cal.spec=N2.cal.spec,
	N2.NO.cal.spec=N2.cal.spec,
	N2.dark.cal.spec=NULL,
	N2.cal.spec=NULL,
	force.write.daily=FALSE,
	lite=TRUE,
	return.AICcTab=FALSE,
	return.specData=FALSE,
    return.fitData=FALSE,
    ncores=1,
    add.name="",
    RawData = NULL,
    CalRefSpecs = NULL,
    Edinburgh_correction = TRUE,
    Serial = NULL,
    ...
  ){

  specialArgs <- list(b=3.5
    ,Scale=function(r) median(abs(r))/0.6745
    ,delta=0
    ,filter.strength_multiplication=1.25 
    ,filter.strength_loess=0.2
    ,fam="gaussian")
  specialArgs[names(special.Args)] <- special.Args
  b <- specialArgs$b
  Scale  <- specialArgs$Scale
  delta <- specialArgs$delta
  filter.strength_multiplication <- specialArgs$filter.strength_multiplication
  filter.strength_loess <- specialArgs$filter.strength_loess

	program.version <- programVersion
	if(use.arima&is.null(arima.Order))arima.Order <- rbind(expand.grid(0:2,0,0,KEEP.OUT.ATTRS = FALSE),expand.grid(0,0,1:2,KEEP.OUT.ATTRS = FALSE))
	if(use.arima&is.null(arima.Order.FP)&corr.fixed.pattern)arima.Order.FP <- rbind(expand.grid(0:2,0,0,KEEP.OUT.ATTRS = FALSE),expand.grid(0,0,1:2,KEEP.OUT.ATTRS = FALSE))
	tau.fix <- length(tau.shift)==1 

	cat("\n************\nevaluating miniDOAS model",DOAS.model,"\n")

	### prepare time range:
	### ******************************************************************************
	timerange <- prepTimeRange(evalperiod,tz.DOAS)

	### get DOAS information:
	### ******************************************************************************
	DOAS.info <- getDOASinfo(DOAS.model,timerange, Serial = Serial)

	### initialize parallelism:
	if(parl <- ncores>1){
		require(snowfall)
		if(!(wasrunning <- sfIsRunning())){
			on.exit(sfStop())
			sfInit(TRUE,ncores)
		}
		cl <- sfGetCluster()
	}

	### Read raw data:
  if (is.null(RawData)) {
  	cat("read raw data\n")
  	RawData <- readDOASdata(DOAS.info,rawdata.dir,skip.check.daily = skip.check.daily, force.write.daily = force.write.daily, ncores = ncores)
  } else {
    cat("raw data supplied - ignoring argument 'evalperiod'...\n")
    # NOTE: get time subset?
  }
	if(DOAS.model=="S1")DOAS.info <- RawData$DOASinfo

	### get DOAS windows:
	### ******************************************************************************
	DOAS.win <- getWindows(DOAS.info, filter.type, timerange, 
	straylight.window, filter.window, fit.window, 
	filter.strength, tau.shift, double.AVG)


	### read (average) reference, calibration & noise spectra
	### ******************************************************************************
	cat("read and process reference/calibration files\n")
  if (is.null(CalRefSpecs)) {
    CalRefSpecs <- getSpecSet(
      spec.dir=reference.dir,
      ref.spec=ref.spec,
      ref.dark.spec=ref.dark.spec,
      dark.spec=dark.spec,
      NH3.cal.spec=NH3.cal.spec,
      SO2.cal.spec=SO2.cal.spec,
      NO.cal.spec=NO.cal.spec,
      N2.NH3.cal.spec=N2.NH3.cal.spec,
      N2.SO2.cal.spec=N2.SO2.cal.spec,
      N2.NO.cal.spec=N2.NO.cal.spec,
      N2.dark.cal.spec=N2.dark.cal.spec,
      DOAS.model=DOAS.model
      )
  } else {
    cat("argument 'CalRefSpecs' supplied...\n")
  }


	# correct cal/ref specs:
  SpecCorr <- correctSpectra(CalRefSpecs,rawData=RawData,correct.dark=correct.dark,correct.linearity=correct.linearity,correct.straylight=correct.straylight,straylight.pix=DOAS.win$pixel3)

	### calibration concentrations in ppb
	### ******************************************************************************
	cat(sprintf("\nNH3 calibration of %1.2f mg/m3 (corresponds to %1.2f ug/m3 in a path of %1.1f m)\n",CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg,CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg * 1000 / path.length * CalRefSpecs$dat.NH3$cuvette$cuvetteLength,path.length))
	cat(sprintf("SO2 calibration of %1.2f mg/m3 (corresponds to %1.2f ug/m3 in a path of %1.1f m)\n",CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg,CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg * 1000 / path.length * CalRefSpecs$dat.SO2$cuvette$cuvetteLength,path.length))
	cat(sprintf("NO calibration of %1.2f mg/m3 (corresponds to %1.2f ug/m3 in a path of %1.1f m)\n",CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg,CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg * 1000 / path.length * CalRefSpecs$dat.NO$cuvette$cuvetteLength,path.length))

  # ### create synthetic theoretical calibration DOAS curves only based on published absorption crosssection and compare that with calibration DOAS curve
  # ### ******************************************************************************
  # f <- DOAS.info$Spectrometer$wavelength
  # sft <- -9
  # x0 <- DOAS.win$pixel1 + sft
  # f2 <- f[x0]
  # f3 <- (f[c(x0[1]-1,x0)]+f[c(x0[1],x0+1)])/2
  # if(DOAS.model!="S1"){
  #   filter.length <- 21 
  # } else {
  #   filter.length <- 41
  # }

  # ### NH3 (Cheng et al., 2006):
  # lit.crossec.NH3 <- read.table(paste0(dirname(ref.spec),"/Literature_absorption_crossections/NH3_Cheng(2006)_298K_140-230nm(0.02nm).txt"), header=FALSE, stringsAsFactors=FALSE)
  # # xb.NH3 <- which(lit.crossec.NH3[,1] >= (f3[1]-df2a) & lit.crossec.NH3[,1] <= (f3[length(f3)]+df2b))
  # f2.NH3 <- lit.crossec.NH3[,1]
  # crossec.NH3 <- lit.crossec.NH3[,2] ### in cm2 / molecule
  # # 
  # f2.NH3_filt <- f2.NH3[-c(1:((filter.length-1)/2),length(f2.NH3)-(0:((filter.length-1)/2-1)))]
  # crossec.NH3_filt <- convolve(crossec.NH3,rep(1/filter.length,filter.length),type="filter")
  # abscrossec.NH3_scaled <- crossec.NH3_filt*6.022 * 10^23/17/10^6/100^2 # cm2/molecule -> m2/ug
  # abscrossec.NH3_cuv <- 1 - abscrossec.NH3_scaled*NH3.cal*dat.NH3$cuvette[[1]] # m2/ug -> [-]
  # cs.Theor <- sapply(seq_along(f2),
  #   function(x,i){
  #     integrate(
  #       function(x){
  #         out <- approx(f2.NH3_filt,abscrossec.NH3_cuv,xout=x)$y
  #         out[is.na(out)] <- 0
  #         return(out)
  #       }
  #       ,x[i],x[i+1])$value/diff(x[i+(0:1)])
  #   },x=f3)
  # NH3.doascurve.theor <- highpass.filter(cs.Theor,filter.type, DOAS.win$filter.strength,b,Scale,delta, filter.strength_multiplication, filter.strength_loess,fam)/(NH3.cal*dat.NH3$cuvette[[1]])
  # regr.NH3 <- suppressWarnings(lmrob(NH3.doascurve.theor[DOAS.win$pixel4 + sft] ~ NH3.doascurve[DOAS.win$pixel4 + sft], setting="KS2014"))
  # theory.NH3.factor <- round(coefficients(regr.NH3)[2],3)

  # plot(f[DOAS.win$pixel1],NH3.doascurve.theor,type="n",ylab="NH3 doascurve",xlab="",panel.first={grid();abline(h=0,col="darkgrey")})
  # lines(f[DOAS.win$pixel1],NH3.doascurve,lwd=2)
  # lines(f[DOAS.win$pixel1],NH3.doascurve.theor,col="blue",lwd=2)
  # lines(f[DOAS.win$pixel1],NH3.doascurve*theory.NH3.factor,col="red",lwd=2)
  # legend("bottomright",c("NH3.doascurve","Chen-Theory",sprintf("NH3.doascurve x %1.3f",theory.NH3.factor)),lwd=2,col=c("black","blue","red"))

  # ind <- seq.int(320)
  # par(mfrow=c(2,1))
  # plot(f[DOAS.win$pixel1][ind],NH3.doascurve[ind]*theory.NH3.factor,type="l",ylab="NH3 doascurve (m2/ug)",xlab="",lwd=2,panel.first={grid();abline(h=0,col="darkgrey")})
  # legend("bottomright","NH3 doascurve (cuvette suisse 4, assuming 192 mg/m3)",lwd=2,col="black",bty="n")
  # plot(f2[ind],NH3.doascurve[ind]*theory.NH3.factor,type="l",ylab="NH3 doascurve (m2/ug)",xlab="",lwd=2,col="red",panel.first={grid();abline(h=0,col="darkgrey")})
  # lines(f2[ind],NH3.doascurve.theor[ind],col="black",lwd=2)
  # legend("bottomright",c("NH3 doascurve (cuvette suisse 4, shifted)","Lit. based NH3 doascurve (Cheng et al. 2006) convolved"),lwd=2,col=c("red","black"),bty="n")


  # ### optionally: use theoretical DOAS curves
  # ### ******************************************************************************
  # if (!use.calibration) {
  #   NH3.doascurve <- NH3.doascurve.theor
  #   #NO.doascurve <- NO.doascurve.theor ### technically possible, but no good literature NO crossection data available
  #   #SO2.doascurve <- SO2.doascurve.theor ### let out SO2 -> focus on NH3
  # }


	################################################################################
	### data processing ############################################################
	################################################################################
  if(DOAS.model=="S1" && RawData$Header[1,"st"] > strptime("20160614",format="%Y%m%d") && RawData$Header[1,"st"] < strptime("20170101",format="%Y%m%d")){
    RawData$Header[,"TECTemp"] <- sapply(strsplit(RawData$Header[,"Klima"],","),"[",4)
  }

	# averaging:
	if(avg.period!=1){
		cat("Average based on",avg.period,"minute intervals...\n")
		dts <- as.numeric(RawData$Header[,"st"] + as.numeric(RawData$Header[,"et"]-RawData$Header[,"st"],units="secs")/2 - timerange[1],units="secs")
		avgIndex <- floor(dts/(avg.period*60))
		diffAvgI <- as.logical(diff(avgIndex))
		NoAvgInt <- tapply(avgIndex,avgIndex,length)
		###
		SpecCorr$I.meas <- t(apply(SpecCorr$I.meas,1,function(x,y)tapply(x,y,mean),y=avgIndex))
		###
		n <- tapply(RawData$Header[,"AccNum"],avgIndex,sum)
		st <- RawData$Header[c(TRUE,diffAvgI),"st"]
		et <- RawData$Header[c(diffAvgI,TRUE),"et"]
		TEC.Temp <- tapply(as.numeric(RawData$Header[,"TECTemp"]),avgIndex,median)
		bT <- strsplit(RawData$Header[,"Klima"],",")
		board.Temp <- round(tapply(as.numeric(sapply(bT,"[",1)),avgIndex,mean),1)
		ambient.Temp <- round(tapply(as.numeric(sapply(bT,"[",2)),avgIndex,mean),1)
		ambient.RH <- round(tapply(as.numeric(sapply(bT,"[",3)),avgIndex,mean),1)
		time.res <- tapply(RawData$Header[,"Expos"]*RawData$Header[,"AccNum"],avgIndex,sum)
		integr.time <- tapply(RawData$Header[,"Expos"],avgIndex,function(x)paste(unique(x),collapse="/"))
		RevPos <- tapply(RawData$Header[,"RevPos"],avgIndex,function(x)paste(unique(x),collapse="/"))
		ShuPos <- tapply(RawData$Header[,"ShuPos"],avgIndex,function(x)paste(unique(x),collapse="/"))
	} else {
		NoAvgInt <- rep(1,length(RawData$Header[,"st"]))
		###
		n <- RawData$Header[,"AccNum"]
		st <- RawData$Header[,"st"]
		et <- RawData$Header[,"et"]
		TEC.Temp <- as.numeric(RawData$Header[,"TECTemp"])
		bT <- strsplit(RawData$Header[,"Klima"],",")
		board.Temp <- as.numeric(sapply(bT,"[",1))
		ambient.Temp <- as.numeric(sapply(bT,"[",2))
		ambient.RH <- as.numeric(sapply(bT,"[",3))
		time.res <- RawData$Header[,"Expos"]*RawData$Header[,"AccNum"]
		integr.time <- RawData$Header[,"Expos"]
		RevPos <- RawData$Header[,"RevPos"]
		ShuPos <- RawData$Header[,"ShuPos"]
	}
    
	# calculate differential spectra:
	DiffSpec <- diffSpecs(SpecCorr, use.ref = use.ref)

	# get calibration doascurves:
	Cal.dc <- getCalCurves(DiffSpec,DOAS.win,CalRefSpecs,...)

  ### process raw-data record-wise
  ### ******************************************************************************
  ### initiate secondary data matrices
  ### ******************************************************************************

  dim1 <- length(DOAS.win$pixel1)
  dim2 <- length(DOAS.win$pixel2)

  files <- ncol(DiffSpec$diffspec)
  cat("\n"); cat("initiate data matrices\n")
  meas.doascurve <- matrix(nrow=dim1, ncol=files)
  best.order <- delta.AICc.zero <- best.tau <- numeric(files)*NA
  aicctab <- vector(mode="list", length=files)
  se <- coeffs <- matrix(nrow=3, ncol=files)
  residual.best <- fitted.doascurve <- matrix(nrow=dim2, ncol=files)
  isna <- logical(files)
  dyn.fixed.pattern <- rep(0, dim2)

  # get correct function:
  ord <- as.character(as.numeric(NROW(arima.Order) > 1) + 1)[use.arima]
  rob <- ".rob"[use.robust]
  fitcurve <- get(paste0("fit.curves.",ifelse(use.arima,"ARIMA","OLS"),ord,rob),mode="function")

  if(any(c("loess","special") %in% DOAS.win$filter.type)){
    highpass.filter2 <- function(dat, DOAS.win) highpass.filter(dat[DOAS.win$pixel1], DOAS.win)
  } else {
    winFUN <- switch(DOAS.win$filter.type,
      "Rect" = winRect,
      "Hann" = winHann,
      "Hamming" = winHamming,
      "Blackman" = winBlackman,
      "BmNuttall" = winBmNuttall,
      "FlatTop" = winFlatTop,
      "Sin" = winSin,
      "Gauss" = winGauss,
      "Kaiser" = winKaiser,
      "BmHarris" = winBmHarris,
      "Tukey" = winTukey,
      "Poisson" = winPoisson,
      "Exp" = winExp,
      "ExpHamming" = winExpHamming,
      "DolphChebyshev" = winDolphChebyshev)
    DOAS.win$filt <- winFUN(DOAS.win$filter.strength,...)
    C_cfilter <- getFromNamespace('C_cfilter', 'stats')
    if(DOAS.win$double){
      highpass.filter2 <- function(dat,DOAS.win){
          dat - rev(.Call(C_cfilter, rev(.Call(C_cfilter, dat, DOAS.win$filt, 2L, FALSE)), DOAS.win$filt, 2L, FALSE))
      }
    } else {
      highpass.filter2 <- function(dat,DOAS.win){
          dat - .Call(C_cfilter, dat, DOAS.win$filt, 2L, FALSE)
      }
    }
  }

  if(return.fitData){

    cat("
        # highpass filtering:
        meas.doascurve[,i] <- highpass.filter2(DiffSpec$diffspec[DOAS.win$pixel1,i],DOAS.win)
        
        ### fit calibration curves to measured DOAS curve
        fitcurve(meas.doascurve[,i], DOAS.win$pixel4, dyn.fixed.pattern, Cal.dc$Xreg[DOAS.win$pixel4,], NULL, DOAS.win$tau.shift, arima.Order, path.length)
      ")

    return(
      list(
        highpass.filter2 = highpass.filter2,
        fitcurve = fitcurve,
        DiffSpec = DiffSpec,
        DOAS.win = DOAS.win,
        dyn.fixed.pattern = dyn.fixed.pattern,
        Cal.dc = Cal.dc,
        arima.Order = arima.Order,
        path.length = path.length
        )
      )
  }

	if(parl){

		cat("Parallel computing doascurve and fit...\n\n")
		pindex <- clusterSplit(cl,seq(files))

		if("special" %in% DOAS.win$filter.type){
      sfLibrary(IDPmisc)
      sfExport("b","Scale","delta","b","filter.strength_multiplication","filter.strength_loess","fam")
    }
    if(use.robust) sfLibrary(robustbase)
		sfExport("highpass.filter2","fitcurve","AICc","fitparallel","forecastArima")
    
		cat("This might take a while...\n\n")
    # p <- clusterApplyLB(cl,pindex,fitparallel,DiffSpec,DOAS.win,meas.doascurve,dyn.fixed.pattern,Cal.dc,arima.Order,path.length,
    #       isna,best.tau,delta.AICc.zero,best.order,aicctab,coeffs,se,fitted.doascurve,residual.best)
    p <- clusterApply(cl,pindex,fitparallel,DiffSpec,DOAS.win,meas.doascurve,dyn.fixed.pattern,Cal.dc,arima.Order,path.length,
          isna,best.tau,delta.AICc.zero,best.order,aicctab,coeffs,se,fitted.doascurve,residual.best)

		for(i in seq_along(p)){
		  meas.doascurve[,pindex[[i]]] <- p[[i]][[1]][,pindex[[i]]]
		  isna[pindex[[i]]] <- p[[i]][[2]][pindex[[i]]]
		  best.tau[pindex[[i]]] <- p[[i]][[3]][pindex[[i]]]
		  delta.AICc.zero[pindex[[i]]] <- p[[i]][[4]][pindex[[i]]]
		  best.order[pindex[[i]]] <- p[[i]][[5]][pindex[[i]]]
		  aicctab[pindex[[i]]] <- p[[i]][[6]][pindex[[i]]]
		  coeffs[,pindex[[i]]] <- p[[i]][[7]][,pindex[[i]]]
		  se[,pindex[[i]]] <- p[[i]][[8]][,pindex[[i]]]
		  fitted.doascurve[,pindex[[i]]] <- p[[i]][[9]][,pindex[[i]]]
		  residual.best[,pindex[[i]]] <- p[[i]][[10]][,pindex[[i]]]
		}
		# if(!corr.fixed.pattern){sfStop();on.exit()}
    if(!wasrunning){
      sfStop();on.exit()
    }
	} else {
		for(i in seq(files)){
			cat("\r",i,"/",files)

			if(median(DiffSpec$diffspec[,i],na.rm=TRUE) > -5){
				# highpass filtering:
				meas.doascurve[,i] <- highpass.filter2(DiffSpec$diffspec[DOAS.win$pixel1,i],DOAS.win)
				
				### fit calibration curves to measured DOAS curve (see Stutz & Platt, 1996, Applied Optics 30), determine best fit after shifting over given tau range, no fixed pattern considered
				### ******************************************************************************#
				fitted <- fitcurve(meas.doascurve[,i], DOAS.win$pixel4, dyn.fixed.pattern, Cal.dc$Xreg[DOAS.win$pixel4,], NULL, DOAS.win$tau.shift, arima.Order, path.length)

				if(is.null(fitted)){
					isna[i] <- TRUE
				} else {
					best.tau[i] <- fitted[[1]]
					delta.AICc.zero[i] <- fitted[[2]]
					best.order[i] <- fitted[[3]]
					aicctab[[i]] <- fitted[[4]]
					coeffs[,i] <- fitted[[5]]
					se[,i] <- fitted[[6]]
					fitted.doascurve[,i] <- fitted[[7]]
					residual.best[,i] <- fitted[[8]]
				}
			} else {
				isna[i] <- TRUE
			}
		}
	}

	# broadband:
	meas.diffspec.broadband <- DiffSpec$diffspec[DOAS.win$pixel1,] - meas.doascurve

  
	## test stats:
    nh3wts <- abs(Cal.dc$NH3.dc[DOAS.win$pixel4])/sum(abs(Cal.dc$NH3.dc[DOAS.win$pixel4]))
    so2wts <- abs(Cal.dc$SO2.dc[DOAS.win$pixel4])/sum(abs(Cal.dc$SO2.dc[DOAS.win$pixel4]))
    nowts <- abs(Cal.dc$NO.dc[DOAS.win$pixel4])/sum(abs(Cal.dc$NO.dc[DOAS.win$pixel4]))
    gof.NH3 <- round((colSums((fitted.doascurve-1)/(meas.doascurve[DOAS.win$pixel4,]-1)*nh3wts)*(colSums(meas.doascurve[DOAS.win$pixel4,]*nh3wts)-1)+1)/colSums(meas.doascurve[DOAS.win$pixel4,]*nh3wts),2)
    gof.SO2 <- round((colSums((fitted.doascurve-1)/(meas.doascurve[DOAS.win$pixel4,]-1)*so2wts)*(colSums(meas.doascurve[DOAS.win$pixel4,]*so2wts)-1)+1)/colSums(meas.doascurve[DOAS.win$pixel4,]*so2wts),2)
    gof.NO <- round((colSums((fitted.doascurve-1)/(meas.doascurve[DOAS.win$pixel4,]-1)*nowts)*(colSums(meas.doascurve[DOAS.win$pixel4,]*nowts)-1)+1)/colSums(meas.doascurve[DOAS.win$pixel4,]*nowts),2)
    gof.perc <- round((colSums((fitted.doascurve-1)/(meas.doascurve[DOAS.win$pixel4,]-1)*(nh3wts+so2wts+nowts)/3)*(colSums(meas.doascurve[DOAS.win$pixel4,]*(nh3wts+so2wts+nowts)/3)-1)+1)/colSums(meas.doascurve[DOAS.win$pixel4,]*(nh3wts+so2wts+nowts)/3),2)

    
    SAE <- apply(residual.best,2,function(x)sum(abs(x)))
    SAT <- apply(fitted.doascurve+residual.best,2,function(x)sum(abs(x)))
    SRT <- apply(fitted.doascurve+residual.best,2,function(x)sum(x))
    SAF <- apply(fitted.doascurve,2,function(x)sum(abs(x)))

    SRE <- apply(residual.best,2,function(x)sum(x))
    # STAT1 <- round(SRE/SAE,2)
    OLS <- if(use.arima)sapply(aicctab,function(x){i <- grep("0/0/0",x[,1]);out <- x[i,2];if(is.null(out))out <- NA;return(out)}) else seq_along(aicctab)*NA
    # STAT2 <- SAE/SAT
    # STAT3 <- SAE/SAF

    cat("\n")
   
    # ### running median residual patterns
    # cat("\n"); cat("creating median residual pattern\n")
    # if (files <= fixed.pattern.length) {
    #   fp.avg <- matrix(rep(apply(residual.best, 1, median, na.rm=TRUE),files), nrow=files, ncol=length(DOAS.win$pixel4), byrow=TRUE)
    # } else {      
    #   if(!(fixed.pattern.length%%2))fixed.pattern.length <- fixed.pattern.length + 1
    #   xin <- 1:ncol(residual.best)
    #   r.b <- apply(residual.best,1,function(x){y<-x[!isna];approx(c(0,xin[!isna],length(x)+1),c(y[1],y,rev(y)[1]),xout=1:length(x))$y})
    #   fp.avg <- apply(r.b, 2, runmed, k=fixed.pattern.length, endrule="constant")
    # }


    ### step 3: optionally subtract averaged fixed residual pattern from measurement data and re-evaluate (using parallel computing on OS other than Windows)
    ### ******************************************************************************
    if (corr.fixed.pattern) {
    	stop("fixed.pattern not yet available.")

      # fp.best.order <- fp.delta.AICc.zero <- fp.best.tau <- numeric(files)*NA
      # fp.aicctab <- vector(mode="list", length=files)
      # fp.se <- fp.coeffs <- matrix(nrow=3, ncol=files)
      # fp.residual.best <- fp.fitted.doascurve <- matrix(nrow=dim2, ncol=files)
    
      # fp.isna <- isna
      # if(use.arima){
      #   if(is.null(dim(arima.Order.FP))){
      #     fitcurve <- fit.curves.ARIMA1
      #   } else {
      #     fitcurve <- if(use.robust) fit.curves.ARIMA2.rob else fit.curves.ARIMA2
      #   }
      # }      
      # cat("\n"); cat("re-evaluate considering fixed residual pattern:\n")
      # if(ncores>1){
      #   sfExport("fp.avg")
      #   cat("Parallel computing fix pattern corrected fit...\nThis might take a while...\n\n")
      #   p <- clusterApply(cl,pindex,fitpWrapper,x4,fp.avg,X,fit.weights,log.meas.diffspec,tau.shift,arima.Order,meas.doascurve,filter.type,filter.strength,m.d,fp.isna,fp.best.tau,fp.delta.AICc.zero,fp.best.order,fp.aicctab,fp.coeffs,fp.se,fp.fitted.doascurve,fp.residual.best,path.length,b,Scale,delta,iter.lims,use.iter,postproc,filter.type_post,filter.strength_post,b,Scale_post,zeta, filter.strength_multiplication, filter.strength_loess,fam,FP=TRUE)
      #   for(i in seq_along(p)){
      #     fp.isna[pindex[[i]]] <- p[[i]][[2]][pindex[[i]]]
      #     fp.best.tau[pindex[[i]]] <- p[[i]][[3]][pindex[[i]]]
      #     fp.delta.AICc.zero[pindex[[i]]] <- p[[i]][[4]][pindex[[i]]]
      #     fp.best.order[pindex[[i]]] <- p[[i]][[5]][pindex[[i]]]
      #     fp.aicctab[pindex[[i]]] <- p[[i]][[6]][pindex[[i]]]
      #     fp.coeffs[,pindex[[i]]] <- p[[i]][[7]][,pindex[[i]]]
      #     fp.se[,pindex[[i]]] <- p[[i]][[8]][,pindex[[i]]]
      #     fp.fitted.doascurve[,pindex[[i]]] <- p[[i]][[9]][,pindex[[i]]]
      #     fp.residual.best[,pindex[[i]]] <- p[[i]][[10]][,pindex[[i]]]
      #   }
      #   sfStop()
      #   on.exit()
      # } else {     
      #   for(i in seq(files)){
      #     cat("\rFP corrected... ",i,"/",files)

      #     if(!isna[i]){
      #         ### fit calibration curves to measured DOAS curve (see Stutz & Platt, 1996, Applied Optics 30), determine best fit after shifting over given tau range, no fixed pattern considered
      #         ### ******************************************************************************#
      #       fitted <- fitcurve(meas.doascurve[,i], x4, fp.avg[i,], X, fit.weights, tau.shift, arima.Order.FP, path.length)

      #       if(is.null(fitted)){
      #         fp.isna[i] <- TRUE
      #       } else {
      #         fp.best.tau[i] <- fitted[[1]]
      #         fp.delta.AICc.zero[i] <- fitted[[2]]
      #         fp.best.order[i] <- fitted[[3]]
      #         fp.aicctab[[i]] <- fitted[[4]]
      #         fp.coeffs[,i] <- fitted[[5]]
      #         fp.se[,i] <- fitted[[6]]
      #         fp.fitted.doascurve[,i] <- fitted[[7]]
      #         fp.residual.best[,i] <- fitted[[8]]
      #       }
      #     }
      #   }
      # cat("\n \n")   
      # }
    }
    if (Edinburgh_correction) {
      NH3cor <- 1.16
    } else {
      NH3cor <- 1
    }

    ### write to result matrix
    ### ******************************************************************************
    if(corr.fixed.pattern){
      # # results <- as.data.frame(matrix(nrow=files, ncol=37, dimnames=list(1:files,c("date time start","date time end","millisec start","millisec end","no of acc.","fpc tau [pixel]","fpc deltaAICc(tau0-tau)","fpc order","fpc NH3 [ug/m3]","fpc NH3 SE [ug/m3]", "fpc SO2 [ug/m3]","fpc SO2 SE [ug/m3]", "fpc NO [ug/m3]","fpc NO SE [ug/m3]", "I.max","I.avg","fpc resid SSE","fp SSE","dark avg offset [counts]","TEC T [degC]","panel T [degC]","ambient T [degC]","ambient RH [perc]","shutter position","revolver position","min I/I0","threshold NH3 [ug/m3]"
      # #   ,"tau [pixel]","deltaAICc(tau0-tau)","order","NH3 [ug/m3]","NH3 SE [ug/m3]","SO2 [ug/m3]","SO2 SE [ug/m3]","NO [ug/m3]","NO SE [ug/m3]","resid SSE"
      # #   ))))#,"NH3 alternative [ug/m3]"))))
      # results <- as.data.frame(matrix(nrow=files, ncol=57, dimnames=list(1:files,c(paste0("date time start (",tz.Output,")"),paste0("date time end (",tz.Output,")"),"millisec start","millisec end","no of acc.","fpc tau [pixel]","fpc deltaAICc(tau0-tau)","fpc order","fpc NH3 [ug/m3]","fpc NH3 SE [ug/m3]", "fpc SO2 [ug/m3]","fpc SO2 SE [ug/m3]", "fpc NO [ug/m3]","fpc NO SE [ug/m3]", "I.max","I.avg","fpc resid SSE","fp SSE","dark avg offset [counts]","TEC T [degC]","panel T [degC]","ambient T [degC]","ambient RH [perc]","shutter position","revolver position","min I/I0"
      #   ,"tau [pixel]","deltaAICc(tau0-tau)","order","NH3 [ug/m3]","NH3 SE [ug/m3]","SO2 [ug/m3]","SO2 SE [ug/m3]","NO [ug/m3]","NO SE [ug/m3]","resid SSE"
      #   ,"SAE","SRE","SRE/SAE","SAF","SAT","SRT","SAE/SAF","SRE/SAF","SRE/SAT"
      #   ,"sign changes 25","sign changes 75","sign changes 125","perc fit","perc NH3","perc SO2","perc NO","dev diffspec","No. Intervals"
      #   ,"NH3-NH3_000","(NH3-NH3_000)/NH3"
      #   ))))#,"NH3 alternative [ug/m3]"))))
      # results[,1] <- format(datetime.start,"%d.%m.%Y %H:%M:%S",tz=tz.Output)
      # results[,2] <- format(datetime.end,"%d.%m.%Y %H:%M:%S",tz=tz.Output)
      # results[,3] <- as.numeric(substring(datetime.start.frac,nchar(datetime.start.frac)-2,nchar(datetime.start.frac)))
      # results[,4] <- as.numeric(substring(datetime.end.frac,nchar(datetime.end.frac)-2,nchar(datetime.end.frac)))
      # results[,5] <- n
      # results[,6] <- fp.best.tau
      # results[,7] <- fp.delta.AICc.zero
      # results[,8] <- fp.best.order
      # results[,9] <- fp.coeffs[1,]
      # results[,10] <- fp.se[1,]
      # results[,11] <- fp.coeffs[2,]  
      # results[,12] <- fp.se[2,]
      # results[,13] <- fp.coeffs[3,]
      # results[,14] <- fp.se[3,]  
      # results[,15] <- apply(I.meas, 2, max, na.rm=TRUE)
      # results[,16] <- apply(I.meas[x2,], 2, mean)
      # results[,17] <- apply(fp.residual.best^2, 2, sum)
      # results[,18] <- apply(fp.avg^2, 1, sum)
      # results[,19] <- dark.offset
      # results[,20] <- TEC.Temp
      # bT <- strsplit(board.Temp,",")
      # if(DOAS.model!="S1"&&!TwoFourTeen)results[,21] <- as.numeric(sapply(bT,"[[",1))
      # if(DOAS.model!="S1"&&!TwoFourTeen)results[,22] <- as.numeric(sapply(bT,"[[",2))
      # if(DOAS.model!="S1"&&!TwoFourTeen)results[,23] <- as.numeric(sapply(bT,"[[",3))
      # if(DOAS.model!="S1")results[,24] <- ShuPos
      # if(DOAS.model!="S1")results[,25] <- RevPos
      # results[,26] <- apply(meas.diffspec[x4,], 2, min)
      # # suppressWarnings(results[,27] <- theoretical.conc(I.attenuation.min[2], I.meas[x1[ind.max.crosssect],], path.length, 1013, 293, 17, NH3.ccd.crossec[ind.max.crosssect]))

      # results[,28] <- best.tau
      # results[,29] <- delta.AICc.zero
      # results[,30] <- best.order
      # results[,31] <- coeffs[1,]
      # results[,32] <- se[1,]
      # results[,33] <- coeffs[2,]  
      # results[,34] <- se[2,]
      # results[,35] <- coeffs[3,]
      # results[,36] <- se[3,]  
      # results[,37] <- apply(residual.best^2, 2, sum)

      # results[,38] <- SAE
      # results[,39] <- SRE
      # results[,40] <- STAT1

      # results[,41] <- SAF
      # results[,42] <- SAT
      # results[,43] <- SRT
      
      # results[,44] <- round(SAE/SAF,2)
      # results[,45] <- round(SRE/SAF,2)
      # results[,46] <- round(SRE/SAT,2)

      # results[,47] <- sc1
      # results[,48] <- sc2
      # results[,49] <- sc3
      # results[,50] <- gof.perc
      # results[,51] <- gof.NH3
      # results[,52] <- gof.SO2
      # results[,53] <- gof.NO
      # results[,54] <- MAX
      # results[,55] <- NoAvgInt
      # results[,56] <- coeffs[1,] - OLS
      # results[,57] <- (coeffs[1,] - OLS)/coeffs[1,]
 

    }else if(lite) {
        # reduce results and simplify names
		results <- as.data.frame(matrix(nrow=files, ncol=13, dimnames=list(1:files,
            c('st', 'et', 'nh3', 'so2', 'no', 'nh3_se', 'so2_se', 'no_se', 'Imax', 'n', 'tau', 'shutter', 'revolver'))))
		results[,1] <- st
		results[,2] <- et
		results[,3] <- coeffs[1,] * NH3cor
		results[,4] <- coeffs[2,]  
		results[,5] <- coeffs[3,]
		results[,6] <- se[1,] * NH3cor
		results[,7] <- se[2,]
		results[,8] <- se[3,]  
		results[,9] <- apply(RawData$RawData, 2, max, na.rm=TRUE)
		results[,10] <- n
		results[,11] <- best.tau
		results[,12] <- ShuPos
		results[,13] <- RevPos
    } else {
		results <- as.data.frame(matrix(nrow=files, ncol=42, dimnames=list(1:files,c(paste0("date time start (",tz.Output,")"),paste0("date time end (",tz.Output,")"),"millisec start","millisec end","no of acc.","tau [pixel]","deltaAICc(tau0-tau)","order","NH3 [ug/m3]","NH3 SE [ug/m3]","SO2 [ug/m3]","SO2 SE [ug/m3]","NO [ug/m3]","NO SE [ug/m3]", "I.max","I.avg(fit)","I.min(fit)","resid SSE","fp SSE","dark avg offset [counts]","TEC T [degC]","panel T [degC]","ambient T [degC]","ambient RH [perc]","shutter position","revolver position","min I/I0"
		        ,"No. intervals","SAE","SRE","SAF","SAT","SRT","SRE/SAE","SAE/SAF","SAE/SAT"
		        ,"perc fit","perc NH3","perc SO2","perc NO","NH3-NH3_000","(NH3-NH3_000)/NH3"
		        ))))
		results[,1] <- format(st,"%d.%m.%Y %H:%M:%S",tz=tz.Output)
		results[,2] <- format(et,"%d.%m.%Y %H:%M:%S",tz=tz.Output)
		results[,3] <- round(as.numeric(st)-floor(as.numeric(st)),3)
		results[,4] <- round(as.numeric(et)-floor(as.numeric(et)),3)
		results[,5] <- n
		results[,6] <- best.tau
		results[,7] <- delta.AICc.zero
		results[,8] <- best.order
		results[,9] <- coeffs[1,] * NH3cor
		results[,10] <- se[1,] * NH3cor
		results[,11] <- coeffs[2,]  
		results[,12] <- se[2,]
		results[,13] <- coeffs[3,]
		results[,14] <- se[3,]  
		results[,15] <- apply(RawData$RawData, 2, max, na.rm=TRUE)
		results[,16] <- apply(RawData$RawData[DOAS.win$pixel2,], 2, mean)
		results[,17] <- apply(RawData$RawData[DOAS.win$pixel2,], 2, min)
		results[,18] <- apply(residual.best^2, 2, sum)
		# results[,19] <- apply(fp.avg^2, 1, sum)
		results[,20] <- apply(RawData$RawData[DOAS.win$pixel3,], 2, mean)
		results[,21] <- TEC.Temp
		results[,22] <- board.Temp
		results[,23] <- ambient.Temp
		results[,24] <- ambient.RH
		results[,25] <- ShuPos
		results[,26] <- RevPos
		results[,27] <- exp(apply(DiffSpec$diffspec[DOAS.win$pixel2,], 2, min))

		results[,28] <- NoAvgInt
		results[,29] <- SAE
		results[,30] <- SRE
		results[,31] <- SAF
		results[,32] <- SAT
		results[,33] <- SRT

		results[,34] <- round(SRE/SAE,2)
		results[,35] <- round(SAE/SAF,2)
		results[,36] <- round(SAE/SAT,2)

		results[,37] <- gof.perc
		results[,38] <- gof.NH3
		results[,39] <- gof.SO2
		results[,40] <- gof.NO
		
		results[,41] <- coeffs[1,] - OLS
		results[,42] <- (coeffs[1,] - OLS)/coeffs[1,]

        colnames(results) <- paste(colnames(results),DOAS.model)
    } 


    if(lite){
      return.AICcTab <- FALSE
      return.specData <- FALSE
    } 

  InputList <- list(
    program.version=program.version,
    DOAS.model=DOAS.model,
    evalperiod=evalperiod,
    tz.DOAS=tz.DOAS,
    path.length=path.length,
    tz.Output=tz.Output,
    filter.type=filter.type,
    use.arima=use.arima,
    use.robust=use.robust,
    filter.window=filter.window, 
    filter.strength=filter.strength,
    fit.window=fit.window,
    straylight.window=straylight.window,
    skip.check.daily=skip.check.daily,
    corr.fixed.pattern=corr.fixed.pattern,
    fixed.pattern.length=fixed.pattern.length,
    correct.linearity=correct.linearity,
    correct.straylight=correct.straylight,
    correct.dark=correct.dark,
    use.ref = use.ref,
    plot.results=plot.results,
    avg.period=avg.period,
    plot.cal=plot.cal,
    tau.shift=tau.shift,
    save.dir=save.dir,
    rawdata.dir=rawdata.dir,
    reference.dir=reference.dir,
    save.results=save.results,
    arima.Order=arima.Order,
    arima.Order.FP=arima.Order.FP,
    special.Args=special.Args,
    ref.spec=ref.spec,
    ref.dark.spec=ref.dark.spec,
    dark.spec=dark.spec,
    NH3.cal.spec=NH3.cal.spec,
    SO2.cal.spec=SO2.cal.spec,
    NO.cal.spec=NO.cal.spec,
    N2.NH3.cal.spec=N2.NH3.cal.spec,
    N2.SO2.cal.spec=N2.SO2.cal.spec,
    N2.NO.cal.spec=N2.NO.cal.spec,
    N2.dark.cal.spec=N2.dark.cal.spec,
    N2.cal.spec=N2.cal.spec,
    force.write.daily=force.write.daily,
    lite=lite,
    return.AICcTab=return.AICcTab,
    return.specData=return.specData,
    ncores=ncores,
    add.name=add.name,
    ...=...
  )


    ### save results
    ### ******************************************************************************
    if(any(c(save.results,plot.cal,!isFALSE(plot.results)))){
      ### create result folder
      ### ******************************************************************************
      cat("create result directory\n")
      now <- format(now1 <- Sys.time(),"%Y%m%d%H%M")
      meas.date <- paste(format(timerange,"%y%m%d%H%M"),collapse="-")
      if(add.name!=""){add.name <- paste0("_",add.name)}
      folder <- paste(c("miniDOAS_",DOAS.model,"_", meas.date,"_Eval",now,add.name),collapse="")
      dir.create(paste(c(save.dir,"/",folder),collapse=""),recursive=TRUE)
      ### create and save log file
      ### ******************************************************************************
      cat("\n"); cat("create and save log file\n")
      cat("to do: logfile!\n")
      time.res2 <- paste(unique(round(time.res/1000, 0)),collapse="/")
      logfile <- c(
        paste0("miniDOAS evaluation program version: ",program.version),
        "*****************************************************************************",
        paste0("calculation performed on: ",format(now1,format="%d.%m.%Y %H:%M")),
        paste0("raw data file path: ",rawdata.dir),
        paste0("evaluated ",ncol(RawData$RawData)," raw data files"),
        paste0("results saved at: ",save.dir,"/",folder),
        "",
        paste0("reference files directory: ",reference.dir),
        if(is.list(ref.spec)){
          paste0("\"in-calculation\" lamp reference spectrum chosen between: ",paste(CalRefSpecs$dat.ref$timerange,collapse=" - ")," (",tz(CalRefSpecs$dat.ref$timerange),")\nfrom the raw data located at: ",ref.spec$dir)
        } else {
          paste0("lamp reference spectrum recorded on: ",paste(CalRefSpecs$dat.ref$timerange,collapse=" - "))
        },
        # if(!any(is.na(ref.dark.spec)))paste(c("\"in-calculation\" dark reference spectrum chosen between: ",paste(ref.dark.spec[[2]],collapse=" - "),"\nfrom the data located at: ",ref.dark.spec[[1]]),collapse=""),
        paste0("NH3 calibration conditions: NH3 conc = ",round(CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg,2)," mg/m3 in a cuvette of ",CalRefSpecs$dat.NH3$cuvette$cuvetteLength," m length"),
        paste0("SO2 calibration conditions: SO2 conc = ",round(CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg,2)," mg/m3 in a cuvette of ",CalRefSpecs$dat.SO2$cuvette$cuvetteLength," m length"),
        paste0("NO calibration conditions: NO conc = ",round(CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg,2)," mg/m3 in a cuvette of ",CalRefSpecs$dat.NO$cuvette$cuvetteLength," m length"),
        "",
        sprintf("total light path = %1.2f m",path.length),
        paste0("spectrometer: ",DOAS.info$Spectrometer$"Spectrometer Name"),
        paste0(round(mean(n),0)," spectra were accumulated on average over ",time.res2," seconds for one data file"),    
        "",
        paste0("data evaluated from ",format(st[1],format="%d.%m.%Y %H:%M:%S")," to ",format(et[length(et)],format="%d.%m.%Y %H:%M:%S")," (",tz(st),")"),
        paste0("spectra were filtered over ",DOAS.win$filter.window[1]," to ",DOAS.win$filter.window[2]," nm"),
        paste0("... using ",if(DOAS.win$filter.type %in% "special"){
            paste0("the 'special' filter by Sintermann et al. (2016) with the following control settings:\n",
              paste(names(DOAS.win$special.control),gsub("\n","",DOAS.win$special.control),sep=": ",collapse=", ")
            )
          } else if(DOAS.win$filter.type %in% "loess"){
            paste0("a loess filter with span: ",DOAS.win$filter.strength," and family: ",DOAS.win$special.control$fam)
          } else {
            paste0("a ",if(DOAS.win$double) "double " else "",sub("Rect","rectangular",DOAS.win$filter.type)," moving average")
          }),
        paste0("spectra were ",ifelse(!tau.fix,"","not "),"allowed to be shifted",
          ifelse(!tau.fix,
            if(all(diff(tau.shift)==1)){
              paste0(" between ",paste(range(tau.shift),collapse=" to ")," pixel")
            }else{
              paste0(" at c(",paste(tau.shift,collapse=","),") pixel")
            },
            ""
          )
          ),
        if(tau.fix)paste0("pixel shift was fixed at ",tau.shift),
        if(corr.fixed.pattern)paste0("a dynamically pre-evaluated fixed residual pattern calculated by runmed over ",fixed.pattern.length," datapoints was subtracted"),
        if(use.arima){
          paste0("the multiple linear fit was performed using",if(use.robust) " robust " else " ","ARIMA with",ifelse(is.null(dim(arima.Order)),"out "," "),"AICc model averaging",if(is.null(dim(arima.Order))){c(", given ARIMA order ",paste(arima.Order,collapse="/"))})
        }else{
          paste0("the multiple linear fit was performed using",if(use.robust) " robust " else " ","OLS")
        },
        paste0("the fit was performed between ",DOAS.win$fit.window[1]," and ",DOAS.win$fit.window[2]," nm"),
        if(corr.fixed.pattern&&use.arima)paste(c("the multiple linear fit of the fix pattern corrected data was performed using ARIMA with",ifelse(is.null(dim(arima.Order.FP)),"out "," "),"AICc model averaging",if(is.null(dim(arima.Order.FP))){c(", given ARIMA order ",paste(arima.Order.FP,collapse="/"))}),collapse=""),
        paste0("(potential) dark-count/straylight logging from ",DOAS.win$straylight.window[1]," to ",DOAS.win$straylight.window[2]," nm"),
        # paste(c( "NH3 concentrations were calculated based on ",ifelse(use.calibration, "calibrated", "published"), " absorption-crosssections"),collapse=""),
        # if(DOAS.model!="S1"){ 
        #   paste(c( "... the respective NH3 DOAS-curves differ roughly by a factor of ", theory.NH3.factor, " (calibration vs. published)\n"),
        #   c( "theoretically, the spectrum I will be attenuated by ",round(theoretical.conc(I.attenuation.min[1], 1, path.length, 1000, 293, 17, NH3.ccd.crossec[ind.max.crosssect]),0)," ug/m3 NH3 down to ", I.attenuation.min[1]*100, "% I at wavelength of ",round(f[x1[ind.max.crosssect]],2)," nm (= location of maximum absorption cross-section)"),collapse="")
        # },
        if(use.arima&&!is.null(dim(arima.Order))){
          paste("*****************************************************************************\nSupplied Orders:\n",
          paste(apply(arima.Order,1,paste,collapse="/"),collapse="\n"),
          "\n*****************************************************************************")
        }

      )

      write(logfile, file=paste0(save.dir,"/",folder,"/",folder,"_log.txt"))    
    } else {
      folder <- NULL
    }

    if(save.results){
      cat("\n"); cat("write results\n")
      write.table(results, file=paste(c(save.dir,"/",folder,"/",folder,"_avg",avg.period,".csv"),collapse=""), sep=";", na="#N/A", dec=".", row.names=FALSE, col.names=TRUE)
    
    }


  ### plot calibration spectra and DOAS curves
  ### ******************************************************************************
  if (plot.cal) {
    cat("\n"); cat("plot calibration spectra and DOAS curves")
    calref.cols <- c("blue", "green", "orange", "cyan")                                                                                                                                                                                        
    f <- DOAS.info$Spectrometer$wavelength
    x1 <- DOAS.win$pixel1
    x2 <- DOAS.win$pixel2
    x4 <- DOAS.win$pixel4
    NH3.cal.ug <- CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg * 1000 * CalRefSpecs$dat.NH3$cuvette$cuvetteLength
    SO2.cal.ug <- CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg * 1000 * CalRefSpecs$dat.SO2$cuvette$cuvetteLength
    NO.cal.ug <- CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg * 1000 * CalRefSpecs$dat.NO$cuvette$cuvetteLength
    NH3.doascurve <- Cal.dc$NH3.dc
    SO2.doascurve <- Cal.dc$SO2.dc
    NO.doascurve <- Cal.dc$NO.dc
    pdf(file=paste(c(save.dir,"/",folder,"/",folder,"_cal.pdf"),collapse=""), width=0.8*7, height=1.25*7)
    par(mfcol=c(3,1))
    plot.calibration.spectra(calref.cols, f[x1], CalRefSpecs$dat.ref$dat.spec[x1], CalRefSpecs$dat.dark$dat.spec[x1], CalRefSpecs$dat.N2.NH3$dat.spec[x1], CalRefSpecs$dat.N2.dark$dat.spec[x1], CalRefSpecs$dat.NH3$dat.spec[x1], CalRefSpecs$dat.SO2$dat.spec[x1], CalRefSpecs$dat.NO$dat.spec[x1])
    plot.calibration.diffspecs(calref.cols, f[x1], exp(DiffSpec$NH3.diffspec[x1]), exp(DiffSpec$SO2.diffspec[x1]), exp(DiffSpec$NO.diffspec[x1]), NH3.doascurve * NH3.cal.ug, SO2.doascurve * SO2.cal.ug, NO.doascurve * NO.cal.ug)
    plot.calibration.DOAScurves(calref.cols, f[x1], f[x2], NH3.doascurve[x4], SO2.doascurve[x4], NO.doascurve[x4])
    # if(DOAS.model!="S1")lines(f[x2], NH3.doascurve.theor[x4], col="gray")
    #lines(f[x2], SO2.doascurve.theor[x4], col="gray")
    graphics.off()
  }   

  ### plot evaluation steps, and optionally merge individual plots into one pdf animation (via LaTex)
  ### ******************************************************************************
  if (!isFALSE(plot.results)) {
      cat("\nplot results evaluation details\n")
      plot.name <- paste0(folder,"_avg",avg.period)
      gwd <- getwd()
      on.exit(setwd(gwd))
      setwd(paste0(save.dir,"/",folder))
      if(isTRUE(plot.results) || !(tolower(plot.results) %in% "html")){
        if(isTRUE(plot.results)) plot.results <- "pdf"
        plotfu <- get(plot.results,mode="function")
        if(plot.results %in% "pdf") pdf(file=paste0(plot.name,".pdf"),onefile=TRUE)
          for (i in 1:nrow(results)) {

            cat(paste("\r saving plots... ",format(x=round(i/nrow(results)*100,1),nsmall=1),"%",sep=""))

            if (corr.fixed.pattern) {
              avg.fixed.pattern <- fp.avg[i,]
            } else {
              avg.fixed.pattern <- rep(0, dim2)
            }

            best.tau <- results[i,6]

            m.diffspec <- DiffSpec$diffspec[DOAS.win$pixel1,i]
            m.diffspec.broadband <- meas.diffspec.broadband[,i]
            m.meas.doascurve <- meas.doascurve[,i]
            f <- DOAS.info$Spectrometer$wavelength
            x1 <- DOAS.win$pixel1
            x2 <- DOAS.win$pixel2
            x4 <- DOAS.win$pixel4

            if (corr.fixed.pattern) {
              m.fitted.doascurve.best <- fp.fitted.doascurve[,i]
              m.residual.best <- fp.residual.best[,i]
            } else {
              m.fitted.doascurve.best <- fitted.doascurve[,i]
              m.residual.best <- residual.best[,i]  
            }
           
            if(!all(is.na(m.meas.doascurve))){
              if(!(plot.results %in% "pdf")) plotfu(file=paste0(plot.name,"_",i,".",plot.results), width=400, height=600)
              par(mfcol=c(3,1),mar=c(2, 4, 2, 1)+0.5,oma=c(2,0,3,0))
              cex.annotations <- 1.25
              
              ### differential spectrum and filtered differential broadband spectrum
              plot(f[x1], x1, type="n", ylim=range(m.diffspec, na.rm=TRUE), xlab="", ylab="[relative units]", main="differential spectrum", cex.axis=cex.annotations, cex.lab=cex.annotations)
              lines(f[x1], m.diffspec, lty=1, lwd=1.5, col="black") 
              lines(f[x1], m.diffspec.broadband, lty=1, lwd=1.5, col="green") 
              legend("bottomright", legend=c("differential","differential broadband"), bty="n", col=c("black","green"), lty=c(1,2), cex=0.75)
              
              ### measued & fitted DOAS curve
              plot(f[x1], x1, type="n", ylim=range(c(m.fitted.doascurve.best,m.meas.doascurve[x4],m.meas.doascurve[x4] - avg.fixed.pattern,m.meas.doascurve[x4+best.tau] - avg.fixed.pattern), na.rm=TRUE), xlab="", ylab="[relative units]", main=paste("DOAS curves &",ifelse(use.robust,"robust",""),ifelse(use.arima,"ARIMA fit","OLS fit")), cex.axis=cex.annotations, cex.lab=cex.annotations)
              abline(h=0,lty=2,col="gray60")
              lines(f[x1 - best.tau], m.meas.doascurve, lty=1.5, lwd=1, col="darkgreen")
              lines(f[x2], (m.meas.doascurve[x4] - avg.fixed.pattern), lty=1, lwd=1, col="gray80")
              lines(f[x2], (m.meas.doascurve[x4+best.tau] - avg.fixed.pattern), lty=1, lwd=1.5, col="black")
              lines(f[x2], m.fitted.doascurve.best, lty=1, lwd=1.5, col="blue")
              lgnd <- if(corr.fixed.pattern) {c("measured","measured-FP+shift","measured-FP","fitted")} else{c("measured","measured+shift","measured","fitted")}
              lgnd.lty <- c(1,1,1,1)
              lgnd.col <- c("darkgreen","black","gray80","blue")
              legend("bottomright", legend=lgnd, bty="n", col=lgnd.col, lty=lgnd.lty, cex=0.75)
              
              ### residual spectrum
              plot(f[x1], x1, type="n", ylim=range(c(m.fitted.doascurve.best,m.meas.doascurve[x4],m.meas.doascurve[x4] - avg.fixed.pattern,m.meas.doascurve[x4+best.tau] - avg.fixed.pattern), na.rm=TRUE), xlab="", ylab="[relative units]", main=sprintf("residual spectrum: SAE = %1.2E, SRE = %1.2E",sum(abs(m.residual.best)),sum(m.residual.best)), cex.axis=cex.annotations, cex.lab=cex.annotations)
              abline(h=0,lty=2,col="gray60")
              if (corr.fixed.pattern) {lines(f[x2 - best.tau], avg.fixed.pattern, lty=1, lwd=1.5, col="gold")}
              lines(f[x2], m.fitted.doascurve.best, lty=2, lwd=1.5, col="blue")
              # lines(f[x2], m.residual, lty=2, lwd=1, col="gray80")
              lines(f[x2], m.residual.best, lty=1, lwd=1.5, col="black")
              
              ### legend & title
              lgnd <- if(corr.fixed.pattern) {c("fixed pattern","fitted DOAS curve","residuals")} else{c("fitted DOAS curve","residuals")}
              lgnd.lty <- if(corr.fixed.pattern) {c(2,2,1)} else{c(2,1)}
              lgnd.col <- if(corr.fixed.pattern) {c("gold","blue","black")} else{c("blue","black")}
              legend("bottomright", legend=lgnd, bty="n", col=lgnd.col, lty=lgnd.lty, cex=0.75)
              msg1 <- paste(results[i,1]," to ",results[i,2],", n=",results[i,5],", tau=",best.tau," px, Imax=",round(results[i,15],0)," counts", sep="")
              msg2 <- paste("NH3=",round(results[i,9],1)," +- ",round(results[i,10],2)," ug/m3, SO2=",round(results[i,11],1)," +- ",round(results[i,12],2)," ug/m3, NO=",round(results[i,13],1)," +- ",round(results[i,14],2)," ug/m3", sep="")
              # mtext(as.expression(substitute(italic(msg), list(msg=msg))), line=-1.25, outer=TRUE, cex=0.5)
              # mtext(substitute(italic(msg), list(msg=msg1)), line=1, outer=TRUE, cex=0.75)
              # mtext(substitute(italic(msg), list(msg=msg2)), line=-0.5, outer=TRUE, cex=0.75)
              mtext(msg1, line=1, outer=TRUE, cex=0.75)
              mtext(msg2, line=-0.5, outer=TRUE, cex=0.75)
              mtext("wavelength [nm]",side=1,outer=TRUE,cex=0.75,adj=0.53)

              if(!(plot.results %in% "pdf")) dev.off()
            }
        }
      } else {
        images_dir <- paste0(save.dir,"/",folder,"/images")
        dir.create(images_dir)
        require(animation)
        saveHTML(
          {for (i in 1:nrow(results)) {

            cat(paste("\r saving plots... ",format(x=round(i/nrow(results)*100,1),nsmall=1),"%",sep=""))

            if (corr.fixed.pattern) {
              avg.fixed.pattern <- fp.avg[i,]
            } else {
              avg.fixed.pattern <- rep(0, dim2)
            }

            best.tau <- results[i,6]

            m.diffspec <- DiffSpec$diffspec[DOAS.win$pixel1,i]
            m.diffspec.broadband <- meas.diffspec.broadband[,i]
            m.meas.doascurve <- meas.doascurve[,i]
            f <- DOAS.info$Spectrometer$wavelength
            x1 <- DOAS.win$pixel1
            x2 <- DOAS.win$pixel2
            x4 <- DOAS.win$pixel4

            if (corr.fixed.pattern) {
              m.fitted.doascurve.best <- fp.fitted.doascurve[,i]
              m.residual.best <- fp.residual.best[,i]
            } else {
              m.fitted.doascurve.best <- fitted.doascurve[,i]
              m.residual.best <- residual.best[,i]  
            }
           
            if(!all(is.na(m.meas.doascurve))){
              png(file=paste0(images_dir,"/",plot.name,"_",i,".png"), width=400, height=600)
              par(mfcol=c(3,1),mar=c(2, 4, 2, 1)+0.5,oma=c(2,0,3,0))
              cex.annotations <- 1.25
              
              ### differential spectrum and filtered differential broadband spectrum
              plot(f[x1], x1, type="n", ylim=range(m.diffspec, na.rm=TRUE), xlab="", ylab="[relative units]", main="differential spectrum", cex.axis=cex.annotations, cex.lab=cex.annotations)
              lines(f[x1], m.diffspec, lty=1, lwd=1.5, col="black") 
              lines(f[x1], m.diffspec.broadband, lty=1, lwd=1.5, col="green") 
              legend("bottomright", legend=c("differential","differential broadband"), bty="n", col=c("black","green"), lty=c(1,2), cex=0.75)
              
              ### measued & fitted DOAS curve
              plot(f[x1], x1, type="n", ylim=range(c(m.fitted.doascurve.best,m.meas.doascurve[x4],m.meas.doascurve[x4] - avg.fixed.pattern,m.meas.doascurve[x4+best.tau] - avg.fixed.pattern), na.rm=TRUE), xlab="", ylab="[relative units]", main=paste("DOAS curves &",ifelse(use.robust,"robust",""),ifelse(use.arima,"ARIMA fit","OLS fit")), cex.axis=cex.annotations, cex.lab=cex.annotations)
              abline(h=0,lty=2,col="gray60")
              lines(f[x1 - best.tau], m.meas.doascurve, lty=1.5, lwd=1, col="darkgreen")
              lines(f[x2], (m.meas.doascurve[x4] - avg.fixed.pattern), lty=1, lwd=1, col="gray80")
              lines(f[x2], (m.meas.doascurve[x4+best.tau] - avg.fixed.pattern), lty=1, lwd=1.5, col="black")
              lines(f[x2], m.fitted.doascurve.best, lty=1, lwd=1.5, col="blue")
              lgnd <- if(corr.fixed.pattern) {c("measured","measured-FP+shift","measured-FP","fitted")} else{c("measured","measured+shift","measured","fitted")}
              lgnd.lty <- c(1,1,1,1)
              lgnd.col <- c("darkgreen","black","gray80","blue")
              legend("bottomright", legend=lgnd, bty="n", col=lgnd.col, lty=lgnd.lty, cex=0.75)
              
              ### residual spectrum
              plot(f[x1], x1, type="n", ylim=range(c(m.fitted.doascurve.best,m.meas.doascurve[x4],m.meas.doascurve[x4] - avg.fixed.pattern,m.meas.doascurve[x4+best.tau] - avg.fixed.pattern), na.rm=TRUE), xlab="", ylab="[relative units]", main=sprintf("residual spectrum: SAE = %1.2E, SRE = %1.2E",sum(abs(m.residual.best)),sum(m.residual.best)), cex.axis=cex.annotations, cex.lab=cex.annotations)
              abline(h=0,lty=2,col="gray60")
              if (corr.fixed.pattern) {lines(f[x2 - best.tau], avg.fixed.pattern, lty=1, lwd=1.5, col="gold")}
              lines(f[x2], m.fitted.doascurve.best, lty=2, lwd=1.5, col="blue")
              # lines(f[x2], m.residual, lty=2, lwd=1, col="gray80")
              lines(f[x2], m.residual.best, lty=1, lwd=1.5, col="black")
              
              ### legend & title
              lgnd <- if(corr.fixed.pattern) {c("fixed pattern","fitted DOAS curve","residuals")} else{c("fitted DOAS curve","residuals")}
              lgnd.lty <- if(corr.fixed.pattern) {c(2,2,1)} else{c(2,1)}
              lgnd.col <- if(corr.fixed.pattern) {c("gold","blue","black")} else{c("blue","black")}
              legend("bottomright", legend=lgnd, bty="n", col=lgnd.col, lty=lgnd.lty, cex=0.75)
              msg1 <- paste(results[i,1]," to ",results[i,2],", n=",results[i,5],", tau=",best.tau," px, Imax=",round(results[i,15],0)," counts", sep="")
              msg2 <- paste("NH3=",round(results[i,9],1)," +- ",round(results[i,10],2)," ug/m3, SO2=",round(results[i,11],1)," +- ",round(results[i,12],2)," ug/m3, NO=",round(results[i,13],1)," +- ",round(results[i,14],2)," ug/m3", sep="")
              # mtext(as.expression(substitute(italic(msg), list(msg=msg))), line=-1.25, outer=TRUE, cex=0.5)
              # mtext(substitute(italic(msg), list(msg=msg1)), line=1, outer=TRUE, cex=0.75)
              # mtext(substitute(italic(msg), list(msg=msg2)), line=-0.5, outer=TRUE, cex=0.75)
              mtext(msg1, line=1, outer=TRUE, cex=0.75)
              mtext(msg2, line=-0.5, outer=TRUE, cex=0.75)
              mtext("wavelength [nm]",side=1,outer=TRUE,cex=0.75,adj=0.53)

              dev.off()
            }
          }
          cat("\n")}
          ,htmlfile=paste0(plot.name,".html"), autobrowse=FALSE, autoplay=FALSE, use.dev=FALSE, ani.type='png', width=400, height=600, interval=0.3, verbose=FALSE, img.name=paste0("Rplot_",plot.name,"_"), title=paste0(plot.name,"_")
        )
      }

      graphics.off()
      setwd(gwd)
      on.exit()
      cat("\nfinished plotting")
  }


  cat("\nfinished evaluation of",DOAS.model,"\n")
  if(!return.AICcTab)aicctab <- NULL
  if(return.specData){
    vars <- list(f=f,x1=x1,x2=x2,x3=x3,x4=x4,dat.ref=dat.ref,dat.ref.dark=dat.ref.dark,dat.dark=dat.dark,dat.N2.NH3=dat.N2.NH3,dat.N2.SO2=dat.N2.SO2,dat.N2.NO=dat.N2.NO,dat.N2.dark=dat.N2.dark,dat.NH3=dat.NH3,dat.SO2=dat.SO2,dat.NO=dat.NO,cal.diffspecs=cal.diffspecs,NH3.doascurve=NH3.doascurve,SO2.doascurve=SO2.doascurve,NO.doascurve=NO.doascurve,NH3.doascurve.theor=NH3.doascurve.theor)
  } else {
    vars <- NULL
  }

  if(lite){
    CalRefSpecs <- NULL 
    RawData <- NULL
  }

  results <- structure(results
      ,AICc.Tables=aicctab
      ,vars=vars
      ,CalRefSpecs=CalRefSpecs
      ,RawData=RawData
      ,inputlist=InputList
      ,callEval=match.call()
      ,class=c("DOASeval","data.frame")
    )
  return(results)
}



# ev. noch plot=T/F einbauen
extractResults <- function(resList,excols,as.xts=TRUE,tz="UTC",Ithresh=NULL,DOASmodel=NULL,tzResults=NULL){ 
  # resList=Results;excols=9:12
  # resList=test
  # as.xts=TRUE;Ithresh=NULL;tz="CET"
  # excols= c("NH3 [ug/m3]","SO2 sigma [ug/m3]","NO sigma [ug/m3]")
  if(as.xts)require(xts)
  if(cD <- class(resList)[1]=="DOASeval")resList <- list(resList)
  if(is.null(DOASmodel)){
    Doasmodels <- sapply(lapply(resList,"attr","inputlist"),"[[","DOAS.model")
  } else {
    Doasmodels <- DOASmodel
  }
  if(is.null(tzResults)){
    tzIn <- sapply(lapply(resList,"attr","inputlist"),"[[","tz.Output")
  } else {
    tzIn <- tzResults
  }
  if(cD)names(resList) <- Doasmodels
  lrL <- length(resList)
  takecols <- excols
  excols <- colnams <- colnams0 <- vector(lrL,mode="list")
  for(i in 1:lrL){
    colnams0[[i]] <- gsub(paste0(" ",Doasmodels[1]),"",colnames(resList[[1]]))
    colnams[[i]] <- gsub("fpcorr ","",colnams0[[i]])
    if(class(takecols)=="character"){
      excols2 <- which(colnams[[i]] %in% takecols)
      isna <- which(!(takecols %in% colnams[[i]][excols2]))
      if(length(isna))stop("Could not find colnames: ",paste(takecols[isna],collapse=", "))
      excols[[i]] <- excols2
    } else {
      excols[[i]] <- takecols
    } 
  }

  index <- vector(lrL,mode="list")
  if(!is.null(Ithresh)){
    Ithresh <- rep(Ithresh,lrL)[seq(lrL)]
    for(i in 1:lrL){
      index[[i]] <- resList[[i]][,15]<Ithresh[i]
    }
  } else {
    for(i in 1:lrL){
      index[[i]] <- logical(length(resList[[i]][,15]))
    }    
  }
  run_names <- names(resList)
  if(is.null(run_names))run_names <- seq(lrL)
  abbr_rn <- abbreviate(run_names)
  out <- vector(lrL,mode="list")
  names(out) <- abbr_rn
  cn <- paste(unlist(mapply(function(x,y)rep(x,length(y)),x=abbr_rn,y=excols,SIMPLIFY=FALSE)),unlist(mapply(function(x,y)x[y],x=colnams0,y=excols,SIMPLIFY=FALSE)),sep=" - ")
  if(as.xts){
    for(i in 1:lrL){
      st <- as.POSIXct(strptime(resList[[i]][,1],"%d.%m.%Y %H:%M:%S",tz=tzIn[i]))
      et <- as.POSIXct(strptime(resList[[i]][,2],"%d.%m.%Y %H:%M:%S",tz=tzIn[i]))
      Time <- st + (et-st)/2
      if(lrL>1){
        Time <- as.POSIXct(round(Time,"mins"))
      }
      out[[i]] <- xts(resList[[i]][,excols[[i]],drop=FALSE],Time,tz=tz)
      out[[i]][index[[i]],] <- NA
    }
    out <- do.call(merge,out)
  } else {
    for(i in 1:lrL){
      st <- as.POSIXct(strptime(resList[[i]][,1],"%d.%m.%Y %H:%M:%S",tz=tzIn[i]))
      et <- as.POSIXct(strptime(resList[[i]][,2],"%d.%m.%Y %H:%M:%S",tz=tzIn[i]))
      Time <- st + (et-st)/2
      if(lrL>1){
        Time <- format(round(Time,"mins"),"%Y%m%d%H%M")
      }
      if(i==1){
        out <- data.frame(Time=Time,resList[[i]][,excols[[i]],drop=FALSE])
        out[index[[i]],-1] <- NA
      } else {
        out <- merge(out,data.frame(Time=Time,resList[[i]][,excols[[i]],drop=FALSE]),by="Time",all=TRUE)
        out[index[[i]],-1] <- NA
      }
    }
    out <- out[,-1]
  }
  colnames(out) <- cn
  out <- structure(out,class=c("DOASres",class(out)))
  return(out)
}

runAppDOAS <- function(ShinyInput){
  require(shiny)
  require(dygraphs)
  require(xts)
  require(lubridate)
  if(is.ibts(ShinyInput)){
		ShinyInput <- as.xts(ShinyInput)
  }
  dat <- cbind(ShinyInput,"none"=NA)
  r10 <- if(all(is.na(dat[,1]))) c(-1,1) else range(dat[,1],na.rm=TRUE)
  dr1 <- diff(r10)
  r1_ini <- round(r10 + c(-1,1)*0.1*ifelse(dr1<1,1,dr1),1)
  r20 <- if(all(is.na(dat[,2]))) c(-1,1) else range(dat[,2],na.rm=TRUE)
  dr2 <- diff(r20)
  r2_ini <- round(r20 + c(-1,1)*0.1*ifelse(dr2<1,1,dr2),1)
  runApp(
    list(
    ui={
      fluidPage(
        
        # Application title
        titlePanel("miniDOAS time series"),
        
        sidebarLayout(position="right",
                      
                      sidebarPanel(width=2,
                        helpText("display variables"),
                        selectInput("var1", 
                                    label = "var1:",
                                    choices = names(dat),
                                    selected = names(dat)[1]
                        ),
                        selectInput("yax1", 
                                    label = "axis1:",
                                    choices = c("left","right"),
                                    selected = "left"
                        ),                         
                        selectInput("var2", 
                                    label = "var2:",
                                    choices = names(dat),
                                    selected = names(dat)[2]
                        ),
                        selectInput("yax2", 
                                    label = "axis2:",
                                    choices = c("left","right"),
                                    selected = "right"
                        ),
                         
                        selectInput("var3", 
                                    label = "var3:",
                                    choices = names(dat),
                                    selected = "none"
                        ),
                        selectInput("yax3", 
                                    label = "axis3:",
                                    choices = c("left","right"),
                                    selected = "right"
                        ),
                        sliderInput("sl1",
                                    label = "Range left:",
                                    value= r1_ini,
                                    min = r1_ini[1],
                                    max = r1_ini[2],
                                    sep="",
                                    pre="",
                                    post="",
                                    step=0.1,
                                    ticks=TRUE,
                                    round=-1
                          ),
                        sliderInput("sl2",
                                    label = "Range right:",
                                    value= r2_ini,
                                    min = r2_ini[1],
                                    max = r2_ini[2],
                                    sep="",
                                    pre="",
                                    post="",
                                    step=0.1,
                                    round=-1
                          )
                        ,checkboxInput("checkzero", label="include zero", value = FALSE)
                        ,downloadButton('downloadPlot', 'Save Plot to png')
                           
                      ),     
                      
                      mainPanel(dygraphOutput("timeseries"))
                      
        )
            
      )
    },
    server={
      
        function(input, output, session) {


              observe({
                ydat <- cbind(dat[,input$var1],dat[,input$var2],dat[,input$var3])
                Times <- as.POSIXct(ydat)
                if(!is.null(input$timeseries_date_window)){
                  start <- strptime(deparse(input$timeseries_date_window[[1]]),"\"%Y-%m-%dT%H:%M:%OSZ\"",tz="UTC")
                  end <- strptime(deparse(input$timeseries_date_window[[2]]),"\"%Y-%m-%dT%H:%M:%OSZ\"",tz="UTC")
                  start <- Times[which.min(abs(start-Times))]
                  end <- Times[which.min(abs(end-Times))]
                  ydat <- window(ydat,start=start,end=end) 
                }
                yax1 <- which(c(input$yax1,input$yax2,input$yax3)=="left")
                yt1 <- c(input$var1,input$var2,input$var3)[yax1]
                yax1 <- yax1[!(yt1 %in% "none")]
                yax2 <- which(c(input$yax1,input$yax2,input$yax3)=="right")
                yt2 <- c(input$var1,input$var2,input$var3)[yax2]
                yax2 <- yax2[!(yt2 %in% "none")]
                drawy1 <- as.logical(length(yax1))
                drawy2 <- as.logical(length(yax2))
                if(drawy1){
                  r10 <- if(all(is.na(ydat[,yax1]))) c(-1,1) else range(ydat[,yax1],na.rm=TRUE)
                  dr1 <- diff(r10)
                  r1 <- round(range(ydat[,yax1],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr1<1,1,dr1),1)
                          if(input$checkzero&r1[1]>0)r1[1]<-0
                          if(input$checkzero&r1[2]<0)r1[2]<-0
                }else{
                  r1 <- c(-1,1) + c(-1,1)*0.04
                }
                if(drawy2){
                  r20 <- if(all(is.na(ydat[,yax2]))) c(-1,1) else range(ydat[,yax2],na.rm=TRUE)
                  dr2 <- diff(r20)
                  r2 <- round(range(ydat[,yax2],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr2<1,1,dr2),1)
                          if(input$checkzero&r2[1]>0)r2[1]<-0
                          if(input$checkzero&r2[2]<0)r2[2]<-0
                }else{
                  r2 <- c(-1,1) + c(-1,1)*0.04
                }
                # print(r1)
                # print(r2)
                updateSliderInput(session, "sl1", label = "Range left:", step=0.1,value = r1,
                  min = r1[1], max = r1[2])
                updateSliderInput(session, "sl2", label = "Range right:", step=0.1,value = r2,
                  min = r2[1], max = r2[2])
              })

             plotReactive <- reactive({
                ydat <- cbind(dat[,input$var1],dat[,input$var2],dat[,input$var3])
                Times <- as.POSIXct(ydat)
                if(!is.null(input$timeseries_date_window)){
                  start <- strptime(deparse(input$timeseries_date_window[[1]]),"\"%Y-%m-%dT%H:%M:%OSZ\"",tz="UTC")
                  end <- strptime(deparse(input$timeseries_date_window[[2]]),"\"%Y-%m-%dT%H:%M:%OSZ\"",tz="UTC")
                  start <- Times[which.min(abs(start-Times))]
                  end <- Times[which.min(abs(end-Times))]
                  ydat <- window(ydat,start=start,end=end) 
                  px <- pretty(xlim1<-c(start,end))
                } else {
                  px <- pretty(xlim1<-c(min(Times),max(Times)))
                }
                yax1 <- which(c(input$yax1,input$yax2,input$yax3)=="left")
                yt1 <- c(input$var1,input$var2,input$var3)[yax1]
                yax1 <- yax1[!(yt1 %in% "none")]
                yax2 <- which(c(input$yax1,input$yax2,input$yax3)=="right")
                yt2 <- c(input$var1,input$var2,input$var3)[yax2]
                yax2 <- yax2[!(yt2 %in% "none")]
                drawy1 <- as.logical(length(yax1))
                drawy2 <- as.logical(length(yax2))
                ploty <- FALSE
                par(oma=c(0,0,0,2),cex.lab=1.25,cex.axis=1.25)
                if(drawy1){
                  ylim1 <- ylim1()
                  if(input$checkzero&ylim1[1]>0)ylim1[1]<-0
                  if(input$checkzero&ylim1[2]<0)ylim1[2]<-0
                  plot(ydat[,1],ann=FALSE,axes=FALSE,auto.grid=FALSE,type="n",main="",ylim=ylim1,xlim=xlim1)
                  axis(1,at=px,labels=strftime(px,"%d.%m %H:%M"))
                  grid(nx=NA,ny=NULL)
                  abline(v=px,lty="dotted",col="lightgrey")
                  for(i in yax1)lines(ydat[,i],col=c("yellowgreen","dodgerblue","orange")[i])
                  axis(2,at=pretty(ylim1),labels=pretty(ylim1))
                  mtext(ylab1(),side=2,line=2.5,cex=1.25)
                  ploty <- TRUE
                }
                if(drawy2){
                  ylim2 <- ylim2()
                  if(input$checkzero&ylim2[1]>0)ylim2[1]<-0
                  if(input$checkzero&ylim2[2]<0)ylim2[2]<-0
                  if(ploty){
                    par(new=TRUE,oma=c(0,0,0,2))
                    plot(ydat[,1],ann=FALSE,axes=FALSE,auto.grid=FALSE,type="n",main="",ylim=ylim2,xlim=xlim1)
                  } else {
                    plot(ydat[,1],ann=FALSE,axes=FALSE,auto.grid=FALSE,type="n",main="",ylim=ylim2,xlim=xlim1)
                    axis(1,at=px,labels=strftime(px,"%d.%m %H:%M"))
                    grid(nx=NA,ny=NULL)
                    abline(v=px,lty="dotted",col="lightgrey")
                  }
                  axis(4,at=pretty(ylim2),labels=pretty(ylim2))
                  mtext(ylab(),side=4,line=2,cex=1.25)
                  for(i in yax2)lines(ydat[,i],col=c("yellowgreen","dodgerblue","orange")[i])     
                }
                legend("top",colnames(ydat),lty=1,col=c("yellowgreen","dodgerblue","orange"),horiz=TRUE,bty="n",cex=1,inset=c(0,-0.075),xpd=T)
              })


               # dat2 <- reactive({
               #  cbind(dat[,input$var1],dat[,input$var2],dat[,input$var3])
               #  })

               ylab1fu <- reactive({
                yt1 <- c(input$var1,input$var2,input$var3)[which(c(input$yax1,input$yax2,input$yax3)=="left")]
                paste(yt1[!(yt1 %in% "none")],collapse=" / ")
              })

               ylabfu <- reactive({
                yt2 <- c(input$var1,input$var2,input$var3)[which(c(input$yax1,input$yax2,input$yax3)=="right")]
                paste(yt2[!(yt2 %in% "none")],collapse=" / ")
              }
                )


               ylab1 <- renderText(ylab1fu())
               ylab <- renderText(ylabfu())


               ylim1 <- reactive({
                yax1 <- which(c(input$yax1,input$yax2,input$yax3)=="left")
                yt1 <- c(input$var1,input$var2,input$var3)[yax1]
                yax1 <- yax1[!(yt1 %in% "none")]
                drawy1 <- as.logical(length(yax1))
                if(drawy1){
                  input$sl1
                } else {
                  NULL
                }
               })
               ylim2 <- reactive({
                yax1 <- which(c(input$yax1,input$yax2,input$yax3)=="right")
                yt1 <- c(input$var1,input$var2,input$var3)[yax1]
                yax1 <- yax1[!(yt1 %in% "none")]
                drawy1 <- as.logical(length(yax1))
                if(drawy1){
                  input$sl2
                } else {
                  NULL
                }
               })

          output$timeseries <- renderDygraph({
            ydat <- cbind(dat[,input$var1],dat[,input$var2],dat[,input$var3])
            # Times <- as.POSIXct(ydat)
            # if(!is.null(input$timeseries_date_window)){
            #   start <- strptime(deparse(input$timeseries_date_window[[1]]),"\"%Y-%m-%dT%H:%M:%OSZ\"",tz="UTC")
            #   end <- strptime(deparse(input$timeseries_date_window[[2]]),"\"%Y-%m-%dT%H:%M:%OSZ\"",tz="UTC")
            #   start <- Times[which.min(abs(start-Times))]
            #   end <- Times[which.min(abs(end-Times))]
            #   ydat <- window(ydat,start=start,end=end) 
            # }
            # yax1 <- which(c(input$yax1,input$yax2,input$yax3)=="left")
            # yt1 <- c(input$var1,input$var2,input$var3)[yax1]
            # yax1 <- yax1[!(yt1 %in% "none")]
            # yax2 <- which(c(input$yax1,input$yax2,input$yax3)=="right")
            # yt2 <- c(input$var1,input$var2,input$var3)[yax2]
            # yax2 <- yax2[!(yt2 %in% "none")]
            # drawy1 <- as.logical(length(yax1))
            # drawy2 <- as.logical(length(yax2))
            # if(drawy1){
            #   r10 <- if(all(is.na(ydat[,yax1]))) c(-1,1) else range(ydat[,yax1],na.rm=TRUE)
            #   dr1 <- diff(r10)
            #   r1 <- round(range(ydat[,yax1],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr1<1,1,dr1),1)
            #           if(input$checkzero&r1[1]>0)r1[1]<-0
            #           if(input$checkzero&r1[2]<0)r1[2]<-0
            # }else{
            #   r1 <- c(-1,1) + c(-1,1)*0.04
            # }
            # if(drawy2){
            #   r20 <- if(all(is.na(ydat[,yax2]))) c(-1,1) else range(ydat[,yax2],na.rm=TRUE)
            #   dr2 <- diff(r20)
            #   r2 <- round(range(ydat[,yax2],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr2<1,1,dr2),1)
            #           if(input$checkzero&r2[1]>0)r2[1]<-0
            #           if(input$checkzero&r2[2]<0)r2[2]<-0
            # }else{
            #   r2 <- c(-1,1) + c(-1,1)*0.04
            # }
          
            dygraph(ydat, main=" ")  %>%
            dySeries(name=input$var1, label=input$var1, axis=ifelse(input$yax1=="left","y","y2")) %>%
            dySeries(name=input$var2, label=input$var2, axis=ifelse(input$yax2=="left","y","y2")) %>%
            dySeries(name=input$var3, label=input$var3, axis=ifelse(input$yax3=="left","y","y2")) %>%
            dyAxis(name="y", label=ylab1(), valueRange = ylim1())  %>% 
            dyAxis(name="y2", label=ylab(), valueRange = ylim2())  %>% 
            # dyAxis(name="y", label=ylab1(), valueRange = r1)  %>% 
            # dyAxis(name="y2", label=ylab(), valueRange = r2)  %>% 
            dyOptions(colors=c("yellowgreen","dodgerblue","orange"), includeZero=input$checkzero, axisLineWidth=2, axisLabelFontSize=1.1*14, axisTickSize=2*3, maxNumberWidth=3, logscale=FALSE, stepPlot=FALSE, useDataTimezone=TRUE) %>%
            dyRoller(rollPeriod=1) %>%
            dyHighlight(highlightCircleSize=3, highlightSeriesBackgroundAlpha=0.3, hideOnMouseOut=TRUE) %>%
            dyRangeSelector(height=30, strokeColor="", retainDateWindow=TRUE) %>%
            dyLegend(width=2*250, hideOnMouseOut=TRUE)
            
          })

          output$downloadPlot <- downloadHandler(
            filename=function(){paste0("Shinysave_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png")},
            content = function(file){
              on.exit(dev.off())
              png(file,width=2*480,height=480)
              plotReactive()      
            }
            )
          
        }
    }
    )
  )
}


