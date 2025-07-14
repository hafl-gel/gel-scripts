

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; prepare time range:
prepTimeRange <- function(tr,tzDOAS="") {
    # check time zone (rather unnecessary, though)
    if (grepl("^GMT[+]",tzDOAS)) {
        warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("GMT+","Etc/GMT-",tzDOAS,fixed=TRUE))    
    } else if (grepl("^UTC[+]",tzDOAS)) {
        warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("UTC+","Etc/GMT-",tzDOAS,fixed=TRUE))     
    } else if (grepl("^GMT[-]",tzDOAS)) {
        warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("GMT-","Etc/GMT+",tzDOAS,fixed=TRUE))     
    } else if (grepl("^UTC[-]",tzDOAS)) {
        warning(tzDOAS," is not an official time zone! Using corresponding time zone ",tzDOAS <- sub("UTC-","Etc/GMT+",tzDOAS,fixed=TRUE))       
    }
    # prepare and return times
    if (inherits(tr,"character")) {
        evalperiod <- rep(tr,2)[1:2]
        timerange <- parse_date_time(evalperiod, c("d.m.Y","d-m-Y","Y-m-d"), tz=tzDOAS, quiet=TRUE)
        if (!is.na(timerange[2])) timerange[2] <- timerange[2] + 24*3600
        timerange2 <- parse_date_time(evalperiod, c("d.m.Y H:M:S","d.m.Y H:M","d-m-Y H:M:S","d-m-Y H:M","Y-m-d H:M:S","Y-m-d H:M"), tz=tzDOAS, quiet=TRUE)
        timerange[is.na(timerange)] <- timerange2[is.na(timerange)]
        timerange
    } else if (inherits(tr,"POSIXlt")) {
        as.POSIXct(tr)
    } else if (!inherits(tr,"POSIXct")) {
        stop("prepTimeRange: argument \"tr\" must be either of class character or POSIXt")
    } else {
        tr
    }
}

### describe spectrometer's linearity (see Ocean Optics calibration sheets)
### ****************************************************************************** 
linearity.func <- function(x, cfs) {
    y <- cfs[1] 
    for (i in seq.int(length(cfs) - 1)) {
        y <- y + cfs[i + 1] * x ^ i  
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
    ser <- qs::qserialize(obj, preset = 'high')
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
        out <- qs::qdeserialize(readBin(con, 'raw', N2, endian = 'big'))
    }
    # close connection
    close(con)
    on.exit()
    # return result
    out
}

process_dailyfiles <- function(folders, path_data, doas_info, RawData, 
    path_dailyfiles = path_data) {
    require(qs)
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
            file.rename(rawdat_all[rawdatSize <= 4000], sub("[.]txt$", ".invalid", rawdat_all[rawdatSize <= 4000]))
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
            out <- list(intensity = lapply(dat[-seq.int(doas_info$rawdata.structure$'Header Lines'), ], as.numeric), header = header)
            # write daily file (files, data, path)
            write_daily(basename(rawdat), out, file.path(path_dailyfiles, i))
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
    lf_daily, path_data, doas_info, rd_list, path_dailyfiles = path_data) {
    # target daily files
    target_daily <- paste0(rawdat_folder, ".bin")
    # check existing daily files
    if (force.write.daily) {
        # "no daily files exist"
        daily_exist <- rep(FALSE, length(target_daily))
    } else {
        # get existing daily files
        daily_exist <- target_daily %in% lf_daily
        # check daily files
        if(!skip.check.daily){
            cat("checking daily files...\n")
            for (i in which(daily_exist)) {
                # read 'old' files
                files_daily <- try(read_daily(file.path(path_dailyfiles, target_daily[i]), files = TRUE), silent = TRUE)
                # read present files
                files_now <- list.files(file.path(path_data, rawdat_folder[i]), pattern = '[.]txt')
                # check for changes
                if (length(setdiff(files_now, files_daily)) == 0){
                    cat(paste0("..daily file ", target_daily[i], " exists and is ok...\n"))
                    # read data from daily file
                    rd_list[[rawdat_folder[i]]] <- read_daily(file.path(path_dailyfiles, target_daily[i]))
                } else {
                    # delete daily file (and re-create it later)
                    unlink(file.path(path_dailyfiles, target_daily[i]))
                    daily_exist[i] <- FALSE
                }
            }   
        } else {
            # read existing daily files
            for (i in which(daily_exist)) {
                cat("reading daily file", target_daily[i], "...\n")
                rd_list[[rawdat_folder[i]]] <- read_daily(file.path(path_dailyfiles, target_daily[i]))
            }     
        }
    }
    # process daily files
    rd_list[!daily_exist] <- process_dailyfiles(rawdat_folder[!daily_exist], path_data, 
        doas_info, rd_list, path_dailyfiles)
    # return
    rd_list[!is.null(rd_list)]
}

readDOASdata <- function(DOASinfo, dataDir, rawdataOnly = FALSE, skip.check.daily = FALSE, 
    force.write.daily = FALSE, timerange = DOASinfo$timerange, timezone = "", ncores = 1,
    path_dailyfiles = dataDir) {

    if(parl <- (is.list(ncores) || ncores > 1)) {
        require(parallel)
        if (cl_was_up <- is.list(ncores)) {     
            cl <- ncores
        } else {
            cl <- parallel::makePSOCKcluster(rep('localhost', ncores))
            on.exit(parallel::stopCluster(cl))
        }
        parallel::clusterCall(cl, library, package = 'lubridate', character.only = TRUE)
        parallel::clusterCall(cl, library, package = 'qs', character.only = TRUE)
        # get progress bar if blsmodelr is available
        if (requireNamespace('bLSmodelR', quietly = TRUE)) {
            cluster_apply_lb <- bLSmodelR:::.clusterApplyLB
        } else {
            cluster_apply_lb <- parallel::clusterApplyLB
        }
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
    lf.raw.dir <- list.files(dataDir, pattern = '^[0-9]{8}$')
    # get existing daily files
    if (!dir.exists(path_dailyfiles)) {
        cat('directory "', path_dailyfiles, '" does not exist - creating it...\n', sep = '')
        dir.create(path_dailyfiles, recursive = TRUE)
    }
    lf_daily <- list.files(path_dailyfiles, pattern = '^[0-9]{8}[.]bin$')

    if (!length(lf.raw.dir)) {
        if (rawdataOnly || force.write.daily) {
            if(!dir.exists(dirname(dataDir))){
                stop("Folder \"",dirname(dataDir),"\" does not exist!")
            } else if(!dir.exists(dataDir)){
                stop("Folder \"",dataDir,"\" does not exist!")
            } else {
                stop("No raw data available for specified timerange!\n",
                "Folder \"",dataDir,"\" does not contain any doas data!")
            }
        } else if (!any(rawdatfolder %in% sub('.bin', '', lf_daily, fixed = TRUE))) {
            stop("No data available for specified timerange!")
        }
    }

    # prepare raw data output
    RawData <- vector(length(rawdatfolder), mode = "list")
    names(RawData) <- rawdatfolder

    # if rawdataOnly else daily files
    if (rawdataOnly) {

        if (parl) {
            # # export functions
            # parallel::clusterExport(cl, list(''))
            # loop over folders
            cat("reading raw data files in parallel...\n")
            raw_data <- cluster_apply_lb(cl, rawdatfolder, read_rawdata, data_dir = dataDir, doas_info = DOASinfo)
            names(raw_data) <- rawdatfolder
        } else {
            # loop over folders
            raw_data <- lapply(rawdatfolder, read_rawdata, data_dir = dataDir, doas_info = DOASinfo)
            names(raw_data) <- rawdatfolder
        }
        
    } else if (parl) {

        # export functions
        parallel::clusterExport(cl, list('process_dailyfiles', 'write_daily', 'read_daily'))

        # call check daily in parallel
        cat("checking daily files in parallel...\n")
        raw_data <- cluster_apply_lb(cl, seq_along(RawData), function(i, check_fu, rdf, rd, ...) {
            unlist(
                check_fu(rdf[[i]], rd_list = rd[[i]], ...),
                recursive = FALSE)
                            }, 
            check_fu = check_dailyfiles, rdf = rawdatfolder, rd = RawData, skip.check.daily, force.write.daily,
            lf_daily, dataDir, DOASinfo, path_dailyfiles)

        # close cluster connection
        if (!cl_was_up) {
            parallel::stopCluster(cl)
            on.exit()
        }

    } else {

        raw_data <- check_dailyfiles(rawdatfolder, skip.check.daily, force.write.daily, 
            lf_daily, dataDir, DOASinfo, RawData, path_dailyfiles)

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

    # fix time zone
    header$st <- lubridate::force_tz(header$st, tz(timerange))
    header$et <- lubridate::force_tz(header$et, tz(timerange))

    # get eval period subset
    takeme <- (header$et >= timerange[1]) & (header$st <= timerange[2])

    # check timerange:
    if (!any(takeme)) stop("No data available for specified time range!\n")

    # subset header
    header <- header[takeme, ]

    # write header
    Header <- data.frame(
        DOASmodel = DOASinfo$DOASmodel,
        Spectrometer = DOASinfo$Spectrometer$"Spectrometer Name",
        header[, 1:8],
        AccTime = round(as.numeric(difftime(header$et, header$st, units = "secs")), 3),
        st = header$st,
        et = header$et,
        row.names = seq_along(header$st),
        stringsAsFactors = FALSE
    )

    # select data in timerange, keep list
    RawData <- unlist(lapply(raw_data, function(x) x[[1]][-1]), recursive = FALSE)[takeme]

    # devide total counts by # of accumulations
    RawData <- mapply('/', RawData, Header[, 'AccNum'], SIMPLIFY = FALSE)

    return(list(RawData=RawData,Header=Header,DOASinfo=DOASinfo))
}

#### raw data processing function
read_rawdata <- function(dir, data_dir = dataDir, doas_info = DOASinfo) {
    cat(paste0("..reading raw data from folder ",dir,"\n"))
    # list files
    rawdat_all <- list.files(file.path(data_dir, dir), pattern = ".txt", full.names = TRUE)
    # read files
    if (length(rawdat_all)) {
        timerange <- doas_info$timerange
        if (doas_info$DOASmodel == "S1" && 
            timerange[1] > strptime("20160101","%Y%m%d") && 
            timerange[1] < strptime("20170101", "%Y%m%d")) {
            rawdat.et <- parse_date_time(sub("^.*/(.*)[.]txt$", "\\1", rawdat_all),
                c("%Y%m%d %H%M%OS", "%y%m%d%H%M%S"), tz = tz(timerange))
        } else {
            rawdat.et <- parse_date_time(sub("^.*_(.*)[.]txt$", "\\1", rawdat_all), 
                c("%Y%m%d %H%M%OS", "%y%m%d%H%M%S"), tz = tz(timerange))
        }
        rawdat.st <- rawdat.et - median(diff(rawdat.et))
        takeme <- (rawdat.et >= timerange[1]) & (rawdat.st <= timerange[2])
        rawdat_all <- rawdat_all[takeme]
        rawdatSize <- file.size(rawdat_all)
        rawdat <- rawdat_all[rawdatSize > 4000]
        if (any(rawdatSize <= 4000)) {
            cat(paste0(sum(rawdatSize <= 4000)," invalid files (< 4000 bytes)!\n"))
            file.rename(rawdat_all[rawdatSize <= 4000], 
                sub("[.]txt$", ".invalid", rawdat_all[rawdatSize <= 4000]))
        }
        if (length(rawdat)) {
            if (doas_info$DOASmodel == "S1" && (as.numeric(dir) - 20160101) >= 0 && (as.numeric(dir) - 20170101) < 0) {
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
                # read files
                firstRow <- c(seq.int(doas_info$rawdata.structure$"Header Lines"),doas_info$Spectrometer$wavelength)
                fls <- lapply(rawdat, function(x) read.table(x, header=FALSE, nrows=doas_info$Spectrometer$"Pixel Number"+doas_info$rawdata.structure$"Header Lines",sep="\n",as.is=TRUE,quote="",comment.char=""))
                # temporary solution to missing Klima entries:
                if(nrow(fls[[1]]) == (length(firstRow) - 1)){
                    fls <- lapply(fls, function(x) rbind(x[1:5,, drop = FALSE], 
                            "Board- / Ambient-T (degC) / ambient-RH (perc)= NA, NA, NA,", 
                            x[6:1051,, drop = FALSE]))
                } else if(any(check_fls <- sapply(fls,nrow) != length(firstRow))){
                    stop("check files:\n",paste(rawdat[check_fls],collapse="\n"))
                }
                dat <- data.frame(
                    wavelength=c(seq.int(doas_info$rawdata.structure$"Header Lines"),doas_info$Spectrometer$wavelength),
                    fls)
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
            # write raw data as list
            out <- list(intensity = lapply(dat[-seq.int(doas_info$rawdata.structure$'Header Lines'), ], as.numeric), header = header)
        } else {
            out <- NULL
        } 
    } else {
        out <- NULL
        cat("No data in folder", dir, "\n")
    }
    out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; average Spectra:
# TODO:
avgSpec <- function(rawdat,type=c("raw","cal","ref","dark"),tracer=c("ambient","NH3","NO","SO2"),saveToFile=FALSE
    ,pathLength=NA,NH3ambient_ug=NA,NOambient_ug=NA,SO2ambient_ug=NA,cuvetteLength=ifelse(tracer[1]=="ambient",NA,0.075),cuvetteConc_mg=NA,Dirname=getwd(),verbose=TRUE) {

    # Calc Sum/n:
    # Specs <- apply(rawdat[[1]][-(1:13),],1,as.numeric)/rawdat[[2]][,"AccNum"]
    Specs <- rawdat[[1]]
    # average specs:
    SpecAvg <- rowMeans(as.data.frame(Specs),na.rm=TRUE)

    type <- type[1]
    tracer <- tracer[1]

    cuvetteGas <- switch(type,
        cal = {
            if (tracer=="ambient")stop("ambient & cal specified!")
            if (is.na(cuvetteConc_mg))stop(
                "Calibration concentration must be specified!\n",
                "Cuvette revolver suisse 1-6:\n",
                "   NH3: 193.4095 mg/m3 (not corrected with 1.16 Edinburgh correction)\n",
                "   NO:  593.9938 mg/m3\n",
                "   SO2: 76.29128 mg/m3\n"
            )
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
        Imax=list(range=range(lapply(Specs,max)),mean=max(SpecAvg))
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
        ,paste0("number of accumulations: ",paste(Info$AccNum$range,collapse=" to ")," (avg: ",sprintf("%1.1f",Info$AccNum$mean),")")
        ,paste0("exposure time (ms): ",paste(Info$Expos$range,collapse=" to ")," (avg: ",sprintf("%1.1f",Info$Expos$mean),")")
        ,paste0("TEC temperature (deg C): ",paste(Info$TEC$range,collapse=" to ")," (avg: ",sprintf("%1.1f",Info$TEC$mean),")")
        ,paste0("ambient temperature (deg C): ",paste(Info$Ambient$range,collapse=" to ")," (avg: ",sprintf("%1.1f",Info$Ambient$mean),")")
        ,paste0("Imax: ",sprintf("%1.1f to %1.1f",Info$Imax$range[1],Info$Imax$range[2])," (avg: ",sprintf("%1.1f",Info$Imax$mean),")")
        ,"-----")

    if (saveToFile) {
        filename <- paste0("miniDOAS_",rawdat$DOASinfo$DOASmodel,"_",paste(tracer,type,"spec",sep="_"),"_",paste(format(c(rawdat[[2]][1,"st"],rev(rawdat[[2]][,"et"])[1]),"%y%m%d%H%M%S"),collapse="-"),".txt")
        write.table(txt,paste(Dirname,filename,sep="/"),quote=FALSE,sep="\n",row.names=FALSE,col.names=FALSE)
        write.table(data.frame(rawdat$DOASinfo$Spectrometer$wavelength,SpecAvg),paste(Dirname,filename,sep="/"),quote=FALSE,sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
        cat("Averaged spectra saved to",paste(Dirname,filename,sep="/"),"\n")
    }

    txt[1] <- paste0("Averaging ",insertme,":\n") 
    txt <- paste(c(txt,"\n"),collapse="\n")

    if (verbose) {
        cat(txt)
    }

    return(list(SpecAvg=SpecAvg,Info=Info,Specs=Specs,DOASinfo=rawdat$DOASinfo,txt=txt))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; get dark/filter/fit windows and filter strength:
getWindows <- function(DOASinfo, filter.type = NULL, timerange = Sys.time(), 
    straylight.window = NULL, filter.window = NULL, fit.window = NULL, filter.strength = NULL, 
    tau.shift = NULL) {

    # get doas info
    if (inherits(DOASinfo, "character")) DOASinfo <- getDOASinfo(DOASinfo, timerange)

    # check filter
    if (is.null(filter.type)) {
        filter.type <- getOption('md.filter.type')
    } else if (!(filter.type %in% names(getOption('md.filter.function.list')))) {
        stop("filter.type = \"", filter.type, "\" is not a valid filter type! -> see names(getOption('md.filter.function.list'))")
    }
    
    # check filter window
    if (is.null(filter.window)) filter.window <- getOption('md.filter.window')

    # check filter strength
    if (is.null(filter.strength)) {
        filter.strength <- getOption('md.filter.strength')
    } else if ((filter.strength %% 2) < 1) {
        filter.strength <- round(filter.strength)
        if ((filter.strength %% 2) < 1) {
            filter.strength <- filter.strength + 1
        }
        warning("filter.strength must be an uneven number and has been rounded to the nearest uneven number.")
    }  

    # check fit.window
    if (is.null(fit.window)) fit.window <- getOption('md.fit.window')

    # check straylight window
    if (is.null(straylight.window)) straylight.window <- c(150, 160)

    # check tau shift
    if (is.null(tau.shift) || is.na(tau.shift)) tau.shift <- 0

    # filter pixel
    pixel_filter <- which(DOASinfo$Spectrometer$wavelength >= filter.window[1] & DOASinfo$Spectrometer$wavelength <= filter.window[2])

    # fitting pixel
    pixel_fit <- max(
        which(DOASinfo$Spectrometer$wavelength >= fit.window[1])[1]
        ,pixel_filter[1]-tau.shift
        ):min(
        rev(which(DOASinfo$Spectrometer$wavelength <= fit.window[2]))[1]
        ,rev(pixel_filter)[1]-tau.shift
    )

    # straylight pixel
    pixel_straylight <- which(DOASinfo$Spectrometer$wavelength >= straylight.window[1] & DOASinfo$Spectrometer$wavelength <= straylight.window[2])

    # check fit window lower bound
    if ((pixel_fit[1] - pixel_filter[1]) < filter.strength) {
        pixel_fitLo <- pixel_filter[1] + filter.strength
        fit.window[1] <- DOASinfo$Spectrometer$wavelength[pixel_fitLo]
    } else {
        pixel_fitLo <- pixel_fit[1]
    }

    # check fit window upper bound
    if ((pixel_filter[length(pixel_filter)] - pixel_fit[length(pixel_fit)]) < filter.strength) {
        pixel_fitHi <- pixel_filter[length(pixel_filter)] - filter.strength
        fit.window[2] <- DOASinfo$Spectrometer$wavelength[pixel_fitHi]
    } else {
        pixel_fitHi <- pixel_fit[length(pixel_fit)]
    }
    
    # rebuild pixel_fit with updated boundaries
    pixel_fit <- seq(pixel_fitLo, pixel_fitHi)

    # find pixel_fit within pixel_filter
    pixel_dc <- match(pixel_fit, pixel_filter)

    list(
        filter.window = filter.window,
        fit.window = fit.window,
        straylight.window = straylight.window,
        filter.type = filter.type,
        filter.strength = filter.strength,
        tau.shift = tau.shift,
        pixel_filter = pixel_filter,
        pixel_fit = pixel_fit,
        pixel_straylight = pixel_straylight,
        pixel_dc = pixel_dc
    )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~ function; get calibration and reference spectra:
# TODO:
getSpec <- function(spec, DOASmodel = NULL, lite = FALSE, SpecName = NULL) {

    if (is.null(SpecName)) {
        if (is.list(spec) && !is.null(spec$spec.name)) {
            SpecName <- spec$spec.name
        } else {
            SpecName <- as.character(substitute(spec))
        }
    }

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
        "N2.dark.cal.spec"="N2",
        'unknown spec.name'
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
        "N2.dark.cal.spec"="dark",
        'unknown spec.name'
    ) 


    if (is.list(spec)) {
        if (all(names(spec) %in% c('SpecAvg', 'Info', 'Specs', 'DOASinfo', 'txt'))) {
            return(spec)
        } else if (all(names(spec)%in%c("dat.spec","ambient","cuvette","timerange"))) {
            dat.spec <- spec$dat.spec
            ambient <- spec$ambient
            cuvette <- spec$cuvette
            timerange <- spec$timerange
            lite <- TRUE    
        } else {

            if (type != 'unknown spec.name') {
                spec$type <- type
            } else if (is.null(spec$type)) {
                stop('list entry "type" is needed if spec.name is not recognized')
            }

            if (type=="cal") {
                if (is.null(spec$cuvetteLength)) spec$cuvetteLength <- 0.075
                if ((!is.null(spec$cuvetteConc_mg) && spec$cuvetteConc_mg <= 0) ||
                    spec$cuvetteLength <= 0) stop('Please provide valid info on cuvette length and concentration')
            }

            if (tracer != 'unknown spec.name') {
                spec$tracer <- tracer
            } else if (is.null(spec$tracer)) {
                stop('list entry "tracer" is needed if spec.name is not recognized')
            }

            if (is.null(spec$rawdat)) {
                if (is.null(spec$tz)) spec$tz <- "Etc/GMT-1"

                if (is.character(spec$dir)) {
                    spec$rawdat <- try(readDOASdata(DOASmodel,timerange=spec$timerange,dataDir=spec$dir,rawdataOnly=TRUE,timezone=spec$tz))
                    if (inherits(spec$rawdat,"try-error")) {
                        stop(paste0(SpecName,": ",unlist(strsplit(spec$rawdat," : "))[2]))
                    }
                } else {
                    stop('list entry "dir" should contain the path to the raw data if "rawdat" is not provided itself')
                }
            }

            AvgSpec <- do.call(avgSpec,spec[names(spec) %in% names(formals(avgSpec))])

            if (lite) {
                dat.spec <- AvgSpec$SpecAvg
                ambient <- AvgSpec$Info[3:6]
                cuvette <- AvgSpec$Info[7:8]
                timerange <- c(min(spec$rawdat$Header[,"st"]),max(spec$rawdat$Header[,"et"])) 
                info.spec <- unlist(strsplit(AvgSpec$txt, split = '\n'))
                info.spec <- info.spec[!(info.spec %in% '')]
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

        if (any(ind <- grepl(SpecName_old,dir.spec,fixed=TRUE))) {
            info.spec <- read.table(dir.spec[which(ind)], sep="\t", stringsAsFactors=FALSE, fill=TRUE, colClasses="character",nrows=10)[,1]
            dat.spec <- read.table(dir.spec[which(ind)], sep="\t", stringsAsFactors=FALSE, fill=TRUE, skip=10,colClasses="numeric")[,1]
            ambient <- list(pathLength=NA,NH3ambient_ug=NA,SO2ambient_ug=NA,NOambient_ug=NA)
            if (any(SpecName==c("NH3.cal.spec","SO2.cal.spec","NO.cal.spec"))) {
                konz.spec <- as.numeric(unlist(strsplit(gsub("^.*): ","",readLines(dir.spec[which(ind)],n=3)[3]),",")))
                cuvette <- list(cuvetteLength=0.075,cuvetteConc_mg=konz.spec[3]*konz.spec[1]/(8.3144598*konz.spec[2])*switch(tracer,NH3=17,SO2=64,NO=30)/10,cuvetteGas=tracer)
            } else {
                cuvette <- list(cuvetteLength=ifelse(tracer=="ambient",NA,0.075),cuvetteConc_mg=NA,cuvetteGas=ifelse(type=="ref",ifelse(tracer=="ambient","none","N2"),NA))
            }
            # get timerange from file names
            timerange <- strptime(unlist(strsplit(sub('miniDOAS.*_([0-9]{12}[-][0-9]{12})[-]Eval.*', '\\1', basename(dir.spec[which(ind)])), split = '[-]')), '%Y%m%d%H%M', tz = 'Etc/GMT-1')
        } else {
            ind <- grepl(SpecName_new,dir.spec,fixed=TRUE)
            if (all(!ind)) {
                if (type %in% 'dark') {
                    ind[grep('dark', dir.spec)[1]] <- TRUE
                }
                if (all(!ind)) {
                    stop('Calibration file for ', tracer, '/', type, ' is missing in provided path:\n',
                    '    "', spec, '"')
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
            # get timerange from file names
            tms <- unlist(strsplit(info.spec[18], split = ' to '))
            timerange <- c(
                fast_strptime(tms[1], 'averaging period: %Y-%m-%d %H:%M:%S', tz = 'Etc/GMT-1', lt = FALSE),
                fast_strptime(tms[2], '%Y-%m-%d %H:%M:%S', tz = 'Etc/GMT-1', lt = FALSE)
                )
        }
        lite <- TRUE
    }

    if (DOASmodel=="S1"&&lite&&(is.na(timerange[1]) || timerange[1] < parse_date_time("20170101","Ymd",tz=tz(timerange[1])))) {
        dat.spec <- dat.spec[seq.int(1021)]
    }

    if (lite) {
        return(
            list(
                dat.spec = dat.spec,
                ambient = ambient,
                cuvette = cuvette,
                timerange = timerange,
                info.spec = info.spec
                )
        )  
    } else {
        return(AvgSpec)
    }
}


### read set of (average) reference, calibration & noise spectra
### ******************************************************************************
# TODO:
getSpecSet <- function(
    spec.dir="",  
    ref.spec=NULL,
    dark.spec=NULL,
    ref.dark.spec=dark.spec,
    NH3.cal.spec=NULL,
    N2.cal.spec=NULL,
    SO2.cal.spec=NULL,
    NO.cal.spec=NULL,
    N2.NH3.cal.spec=N2.cal.spec,
    N2.SO2.cal.spec=N2.cal.spec,
    N2.NO.cal.spec=N2.cal.spec,
    N2.dark.cal.spec=dark.spec,
    DOAS.model=NULL) {
    # check if serialized calref file is provided
    if (is.character(spec.dir) && isFALSE(file.info(spec.dir)$isdir)) {
        if (file.exists(spec.dir)) {
            spec.dir <- switch(
                sub('.*[.]([a-z]+)$', '\\1', spec.dir)
                , qs = qs::qread(spec.dir)
                , rds = readRDS(spec.dir)
                , stop('unknown calref file format: *.',
                    sub('.*[.]([a-z]+)$', '\\1', spec.dir))
            )
        } else {
            stop('path ', spec.dir, ' is invalid')
        }
    }
    # convert if calref is provided
    if (inherits(spec.dir, 'calref')) {
        out <- convert_calref(spec.dir)
        # add dark spectra
        if (!is.null(ref.dark.spec)) {
                ### lamp reference dark spectrum
                out$dat.ref.dark <- getSpec(ref.dark.spec, DOASmodel = out$DOAS.model, 
                    lite = TRUE)
        }
        if (!is.null(dark.spec)) {
                ### actual dark spectrum
                out$dat.dark <- getSpec(dark.spec, DOASmodel = out$DOAS.model, 
                    lite = TRUE)
        }
        if (!is.null(N2.dark.cal.spec)) {
                ### cuvette calibration dark spectrum
                out$dat.N2.dark <- getSpec(N2.dark.cal.spec, DOASmodel = out$DOAS.model, 
                    lite = TRUE)
        }
    } else {
        if (is.null(DOAS.model)) 
            stop('Please specify DOAS model (argument DOAS.model)!')
        if (is.null(ref.dark.spec)) ref.dark.spec <- spec.dir
        if (is.null(dark.spec)) dark.spec <- spec.dir
        if (is.null(NH3.cal.spec)) NH3.cal.spec <- spec.dir
        if (is.null(SO2.cal.spec)) SO2.cal.spec <- spec.dir
        if (is.null(NO.cal.spec)) NO.cal.spec <- spec.dir
        if (is.null(N2.dark.cal.spec)) N2.dark.cal.spec <- spec.dir
        if (is.null(N2.NH3.cal.spec)) N2.NH3.cal.spec <- spec.dir
        if (is.null(N2.SO2.cal.spec)) N2.SO2.cal.spec <- spec.dir
        if (is.null(N2.NO.cal.spec)) N2.NO.cal.spec <- spec.dir
        out <- list(
            DOAS.model = DOAS.model,
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
    }
    # add reference spectrum
    out$dat.ref <- if (!is.null(ref.spec)) {
        getSpec(ref.spec, DOASmodel = out$DOAS.model, lite = TRUE)
    } else {
        warning('no reference spectrum provided')
        NULL
    }
    out
}

### process spectra data
### ******************************************************************************
process_spectra <- function(specSet, rawData = NULL, correct.dark = TRUE, 
    correct.linearity = TRUE, correct.straylight = c("avg", "linear", "none"), 
    straylight.pix = NULL, return.subset.list = FALSE) {

    # get doas model
    doas_model <- specSet$DOAS.model

    # subset data?
    if (return.subset.list) {
        # get windows
        wins <- getWindows(specSet$DOAS.model)
        # get user defined straylight pixels
        if (!is.null(straylight.pix)) {
            wins$pixel_straylight <- straylight.pix
        }
        # re-define straylight.pix
        straylight.pix <- seq_along(wins$pixel_straylight)
        # remove doas model
        specSet$DOAS.model <- NULL
        # get subset index
        s_ind <- c(wins$pixel_straylight, wins$pixel_filter)
        # subset specSet
        specSet <- lapply(specSet, function(x) {
            out <- x
            out[[1]] <- out[[1]][s_ind]
            out
        })
        # subset rawData
        rawData[[1]] <- lapply(rawData[[1]], '[', s_ind)
    } else {
        if (is.null(straylight.pix)) {
            straylight.pix <- getWindows(specSet$DOAS.model)$pixel_straylight
        }
        # define s_ind for prediction of linear straylight
        s_ind <- seq_along(rawData[[1]][[1]])
    }

    ### dark-corrected reference spectra
    ### ******************************************************************************
    if (correct.dark) {
        out <- list(
            I.N2.NH3 = specSet$dat.N2.NH3[[1]] - specSet$dat.N2.dark[[1]],
            I.N2.SO2 = specSet$dat.N2.SO2[[1]] - specSet$dat.N2.dark[[1]],
            I.N2.NO = specSet$dat.N2.NO[[1]] - specSet$dat.N2.dark[[1]],
            I.NH3 = specSet$dat.NH3[[1]] - specSet$dat.N2.dark[[1]],
            I.SO2 = specSet$dat.SO2[[1]] - specSet$dat.N2.dark[[1]],
            I.NO = specSet$dat.NO[[1]] - specSet$dat.N2.dark[[1]],
            I.ref = specSet$dat.ref[[1]] - specSet$dat.ref.dark[[1]],
            I.meas = lapply(rawData[[1]], '-', specSet$dat.dark[[1]])
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
    if (correct.linearity) {
        if (is.null(rawData)) {
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
        out$I.meas <- lapply(out$I.meas, function(x) x / linearity.func(x, lin.coef))
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
            out$I.meas <- lapply(out$I.meas, function(x) x - mean(x[straylight.pix]))
        },
        {
            # NOTE: possibly introduced error with s_ind...
            require(MASS)
            out$I.N2.NH3 <- out$I.N2.NH3 - predict(rlm(out$I.N2.NH3[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.N2.SO2 <- out$I.N2.SO2 - predict(rlm(out$I.N2.SO2[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.N2.NO <- out$I.N2.NO - predict(rlm(out$I.N2.NO[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.NH3 <- out$I.NH3 - predict(rlm(out$I.NH3[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.SO2 <- out$I.SO2 - predict(rlm(out$I.SO2[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.NO <- out$I.NO - predict(rlm(out$I.NO[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.ref <- out$I.ref - predict(rlm(out$I.ref[straylight.pix]~straylight.pix,psi=psi.huber,maxit=1000),newdata=list(straylight.pix=s_ind))
            out$I.meas <- lapply(out$I.meas, function(x) {
                x - predict(
                    rlm(x[straylight.pix] ~ straylight.pix, psi = psi.huber, maxit = 1000),
                    newdata = list(straylight.pix = s_ind)
                )
            })
        },
        stop("argument correct.straylight must be either \"none\", \"avg\" or \"linear\"")
    )

    # remove straylight data
    if (return.subset.list) {
        out$I.N2.NH3 <- out$I.N2.NH3[-straylight.pix]
        out$I.N2.SO2 <- out$I.N2.SO2[-straylight.pix]
        out$I.N2.NO <- out$I.N2.NO[-straylight.pix]
        out$I.NH3 <- out$I.NH3[-straylight.pix]
        out$I.SO2 <- out$I.SO2[-straylight.pix]
        out$I.NO <- out$I.NO[-straylight.pix]
        out$I.ref <- out$I.ref[-straylight.pix]      
        out$I.meas <- lapply(out$I.meas, '[', -straylight.pix)
    }
    # remove empty I.meas entry if rawData is null
    if (!length(out$I.meas)) out$I.meas <- NULL 

    return(out)
}

### calculate differential spectra
### ******************************************************************************
diffSpecs <- function(Iset,use.ref=TRUE) {

    out <- list(
        NH3.diffspec = log(Iset$I.NH3 / Iset$I.N2.NH3),
        SO2.diffspec = log(Iset$I.SO2 / Iset$I.N2.SO2),
        NO.diffspec = log(Iset$I.NO / Iset$I.N2.NO)
    )

    if (!is.null(Iset$I.meas)) {
        if (use.ref) {
            suppressWarnings(out$diffspec <- lapply(Iset$I.meas, function(x) log(x / Iset$I.ref)))
        } else {
            suppressWarnings(out$diffspec <- lapply(Iset$I.meas, log))
        }
        out$diffspec[is.na(out$diffspec)] <- log(.Machine$double.eps)
    }

    return(out)
}

### get calibration related stuff
### ******************************************************************************
# TODO: tidy up later
getCalCurves <- function(diffspec, DOAS.win, calrefspec, warn = TRUE, ...) {

    # correct filter window for subsetted differential spectra
    if (length(diffspec$NH3.diffspec) < 800) {
        DOAS.win$pixel_filter <- seq_along(diffspec$NH3.diffspec)
    }

    # NH3:
    if (anyNA(c(calrefspec$dat.NH3$ambient$pathLength,calrefspec$dat.N2.NH3$ambient$pathLength))) {
        if (warn)warning("NH3 calibration spectrum: unknown path lengths might be different!")
        NH3_total <- calrefspec$dat.NH3$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NH3$cuvette$cuvetteLength
    } else {
        if (anyNA(c(calrefspec$dat.NH3$ambient$NH3ambient_ug,calrefspec$dat.N2.NH3$ambient$NH3ambient_ug))) {
            if (warn)warning("NH3 calibration spectrum: unknown ambient concentrations!")
            NH3_total <- calrefspec$dat.NH3$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NH3$cuvette$cuvetteLength
        } else {
            NH3_total <- calrefspec$dat.NH3$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NH3$cuvette$cuvetteLength + 
                calrefspec$dat.NH3$ambient$pathLength*calrefspec$dat.NH3$ambient$NH3ambient_ug -
                calrefspec$dat.N2.NH3$ambient$pathLength*calrefspec$dat.N2.NH3$ambient$NH3ambient_ug
        }
    }
    NH3.dc <- highpass.filter(diffspec$NH3.diffspec,DOAS.win,...)/NH3_total


    # SO2:
    if (anyNA(c(calrefspec$dat.SO2$ambient$pathLength,calrefspec$dat.N2.SO2$ambient$pathLength))) {
        if (warn)warning("SO2 calibration spectrum: unknown path lengths might be different!")
        SO2_total <- calrefspec$dat.SO2$cuvette$cuvetteConc_mg*1000*calrefspec$dat.SO2$cuvette$cuvetteLength
    } else {
        if (anyNA(c(calrefspec$dat.SO2$ambient$SO2ambient_ug,calrefspec$dat.N2.SO2$ambient$SO2ambient_ug))) {
            if (warn)warning("SO2 calibration spectrum: unknown ambient concentrations!")
            SO2_total <- calrefspec$dat.SO2$cuvette$cuvetteConc_mg*1000*calrefspec$dat.SO2$cuvette$cuvetteLength
        } else {
            SO2_total <- calrefspec$dat.SO2$cuvette$cuvetteConc_mg*1000*calrefspec$dat.SO2$cuvette$cuvetteLength + 
                calrefspec$dat.SO2$ambient$pathLength*calrefspec$dat.SO2$ambient$SO2ambient_ug -
                calrefspec$dat.N2.SO2$ambient$pathLength*calrefspec$dat.N2.SO2$ambient$SO2ambient_ug
        }
    }
    SO2.dc <- highpass.filter(diffspec$SO2.diffspec,DOAS.win,...)/SO2_total

    # NO:
    if (anyNA(c(calrefspec$dat.NO$ambient$pathLength,calrefspec$dat.N2.NO$ambient$pathLength))) {
        if (warn)warning("NO calibration spectrum: unknown path lengths might be different!")
        NO_total <- calrefspec$dat.NO$cuvette$cuvetteConc_mg*1000*calrefspec$dat.NO$cuvette$cuvetteLength
    } else {
        if (anyNA(c(calrefspec$dat.NO$ambient$NOambient_ug,calrefspec$dat.N2.NO$ambient$NOambient_ug))) {
            if (warn)warning("NO calibration spectrum: unknown ambient concentrations!")
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

fitConc <- function(meas.doascurve, DOAS.win, path.length, Cal.dc, robust=FALSE) {
    rob <- ".rob"[robust]
    fitcurve <- get(paste0("fit.curves",rob),mode="function")
    return(fitcurve(meas.doascurve,DOAS.win$pixel_dc,Cal.dc$Xreg[DOAS.win$pixel_dc,],DOAS.win$tau.shift,path.length))
}




evalOffline <- function(
    DOAS.model=NULL,
    evalperiod=NULL,
    tz.DOAS='Etc/GMT-1',
    path.length=NULL,
    tz.Output='Etc/GMT-1',
    filter.type= NULL,
    use.arima=FALSE,
    use.robust=FALSE,
    filter.window=NULL, 
    filter.strength=NULL,
    fit.window=NULL,
    straylight.window=NULL,
    skip.check.daily=FALSE,
    correct.linearity=TRUE,
    correct.straylight=c("avg","linear","none"),
    correct.dark=TRUE,
    use.ref = TRUE,
    avg.period=1,
    tau.shift=NULL,
    rawdata.dir=NULL,
    reference.dir=NULL,
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
    ncores=1,
    n_chunks = 1e2,
    add.name="",
    RawData = NULL,
    CalRefSpecs = NULL,
    Edinburgh_correction = TRUE,
    Serial = NULL,
    ...
    ) {


    ### verbose
    cat("\n************\nevaluating miniDOAS model", DOAS.model, "\n")

    ### initialize parallelism:
    if (parl <- (is.list(ncores) || ncores > 1)) {
        require(parallel)
        if (cl_was_up <- is.list(ncores)) {
            cl <- ncores
        } else {
            cl <- parallel::makePSOCKcluster(rep('localhost', ncores))
            on.exit(parallel::stopCluster(cl))
        }
    }

    ### read raw data:
    if (is.null(RawData)) {
        ### prepare time range:
        timerange <- prepTimeRange(evalperiod,tz.DOAS)

        ### get DOAS information:
        DOAS.info <- getDOASinfo(DOAS.model,timerange, Serial = Serial)

        cat("read raw data\n")
        RawData <- readDOASdata(DOAS.info, rawdata.dir, skip.check.daily = skip.check.daily, 
            force.write.daily = force.write.daily, ncores = ncores)
    } else {
        cat("raw data supplied - ignoring arguments 'evalperiod', 'DOAS.model', etc. ...\n")
        # NOTE: get time subset?
        ### get DOAS information:
        DOAS.info <- RawData[['DOASinfo']]
        ### prepare time range:
        timerange <- DOAS.info[['timerange']]
        # overwrite arguments provided with raw data
        tz.DOAS <- tzone(timerange)
        DOAS.model <- DOAS.info[['DOASmodel']]
    }

    ### set tzone output
    RawData$Header$st <- lubridate::with_tz(RawData$Header$st, tz.Output)
    RawData$Header$et <- lubridate::with_tz(RawData$Header$et, tz.Output)

    ### get DOAS windows:
    DOAS.win <- getWindows(DOAS.info, filter.type, timerange, 
        straylight.window, filter.window, fit.window, 
        filter.strength, tau.shift)

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
        if (is.character(CalRefSpecs)) {
            CalRefSpec <- switch(
                sub('.*[.]([a-z]+)$', '\\1', CalRefSpecs)
                , qs = {
                    if (!require(qs)) {
                        stop("package 'qs' needs to be installed")
                    }
                    qs::qread(CalRefSpecs)
                }
            )
        }
        # need to convert it?
        if (inherits(CalRefSpecs, 'calref')) {
            CalRefSpecs <- convert_calref(CalRefSpecs, ref.spec = ref.spec)
        } else if (!is.null(ref.spec)) {
            cat('ref.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.ref <- getSpec(ref.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(dark.spec)) {
            cat('dark.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.dark <- getSpec(dark.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(ref.dark.spec)) {
            cat('ref.dark.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.ref.dark <- getSpec(ref.dark.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(N2.dark.cal.spec)) {
            cat('N2.dark.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.N2.dark <- getSpec(N2.dark.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(NH3.cal.spec)) {
            cat('NH3.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.NH3 <- getSpec(NH3.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(SO2.cal.spec)) {
            cat('SO2.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.SO2 <- getSpec(SO2.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(NO.cal.spec)) {
            cat('NO.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.NO <- getSpec(NO.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(N2.NH3.cal.spec)) {
            cat('N2.NH3.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.N2.NH3 <- getSpec(N2.NH3.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(N2.SO2.cal.spec)) {
            cat('N2.SO2.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.N2.SO2 <- getSpec(N2.SO2.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
        if (!is.null(N2.NO.cal.spec)) {
            cat('N2.NO.cal.spec will be replaced with corresponding input\n')
            CalRefSpecs$dat.N2.NO <- getSpec(N2.NO.cal.spec, DOASmodel = DOAS.model, lite = TRUE)
        }
    }

    if (is.null(CalRefSpecs$dat.ref)) stop('Missing reference spectrum...')
    # if (is.null(CalRefSpecs$dat.dark)) stop('Missing dark spectrum...')
    if (is.null(CalRefSpecs$dat.NH3)) stop('Missing NH3 spectrum...')
    if (is.null(CalRefSpecs$dat.SO2)) stop('Missing SO2 spectrum...')
    if (is.null(CalRefSpecs$dat.NO)) stop('Missing NO spectrum...')

    ### calibration concentrations in ppb
    cat(
        sprintf(
            "\nNH3 calibration of %1.2f mg/m3 (corresponds to %1.2f ug/m3 in a path of %1.1f m)\n",
            CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg,
            CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg * 1000 / path.length * CalRefSpecs$dat.NH3$cuvette$cuvetteLength,
            path.length
            )
    )
    cat(
        sprintf(
            "SO2 calibration of %1.2f mg/m3 (corresponds to %1.2f ug/m3 in a path of %1.1f m)\n",
            CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg,
            CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg * 1000 / path.length * CalRefSpecs$dat.SO2$cuvette$cuvetteLength,
            path.length
            )
    )
    cat(
        sprintf(
            "NO calibration of %1.2f mg/m3 (corresponds to %1.2f ug/m3 in a path of %1.1f m)\n",
            CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg,
            CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg * 1000 / path.length * CalRefSpecs$dat.NO$cuvette$cuvetteLength,
            path.length
            )
    )

    ################################################################################
    ### data processing ############################################################
    ################################################################################
    
    cat('\nProcessing spectra...\n')
    # correct cal/ref specs:
    SpecCorr <- process_spectra(CalRefSpecs, rawData = RawData, correct.dark = correct.dark, 
        correct.linearity = correct.linearity, correct.straylight = correct.straylight, 
        straylight.pix = DOAS.win$pixel_straylight, return.subset.list = TRUE)

    # Fix S1
    if (DOAS.model=="S1" && RawData$Header[1,"st"] > strptime("20160614",format="%Y%m%d") && RawData$Header[1,"st"] < strptime("20170101",format="%Y%m%d")) {
        RawData$Header[,"TECTemp"] <- sapply(strsplit(RawData$Header[,"Klima"],","),"[",4)
    }

    # calculate differential spectra:
    DiffSpec <- diffSpecs(SpecCorr, use.ref = use.ref)

    # get calibration doascurves:
    Cal.dc <- getCalCurves(DiffSpec, DOAS.win, CalRefSpecs, ...)

    # get fitting function:
    if (use.robust) {
        fitcurve <- fit.curves.rob
    } else {
        fitcurve <- fit.curves
    }

    # create a 'lighter' highpass filter
    winFUN <- getOption('md.filter.function.list')[[DOAS.win$filter.type]]
    DOAS.win$filt <- winFUN(DOAS.win$filter.strength,...)
    # new double filter
    highpass.filter2 <- function(dat, filt) {
        C_cfilter <- getFromNamespace('C_cfilter', 'stats')
        dat - (
            .Call(C_cfilter, dat, filt, 2L, FALSE) +
                rev(.Call(C_cfilter, rev(dat), filt, 2L, FALSE))
            ) / 2
    }


    ### ******************************************************************************
    ### process raw-data record-wise
    ### ******************************************************************************

    if (parl) {

        # verbose
        cat("\nParallel computing doascurve and fit...\n\n")

        # # call libraries for robust fit
        # if (use.robust) {
        #     parallel::clusterCall(cl, library, package = 'robustbase', character.only = TRUE)
        #     parallel::clusterCall(cl, library, package = 'MASS', character.only = TRUE)
        # }

        # get calibration curves
        xreg <- Cal.dc$Xreg[DOAS.win$pixel_dc, ]

        # export functions & objects
        hpf2 <- compiler::cmpfun(highpass.filter2)
        diff_specs <- DiffSpec$diffspec
        # FIXME: fix parallel comp, move functions outside, remove process_fit, change fitparallel
        fit_parallel <- function(index, doas_win, path_length) {
            # xreg & diff_specs need to be exported
            lapply(index, function(i) {
                # check if light
                if (isTRUE(median(diff_specs[[i]], na.rm = TRUE) > -5)) {
                    # highpass filter and fit curve to calibration
                    fitcurve(
                        hpf2(diff_specs[[i]], doas_win$filt), 
                        doas_win$pixel_dc, xreg, doas_win$tau.shift, path_length)
                } else {
                    # return NAs if too few light
                    as.list(c(
                            rep(NA_real_, 6),
                            NA_integer_
                            ))
                }
            })
        }
        parallel::clusterExport(cl, 
            list('hpf2', 'fitcurve','fit_parallel', 'diff_specs', 'xreg'), 
            envir = environment())

        # parallel calculation
        cat("This might take a while...\n\n")
        if (requireNamespace('bLSmodelR')) {
            # split diffspec indices into 
            j_frac <- ceiling(length(diff_specs) / n_chunks)
            lcl <- ceiling(j_frac / length(cl))
            p_index <- parallel::splitIndices(length(DiffSpec$diffspec), length(cl) * lcl)
            out <- bLSmodelR:::.clusterApplyLB(cl, p_index,
                fit_parallel, doas_win = DOAS.win, path_length = path.length)
        } else {
            # split diffspec indices
            p_index <- parallel::splitIndices(length(DiffSpec$diffspec), length(cl))
            # run in chunks
            out <- parallel::clusterApply(cl, p_index, fit_parallel, 
                doas_win = DOAS.win, path_length = path.length)
        }

        if (!cl_was_up) {
            parallel::stopCluster(cl)
            on.exit()
        }

    } else {

        ### sequential
        ### ******************************************************************************#

        # get # of files to process
        n_files <- length(DiffSpec$diffspec)

        # prepare xreg
        xreg <- Cal.dc$Xreg[DOAS.win$pixel_dc, ]

        # loop over files
        out <- lapply(seq.int(n_files), function(i) {
            # verbose
            cat("\r", i, "/", n_files)
            # check if light
            if (isTRUE(median(DiffSpec$diffspec[[i]], na.rm = TRUE) > -5)) {
                # highpass filter and fit curve to calibration
                fitcurve(
                    highpass.filter2(DiffSpec$diffspec[[i]], DOAS.win$filt), 
                    DOAS.win$pixel_dc, xreg, DOAS.win$tau.shift, path.length, 
                    return_resid = !lite)
            } else {
                # return NAs if too few light
                out <- as.list(c(
                        rep(NA_real_, 6),
                        NA_integer_
                        ))
                if (!lite) {
                    out <- c(out, list(rep(NA_real_, length(DOAS.win$pixel_dc))))
                }
                out
            }
            })
        cat("\n")

    }

    if (lite) {
        # rbind results
        out <- rbindlist(out)
    } else {
        # residuals
        residuals <- sapply(out, '[[', 8)
        colnames(residuals) <- format(RawData$Header[, 'st'])
        # rbind results
        out <- rbindlist(lapply(out, '[', -8))
    }



    # name them
    setnames(out, c('nh3', 'so2', 'no', 'nh3_se', 'so2_se', 'no_se', 'tau'))

    # correct NH3 calibration with Edinburgh correction?
    if (Edinburgh_correction) {
        out[, ':='(
            nh3 = nh3 * 1.16,
            nh3_se = nh3_se * 1.16
            )]
    }

    # gather residual results
    results <- cbind(
        st = RawData$Header[, "st"],
        et = RawData$Header[, "et"],
        out,
        i_max = sapply(RawData$RawData, max, na.rm = TRUE),
        n = RawData$Header[,"AccNum"],
        shutter = RawData$Header[, "ShuPos"],
        revolver = RawData$Header[, "RevPos"],
        tec = as.numeric(RawData$Header[, "TECTemp"]),
        i_fit_min = sapply(RawData$RawData, function(x) min(x[DOAS.win$pixel_fit], na.rm = TRUE)),
        i_fit_max = sapply(RawData$RawData, function(x) max(x[DOAS.win$pixel_fit], na.rm = TRUE))
    )

    # verbose
    cat("\nfinished evaluation of",DOAS.model,"\n")

    # return results
    if (lite) {
        results
    } else {
        structure(
            results
            ,CalRefSpecs = CalRefSpecs
            ,RawData = RawData
            ,residuals = residuals
            ,DOASwin = DOAS.win
            ,callEval = match.call()
            ,class = c("DOASeval", 'data.table', "data.frame")
        )
    }
}

# add resid method
residuals.DOASeval <- function(object, ...) {
    structure(attr(object, 'residuals')
        , DOASinfo = attr(object, 'RawData')[['DOASinfo']]
        , DOASwin = attr(object, 'DOASwin')
        , class = 'DOAS.resid'
        )
}


runAppDOAS <- function(ShinyInput) {
    require(shiny)
    require(dygraphs)
    require(xts)
    require(lubridate)
    if (is.ibts(ShinyInput)) {
        ShinyInput <- as.xts(ShinyInput)
    }
    dat <- cbind(ShinyInput,"none"=NA)
    r10 <- if (all(is.na(dat[,1]))) c(-1,1) else range(dat[,1],na.rm=TRUE)
    dr1 <- diff(r10)
    r1_ini <- round(r10 + c(-1,1)*0.1*ifelse(dr1<1,1,dr1),1)
    r20 <- if (all(is.na(dat[,2]))) c(-1,1) else range(dat[,2],na.rm=TRUE)
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
                        if (!is.null(input$timeseries_date_window)) {
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
                        if (drawy1) {
                            r10 <- if (all(is.na(ydat[,yax1]))) c(-1,1) else range(ydat[,yax1],na.rm=TRUE)
                            dr1 <- diff(r10)
                            r1 <- round(range(ydat[,yax1],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr1<1,1,dr1),1)
                            if (input$checkzero&r1[1]>0)r1[1]<-0
                            if (input$checkzero&r1[2]<0)r1[2]<-0
                        }else{
                            r1 <- c(-1,1) + c(-1,1)*0.04
                        }
                        if (drawy2) {
                            r20 <- if (all(is.na(ydat[,yax2]))) c(-1,1) else range(ydat[,yax2],na.rm=TRUE)
                            dr2 <- diff(r20)
                            r2 <- round(range(ydat[,yax2],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr2<1,1,dr2),1)
                            if (input$checkzero&r2[1]>0)r2[1]<-0
                            if (input$checkzero&r2[2]<0)r2[2]<-0
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
                        if (!is.null(input$timeseries_date_window)) {
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
                        if (drawy1) {
                            ylim1 <- ylim1()
                            if (input$checkzero&ylim1[1]>0)ylim1[1]<-0
                            if (input$checkzero&ylim1[2]<0)ylim1[2]<-0
                            plot(ydat[,1],ann=FALSE,axes=FALSE,auto.grid=FALSE,type="n",main="",ylim=ylim1,xlim=xlim1)
                            axis(1,at=px,labels=strftime(px,"%d.%m %H:%M"))
                            grid(nx=NA,ny=NULL)
                            abline(v=px,lty="dotted",col="lightgrey")
                            for(i in yax1)lines(ydat[,i],col=c("yellowgreen","dodgerblue","orange")[i])
                            axis(2,at=pretty(ylim1),labels=pretty(ylim1))
                            mtext(ylab1(),side=2,line=2.5,cex=1.25)
                            ploty <- TRUE
                        }
                        if (drawy2) {
                            ylim2 <- ylim2()
                            if (input$checkzero&ylim2[1]>0)ylim2[1]<-0
                            if (input$checkzero&ylim2[2]<0)ylim2[2]<-0
                            if (ploty) {
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
                        if (drawy1) {
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
                        if (drawy1) {
                            input$sl2
                        } else {
                            NULL
                        }
                    })

                    output$timeseries <- renderDygraph({
                        ydat <- cbind(dat[,input$var1],dat[,input$var2],dat[,input$var3])
                        # Times <- as.POSIXct(ydat)
                        # if (!is.null(input$timeseries_date_window)) {
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
                        # if (drawy1) {
                        #   r10 <- if (all(is.na(ydat[,yax1]))) c(-1,1) else range(ydat[,yax1],na.rm=TRUE)
                        #   dr1 <- diff(r10)
                        #   r1 <- round(range(ydat[,yax1],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr1<1,1,dr1),1)
                        #           if (input$checkzero&r1[1]>0)r1[1]<-0
                        #           if (input$checkzero&r1[2]<0)r1[2]<-0
                        # }else{
                        #   r1 <- c(-1,1) + c(-1,1)*0.04
                        # }
                        # if (drawy2) {
                        #   r20 <- if (all(is.na(ydat[,yax2]))) c(-1,1) else range(ydat[,yax2],na.rm=TRUE)
                        #   dr2 <- diff(r20)
                        #   r2 <- round(range(ydat[,yax2],na.rm=TRUE) + c(-1,1)*0.1*ifelse(dr2<1,1,dr2),1)
                        #           if (input$checkzero&r2[1]>0)r2[1]<-0
                        #           if (input$checkzero&r2[2]<0)r2[2]<-0
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
                        filename=function() {paste0("Shinysave_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png")},
                        content = function(file) {
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


inspectEvaluation <- function(rawdat,CalRefSpecs, path.length, index = 1,
    filter.type = NULL, fit.type = "OLS", robust = TRUE, straylight.window = NULL, filter.window = NULL, fit.window = NULL, 
    filter.strength = NULL, tau.shift = NULL, correct.dark = TRUE, correct.linearity = TRUE, 
    correct.straylight = c("avg", "linear", "none"), use.ref = TRUE,
    Edinburgh_correction = TRUE) {

    require(shiny)
    require(shinyWidgets)
    require(lubridate)
    require(robustbase)
    require(MASS)

    # windows first time:
    DOASwindows <- getWindows(rawdat$DOASinfo,filter.type = filter.type, straylight.window = straylight.window, 
        filter.window = filter.window, fit.window = fit.window, filter.strength = filter.strength,
        tau.shift = tau.shift)

    # correct cal/ref specs:
    SpecCorr <- process_spectra(CalRefSpecs,NULL, correct.dark = correct.dark, correct.linearity = correct.linearity, 
        correct.straylight = correct.straylight, straylight.pix=DOASwindows$pixel_straylight)

    # diffspec:
    DiffSpec <- diffSpecs(SpecCorr,use.ref=use.ref)

    # cal doascurves:
    Cal.dc <- getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE)

    wavelength <- rawdat$DOASinfo$Spectrometer$wavelength[DOASwindows$pixel_filter]

    index.max <- length(rawdat$RawData)

    frmls <- formals(getOption('md.filter.function.list')[[DOASwindows$filter.type]])
    n <- DOASwindows$filter.strength
    p1 <- eval(frmls$p1)
    p2 <- eval(frmls$p2)
    p3 <- eval(frmls$p3) 
    p4 <- eval(frmls$p4)

    # input <- c(list(
    #     index = 1,
    #     robust = robust,
    #     correct.dark = correct.dark,
    #     correct.linearity = correct.linearity,
    #     correct.straylight = correct.straylight,
    #     use.ref = use.ref
    #     ),
    #     DOASwindows
    #     )

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
                            ,numericInput("filter.strength", 
                              label = "filter.strength:",
                              step = 2,
                              value = DOASwindows$filter.strength
                            )
                            ,selectInput("filter.type", 
                                label = "filter.type:",
                                choices = names(getOption('md.filter.function.list')),
                                selected = DOASwindows$filter.type
                            )
                            ,numericRangeInput("filter.window", 
                                label = "filter.window:", 
                                min = 190, 
                                max = 240,
                                step = 0.2, 
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

                    # windows:
                    DOASwindows_reactive <- reactive({
                        if (is.null(input$fit.window)) {
                            fit_window <- DOASwindows$fit.window
                        } else {
                            fit_window <- input$fit.window
                        }
                        if (is.null(input$filter.strength)) {
                            filter.strength <- DOASwindows$filter.strength
                        } else {
                            filter.strength <- input$filter.strength
                        }
                        if (filter.strength < 3) {
                                filter_strength <- getOption('md.filter.strength')
                        } else {
                            filter_strength <- filter.strength
                        }
                        getWindows(rawdat$DOASinfo, filter.type = input$filter.type, straylight.window = straylight.window, 
                            filter.window = input$filter.window, fit.window = fit_window,
                            filter.strength = filter_strength, tau.shift = input$tau.shift)
                    })

                    output$filterArgs <- renderUI({
                        DOASwindows <- DOASwindows_reactive()
                            if (is.null(input$filter.strength)) {
                                filter.strength <- DOASwindows$filter.strength
                            } else {
                                filter.strength <- input$filter.strength
                            }
                            frmls <- formals(getOption('md.filter.function.list')[[DOASwindows$filter.type]])
                            n <- filter.strength
                            p1 <- eval(frmls$p1)
                            p2 <- eval(frmls$p2)
                            p3 <- eval(frmls$p3)
                            p4 <- eval(frmls$p4)

                            tagList(
                                if (!is.null(p1))numericInput("p1", 
                                    label = "p1:",
                                    step = 0.1,
                                    value = p1
                                )
                                ,if (!is.null(p2))numericInput("p2", 
                                    label = "p2:",
                                    step = 0.1,
                                    value = p2
                                )
                                ,if (!is.null(p3))numericInput("p3", 
                                    label = "p3:",
                                    step = 0.1,
                                    value = p3
                                )
                                ,if (!is.null(p4))numericInput("p4", 
                                    label = "p4:",
                                    step = 0.1,
                                    value = p4
                                )
                            )
                        }            
                    )

                    output$fit_window <- renderUI({
                        DOASwindows <- DOASwindows_reactive()
                        numericRangeInput("fit.window", 
                            label = "fit.window:", 
                            min = DOASwindows$filter.window[1], 
                            max = DOASwindows$filter.window[2],
                            step = 0.2, 
                            value = DOASwindows$fit.window
                        )
                    })


                    plot_reactive <- reactive({
                        i <- input$index
                        RawDat <- rawdat
                        RawDat[[1]] <- list(rawdat[[1]][[i]])
                        DOASwindows <- DOASwindows_reactive()
                        SpecCorr <- process_spectra(CalRefSpecs,RawDat,correct.dark = input$correct.dark, 
                            correct.linearity = input$correct.linearity, correct.straylight = input$correct.straylight, 
                            straylight.pix=DOASwindows$pixel_straylight)
                        DiffSpec <- diffSpecs(SpecCorr,use.ref=input$use.ref)
                        Cal.dc <- getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE,input$p1,input$p2,input$p3,input$p4)
                        wavelength <- rawdat$DOASinfo$Spectrometer$wavelength[DOASwindows$pixel_filter]


                        # plot(1,main="diffspec")
                        # plot(1,main="doascurve")
                        # plot(1,main="darkspec")
                        # reactive:
                        xlim <- range(wavelength)
                        # main <- rawdat$Header[i,]

                        # I.meas + I.ref
                        ylim <- range(SpecCorr$I.meas[[1]][DOASwindows$pixel_filter],SpecCorr$I.ref[DOASwindows$pixel_filter],na.rm=TRUE)
                        if (!all(is.finite(ylim))) {
                            ylim <- c(0.01,2)
                        } else if (any(ylimBelow <- ylim<0)) {
                            ylim[ylimBelow] <- c(0.01,2)[ylimBelow]
                        } 
                        # windows(width=10,height=7)
                        plot(1,1,xlim=xlim,ylim=ylim,log="y",type="n",ylab="counts",xlab="",main="spectra")
                        if (input$use.ref)lines(wavelength,SpecCorr$I.ref[DOASwindows$pixel_filter],col="darkgrey",lwd=2)
                        lines(wavelength,SpecCorr$I.meas[[1]][DOASwindows$pixel_filter],col="black",lwd=2)
                        legend("bottomright",c("meas.","ref."),lwd=2,col=c("black","darkgrey"),bty="n")

                        # log(I.meas/I.ref)
                        ylim <- range(DiffSpec$diffspec[[1]][DOASwindows$pixel_filter],na.rm=TRUE)
                        if (!all(is.finite(ylim))) ylim <- c(0,1)
                        meas.dc <- highpass.filter(DiffSpec$diffspec[[1]],DOASwindows,input$p1,input$p2,input$p3,input$p4)
                        # isna <- is.na(meas.dc)
                        # windows(width=10,height=7)
                        plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.meas/I.ref)",xlab="",main="diffspec")
                        lines(wavelength,DiffSpec$diffspec[[1]][DOASwindows$pixel_filter] - meas.dc,lwd=2,col="orange")
                        lines(wavelength,DiffSpec$diffspec[[1]][DOASwindows$pixel_filter],col="black")
                        fit <- fitConc(meas.dc, DOASwindows, path.length, Cal.dc, robust=input$robust)

                        fit.NH3 <- fit[[1]]*path.length*Cal.dc$NH3.dc
                        fit.SO2 <- fit[[2]]*path.length*Cal.dc$SO2.dc
                        fit.NO <- fit[[3]]*path.length*Cal.dc$NO.dc

                        if (input$tau.shift>0) {
                            meas.dc <- c(meas.dc[-seq.int(input$tau.shift)],rep(NA,abs(input$tau.shift)))
                        } else if (input$tau.shift<0) {
                            tau.shift <- abs(input$tau.shift)
                            meas.dc <- c(rep(NA,tau.shift),meas.dc[-(length(meas.dc) - seq.int(tau.shift) + 1)])
                        }

                        # doascurve
                        ylim <- range(meas.dc,na.rm=TRUE)
                        if (!all(is.finite(ylim))) ylim <- c(0,1)
                        # windows(width=10,height=7)
                        plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)},main="doascurve")
                        lines(wavelength,fit.SO2,col="#00bb00aa")
                        lines(wavelength,fit.SO2+fit.NO,col="#ff0000dd")
                        lines(wavelength,fit.SO2+fit.NO+fit.NH3,col="blue")
                        lines(wavelength[DOASwindows$pixel_dc],fit.SO2[DOASwindows$pixel_dc],lwd=2,col="#00bb00aa")
                        lines(wavelength[DOASwindows$pixel_dc],fit.SO2[DOASwindows$pixel_dc]+fit.NO[DOASwindows$pixel_dc],lwd=2,col="#ff0000dd")
                        lines(wavelength[DOASwindows$pixel_dc],fit.SO2[DOASwindows$pixel_dc]+fit.NO[DOASwindows$pixel_dc]+fit.NH3[DOASwindows$pixel_dc],lwd=2,col="blue")
                        abline(v=wavelength[range(DOASwindows$pixel_dc)],lty=3)
                        legend("topright",c("SO2","SO2 + NO","SO2 + NO + NH3"),lty=1,col=c("#00bb00aa","#ff0000dd","blue"),cex=0.7)

                        # residual
                        # windows(width=10,height=7)
                        plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="residuals [-]",xlab="",panel.first={grid();abline(h=0)},col="darkgrey",main="residuals")
                        lines(wavelength,meas.dc-fit.SO2-fit.NO-fit.NH3) 
                        lines(wavelength[DOASwindows$pixel_dc],meas.dc[DOASwindows$pixel_dc]-fit.SO2[DOASwindows$pixel_dc]-fit.NO[DOASwindows$pixel_dc]-fit.NH3[DOASwindows$pixel_dc],lwd=2) 
                        abline(v=wavelength[range(DOASwindows$pixel_dc)],lty=3)

                        # cal diffspec
                        ylim <- range(DiffSpec$NH3.diffspec[DOASwindows$pixel_filter],na.rm=TRUE)
                        plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.NH3/I.N2)",xlab="",panel.first={grid();abline(h=0)})
                        # #
                        # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel_filter] - Cal.dc$Xreg[,3]*CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NO$cuvette$cuvetteLength,lwd=2,col="red")
                        # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel_filter])
                        # #
                        # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel_filter] - Cal.dc$Xreg[,2]*CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.SO2$cuvette$cuvetteLength,lwd=2,col="green")
                        # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel_filter])
                        # #
                        lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel_filter] - Cal.dc$Xreg[,1]*CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NH3$cuvette$cuvetteLength,lwd=2,col="blue")
                        lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel_filter])

                        # cal doascurve
                        ylim <- range(Cal.dc$Xreg,na.rm=TRUE)
                        plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="cal. DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)})
                        lines(wavelength,Cal.dc$Xreg[,3],lwd=2,col="red")
                        lines(wavelength,Cal.dc$Xreg[,2],lwd=2,col="green")
                        lines(wavelength,Cal.dc$Xreg[,1],lwd=2,col="blue")

                        if (Edinburgh_correction) {
                            fit[[1]] <- fit[[1]] * 1.16
                        }

                        msg1 <- sprintf("index %i/%i:  %s  --  NH3: %1.1f +/- %1.1f  --  SO2: %1.1f +/- %1.1f  --  NO: %1.1f +/- %1.1f",
                            i,index.max,format(rawdat$Header[i,"st"]),
                            fit[[1]],fit[[4]],fit[[2]],fit[[5]],fit[[3]],fit[[6]])
                        # mtext(as.expression(substitute(italic(msg), list(msg=msg))), line=-1.25, outer=TRUE, cex=0.5)
                        # mtext(substitute(italic(msg), list(msg=msg1)), line=1, outer=TRUE, cex=0.75)
                        # mtext(substitute(italic(msg), list(msg=msg2)), line=-0.5, outer=TRUE, cex=0.75)
                        mtext(msg1, line=0, outer=TRUE)#, cex=0.75)

                    })

                    output$specs <- renderPlot({
                        par(mfrow=c(3,2),oma=c(0,0,1.5,0))
                        plot_reactive()
                    })
                    output$downloadPlot <- downloadHandler(
                        filename=function() {paste0("shinyInspect_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png")},
                        content = function(file) {
                            png(file,width=1200,height=1200)
                            on.exit(dev.off())

                            par(mfrow=c(3,2),oma=c(0,0,1.5,0))
                            i <- input$index
                            RawDat <- rawdat
                            RawDat[[1]] <- rawdat[[1]][[i]]
                            DOASwindows <- DOASwindows_reactive()
                            SpecCorr <- process_spectra(CalRefSpecs,RawDat,correct.dark = input$correct.dark, correct.linearity = input$correct.linearity, 
                                correct.straylight = input$correct.straylight, straylight.pix=DOASwindows$pixel_straylight)
                            DiffSpec <- diffSpecs(SpecCorr,use.ref=input$use.ref)
                            Cal.dc <- getCalCurves(DiffSpec,DOASwindows,CalRefSpecs,warn=FALSE,input$p1,input$p2,input$p3,input$p4)
                            wavelength <- rawdat$DOASinfo$Spectrometer$wavelength[DOASwindows$pixel_filter]


                            # plot(1,main="diffspec")
                            # plot(1,main="doascurve")
                            # plot(1,main="darkspec")
                            # reactive:
                            xlim <- range(wavelength)
                            # main <- rawdat$Header[i,]

                            # I.meas + I.ref
                            ylim <- range(SpecCorr$I.meas[DOASwindows$pixel_filter,],SpecCorr$I.ref[DOASwindows$pixel_filter],na.rm=TRUE)
                            if (!all(is.finite(ylim))) {
                                ylim <- c(0.01,2)
                            } else if (any(ylimBelow <- ylim<0)) {
                                ylim[ylimBelow] <- c(0.01,2)[ylimBelow]
                            } 
                            # windows(width=10,height=7)
                            plot(1,1,xlim=xlim,ylim=ylim,log="y",type="n",ylab="counts",xlab="",main="spectra")
                            if (input$use.ref)lines(wavelength,SpecCorr$I.ref[DOASwindows$pixel_filter],col="darkgrey",lwd=2)
                            lines(wavelength,SpecCorr$I.meas[DOASwindows$pixel_filter,],col="black",lwd=2)
                            legend("bottomright",c("meas.","ref."),lwd=2,col=c("black","darkgrey"),bty="n")

                            # log(I.meas/I.ref)
                            ylim <- range(DiffSpec$diffspec[DOASwindows$pixel_filter,],na.rm=TRUE)
                            if (!all(is.finite(ylim))) ylim <- c(0,1)
                            meas.dc <- highpass.filter(DiffSpec$diffspec,DOASwindows,input$p1,input$p2,input$p3,input$p4)
                            # isna <- is.na(meas.dc)
                            # windows(width=10,height=7)
                            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.meas/I.ref)",xlab="",main="diffspec")
                            lines(wavelength,DiffSpec$diffspec[DOASwindows$pixel_filter,] - meas.dc,lwd=2,col="orange")
                            lines(wavelength,DiffSpec$diffspec[DOASwindows$pixel_filter,],col="black")
                            fit <- fitConc(meas.dc, DOASwindows, path.length, Cal.dc, robust=input$robust)

                            cfs <- fit[[5]]
                            fit.SO2 <- cfs[2]*path.length*Cal.dc$SO2.dc
                            fit.NO <- cfs[3]*path.length*Cal.dc$NO.dc
                            fit.NH3 <- cfs[1]*path.length*Cal.dc$NH3.dc

                            # fit.SO2[isna] <- NA
                            # fit.NO[isna] <- NA
                            # fit.NH3[isna] <- NA

                            if (input$tau.shift>0) {
                                meas.dc <- c(meas.dc[-seq.int(input$tau.shift)],rep(NA,abs(input$tau.shift)))
                            } else if (input$tau.shift<0) {
                                tau.shift <- abs(input$tau.shift)
                                meas.dc <- c(rep(NA,tau.shift),meas.dc[-(length(meas.dc) - seq.int(tau.shift) + 1)])
                            }

                            # doascurve
                            ylim <- range(meas.dc,na.rm=TRUE)
                            if (!all(is.finite(ylim))) ylim <- c(0,1)
                            # windows(width=10,height=7)
                            plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)},main="doascurve")
                            lines(wavelength,fit.SO2,col="#00bb00aa")
                            lines(wavelength,fit.SO2+fit.NO,col="#ff0000dd")
                            lines(wavelength,fit.SO2+fit.NO+fit.NH3,col="blue")
                            lines(wavelength[DOASwindows$pixel_dc],fit.SO2[DOASwindows$pixel_dc],lwd=2,col="#00bb00aa")
                            lines(wavelength[DOASwindows$pixel_dc],fit.SO2[DOASwindows$pixel_dc]+fit.NO[DOASwindows$pixel_dc],lwd=2,col="#ff0000dd")
                            lines(wavelength[DOASwindows$pixel_dc],fit.SO2[DOASwindows$pixel_dc]+fit.NO[DOASwindows$pixel_dc]+fit.NH3[DOASwindows$pixel_dc],lwd=2,col="blue")
                            abline(v=wavelength[range(DOASwindows$pixel_dc)],lty=3)
                            legend("topright",c("SO2","SO2 + NO","SO2 + NO + NH3"),lty=1,col=c("#00bb00aa","#ff0000dd","blue"),cex=0.7)

                            # residual
                            # windows(width=10,height=7)
                            plot(wavelength,meas.dc,xlim=xlim,ylim=ylim,type="l",ylab="residuals [-]",xlab="",panel.first={grid();abline(h=0)},col="darkgrey",main="residuals")
                            lines(wavelength,meas.dc-fit.SO2-fit.NO-fit.NH3) 
                            lines(wavelength[DOASwindows$pixel_dc],meas.dc[DOASwindows$pixel_dc]-fit.SO2[DOASwindows$pixel_dc]-fit.NO[DOASwindows$pixel_dc]-fit.NH3[DOASwindows$pixel_dc],lwd=2) 
                            abline(v=wavelength[range(DOASwindows$pixel_dc)],lty=3)

                            # cal diffspec
                            ylim <- range(DiffSpec$NH3.diffspec[DOASwindows$pixel_filter],na.rm=TRUE)
                            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="log(I.NH3/I.N2)",xlab="",panel.first={grid();abline(h=0)})
                            # #
                            # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel_filter] - Cal.dc$Xreg[,3]*CalRefSpecs$dat.NO$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NO$cuvette$cuvetteLength,lwd=2,col="red")
                            # lines(wavelength,DiffSpec$NO.diffspec[DOASwindows$pixel_filter])
                            # #
                            # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel_filter] - Cal.dc$Xreg[,2]*CalRefSpecs$dat.SO2$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.SO2$cuvette$cuvetteLength,lwd=2,col="green")
                            # lines(wavelength,DiffSpec$SO2.diffspec[DOASwindows$pixel_filter])
                            # #
                            lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel_filter] - Cal.dc$Xreg[,1]*CalRefSpecs$dat.NH3$cuvette$cuvetteConc_mg*1000*CalRefSpecs$dat.NH3$cuvette$cuvetteLength,lwd=2,col="blue")
                            lines(wavelength,DiffSpec$NH3.diffspec[DOASwindows$pixel_filter])

                            # cal doascurve
                            ylim <- range(Cal.dc$Xreg,na.rm=TRUE)
                            plot(1,1,xlim=xlim,ylim=ylim,type="n",ylab="cal. DOAS curve [-]",xlab="",panel.first={grid();abline(h=0)})
                            lines(wavelength,Cal.dc$Xreg[,3],lwd=2,col="red")
                            lines(wavelength,Cal.dc$Xreg[,2],lwd=2,col="green")
                            lines(wavelength,Cal.dc$Xreg[,1],lwd=2,col="blue")

                            msg1 <- sprintf("index %i/%i:  %s  --  NH3: %1.1f +/- %1.1f  --  SO2: %1.1f +/- %1.1f  --  NO: %1.1f +/- %1.1f",i,index.max,format(rawdat$Header[i,"st"]),cfs[1],fit[[6]][1],cfs[2],fit[[6]][2],cfs[3],fit[[6]][3])
                            # mtext(as.expression(substitute(italic(msg), list(msg=msg))), line=-1.25, outer=TRUE, cex=0.5)
                            # mtext(substitute(italic(msg), list(msg=msg1)), line=1, outer=TRUE, cex=0.75)
                            # mtext(substitute(italic(msg), list(msg=msg2)), line=-0.5, outer=TRUE, cex=0.75)
                            mtext(msg1, line=0, outer=TRUE)#, cex=0.75)

                        },
                        contentType="image/png"
                    )
                }
            } 
        )
    )
}


# ### create synthetic theoretical calibration DOAS curves only based on published absorption crosssection and compare that with calibration DOAS curve
# ### ******************************************************************************
# f <- DOAS.info$Spectrometer$wavelength
# sft <- -9
# x0 <- DOAS.win$pixel_filter + sft
# f2 <- f[x0]
# f3 <- (f[c(x0[1]-1,x0)]+f[c(x0[1],x0+1)])/2
# if (DOAS.model!="S1") {
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
#   function(x,i) {
#     integrate(
#       function(x) {
#         out <- approx(f2.NH3_filt,abscrossec.NH3_cuv,xout=x)$y
#         out[is.na(out)] <- 0
#         return(out)
#       }
#       ,x[i],x[i+1])$value/diff(x[i+(0:1)])
#   },x=f3)
# NH3.doascurve.theor <- highpass.filter(cs.Theor,filter.type, DOAS.win$filter.strength,b,Scale,delta, filter.strength_multiplication, filter.strength_loess,fam)/(NH3.cal*dat.NH3$cuvette[[1]])
# regr.NH3 <- suppressWarnings(lmrob(NH3.doascurve.theor[DOAS.win$pixel_dc + sft] ~ NH3.doascurve[DOAS.win$pixel_dc + sft], setting="KS2014"))
# theory.NH3.factor <- round(coefficients(regr.NH3)[2],3)

# plot(f[DOAS.win$pixel_filter],NH3.doascurve.theor,type="n",ylab="NH3 doascurve",xlab="",panel.first={grid();abline(h=0,col="darkgrey")})
# lines(f[DOAS.win$pixel_filter],NH3.doascurve,lwd=2)
# lines(f[DOAS.win$pixel_filter],NH3.doascurve.theor,col="blue",lwd=2)
# lines(f[DOAS.win$pixel_filter],NH3.doascurve*theory.NH3.factor,col="red",lwd=2)
# legend("bottomright",c("NH3.doascurve","Chen-Theory",sprintf("NH3.doascurve x %1.3f",theory.NH3.factor)),lwd=2,col=c("black","blue","red"))

# ind <- seq.int(320)
# par(mfrow=c(2,1))
# plot(f[DOAS.win$pixel_filter][ind],NH3.doascurve[ind]*theory.NH3.factor,type="l",ylab="NH3 doascurve (m2/ug)",xlab="",lwd=2,panel.first={grid();abline(h=0,col="darkgrey")})
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
