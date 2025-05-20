

### function for reading GasFinder3.0 data:
read_gasfinder <- function(Folder = "", From = NULL, To = NULL, tz = "Etc/GMT-1", 
    DTA_only = TRUE, asCharacter = FALSE, time_add = 0) {

	# parse time to add
	time_add <- parse_time_diff(time_add)

	if((From_null <- is.null(From)) & (To_null <- is.null(To))){
		if(R.utils::isDirectory(Folder[1])){
			cat("Please choose raw data files...\n")
			Files <- tcltk::tclvalue(tcltk::tkgetOpenFile(initialdir=Folder,title="Select GF3 Data...",filetypes="{GF3_Data {*.LOG}}",multiple=TRUE))
			if(Files=="")stop("No files selected.")
			Files <- gsub("[{]|[}]","",unlist(strsplit(Files,"} {",fixed=TRUE)))
		} else {
			Files <- Folder
			# check existing
			fExist <- file.exists(Files)
			if(!all(fExist)){
				Files <- Folder[fExist]
				warning("The following specified raw data files do not exist:",paste("\n\t->",Folder[!fExist]))
			}
		}
	} else {
		if(any(From_null,To_null)){
			stop("Please specify either both arguments 'From' and 'To' or none.")
		}
		From <- parse_timerange(From,tz=tz)
		From_date <- trunc(From + time_add, 'days')
		To <- parse_timerange(To,tz=tz)
		To_date <- trunc(To + time_add, 'days')
		Year_From <- year(From + time_add)
		Year_To <- year(To + time_add)
		Years <- as.character(unique(c(Year_From,Year_To)))
		if(basename(Folder) %in% Years){
			if(length(Years) > 1){
				warning("Meas. period spans several years, however raw data folder for year '",basename(Folder),"' has been specified.")
			}
			Folder <- dirname(Folder)
		}
		Files <- character(0)
		# browser()
		for(Yr in Years){
			Folder_Yr <- paste(Folder,Yr,sep="/")
			if(!file.exists(Folder_Yr)){
				warning("Directory '",Folder_Yr,"' does not exist!\nReading files from directory '",Folder,"'")
				Folder_Yr <- Folder
			}
			Dir_Files <- dir(Folder_Yr,pattern=".LOG")
			Dir_Files_Date <- parse_date_time(Dir_Files,orders = "Ymd",tz=tz)
			Files <- c(Files,paste(Folder_Yr,Dir_Files[Dir_Files_Date>=From_date & Dir_Files_Date<=To_date],sep="/"))
		}
	}

	### use fread from data.table
	DataList <- vector("list",length(Files))
	for(i in seq_along(Files)){
		DataList[[i]] <- fread(cmd = paste0("grep -e '$GFDTA' '",path.expand(Files[i]),"'"),fill=TRUE,blank.lines.skip = TRUE,header=FALSE)
	}
	Data <- na.omit(rbindlist(DataList), cols = c("V2","V6"))
	rm(DataList)

	Data[,V8 := gsub("[*].*$","",V8)][,V1 := NULL]
	setnames(Data,c("CH4 (ppm-m)","R2","distance","received_power","date_time","ser_number","status_code"))

	dateTime <- Data[,parse_date_time(date_time,orders="%Y/%m/%d %H:%M:%S",tz=tz)] + time_add
	st_dateTime <- dateTime[1]
	TimeDiff <- diff(dateTime)
	# correct negative time differences!
	if(any(TimeDiff < 0)){
		TD <- which(TimeDiff < 0)
		for(i in TD){
			warning("Internal clock changed back in time -> Negative time differences occured.\n\t-> Shifting time after ",format(dateTime[i])," by ",-TimeDiff[i] + 1," seconds!")
			dateTime[(i + 1) : length(dateTime)] <- dateTime[(i + 1) : length(dateTime)] - TimeDiff[i] + 1
		}
		TimeDiff <- diff(dateTime)
	} 
	dT0 <- which(TimeDiff==0)
	if(length(dT0)){
		dT1 <- c(as.numeric(TimeDiff!=1),1)
		cs <- mapply(function(x,y,z,dt1,dtT){
			xx <- which(cumsum(dt1[x:y])!=0)[1] - 1
			if(xx==0){
				out <- dtT[z] + 0.5
				names(out) <- z + 1
			} else {
				out <- dtT[z] + (1:xx)*xx/(xx+1)
				names(out) <- z + 1:xx
			}
			out
		},x=dT0+1,y=c(dT0[-1],length(dateTime)),z=dT0,MoreArgs=list(dt1=dT1,dtT=as.numeric(dateTime - st_dateTime,units="secs")))

		csu <- unlist(cs)
		ind <- as.numeric(names(csu))
		dateTime[ind] <- st_dateTime + csu
		TimeDiff <- diff(dateTime) 
	}
	takeme <- c(TRUE,TimeDiff!=0)
	stTime <- dateTime[takeme]
	# browser()
	nms <- names(Data)[!(names(Data) %in% "date_time")]
	Out <- as.ibts(Data[(takeme),nms,with=FALSE],st=stTime,et=stTime + c(pmin(TimeDiff[TimeDiff!=0],3),median(TimeDiff[TimeDiff!=0])))
	if(!asCharacter){
		for(i in 1:4)Out[,i] <- as.numeric(Out[,i])
		colClasses(Out)[1:4] <- c("avg","min","avg","min")
	}

  if(!any(From_null,To_null)){
    Out <- Out[deparse_timerange(From,To)]
  }

	if(!DTA_only){
		### use fread from data.table
		DataList <- vector("list",length(Files))
		for(i in seq_along(Files)){
			DataList[[i]] <- fread(file = Files[i],fill=TRUE,blank.lines.skip = TRUE,sep="",header=FALSE)
		}
		Data <- rbindlist(DataList)
		rm(DataList)

		#### LCA
    LCA <- fread(text=unlist(Data[grepl("$GFLCA",V1,fixed = TRUE)]),header=FALSE)
		LCA[,V18 := gsub("[*].*$","",V18)][,V1 := NULL]
		setnames(LCA,c("LinReg_Enabled","LD_Status","Inner_TEC_Therm","Local_TEC_Therm","Outer_TEC_Therm","Inner_TEC_SetPt","Ramp_DC_Offset",
			"Curr_Peak","Targ_Peak","Ref_Gain","Path_Gain","Ref_Qual","Ref_R2","Phase_Avg","Ref_Light","Status","date_time"))
		LCA[,date_time := parse_date_time(date_time,orders="%Y/%m/%d %H:%M:%S",tz=tz) + time_add]
		#### DBG
    DBG <- fread(text=unlist(Data[grepl("$GFDBG",V1,fixed = TRUE)]),header=FALSE)
		DBG[,V9 := gsub("[*].*$","",V9)][,V1 := NULL]
		setnames(DBG,c("Internal_temperature","offset","ref_cell_quality","ref_R2","ref_status","supply_voltage","system_status","date_time"))
		DBG[,date_time := parse_date_time(date_time,orders="%Y/%m/%d %H:%M:%S",tz=tz) + time_add]
		# stop("To Do!")
		Out <- structure(Out,LCA=LCA,DBG=DBG)
	}

  Out
}


################################## Calibration of TDL analyzers
correct.pressure <- function(P,a,b){
	a / (
		4 * P ^ 2 / b - 
		sqrt(2) * P ^ 2 / b * 
		((b / P ^ 2 + 2) * sqrt(2 + 2 * b / P ^ 2) / (1 + b / P ^ 2))
		)
}
correct.temperature <- function(T,a0,a1,a2,a3,a4,a5,a6,a7){
	a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4 + a5*T^5 + a6*T^6 + a7*T^7
}

# #### CH4OP-1003
# P_1003 <- data.frame(p=(7:12)*10,C=c(0.5865,0.7015,0.8325,0.9794,1.1419,1.3202))
# T_1003 <- data.frame(T=seq(-50,150,10),C=c(0.9693,0.9704,0.9728,0.9763,0.9809,0.9865,0.9928,1.0000,1.0079,1.0165,1.0257,1.0355,1.0459,1.0568,1.0683,1.0802,1.0927,1.1056,1.1189,1.1327,1.1470))

# #### CH4OP-1049
# P_1049 <- data.frame(p=(7:12)*10,C=c(0.9933,0.9806,0.9832,0.9974,1.0210,1.0526))
# T_1049 <- data.frame(T=seq(-50,150,10),C=c(0.7137,0.7488,0.7859,0.8249,0.8658,0.9086,0.9534,1.0000,1.0486,1.0991,1.1515,1.2060,1.2624,1.3208,1.3813,1.4438,1.5084,1.5751,1.6438,1.7148,1.7878))

# #### CH4OP-30015
# P_30015 <- data.frame(p=(7:12)*10,C=c(0.9133,0.9274,0.9551,0.9941,1.0426,1.0996))
# T_30015 <- data.frame(T=seq(-50,150,10),C=c(0.7580,0.7875,0.8188,0.8517,0.8863,0.9226,0.9605,1.0000,1.0411,1.0839,1.1283,1.1743,1.2220,1.2713,1.3223,1.3749,1.4292,1.4852,1.5429,1.6023,1.6635))

# #### CH4OP-30016
# P_30016 <- data.frame(p=(7:12)*10,C=c(0.8573,0.8896,0.9351,0.9917,1.0584,1.1342))
# T_30016 <- data.frame(T=seq(-50,150,10),C=c(0.7909,0.8161,0.8430,0.8714,0.9013,0.9328,0.9657,1,1.0358,1.0730,1.1116,1.1516,1.1931,1.2359,1.2802,1.3259,1.3730,1.4216,1.4716,1.5230,1.5759))

# #### CH4OP-30017
# P_30017 <- data.frame(p=(7:12)*10,C=c(0.9133,0.9274,0.9551,0.9941,1.0426,1.0996))
# T_30017 <- data.frame(T=seq(-50,150,10),C=c(0.7580,0.7875,0.8188,0.8517,0.8863,0.9226,0.9605,1,1.0411,1.0839,1.1283,1.1743,1.2220,1.2713,1.3223,1.3749,1.4292,1.4852,1.5429,1.6023,1.6635))

# #### CH4OP-30018
# P_30018 <- data.frame(p=(7:12)*10,C=c(0.9261,0.9359,0.9597,0.9946,1.0391,1.0919))
# T_30018 <- data.frame(T=seq(-50,150,10),C=c(0.7507,0.7812,0.8134,0.8473,0.8830,0.9203,0.9593,1,1.0423,1.0864,1.1321,1.1794,1.2285,1.2793,1.3318,1.3860,1.4419,1.4996,1.5591,1.6203,1.6834))

# #### CH4OP-30025
# P_30025 <- data.frame(p=(7:12)*10,C=c(0.8863,0.9092,0.9455,0.9930,1.0501,1.1161))
# T_30025 <- data.frame(T=seq(-50,150,10),C=c(0.7737,0.8012,0.8303,0.8611,0.8935,0.9274,0.9630,1.0000,1.0386,1.0787,1.1203,1.1634,1.2081,1.2543,1.3021,1.3514,1.4022,1.4546,1.5086,1.5642,1.6213))

# #### CH4OP-30026
# P_30026 <- data.frame(p=(7:12)*10,C=c(0.9261,0.9359,0.9597,0.9946,1.0391,1.0919))
# T_30026 <- data.frame(T=seq(-50,150,10),C=c(0.7507,0.7812,0.8134,0.8473,0.8830,0.9203,0.9593,1.0000,1.0423,1.0864,1.1321,1.1794,1.2285,1.2793,1.3318,1.3860,1.4419,1.4996,1.5591,1.6203,1.6834))

if (FALSE) {
    gf.calibration = list(
        "CH4OP.1003" = list(p = c(-0.09970423838,2567.14526277485),T = c(9.864418e-01,5.977530e-04,4.289765e-06,-1.373450e-08,5.874422e-11,6.124736e-14,-2.473273e-15,8.490556e-18))
        ,"CH4OP.1049" = list(p = c(-0.3362784861,33263.7220313790),T = c(9.086332e-01,4.378029e-03,9.536197e-06,3.367833e-10,2.385533e-11,-6.367599e-14,-1.518728e-16,7.363739e-19))
        ,"CH4OP.30015" = list(p = c(-0.3130739872,21583.3936239178),T = c(9.225870983e-01,3.707544541e-03,8.155413500e-06,-1.531357343e-09,2.666445763e-11,-1.327184828e-13,8.755430303e-17,1.190407513e-18))
        ,"CH4OP.30016" = list(p = c(-0.2886017923,16039.7106450116),T = c(9.327787473e-01,3.217314336e-03,7.313661791e-06,-3.882584765e-09,2.951756048e-11,-8.829885883e-14,4.191947142e-17,2.047602667e-19))
        ,"CH4OP.30017" = list(p = c(-0.3130739872,21583.3936239178),T = c(9.225870983e-01,3.707544541e-03,8.155413500e-06,-1.531357343e-09,2.666445763e-11,-1.327184828e-13,8.755430303e-17,1.190407513e-18))
        ,"CH4OP.30018" = list(p = c(-0.3177601163,23100.0009584746),T = c(9.202979911e-01,3.817143728e-03,8.401503593e-06,-1.683137404e-09,1.661059500e-11,9.442030054e-14,-1.414050661e-15,4.508541276e-18))
        ,"CH4OP.30025" = list(p = c(-0.3020502714,18707.3106621851),T = c(9.274455217e-01,3.474331304e-03,7.764144758e-06,-3.589742775e-09,2.733007984e-11,7.828519152e-14,-1.551044182e-15,4.698039907e-18))
        ,"CH4OP.30026" = list(p = c(-0.3177601163,23100.0009584746),T = c(9.202979911e-01,3.817143728e-03,8.401503593e-06,-1.683137404e-09,1.661059500e-11,9.442030054e-14,-1.414050661e-15,4.508541276e-18))
        )
    save(gf.calibration, file = 'data/gf.calibration.rds')
}

correct_pT <- function(instrument,p_hPa,T_deg,p_hPa.inst=1013.25,T_deg.inst=20,units=c("mgm3","ppm")){
	# get instrument's calibration
	CalPars <- gel::gf.calibration[[switch(which(sapply(c("1003","1049","15","16","17","18","25","26"),function(x)grepl(x,instrument)))
		,"CH4OP.1003"
		,"CH4OP.1049"
		,"CH4OP.30015"
		,"CH4OP.30016"
		,"CH4OP.30017"
		,"CH4OP.30018"
		,"CH4OP.30025"
		,"CH4OP.30026"
		,stop("Could not relate argument 'instrument' to any available instrument ('1003','1049','30015','30016','30017','30018','30025','30026')")
		)]]
	# convert raw data to ppm @ ambient T and p
	C_out <- correct.pressure(p_hPa/10,CalPars$p[1],CalPars$p[2])/
		correct.pressure(p_hPa.inst/10,CalPars$p[1],CalPars$p[2])*
		correct.temperature(T_deg,CalPars$T[1],CalPars$T[2],CalPars$T[3],CalPars$T[4],CalPars$T[5],CalPars$T[6],CalPars$T[7],CalPars$T[8])/
		correct.temperature(T_deg.inst,CalPars$T[1],CalPars$T[2],CalPars$T[3],CalPars$T[4],CalPars$T[5],CalPars$T[6],CalPars$T[7],CalPars$T[8])
	# convert ppm to chosen unit
	switch(pmatch(units[1],c("ppm","ppb","mgm3","ugm3"))
		# "ppm"
		, C_out
		# "ppb"
		, C_out*1000
		# "mgm3"
		, C_out*16.04*p_hPa/10/(8.3144598*(T_deg + 273.15))
		# "ugm3"
		, C_out*16.04*p_hPa*100/(8.3144598*(T_deg + 273.15))
		# else
		, stop("Argument 'units' was not recognized. Use ('ppm','ppb','mgm3','ugm3').")
		)
}

# wrapper function for ppm
gf2ppm <- function(data, 
    single_path,
    p_hPa = '(p|P)[_.-]?(air|Air|hpa|hPa)',
    t_deg = '(t|T)[_.-]?(air|Air|deg|Deg)',
    p_hPa.inst = 1013.25, 
    T_deg.inst = 20, 
    instrument = data[1, 'ser_number']
    ) {
    # check p_hPa
    if (is.character(p_hPa)) {
        p_name <- grep(p_hPa, names(data), value = TRUE)[1]
        p_hPa <- data[[p_name]]
    }
    # check t_dep
    if (is.character(t_deg)) {
        t_name <- grep(t_deg, names(data), value = TRUE)[1]
        t_deg <- data[[t_name]]
    }
    # convert
    data[, 'CH4 (ppm-m)'] * 
        correct_pT(instrument, p_hPa, t_deg, p_hPa.inst, T_deg.inst, units = 'ppm') /
        single_path
}

# wrapper function for mg/m3
gf2mgm3 <- function(data, 
    single_path,
    p_hPa = '(p|P)[_.-]?(air|Air|hpa|hPa)',
    t_deg = '(t|T)[_.-]?(air|Air|deg|Deg)',
    p_hPa.inst = 1013.25, 
    T_deg.inst = 20, 
    instrument = data[1, 'ser_number']
    ) {
    # check p_hPa
    if (is.character(p_hPa)) {
        p_name <- grep(p_hPa, names(data), value = TRUE)[1]
        p_hPa <- data[[p_name]]
    }
    # check t_dep
    if (is.character(t_deg)) {
        t_name <- grep(t_deg, names(data), value = TRUE)[1]
        t_deg <- data[[t_name]]
    }
    # convert
    data[, 'CH4 (ppm-m)'] * 
        correct_pT(instrument, p_hPa, t_deg, p_hPa.inst, T_deg.inst, units = 'mgm3') /
        single_path
}

