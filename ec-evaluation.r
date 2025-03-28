## 13.04.2018 - Christoph - 
##     Funktion 'readWindMaster_ascii' abgeändert. Korrupte Files werden führen nun nicht mehr zu einem Fehler.
##     Funktion 'evalReddy': Argument 'create_graphs' hinzugefügt (FALSE -> kein plotten).
##                           ca. L814: zu kurze Messfiles werden abgefangen und führen nicht mehr zu einem Fehler.
# *************************************************************************************
## 16.01.2018 - Marcel - Funktion zur Berechnung der Pfadlänge GF hinzugefügt
# *************************************************************************************

# REddy.version <- "v20220218"
library(data.table)
library(ibts)
library(Rcpp)

# library(latticeExtra)
# library(gridExtra)
# library(grid)

# *************************************************************************************
# new formula database
# readWindMaster_ascii
# readART2012_binary
# readART2012_ascii
# readART2012_rds
# convertARTbinary_2012
# "%w/o%"
# H.flags
# trend
# rotate_detrend

# ************************************
# Reading raw data
# ************************************

# readART2012_binary
# readART2012_ascii
# readLicor_ascii (Jesper's Licor EC measurements)


# FilePath <- "~/Y-Drive/Jesper/Foulum_MiniDOAS/EC_raw/2019-08-15T200000_AIU-1273.ghg"
# FilePath <- "~/Y-Drive/Jesper/Foulum_MiniDOAS/EC_raw/2019-08-15T200000_AIU-1273.data"
read4hLicor <- function(FilePath){
	filename <- gsub("ghg$","data",basename(FilePath))
	con <- unz(FilePath, filename)
	read.table(con, skip = 7, header = TRUE, sep = "\t", check.names = FALSE)
}

# Folder <- "~/Y-Drive/Jesper/Foulum_MiniDOAS/EC_raw"
convert2dailyLicor <- function(Folder){
	ghg <- dir(Folder, pattern = ".ghg")
	dates <- gsub("([0-9]{4}[-][0-9]{2}[-][0-9]{2})T.*","\\1",ghg)
	times <- gsub("([0-9]{4}[-][0-9]{2}[-][0-9]{2})T([0-9]{6}).*","\\2",ghg)
	ending <- gsub("([0-9]{4}[-][0-9]{2}[-][0-9]{2})T([0-9]{6})(.*)[.]ghg","\\3",ghg)
	uend <- unique(ending)
	udates <- unique(dates)
	for(u in uend){
		for(d in udates){
			cat("Gathering data for",d,"\n")
			ghgList <- lapply(file.path(Folder,ghg[dates == d][order(times[dates == d])]), read4hLicor)
			cat("Writing daily file.\n")
			fwrite(rbindlist(ghgList), file = file.path(Folder,paste0(d,u,".txt")))
		}
	}
	invisible(NULL)
}
# convert2dai1lyLicor(Folder)

# FilePath <- "~/Y-Drive/Jesper/Foulum_MiniDOAS/EC_raw/2019-08-15_AIU-1273.txt"
readLicor_ascii <- function(FilePath){
	Data <- fread(FilePath)
	Data[, st := fast_strptime(paste(Date[1],gsub("[:]([0-9]{3})$",".\\1",
		Time[1])),"%Y-%m-%d %H:%M:%OS", lt = FALSE) - Seconds[1] - Nanoseconds[1]*1E-9 + Seconds + Nanoseconds*1E-9]
	### check delta-t
	Data[,dt := c(round(as.numeric(diff(st))*1000),100)]
	Data <- Data[dt>50]
 	Data[dt > 150, dt := 100]
	### set Output names and order
	setnames(Data,
		c("Aux 1 - U (m/s)","Aux 2 - V (m/s)","Aux 3 - W (m/s)","Aux 4 - Ts (K)","dt"), 
		c("u_analog", "v_analog", "w_analog", "T_analog","delta-t"))
	# Conversion from aux (analog) to m/s:
	# u = u_analog * 24 - 60;
	# v = v_analog * 24 - 60;
	# w = w_analog * 2 - 5;
	# T = T_analog * 20 - 40;
	Data[,":="(
		u = u_analog * 24 - 60,
		v = v_analog * 24 - 60,
		w = w_analog * 2 - 5,
		T = T_analog * 20 - 40 + 273.15
		)]
	nms <- names(Data)
	setcolorder(Data,c("st","delta-t","u","v","w","T",nms %w/o% c("st","delta-t","u","v","w","T")))

	as.data.frame(Data)
}

readWindMaster_ascii <- function(FilePath,skip=","){
	### get Date
	bn <- basename(FilePath)
	Date <- gsub(".*_([0-9]{6})_.*","\\1",bn)
	### read File
	suppressWarnings(out <- try(fread(FilePath,sep=",",fill=TRUE,blank.lines.skip = TRUE),silent=TRUE))
	if(inherits(out,"try-error")){
		cat("\nRaw data file is corrupt. Trying to read it anyway.\n")
		msg <- conditionMessage(attr(out,"condition"))
		lin <- as.numeric(gsub(".*line ([0-9]*) contains.*","\\1",msg))
		out1 <- fread(FilePath,sep=",",fill=TRUE,nrows=lin-1)
		out2 <- fread(FilePath,sep=",",fill=TRUE,skip=lin-1)
		suppressWarnings(out2[,V2 := as.character(V2)])
		suppressWarnings(out2[,V3 := as.numeric(V3)])
		suppressWarnings(out2[,V4 := as.numeric(V4)])
		suppressWarnings(out2[,V5 := as.numeric(V5)])
		suppressWarnings(out2[,V6 := as.numeric(V6)])
		suppressWarnings(out2[,V7 := as.numeric(V7)])
		out2 <- na.omit(out2,cols=2:7)
		out2 <- out2[!(V2 %chin% ""),1:7]
		out <- rbind(out1,out2)
	}
 	out[,V1 := NULL]
 	out <- na.omit(out)
	### set times
	out[,st :=fast_strptime(paste0(Date,V2),lt=FALSE,format="%y%m%d%H:%M:%S",tz="Etc/GMT-1")+V3][,c("V2","V3"):=NULL]
	### check delta-t
	out[,dt := c(round(as.numeric(diff(st))*1000),100)]
	# browser()
	out[min(which(dt < -10000),.N):.N,":="(
		st = st[1]
		,dt = 100
		)]
	out <- out[dt>50]
 	out[dt > 150, dt := 100]
	### set Output names and order
	setnames(out,c("u", "v", "w", "T", "st","delta-t"))
	setcolorder(out,c("st","delta-t","u","v","w","T"))
	### remove 999.99 entries
	out <- out[!(u%in%999.99|v%in%999.99|w%in%999.99|T%in%999.99),]
	### change units from °C to K
	out[,T := T + 273.15]
	as.data.frame(out)
}

read_ht_merged <- function(FilePath) {
    out <- fread(file = FilePath, showProgress = FALSE)
	### check delta-t
	out[,dt := c(round(as.numeric(diff(Time))*1000),100)]
    dt_med <- out[, median(dt)]
 	out[dt > dt_med * 1.5, dt := dt_med]
    as.data.frame(out[dt > dt_med / 2, 
        .(st = with_tz(Time, 'Etc/GMT-1'), dt, u, v, w, T, 
            nh3_ugm3, nh3_ppb, temp_amb, press_amb, oss)
        ])
}


readART2012_rds <- function(FilePath){
	Data <- readRDS(FilePath)
	Start <- parse_date_time2(paste(Data[1,1:2],collapse=" "),"%d.%m.%Y %H:%M:%S",tz="Etc/GMT-1")
	cn <- as.character(Data[2,])
	Data <- data.frame(lapply(Data[-(1:2),],as.double))
	names(Data) <- cn
	# delta-t kommt von LabVIEW -> Buffer Problem!
	# csum <- cumsum(Data[,1])/1000
	out <- cbind(
		st=Start + seq.int(0,(nrow(Data) - 1)/10,0.1)
		,Data
		)
	out
}

convertARTbinary_2012 <- function(Pfad,FilePattern="PSX%y%m%d_0.dat",asRDS=TRUE,asCSV=FALSE){	

	if(file.info(Pfad)$isdir){
		Files <- choose.files(default=dir(Pfad,full.names=TRUE)[1])
	} else {
		Files <- Pfad
	}

	for(File in Files){
		### Read Header:
		header.length <- read.table(File, header=FALSE, sep=";", stringsAsFactors=FALSE, nrows=1, flush=TRUE) # read only the first row to extract the header length (recorded in the first row as ascii)
		header.compl <- read.table(File, header=FALSE, sep=";", stringsAsFactors=FALSE, row.names=1,       # read entire header, except 1st line
				skip=1, nrows=as.numeric(header.length[2])-1, flush=TRUE, fill=TRUE)

		running.bytes <- as.numeric(header.compl[3,])
	  
		channels <- length(running.bytes)                                                              # number of data channels
		total.bytes <- sum(running.bytes)
		signed <- as.numeric(header.compl[5,]) == 1

		##### read data
		con.binary <- file(File,open="rb")
		on.exit(close(con.binary))
		Nirvana <- readBin(con.binary,"raw",n=as.numeric(header.length[1]),size=1)
		rm(Nirvana)
		File_Seconds <- as.numeric(fast_strptime(paste(header.compl[1,3], header.compl[1,4]),format="%d.%m.%Y %H:%M:%S", tz="UTC",lt=FALSE) - fast_strptime(paste(header.compl[1,1], header.compl[1,2]),format="%d.%m.%Y %H:%M:%S", tz="UTC",lt=FALSE),units="secs")
		if(is.na(File_Seconds)) File_Seconds <- 24*3600
		max.datapoints.estimate <- 11*File_Seconds
		output <- matrix(nrow=max.datapoints.estimate, ncol=channels)
		KeepOnReading <- TRUE
		starts <- cumsum(c(1,running.bytes[-length(running.bytes)]))
		ends <- cumsum(running.bytes)
		max.running.bytes <- max(running.bytes)
		running.factor <- 2^(((max.running.bytes-1):0)*8)
		j <- 0
		cat("Reading binary data...\n")
		pb <- tcltk::tkProgressBar("Reading binary data...","progress",0,max.datapoints.estimate,0)
		on.exit(close(pb),add=TRUE)
		while(KeepOnReading){
			dat.raw <- readBin(con.binary, "raw", n=total.bytes)
			if(length(dat.raw)==0){
				KeepOnReading <- FALSE
			} else {
				j <- j + 1
				if(j==max.datapoints.estimate){
					cat("data file has more data than estimated. ")
					cat("extending output matrix...\n")
					close(pb)
					pb <- tcltk::tkProgressBar("Reading binary data...","progress",max.datapoints.estimate,max.datapoints.estimate + File_Seconds,max.datapoints.estimate)
					max.datapoints.estimate <- max.datapoints.estimate + File_Seconds
					output <- rbind(output,matrix(nrow=File_Seconds, ncol=channels))
				}
				# cat("\r",j)
				tcltk::setTkProgressBar(pb,j)
				for (i in 1:channels) {                                                                    # read the number of bytes for one complete "row" of the dataset                                 
					if (signed[i] & as.integer(dat.raw)[starts[i]] >= 128) {# consider sign-byte (readBin() itself can only consider signed bytes in 4-byte or 2-byte system...)
						output[j,i] <- -sum(running.factor[max.running.bytes + 1 - (running.bytes[i]:1)] * as.integer(!dat.raw[starts[i]:ends[i]])) + 1														
					} else {
						output[j,i] <- sum(running.factor[max.running.bytes + 1 - (running.bytes[i]:1)] * as.integer(dat.raw[starts[i]:ends[i]]))
					}
				}
			}
		}
		close(con.binary)
		close(pb)
		on.exit()
		cat("done.\n")
		output <- output[1:j,]

		output <- sweep(output,2,as.numeric(header.compl[4,]),"/")
		colnames(output) <- names(header.compl)
		mode(output) <- "character"
		out <- rbind(header.compl[1:2,],output)
		if(asCSV)write.table(out,file=sub(".dat$",".csv",File),sep=";",row.names=FALSE,col.names=FALSE)
		if(asRDS)saveRDS(out,file=sub(".dat$",".rds",File))
	}
	invisible(NULL)
}


# misc
"%w/o%" <- function(x,y){
	x[!(x %in% y)]
}


# hard flag function
H.flags <- function(input, time, d_t, limits, wind = 500, hflg.met = "norepl"){ 
	# *****************************************************************************
	# ***************** Detects values outside of physical range... ***************
	# ********************** replaces them by running mean ************************ 
	# ********************* author: Raphael Felber, 16.08.2012 ********************
	# ********************* adpted: Raphael Felber, 25.02.2014 ********************
	# *****************************************************************************

	nms <- names(input)
	missing <- !(nms %in% colnames(limits))

	if(length(nms[missing]) > 0) stop(paste0("physical limits missing for: ", paste(nms[missing], collapse=", ")), call. = TRUE, domain = NULL)
	limits <- limits[,nms]
	
	### replace values outside of limits by NAs
	dat <- mapply(function(x,lo,hi){
		x[x < lo | x > hi] <- NA
		x
	},x=input,lo=limits["lower",],hi=limits["top",],SIMPLIFY = FALSE)

	hflgs <- unlist(lapply(dat,anyNA))
	if(any(hflgs))cat("Hard Flag: flagged values in columns ",paste(names(hflgs)[hflgs],sep=", "),"\n")
	
	if(pmatch(hflg.met, "replace", nomatch = 0) && any(hflgs)){

		# cat("Replacing NA values by window mean...\n")
		cat("Replacing NA values by window median...\n")
		# since we're only intrested in specific time windows, find first NAs:
		isna <- lapply(dat[hflgs],function(x)which(is.na(x)))
		# replace if all NA? (Why could this happen?)
		l <- length(dat[hflgs][[1]])
		replAll <- (l - lengths(isna)) < 2
		if(any(replAll)) dat[hflgs][replAll] <- rep(-99999,l)
		
		# replace wind by seconds
		wind <- parse_time_diff(wind)

		### create matrix with running means of size wind (default = 500)
		# original code
		# run.m <- lapply(dat_r,function(x)rollapply(x, wind, mean, by.column=TRUE, na.rm = TRUE, fill="extend"))
		# using much! faster function caTools::runmean
		# run.m <- lapply(dat_r,function(x)caTools::runmean(x, wind, alg = "C", endrule = "mean"))
		# using correct windows, bit slower though
		st1 <- as.numeric(time)
		st2 <- c(st1[-1], st1[length(st1)] + d_t / 1000)
		for(i in seq_along(isna)){
			x1 <- st1[isna[[i]]] - wind/2 + d_t/2000
			x2 <- st1[isna[[i]]] + wind/2 - d_t/2000
            # NOTE: fixme cutIntervals & getIntervals sind continuous!!! -> add contin. check in ibts!!!
            ind <- find_window(time, x1, x2)
			# dat[hflgs][[i]][isna[[i]]] <- sapply(ind, function(x) mean(dat[hflgs][[i]][x], na.rm = TRUE))
			dat[hflgs][[i]][isna[[i]]] <- sapply(ind, function(x) median(dat[hflgs][[i]][x], na.rm = TRUE))
		}	
		cat("number of replaced values\n*~~~~*\n", names(dat[hflgs]),"\n", lengths(isna),"\n*~~~~*\n")
	}
							
	dat
	
}

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List find_window(NumericVector x, NumericVector y1, NumericVector y2)
	{
	Rcpp::List Out(y1);
	int lenx = x.size();
	int leny = y1.size();
	int run = 0;
	int j = 0;
	for (int i = 0; i < leny; i++) {
        // pass j to runner & set j to 0
        run = j;
        j = 0;
        LogicalVector LogVec = rep(false, lenx);
		if ((x[lenx - 1] < y1[i]) || (x[run] > y2[i])) {
			Out[i] = LogVec;
		} else {
            // goto y1[i]
			while ((x[run] < y1[i]) && (run < lenx)) {
				run += 1;
			}
            // goto y2[i]
			while ((run < lenx) && (x[run] < y2[i])) {
                if (i < leny && j == 0 && x[run] >= y1[i + 1]) {
                    j = run;
                }
                LogVec[run] = true;
				run += 1; 
			}
            // check j
            if (j == 0) {
                j = run;
            } else if (j == lenx) {
                j -= 1;
            }
            Out[i] = LogVec;
		}
	}
	return(Out);
}
')

trend <- function(y,method=c("blockAVG","linear","linear_robust","ma_360"),Hz_ts=10){
	n <- length(y)
	method <- method[1]
	switch(method,
		"blockAVG"={
			my <- .Internal(mean(y))
			fitted <- rep(my,n) 
			list(
				coefficients=c(intercept=my,slope=0)
				,fitted=fitted
				,residuals=y - fitted
				) 
		}
		,"linear"={
			my <- .Internal(mean(y))
			mx <- (n+1)/2
			xstr <- (x <- seq.int(n)) - mx
			b <- sum(xstr*(y-my))/(n-1)/(sum(xstr^2)/(n-1))
			a <- my - mx*b
			fitted <- a + b*x
			list(
				coefficients=c(intercept=a,slope=b)
				,fitted=fitted
				,residuals=y - fitted
				) 
		}
		,"linear_robust"={
			mod <- MASS::rlm(y~seq.int(n))
			cfs <- mod$coefficients
			a <- cfs[1]
			b <- cfs[2]
			fitted <- mod$fitted
			if(!mod$converged){
				stop("robust linear regression did not converge!")
			}
			list(
				coefficients=c(intercept=a,slope=b)
				,fitted=fitted
				,residuals=y - fitted
				) 
		}
		,{
            if (!grepl('^ma_', method)) stop('detrending method not valid! ',
                'Should be one of "blockAVG", "linear", "linear_robust" or "ma_xxx" where',
                ' xxx represents the moving average seconds')
			ma <- round(as.numeric(gsub("ma_","",method))*Hz_ts)
			fitted <- caTools::runmean(y,ma,"C","mean")
			list(
				coefficients=c(intercept=NA,slope=NA)
				,fitted=fitted
				,residuals=y - fitted
				) 
		}
	)
}


rotate_detrend <- function(u, v, w, T, phi = NULL, 
    method = c(u = "blockAVG", v = "blockAVG", w = "blockAVG", T = "linear"), 
    c.system = "Windmaster", Hz_ts = 10) {
	thetam <- atan2(.Internal(mean(v)), .Internal(mean(u)))
	u1 <- u * cos(thetam) + v * sin(thetam)
	if (is.null(phi)) {
		phi <- atan2(.Internal(mean(w)), .Internal(mean(u1)))
	}
	ud <- trend(u1 * cos(phi) + w * sin(phi), method["u"], Hz_ts)
	vd <- trend(-u * sin(thetam) + v * cos(thetam), method["v"], Hz_ts)
	wd <- trend(-u1 * sin(phi) + w * cos(phi), method["w"], Hz_ts)		
	Td <- trend(T, method["T"], Hz_ts)
	out <- list(
		uprot = ud$residuals
		, vprot = vd$residuals
		, wprot = wd$residuals
		, Tdet = Td$residuals
		, wd = switch(c.system
			, "Windmaster" = (180 - thetam / pi * 180) %% 360 
            , "HS-Wauwilermoos" = (150 - thetam / pi * 180) %% 360
			, "ART.EC1" = ((180 / pi) * -thetam + 150 + 147) %% 360
			# new (coordinate) system
			, "Licor" = (270 - thetam / pi * 180) %% 360
			# # old (coordinate) system
			# ,"Licor" = (180 + thetam / pi * 180) %% 360 
			)
		, phi = phi
		, umrot = ud$fitted
		, vmrot = vd$fitted
		, wmrot = wd$fitted
		, Tmdet = Td$fitted
		)	
	return(out)
}

cov_fun <- function(x1,x2){
	n <- length(x1)
	x1 <- fft(x1)/n
	x2 <- fft(x2)/n
	if(n%%2){
		# n ungerade:
		Re(fft(Conj(x2) * x1, inverse=TRUE))[c(((n+1)/2+1):n,1:((n+1)/2))]*n/(n-1)
	} else {
		# n gerade:
		Re(fft(Conj(x2) * x1, inverse=TRUE))[c((n/2+1):n,1:(n/2))]*n/(n-1)
	}
}

cospec_fun <- function(x1,x2,lag=0){
	n <- length(x1)
	x1 <- fft(x1)/n
	x2 <- fft(shift(x2,-lag))/n
	Re(Conj(x2) * x1)[seq(n/2)+1]*n/(n-1)*2
}


shift <- function(x, lag){
	if (lag != 0) {
		if (lag > 0) {
			indr <- seq.int(lag)
			indl <- seq.int(length(x) - lag) + lag
			x <- x[c(indl,indr)]
		} else {
			indr <- seq.int(length(x) + lag)
			indl <- length(x) - seq.int(-lag,1) + 1
			x <- x[c(indl,indr)]
		}
	}
	x
}

find_dynlag <- function(x,dyn){
	n <- length(x)
	if(n%%2){
		# ungerade
		m <- (n+1)/2
	} else {
		# gerade
		m <- n/2 + 1
	}
	ind <- seq(dyn[1],dyn[2]) + m
	# find max:
	maxis <- ind[which.max(abs(x[ind]))]
	c(index = maxis, tau = maxis - m)
}

reduce_cospec <- function(cospec,freq,length.out=100){
	log_freq <- log(freq)
	log_cuts <- seq(min(log_freq),max(log_freq),length.out=length.out+1)
	ind_cuts <- findInterval(log_freq,log_cuts,rightmost.closed=TRUE)
	freq_out <- exp(log_cuts[-1] - diff(log_cuts) / 2)
	list(
        cospec = tapply(cospec * freq,ind_cuts,mean),
        freq = freq_out[unique(ind_cuts)]
    )
}

damp_hac5 <- function(ogive, freq, freq.limits, ogive_ref){
	require(deming)
	# select undamped frequqency region:
	undamped_ind <- freq < (1/freq.limits[1]) & freq > (1/freq.limits[2])
	undamped_ogv <- ogive[undamped_ind]
	undamped_ogvref <- ogive_ref[undamped_ind]

	if(sign(sum(undamped_ogvref))!=sign(sum(undamped_ogv))){
		undamped_ogvref <- -undamped_ogvref
		ogive_ref <- -ogive_ref
	}	
	# linear regression + prediction:
	# mod <- deming::deming(undamped_ogv~undamped_ogvref,weights=1/freq[undamped_ind])
	mod <- try(deming::deming(undamped_ogv ~ undamped_ogvref), silent = TRUE)
    if (inherits(mod, 'try-error')) {
        cfs <- c(0, 1)
    } else {
        cfs <- coef(mod)
    }
	pred_ogv_deming <- cfs[2]*ogive_ref + cfs[1]

	# robust linear regression + prediction:
	mod2 <- try(deming::pbreg(undamped_ogv ~ undamped_ogvref,method=3,eps=min(abs(undamped_ogv))*1E-8), silent = TRUE)
    if (inherits(mod2, 'try-error')) {
        cfs2 <- c(0, 1)
    } else {
        cfs2 <- coef(mod2)
    }
	pred_ogv_pbreg <- cfs2[2]*ogive_ref + cfs2[1]
	
	list(dampf_pbreg=ogive[1]/(ogive[1] - cfs2[1]),dampf_deming=ogive[1]/(ogive[1] - cfs[1]),freq.limits=freq.limits,ogive=ogive,ogive_ref_pbreg=pred_ogv_pbreg,ogive_ref_deming=pred_ogv_deming)	
}

# -cfs2 = og1/dp - og1 = (1/dp - 1)*og1


split_index <- function(index_length,n_subint){
	nm <- index_length/n_subint
	fn <- floor(nm)
	res <- index_length - n_subint*fn 
	len <- rep(fn,n_subint) + c(rep(1,res),rep(0,n_subint-res))
	mapply(function(x,y)seq.int(x)+y,x=len,y=cumsum(c(0,len[-length(len)])),SIMPLIFY=FALSE)
}

wind_statistics <- function(wind,z_canopy,z_sonic){
	Cov_sonic <- cov(list2DF(wind[1:4]))
	Var_sonic <- diag(Cov_sonic)
	names(Var_sonic) <- c('var_u', 'var_v', 'var_w', 'var_T')
	Cov_sonic <- Cov_sonic[cbind(c("uprot","uprot","uprot","vprot","vprot","wprot"),c("vprot","wprot","Tdet","wprot","Tdet","Tdet"))]
	names(Cov_sonic) <- c('cov_uv', 'cov_uw', 'cov_uT', 'cov_vw', 'cov_vT', 'cov_wT')
	suppressWarnings(Ustar <- c(sqrt(-Cov_sonic["cov_uw"]),use.names = FALSE))
	T_K <- mean(wind$Tmdet + wind$Tdet)
	U <- mean(wind$umrot + wind$uprot)
	L <- c(-Ustar ^ 3 * T_K / (0.4 * 9.80620 * Cov_sonic["cov_wT"]), use.names = FALSE)
	if (!is.na(z_canopy)) {
		d <- 2/3 * z_canopy
		suppressWarnings(z0 <- optimize(function(x, ustar, L, z, d, U) 
                abs(U - calcU(ustar, x, L, z - d)), c(0, z_sonic * 1.1), 
                ustar = Ustar, L = L, U = U, z = z_sonic, d = d)$minimum
        )
		if (z0 >= z_sonic * 1.09) z0 <- NA
	} else {
		d <- NA
		z0 <- NA
	}
	c(Var_sonic, Cov_sonic, Ustar = Ustar, L = L, z_sonic = z_sonic, z_canopy = z_canopy, 
        d = d, Zo = z0, WD = wind$wd, phi = wind$phi, U_sonic = U, T_sonic = T_K
    )
}


calcU <- function (ustar, Zo, L, z, kv = 0.4){
    zL <- z/L
    ZoL <- Zo/L
    psiMz <- ifelse(zL < 0, {
        x <- (1 - 16 * zL)^(1/4)
        -2 * log((1 + x)/2) - log((1 + x^2)/2) + 2 * atan(x) - 
            pi/2
    }, 4.8 * zL)
    psiMZo <- ifelse(zL < 0, {
        x <- (1 - 16 * ZoL)^(1/4)
        -2 * log((1 + x)/2) - log((1 + x^2)/2) + 2 * atan(x) - 
            pi/2
    }, 4.8 * ZoL)
    ustar/kv * (log(z/Zo) + psiMz - psiMZo)
}

plot.tseries <- function(dat,wind,scal,selection,color,units){
	msg <- paste(c(format(dat[1,1],"%d.%m.%Y")," - time series"),collapse="")
	tsbeginning <- dat[1,1]
	dat[,c("u","v","w","T")] <- mapply("+",wind[c("umrot","vmrot","wmrot","Tmdet")],wind[c("uprot","vprot","wprot","Tdet")])
	dat2 <- reshape2::melt(dat[,c("st",selection)],id="st")
	### get trends:
	dat3 <- list2DF(wind[c("umrot","vmrot","wmrot","Tmdet")])
	names(dat3) <- c("u","v","w","T")
	if(!is.null(scal)){
        # add scalars
		dat3 <- cbind(dat3, lapply(scal, \(x) {
                if (is.null(isna <- na.action(x$residuals))) {
                    x$fitted
                } else {
                    out <- rep(NA_real_, nrow(dat3))
                    out[-isna] <- x$fitted
                    out
                }
            })
        )
        # add residuals
        res <- setdiff(selection, names(dat3))
        if (length(res) > 0) {
            dat3 <- cbind(dat3, dat[, res])
        }
	}
	dat3 <- dat3[,selection]
	### melt and add trends:
	dat4 <- reshape2::melt(dat3, id = NULL, value.name = "trend")
	dat2[, "trend"] <- dat4[, "trend"]
	myxscale.component <- function(...) {
		ans <- xscale.components.default(...)
		ans$top <- ans$bottom
		ans
	}
	myyscale.component <- function(...) {
		ans <- yscale.components.default(...)
		ans$right <- ans$left
		ans$right$labels$labels <- NULL
		ans
	}
	ylab <- paste0(selection," (",units,")")
	ylab_r <- ylab_l <- rep(" ",length(selection))
	ylab_l[seq.int(1,length(selection),2)] <- ylab[seq.int(1,length(selection),2)]
	ylab_r[seq.int(2,length(selection),2)] <- ylab[seq.int(2,length(selection),2)]
	xyplot(value ~ st | variable, data=dat2, groups=variable, aspect=0.2, type="l",
		xlab=list("time of the day",cex=1.25), ylab.right=list(ylab_r,cex=1.25), ylab=list(ylab_l,cex=1.25), main=list(msg,cex=1.5),
		scales=list(
				x=list(cex=1.25, tck=c(-0.75,-0.75), format="%H:%M"),
				y=list(relation="free", cex=1.25, tck=c(-0.75,-0.75), rot=0)),
		xscale.component=myxscale.component, yscale.component=myyscale.component,
		strip=FALSE, layout=c(1, length(selection)), between=list(x=0,y=1), subscripts=TRUE, lwd=rep(1,length(color)), lty=rep(1,length(color)), col=color,
		panel=function(x, y, ...) {
			#panel.grid(h=-1, v=-1, lty=3, col="gray80")
			y2 <- dat2[list(...)$subscripts,"trend"]
			panel.xyplot(x,y,...)
			# panel.xyplot(x,y2,type="l",lwd=1.5, lty=3, col="gray30")
			# panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightblue")
			panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightgrey")
		}
	)  	
}

plot_covfunc <- function(cov_func,avg_t,dynLag,fixLag, ylab=NULL, xlim = NULL, cx=1.5, cxmt=1.25, cl="black"){
	midP <- dynLag[1] - dynLag[2]
	n <- length(cov_func)
	fix_cov <- cov_func[midP+fixLag]
	dyn_cov <- cov_func[dynLag[1]]

	Hz <- n/avg_t
	if(is.null(xlim)){
		xlim <- c(-avg_t/2,avg_t/2)
		x <- seq(-avg_t/2,avg_t/2,length.out=n)
	} else {
		xlim2 <- round(xlim*Hz)
		x_hi <- intersect(seq(xlim2[1],xlim2[2]),seq(n)-midP)
		x <- x_hi/Hz
		cov_func <- cov_func[midP + x_hi]
	}
	plot(0, type="n", xlim=xlim, ylim=range(cov_func,na.rm=TRUE), xlab=expression(paste(italic(tau)," (s)")), ylab=ylab, cex.axis=cx, cex.lab=cx)
	abline(h=0, lwd=1.5, lty=2, col="gray60"); abline(v=0, lwd=1.5, lty=2, col="gray60")
	# abline(h=fix_cov,lty=3,col="lightgrey")
	# abline(h=dyn_cov,lty=4,col="lightgrey")
	abline(v=fixLag/Hz,lty=3,col="lightgrey")
	abline(v=dynLag[2]/Hz,lty=4,col="lightgrey")
	lines(x, cov_func, lwd=2, col=cl)
	if(fixLag==dynLag[2]){
		mtext(substitute(paste(italic(tau)==a,"s: ",y),list(y=sprintf("%1.8f",dyn_cov),a=sprintf("%1.1f",fixLag/Hz))),side=3,cex=cxmt)
	} else {
		mtext(substitute(paste(italic(tau)[fix]==a,"s: ",x," / ",italic(tau)[dyn]==b,"s: ",y),list(x=sprintf("%1.8f",fix_cov),y=sprintf("%1.8f",dyn_cov),a=sprintf("%1.1f",fixLag/Hz),b=sprintf("%1.1f",dynLag[2]/Hz))),side=3,cex=cxmt)
	}
}

plot_cospec_ogive <- function(ogive,cospec,freq,ylab=NULL,xlim=NULL,cx=1.5,col="lightblue",nred=floor(sqrt(sqrt(length(ogive)))*3)){
	# reduced cospec 1:
	cospec_reduced0 <- reduce_cospec(cospec,freq,nred*10)
	cospec_f <- cospec_reduced0$cospec
	# reduced cospec 2:
	cospec_reduced <- reduce_cospec(cospec,freq,nred)
	cospec_rm <- cospec_reduced$cospec

	rCo <- range(cospec_rm,na.rm=TRUE)
	if(is.null(xlim))xlim <- rev(range(freq))
	ylim <- c(min(ogive,0),max(0,max(ogive)))
	pxlim <- pretty(log10(xlim),n=ceiling(abs(diff(log10(xlim)))))
	pxlims <- rep(pxlim,each=9) - 1 + log10(seq(9, 1))
	pylim <- pretty(ylim)
	prCo <- pretty(rCo) 
	y_cf <- (cospec_f - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	y_crm <- (cospec_rm - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	py2 <- (prCo - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	
	plot(1,xlim=xlim,ylim=ylim, cex.axis=cx, cex.lab=cx,type="n",log="x",xaxt="n",yaxt="n",xlab="frequency [Hz]",ylab="",panel.first=abline(h=0,col=col,lty=2))
	abline(h=(0 - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1],lty=2,col="darkgrey")
	axis(1,at=10^pxlims,labels=FALSE,tck=-0.01, cex.axis=cx, cex.lab=cx)
	axis(1,at=10^pxlim,labels=10 ^ pxlim, cex.axis=cx, cex.lab=cx)
	axis(2,at=pylim,labels=pylim, cex.axis=cx, cex.lab=cx, col=col,col.axis=col,lwd=2,font=2)
	title(ylab=ylab,col.lab=col, cex.lab=cx,font.lab=2)
	axis(4,at=py2,labels=prCo, cex.axis=cx, cex.lab=cx)
    at3 <- 1 / c(0.1, 1, 10, 60, 600)
    lab3 <- c('100ms', '1s', '10s', '1min', '10min')
    sub3 <- 1 / c(seq(0.01, 0.09, by = 0.01), seq(0.2, 0.9, by = 0.1),
        seq(2, 9), seq(20, 50, by = 10), seq(120, 540, by = 60)) 
	axis(3,at=at3,labels=lab3, cex.axis = cx, cex.lab = cx)
	axis(3,at=sub3,labels=NA, cex.axis = cx, cex.lab = cx, tck = -0.01)
	lines(cospec_reduced0$freq,y_cf,col="lightgrey")
	lines(cospec_reduced$freq,y_crm,type="b",col="black",lwd=2)
	lines(freq,ogive,col=col,lwd=2)
}

plot_damping <- function(ogive_damp,freq,ylab=NULL,xlim=NULL,ylim=NULL,cx=1.5,cx.leg=1.5,col="lightblue",main=NULL){
	ogive <- ogive_damp$ogive
	ogive_ref_pbreg <- ogive_damp$ogive_ref_pbreg
	ogive_ref_deming <- ogive_damp$ogive_ref_deming
	dampf_pbreg <- ogive_damp$dampf_pbreg
	dampf_deming <- ogive_damp$dampf_deming
	freq.limits <- ogive_damp$freq.limits

	if(is.null(xlim))xlim <- rev(range(freq))
	if(is.null(ylim))ylim <- c(min(0,ogive,ogive_ref_pbreg,ogive_ref_deming),max(0,max(ogive,ogive_ref_pbreg,ogive_ref_deming)))
	pxlim <- pretty(log10(xlim),n=ceiling(abs(diff(log10(xlim)))))
	pxlims <- rep(pxlim,each=9) - 1 + log10(seq(9, 1))
	pylim <- pretty(ylim)

	plot(1,xlim=xlim,ylim=ylim,main=main, cex.axis=cx, cex.lab=cx,type="n",log="x",xaxt="n",xlab="frequency [Hz]",ylab=ylab,panel.first=abline(h=0,col=col,lty=2))
	abline(h=0,lty=2,col="darkgrey")
	axis(1,at=10^pxlims,labels=FALSE,tck=-0.01, cex.axis=cx, cex.lab=cx)
	axis(1,at=10^pxlim,labels=10 ^ pxlim, cex.axis=cx, cex.lab=cx)
    at3 <- 1 / c(0.1, 1, 10, 60, 600)
    lab3 <- c('100ms', '1s', '10s', '1min', '10min')
    sub3 <- 1 / c(seq(0.01, 0.09, by = 0.01), seq(0.2, 0.9, by = 0.1),
        seq(2, 9), seq(20, 50, by = 10), seq(120, 540, by = 60)) 
	axis(3,at=at3,labels=lab3, cex.axis = cx, cex.lab = cx)
	axis(3,at=sub3,labels=NA, cex.axis = cx, cex.lab = cx, tck = -0.01)
	lines(freq,ogive_ref_deming,col="darkgrey",lwd=2)
	lines(freq,ogive_ref_pbreg,col="lightgrey",lwd=2)
	lines(freq,ogive,col=col,lwd=2)
	abline(h=(1 - 1/dampf_deming)*ogive[1],col="darkgrey",lty=2)
	abline(h=(1 - 1/dampf_pbreg)*ogive[1],col="lightgrey",lty=2)
	abline(v=1/freq.limits,lty=4,col="black")
	pos <- if(diff(abs(ylim))>0) "topleft" else "bottomleft"
	legend(pos,legend=c("damped","scaled reference (pbreg)","scaled reference (deming)",sprintf("damping (pbreg): %1.1f%%",(1 - dampf_pbreg)*100),sprintf("damping (deming): %1.1f%%",(1 - dampf_deming)*100)),col=c(col,"lightgrey","darkgrey",NA,NA),bty="n",lty=1,lwd=2,cex=cx.leg)
}

# wpl function
correct_wpl <- function(ec_dat, fluxes = c('nh3_ugm3', 'co2_mmolm3'), 
    dynamic_lag = TRUE, return_extra_cols = FALSE, water = 'h2o_mmolm3', 
    temperature = 'li_temp_amb', pressure = 'li_press_amb',
    coefs = list(ht8700 = c(A = 1, B = 1, C = 1), licor = c(A = 1, B = 1, C = 1))) {
    # switch to data.table
    if (is_ibts <- is.ibts(ec_dat)) {
        ec_dat <- as.data.table(ec_dat)
    }
    out <- copy(ec_dat)
    dynfix <- ifelse(dynamic_lag, 'dyn_', 'fix_')
    cat('using', ifelse(dynamic_lag, 'dynamic', 'fixed'), 'lag fluxes\n')
    # get pressure in Pa
    f_press <- switch(pressure
        , li_press_amb = 1e3
        , ht_press_amb = 1e2
        , stop('unknown pressure variable: ', pressure)
    )
    out[, p_pa := get(pressure) * f_press]
    # get temperature in K
    o_temp <- switch(temperature
        , li_temp_amb = 
        , ht_temp_amb = 273.15
        , T_sonic = 0
        , stop('unknown temperature variable: ', temperature)
    )
    out[, t_k := get(temperature) + o_temp]
    # dry air molar mass (g/mol)
    M_dry <- 28.965
    # molar gas constant
    R <- 8.31446
    # get dry air density in g/m3 (p*V=m/M*R*T->m/V=p*M/R/T)
    out[, rho_dry := p_pa * M_dry / R / t_k]
    # h2o molar mass (g/mol)
    M_h2o <- 18.015
    # get h2o density/flux in g/m3
    f_h2o <- switch(water
        , h2o_mmolm3 = M_h2o * 1e-3
        , stop('only h2o_mmolm3 accepted by now')
    )
    # density
    out[, rho_h2o := get(paste0('avg_', water)) * f_h2o]
    # flux
    out[, flux_h2o := get(paste0('flux_', dynfix, 'wx', water)) * f_h2o]
    # get nh3 density/flux in g/m3
    if (nh3_ppb <- 'nh3_ppb' %in% fluxes) {
        stop('nh3_ppb is not implemented yet...')
    }
    if (nh3_ugm3 <- 'nh3_ugm3' %in% fluxes) {
        f_nh3 <- 1e-6
        # density
        out[, rho_nh3 := avg_nh3_ugm3 * f_nh3]
        # raw flux
        out[, flux_nh3_raw := get(paste0('flux_', dynfix, 'wxnh3_ugm3')) * f_nh3]
        # corrected flux
        out[, flux_nh3_gm2s := coefs$ht8700['A'] * (
            # first term
            flux_nh3_raw +
            # second term
            coefs$ht8700['B'] * rho_nh3 / rho_dry * flux_h2o +
            # third term
            coefs$ht8700['C'] * (1 + M_dry / M_h2o * rho_h2o / rho_dry) * 
                rho_nh3 / t_k * cov_wT
            )]
        # get wpl factor
        out[, flux_nh3_wpl_factor := flux_nh3_gm2s / flux_nh3_raw]
        # convert back to ug/m3
        out[, flux_wpl_nh3 := flux_nh3_gm2s / f_nh3]
    }
    # get co2 density/flux in g/m3
    if (co2_mmolm3 <- 'co2_mmolm3' %in% fluxes) {
        M_co2 <- 44.01
        f_co2 <- M_co2 * 1e-3
        # density
        out[, rho_co2 := avg_co2_mmolm3 * f_co2]
        # raw flux
        out[, flux_co2_raw := get(paste0('flux_', dynfix, 'wxco2_mmolm3')) * f_co2]
        # corrected flux
        out[, flux_co2_gm2s := coefs$licor['A'] * (
            # first term
            flux_co2_raw +
            # second term
            coefs$licor['B'] * rho_co2 / rho_dry * flux_h2o +
            # third term
            coefs$licor['C'] * (1 + M_dry / M_h2o * rho_h2o / rho_dry) * 
                rho_co2 / t_k * cov_wT
            )]
        # get wpl factor
        out[, flux_co2_wpl_factor := flux_co2_gm2s / flux_co2_raw]
        # convert back to mmol/m3
        out[, flux_wpl_co2 := flux_co2_gm2s / f_co2]
    }
    if (!return_extra_cols) {
        out <- out[, .SD, .SDcols = c(
            names(ec_dat)
            , if (nh3_ugm3) 'flux_wpl_nh3'
            , if (co2_mmolm3) 'flux_wpl_co2'
            )]
    }
    # switch back to ibts
    if (is_ibts) {
        out <- as.ibts(out)
    }
    out
}

# helper function for process_ec_fluxes()
fix_defaults <- function(x, vars) {
    nms <- names(x)
    if (is.null(nms)) {
        if (length(x) != 1) {
            stop(deparse(substitute(x)), 
                ' - provide either a named vector or a single value!')
        }
        nms <- ''
    }
    # find vars not in x
    missing_vars <- vars[!(vars %in% nms)]
    if (length(missing_vars)) {
        # get formals
        frms <- formals(sys.function(sys.parent(1L)))
        # get formal name
        fname <- as.character(substitute(x))
        # get default values
        default <- eval(frms[[fname]])
        # extend/replace x with default
        if (nms[1] == '') {
            # replace
            x <- default[missing_vars]
        } else {
            # extend
            x <- c(x, default[missing_vars])
        }
    }
    # return
    x[vars]
}

## main flux processing function
process_ec_fluxes <- function(
		sonic_directory
        , ht_directory = NULL
        , licor_directory = NULL
		, start_time = NULL
		, end_time = NULL
		, avg_period = '30mins'
		, tz_user = 'UTC'
		, dev_north = NULL
        , declination = NULL
		, z_ec = NULL
		, z_canopy = NULL
        # detrending -> valid entries are blockAVG,linear,linear_robust,ma_xx (xx = time in seconds)
        , detrending = c(u = 'linear', v = 'linear', w = 'linear', T = 'linear', nh3_ppb = 'linear', nh3_ugm3 = 'linear', h2o_mmolm3 = 'linear', co2_mmolm3 = 'linear')
        , hard_limits = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, nh3_ppb = TRUE, nh3_ugm3 = TRUE, h2o_mmolm3 = TRUE, co2_mmolm3 = TRUE)
        , hard_limits_lower = c(u = -30, v = -30, w = -10, T = 243, nh3_ppb = -100, nh3_ugm3 = -100, h2o_mmolm3 = -100, co2_mmolm3 = -100)
        , hard_limits_upper = c(u = 30, v = 30, w = 10, T = 333, nh3_ppb = 5000, nh3_ugm3 = 5000, h2o_mmolm3 = 5000, co2_mmolm3 = 5000)
        , hard_limits_window = '5mins'
        , hard_limits_replace = FALSE
		, covariances = c('uxw', 'wxT', 'wxnh3_ugm3', 'wxh2o_mmolm3', 'wxco2_mmolm3')
        # fix lag in seconds
		, lag_fix = c(uxw = 0, wxT = 0, wxnh3_ppb = -0.4, wxnh3_ugm3 = -0.4, wxh2o_mmolm3 = -0.2, wxco2_mmolm3 = -0.2)
        # dyn lag in seconds around lag_fix
		, lag_dyn = c(uxw = 0.2, wxT = 0.2, wxnh3_ppb = 1.5, wxnh3_ugm3 = 1.5, wxh2o_mmolm3 = 1.5, wxco2_mmolm3 = 1.5)
        # re_rmse window
        , gamma_time_window = c(2.5, 5)
		, damping_reference = c(wxnh3_ppb = 'wxT', wxnh3_ugm3 = 'wxT', wxh2o_mmolm3 = 'wxT', wxco2_mmolm3 = 'wxT')
        # lower & upper bounds of fitting ogives (in seconds)
		, damping_lower = c(wxnh3_ppb = 2, wxnh3_ugm3 = 2, wxh2o_mmolm3 = 2, wxco2_mmolm3 = 2)
		, damping_upper = c(wxnh3_ppb = 20, wxnh3_ugm3 = 20, wxh2o_mmolm3 = 20, wxco2_mmolm3 = 20)
        , subintervals = TRUE
        , subint_n = 5
        , subint_detrending = c(u = 'linear', v = 'linear', w = 'linear', T = 'linear', nh3_ppb = 'linear', nh3_ugm3 = 'linear', h2o_mmolm3 = 'linear', co2_mmolm3 = 'linear')
        , oss_threshold = 0
        , na_alarm_code = c(1:3, 5:8, 11, 13)
        , thresh_period = 0.75
		, create_graphs = TRUE
		, graphs_directory = NULL
		, add_name = ''
        , plot_timeseries = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, ht_oss = TRUE, nh3_ppb = FALSE, nh3_ugm3 = TRUE, li_co2ss = TRUE, h2o_mmolm3 = TRUE, co2_mmolm3 = TRUE)
        , plotting_var_units = c(u = 'm/s', v = 'm/s', w = 'm/s', T = 'K', ht_oss = '-', nh3_ppb = 'ppb', nh3_ugm3 = 'ug/m3', li_co2ss = '-', h2o_mmolm3 = 'mmol/m3', co2_mmolm3 = 'mmol/m3')
        , plotting_var_colors = c(u = 'gray20', v = 'gray20', w = 'gray20', T = 'orange', ht_oss = 'grey', nh3_ppb = 'indianred', nh3_ugm3 = 'indianred', li_co2ss = 'grey', h2o_mmolm3 = '#8FC1E6', co2_mmolm3 = 'seagreen4')
        , plotting_covar_units = c(uxw = 'm2/s2', wxT = 'K*m/s', wxnh3_ppb = 'ppb*m/s', wxnh3_ugm3 = 'ug/m2/s', wxh2o_mmolm3 = 'mmol/m2/s', wxco2_mmolm3 = 'mmol/m2/s')
        , plotting_covar_colors = c(uxw = 'gray70', wxT = 'orange', wxnh3_ppb = 'indianred', wxnh3_ugm3 = 'indianred', wxh2o_mmolm3 = '#8FC1E6', wxco2_mmolm3 = 'seagreen4')
		, ogives_out = FALSE
        , as_ibts = TRUE
	){

	script.start <- Sys.time()
	################################################################################
	# ----------------------------- read config file -------------------------------

    # check input
    if (create_graphs) {
        # TODO: eventually remove one of either arguments
        if (is.null(graphs_directory)) {
            stop('argument "create_graphs" is TRUE, but "graphs_directory" is not specified!')
        }
        library(reshape2)
        library(lattice)
    }
    if (!is.null(graphs_directory) && (
            !is.character(graphs_directory) || !dir.exists(graphs_directory)
            )) {
        stop('argument "graphs_directory": directory "', graphs_directory, 
            '" does not exist!')
    }

    # check sonic function
    if (!exists('read_sonic', mode = "function")) {
        stop("function 'read_sonic' doesn't exist -> please source gel script",
        " read-sonic-data.r")
    }

    # check data for nh3, co2 & h2o
    if (is.data.table(sonic_directory)) {
        if (any(grepl('nh3', names(sonic_directory)))) {
            ht_directory <- 'data provided with sonic'
        }
        if (any(grepl('h2o', names(sonic_directory)))) {
            licor_directory <- 'data provided with sonic'
        }
    }

    # get scalars and fix missing instruments
    scalars <- sub('wx', '', grep('wx[^T].+', covariances, value = TRUE))
    # fix missing ht8700
    if (ht_null <- is.null(ht_directory)) {
        scalars <- grep('nh3', scalars, value = TRUE, invert = TRUE)
        covariances <- grep('nh3', covariances, value = TRUE, invert = TRUE)
    } else if (!exists('read_ht8700', mode = "function")) {
        stop("function 'read_ht8700' doesn't exist -> please source gel script",
        " ht8700-functions.r")
    }
    # fix missing licor
    if (licor_null <- is.null(licor_directory)) {
        scalars <- grep('h2o|co2', scalars, value = TRUE, invert = TRUE)
        covariances <- grep('h2o|co2', covariances, value = TRUE, invert = TRUE)
    } else if (!exists('read_licor', mode = "function")) {
        stop("function 'read_licor' doesn't exist -> please source gel script",
        " licor-functions.r")
    }

    # prepare variables etc.
    variables <- c('u', 'v', 'w', 'T', scalars)
    scalar_covariances <- grepl('wx[^T].+', covariances)
    names(scalar_covariances) <- covariances
    covariances_plotnames <- make.names(sub("x", "", covariances))
    names(covariances_plotnames) <- covariances
    covariances_variables <- strsplit(covariances, "x")
    hl_method <- if (hard_limits_replace) 'repl' else 'norepl'
    if (is.null(z_ec) || !is.numeric(z_ec)) {
        stop('argument "z_ec" must be provided as numeric value (height in m a.g.l)!')
    }
    if (is.null(z_canopy) || !is.numeric(z_canopy)) {
        stop('argument "z_canopy" must be provided as numeric value ',
            '(height of canopy in meters)!')
    }
    if (is.null(dev_north) || !is.numeric(dev_north)) {
        stop('argument "dev_north" must be provided as numeric value!')
    }
    if (is.null(declination)) {
        stop('argument "declination" must be provided!',
        ' -> it is also possible to provide a list with lon/lat entries...')
    } else if (is.list(declination) || length(declination) == 2) {
        if (!require(oce)) {
            stop('R package "oce" must be installed when lat/lon is provided')
        }
        mag_dec <- \(x) oce::magneticField(declination$lon, declination$lat, x)$declination
    } else if (!is.numeric(declination) || length(declination) != 1) {
        stop('argument "declination" must be a single numeric value!')
    } else {
        mag_dec <- \(x) declination
    }

    # fix input (vectors of default values)
    detrending <- fix_defaults(detrending, variables)
    hard_limits <- fix_defaults(hard_limits, variables)
    hard_limits_lower <- fix_defaults(hard_limits_lower, variables)
    hard_limits_upper <- fix_defaults(hard_limits_upper, variables)
    lag_fix <- fix_defaults(lag_fix, covariances)
    lag_dyn <- fix_defaults(lag_dyn, covariances)
    damping_reference <- fix_defaults(damping_reference, covariances[scalar_covariances])
    damping_lower <- fix_defaults(damping_lower, covariances[scalar_covariances])
    damping_upper <- fix_defaults(damping_upper, covariances[scalar_covariances])
    subint_detrending <- fix_defaults(subint_detrending, variables)
    ts_vars <- variables
    if (!ht_null) {
        # add ht8700 quality parameters (oss) to plotting
        ts_vars <- c(ts_vars, 'ht_oss')
    }
    if (!licor_null) {
        # add licor quality parameters (co2ss) to plotting
        ts_vars <- c(ts_vars, 'li_co2ss')
    }
    # extend & sort ts_vars
    ts_vars <- union(names(plot_timeseries), ts_vars)
    plot_timeseries <- fix_defaults(plot_timeseries, ts_vars)
    # only take TRUE
    plot_timeseries <- plot_timeseries[plot_timeseries]
    plotting_var_colors <- fix_defaults(plotting_var_colors, ts_vars)[names(plot_timeseries)]
    plotting_var_units <- fix_defaults(plotting_var_units, ts_vars)[names(plot_timeseries)]
    plotting_covar_units <- fix_defaults(plotting_covar_units, covariances)
    plotting_covar_colors <- fix_defaults(plotting_covar_colors, covariances)

    lim_range <- rbind(lower = hard_limits_lower, top = hard_limits_upper)
    damp_region <- mapply(c, damping_lower, damping_upper, SIMPLIFY = FALSE)

    # get flux variables and lag times:                                      
    # ------------------------------------------------------------------------------ 
    flux_variables <- unique(unlist(strsplit(covariances, split = "x")))

    # ------------------------------------------------------------------------------

    # get data
    # ------------------------------------------------------------------------------

    cat("\n************************************************************\n")
    cat("HT8700/LI-COR EC flux processing\n")

    # mandatory sonic data
    if (sonic_has_data <- inherits(sonic_directory, 'data.table')) {
        cat('Sonic Anemometer: raw data provided...\n')
        setDT(sonic_directory)
        # set time zone to UTC
        sonic_directory[, Time := with_tz(Time, 'UTC')]
        # get start & end
        start_sonic <- sonic_directory[, Time[1]]
        end_sonic <- sonic_directory[, Time[.N]]
    } else {
        # get files
        sonic_files <- dir(sonic_directory, pattern = '^(py_)?fnf_.*_sonic_.*')
        if (sonic_old_format <- length(sonic_files) == 0) {
            cat('Implement old CET format sonic!\n')
            browser()
            # old loggerbox format
            sonic_files <- dir(sonic_directory, pattern = '^data_sonic-.')
        }
        # check files available
        if (length(sonic_files) == 0) {
            stop('No sonic data available in "', sonic_directory, '"!')
        }
        # prioritize py_ files
        # get py_
        i_py <- grepl('^py_', sonic_files)
        # remove duplicated non-py
        sonic_files <- c(
            setdiff(sonic_files[!i_py], sub('^py_', '', sonic_files[i_py])),
            sonic_files[i_py]
        )
        # prioritize gz files
        i_gz <- grepl('\\.gz$', sonic_files)
        # remove duplicated non-py
        sonic_files <- c(
            setdiff(sonic_files[!i_gz], sub('\\.gz$', '.csv', sonic_files[i_gz])),
            sonic_files[i_gz]
        )
        # get date
        sonic_dates <- sub('^(py_)?fnf_01_sonic_', '', sonic_files)
        # sort by date
        sonic_files <- sonic_files[order(sonic_dates)]
        # sort files by date (& time)
        sonic_files <- sort(sonic_files)
        # get start
        if (sonic_old_format) {
            sonic_pattern <- c('.*_(\\d{8})_(\\d{6}).*', '\\1\\2', 'Etc/GMT-1', 
                '\\1000000')
        } else {
            sonic_pattern <- c('.*_(\\d{4})_(\\d{2})_(\\d{2}).csv', '\\1\\2\\3000000',
                'UTC')
            sonic_pattern[4] <- sonic_pattern[2]
        }
        start_sonic <- strptime(sub(sonic_pattern[1], sonic_pattern[2], sonic_files[1]),
            '%Y%m%d%H%M%S', tz = sonic_pattern[3])
        # get end (last date + 24h)
        end_sonic <- strptime(sub(sonic_pattern[1], sonic_pattern[4], tail(sonic_files, 1)),
            '%Y%m%d%H%M%S', tz = sonic_pattern[3]) + 24 * 3600
    }

    # optional ht8700 data
    ht_with_sonic <- FALSE
    if (ht_provided <- !is.null(ht_directory)) {
        if (ht_has_data <- inherits(ht_directory, 'data.table')) {
            cat('HT8700: raw data provided in input...\n')
            setDT(ht_directory)
            # set time zone to UTC
            ht_directory[, Time := with_tz(Time, 'UTC')]
            # get start & end
            start_ht <- ht_directory[, Time[1]]
            end_ht <- ht_directory[, Time[.N]]
        } else if (!(ht_with_sonic <- ht_directory == 'data provided with sonic')) {
            cat('HT8700: data provided...\n')
            # get files
            ht_files <- dir(ht_directory, pattern = '^(py_)?fnf_.*_ht8700_.*')
            if (ht_old_format <- length(ht_files) == 0) {
                cat('Implement old CET format ht8700!\n')
                # old loggerbox format
                ht_files <- dir(ht_directory, pattern = '^ht8700_sonic-.')
            }
            # check files available
            if (length(ht_files) == 0) {
                stop('No ht8700 data available in "', ht_directory, '"!')
            }
            # prioritize py_ files
            # get py_
            i_py <- grepl('^py_', ht_files)
            # remove duplicated non-py
            ht_files <- c(
                setdiff(ht_files[!i_py], sub('^py_', '', ht_files[i_py])),
                ht_files[i_py]
            )
            # prioritize gz files
            i_gz <- grepl('\\.gz$', ht_files)
            # remove duplicated non-py
            ht_files <- c(
                setdiff(ht_files[!i_gz], sub('\\.gz$', '.csv', ht_files[i_gz])),
                ht_files[i_gz]
            )
            # get date
            ht_dates <- sub('^(py_)?fnf_01_ht8700_', '', ht_files)
            # sort by date
            ht_files <- ht_files[order(ht_dates)]
            # get start
            if (ht_old_format) {
                ht_pattern <- c('.*_(\\d{8})_(\\d{6}).*', '\\1\\2', 'Etc/GMT-1', '\\1000000')
            } else {
                ht_pattern <- c('.*_(\\d{4})_(\\d{2})_(\\d{2}).csv', '\\1\\2\\3000000', 'UTC')
                ht_pattern[4] <- ht_pattern[2]
            }
            start_ht <- strptime(sub(ht_pattern[1], ht_pattern[2], ht_files[1]),
                '%Y%m%d%H%M%S', tz = ht_pattern[3])
            # get end (last date + 24h)
            end_ht <- strptime(sub(ht_pattern[1], ht_pattern[4], tail(ht_files, 1)),
                '%Y%m%d%H%M%S', tz = ht_pattern[3]) + 24 * 3600
        } else if (ht_with_sonic) {
            cat('HT8700: raw data provided with sonic input...\n')
        }
    }

    # optional licor data
    licor_with_sonic <- FALSE
    if (licor_provided <- !is.null(licor_directory)) {
        if (licor_has_data <- inherits(licor_directory, 'data.table')) {
            cat('LI-7500: raw data provided in input...\n')
            licor_files <- NULL
            setDT(licor_directory)
            # set time zone to UTC
            licor_directory[, Time := with_tz(Time, 'UTC')]
            # get start & end
            start_licor <- licor_directory[, Time[1]]
            end_licor <- licor_directory[, Time[.N]]
        } else if (!(licor_with_sonic <- licor_directory == 'data provided with sonic')) {
            cat('LI-7500: data provided...\n')
            licor_files <- dir(licor_directory, pattern = '^(py_)?fnf_.*_licor_.*')
            if (length(licor_files) == 0) {
                stop('No licor data available in "', licor_directory, '"!')
            }
            # prioritize py_ files
            # get py_
            i_py <- grepl('^py_', licor_files)
            # remove duplicated non-py
            licor_files <- c(
                setdiff(licor_files[!i_py], sub('^py_', '', licor_files[i_py])),
                licor_files[i_py]
            )
            # prioritize gz files
            i_gz <- grepl('\\.gz$', licor_files)
            # remove duplicated non-py
            licor_files <- c(
                setdiff(licor_files[!i_gz], sub('\\.gz$', '.csv', licor_files[i_gz])),
                licor_files[i_gz]
            )
            # get date
            licor_dates <- sub('^(py_)?fnf_01_licor_', '', licor_files)
            # sort by date
            licor_files <- licor_files[order(licor_dates)]
            # get start & end
            licor_pattern <- c('.*_(\\d{4})_(\\d{2})_(\\d{2}).csv', '\\1\\2\\3000000',
                'UTC')
            start_licor <- strptime(sub(licor_pattern[1], licor_pattern[2], licor_files[1]),
                '%Y%m%d%H%M%S', tz = licor_pattern[3])
            # get end (last date + 24h)
            end_licor <- strptime(sub(licor_pattern[1], licor_pattern[2], tail(licor_files, 1)),
                '%Y%m%d%H%M%S', tz = licor_pattern[3]) + 24 * 3600
        } else if (licor_with_sonic) {
            cat('LI-7500: raw data provided with sonic input...\n')
        }
    } else {
        # check if ht available
        if (!ht_provided) {
            stop('neither ht8700 nor licor data or directory has been provided -> cannot process fluxes without concentration data!')
        }
    }

    cat("************************************************************\n")

    # parse time diff
    avg_secs <- parse_time_diff(avg_period)

    # check times
    if (is.null(start_time) && is.null(end_time)) {
        start_time <- 'first'
        end_time <- 'last'
    } 

    # fix start_time (only use sonic)
    if (length(start_time) == 1 && start_time == 'first') {
        start_time <- with_tz(start_sonic, tz_user)
    } else if (!is.null(start_time)) {
        start_time <- parse_date_time3(start_time, tz = tz_user)
    }

    # fix end_time (only use sonic)
    if (is.null(end_time)) {
        end_time <- start_time + avg_secs
    } else if (length(end_time) == 1 && end_time == 'last') {
        end_time <- with_tz(end_sonic, tz_user)
    } else {
        end_time <- parse_date_time3(end_time, tz = tz_user)
    }

    # fix end_time != NULL & start_time == NULL
    if (is.null(start_time)) {
        start_time <- end_time - avg_secs
    }

    # check & extend start/end times
    if (length(start_time) != length(end_time)) {
        stop('arguments "start_time" and "end_time" must have equal lengths!')
    } else if (length(start_time) == 1) {
        start_time <- seq(start_time, end_time, by = avg_secs)
        end_time <- start_time + avg_secs
    }

    # get vector of dates
    start_dates <- format(date(start_time), '%Y-%m-%d')
    # subtract 1 second to convert date of 00:00:00 to 23:59:59
    end_dates <- format(date(end_time - 1), '%Y-%m-%d')
    dates <- unique(c(start_dates, end_dates))

    # get start/end per date in UTC & CET
    times_list <- sapply(dates, \(x) {
        eind <- end_dates == x
        times <- unique(c(start_time[eind], end_time[eind] - 1))
        UTC_times <- with_tz(times, 'UTC')
        CET_times <- with_tz(times, 'Etc/GMT-1')
        list(
            UTC = unique(date(UTC_times)),
            CET = unique(date(CET_times)),
            UTC_times = list(
                st = with_tz(start_time[eind], 'UTC'),
                et = with_tz(end_time[eind], 'UTC')
            ),
            CET_times = list(
                st = with_tz(start_time[eind], 'Etc/GMT-1'),
                et = with_tz(end_time[eind], 'Etc/GMT-1')
            )
        )
    }, simplify = FALSE)

    # create output list
    result_list <- vector(mode = 'list', length(dates))
    names(result_list) <- dates

    # prepare ogive output
    if (ogives_out) {
        e_ogive <- new.env()
        e_ogive$Cospec_dyn_Out <- e_ogive$Cospec_fix_Out <- e_ogive$Covars_Out <- e_ogive$Ogive_fix_Out <- e_ogive$Ogive_dyn_Out <- setNames(vector("list", length(dates)), dates)
    }

    # copy scalars etc. to fix missing
    input_scalars <- scalars
    input_covariances <- covariances
    input_covariances_variables <- covariances_variables
    input_covariances_plotnames <- covariances_plotnames
    input_scalar_covariances <- scalar_covariances
    input_flux_variables <- flux_variables
    input_plot_timeseries <- plot_timeseries
    input_damping_reference <- damping_reference
    input_damp_region <- damp_region

    # scalar covariances
    scalar_covariances_only <- covariances[scalar_covariances]

    # create result folder (folder name includes first input-filename) and set this directory as working directory                                      
    # ------------------------------------------------------------------------------ 
    if (is.null(graphs_directory) || isFALSE(graphs_directory)) {
        create_graphs <- FALSE
    } else {					
        # tstamp <- format(Sys.time(), "%Y%m%d_%H%M")
        if (add_name != "") {
            # folder <- paste0("ec-fluxes-", add_name, "-eval", tstamp, "-avg", avg_period)
            folder <- paste0("ec-fluxes-", add_name, "-avg", avg_period)
        } else {
            # folder <- paste0("ec_fluxes-eval", tstamp, "-avg", avg_period)
            folder <- paste0("ec_fluxes-avg", avg_period)
        }
        path_folder <- file.path(graphs_directory, folder)
    }


    # loop over dates
    for (day in seq_along(dates)) {

        # day <- dates[1]
        cat(
            '\n\n-------------------------------------------',
            '\nReading raw data for', 
            dates[[day]],
            '\n-------------------------------------------\n\n'
        )

        # prepare dates & times
        dates_utc <- times_list[[day]]$UTC
        dates_cet <- times_list[[day]]$CET
        times_utc <- times_list[[day]]$UTC_times
        times_cet <- times_list[[day]]$CET_times

        # reading data file:                                      
        # ------------------------------------------------------------------------------

        # get sonic data
        if (sonic_has_data) {
            cat('~~~\nSubsetting sonic data - ')
            # copy
            sonic <- sonic_directory[Time >= times_utc$st[1] &
                Time < tail(times_utc$et, 1), ]
        } else {
            cat('~~~\nReading sonic files - ')
            # get data from old files
            if (day > 1 && !is.null(sonic)) {
                # copy old data
                old_sonic <- sonic[Time >= times_utc$st[1], ]
                # remove old date
                dates_sonic <- dates_utc[!(dates_utc %in% dates_sonic)]
            } else {
                dates_sonic <- dates_utc
                old_sonic <- NULL
            }
            # get files
            sonic_read <- grep(
                gsub('-', '_', dates_sonic, fixed = TRUE)
                , sonic_files, value = TRUE)
            if (length(sonic_read)) {
                # read new sonic files
                capture.output(
                    sonic <- rbindlist(lapply(file.path(sonic_directory, sonic_read), 
                            read_sonic))
                )
            } else {
                sonic <- NULL
            }
            # no UTC conversion needed, since this case is excluded up to now
            # # convert to UTC
            # sonic[, Time := with_tz(Time, 'UTC')]
            # bind together
            sonic <- rbind(sonic, old_sonic)
            if (is.null(sonic)) {
                cat('No sonic data for current day. -> skipping day.\n')
                next
            }
        }
        cat('done\n')

        # get ht8700 data
        if (ht_provided && ht_has_data) {
            cat('Subsetting HT8700 data - ')
            # copy
            ht <- ht_directory[Time >= times_utc$st[1] &
                Time < tail(times_utc$et, 1), ]
            cat('done\n')
        } else if (ht_provided && !ht_with_sonic) {
            cat('Reading HT8700 files - ')
            # get data from old files
            if (day > 1 && !is.null(ht)) {
                # copy old data
                old_ht <- ht[Time >= times_utc$st[1], ]
                # remove old date
                dates_ht <- dates_utc[!(dates_utc %in% dates_ht)]
            } else {
                dates_ht <- dates_utc
                old_ht <- NULL
            }
            # get files
            ht_read <- grep(
                gsub('-', '_', dates_ht, fixed = TRUE)
                , ht_files, value = TRUE)
            if (length(ht_read)) {
                # read new ht files
                capture.output(
                    ht <- rbindlist(lapply(file.path(ht_directory, ht_read), read_ht8700))
                )
            } else {
                ht <- NULL
            }
            # no UTC conversion needed, since this case is excluded up to now
            # # convert to UTC
            # ht[, Time := with_tz(Time, 'UTC')]
            # bind together
            ht <- rbind(ht, old_ht[, names(ht), with = FALSE])
            if (is.null(ht)) {
                cat('No HT8700 data for current day.\n')
            } else {
                cat('done\n')
            }
        } else {
            # no ht data
            ht <- NULL
        }

        # check HT8700 quality: alarm codes & OSS
        if (ht_provided && ht_with_sonic) {
            cat('Checking oss and alarm codes - ')
            # add alarm code
            if ('ht_alarm_code' %in% names(sonic)) {
                # rename all ht_* variables
                ht_vars <- grep('ht_', names(sonic), value = TRUE)
                setnames(sonic, ht_vars, sub('ht_', '', ht_vars, fixed = TRUE))
            } else if (!('alarm_code' %in% names(sonic))) {
                sonic[, alarm_code := get_alarms(.SD)]
            }
            # create regex pattern
            na_alarm_pattern <- paste(paste0('\\b', na_alarm_code, '\\b'), collapse = '|')
            # check alarms and set nh3 NA
            nh3_vars <- grep('nh3', names(sonic), value = TRUE)
            na0 <- sonic[, sum(is.na(get(nh3_vars[1])))]
            sonic[grepl(na_alarm_pattern, alarm_code),
                (nh3_vars) := NA_real_]
            na1 <- sonic[, sum(is.na(get(nh3_vars[1])))]
            # check oss
            sonic[oss < oss_threshold,
                (nh3_vars) := NA_real_]
            cat('done\nBad alarms:', na1 - na0, '\nValues below OSS:', 
                sonic[, sum(oss < oss_threshold)], '\n')
        } else if (ht_provided && !is.null(ht)) {
            cat('Checking oss and alarm codes - ')
            # add alarm code
            if (!('alarm_code' %in% names(ht))) {
                ht[, alarm_code := get_alarms(.SD)]
            }
            # create regex pattern
            na_alarm_pattern <- paste(paste0('\\b', na_alarm_code, '\\b'), collapse = '|')
            # check alarms and set nh3 NA
            nh3_vars <- grep('nh3', names(ht), value = TRUE)
            na0 <- ht[, sum(is.na(get(nh3_vars[1])))]
            ht[grepl(na_alarm_pattern, alarm_code),
                (nh3_vars) := NA_real_]
            na1 <- ht[, sum(is.na(get(nh3_vars[1])))]
            # check oss
            ht[oss < oss_threshold,
                (nh3_vars) := NA_real_]
            cat('done\nBad alarms:', na1 - na0, '\nValues below OSS:', 
                ht[, sum(oss < oss_threshold)], '\n')
        }

        # get licor data
        if (licor_provided && licor_has_data) {
            cat('Subsetting LI-7500 data - ')
            # copy
            licor <- licor_directory[Time >= times_utc$st[1] &
                Time < tail(times_utc$et, 1), ]
            cat('done\n')
        } else if (licor_provided && !licor_with_sonic) {
            cat('Reading LI-7500 files - ')
            # get data from old files
            if (day > 1 && !is.null(licor)) {
                # copy old data
                old_licor <- licor[Time >= times_utc$st[1], ]
                # remove old date
                dates_licor <- dates_utc[!(dates_utc %in% dates_licor)]
            } else {
                dates_licor <- dates_utc
                old_licor <- NULL
            }
            # get files
            licor_read <- grep(
                gsub('-', '_', dates_licor, fixed = TRUE)
                , licor_files, value = TRUE)
            if (length(licor_read)) {
                # read new licor files
                capture.output(
                    licor <- rbindlist(lapply(file.path(licor_directory, licor_read), 
                            read_licor))
                )
            } else {
                licor <- NULL
            }
            # no UTC conversion needed, since this case is excluded up to now
            # # convert to UTC
            # licor[, Time := with_tz(Time, 'UTC')]
            # bind together
            licor <- rbind(licor, old_licor[, names(licor), with = FALSE])
            if (is.null(licor)) {
                cat('No LI-7500 data for current day.\n')
            } else {
                cat("done\n")
            }
        } else {
            # no licor data
            licor <- NULL
        }

        cat('Merging files - ')

        daily_data <- merge_data(sonic, ht, licor)

        cat('done\n~~~\n')

        # check if empty
        if (nrow(daily_data) == 0) {
            cat('No valid data available...\n')
            next
        }
        
        # check variables and covars and subset by columns
        if (!all(cn_check <- variables %in% names(daily_data))) {
            stop(
                "cannot find variable(s): ", 
                paste(variables[!cn_check], collapse = ", "),
                "\nAvailable column names are: ", 
                paste(names(daily_data), collapse = ", ")
            )
        }
        
        # define coordinate system
        if (daily_data[, sonic[1] == 'HS']) {
            coord.system <-  'HS-Wauwilermoos'
        } else {
            coord.system <-  'Windmaster'
        }

        # get intervals:                                      
        # ------------------------------------------------------------------------------ 
        st_interval <- times_utc$st
        et_interval <- times_utc$et

        # define bins & subset again & just make sure sonic has no missing data
        daily_data <- daily_data[, bin := getIntervals(Time, st_interval, et_interval)][
            bin != 0L & is.finite(u) & is.finite(v) & is.finite(w) & is.finite(T) &
                !is.na(Time)]

        # check again if empty
        if (nrow(daily_data) == 0) {
            cat('No valid data available...\n')
            # advance to next day
            next
        }
        
        # be verbose
        daily_data[, {
            cat("\n~~~\nCalculation will include",
                uniqueN(bin)
                ,"intervals between", 
                format(Time[1], format = "%Y-%m-%d %H:%M"), "and", 
                format(Time[.N], 
                    format = "%Y-%m-%d %H:%M", usetz = TRUE), 
                "on a", avg_period, "min basis\n~~~\n\n")
        }]

        # --------------------- read fix lags from lag_lookuptable ---------------------
        # TODO: only if necessary/wanted

        # loop over individual intervals:
        # ------------------------------------------------------------------------------ 
        .Hz <- daily_data[, {
            d_t <- diff(as.numeric(Time))
            round(1 / median(d_t, na.rm = TRUE), -1)
        }]
        if (.Hz < 10) {
            cat('Frequency is lower than 10 Hz! Skipping processing fluxes!\n')
            next
        }
        n_period <- avg_secs * .Hz
        n_threshold <- thresh_period * n_period

        # fix and dyn lags:                                      
        # ------------------------------------------------------------------------------ 
        input_fix_lag <- round(lag_fix * .Hz)
        input_dyn_lag <- rbind(
            lower = round((lag_fix - lag_dyn) * .Hz)
            , upper = round((lag_fix + lag_dyn) * .Hz)
        )

        # loop over bins
        result_list[[day]] <- daily_data[, {

            cat("\n\n~~~~~~~~\n", dates[[day]], ": interval ", .GRP, " of ", .NGRP, "\n", sep = '')
            # get subset:
            # --------------------------------------------------------------------------
            if (.N > n_threshold) {

                # ugly copy because of inability to change .SD values
                SD <- copy(.SD)

                # times, frequencies etc.
                # --------------------------------------------------------------------------
                Int_Start <- st_interval[.BY[[1]]]
                freq <- .Hz * seq(floor(n_period / 2)) / floor(n_period / 2)
                # TODO (maybe for later): include NA where d_t>2*mean, remove entries where d_t < 0.5*mean

                # raw data quality control I, i.e. hard limits = physical range
                # --------------------------------------------------------------------------
                cat("~~~\nchecking hard limits...\n")
                if (any(hard_limits)) {
                    SD[, variables[hard_limits] := H.flags(mget(variables[hard_limits]), Time, .Hz, lim_range[, variables[hard_limits], drop = FALSE], hard_limits_window, hl_method)]
                } 

                # check NA values in scalars
                scalars <- input_scalars
                covariances <- input_covariances
                covariances_variables <- input_covariances_variables
                covariances_plotnames <- input_covariances_plotnames
                scalar_covariances <- input_scalar_covariances
                flux_variables <- input_flux_variables
                plot_timeseries <- input_plot_timeseries
                fix_lag <- input_fix_lag
                dyn_lag <- input_dyn_lag
                damping_reference <- input_damping_reference
                damp_region <- input_damp_region
                if (SD[, anyNA(.SD), .SDcols = scalars]) {
                    # loop over scalars & check
                    for (s in scalars) {
                        x <- SD[, unlist(mget(s, ifnotfound = NA))]
                        if (sum(is.finite(x)) < n_threshold) {
                            cat('=> ', s, ': number of valid measurements is below',
                                ' threshold -> exclude from current interval\n', sep = '')
                            scalars <- scalars[!(scalars %in% s)]
                            flux_variables <- flux_variables[!(flux_variables %in% s)]
                            plot_timeseries <- plot_timeseries[!(names(plot_timeseries) %in% s)]
                            sind <- grep(s, covariances)
                            if (length(sind)) {
                                damping_reference <- damping_reference[!(names(damping_reference) %in% covariances[sind])]
                                damp_region <- damp_region[!(names(damp_region) %in% covariances[sind])]
                                covariances <- covariances[-sind]
                                fix_lag <- fix_lag[-sind]
                                dyn_lag <- dyn_lag[, -sind, drop = FALSE]
                                covariances_variables <- covariances_variables[-sind]
                                covariances_plotnames <- covariances_plotnames[-sind]
                                scalar_covariances <- scalar_covariances[-sind]
                            }
                        }
                    }
                }


                # calculate wind direction, rotate u, v, w, possibly detrend T (+ u,v,w)
                # -------------------------------------------------------------------------- 
                cat("~~~\nrotating and deterending sonic data...\n")
                cat("     -> coordinate system:", coord.system, "\n")
                wind <- SD[, I(rotate_detrend(u, v, w, T, c.system = coord.system,
                        method = detrending[c("u", "v", "w", "T")], Hz_ts = .Hz))]
                SD[, c("u", "v", "w", "T") := wind[c("uprot", "vprot", "wprot", "Tdet")]]

                # calculate some turbulence parameters and collect some wind parameters
                # -------------------------------------------------------------------------- 
                wind_stats <- wind_statistics(wind, z_canopy, z_ec)
                ### correct for sonic north deviation and magn. declination
                current_declination <- mag_dec(Time[1])
                wind_stats["WD"] <- (wind_stats["WD"] + dev_north + current_declination) %% 360

                # switch to list with different lengths
                # only for scalar fluxes!
                # -> check na.action for omitted indices
                detrended_scalars <- NULL
                scalar_means <- setNames(
                    rep(NA_real_, length(scalars)),
                    scalars
                )
                scalar_sd <- setNames(
                    rep(NA_real_, length(scalars)),
                    scalars
                )
                if (length(scalars)) {
                    scalar_list <- SD[, I(lapply(.SD, na.omit)), .SDcols = scalars]

                # detrend scalars
                # -------------------------------------------------------------------------- 
                    cat("~~~\ndetrending scalars...\n")
                    detrended_scalars <- mapply(trend, y = scalar_list, method = 
                        detrending[scalars], MoreArgs = list(Hz_ts = .Hz), SIMPLIFY = FALSE)

                # calculate scalar averages and sd:
                # -------------------------------------------------------------------------- 
                    scalar_means[scalars] <- sapply(scalar_list, mean)
                    scalar_sd[scalars] <- sapply(scalar_list, sd)
                } else {
                    cat('No scalars available -> skipping interval!\n')
                    next
                }

                # start of flux relevant data manipulation
                # -------------------------------------------------------------------------- 
                # -------------------------------------------------------------------------- 
                cat("~~~\nstarting flux evaluation...\n")

                # get fft (keep list format)
                FFTs <- SD[, I(lapply(.SD, \(x) fft(na.omit(x)) / .N)), .SDcols = flux_variables]

                # ### testing auto-covariances
                # AutoCovars <- lapply(names(FFTs), \(i) {
                #     # check scalars
                #     if (i %in% scalars && !is.null(na_ind <- na.action(FFTs[[i]]))) {
                #         # fix lengths
                #         FFTs[[i]] <- FFTs[[i]][-na_ind]
                #         # get N
                #         N <- length(FFTs[[i]])
                #     } else {
                #         N <- .N
                #     }
                #     # get Re
                #     re <- Re(fft(Conj(FFTs[[i]]) * FFTs[[i]], inverse = TRUE))
                #     # get missing
                #     # subset
                #     if (N %% 2) {
                #         out <- re[c(((N + 1) / 2 + 1):N, 1:((N + 1) / 2))] * N / (N - 1)
                #     } else {
                #         out <- re[c((N / 2 + 1):N, 1:(N / 2))] * N / (N - 1)
                #     }
                #     n_missing <- n_period - N
                #     if (sign(n_missing) >= 0) {
                #         if (n_missing %% 2) {
                #             n1 <- rep(NA_real_, (n_missing + 1) / 2)
                #             n2 <- rep(NA_real_, (n_missing - 1) / 2)
                #         } else {
                #             n1 <- n2 <- rep(NA_real_, n_missing / 2)
                #         }
                #         c(n1, out, n2)
                #     } else {
                #         if (-n_missing %% 2) {
                #             ind <- ((1 - n_missing) / 2 + 1):(N - (-n_missing - 1) / 2)
                #         } else {
                #             ind <- (1 - n_missing / 2):(N + n_missing / 2)
                #         }
                #         out[ind]
                #     }
                # })
                # names(AutoCovars) <- paste(flux_variables, flux_variables, sep = 'x')

                # calculate covariances with fix lag time:
                # -------------------------------------------------------------------------- 
                cat("\t- covariances\n")
                Covars <- lapply(covariances_variables, function(i, ffts) {
                    # check scalars
                    if (i[2] %in% scalars && !is.null(na_ind <- na.action(ffts[[i[2]]]))) {
                        # fix lengths
                        ffts[[i[1]]] <- ffts[[i[1]]][-na_ind]
                        # get N
                        N <- length(ffts[[i[1]]])
                    } else {
                        N <- .N
                    }
                    # get Re
                    re <- Re(fft(Conj(ffts[[i[2]]]) * ffts[[i[1]]], inverse = TRUE))
                    # get missing
                    # subset
                    if (N %% 2) {
                        out <- re[c(((N + 1) / 2 + 1):N, 1:((N + 1) / 2))] * N / (N - 1)
                    } else {
                        out <- re[c((N / 2 + 1):N, 1:(N / 2))] * N / (N - 1)
                    }
                    n_missing <- n_period - N
                    if (sign(n_missing) >= 0) {
                        if (n_missing %% 2) {
                            n1 <- rep(NA_real_, (n_missing + 1) / 2)
                            n2 <- rep(NA_real_, (n_missing - 1) / 2)
                        } else {
                            n1 <- n2 <- rep(NA_real_, n_missing / 2)
                        }
                        c(n1, out, n2)
                    } else {
                        if (-n_missing %% 2) {
                            ind <- ((1 - n_missing) / 2 + 1):(N - (-n_missing - 1) / 2)
                        } else {
                            ind <- (1 - n_missing / 2):(N + n_missing / 2)
                        }
                        out[ind]
                    }
                }, ffts = FFTs)
                names(Covars) <- covariances

                # find maximum in dynamic lag time range:
                # -------------------------------------------------------------------------- 
                cat("\t- dyn lag\n")
                dyn_lag_max <- sapply(covariances, function(i, x, lag) {
                    find_dynlag(x[[i]], lag[, i])
                }, x = Covars, lag = dyn_lag)

                # covariance function's standard deviation and mean values left and right of fix lag
                # ------------------------------------------------------------------------
                # RE_RMSE (Eq. 9 in Langford et al. 2015)
                # -> ranges lo/hi (-/+180 to -/+150 secs? => define range as argument)
                # gamma_time_window <- c(2.5, 5)
                n <- length(Covars[[1]])
                m <- ifelse(n %% 2, (n + 1) / 2, n / 2 + 1)
                lo_range <- m - rev(gamma_time_window) * 60 * .Hz
                hi_range <- m + gamma_time_window * 60 * .Hz
                # -> sd_cov_low
                sd_cov_lo <- sapply(Covars, \(x) sd(x[lo_range]))
                # -> avg_cov_low
                avg_cov_lo <- sapply(Covars, \(x) mean(x[lo_range]))
                # -> sd_cov_hi
                sd_cov_hi <- sapply(Covars, \(x) sd(x[hi_range]))
                # -> avg_cov_hi
                avg_cov_hi <- sapply(Covars, \(x) mean(x[hi_range]))
                re_rmse <- sqrt(0.5 * (sd_cov_lo ^ 2 + avg_cov_lo ^ 2 +
                        sd_cov_hi ^ 2 + avg_cov_hi ^ 2))
                # ## auto-covariances
                # sd_acov_lo <- sapply(AutoCovars, \(x) sd(x[lo_range]))
                # # avg_acov_lo <- sapply(AutoCovars, \(x) mean(x[lo_range]))
                # sd_acov_hi <- sapply(AutoCovars, \(x) sd(x[hi_range]))
                # # avg_acov_hi <- sapply(AutoCovars, \(x) mean(x[hi_range]))
                # # re_rmse_acov <- sqrt(0.5 * (sd_acov_lo ^ 2 + avg_acov_lo ^ 2 +
                # #         sd_acov_hi ^ 2 + avg_acov_hi ^ 2))
                # re_rmse_acov <- sqrt(0.5 * (sd_acov_lo ^ 2 + sd_acov_hi ^ 2))

                # covariance function values +/- tau.off.sec of fix lag
                # ------------------------------------------------------------------------
                # not implemented...

                # cospectra for fixed & dynamic lags
                # ------------------------------------------------------------------------ 			
                cat("\t- co-spectra\n")
                Cospec_fix <- mapply(function(i, lag, ffts) {
                        # check scalars
                        N <- .N
                        if (i[2] %in% scalars) {
                            if (!is.null(na_ind <- na.action(scalar_list[[i[2]]]))) {
                                # fix lengths
                                ffts[[i[1]]] <- ffts[[i[1]]][-na_ind]
                                # get N
                                N <- length(ffts[[i[1]]])
                            }
                            x <- scalar_list[[i[2]]]
                        } else {
                            x <- SD[, get(i[2])]
                        }
                        xs <- fft(shift(x, -lag)) / N
                        re <- Re(Conj(xs) * ffts[[i[1]]])[seq(N / 2) + 1] * N / 
                            (N - 1) * 2
                        # get missing
                        n_missing <- n_period / 2 - length(re)
                        if (sign(n_missing) >= 0) {
                            c(re, rep(NA_real_, n_missing))
                        } else {
                            re[seq_len(n_period / 2)]
                        }
                    }, i = covariances_variables, lag = fix_lag, 
                    MoreArgs = list(ffts = FFTs), SIMPLIFY = FALSE)
                names(Cospec_fix) <- covariances
                Cospec_dyn <- mapply(function(i, lag, ffts) {
                        # check scalars
                        N <- .N
                        if (i[2] %in% scalars) {
                            if (!is.null(na_ind <- na.action(scalar_list[[i[2]]]))) {
                                # fix lengths
                                ffts[[i[1]]] <- ffts[[i[1]]][-na_ind]
                                # get N
                                N <- length(ffts[[i[1]]])
                            }
                            x <- scalar_list[[i[2]]]
                        } else {
                            x <- SD[, get(i[2])]
                        }
                        xs <- fft(shift(x, -lag)) / N
                        re <- Re(Conj(xs) * ffts[[i[1]]])[seq(N / 2) + 1] * N / 
                            (N - 1) * 2
                        # get missing
                        n_missing <- n_period / 2 - length(re)
                        if (sign(n_missing) >= 0) {
                            c(re, rep(NA_real_, n_missing))
                        } else {
                            re[seq_len(n_period / 2)]
                        }
                    }, i = covariances_variables, lag = dyn_lag_max[2, ],
                    MoreArgs = list(ffts = FFTs), SIMPLIFY = FALSE)
                names(Cospec_dyn) <- covariances
        
                # ogives for fixed & dynamic lags 
                # ------------------------------------------------------------------------ 
                Ogive_fix <- lapply(Cospec_fix, function(x) {
                    x[is.na(x)] <- 0
                    rev(cumsum(rev(x)))
                })
                names(Ogive_fix) <- covariances
                Ogive_dyn <- lapply(Cospec_dyn, function(x) {
                    x[is.na(x)] <- 0
                    rev(cumsum(rev(x)))
                })
                names(Ogive_dyn) <- covariances
                if (ogives_out) {
                    int_start <- format(Int_Start)
                    e_ogive$Cospec_fix_Out[[dates[day]]][[int_start]] <- c(list(freq = freq), Cospec_fix)
                    e_ogive$Cospec_dyn_Out[[dates[day]]][[int_start]] <- c(list(freq = freq), Cospec_dyn)
                    e_ogive$Ogive_fix_Out[[dates[day]]][[int_start]] <- c(list(freq = freq), Ogive_fix)
                    e_ogive$Ogive_dyn_Out[[dates[day]]][[int_start]] <- c(list(freq = freq), Ogive_dyn)
                    e_ogive$Covars_Out[[dates[day]]][[int_start]] <- c(list(Hz = .Hz), Covars)
                }


                # empirical damping estimation, dyn and fix should have best reference (dyn/fix)...
                # ------------------------------------------------------------------------
                if (any(scalar_covariances)) {
                    cat("\t- damping\n")
                    Damping_fix <- mapply(damp_hac5, ogive = Ogive_fix[scalar_covariances], 
                        ogive_ref = Ogive_fix[damping_reference], freq.limits = damp_region,
                        MoreArgs = list(freq = freq), SIMPLIFY = FALSE)
                    Damping_dyn <- mapply(damp_hac5, ogive = Ogive_dyn[scalar_covariances], 
                        ogive_ref = Ogive_dyn[damping_reference], freq.limits = damp_region,
                        MoreArgs = list(freq = freq), SIMPLIFY = FALSE)
                } else {
                    Damping_dyn <- Damping_fix <- NULL
                }

                ##### ~~~> sub-intervals switched off
                # if (FALSE) {
                #     # sub-int: sub-interval calculations (switch do data.frame since code already exists (I'm lazy))
                #     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                #     cat("~~~\nsub-intervals (subint_n = ", subint_n, ")\n", sep = '')
                #     sub_indices <- split_index(.N, subint_n)
                #     sub_Data <- lapply(sub_indices, function(i, x) x[i, ], x = as.data.frame(SD))

                #     # sub-int: calculate wind direction, rotate u, v, w, possibly detrend T (+ u,v,w)
                #     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                #     sub_wind <- lapply(sub_Data, function(x) 
                #         rotate_detrend(x[, "u"], x[, "v"], x[, "w"], x[, "T"], 
                #             method = subint_detrending[c("u", "v", "w", "T")], Hz_ts = .Hz))
                #     for(sub in seq_along(sub_wind)){
                #         sub_Data[[sub]][,c("u","v","w","T")] <- sub_wind[[sub]][c("uprot","vprot","wprot","Tdet")]
                #     }

                #     # sub-int: detrend scalars
                #     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                #     if (length(scalars)) {
                #         sub_detrended_scalars <- lapply(sub_Data, 
                #             function(x) mapply(trend, y = x[, scalars, drop = FALSE], 
                #                 method = subint_detrending[scalars], 
                #                 MoreArgs = list(Hz_ts = .Hz), SIMPLIFY = FALSE))
                #         for (sub in seq_along(sub_wind)) {
                #             sub_Data[[sub]][, scalars] <- lapply(sub_detrended_scalars[[sub]], "[[", "residuals")
                #         }
                #     }

                #     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                    
                #     # extract covariances etc.:
                #     # ------------------------------------------------------------------------ 
                #     # sub-int: calculate covariances of sub-intervals
                #     # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                #     sub_covs <- sapply(sub_wind, function(x) cov(list2DF(x[1:4]))[
                #         cbind(
                #             c("uprot", "vprot", "wprot", "Tdet", "uprot", "uprot", "uprot", "vprot", "vprot", "wprot"), 
                #             c("uprot", "vprot", "wprot", "Tdet", "vprot", "wprot", "Tdet", "wprot", "Tdet", "Tdet")
                #             )]
                #     )
                #     sub_cov_means <- apply(sub_covs, 1, mean)
                #     sub_cov_sd <- apply(sub_covs, 1, sd)
                #     names(sub_cov_means) <- paste0("mean(sub.int.", c("<u'u'>", "<v'v'>", "<w'w'>", "<T'T'>", "<u'v'>", "<u'w'>", "<u'T'>", "<v'w'>", "<v'T'>", "<w'T'>"), ")")
                #     names(sub_cov_sd) <- paste0("sd(sub.int.", c("<u'u'>", "<v'v'>", "<w'w'>", "<T'T'>", "<u'v'>", "<u'w'>", "<u'T'>", "<v'w'>", "<v'T'>", "<w'T'>"), ")")
                # }

                ### some output and namings
                fix_lag_out <- fix_lag
                dyn_lag_out <- dyn_lag_max["tau", ]
                flux_fix_lag <- sapply(Ogive_fix, "[", 1)
                flux_dyn_lag <- sapply(Ogive_dyn, "[", 1)
                fix_damping_pbreg <- sapply(Damping_fix, "[[", 1)
                dyn_damping_pbreg <- sapply(Damping_dyn, "[[", 1)
                fix_damping_deming <- sapply(Damping_fix, "[[", 2)
                dyn_damping_deming <- sapply(Damping_dyn, "[[", 2)
                if (!is.null(Damping_fix)) {
                    names(fix_damping_pbreg) <- names(Damping_fix)
                    names(fix_damping_deming) <- names(Damping_fix)
                }
                if (!is.null(Damping_dyn)) {
                    names(dyn_damping_pbreg) <- names(Damping_dyn)
                    names(dyn_damping_deming) <- names(Damping_dyn)
                }

                # write results:
                # -------------------------------------------------------------------------- 

                out <- c(
                    list(
                        st = Int_Start
                        , et = et_interval[.BY[[1]]]
                        , n_values = .N
                        #, SubInts =  subint_n
                    )
                    , if (ht_provided) {
                        list(
                            ht_temp_amb = mean(ht_temp_amb, na.rm = TRUE)
                            , ht_press_amb = mean(ht_press_amb, na.rm = TRUE)
                            , ht_oss = mean(ht_oss, na.rm = TRUE)
                            , ht_alarm_codes = paste(unique(unlist(strsplit(
                                            unique(ht_alarm_code[!is.na(ht_alarm_code)])
                                            , split = ','))), collapse = ',')
                        )
                    }
                    , if (licor_provided) {
                        list(
                            li_temp_amb = mean(li_temp_amb, na.rm = TRUE)
                            , li_press_amb = mean(li_press_amb, na.rm = TRUE)
                            , li_co2ss = mean(li_co2ss, na.rm = TRUE)
                        )
                    }
                    , as.list(c(
                        wind_stats
                        # scalar avg & sd
                        , setNames(scalar_means[input_scalars], paste0('avg_', input_scalars))
                        , setNames(scalar_sd[input_scalars], paste0('sd_', input_scalars))
                        # fluxes
                        , setNames(
                            flux_fix_lag[input_covariances],
                            paste0('flux_fix_', input_covariances)
                        )
                        , setNames(
                            flux_dyn_lag[input_covariances],
                            paste0('flux_dyn_', input_covariances)
                        )
                        # , sub_cov_means
                        # , sub_cov_sd
                        # fix/dyn lag as seconds
                        , setNames(
                            fix_lag_out[input_covariances] / .Hz, 
                            paste0('lag_fix_', input_covariances)
                        )
                        , setNames(
                            dyn_lag_out[input_covariances] / .Hz, 
                            paste0('lag_dyn_', input_covariances)
                        )
                        , setNames(
                            re_rmse[input_covariances],
                            paste0('re_rmse_', input_covariances)
                        )
                        , if (!is.null(Damping_fix)) {
                            c(
                                setNames(
                                    fix_damping_pbreg[scalar_covariances_only],
                                    paste0('damping_fix_pbreg_', scalar_covariances_only)
                                ),
                                setNames(
                                    fix_damping_deming[scalar_covariances_only],
                                    paste0('damping_fix_deming_', scalar_covariances_only)
                                )
                            )
                        }
                        , if (!is.null(Damping_dyn)) {
                            c(
                                setNames(
                                    dyn_damping_pbreg[scalar_covariances_only],
                                    paste0('damping_dyn_pbreg_', scalar_covariances_only)
                                ),
                                setNames(
                                    dyn_damping_deming[scalar_covariances_only],
                                    paste0('damping_dyn_deming_', scalar_covariances_only)
                                )
                            )
                        }
                        ))
                )				


                if (create_graphs) {
                    if (!dir.exists(path_folder)) {
                        dir.create(path_folder, recursive = FALSE)
                    }
                    # get end of current interval time in UTC
                    soi_utc <- st_interval[.GRP]
                    eoi_utc <- et_interval[.GRP]
                    # get date in correct format
                    date_formatted <- gsub('-', '', dates[day], fixed = TRUE)
                    # plotting:
                    # -------------------------------------------------------------------------- 
                    # -------------------------------------------------------------------------- 
                    cat("~~~\nplotting timeseries and fluxes\n")
                    # time series:
                    # -------------------------------------------------------------------------- 
                    # plot and save (rotated) data time series with raw-data trends...
                    # ------------------------------------------------------------------------
                    ## TODO: -> fix tz!!! -> use UTC but indicate in name!!!
                    time2 <- format(soi_utc, format = "%H%M")
                    plotname <- paste("timeseries", date_formatted, time2, sep="-") 
                    ts_vars <- names(plot_timeseries)[plot_timeseries]
                    jpeg(file = paste0(path_folder, '/', plotname, ".jpg"), width = 600, height = (sum(plot_timeseries)) * 100, quality = 60)
                        ts_plot <- plot.tseries(
                            cbind(st = Time, as.data.frame(SD)),
                            wind, detrended_scalars, ts_vars,
                            plotting_var_colors, plotting_var_units)
                        print(ts_plot)
                    dev.off()

                    # plot and save flux evaluation...
                    # ------------------------------------------------------------------------
                    for(i in covariances){
                        # i <- "w'TDL CH4'"
                        plotname <- paste("plots", date_formatted, time2, covariances_plotnames[i], sep = "-")
                        # fix ylab
                        ylab <- sub('(.+)x(.+)', "<\\1'\\2'>", i)
                        jpeg(file=paste0(path_folder, '/', plotname,".jpg"),width=1350, height=900, quality=60)
                            par(mfrow=c(2,3))
                            # ----------------------- Covariance -------------------------------------
                            plot_covfunc(Covars[[i]],n_period / .Hz,dyn_lag_max[,i],fix_lag[i],ylab=ylab, xlim = c(-50,50), cx=1.5, cxmt=1.25, cl=plotting_covar_colors[i])
                            # ---------------------- Co-Spec/Ogive fix lag -----------------------------------
                            plot_cospec_ogive(Ogive_fix[[i]],Cospec_fix[[i]],freq,ylab=paste0("ogive (fix lag) of ",ylab),cx=1.5,col=plotting_covar_colors[i])
                            # ---------------------- Co-Spec/Ogive dyn lag -----------------------------------
                            plot_cospec_ogive(Ogive_dyn[[i]],Cospec_dyn[[i]],freq,ylab=paste0("ogive (dyn lag) of ",ylab),cx=1.5,col=plotting_covar_colors[i])
                            # ---------------------- empirical damping -----------------------------------
                            if(scalar_covariances[i]){
                                plot_damping(Damping_fix[[i]],freq,ylab=paste0("ogive (fix lag) of ",i),cx=1.5,col=plotting_covar_colors[i])
                                plot_damping(Damping_dyn[[i]],freq,ylab=paste0("ogive (dyn lag) of ",i),cx=1.5,col=plotting_covar_colors[i])
                            }
                            title(paste0(ylab, " flux ", 
                                format(soi_utc, format = "(%H:%M:%S"), " - ", format(eoi_utc, format = "%H:%M:%S)"), 
                                if (scalar_covariances[i]) {
                                    reflab <- sub('(.+)x(.+)', "<\\1'\\2'>", damping_reference[i])
                                    paste0(" - damping reference flux: ", reflab)
                                }
                            ), outer = TRUE, line = -1)
                        dev.off()	
                    }
                } # end plotting
                # return out
                out
            } else {
                cat("less than ", round(thresh_period * 100),"% of data points in raw data - skipping interval.\n", sep = '')
                # return NULL
                NULL
            }
        }, by = bin]

    } # end for loop over days

    # rbind list to data.table
    results <- rbindlist(result_list, fill = TRUE)

	cat("************************************************************\n") 
	cat("operation finished @", format(Sys.time(), "%d.%m.%Y %H:%M:%S"), 
        "time elapsed: ", difftime(Sys.time(), script.start, unit = "mins"),
        "minutes\n")
	cat("************************************************************\n")  

    if (nrow(results) == 0) {
        return(NULL)
    }

    # remove bin column
    results[, bin := NULL]

    # fix time zone
    results[, ':='(
        st = with_tz(st, tz_user),
        et = with_tz(et, tz_user)
    )]

    # convert to ibts
    if (as_ibts) {
        results <- as.ibts(results)
    }

	# #################################### END VERSION HISTORY #################################### #

	if (ogives_out) {
		structure(
            results, 
            covars = e_ogive$Covars_Out,
            cospec_fix = e_ogive$Cospec_fix_Out, 
            cospec_dyn = e_ogive$Cospec_dyn_Out,
            ogv_fix = e_ogive$Ogive_fix_Out, 
            ogv_dyn = e_ogive$Ogive_dyn_Out
        )
	} else {
        results
	}
}

evalREddy <- function(
		config_file = NULL
		,file_directory = NULL
		,save_directory = NULL
		,add_name = NULL
		,start_time = NULL
		,end_time = NULL
		,avg_period = NULL                                                                                     
		,variables = NULL
		,covariances = NULL
		,ogives_out = FALSE
		,z_sonic = NULL
		,dev_north = NULL
		,tz_sonic = NULL
		,z_canopy = NULL
		,create_graphs=TRUE
        ,as_ibts = TRUE
	){

    library(tcltk)
    library(xlsx)
    # library(Rcpp)
    library(reshape2)
    library(lattice)
    # library(caTools)

	script.start <- Sys.time()
	################################################################################
	# ----------------------------- read config file -------------------------------
	doCalc <- TRUE
	if(is.null(config_file)){
		config_file <- choose.files(caption = "select config file",filters=matrix(c("Excel file (*.xlsx)","*.xlsx"),nrow=1))
		if(length(config_file)==0){
			doCalc <- FALSE
		}
	}
	if(doCalc){
		# check number of columns:
		check_cols <- read.xlsx(config_file,sheetName="REddy_config",rowIndex=c(20,32),check.names=FALSE,stringsAsFactors=FALSE,header=FALSE)
		cols <- 1:(max(rowSums(!is.na(check_cols[,-(1:6)])))+6)
		config <- read.xlsx2(config_file,sheetName="REddy_config",check.names=FALSE,stringsAsFactors=FALSE,colIndex=cols)
		config_vars <- config[,"var-name in R-script"]
		if(is.null(file_directory)) file_directory <- config[config_vars %in% "file_directory",7]
		if(is.null(save_directory)) save_directory <- config[config_vars %in% "save_directory",7]
		if(is.null(add_name)) add_name <- config[config_vars %in% "add_name",7]
		if(is.null(avg_period)) avg_period <- as.numeric(config[config_vars %in% "avg_period",7])
		if(is.null(start_time)) start_time <- config[config_vars %in% "start_time",7]
		if(is.null(end_time)) end_time <- config[config_vars %in% "end_time",7]
		if(is.null(tz_sonic)) tz_sonic <- config[config_vars %in% "tz_sonic",7]
		if(is.null(z_sonic)) z_sonic <- as.numeric(config[config_vars %in% "z_sonic",7])
		if(is.null(dev_north)) dev_north <- as.numeric(config[config_vars %in% "dev_north",7])
		if(is.null(z_canopy)) suppressWarnings(z_canopy <- as.numeric(config[config_vars %in% "z_canopy",7]))
		if(is.null(create_graphs)) create_graphs <- as.logical(config[config_vars %in% "create_graphs",7])
		sonic_type <- config[config_vars %in% "sonic_type",7]
		rotation_method <- as.character(config[config_vars %in% "rotation_method",7])
		rawdata_function <- config[config_vars %in% "rawdata_function",7]
		fixlag_lookup_file <- as.character(config[config_vars %in% "fixlag_lookup_file",7])
		n_stationarity_subint <- as.numeric(config[config_vars %in% "n_stationarity_subint",7])



		if(v_flag <- !is.null(variables)){
			variables_temp <- variables
		}
		variables <- as.character(config[config_vars %in% "variables",-(1:6)]) %w/o% ""
		variables_units <- as.character(config[config_vars %in% "variables_units",seq.int(len_variables)+6])
		variables_cols <- as.character(config[config_vars %in% "variables_cols",seq.int(len_variables)+6])
		variables_timeseries <- as.logical(config[config_vars %in% "variables_timeseries",seq.int(len_variables)+6])
		hard_flag <- as.logical(config[config_vars %in% "hard_flag",seq.int(len_variables)+6])
		detrending_method <- as.character(config[config_vars %in% "detrending_method",seq.int(len_variables)+6])
		subint_detrending_method <- as.character(config[config_vars %in% "subint_detrending_method",seq.int(len_variables)+6])
		hf_lo <- suppressWarnings(as.numeric(config[config_vars %in% "hf_lo",seq.int(len_variables)+6]))
		hf_hi <- suppressWarnings(as.numeric(config[config_vars %in% "hf_hi",seq.int(len_variables)+6]))
		hf_window <- as.numeric(config[config_vars %in% "hf_window",7])
		hf_method <- as.character(config[config_vars %in% "hf_method",7])
		names(subint_detrending_method) <-names(detrending_method) <- names(variables_units) <- names(variables_cols) <- names(variables_timeseries) <- names(hard_flag) <- names(hf_lo) <- names(hf_hi) <- variables
		if(v_flag){
			variables <- variables_temp
			variables_units <- variables_units[variables]
			variables_cols <- variables_cols[variables]
			variables_timeseries <- variables_timeseries[variables]
			hard_flag <- hard_flag[variables]
			hf_lo <- hf_lo[variables]
			hf_hi <- hf_hi[variables]
			detrending_method <- detrending_method[variables]
			subint_detrending_method <- subint_detrending_method[variables]
		}
		scalars <- variables %w/o% c("u","v","w","T")
		lim_range <- rbind(lower=hf_lo,top=hf_hi)

		if(cov_flag <- !is.null(covariances)){
			covariances_temp <- covariances
		}
		covariances <- as.character(config[config_vars %in% "covariances",-(1:6)]) %w/o% ""
		len_covariances <- length(covariances)
		covariances_cols <- as.character(config[config_vars %in% "covariances_cols",seq.int(len_covariances)+6])
		lag_fix <- as.numeric(config[config_vars %in% "lag_fix",seq.int(len_covariances)+6])
		lag_dev <- as.numeric(config[config_vars %in% "lag_dev",seq.int(len_covariances)+6])
		damp_ref <- as.character(config[config_vars %in% "damp_ref",seq.int(len_covariances)+6])
		damp_lower <- suppressWarnings(as.numeric(config[config_vars %in% "damp_lower",seq.int(len_covariances)+6]))
		damp_upper <- suppressWarnings(as.numeric(config[config_vars %in% "damp_upper",seq.int(len_covariances)+6]))
		damp_region <- mapply(c,damp_lower,damp_upper,SIMPLIFY=FALSE)
		names(damp_ref) <- names(covariances_cols) <- names(lag_fix) <- names(lag_dev) <- names(damp_region) <- covariances
		if(cov_flag){
			covariances <- covariances_temp
			covariances_cols <- covariances_cols[covariances]
			lag_fix <- lag_fix[covariances]
			lag_dev <- lag_dev[covariances]
			damp_region <- damp_region[covariances]
		}
		covariances_variables <- strsplit(covariances,"'")
		covariances_plotnames <- make.names(gsub("'","",covariances))
		scalar_covariances <- as.numeric(grepl("(^u'|'u'$)",covariances)) + as.numeric(grepl("(^w'|'w'$)",covariances)) + as.numeric(grepl("(^T'|'T'$)",covariances)) < 2
		reference_covariances <- rep(NA_character_,length(covariances))
		ref_cov <- strsplit(covariances[scalar_covariances],"'")
		ref_scalar <- lapply(ref_cov,function(x)!grepl("(^u$|^w$|^T$)",x))
		damp_ref[damp_ref %in% ""] <- "T"
		damp_ref <- paste0(damp_ref,"'")
		reference_covariances[scalar_covariances] <- mapply(function(x,y,z){
			x[y] <- z 
			paste(x,collapse="'")
			},ref_cov,ref_scalar,damp_ref[scalar_covariances])
		names(covariances_variables) <- names(covariances_plotnames) <- names(scalar_covariances) <- names(reference_covariances) <- covariances

		# ------------------------------------------------------------------------------

		# get rawdata function:                                      
		# ------------------------------------------------------------------------------ 
		if(!exists(rawdata_function,mode="function")){
			doCalc <- FALSE
			stop("rawdata_function doesn't exist")
		}
		# start reading and processing single data files:                                      
		# ------------------------------------------------------------------------------ 
		if(doCalc){


			# start of code
			# ==============================================================================
			# ==============================================================================
			# request data directory/filename
			# ------------------------------------------------------------------------------
			if(file.info(file_directory)$isdir){
				dfname <- tclvalue(tkgetOpenFile(initialdir=file_directory,title="select raw data file",multiple=FALSE))
			} else {
				dfname <- file_directory
			}
			now <- Sys.time()
			now1 <- format(now,"%d.%m.%Y %H:%M:%S")
			now0 <- format(now,"%Y%m%d_%H%M")
			filename <- sub("([.]{1}.*)","",basename(dfname[1]))


			# create result folder (folder name includes first input-filename) and set this directory as working directory                                      
			# ------------------------------------------------------------------------------ 
			former_wd <- getwd()
			if(isFALSE(save_directory)){
				create_graphs <- FALSE
			} else {					
				if(add_name!=""){
					folder <- paste0("REddy-",filename,"-",add_name,"-eval",now0,"-avg",avg_period)
				} else {
					folder <- paste0("REddy-",filename,"-eval",now0,"-avg",avg_period)
				}
                path_folder <- file.path(save_directory, folder)
				dir.create(path_folder, recursive = FALSE)
				# on.exit(setwd(former_wd))
				# setwd(paste0(save_directory,"/",folder))
			}

			read_rawdata <- get(rawdata_function,mode="function")

			# reading data file:                                      
			# ------------------------------------------------------------------------------
			cat("\n************************************************************\n")
            cat("EC flux processing\n")
			cat("************************************************************\n")
            cat('\nReading raw data from file "', dfname, '"\n', sep = '')
			# browser() 
			Data <- read_rawdata(dfname)
            # check if empty
            if (nrow(Data) == 0) {
                cat('Empty file...\n')
                return(NULL)
            }
			if(!all(cn_check <- variables %in% names(Data))){
				stop("cannot find variable(s): ",paste(variables[!cn_check],collapse=" ,"),"\nAvailable column names are: ",paste(names(Data),collapse=" ,"))
			}
			# check if delta_t exists:
			if(names(Data)[2] %in% c("u","v","w","T")){
				delta_t <- as.numeric(diff(Data[,1]))*1000
				Data <- cbind(Data[,1],delta_t=c(median(delta_t),delta_t),Data[,-1])
			}
			Data <- Data[,c(1:2,which(names(Data) %in% variables))]
            
			# correct time zone:                                      
			# ------------------------------------------------------------------------------ 
			if(!identical(tz_sonic,tz(Data[,1]))){
				tz(Data[,1]) <- tz_sonic
			}
			# enforce POSIXct
			if(!is.POSIXct(Data[,1]))Data[,1] <- as.POSIXct(Data[,1])

			# browser()
			# get intervals:                                      
			# ------------------------------------------------------------------------------ 
			Data_Day <- format(Data[1,1],"%Y.%m.%d")
			# start
			if(is.character(start_time)){
				if(grepl("[:]",start_time[1])){
					Start <- ibts::parse_date_time3(paste(Data_Day,start_time),tz=tz_sonic)
				} else {
					Start <- switch(start_time
						,"full_interval" =
						, "interval" = {
							is <- ibts::parse_date_time3(paste(Data_Day,"00:00:00"),tz=tz_sonic) + seq(0,24*3600,avg_period*60)
							is[which(is >= Data[1,1])[1]]
						}
						,"start" = 
						,"beginning" = Data[1,1]
						)
				}
			} else {
				# check POSIX
				if(is.POSIXlt(start_time)){
					Start <- as.POSIXct(start_time)
				} else if(is.POSIXct(start_time)){
					Start <- start_time
				} else {
					stop("invalid argument type: 'start_time'")
				}
			}
			# end
			if(is.character(end_time)){
				if(grepl("[:]",end_time[1])){
					End <- ibts::parse_date_time3(paste(Data_Day,end_time),tz=tz_sonic)
					if(any(add1day <- end_time %in% c("00:00:00","00:00"))){
						End[add1day] <- End[add1day] + as.difftime(1,units="days")
					}
				} else {
					End <- switch(end_time
						,"full_interval" =
						, "interval" = {
							is <- ibts::parse_date_time3(paste(Data_Day,"00:00:00"),tz=tz_sonic) + seq(0,24*3600,avg_period*60)
							is[rev(which(is <= Data[nrow(Data),1]))[1]]
						}
						,"all" = 
						,"end" = Data[nrow(Data),1]
						)
				}
			} else {
				# check POSIX
				if(is.POSIXlt(end_time)){
					End <- as.POSIXct(end_time)
				} else if(is.POSIXct(end_time)){
					End <- end_time
				} else {
					stop("invalid argument type: 'end_time'")
				}
			}

			if(length(Start)==1){
				if(End-avg_period*60 < Start){
					warning(paste0("File '",filename,"' (",format(Data_Day),"): Too few data to evaluate! Returning 'NULL'"))
					return(NULL)
				}
				st_interval <- seq(Start,End-avg_period*60,avg_period*60)
				et_interval <- st_interval + avg_period*60
				cat("\n~~~\nCalculation will include all intervals between",format(st_interval[1],format="%Y-%m-%d %H:%M"),"and",format(et_interval[length(et_interval)],format="%Y-%m-%d %H:%M",usetz=TRUE),"on a",avg_period,"min basis\n~~~\n\n")
			} else {
				st_interval <- Start
				et_interval <- End
				cat("\n~~~\nCalculation will include the following intervals:\n",paste(format(st_interval,format="   %Y-%m-%d %H:%M"),"to",format(et_interval,format="%Y-%m-%d %H:%M",usetz=TRUE)), '\n~~~\n')
			}

			# prepare results:                                      
			# ------------------------------------------------------------------------------ 
			results <- cbind(
				data.frame(
					st=st_interval
					,et=et_interval
					,Data_Start = with_tz(NA_POSIXct_, tz_sonic)
					,Data_End = with_tz(NA_POSIXct_, tz_sonic)
					,tzone = tz_sonic
					,"Int_Length (min)" = avg_period
					,stringsAsFactors=FALSE,check.names=FALSE
					)
				,
				matrix(NA,ncol=45 + 4*length(covariances) + 4*sum(scalar_covariances) + 3*length(scalars),nrow=length(st_interval))
				)
			checkNames <- TRUE

			# get flux variables and lag times:                                      
			# ------------------------------------------------------------------------------ 
			flux_variables <- unique(unlist(strsplit(covariances,"'")))
			# --------------------- read fix lags from lag_lookuptable ---------------------
			cov_names <- make.names(gsub("'","",covariances))
			names(cov_names) <- covariances
			if(fixlag_lookup_file!=""){
				df.lag.fix.sec <- read.table(fixlag_lookup_file, header=TRUE, sep=";", as.is=TRUE, fill=TRUE)
				table_time <- fast_strptime(df.lag.fix.sec[,1], "%d.%m.%Y",tz="Etc/GMT-1",lt=FALSE)
				rownames(df.lag.fix.sec) <- format(table_time, "%y%m%d")
				cov_exist <- cov_names %in% names(df.lag.fix.sec)
				# select for specified time range:
				table_ind <- which(table_time >= trunc(Start[1],"days") & table_time <= trunc(End[1],"days"))
				df.lag.fix.sec <- df.lag.fix.sec[table_ind,,drop=FALSE]
				# take existing fix lags from look up table:
				df.lag.fix.sec <- -df.lag.fix.sec[,cov_names[cov_exist],drop=FALSE]
				for(i in seq.int(ncol(df.lag.fix.sec))){
					df.lag.fix.sec[is.na(df.lag.fix.sec[,i]),i] <- as.numeric(lag_fix[covariances[cov_exist][i]])
					# what if this is not numeric!!!
				}
				names(df.lag.fix.sec) <- covariances[cov_exist]
			} else {
				cov_exist <- rep(TRUE,length(covariances))
				df.lag.fix.sec <- as.data.frame(as.list(lag_fix),check.names=FALSE)
				rownames(df.lag.fix.sec) <- format(Data[1,1],"%y%m%d")
			}
			# add non-existing fix lags to look up table:
			if(any(!cov_exist)){
				for(i in covariances[!cov_exist]){
					lag_add <- lag_fix[i]
					if(lag_add %in% names(df.lag.fix.sec)){
						lag_value <- df.lag.fix.sec[,lag_add]
					} else {
						lag_value <- as.numeric(lag_fix[i])
					}
					df.lag.fix.sec[,i] <- lag_value
				}
			}
			# sort (not really necessary, though):
			df.lag.fix.sec <- df.lag.fix.sec[,covariances]

			# find data in individual intervals:
			# ------------------------------------------------------------------------------ 
			if(any(st_interval[-1] - et_interval[-length(et_interval)] < 0))stop("Intervals need to be in strictly increasing order!")
			ind_interval <- getIntervals(Data[, 1], st_interval, et_interval)
			unique_ind <- unique(ind_interval) %w/o% 0

			if(ogives_out){
				Cospec_dyn_Out <- Cospec_fix_Out <- Covars_Out <- Ogive_fix_Out <- Ogive_dyn_Out <- vector("list",length(st_interval))
				names(Cospec_fix_Out) <- names(Cospec_dyn_Out) <- names(Covars_Out) <- names(Ogive_fix_Out) <- names(Ogive_dyn_Out) <- format(st_interval,format="%H%M")
			}

			# browser()

			# loop over individual intervals:
			# ------------------------------------------------------------------------------ 
			for(intval in unique_ind){
				
				# intval <- 1L
				# intval <- 65L
				# intval <- 1L
				# intval <- 2L
				# browser()
				cat("~~~~~~~~\ninterval",which(intval == unique_ind),"of",length(unique_ind),"\n")
				# get subset:
				# --------------------------------------------------------------------------
				index <- ind_interval %in% intval
				Int_length <- sum(index)
				if(Int_length > 10){
					Data_int <- Data[index,]

					# times, frequencies etc.
					# --------------------------------------------------------------------------
					Time <- Data_int[,1]
					Int_Start <- Time[1]
					Int_End <- Time[length(Time)] + Data_int[length(Time),2]/1000
					d_t <- diff(as.numeric(Time))
					Hz <- round(1/summary(d_t)[['Median']])
                    if (Hz < 10) {
                        cat('Frequency is lower than 10 Hz! Skipping interval.\n')
                        next
                    }
					freq <- Hz * seq(floor(Int_length / 2)) / floor(Int_length / 2)
					# for later: include NA where d_t>2*mean, remove entries where d_t < 0.5*mean

					# fix and dyn lags:                                      
					# ------------------------------------------------------------------------------ 
					fix_raw <- unlist(df.lag.fix.sec[format(Data[1,1],"%y%m%d"),])
					fix_lag <- round(fix_raw*Hz)
					dyn_lag <- rbind(
						lower=round((fix_raw - lag_dev)*Hz)
						,upper=round((fix_raw + lag_dev)*Hz)
						)
					# # check dyn lag <= 0 for scalars:
					# dyn_lag[2,scalar_covariances] <- pmin(dyn_lag[2,scalar_covariances],0)

					# raw data quality control I, i.e. hard flags = physical range
					# --------------------------------------------------------------------------
					cat("~~~\nchecking hard flags...\n")
					if(any(hard_flag)){
						Data_int[,variables[hard_flag]] <- H.flags(Data_int[,variables[hard_flag],drop=FALSE],Time,Hz,lim_range[,variables[hard_flag],drop=FALSE],hf_window,hf_method)
					} 

					# check for remaining NA values
					if(anyNA(Data_int[,c("u","v","w","T", scalars)])){
						cat("NA values in data. Skipping interval...")
						next
					}

					# calculate wind direction, rotate u, v, w, possibly detrend T (+ u,v,w)
					# -------------------------------------------------------------------------- 
					cat("~~~\nrotating and deterending sonic data...\n")
					wind <- rotate_detrend(Data_int[,"u"],Data_int[,"v"],Data_int[,"w"],Data_int[,"T"],phi=if(rotation_method %in% "two_axis") NULL else 0,method=detrending_method[c("u","v","w","T")],c.system=sonic_type,Hz_ts=Hz)
					Data_int[,c("u","v","w","T")] <- wind[c("uprot","vprot","wprot","Tdet")]

					# detrend scalars
					# -------------------------------------------------------------------------- 
					if(length(scalars)){
						cat("~~~\ndetrending scalars...\n")
						detrended_scalars <- mapply(trend,y=Data_int[,scalars, drop = FALSE],method=detrending_method[scalars],MoreArgs=list(Hz_ts=Hz),SIMPLIFY=FALSE)
						Data_int[,scalars] <- lapply(detrended_scalars,"[[","residuals")
					} else {
						detrended_scalars <- NULL
					}

					# start of flux relevant data manipulation
					# -------------------------------------------------------------------------- 
					# -------------------------------------------------------------------------- 
					cat("~~~\nstarting flux evaluation...\n")

					cov_vars <- unique(unlist(covariances_variables))
					FFTs <- lapply(Data_int[,cov_vars],function(x)fft(x)/Int_length)

					# calculate covariances with fix lag time:
					# -------------------------------------------------------------------------- 
					cat("\t- covariances\n")
					if(Int_length%%2){
						Covars <- lapply(covariances_variables,function(i,x){
							Re(fft(Conj(FFTs[[i[2]]]) * FFTs[[i[1]]], inverse=TRUE))[c(((Int_length+1)/2+1):Int_length,1:((Int_length+1)/2))]*Int_length/(Int_length-1)
							},x=FFTs)
					} else {
						Covars <- lapply(covariances_variables,function(i,x){
							Re(fft(Conj(FFTs[[i[2]]]) * FFTs[[i[1]]], inverse=TRUE))[c((Int_length/2+1):Int_length,1:(Int_length/2))]*Int_length/(Int_length-1)
							},x=FFTs)
					}
					names(Covars) <- covariances

					# find maximum in dynamic lag time range:
					# -------------------------------------------------------------------------- 
					cat("\t- dyn lag\n")
					# browser()
					dyn_lag_max <- sapply(covariances,function(i,x,lag){
						find_dynlag(x[[i]],lag[,i])
					},x=Covars,lag=dyn_lag)

					# covariance function's standard deviation and mean values left and right of fix lag
					# ------------------------------------------------------------------------
					# not implemented...

					# covariance function values +/- tau.off.sec of fix lag
					# ------------------------------------------------------------------------
					# not implemented...

					# cospectra for fixed & dynamic lags
					# ------------------------------------------------------------------------ 			
					cat("\t- co-spectra\n")
					Cospec_fix <- mapply(function(i,lag,x1,x2){
							xs <- fft(shift(x2[,i[2]],-lag))/Int_length
							Re(Conj(xs) * x1[[i[1]]])[seq(Int_length/2)+1]*Int_length/(Int_length-1)*2
						},i=covariances_variables,lag=fix_lag,MoreArgs=list(x1=FFTs,x2=Data_int),SIMPLIFY=FALSE)
					names(Cospec_fix) <- covariances
					# browser()
					Cospec_dyn <- mapply(function(i,lag,x1,x2){
							xs <- fft(shift(x2[,i[2]],-lag))/Int_length
							Re(Conj(xs) * x1[[i[1]]])[seq(Int_length/2)+1]*Int_length/(Int_length-1)*2
						},i=covariances_variables,lag=dyn_lag_max[2,],MoreArgs=list(x1=FFTs,x2=Data_int),SIMPLIFY=FALSE)
					names(Cospec_dyn) <- covariances
			
					
					#*** hac5: *** ogives for fixed & dynamic lags 
					# ------------------------------------------------------------------------ 
					Ogive_fix <- lapply(Cospec_fix,function(x)rev(cumsum(rev(x))))
					Ogive_dyn <- lapply(Cospec_dyn,function(x)rev(cumsum(rev(x))))		
					if(ogives_out){
						Cospec_fix_Out[[intval]] <- c(list(freq = freq), Cospec_fix)
						Cospec_dyn_Out[[intval]] <- c(list(freq = freq), Cospec_dyn)
						Ogive_fix_Out[[intval]] <- c(list(freq = freq), Ogive_fix)
						Ogive_dyn_Out[[intval]] <- c(list(freq = freq), Ogive_dyn)
						Covars_Out[[intval]] <- c(list(Hz = Hz), Covars)
					}


					# empirical damping (hac5), dyn and fix should have best reference (dyn/fix)...
					# ------------------------------------------------------------------------
					if(any(scalar_covariances)){
						cat("\t- damping\n")
						Damping_fix <- mapply(damp_hac5,ogive=Ogive_fix[scalar_covariances],ogive_ref=Ogive_fix[reference_covariances[scalar_covariances]],freq.limits = damp_region[scalar_covariances],MoreArgs = list(freq = freq),SIMPLIFY=FALSE)
						Damping_dyn <- mapply(damp_hac5,ogive=Ogive_dyn[scalar_covariances],ogive_ref=Ogive_dyn[reference_covariances[scalar_covariances]],freq.limits = damp_region[scalar_covariances],MoreArgs = list(freq = freq),SIMPLIFY=FALSE)
					} else {
						Damping_dyn <- Damping_fix <- NULL
					}

					# sub-int: sub-interval calculations
					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
					cat("~~~\nsub-intervals\n")
					sub_indices <- split_index(Int_length,n_stationarity_subint)
					sub_Data <- lapply(sub_indices,function(i,x)x[i,],x=Data[index,])

					# sub-int: calculate wind direction, rotate u, v, w, possibly detrend T (+ u,v,w)
					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
					sub_wind <- lapply(sub_Data,function(x)rotate_detrend(x[,"u"],x[,"v"],x[,"w"],x[,"T"],phi=if(rotation_method %in% "two_axis") NULL else 0,method=subint_detrending_method[c("u","v","w","T")],c.system=sonic_type,Hz_ts=Hz))
					for(sub in seq_along(sub_wind)){
						sub_Data[[sub]][,c("u","v","w","T")] <- sub_wind[[sub]][c("uprot","vprot","wprot","Tdet")]
					}

					# sub-int: detrend scalars
					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
					if(length(scalars)){
						sub_detrended_scalars <- lapply(sub_Data,function(x)mapply(trend,y=x[,scalars, drop = FALSE],method=subint_detrending_method[scalars],MoreArgs=list(Hz_ts=Hz),SIMPLIFY=FALSE))
						Data_int[,scalars] <- lapply(detrended_scalars,"[[","residuals")
						for(sub in seq_along(sub_wind)){
							sub_Data[[sub]][,scalars] <- lapply(sub_detrended_scalars[[sub]],"[[","residuals")
						}
					}

					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
					
					# extract covariances etc.:
					# ------------------------------------------------------------------------ 

					# calculate some turbulence parameters and collect some wind parameters
					# -------------------------------------------------------------------------- 
					wind_stats <- wind_statistics(wind,z_canopy,z_sonic)
					### correct for soni north deviation 
					wind_stats["WD"] <- (wind_stats["WD"] + dev_north) %% 360

					# calculate scalar averages and sd:
					# -------------------------------------------------------------------------- 
					if(length(scalars)){
						scalar_means <- sapply(detrended_scalars,function(x)mean(x[["fitted"]]+x[["residuals"]]))
						scalar_means_trend <- sapply(detrended_scalars,function(x)mean(x[["fitted"]]))
						scalar_sd <- sapply(detrended_scalars,function(x)sd(x[["residuals"]]))
						names(scalar_means) <- paste0("mean(",names(scalar_means),")") 
						names(scalar_means_trend) <- paste0("mean(trend.",names(scalar_means_trend),")") 
						names(scalar_sd) <- paste0("sd(",names(scalar_sd),")") 
					} else {
						scalar_sd <- scalar_means_trend <- scalar_means <- NULL
					}

					# sub-int: calculate covariances of sub-intervals
					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
					sub_covs <- sapply(sub_wind,function(x)cov(as.data.frame(x[1:4]))[cbind(c("uprot","vprot","wprot","Tdet","uprot","uprot","uprot","vprot","vprot","wprot"),c("uprot","vprot","wprot","Tdet","vprot","wprot","Tdet","wprot","Tdet","Tdet"))])
					sub_cov_means <- apply(sub_covs,1,mean)
					sub_cov_sd <- apply(sub_covs,1,sd)
					names(sub_cov_means) <- paste0("mean(sub.int.",c("<u'u'>","<v'v'>","<w'w'>","<T'T'>","<u'v'>","<u'w'>","<u'T'>","<v'w'>","<v'T'>","<w'T'>"),")")
					names(sub_cov_sd) <- paste0("sd(sub.int.",c("<u'u'>","<v'v'>","<w'w'>","<T'T'>","<u'v'>","<u'w'>","<u'T'>","<v'w'>","<v'T'>","<w'T'>"),")")

					### some output and namings
					fix_lag_out <- fix_lag
					names(fix_lag_out) <- paste("fix.lag",names(fix_lag_out))
					dyn_lag_out <- dyn_lag_max["tau",]
					names(dyn_lag_out) <- paste("dyn.lag",names(dyn_lag_out))
					flux_fix_lag <- sapply(Ogive_fix,"[",1)
					names(flux_fix_lag) <- paste("flux.fix.lag",names(flux_fix_lag))
					flux_dyn_lag <- sapply(Ogive_dyn,"[",1)
					names(flux_dyn_lag) <- paste("flux.dyn.lag",names(flux_dyn_lag))
					fix_damping_pbreg <- sapply(Damping_fix,"[[",1)
					if(!is.null(Damping_fix))names(fix_damping_pbreg) <- paste("damping.pbreg.fix.lag",names(Damping_fix))
					dyn_damping_pbreg <- sapply(Damping_dyn,"[[",1)
					if(!is.null(Damping_dyn))names(dyn_damping_pbreg) <- paste("damping.pbreg.dyn.lag",names(Damping_dyn))
					fix_damping_deming <- sapply(Damping_fix,"[[",2)
					if(!is.null(Damping_fix))names(fix_damping_deming) <- paste("damping.deming.fix.lag",names(Damping_fix))
					dyn_damping_deming <- sapply(Damping_dyn,"[[",2)
					if(!is.null(Damping_dyn))names(dyn_damping_deming) <- paste("damping.deming.dyn.lag",names(Damping_dyn))

					
					# write results:
					# -------------------------------------------------------------------------- 
                    # browser()
					if(checkNames){
						# browser()
						res_temp <- cbind(data.frame(
								Data_Start = Int_Start
								,Data_End = Int_End
								,tzone = tz_sonic
								,"Int_Length (min)" = round(Int_Time/60,4)
								,N = Int_length
								,SubInts = n_stationarity_subint
								,check.names=FALSE,stringsAsFactors=FALSE
								)
							,as.list(c(
								wind_stats
								,scalar_means
								,scalar_means_trend
								,scalar_sd
								,sub_cov_means
								,sub_cov_sd
								,fix_lag_out
								,dyn_lag_out
								,flux_fix_lag
								,flux_dyn_lag
								,fix_damping_pbreg
								,dyn_damping_pbreg
								,fix_damping_deming
								,dyn_damping_deming
								))
							)
						names(results)[-(1:2)] <- names(res_temp)
						checkNames <- FALSE
						results[intval,-(1:2)] <- res_temp
					} else {
						results[intval,-(1:2)] <- cbind(data.frame(
								Data_Start = Int_Start
								,Data_End = Int_End
								,tzone = tz_sonic
								,"Int_Length (min)" = round(Int_Time/60,4)
								,N = Int_length
								,SubInts = n_stationarity_subint
								,check.names=FALSE,stringsAsFactors=FALSE
								)
							,as.list(c(
								wind_stats
								,scalar_means
								,scalar_means_trend
								,scalar_sd
								,sub_cov_means
								,sub_cov_sd
								,fix_lag_out
								,dyn_lag_out
								,flux_fix_lag
								,flux_dyn_lag
								,fix_damping_pbreg
								,dyn_damping_pbreg
								,fix_damping_deming
								,dyn_damping_deming
								))
							)				
					}


					if(create_graphs){
						# plotting:
						# -------------------------------------------------------------------------- 
						# -------------------------------------------------------------------------- 
						cat("~~~\nplotting timeseries and fluxes\n")
						# time series:
						# -------------------------------------------------------------------------- 
			    		# plot and save (rotated) data time series with raw-data trends...
						# ------------------------------------------------------------------------
						time2 <- format(Int_End,format="%H%M%S")
						plotname <- paste("timeseries", filename, time2, sep="-") 
						jpeg(file=paste0(path_folder, '/', plotname,".jpg"),width=600, height=(sum(variables_timeseries))*100, quality=60)
							ts_plot <- plot.tseries(Data[index,],wind,detrended_scalars,names(variables_timeseries)[variables_timeseries],variables_cols,variables_units)
							print(ts_plot)
						dev.off()

			    		# plot and save flux evaluation...
						# ------------------------------------------------------------------------
						for(i in covariances){
							# i <- "w'TDL CH4'"
							plotname <- paste("plots",filename,time2,covariances_plotnames[i],sep="-")
							jpeg(file=paste0(path_folder, '/', plotname,".jpg"),width=1350, height=900, quality=60)
								par(mfrow=c(2,3))
								# ----------------------- Covariance -------------------------------------
								plot_covfunc(Covars[[i]],Int_Time,dyn_lag_max[,i],fix_lag[i],ylab=i, xlim = c(-50,50), cx=1.5, cxmt=1.25, cl=covariances_cols[i])
								# ---------------------- Co-Spec/Ogive fix lag -----------------------------------
								plot_cospec_ogive(Ogive_fix[[i]],Cospec_fix[[i]],freq,ylab=paste0("ogive (fix lag) of ",i),cx=1.5,col=covariances_cols[i])
								# ---------------------- Co-Spec/Ogive dyn lag -----------------------------------
								plot_cospec_ogive(Ogive_dyn[[i]],Cospec_dyn[[i]],freq,ylab=paste0("ogive (dyn lag) of ",i),cx=1.5,col=covariances_cols[i])
								# ---------------------- empirical damping -----------------------------------
								if(scalar_covariances[i]){
									plot_damping(Damping_fix[[i]],freq,ylab=paste0("ogive (fix lag) of ",i),cx=1.5,col=covariances_cols[i])
									plot_damping(Damping_dyn[[i]],freq,ylab=paste0("ogive (dyn lag) of ",i),cx=1.5,col=covariances_cols[i])
								}
								title(paste0(i," flux ",format(Int_Start,format="(%H:%M:%S")," - ",format(Int_End,format="%H:%M:%S)"),if(scalar_covariances[i])paste0(" - damping reference flux: ",reference_covariances[i])),outer=TRUE,line=-1)
							dev.off()	
						}
					} # end plotting
				} else {
					cat("less than 10 data points in raw data - skipping interval.\n")
				}
			} # end for loop intval	

		} # end inner doCalc
		# setwd(former_wd)	
	} # end outer doCalc

    if (as_ibts) {
        results <- as.ibts(results)
    }


	# #################################### END VERSION HISTORY #################################### #
	cat("************************************************************\n") 
	cat("operation finished @", format(Sys.time(), "%d.%m.%Y %H:%M:%S"),"time elapsed: ", difftime(Sys.time(), script.start, unit="mins"),"minutes\n")
	cat("************************************************************\n")  
	if(ogives_out){
		return(list(
                results = results, 
                covars = Covars_Out,
                cospec_fix = Cospec_fix_Out, 
                cospec_dyn = Cospec_dyn_Out,
                ogv_fix = Ogive_fix_Out, 
                ogv_dyn = Ogive_dyn_Out
                ))
	} else {
		return(results)
	}
}


check_damping <- function(i,flux=c("CO2","H2O","CH4")[1],Data=Res1,ref="w'T'"){
	Ogives <- Data[[2]][[i]]
	# Ogives <- Res_2hrs[[2]][[1]]
	# Ogives <- Res_30mins[[2]][[2]]
	freq <- Ogives$freq


	switch(flux
		,CO2 = {
			ylab_wC <- "w'CO2 Analog'"
			f.lim <- c(2,15)
		}
		,H2O = {
			ylab_wC <- "w'H2O Analog'"
			f.lim <- c(2,15)
		}
		,CH4 = {
			ylab_wC <- "w'TDL CH4'"
			f.lim <- c(2,15)
		})

	wC <- Ogives[[ylab_wC]]
	wT <- Ogives[[ref]]
	m <- 3/4
	ylab_wT <- ref


	#### wT:

	# wts <- 1/(1 + abs(1/5 - freq)^2)
	# wts <- 1/(1 + abs(1/120 - freq)^2)
	# wts <- 1/(1 + abs(1/600 - freq)^2)
	# wts0 <- 1/(1 + abs(1/600 - freq))
	wts1 <- 1/(1 + abs(1/600 - freq))/sqrt(freq)
	# wts2 <- 1/(1 + abs(1/60 - freq))/sqrt(freq)
	# wts <- wts1 + 1/(1 + abs(10 - freq))
	# wts <- wts1 + 3/(1 + abs(10 - freq))
	# wts <- wts1 + 3/(1 + abs(10 - freq)^(2))
	wts <- wts1 + 3/(1 + abs(5 - freq)^(2))
	# wts <- wts*0+1
	ind <- seq_along(freq)
	opt_wT <- optim(c(A0=wT[1],fx=0.001,m=3/5,mu=2,offset=0),fn=opt_Ogv,f=freq[ind],Ogv=wT[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
	if(opt_wT$convergence != 0) opt_wT <- NA
	fmax <- c(1/50,2*3600)
	f2 <- seq(1/fmax[2],1/fmax[1],2/fmax[2])
	wT_mod <- Ogive_Model_pred(f2,opt_wT$par,m)
	wT_fit <- Ogive_Model(freq,opt_wT$par,m)



	#### wC:

	# wts <- 1/(1 + abs(1/600 - freq))
	# wts <- wts*0+1
	ind <- seq_along(freq)
	# opt_wC <- optim(c(A0=wC[1],fx=0.001,m=3/5,mu=2,offset=0),fn=opt_Ogv,f=freq[ind],Ogv=wC[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
	opt_wC <- optim(c(A0=wC[1],fx=opt_wT$par[["fx"]],m=3/5,mu=opt_wT$par[["mu"]],offset=opt_wT$par[["offset"]]),fn=opt_Ogv,f=freq[ind],Ogv=wC[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
	if(opt_wC$convergence != 0) opt_wC <- NA
	wC_mod <- Ogive_Model_pred(f2,opt_wC$par,m)
	wC_fit <- Ogive_Model(freq,opt_wC$par,m)
	# add_og <- sum(Cospec_Model(rev(add_f),opt_wC$par,m)*df)

	x11(width=10,height=5)
	par(mfrow=c(1,2))
	plot(freq,wT,log="x",type="l",panel.first={grid();abline(h=0)},xlim=1/fmax,ylim=range(wT_mod,wT,na.rm=TRUE),ylab=ylab_wT)
	lines(freq,wT_fit,col="red")
	lines(f2,wT_mod,col="blue")
	plot(freq,wC,log="x",type="l",panel.first={grid();abline(h=0)},xlim=1/fmax,ylim=range(wC_mod,wC,na.rm=TRUE),ylab=ylab_wC)
	lines(freq,wC_fit,col="red")
	lines(f2,wC_mod,col="blue")


	#### damping wC:

	dmp_1 <- damp_hac5(wC,freq,f.lim,wT)
	dmp_2 <- damp_hac5(wC,freq,f.lim,wT_fit)
	dmp_3 <- damp_hac5(wC_fit,freq,f.lim,wT_fit)
	dmp_4 <- damp_hac5(wC_mod,f2,f.lim,wT_mod)

	ylim <- range(dmp_1$ogive
		,dmp_1$ogive_ref_pbreg,dmp_1$ogive_ref_deming
		,dmp_2$ogive_ref_pbreg,dmp_2$ogive_ref_deming
		,dmp_3$ogive_ref_pbreg,dmp_3$ogive_ref_deming
		,dmp_4$ogive_ref_pbreg,dmp_4$ogive_ref_deming
		)

	x11(width=10,height=10)
	par(mfrow=c(2,2))
	plot_damping(dmp_1,freq,cx.leg=0.7,ylim=ylim)
	plot_damping(dmp_2,freq,cx.leg=0.7,ylim=ylim)
	plot_damping(dmp_3,freq,cx.leg=0.7,col="indianred",ylim=ylim)
	plot_damping(dmp_4,f2,cx.leg=0.7,col="indianred",ylim=ylim)
	invisible(NULL)
}


model_damping <- function(i,flux=c("CO2","H2O","CH4","uw","wT")[1],Data=NULL,ref="w'T'"
	,wtfu=function(freq){
		wts1 <- 1/(1 + abs(1/600 - freq))/sqrt(freq)
		wts1 + 3/(1 + abs(5 - freq)^(2))
	}
	,dyn=FALSE,wtfu_damped=wtfu,f.lim=c(2,15)
	){

	Ogives <- Data[[as.numeric(dyn)+2]][[i]]
	freq <- Ogives$freq

	m <- 3/4
	switch(flux
		,CO2 = {
			ylab_wC <- "w'CO2 Analog'"
			# f.lim <- c(2,15)
		}
		,H2O = {
			ylab_wC <- "w'H2O Analog'"
			# f.lim <- c(2,15)
		}
		,CH4 = {
			ylab_wC <- "w'TDL CH4'"
			# f.lim <- c(2,15)
		}
		,uw = {
			ylab_wC <- "w'u'"
			m <- 3/5
			ref <- NULL
		}
		,wT = {
			ylab_wC <- "w'T'"
			ref <- NULL
		}
		)

	wC <- Ogives[[ylab_wC]]
	main <- paste0(ylab_wC,format(Data$results[i,1],format=": %d.%m.%Y %H:%M - "),format(Data$results[i,2],format="%H:%M")," (i = ",i,")")

	if(is.null(ref)){
		#### wC:

		# # wts <- 1/(1 + abs(1/5 - freq)^2)
		# # wts <- 1/(1 + abs(1/120 - freq)^2)
		# # wts <- 1/(1 + abs(1/600 - freq)^2)
		# # wts0 <- 1/(1 + abs(1/600 - freq))
		# wts1 <- 1/(1 + abs(1/600 - freq))/sqrt(freq)
		# # wts2 <- 1/(1 + abs(1/60 - freq))/sqrt(freq)
		# # wts <- wts1 + 1/(1 + abs(10 - freq))
		# # wts <- wts1 + 3/(1 + abs(10 - freq))
		# # wts <- wts1 + 3/(1 + abs(10 - freq)^(2))
		# wts <- wts1 + 3/(1 + abs(5 - freq)^(2))
		# wts <- wts*0+1
		wts <- wtfu(freq)
		ind <- seq_along(freq)
		# opt_wC <- optim(c(A0=wC[1],fx=0.001,m=3/5,mu=2,offset=0),fn=opt_Ogv,f=freq[ind],Ogv=wC[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
		opt_wC <- optim(c(A0=wC[1],fx=0.001,mu=2,offset=0),fn=opt_Ogv,f=freq[ind],Ogv=wC[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
		if(opt_wC$convergence != 0) opt_wC <- NA
		wC_mod <- Ogive_Model_pred(f2,opt_wC$par,m)
		wC_fit <- Ogive_Model(freq,opt_wC$par,m)
		# add_og <- sum(Cospec_Model(rev(add_f),opt_wC$par,m)*df)
		x11()
		plot(freq,wC,log="x",type="l",panel.first={grid();abline(h=0)},xlim=1/fmax,ylim=range(wC_mod,wC,na.rm=TRUE),ylab=ylab_wC,main=main)
		lines(freq,wC_fit,col="red")
		lines(f2,wC_mod,col="blue")
		# out <- list(freq=f2,model=wC_mod,opt_flux=opt_wC)
		out <- list(flux=wC[1],model=wC_mod[1],opt_flux=opt_wC)
	} else {
		#### wT:
		wT <- Ogives[[ref]]

		ylab_wT <- ref

		# wts <- 1/(1 + abs(1/5 - freq)^2)
		# # wts <- 1/(1 + abs(1/120 - freq)^2)
		# # wts <- 1/(1 + abs(1/600 - freq)^2)
		# # wts0 <- 1/(1 + abs(1/600 - freq))
		# wts1 <- 1/(1 + abs(1/600 - freq))/sqrt(freq)
		# # wts2 <- 1/(1 + abs(1/60 - freq))/sqrt(freq)
		# # wts <- wts1 + 1/(1 + abs(10 - freq))
		# # wts <- wts1 + 3/(1 + abs(10 - freq))
		# # wts <- wts1 + 3/(1 + abs(10 - freq)^(2))
		# wts <- wts1 + 3/(1 + abs(5 - freq)^(2))
		# # wts <- wts*0+1
		wts <- wtfu(freq)
		ind <- seq_along(freq)
		opt_wT <- optim(c(A0=wT[1],fx=0.001,mu=2,offset=0),fn=opt_Ogv,f=freq[ind],Ogv=wT[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
		if(opt_wT$convergence != 0) opt_wT <- NA
		fmax <- c(1/50,2*3600)
		f2 <- seq(1/fmax[2],1/fmax[1],2/fmax[2])
		wT_mod <- Ogive_Model_pred(f2,opt_wT$par,m)
		wT_fit <- Ogive_Model(freq,opt_wT$par,m)

		#### wC:
		wts <- wtfu_damped(freq)
		# wts <- 1/(1 + abs(1/600 - freq))
		# wts <- wts*0+1
		ind <- seq_along(freq)
		# opt_wC <- optim(c(A0=wC[1],fx=0.001,m=3/5,mu=2,offset=0),fn=opt_Ogv,f=freq[ind],Ogv=wC[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
		opt_wC <- optim(c(A0=wC[1],fx=opt_wT$par[["fx"]],mu=opt_wT$par[["mu"]],offset=opt_wT$par[["offset"]]),fn=opt_Ogv,f=freq[ind],Ogv=wC[ind],m=m,control=list(maxit=5E4),wts=wts[ind])
		if(opt_wC$convergence != 0) opt_wC <- NA
		wC_mod <- Ogive_Model_pred(f2,opt_wC$par,m)
		wC_fit <- Ogive_Model(freq,opt_wC$par,m)

		x11(width=10,height=5)
		par(mfrow=c(1,2))
		plot(freq,wT,log="x",type="l",panel.first={grid();abline(h=0)},xlim=1/fmax,ylim=range(wT_mod,wT,na.rm=TRUE),ylab=ylab_wT)
		lines(freq,wT_fit,col="red")
		lines(f2,wT_mod,col="blue")
		plot(freq,wC,log="x",type="l",panel.first={grid();abline(h=0)},xlim=1/fmax,ylim=range(wC_mod,wC,na.rm=TRUE),ylab=ylab_wC)
		lines(freq,wC_fit,col="red")
		lines(f2,wC_mod,col="blue")
		title(main,outer=TRUE,line=-1)

		#### damping wC:

		dmp_1 <- damp_hac5(wC,freq,f.lim,wT)
		# dmp_2 <- damp_hac5(wC,freq,f.lim,wT_fit)
		# dmp_3 <- damp_hac5(wC_fit,freq,f.lim,wT_fit)
		dmp_4 <- damp_hac5(wC_mod,f2,f.lim,wT_mod)

		ylim1 <- range(dmp_1$ogive,dmp_1$ogive_ref_pbreg,dmp_1$ogive_ref_deming)
			# ,dmp_2$ogive_ref_pbreg,dmp_2$ogive_ref_deming
			# ,dmp_3$ogive_ref_pbreg,dmp_3$ogive_ref_deming
		ylim2 <- range(dmp_4$ogive,dmp_4$ogive_ref_pbreg,dmp_4$ogive_ref_deming)
			

		x11(width=10,height=5)
		par(mfrow=c(1,2))
		plot_damping(dmp_1,freq,cx.leg=0.7,ylim=ylim1)
		# plot_damping(dmp_2,freq,cx.leg=0.7,ylim=ylim)
		# plot_damping(dmp_3,freq,cx.leg=0.7,col="indianred",ylim=ylim)
		plot_damping(dmp_4,f2,cx.leg=0.7,col="indianred",ylim=ylim2)
		title(main,outer=TRUE,line=-1)
		# out <- list(freq=f2,model=wC_mod,ref_model=wT_mod,damping=dmp_4,opt_flux=opt_wC,opt_ref=opt_wT)
		out <- list(flux=wC[1],model=wC_mod[1],damping_pbreg=dmp_4$dampf_pbreg,damping_deming=dmp_4$dampf_deming,opt_flux=opt_wC,opt_ref=opt_wT)
	}
	invisible(out)
}

Cospec_Model <- function(f,x,m=NULL){
	if(is.null(m)) m <- x["m"]
	mu2 <- 2*x["mu"]
	x["A0"]/x["fx"] / (1 + m*(f/x["fx"])^mu2)^((m+1)/m/mu2)
}

# Ogive_Model <- function(f,x,m=NULL,f1=100){
# 	df <- f[2]-f[1]
# 	add_f <- seq(f[length(f)],f1,df)
# 	add_og <- sum(Cospec_Model(rev(add_f),x,m)*df)
# 	x["offset"] + rev(cumsum(Cospec_Model(rev(f),x,m)*df)) + add_og
# }
Ogive_Model <- function(f,x,m=NULL){
	x["offset"] + rev(cumsum(Cospec_Model(rev(f),x,m)*(f[2]-f[1])))
}

Ogive_Model_pred <- function(f,x,m=NULL){
	rev(cumsum(Cospec_Model(rev(f),x,m)*(f[2]-f[1])))
}




opt_Ogv <- function(x,f,Ogv,m=NULL,wts=1){
	out <- try(sum(abs(Ogive_Model(f,x,m) - Ogv)*wts)/sum(wts))
	if(inherits(out,"try-error")){
		return(1E9)
	} else {
		return(out)
	}
}



################################ Funktionen von Marcel ####################################################

################## Funktion zur Berechnung der Pfadlängen ##############################################
Path.Length.SR <- function(Sensor,Reflektor){
	pathLength <- sqrt((Sensor[1]-Reflektor[1])^2+(Sensor[2]-Reflektor[2])^2) 
	return(pathLength)
}
