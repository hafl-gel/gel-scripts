## 13.04.2018 - Christoph - 
##     Funktion 'readWindMaster_ascii' abgeändert. Korrupte Files werden führen nun nicht mehr zu einem Fehler.
##     Funktion 'evalReddy': Argument 'create_graphs' hinzugefügt (FALSE -> kein plotten).
##                           ca. L814: zu kurze Messfiles werden abgefangen und führen nicht mehr zu einem Fehler.
# *************************************************************************************
## 16.01.2018 - Marcel - Funktion zur Berechnung der Pfadlänge GF hinzugefügt
# *************************************************************************************

REddy.version <- "v20190905"
library(tcltk)
library(xlsx)
library(lubridate)
library(Rcpp)
library(reshape2)
library(lattice)
library(caTools)
library(data.table)

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

# find time in intervals (from ibts: findI_st)
cppFunction("
	IntegerVector findInterval2(NumericVector x,NumericVector y1,NumericVector y2){
		// case: st_closed = TRUE
		int lenx = x.size();
		IntegerVector Out = rep(0,lenx);
		int leny = y1.size() - 1;
		int run = 0;

		for(int i = 0; (i < lenx) & (run <= leny); i++){
			if(x[i] >= y2[leny]){
				break;
			}
			while((y2[run] <= x[i]) & (run <= leny)){
				run += 1;
			}
			if(y1[run] <= x[i]){
				Out[i] = run + 1;
			}
		}
		return Out;
	}
	")

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
	if(any(hflgs))cat("Hard Flag: NA flagged values in columns ",paste(names(hflgs)[hflgs],sep=", "),"\n")
	
	if(pmatch(hflg.met, "replace", nomatch = 0) && any(hflgs)){

		cat("Replacing NA values by running mean...\n")
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
		st2 <- c(st1[-1],st1[length(st1)] + d_t)
		for(i in seq_along(isna)){
			x1 <- st1[isna[[i]]] - wind/2 + d_t/2
			x2 <- st1[isna[[i]]] + wind/2 - d_t/2

			ind <- cutIntervals(st1,st2,x1,x2)

			dat[hflgs][[i]][isna[[i]]] <- sapply(ind,function(x)mean(dat[hflgs][[i]][x[,1]],na.rm=TRUE))
		}	
		cat("number of replaced values\n*~~~~*\n", names(dat[hflgs]),"\n", lengths(isna),"\n*~~~~*\n")
	}
							
	dat
	
}

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


rotate_detrend <- function(u,v,w,T,phi=NULL,method=c(u="blockAVG",v="blockAVG",w="blockAVG",T="linear"),c.system="Windmaster",Hz_ts=10){
	thetam <- atan2(.Internal(mean(v)),.Internal(mean(u)))
	u1 <- u*cos(thetam) + v*sin(thetam)
	if(is.null(phi)){
		phi <- atan2(.Internal(mean(w)),.Internal(mean(u1)))
	}
	ud <- trend(u1*cos(phi) + w*sin(phi),method["u"],Hz_ts)
	vd <- trend(-u*sin(thetam) + v*cos(thetam),method["v"],Hz_ts)
	wd <- trend(-u1*sin(phi) + w*cos(phi),method["w"],Hz_ts)		
	Td <- trend(T,method["T"],Hz_ts)
	# Yamartino 1984:
	n <- length(u)
	theta_i <- atan2(v,u)
	delta_i <- abs(theta_i - thetam)
	delta_i <- ifelse((thetam - pi) < theta_i & theta_i < thetam,-delta_i,delta_i)
	sd_wd <- sqrt(sum(delta_i^2)/n - (sum(delta_i)/n)^2) / pi * 180

	out <- list(
		uprot=ud$residuals
		,vprot=vd$residuals
		,wprot=wd$residuals
		,Tdet=Td$residuals
		# ,wd=if(c.system=="Windmaster") (180 - thetam / pi * 180)%%360 else if(c.system=="ART.EC1") ((180 / pi) * -thetam + 150 + 147)%%360
		,wd=switch(c.system
			,"Windmaster" = (180 - thetam / pi * 180) %% 360 
			,"ART.EC1" = ((180 / pi) * -thetam + 150 + 147) %% 360
			# new (coordinate) system
			,"Licor" = (270 - thetam/pi*180) %% 360
			# # old (coordinate) system
			# ,"Licor" = (180 + thetam / pi * 180) %% 360 
			)
		,sd_wd=sd_wd
		,phi=phi
		,umrot=ud$fitted
		,vmrot=vd$fitted
		,wmrot=wd$fitted
		,Tmdet=Td$fitted
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
	freq_out <- exp(log_cuts[-1] - diff(log_cuts))
	list(cospec=tapply(cospec,ind_cuts,sum),freq=freq_out[unique(ind_cuts)],d_freq=diff(exp(log_cuts))[unique(ind_cuts)])
}

damp_hac5 <- function(ogive, freq, freq.limits, ogive_ref){
	# browser()
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
	# mod <- lm(undamped_ogv~undamped_ogvref,weights=1/freq[undamped_ind])
	# pred_ogv_lm <- predict(mod,newdata=list(undamped_ogvref=ogive_ref))
	mod <- deming::deming(undamped_ogv~undamped_ogvref,weights=1/freq[undamped_ind])
	cfs <- coef(mod)
	pred_ogv_deming <- cfs[2]*ogive_ref + cfs[1]

	# robust linear regression + prediction:
	# mod2 <- MASS::rlm(undamped_ogv~undamped_ogvref,weights=1/freq[undamped_ind])
	# pred_ogv_rlm <- predict(mod2,newdata=list(undamped_ogvref=ogive_ref))

	mod2 <- deming::pbreg(undamped_ogv~undamped_ogvref,method=3,eps=min(abs(undamped_ogv))*1E-8)
	cfs2 <- coef(mod2)
	pred_ogv_pbreg <- cfs2[2]*ogive_ref + cfs2[1]
	
	list(dampf_pbreg=ogive[1]/(ogive[1] - cfs2[1]),dampf_deming=ogive[1]/(ogive[1] - coef(mod)[1]),freq.limits=freq.limits,ogive=ogive,ogive_ref_pbreg=pred_ogv_pbreg,ogive_ref_deming=pred_ogv_deming)	
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
	Cov_sonic <- cov(as.data.frame(wind[1:4]))
	Var_sonic <- diag(Cov_sonic)
	names(Var_sonic) <- c("<u'u'>","<v'v'>","<w'w'>","<T'T'>")
	Cov_sonic <- Cov_sonic[cbind(c("uprot","uprot","uprot","vprot","vprot","wprot"),c("vprot","wprot","Tdet","wprot","Tdet","Tdet"))]
	names(Cov_sonic) <- c("<u'v'>","<u'w'>","<u'T'>","<v'w'>","<v'T'>","<w'T'>")
	Ustar <- c(sqrt(-Cov_sonic["<u'w'>"]),use.names = FALSE)
	T_K <- mean(wind$Tmdet + wind$Tdet)
	U <- mean(wind$umrot + wind$uprot)
	L <- c(-Ustar^3 * T_K / (0.4 * 9.80620 * Cov_sonic["<w'T'>"]),use.names = FALSE)
	if(!is.na(z_canopy)){
		d <- 2/3 * z_canopy
		suppressWarnings(z0 <- optimize(function(x,ustar,L,z,d,U)abs(U - calcU(ustar, x, L, z-d)),c(0,z_sonic*1.1),ustar=Ustar,L=L,U=U,z=z_sonic,d=d)$minimum)
		if(z0>=z_sonic)z0 <- NA
	} else {
		d <- NA
		z0 <- NA
	}
	c(Var_sonic,Cov_sonic,Ustar=Ustar,L=L,z_sonic=z_sonic,z_canopy=z_canopy,d=d,z0=z0,WD=wind$wd,"sd(WD)"=wind$sd_wd,phi=wind$phi,U_sonic=U,T_sonic=T_K,U_trend=mean(wind$umrot),T_trend=mean(wind$Tmdet))
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
	dat3 <- data.frame(wind[c("umrot","vmrot","wmrot","Tmdet")])
	names(dat3) <- c("u","v","w","T")
	if(!is.null(scal)){
		dat3 <- cbind(dat3,lapply(scal,"[[","fitted"))
	}
	dat3 <- dat3[,selection]
	### melt and add trends:
	dat4 <- reshape2::melt(dat3,id=NULL,value.name="trend")
	dat2[,"trend"] <- dat4[,"trend"]

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
		strip=FALSE, layout=c(1,ncol(dat)-2), between=list(x=0,y=1), subscripts=TRUE, lwd=rep(1,length(color)), lty=rep(1,length(color)), col=color,
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
	freq_r0 <- cospec_reduced0$freq
	cospec_f <- cospec_reduced0$cospec/cospec_reduced0$d_freq*freq_r0
	# reduced cospec 2:
	cospec_reduced <- reduce_cospec(cospec,freq,nred)
	freq_rm <- cospec_reduced$freq
	cospec_rm <- cospec_reduced$cospec/cospec_reduced$d_freq*freq_rm

	rCo <- range(cospec_rm,na.rm=TRUE)
	if(is.null(xlim))xlim <- rev(range(freq))
	ylim <- c(min(ogive,0),max(0,max(ogive)))
	pxlim <- pretty(log10(xlim),n=ceiling(abs(diff(log10(xlim)))))
	pxlims <- rep(pxlim,each=9) - log10(seq(9, 1))
	pylim <- pretty(ylim)
	prCo <- pretty(rCo) 
	y_cf <- (cospec_f - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	y_crm <- (cospec_rm - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	py2 <- (prCo - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	
	plot(1,xlim=xlim,ylim=ylim, cex.axis=cx, cex.lab=cx,type="n",log="x",xaxt="n",yaxt="n",xlab="log of frequency [Hz]",ylab="",panel.first=abline(h=0,col=col,lty=2))
	abline(h=(0 - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1],lty=2,col="darkgrey")
	axis(1,at=10^pxlims,labels=FALSE,tck=-0.01, cex.axis=cx, cex.lab=cx)
	axis(1,at=10^pxlim,labels=pxlim, cex.axis=cx, cex.lab=cx)
	axis(2,at=pylim,labels=pylim, cex.axis=cx, cex.lab=cx, col=col,col.axis=col,lwd=2,font=2)
	title(ylab=ylab,col.lab=col, cex.lab=cx,font.lab=2)
	axis(4,at=py2,labels=prCo, cex.axis=cx, cex.lab=cx)
	axis(3,at=1/c(0.02,0.05,0.1,1,10,60,120,300,600,1800,3600,7200),labels=c("20ms","50ms","100ms","1s","10s","1min","2min","5min","10min","30min","1hr","2hrs"), cex.axis=cx, cex.lab=cx)
	lines(freq_r0,y_cf,col="lightgrey")
	# # ref cospec:
	# if(!is.null(cospec_ref)){
	# 	cospec_refred <- reduce_cospec(cospec_ref,freq,nred)
	# 	cospec_rf <- cospec_refred$cospec/cospec_refred$d_freq*freq_rm
	# 	y_crf <- (cospec_rf - min(prCo))/diff(range(prCo))*diff(ylim) + ylim[1]
	# 	lines(freq_rm,y_crf,type="b",col="darkgrey",lwd=2)
	# }
	lines(freq_rm,y_crm,type="b",col="black",lwd=2)
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
	pxlims <- rep(pxlim,each=9) - log10(seq(9, 1))
	pylim <- pretty(ylim)

	plot(1,xlim=xlim,ylim=ylim,main=main, cex.axis=cx, cex.lab=cx,type="n",log="x",xaxt="n",xlab="log of frequency [Hz]",ylab=ylab,panel.first=abline(h=0,col=col,lty=2))
	abline(h=0,lty=2,col="darkgrey")
	axis(1,at=10^pxlims,labels=FALSE,tck=-0.01, cex.axis=cx, cex.lab=cx)
	axis(1,at=10^pxlim,labels=pxlim, cex.axis=cx, cex.lab=cx)
	axis(3,at=1/c(0.02,0.05,0.1,1,10,60,120,300,600,1800,3600,7200),labels=c("20ms","50ms","100ms","1s","10s","1min","2min","5min","10min","30min","1hr","2hrs"), cex.axis=cx, cex.lab=cx)
	lines(freq,ogive_ref_deming,col="darkgrey",lwd=2)
	lines(freq,ogive_ref_pbreg,col="lightgrey",lwd=2)
	lines(freq,ogive,col=col,lwd=2)
	abline(h=(1 - 1/dampf_deming)*ogive[1],col="darkgrey",lty=2)
	abline(h=(1 - 1/dampf_pbreg)*ogive[1],col="lightgrey",lty=2)
	abline(v=1/freq.limits,lty=4,col="black")
	pos <- if(diff(abs(ylim))>0) "topleft" else "bottomleft"
	legend(pos,legend=c("damped","scaled reference (pbreg)","scaled reference (deming)",sprintf("damping (pbreg): %1.1f%%",(1 - dampf_pbreg)*100),sprintf("damping (deming): %1.1f%%",(1 - dampf_deming)*100)),col=c(col,"lightgrey","darkgrey",NA,NA),bty="n",lty=1,lwd=2,cex=cx.leg)
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
	){
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
		len_variables <- length(variables)
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
				dir.create(paste0(save_directory,"/",folder),recursive=FALSE)
				on.exit(setwd(former_wd))
				setwd(paste0(save_directory,"/",folder))
			}

			read_rawdata <- get(rawdata_function,mode="function")

			# reading data file:                                      
			# ------------------------------------------------------------------------------
			cat("Reading raw data...\n")
			# browser() 
			Data <- read_rawdata(dfname)
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
				cat("Calculation will include all intervals between",format(st_interval[1],format="%Y-%m-%d %H:%M"),"and",format(et_interval[length(et_interval)],format="%Y-%m-%d %H:%M",usetz=TRUE),"on a",avg_period,"min basis\n")
			} else {
				st_interval <- Start
				et_interval <- End
				cat("Calculation will include the following intervals:\n",paste(format(st_interval,format="   %Y-%m-%d %H:%M"),"to",format(et_interval,format="%Y-%m-%d %H:%M",usetz=TRUE)))
			}

			# prepare results:                                      
			# ------------------------------------------------------------------------------ 
			POSIxct_dummy <- rep(now(),length(st_interval))
			POSIxct_dummy[] <- NA
			results <- cbind(
				data.frame(
					Interval_Start=st_interval
					,Interval_End=et_interval
					,Data_Start = POSIxct_dummy
					,Data_End = POSIxct_dummy
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
			ind_interval <- findInterval2(as.numeric(Data[,1]),as.numeric(st_interval),as.numeric(et_interval))
			unique_ind <- unique(ind_interval) %w/o% 0

			if(ogives_out){
				Ogive_fix_Out <- Ogive_dyn_Out <- vector("list",length(et_interval))
				names(Ogive_fix_Out) <- names(Ogive_dyn_Out) <- format(et_interval,format="%H%M")
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
					Int_Time <- as.numeric(Int_End - Int_Start,units="secs")
					d_t <- diff(as.numeric(Time))
					Hz <- 1/summary(d_t)
					freq <- seq(floor(Int_length/2))*2/Int_Time
					# for later: include NA where d_t>2*mean, remove entries where d_t < 0.5*mean

					# fix and dyn lags:                                      
					# ------------------------------------------------------------------------------ 
					fix_raw <- unlist(df.lag.fix.sec[format(Data[1,1],"%y%m%d"),])
					fix_lag <- round(fix_raw*Hz["Mean"])
					dyn_lag <- rbind(
						lower=round((fix_raw - lag_dev)*Hz["Mean"])
						,upper=round((fix_raw + lag_dev)*Hz["Mean"])
						)
					# # check dyn lag <= 0 for scalars:
					# dyn_lag[2,scalar_covariances] <- pmin(dyn_lag[2,scalar_covariances],0)

					# raw data quality control I, i.e. hard flags = physical range
					# --------------------------------------------------------------------------
					cat("~~~\nchecking hard flags...\n")
					if(any(hard_flag)){
						Data_int[,variables[hard_flag]] <- H.flags(Data_int[,variables[hard_flag],drop=FALSE],Time,Hz["Mean"],lim_range[,variables[hard_flag],drop=FALSE],hf_window,hf_method)
					} 

					# check for remaining NA values
					if(anyNA(Data_int[,c("u","v","w","T", scalars)])){
						cat("NA values in data. Skipping interval...")
						next
					}

					# calculate wind direction, rotate u, v, w, possibly detrend T (+ u,v,w)
					# -------------------------------------------------------------------------- 
					cat("~~~\nrotating and deterending sonic data...\n")
					wind <- rotate_detrend(Data_int[,"u"],Data_int[,"v"],Data_int[,"w"],Data_int[,"T"],phi=if(rotation_method %in% "two_axis") NULL else 0,method=detrending_method[c("u","v","w","T")],c.system=sonic_type,Hz_ts=Hz["Mean"])
					Data_int[,c("u","v","w","T")] <- wind[c("uprot","vprot","wprot","Tdet")]

					# detrend scalars
					# -------------------------------------------------------------------------- 
					if(length(scalars)){
						cat("~~~\ndetrending scalars...\n")
						detrended_scalars <- mapply(trend,y=Data_int[,scalars, drop = FALSE],method=detrending_method[scalars],MoreArgs=list(Hz_ts=Hz["Mean"]),SIMPLIFY=FALSE)
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
						Ogive_fix_Out[[intval]] <- c(list(freq=freq),Ogive_fix)
						Ogive_dyn_Out[[intval]] <- c(list(freq=freq),Ogive_dyn)
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
					sub_wind <- lapply(sub_Data,function(x)rotate_detrend(x[,"u"],x[,"v"],x[,"w"],x[,"T"],phi=if(rotation_method %in% "two_axis") NULL else 0,method=subint_detrending_method[c("u","v","w","T")],c.system=sonic_type,Hz_ts=Hz["Mean"]))
					for(sub in seq_along(sub_wind)){
						sub_Data[[sub]][,c("u","v","w","T")] <- sub_wind[[sub]][c("uprot","vprot","wprot","Tdet")]
					}

					# sub-int: detrend scalars
					# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
					if(length(scalars)){
						sub_detrended_scalars <- lapply(sub_Data,function(x)mapply(trend,y=x[,scalars, drop = FALSE],method=subint_detrending_method[scalars],MoreArgs=list(Hz_ts=Hz["Mean"]),SIMPLIFY=FALSE))
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
						jpeg(file=paste0(plotname,".jpg"),width=600, height=(sum(variables_timeseries))*100, quality=60)
							ts_plot <- plot.tseries(Data[index,],wind,detrended_scalars,names(variables_timeseries)[variables_timeseries],variables_cols,variables_units)
							print(ts_plot)
						dev.off()

			    		# plot and save flux evaluation...
						# ------------------------------------------------------------------------
						for(i in covariances){
							# i <- "w'TDL CH4'"
							plotname <- paste("plots",filename,time2,covariances_plotnames[i],sep="-")
							jpeg(file=paste0(plotname,".jpg"),width=1350, height=900, quality=60)
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

			if(!isFALSE(save_directory)){
				# write results
				# ------------------------------------------------------------------------
				write.table(results,file=paste0("_Reddy-",filename,".csv"),sep=";",row.names=FALSE)
			}
		} # end inner doCalc
		setwd(former_wd)	
	} # end outer doCalc


	# #################################### END VERSION HISTORY #################################### #
	cat("************************************************************\n") 
	cat("operation finished @", format(Sys.time(), "%d.%m.%Y %H:%M:%S"),"time elapsed: ", difftime(Sys.time(), script.start, unit="mins"),"minutes\n")
	cat("************************************************************\n")  
	if(ogives_out){
		return(list(results=results,ogv_fix=Ogive_fix_Out,ogv_dyn=Ogive_dyn_Out))
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
