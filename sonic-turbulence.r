

library(data.table)
library(ibts)


#### To Do:
# - neue readWindMaster Routine
# - Plotting & anderes Gheu aufräumen

readWindMaster_ascii <- function(FilePath, tz = "Etc/GMT-1"){
	### get Date
	bn <- basename(FilePath)
	if(!grepl("^data_", bn)){
		# run old script
		return(readWindMaster_old_ascii(FilePath, tz))
	}
    if (grepl('[.]gz$', bn)) {
        require(R.utils)
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
        if (grep('File is empty', out)) {
            cat('empty\n')
        } else {
            cat('error reading file\n')
            cat(out)
        }
		return(NULL)
	} else {
        out <- out[grepl('^[\x01-\x1A]', V9)]
        # check which columns to convert columns if necessary
        vnums <- paste0('V', c(3, 4, 5, 7))
        is.char <- out[, sapply(.SD, is.character), .SDcols = vnums]
        # convert to numeric
        if (any(is.char)) {
            out[, vnums[is.char] := {
                lapply(.SD, as.numeric)
            }, .SDcols = vnums[is.char]]
        }
        # remove NA lines that come from conversion
        out <- na.omit(out)
        # be verbose and print sonic names:
        sonic <- out[, sub('[\x01-\x1A]', '', V2[1])]
        cat(paste0("data recorded by sonic-", tolower(sonic), "\n"))
        sonic_file <- sub("data_(.*)_[0-9]{8}_[0-9]{6}([.]gz)?$", "\\1", bn)
        if(sonic != toupper(sub("sonic-", "", sonic_file))){
            warning(paste0("Sonic name '", sonic, "', and hostname '", sonic_file, "' don't match!"), call. = FALSE)
        }
        # check units
        if(out[, V6[1]] != "M"){
            stop("Units of recorded data not compatible with evaluation script! Column 6 should contain 'M' for m/s!")
        }
        # only call dt once
        out[, c(
            # remove them
            'V1', 'V2', 'V6', 'V8', 'V9',
            # add them
            'sonic', 'Time', 'Hz',
            # replace °C by K
            'V7'
            ) := {
            # set times correctly
            st.dec <- fast_strptime(paste(Date, V1), lt = FALSE, format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")
            # get start time
            start_time <- as.POSIXct(st.dec[1])
            dt <- as.numeric(st.dec - start_time, units = "secs")
            # correct new day
            ind <- trunc(dt) <= 0L
            # exclude first cuple of 20Hz data since trunc 0
            ind[1:min(.N, 30)] <- FALSE
            dt[ind] <- dt[ind] + 24 * 3600
            # return list
            list(
                # remove them
                NULL, NULL, NULL, NULL, NULL,
                # sonic
                sonic,
                # Time
                start_time + dt,
                # Hz
                round(median(tabulate(trunc(dt))), -1),
                # T
                V7 + 273.15
                )
        }]
        ### set Output names and order
        setnames(out, c("u", "v", "w", "T", "sonic", "Time","Hz"))
        setcolorder(out,c("Time","Hz","u","v","w","T", "sonic"))
        # check if old sonic
        if (out[, sonic[1] %in% c('C', 'D')]) setattr(out, 'OldDevice', TRUE)
        # return
        out
    }
}

# readWindMaster_ascii("~/repos/3_Scripts/5_shellSonic/test_data/data_sonicb_20210120_170910")


readWindMaster_old_ascii <- function(FilePath, tz = "Etc/GMT-1"){
	### get Date
	bn <- basename(FilePath)
	Date <- gsub("^..._([0-9]{6})_.*","\\1",bn)
	### read File
	# browser()
	suppressWarnings(out <- fread(cmd=paste0("grep -v -e ',,' -e [A-Za-z] '",path.expand(FilePath),"'"),fill = TRUE,blank.lines.skip = TRUE))
	if(nrow(out) == 0){
		cat("File empty:",path.expand(FilePath),"\n")
		return(NULL)
	}
	# remove first (empty) column
    if ('V1' %in% names(out)) out[, V1 := NULL]
    # remove NAs
 	out <- na.omit(out)
	### set times
	out[,st.dec := fast_strptime(paste0(Date,V2),lt = FALSE,format = "%y%m%d%H:%M:%S",tz = "Etc/GMT-1")+V3][,c("V2","V3"):=NULL]
	# get start time
	start_time <- out[, as.POSIXct(trunc(st.dec[1]))]
	out[, dt := as.numeric(st.dec - start_time, units = "secs")]
	# correct new day
	out[(seq_len(.N) > 30) & trunc(dt) == 0, dt := dt + 24 * 3600]
	out[trunc(dt) < 0, dt := dt + 24 * 3600]
	# add st column
	out[, st := start_time + dt]
	# add Hz column
	Hz <- out[, .N, by = trunc(dt)][, round(median(N), -1)]
	out[, Hz := Hz]
	# remove columns
	out[, c("st.dec", "dt") := NULL]
	### set Output names and order
	setnames(out,c("u", "v", "w", "T", "Time","Hz"))
	setcolorder(out,c("Time","Hz","u","v","w","T"))
	### remove 999.99 entries
	out <- out[!(u%in%999.99|v%in%999.99|w%in%999.99|T%in%999.99),]
	### change units from °C to K
	out[,T := T + 273.15]
    setattr(out, "OldDevice", TRUE)
	out
}

# readWindMaster_old_ascii("~/LFE/01_Projekte/04_FerARA/Guerbetal/Daten/Sonic/Sonic1/HS1_190408_0000")
# readWindMaster_ascii("~/LFE/01_Projekte/04_FerARA/Guerbetal/Daten/Sonic/Sonic1/HS1_190408_0000")

readSonicEVS_csv <- function(FilePath, tz = "Etc/GMT-1"){
	### read File
	suppressWarnings(out <- fread(cmd=paste0("grep -v -e ',,' -e [A-Za-z] '", 
		path.expand(FilePath), "'"), fill = TRUE, blank.lines.skip = TRUE))
	### remove rows with NA
 	out <- na.omit(out)
	### set times
	out[,st.dec := fast_strptime(paste(V1, V2), lt = FALSE, format = "%d.%m.%Y %H.%M.%OS",
		tz = tz)][,c("V1","V2"):=NULL]
	# get start time
	start_time <- out[, as.POSIXct(trunc(st.dec[1]))]
	out[, dt := trunc(as.numeric(st.dec - start_time, units = "secs"))]
	# correct new day
	out[(seq_len(.N) > 30) & dt == 0, dt := dt + 24 * 3600]
	out[dt < 0, dt := dt + 24 * 3600 - 1]
	# add st column
	out[, st := start_time + dt]
	# add Hz column
	Hz <- out[, .N, by = trunc(dt)][, round(median(N), -1)]
	out[, Hz := Hz]
	# remove columns
	out[, c("st.dec", "dt") := NULL]
	### set Output names and order
	setnames(out,c("u", "v", "w", "T", "Time","Hz"))
	setcolorder(out,c("Time","Hz","u","v","w","T"))
	### remove 999.99 entries
	out <- out[!(u%in%999.99|v%in%999.99|w%in%999.99|T%in%999.99),]
	### change units from °C to K
	out[,T := T + 273.15]

	out
}


avgStats <- function(FilePath, granularity = "30mins", st_to = NULL, 
	et_to = NULL, as_ibts = FALSE, tz = "Etc/GMT-1", 
	rawdata_function = readSonicEVS_csv, floor_to = granularity, 
	ceiling_to = granularity){
	xx <- rawdata_function(FilePath, tz = tz)
	if(!is.null(xx)){
		# prepare times
		gran <- parse_time_diff(granularity)
		if(is.null(st_to)){
			st_to <- xx[, floor_time(Time[1], floor_to)]
		} else {
			st_to <- parse_date_time3(st_to, tz = tz)
		}
		if(is.null(et_to)){
			et_to <- xx[, ceiling_time(Time[.N], ceiling_to)]
		} else {
			et_to <- parse_date_time3(et_to, tz = tz)
		}
		sttime_to <- seq(st_to, et_to - gran, by = gran)
		ettime_to <- seq(st_to + gran, et_to, by = gran)
		# get time index	
		xx[,timeIndex := findI_st(as.numeric(Time), as.numeric(sttime_to), as.numeric(ettime_to))]
		# exclud time outside time range & return calculation
		xx[timeIndex > 0,{
			avg_u <- mean(u)
			res_u <- u - avg_u
			avg_v <- mean(v)
			res_v <- v - avg_v
			avg_w <- mean(w)
			res_w <- w - avg_w
			var_u <- mean(res_u ^ 2)
			var_v <- mean(res_v ^ 2)
			var_w <- mean(res_w ^ 2)
			avg_U <- sqrt(avg_u ^ 2 + avg_v ^ 2 + avg_w ^ 2)
			U <- sqrt(u ^ 2 + v ^ 2 + w ^ 2)
			eddy_U <- U - mean(U)
			U_scalar <- mean(U)
			list(
				st = sttime_to[timeIndex],
				et = ettime_to[timeIndex],
				n = .N,
				u = avg_u,
				v = avg_v,
				w = avg_w,
				u_sd = sqrt(var_u),
				v_sd = sqrt(var_v),
				w_sd = sqrt(var_w),
				U_scalar = U_scalar,
				Us_horiz = mean(sqrt(u ^ 2 + v ^ 2)),
				U_vector = avg_U,
				Uv_horiz = sqrt(avg_u ^ 2 + avg_v ^ 2),
				I_cfd = sqrt(var_u + var_v + var_w) / 3 / avg_U,
				I_ams = sqrt(mean(eddy_U ^ 2)) / U_scalar,
				T = mean(T),
				T_sd = sd(T)
				)
		}, by = timeIndex][, timeIndex := NULL]
	} else {
		NULL
	}
}


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
	
	if(hflg.met %in% "replace" && any(hflgs)){

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

detrend <- function(u,v,w,T,method=c(u="blockAVG",v="blockAVG",w="blockAVG",T="linear"),Hz_ts=10){
	ud <- trend(u,method["u"],Hz_ts)
	vd <- trend(v,method["v"],Hz_ts)
	wd <- trend(w,method["w"],Hz_ts)		
	Td <- trend(T,method["T"],Hz_ts)

	list(
		up=ud$residuals
		,vp=vd$residuals
		,wp=wd$residuals
		,Tp=Td$residuals
		,um=ud$fitted
		,vm=vd$fitted
		,wm=wd$fitted
		,Tm=Td$fitted
		)	
}

rotate_twoaxis <- function(u,v,w,T,phi=NULL,c.system="Windmaster"){
	thetam <- atan2(.Internal(mean(v)),.Internal(mean(u)))
	# Yamartino 1984:
	n <- length(u)
	theta_i <- atan2(v,u)
    delta_i <- ((theta_i - thetam) + pi) %% (2 * pi) - pi
    sd_wd <- sd(delta_i) / pi * 180
	# delta_i <- abs(theta_i - thetam)
	# delta_i <- ifelse((thetam - pi) < theta_i & theta_i < thetam,-delta_i,delta_i)
	# sd_wd <- sqrt(sum(delta_i^2)/n - (sum(delta_i)/n)^2) / pi * 180
	# rotate
	u1 <- u*cos(thetam) + v*sin(thetam)
	if(is.null(phi)){
		phi <- atan2(.Internal(mean(w)),.Internal(mean(u1)))
	}

	list(
		phi = phi
		,wd = if(tolower(c.system) %in% "windmaster") (180 - thetam / pi * 180)%%360 else if(tolower(c.system) %in% "art.ec1") ((180 / pi) * -thetam + 150 + 147)%%360
		,sd_wd = sd_wd
		,uprot = u1*cos(phi) + w*sin(phi)
		,vprot = -u*sin(thetam) + v*cos(thetam)
		,wprot = -u1*sin(phi) + w*cos(phi)
		)	

}


pf_transf <- function(u,v,w,P,cw){
	list(
		u = P["up","u"]*u + P["up","v"]*v + P["up","w"]*(w - cw)
		,v = P["vp","u"]*u + P["vp","v"]*v + P["vp","w"]*(w - cw)
		,w = P["wp","u"]*u + P["wp","v"]*v + P["wp","w"]*(w - cw)
		)
}
planar_fit <- function(u,v,w,FUN=MASS::rlm,method=c("Wilczak2001","vanDik2004"),...){
	if(method[1] %in% "Wilczak2001"){
		# Wilczak 2001 equation 39, wm = b0 + b1*um + b2*vm:	
		mod <- FUN(w ~ u + v,...)
		b <- coef(mod)
		names(b) <- c("b0","b1","b2")
		cw <- b["b0"]
		b1 <- b["b1"]
		b2 <- b["b2"]
	} else {
		# van Dik 2004, wm = b1*um + b2*vm:	
		mod <- FUN(w ~ u + v - 1,...)
		b <- coef(mod)
		names(b) <- c("b1","b2")
		cw <- 0
		b1 <- b["b1"]
		b2 <- b["b2"]		
	}

	# equation 42, p3j:
	p31 <- -b1/sqrt(b1^2 + b2^2 + 1)
	p32 <- -b2/sqrt(b1^2 + b2^2 + 1)
	p33 <- 1/sqrt(b1^2 + b2^2 + 1)

	# equation 44, angles:
	# tan_b <- -p32/p33
	sin_b <- -p32/sqrt(p32^2 + p33^2)
	cos_b <- p33/sqrt(p32^2 + p33^2)
	sin_a <- p31
	cos_a <- sqrt(p32^2 + p33^2)

	# equation 2, D + C:
	Cmat <- cbind(c(1,0,0),c(0,cos_b,sin_b),c(0,-sin_b,cos_b))
	Dmat <- cbind(c(cos_a,0,-sin_a),c(0,1,0),c(sin_a,0,cos_a))

	# equation 36, P = t(D)%*%t(C):
	Pmat <- t(Dmat)%*%t(Cmat)
	colnames(Pmat) <- c("u","v","w")
	rownames(Pmat) <- c("up","vp","wp")

	list(Pmat=Pmat,cw=cw,mod=mod,alpha=asin(sin_a)/pi*180,beta=asin(sin_b)/pi*180)
}

twoaxis_rotation_matrix <- function(phi, theta = 0){
	rbind(
		up = c(u = cos(theta) * cos(phi), v = sin(theta) * cos(phi), w = sin(phi)),
		vp = c(u = -sin(theta), v = cos(theta), w = 0),
		wp = c(u = cos(theta) * -sin(phi), v = sin(theta) * -sin(phi), w = cos(phi))
		)
}

pf_rotation_matrix <- function(alpha, beta){
		cos_b <- cos(beta / 180 * pi)
		sin_b <- sin(beta / 180 * pi)
		cos_a <- cos(alpha / 180 * pi)
		sin_a <- sin(alpha / 180 * pi)
		Cmat <- cbind(c(1,0,0),c(0,cos_b,sin_b),c(0,-sin_b,cos_b))
		Dmat <- cbind(c(cos_a,0,-sin_a),c(0,1,0),c(sin_a,0,cos_a))

		# equation 36, P = t(D)%*%t(C):
		Pmat <- t(Dmat)%*%t(Cmat)
		colnames(Pmat) <- c("u","v","w")
		rownames(Pmat) <- c("up","vp","wp")
		Pmat
}

split_index <- function(index_length,n_subint){
	nm <- index_length/n_subint
	fn <- floor(nm)
	res <- index_length - n_subint*fn 
	# len <- rep(fn,n_subint) + c(rep(1,res),rep(0,n_subint-res))
	# mapply(function(x,y)seq.int(x)+y,x=len,y=cumsum(c(0,len[-length(len)])),SIMPLIFY=FALSE)
	len <- rep(fn,n_subint) + c(rep(1,res),rep(0,n_subint-res))
	rep(letters[seq_along(len)],len)
}

wind_statistics <- function(wind,z_canopy,z_sonic){
	Cov_sonic <- cov(as.data.frame(wind[c("up","vp","wp","Tp")]))
	Var_sonic <- diag(Cov_sonic)
	names(Var_sonic) <- c("<u'u'>","<v'v'>","<w'w'>","<T'T'>")
	Cov_sonic <- Cov_sonic[cbind(c("up","up","up","vp","vp","wp"),c("vp","wp","Tp","wp","Tp","Tp"))]
	names(Cov_sonic) <- c("<u'v'>","<u'w'>","<u'T'>","<v'w'>","<v'T'>","<w'T'>")
	suppressWarnings(Ustar <- c(sqrt(-Cov_sonic["<u'w'>"]),use.names = FALSE))
	T_K <- mean(wind$Tm + wind$Tp)
	U <- mean(wind$um + wind$up)
	L <- c(-Ustar^3 * T_K / (0.4 * 9.80620 * Cov_sonic["<w'T'>"]),use.names = FALSE)
	if(!is.na(z_canopy)){
		d <- 2/3 * z_canopy
		suppressWarnings(z0 <- optimize(function(x,ustar,L,z,d,U)abs(U - calcU(ustar, x, L, z-d)),c(0,z_sonic*1.1),ustar=Ustar,L=L,U=U,z=z_sonic,d=d)$minimum)
		if(z0>=z_sonic)z0 <- NA
	} else {
		d <- NA
		z0 <- NA
	}
	as.list(c(Var_sonic,Cov_sonic,Ustar=Ustar,L=L,d=d,z0=z0,sUu=sqrt(Var_sonic["<u'u'>"])/Ustar,sVu=sqrt(Var_sonic["<v'v'>"])/Ustar,sWu=sqrt(Var_sonic["<w'w'>"])/Ustar,U_sonic=U,T_sonic=T_K,U_trend=mean(wind$um),T_trend=mean(wind$Tm)))
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


if(FALSE){	
	file_directory = "U:/Documents/13_Campaigns/EVEMBI/Waedenswil/Daten/Sonic/Sonic1"
	# start_time = "22.02.2018 10:00"
	# end_time = "22.02.2018 20:00"
	avg_period = "10mins"
	start_time = c("22.02.2018 10:00","22.02.2018 12:00","22.02.2018 19:00")
	end_time = c("22.02.2018 10:30","22.02.2018 13:00","22.02.2018 20:00")
	# avg_period = NULL
	tz_sonic = "Etc/GMT-1"
	z_canopy = 0.06
	detrending_method = c(u="linear",v="linear",w="linear",T="linear") # blockAVG/linear/linear_robust/ma_xx (xx = time in seconds)
	n_stationarity_subint = 5
	subint_detrending_method = c(u="blockAVG",v="blockAVG",w="blockAVG",T="blockAVG") # blockAVG/linear/linear_robust/ma_xx (xx = time in seconds)
	hard_flag = list(
		low = c(u=-30,v=-30,w=-10,T=243.15)
		,high = c(u=30,v=30,w=10,T=323.15)
		,window = "1mins" # data points
		,method = "replace" # replace/no replace
		)
	z_sonic = 1.41
	dev_north = 334
	rotation_method = c("two axis","planar fit")
	rotation_args = list(coord_system="WindMaster",phi = NULL,pf_avg_time="30mins",pf_FUN=MASS::rlm,pf_method=c("Wilczak2001","vanDik2004"),pf_N_thresh=10,pf_wd_sectors = c(90,270),pf_U_thresh = 0.5)
	correct_rawdata=c("Gill","none","Nakai2012")
	coord_system="WindMaster"
	create_graphs = FALSE
	write_csv = FALSE
	save_directory = paste0(dirname(file_directory),"/Output_evalSonic")
	add_name = ""
	rawdata_function = readWindMaster_ascii
	as_ibts = TRUE
	variables = c("u","v","w","T")
	covariances = c("u'w'","w'T'")
	data_threshold = 0.9
	lite=TRUE
}

###################################################### Start Evaluation:

evalSonic <- function(
		# config_file = NULL
		file_directory = NULL
		,start_time = NULL
		,end_time = NULL
		,add_time = 0
		,avg_period = NULL
		,tz_sonic = "Etc/GMT-1"
		,z_canopy = NULL
		,detrending_method = c(u="linear",v="linear",w="linear",T="linear") # blockAVG/linear/linear_robust/ma_xx (xx = time in seconds)
		,n_stationarity_subint = 5
		,subint_detrending_method = c(u="blockAVG",v="blockAVG",w="blockAVG",T="blockAVG") # blockAVG/linear/linear_robust/ma_xx (xx = time in seconds)
		,hard_flag = list(
			low = c(u=-30,v=-30,w=-10,T=243.15)
			,high = c(u=30,v=30,w=10,T=323.15)
			,window = "1mins" # seconds OR time string
			,method = "replace" # replace/no replace
			)
		,z_sonic = NULL
		,dev_north = 0
		,rotation_method = c("two axis","planar fit")
		,rotation_args = list(coord_system="WindMaster",phi = NULL,pf_avg_time=avg_period,pf_FUN=MASS::rlm,pf_method=c("Wilczak2001","vanDik2004"),pf_N_thresh=10,pf_wd_sectors = c(0,360),pf_U_thresh = 0)
		,correct_rawdata = c("none", "Gill", "Nakai2012")
		,data_threshold = 0.9
		,create_graphs = FALSE
		,write_csv = FALSE
		,save_directory = paste0(dirname(file_directory),"/Output_evalSonic")
		,add_name = ""
		,rawdata_function = readWindMaster_ascii
		,as_ibts = TRUE
		,variables = c("u","v","w","T")
		,covariances = c("u'w'","w'T'")
		,lite = TRUE
		,asDT = TRUE
		,file_pattern_EVS = "([0-9]{4}-[0-9]{2}-[0-9]{2})_Anemometer_Wetterstation.csv"
		# ,ogives_out = FALSE
	){

	script.start <- Sys.time()

	########################### quick checks:
	if(any(is.null(z_canopy),is.null(z_sonic),is.null(file_directory))){
		stop(paste0("Please provide following argument(s): ",paste(c("\n\t-> 'z_canopy'","\n\t-> 'z_sonic'","\n\t-> 'file_directory'")[which(c(is.null(z_canopy),is.null(z_sonic),is.null(file_directory)))],collapse="")))
	}
	########################### daily files:
	if(deparse(substitute(rawdata_function)) %in% c("readWindMaster_ascii", "readWindMaster_old_ascii")){
		dailyfiles <- dir(file_directory, pattern="(^...|^data_.*)_[0-9]{6,8}_[0-9]{4,6}([.]gz)?$")
		if(length(dailyfiles) == 0){
			stop("Directory '",file_directory,"' doesn't have any sonic files with correct name formatting!")
		}
		if(grepl("^data", dailyfiles[1])){
			starttimes <- fast_strptime(gsub("^data_.*_([0-9]{8}_[0-9]{6})([.]gz)?", "\\1", dailyfiles),format="%Y%m%d_%H%M%S", tz=tz_sonic, lt=FALSE)
		} else {
			starttimes <- fast_strptime(gsub("^..._","",dailyfiles),format="%y%m%d_%H%M", tz=tz_sonic, lt=FALSE)
		}
	} else {
		dailyfiles <- dir(file_directory, pattern = file_pattern_EVS, recursive = TRUE)
		if(length(dailyfiles) == 0){
			stop("Directory '",file_directory,"' doesn't have any sonic files matching provided 'file_pattern_EVS'!")
		}
		starttimes <- fast_strptime(gsub(paste0(".*/",file_pattern_EVS),"\\1",dailyfiles),format="%Y-%m-%d", tz=tz_sonic, lt=FALSE)
	}
	


	########################### time range:
	if(is.null(start_time)){
		start_time <- as.POSIXct(trunc(starttimes[1],"day"))
		cat("argument 'start_time' is not specified. Setting 'start_time' to beginning of first day with measurement data: ",format(start_time), "\n")
	} else {
		start_time <- parse_date_time3(start_time,tz=tz_sonic)
	}
	if(is.null(end_time)){
		end_time <- as.POSIXct(trunc(starttimes[length(starttimes)],"day") + as.difftime(1,units="days"))
		cat("argument 'end_time' is not specified. Setting 'end_time' to end of last day with measurement data: ",format(end_time), "\n")
	} else {
		end_time <- parse_date_time3(end_time,tz=tz_sonic)
	}
	
	if(length(start_time) > 1){
		if(length(start_time) != length(end_time)){
			stop("'start_time' vector length doesn't match length of 'end_time' vector!")
		}
		if(is.null(avg_period)){
			avg_period <- as.numeric(start_time - end_time,units="secs")
		} else {
			avg_period <- parse_time_diff(avg_period)
		}
	} else {
		if(is.null(avg_period))stop("Please provide an averaging time! (argument: 'avg_period')")
		avg_period <- parse_time_diff(avg_period)
	}
	if(avg_period > 2*3600)warning("averaging periods are longer than 2h!")
	difftimes <- as.numeric(diff(c(starttimes,trunc(starttimes[length(starttimes)],"day") + as.difftime(1,units="days"))),units="days")
	difftimes[difftimes > 1] <- 1
	endtimes <- starttimes + as.difftime(difftimes,units="days")
	filerange <- unique(unlist(lapply(cutIntervals(starttimes,endtimes,start_time,end_time),function(x)x[,1])))

	########################### pick correct files:
	if(length(filerange)==0){
        cat("***********\nNo sonic data available between",paste(start_time, end_time, sep =" and "),"\n***********\n")
    }
	read_Files <- dailyfiles[filerange]

	########################### complete argument 'rotation_args'
	ra_list = list(coord_system="WindMaster",phi = NULL,pf_avg_time=avg_period,pf_FUN=MASS::rlm,pf_method=c("Wilczak2001","vanDik2004"),pf_N_thresh=10,pf_wd_sectors = c(0,360),pf_U_thresh = 0)
	rotation_args <- c(rotation_args,ra_list[!(names(ra_list) %in% names(rotation_args))])

	########################### read in data:
	cat("~~~\nReading Raw Files:\n")
	
	dfl <- vector(mode="list",length=length(read_Files))
	for(i in seq_along(read_Files)){
		cat(paste0("\t",i," / ",length(read_Files),": ",read_Files[i],"\n"))
		dfl[[i]] <- rawdata_function(paste(file_directory,read_Files[i],sep="/"))
	}

    OldDevice <- isTRUE(attr(dfl[[1]], "OldDevice"))
	Data <- rbindlist(dfl)

	### free memory
	rm(dfl)
	{xalt <- matrix(0,2,3)
	xneu <- gc()
	while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}}

	########################### add time shift
	Data[,Time := Time + parse_time_diff(add_time)]

	########################### extract data within specified time range:
	Hz <- Data[, Hz[1]]
	Data[, Hz := NULL]
	Data <- Data[findI_st(Time,start_time,end_time) > 0,]

    # return null if no data within time range
    if (nrow(Data) == 0L) {
        cat('No data within given timerange\n')
        return(NULL)
    }

	# raw data quality control I, i.e. hard flags = physical range
	# --------------------------------------------------------------------------
	cat("~~~\nchecking hard flags...\n")
	if(!isFALSE(hard_flag)){
		lim_range <- rbind(lower=hard_flag$low,top=hard_flag$high)
		hf_vars <- colnames(lim_range)
		Data[,(hf_vars) := {
			# browser()
			H.flags(mget(hf_vars),Time,1/Hz,lim_range,hard_flag$window,hard_flag$method)
		}]
	}
	# drop NA's
	Data <- na.omit(Data)

	########################### correct wind speed for angle of attack:
	# Data[,c("u0","v0","w0") := .(u,v,w)]
    if(missing(correct_rawdata) && OldDevice){
        correct_rawdata <- "Gill"
    }
	switch(pmatch(correct_rawdata[1],c("Gill","Nakai2012","none")),
		{
			cat("~~~\nCorrecting for Gill software bug. Data is corrected as proposed by Gill instruments...\n")
			Data[,w:=w*ifelse(w<0,1.289,1.166)]
			# change data names
		},
		{
			cat("~~~\nAOA correction: data is corrected as proposed by Nakai and Shimoyama 2012...\n")
			cat("Sourcing C++ Code. This might take a while...\n")
			require(Rcpp)
			RcppText <- "
				#include <Rcpp.h>

				// #include <aoa2012.h>

				double sinerr(double x, double wd)
				{
				// Sine correction function phi_sr(alpha, gamma)
				// Eqs. (10), (13), and (14) of Nakai and Shimoyama (submitted)

					double a1, a2, a3, a4, a5;
					double b1, b2, b3, b4, b5;
					double f_aoa, x_orig, wd_orig;
					
					a1 = -3.19818998552857E-10;
					a2 = -2.69824417931343E-8;
					a3 = 4.16728613218081E-6;
					a4 = 4.85252964763967E-4;
					a5 = 1.67354200080193E-2;
					b1 = 5.92731123831391E-10;
					b2 = 1.44129103378194E-7;
					b3 = 1.20670183305798E-5;
					b4 = 3.92584527104954E-4;
					b5 = 3.82901759130896E-3;
					
					x_orig = x;
					wd_orig = wd; 
					
					if (x > 0) {
						x *= -1;
						wd += 180;
					}
					
					f_aoa = a1 * pow(x, 5) + a2 * pow(x, 4) + a3 * pow(x, 3) + a4 * pow(x, 2) + a5 * x + 1;
					f_aoa -= sin(3 * wd * PI/180) * (b1 * pow(x, 5) + b2 * pow(x, 4) + b3 * pow(x, 3) + b4 * pow(x, 2) + b5 * x);
					
					x = x_orig;
					wd = wd_orig;
					
					return f_aoa;
				}

				double coserr(double x, double wd)
				{
				// Cosine correction function phi_cr(alpha, gamma)
				// Eqs. (11), (12), (15), and (16) of Nakai and Shimoyama (submitted)

					double c1, c2, c3, c4, c5;
					double d1, d2, d3, d4, d5;
					double f_aoa, x_orig, wd_orig;
					
					c1 = -1.20804470033571E-9;
					c2 = -1.58051314507891E-7;
					c3 = -4.95504975706944E-6;
					c4 = 1.60799801968464E-5;
					c5 = 1.28143810766839E-3;
					d1 = 2.2715401644872E-9;
					d2 = 3.85646200219364E-7;
					d3 = 2.03402753902096E-5;
					d4 = 3.94248403622007E-4;
					d5 = 9.18428193641156E-4;

					x_orig = x;
					wd_orig = wd; 
					
					if (x > 0) {
						x *= -1;
						wd += 180;
					}
					
					if (x < -70)	x = -70;
					
					f_aoa = c1 * pow(x, 5) + c2 * pow(x, 4) + c3 * pow(x, 3) + c4 * pow(x, 2) + c5 * x + 1;
					f_aoa += sin(3 * wd * PI/180) * (d1 * pow(x, 5) + d2 * pow(x, 4) + d3 * pow(x, 3) + d4 * pow(x, 2) + d5 * x);
					
					x = x_orig;
					wd = wd_orig;
					
					return f_aoa;
				}

				double gx(double x, double wd, double a)
				{
				// Equation to solve
					
					return atan(a * coserr(x, wd) / sinerr(x, wd)) * 180 / PI;
					
				}

				double Steffensen(double x, double wd, double a)
				{
				// Appendix A of Nakai et al.(2006)

					double x0, x1, x2, x3;
					double key;
					
					x0 = x;
					while (1) {
						x1 = gx(x0, wd, a);
						x2 = gx(x1, wd, a);
						x3 = x2;	// output value if `break' out of the loop
						key = x2 - 2 * x1 + x0;
						if (fabs(key) < 0.01) break;
						x3 = x0 - pow(x1 - x0, 2) / key;
						x0 = x3;
					}
					
					return x3;
				}

				// [[Rcpp::export]]
				Rcpp::List angleAttack(std::vector<double> &uIn, std::vector<double> &vIn, std::vector<double> &wIn) {
					const int N = uIn.size();
					Rcpp::NumericVector Uout(N);
					Rcpp::NumericVector Vout(N);
					Rcpp::NumericVector Wout(N);
					double aoa, sin_err, cos_err, wd, ws;
					int i;	
					for(int i = 0; i < N; i++){
						if (wIn[i] == 0){
							aoa = 0.0;
						} else {
							ws = sqrt(pow(uIn[i], 2) + pow(vIn[i], 2));
							
							if (ws == 0) {
								if (wIn[i] >= 0) aoa = 90;
								if (wIn[i] < 0)  aoa = -90;
							} else {
								aoa = atan(wIn[i]/ws) * 180 / PI;
							}
							
							// Wind direction
							if (ws == 0) {
								wd = 0;
							} else if (vIn[i] >= 0) {
								wd = 180 - acos(uIn[i] / ws) * 180 / PI;
							} else {
								wd = 180 + acos(uIn[i] / ws) * 180 / PI;
							}
							
							// Steffensen's method --- find out true AoA
							if (ws != 0) aoa = Steffensen(aoa, wd, wIn[i]/ws);
						}
						// sine error calculation
						sin_err = sinerr(aoa, wd);	
						// cosine error calculation
						cos_err = coserr(aoa, wd);
						Uout[i] = uIn[i] / cos_err;
						Vout[i] = vIn[i] / cos_err;
						Wout[i] = wIn[i] / sin_err;
					}
					return Rcpp::List::create(
						Rcpp::_[\"uCorr\"] = Uout,
						Rcpp::_[\"vCorr\"] = Vout,
						Rcpp::_[\"wCorr\"] = Wout);
				}"

			sourceCpp(code=RcppText)
			# Attack!
			cat("Correcting for Angle of Attack...\n")
			Data[,c("u","v","w"):=angleAttack(u,v,w)]
		},
		cat("~~~\nno AOA correction is done...\n")
	)

	########################### define period bins:
	Data[,BinPeriod := floor(as.numeric(Time-start_time[1],units="secs")/avg_period)]

	### remove incomplete intervals
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Data <- Data[,keep := .N/(avg_period * Hz) > data_threshold,by = BinPeriod][(keep)][,keep := NULL]

	# calculate wind direction, rotate u, v, w, possibly detrend T (+ u,v,w)
	# --------------------------------------------------------------------------
	c.system <- tolower(rotation_args$"coord_system")
	# planar fit?
	if(rotation_method[1] %in% "planar fit"){
		cat(paste0("~~~\nApply planar fit method (",rotation_args$pf_method,") based on"),if(is.character(rotation_args$pf_avg_time))rotation_args$pf_avg_time else paste(parse_time_diff(rotation_args$pf_avg_time),"secs"),"intervals\n")
		pf_avg_time <- parse_time_diff(rotation_args$pf_avg_time)
		# PF periods:
		Data[,PFPeriod := floor(as.numeric(Time-start_time[1],units="secs")/pf_avg_time)]
		### tag 'complete' intervals
		Data[,keepPF := .N/(pf_avg_time * Hz) > data_threshold,by = PFPeriod]
		# get Uavg & provis. WD
		Data[,c("um","vm","wm","Uavg","WDprov") := {
			um <- mean(u)
			vm <- mean(v)
			wm <- mean(w)
			.(um,vm,wm,sqrt(um^2 + vm^2 + wm^2),if(c.system %in% "windmaster") (180 - atan2(vm,um) / pi * 180)%%360 else if(c.system %in% "art.ec1") ((180 / pi) * -atan2(vm,um) + 150 + 147)%%360)
		},by = PFPeriod]
		Data[,WDprov := (WDprov + dev_north) %% 360]
		# fit to plane
		# browser()
		wd_sec <- unique(c(0,sort(rotation_args$pf_wd_sectors),360))
		no_sec <- length(wd_sec) - 2
		wd_sort <- unique(sort(rotation_args$pf_wd_sectors))
		wd_sec <- c(wd_sort[length(wd_sort)],wd_sort)
		no_sec <- length(wd_sort)
		U_thres <- rotation_args$pf_U_thresh
		if(is.null(U_thres)){
			U_thres <- 0
		}
		Data[,c("alpha","beta","w_bias") := NA_real_]
		cat(sprintf("\tU threshold = %1.2f m/s (%1.1f%% of data set)\n",U_thres,Data[Uavg > U_thres,.N]/nrow(Data)*100))
		# browser()
		if(no_sec > 1){
			cat("\tFollowing wind sectors have been defined:\n")
			for(i in seq.int(no_sec)){
				if(i == 1){
					cat(paste0("\tSector ",i,": (",wd_sec[no_sec + 1],"\u00B0 to ",wd_sec[2],"\u00B0] (",Data[(WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1]) & Uavg > U_thres,length(unique(PFPeriod))]," intervals",sprintf(" %1.1f%% of wind speed subset)",Data[(WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1]) & Uavg > U_thres,.N]/Data[Uavg > U_thres,.N]*100),"\n"))
					if(Data[(WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1]) & Uavg > U_thres,.N] < rotation_args$pf_N_thresh){
						if(Data[(WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1]) & Uavg > U_thres,.N] == 0){
							cat("\t\t|-> No data in this sector.\n")
						} else {
							cat(paste0("\t\t|-> Not enough data (",Data[(WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1]) & Uavg > U_thres,.N],") for planar fit with threshold of ",rotation_args$pf_N_thresh," (applying two axis rotation\n"))
						}
					} else {
						# do subset fit
						PF <- Data[(WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1]) & Uavg > U_thres,.(um=um[1],vm=vm[1],wm=wm[1]),by = PFPeriod][,{
							list2env(planar_fit(um,vm,wm,FUN=rotation_args$pf_FUN,method=rotation_args$pf_method)[c("Pmat","cw","alpha","beta")])
						}]
						# transform:
						Data[WDprov <= wd_sec[2] | WDprov > wd_sec[no_sec + 1],c("u","v","w","alpha","beta","w_bias") := {
							# browser()
							c(pf_transf(u,v,w,with(PF,Pmat),with(PF,cw)),list(with(PF,alpha),with(PF,beta),with(PF,cw)))
						}]
					}
				} else {
					cat(paste0("\tSector ",i,": (",wd_sec[i],"\u00B0 to ",wd_sec[i + 1],"\u00B0] (",Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > U_thres,length(unique(PFPeriod))]," intervals",sprintf(" %1.1f%% of wind speed subset)",Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > U_thres,.N]/Data[Uavg > U_thres,.N]*100),"\n"))
					if(Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > U_thres,.N] < rotation_args$pf_N_thresh){
						if(Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > U_thres,.N] == 0){
							cat("\t\t|-> No data in this sector.\n")
						} else {
							cat(paste0("\t\t|-> Not enough data (",Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > U_thres,.N],") for planar fit with threshold of ",rotation_args$pf_N_thresh," (applying two axis rotation for this sector)\n"))
						}
					} else {
						# do subset fit
						PF <- Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > U_thres,.(um=um[1],vm=vm[1],wm=wm[1]),by = PFPeriod][,{
							list2env(planar_fit(um,vm,wm,FUN=rotation_args$pf_FUN,method=rotation_args$pf_method)[c("Pmat","cw","alpha","beta")])
						}]
						# transform:
						Data[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1],c("u","v","w","alpha","beta","w_bias") := {
							# browser()
							c(pf_transf(u,v,w,with(PF,Pmat),with(PF,cw)),list(with(PF,alpha),with(PF,beta),with(PF,cw)))
						}]
					}
				}
			}
			# browser()
		} else {
			cat("No wind sectors have been defined.  (",Data[Uavg > U_thres,length(unique(PFPeriod))]," intervals available)\n")
			if(Data[Uavg > U_thres,.N] < rotation_args$pf_N_thresh){
				cat(paste0("\t-> Not enough data (",Data[Uavg > U_thres,.N],") for planar fit with threshold of ",rotation_args$pf_N_thresh," (applying two axis rotation)\n"))
			} else {
				# do subset fit
				PF <- Data[Uavg > U_thres,.(um=um[1],vm=vm[1],wm=wm[1]),by = PFPeriod][,{
					list2env(planar_fit(um,vm,wm,FUN=rotation_args$pf_FUN,method=rotation_args$pf_method)[c("Pmat","cw","alpha","beta")])
				}]
				# transform:
				Data[,c("u","v","w","alpha","beta","w_bias") := .(pf_transf(u,v,w,with(PF,Pmat),with(PF,cw)),with(PF,alpha),with(PF,beta),with(PF,cw))]
			}
		}
		### subtract w_avg (is done by detrending) and remove keepPF column:
		Data[,c("keepPF","WDprov","Uavg") := NULL]
	}

	# rotate coordinate system
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	cat("~~~\nrotating and detrending data...\n")
	Data[,
		c("phi","wd","sd_wd","uprot","vprot","wprot") := rotate_twoaxis(u, v, w, T, phi = if(rotation_method[1] %in% "two axis") rotation_args$phi else if(is.na(alpha[1])) NULL else 0, c.system = tolower(rotation_args$coord_system))
		,by = BinPeriod]
	### correct for sonic north deviation
	Data[,wd := (wd + dev_north) %% 360]

	### Write to output & free memory
	Out <- Data[,{
		# browser()
		.(
			start_interval = start_time[1] + .BY[[1]]*avg_period
			,end_interval = start_time[1] + (.BY[[1]] + 1)*avg_period
			,n_data = .N
            ,Hz = Hz
			,z_sonic = z_sonic
			,z_canopy = NA_real_
			,d = NA_real_
			,WD = wd[1]
			,sd_WD = sd_wd[1]
			,Ustar = NA_real_
			,L = NA_real_
			,z0 = NA_real_
			,sUu = NA_real_
			,sVu = NA_real_
			,sWu = NA_real_
			,phi = phi[1]
			,alpha = if(exists("w_bias")) alpha[1] else NA_real_
			,beta = if(exists("w_bias")) beta[1] else NA_real_
			,w_bias = if(exists("w_bias")) w_bias[1] else 0
			,T_sonic = NA_real_
			,U_sonic = NA_real_
			,up.up = NA_real_
			,vp.vp = NA_real_
			,wp.wp = NA_real_
			,Tp.Tp = NA_real_
			,up.vp = NA_real_
			,up.wp = NA_real_
			,up.Tp = NA_real_
			,vp.wp = NA_real_
			,vp.Tp = NA_real_
			,wp.Tp = NA_real_
			,sub_z_canopy = NA_real_
			,sub_d = NA_real_
			,sub_WD = NA_real_
			,sub_sd_WD = NA_real_
			,sub_Ustar = NA_real_
			,sub_L = NA_real_
			,sub_z0 = NA_real_
			,sub_sUu = NA_real_
			,sub_sVu = NA_real_
			,sub_sWu = NA_real_
			,sub_phi = NA_real_
			,sub_up.up = NA_real_
			,sub_vp.vp = NA_real_
			,sub_wp.wp = NA_real_
			,sub_Tp.Tp = NA_real_
			,sub_up.vp = NA_real_
			,sub_up.wp = NA_real_
			,sub_up.Tp = NA_real_
			,sub_vp.wp = NA_real_
			,sub_vp.Tp = NA_real_
			,sub_wp.Tp = NA_real_
			)
	},by=BinPeriod][,BinPeriod := NULL]
	Data[,c("phi","wd","sd_wd") := NULL]

	# detrend variables
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Data[,c("up","vp","wp","Tp","um","vm","wm","Tm") := {
		detrend(uprot,vprot,wprot,T,method=detrending_method[c("u","v","w","T")],Hz_ts=Hz)
	},by = BinPeriod]

	# browser()
	if(is.ibts(z_canopy)){
		tzone(z_canopy) <- tz_sonic
		z_canopy <- pool(z_canopy,st.to = Out[,format(start_interval)],et.to = Out[,format(end_interval)])[,1,drop=TRUE,keepAtts=FALSE]
	} else {
		z_canopy <- rep(z_canopy,nrow(Out))
	}
	names(z_canopy) <- as.character(Data[,unique(BinPeriod)])
	# calculate wind statistics etc.
	# -------------------------------------------------------------------------- 
	cat("~~~\ncalculating wind statistics\n")
	Out[,c("z_canopy","d","Ustar","L","z0","sUu","sVu","sWu","T_sonic","U_sonic","up.up","vp.vp","wp.wp","Tp.Tp","up.vp","up.wp","up.Tp","vp.wp","vp.Tp","wp.Tp")] <- Data[,{
		# browser()
		WS <- wind_statistics(.(up=up,vp=vp,wp=wp,Tp=Tp,um=um,Tm=Tm),z_canopy[as.character(.BY[[1]])],z_sonic)
		.(
			z_canopy = z_canopy[as.character(.BY[[1]])]
			,d = WS$d
			,Ustar = WS$Ustar
			,L = WS$L
			,z0 = WS$z0
			,sUu = WS$sUu
			,sVu = WS$sVu
			,sWu = WS$sWu
			,T_sonic = WS$T_sonic
			,U_sonic = WS$U_sonic
			,up.up = WS$"<u'u'>"
			,vp.vp = WS$"<v'v'>"
			,wp.wp = WS$"<w'w'>"
			,Tp.Tp = WS$"<T'T'>"
			,up.vp = WS$"<u'v'>"
			,up.wp = WS$"<u'w'>"
			,up.Tp = WS$"<u'T'>"
			,vp.wp = WS$"<v'w'>"
			,vp.Tp = WS$"<v'T'>"
			,wp.Tp = WS$"<w'T'>"
			)
	},by=BinPeriod][,BinPeriod := NULL]



	if(create_graphs){
		if(!dir.exists(save_directory))dir.create(save_directory)
		### time series
		cat("~~~\nplotting...\n")
		require(lattice)
		covariances_cols <- c("u'w'" = blues9[7], "w'T'" = "orange")
		covariances_variables <- list("u'w'" = c("u","w"),"w'T'" = c("w","T"))
		Data[,{
			Int_End <- Time[.N]
			Int_Start <- Time[1]
			Int_Time <- as.numeric(Int_End - Int_Start,units="secs")
			filename <- format(Int_End,format="%Y%m%d") 
			time2 <- format(Int_End,format="%H%M%S")
			plotname <- paste("timeseries", filename, time2, sep="-") 
			jpeg(file=file.path(save_directory,paste0(plotname,".jpg")),width=1350, height=900, quality=60)
				ts_plot <- plot.tseries(.SD)
				print(ts_plot)
			dev.off()

			### cov
			# calculate covariances with fix lag time:
			# -------------------------------------------------------------------------- 
			freq <- seq(floor(.N/2))*2/.N
			FFTs <- lapply(list(u=up,w=wp,T=Tp),function(x)fft(x)/.N)

			if(.N%%2){
				Covars <- lapply(covariances_variables,function(i,x){
					Re(fft(Conj(FFTs[[i[2]]]) * FFTs[[i[1]]], inverse=TRUE))[c(((.N+1)/2+1):.N,1:((.N+1)/2))]*.N/(.N-1)
					},x=FFTs)
			} else {
				Covars <- lapply(covariances_variables,function(i,x){
					Re(fft(Conj(FFTs[[i[2]]]) * FFTs[[i[1]]], inverse=TRUE))[c((.N/2+1):.N,1:(.N/2))]*.N/(.N-1)
					},x=FFTs)
			}
			### ogive
			# cospectra for fixed & dynamic lags
			# ------------------------------------------------------------------------ 			
			Cospec_fix <- mapply(function(i,x1,x2){
					xs <- fft(x2[,i[2]])/.N
					Re(Conj(xs) * x1[[i[1]]])[seq(.N/2)+1]*.N/(.N-1)*2
				},i=covariances_variables,MoreArgs=list(x1=FFTs,x2=data.frame(u=up,w=wp,T=Tp)),SIMPLIFY=FALSE)

			
			#*** hac5: *** ogives for fixed & dynamic lags 
			# ------------------------------------------------------------------------ 
			Ogive_fix <- lapply(Cospec_fix,function(x)rev(cumsum(rev(x))))

			# plot and save flux evaluation...
			# ------------------------------------------------------------------------
			for(i in names(covariances_variables)){
				# i <- "w'TDL CH4'"
				plotname <- paste("plots",filename,time2,i,sep="-")
				jpeg(file=file.path(save_directory,paste0(plotname,".jpg")),width=1350, height=600, quality=60)
					par(mfrow=c(1,2))
					# ----------------------- Covariance -------------------------------------
					plot_covfunc(Covars[[i]],Int_Time,ylab=i, xlim = c(-50,50), cx=1.5, cxmt=1.25, cl=covariances_cols[i])
					# ---------------------- Co-Spec/Ogive fix lag -----------------------------------
					plot_cospec_ogive(Ogive_fix[[i]],Cospec_fix[[i]],freq,ylab=paste0("ogive (fix lag) of ",i),cx=1.5,col=covariances_cols[i])
					title(paste0(i," flux ",format(Int_Start,format="(%H:%M:%S")," - ",format(Int_End,format="%H:%M:%S)"))
						,outer=TRUE,line=-1)
				dev.off()	
			}





		}, by = BinPeriod]
	}


	### Free memory:
	Data[,c("up","vp","wp","Tp","um","Tm","uprot","vprot","wprot","vm","wm") := NULL]

	# sub-int: sub-interval calculations
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	cat("~~~\nsub-intervals\n")
	Data[,BinSub := paste0(BinPeriod,split_index(.N,n_stationarity_subint)),by = BinPeriod]

	# sub-int: rotate coordinate system
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	cat("\t |-> rotating and detrending data...\n")
	Data[,
		c("sub_phi","sub_wd","sub_sd_wd","sub_uprot","sub_vprot","sub_wprot") := rotate_twoaxis(u, v, w, T, phi = if(rotation_method[1] %in% "two axis") rotation_args$phi else if(is.na(alpha[1])) NULL else 0, c.system = tolower(rotation_args$coord_system))
		,by = BinSub]
	Out[,c("sub_WD","sub_sd_WD","sub_phi")] <- Data[,.(mean(sub_wd),mean(sub_sd_wd),mean(sub_phi)),by=BinPeriod][,BinPeriod := NULL]
	Data[,c("sub_phi","sub_wd","sub_sd_wd") := NULL]
	
	# sub-int: detrend variables
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	Data[,c("sub_up","sub_vp","sub_wp","sub_Tp","sub_um","sub_vm","sub_wm","sub_Tm") := {
		detrend(sub_uprot,sub_vprot,sub_wprot,T,method=detrending_method[c("u","v","w","T")],Hz_ts=Hz)
	},by = BinSub][,c("sub_uprot","sub_vprot","sub_wprot","sub_vm","sub_wm") := NULL]

	# browser()
	# sub-int: wind statistics
	cat("\t |-> calculating wind statistics\n")
	Out[,paste0("sub_",c("z_canopy","d","Ustar","L","z0","sUu","sVu","sWu","up.up","vp.vp","wp.wp","Tp.Tp","up.vp","up.wp","up.Tp","vp.wp","vp.Tp","wp.Tp"))] <- Data[,{
		# browser()
		WS <- wind_statistics(.(up=sub_up,vp=sub_vp,wp=sub_wp,Tp=sub_Tp,um=sub_um,Tm=sub_Tm),z_canopy[as.character(BinPeriod[1])],z_sonic)
		.(
			z_canopy = z_canopy[as.character(BinPeriod[1])]
			,d = WS$d
			,Ustar = WS$Ustar
			,L = WS$L
			,z0 = WS$z0
			,sUu = WS$sUu
			,sVu = WS$sVu
			,sWu = WS$sWu
			,up.up = WS$"<u'u'>"
			,vp.vp = WS$"<v'v'>"
			,wp.wp = WS$"<w'w'>"
			,Tp.Tp = WS$"<T'T'>"
			,up.vp = WS$"<u'v'>"
			,up.wp = WS$"<u'w'>"
			,up.Tp = WS$"<u'T'>"
			,vp.wp = WS$"<v'w'>"
			,vp.Tp = WS$"<v'T'>"
			,wp.Tp = WS$"<w'T'>"
			,BinPeriod = BinPeriod[1]
			)
	},by=BinSub][,{
		# browser()
		as.list(colMeans(.SD[,-1]))
	},by = BinPeriod][,BinPeriod := NULL]


	rm(Data)
	for(i in 1:20)gc()


	if(write_csv){
		if(!dir.exists(save_directory))dir.create(save_directory)
		Out_csv <- Out[,c("st","et","tz") := .(format(start_interval),format(end_interval),tz_sonic)]
		Out_csv[,c("start_interval","end_interval") := .(st,et)][,c("st","et") := NULL]
		setcolorder(Out_csv,names(Out_csv)[c(1:2,length(Out_csv),3:(length(Out_csv)-1))])
		fwrite(Out_csv,file=paste0(save_directory,"/EvalSonic_",format(start_time[1],format="%Y-%m-%d"),"_to_",format(end_time[length(end_time)],format="%Y-%m-%d"),format(Sys.time(),format="_%y%m%d%H%M"),add_name,".csv"))
	}

	# #################################### END VERSION HISTORY #################################### #
	cat("************************************************************\n") 
	cat("operation finished @", format(Sys.time(), "%d.%m.%Y %H:%M:%S"),"time elapsed: ", difftime(Sys.time(), script.start, unit="mins"),"minutes\n")
	cat("************************************************************\n")  

	if(!asDT){
		setDF(Out)
	}
	
	Out

}

###### plotting functions


plot.tseries <- function(DT, selection = c("u", "v", "w", "T"), units = c("m/s", "m/s", "m/s", "deg C"), color = c(blues9[c(7,7,7)], "orange")){
	# browser()
	dat <- setNames(as.data.frame(DT[,c("Time", "uprot", "vprot", "wprot", "T")]), c("st","u", "v", "w", "T"))
	msg <- DT[,paste(c(format(Time[1],"%d.%m.%Y")," - time series"),collapse="")]
	dat2 <- melt(dat,id="st")
	### get trends:
	dat3 <- setNames(as.data.frame(DT[,c("Time", "um", "vm", "wm", "Tm")]), c("st","u", "v", "w", "T"))
	### melt and add trends:
	dat4 <- melt(dat3,id="st",value.name="trend")
	dat2[,"trend"] <- dat4[,"trend"]

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
		strip=FALSE, layout=c(1,ncol(dat)-1), between=list(x=0,y=1), subscripts=TRUE, lwd=rep(1,length(color)), lty=rep(1,length(color)), col=color,
		panel=function(x, y, ...) {
			panel.grid(h=-1, v=-1, lty=3, col="gray80")
			y2 <- dat2[list(...)$subscripts,"trend"]
			panel.xyplot(x,y,...)
			panel.xyplot(x,y2,type="l",lwd=1.5, lty=3, col="gray60")
			# panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightblue")
			# panel.xyplot(x,y2,type="l",lwd=2, lty=2, col=blues9[5])
			# panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightgrey")
		}
	)  	
}

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

plot_covfunc <- function(cov_func,avg_t, ylab=NULL, xlim = NULL, cx=1.5, cxmt=1.25, cl="black"){
	# browser()
	n <- length(cov_func)
	midP <- floor((n + 1)/2)
	fix_cov <- cov_func[midP]

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
	abline(v=0,lty=3,col="lightgrey")
	lines(x, cov_func, lwd=2, col=cl)
	mtext(substitute(paste(italic(tau)==a,"s: ",y),list(y=sprintf("%1.8f",fix_cov),a=sprintf("%1.1f",0))),side=3,cex=cxmt)
	
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
	pxlims <- rep(pxlim,each=9) + log10(seq(9))
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

reduce_cospec <- function(cospec,freq,length.out=100){
	log_freq <- log(freq)
	log_cuts <- seq(min(log_freq),max(log_freq),length.out=length.out+1)
	ind_cuts <- findInterval(log_freq,log_cuts,rightmost.closed=TRUE)
	freq_out <- exp(log_cuts[-1] - diff(log_cuts))
	list(cospec=tapply(cospec,ind_cuts,sum),freq=freq_out[unique(ind_cuts)],d_freq=diff(exp(log_cuts))[unique(ind_cuts)])
}
