
library(data.table)
library(ibts)
library(Rcpp)

# check (and optionally replace) 'hard' data limits
check_limits <- function(input, time, d_t, limits, wind = 500, hflg.met = "norepl"){ 
	nms <- names(input)
	missing <- !(nms %in% colnames(limits))
	if (length(nms[missing]) > 0) {
        stop(paste0("physical limits missing for: ", paste(nms[missing], collapse=", ")),
            call. = TRUE, domain = NULL)
    }
	limits <- limits[, nms]
	### replace values outside of limits by NAs
	dat <- mapply(function(x, lo, hi) {
		x[x < lo | x > hi] <- NA
		x
	}, x = input, lo = limits["lower", ], hi = limits["top", ], SIMPLIFY = FALSE)
	hflgs <- unlist(lapply(dat, \(x) sum(is.na(x))))
	if (pmatch(hflg.met, "replace", nomatch = 0) && any(hflgs > 0)) {
		cat("Replacing flagged values by window median...\n")
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
	} else {
		cat(
            paste0(
                'number of values outside of hard limits (set to NA)\n*~~~~*\n',
                paste(
                    sprintf('%s: %i', names(dat[hflgs > 0]), hflgs[hflgs > 0]), 
                    collapse = '\n'),
                "\n*~~~~*\n"
            )
        )
    }
	dat
}

# C++ helper function used in hard limits function
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

# filter functions for higph-pass filtering time series
# (copied from minidoas scripts)
filter_list <- list(
    # Blackman-Harris:
    "BmHarris" = function(n, ...){
      N <- n - 1
      x <- 0.35875 - 
        0.48829 * cos(2 * pi * seq.int(0, N) / N) + 
        0.14128 * cos(4 * pi * seq.int(0, N) / N) - 
        0.01168 * cos(6 * pi * seq.int(0, N) / N)
      x / sum(x)
    },
    # Hann:
    "Hann" = function(n, ...){
      N <- n - 1 
      sin(pi * seq.int(0, N) / N) ^ 2 / (0.5 * N)
    },
    # Hamming:
    "Hamming" = function(n, p1 = 25 / 46, p2 = 21 / 46, ...){
      N <- n - 1
      x <- p1 - p2 * cos(2 * pi * seq.int(0, N) / N)
      x / sum(x)
    },
    # Blackman:
    "Blackman" = function(n, p1 = 0.42, p2 = 0.5, p3 = 0.08, ...){
      N <- n - 1
      x <- p1 - 
        p2 * cos(2 * pi * seq.int(0, N) / N) + 
        p3 * cos(4 * pi * seq.int(0, N) / N)
      x/sum(x)
    },
    # Blackman-Nuttall:
    "BmNuttall" = function(n, ...){
      N <- n - 1
      x <- 0.3635819 - 
        0.4891775 * cos(2 * pi * seq.int(0, N) / N) + 
        0.1365995 * cos(4 * pi * seq.int(0, N) / N) - 
        0.0106411 * cos(6 * pi * seq.int(0, N) / N)
      x / sum(x)
    },
    # Rectangular (Moving Average):
    "Rect" = function(n, ...){
      rep(1 / n, n)
    },
    # Exponential:
    "Exp" = function(n, p1 = n / 2, ...){
      N <- (n - 1) / 2
      x <- exp(p1 * sqrt(1 - (seq.int(-N, N) / N) ^ 2)) / exp(p1)
      x / sum(x)
    },
    # Poisson:
    "Poisson" = function(n, p1 = n / 2 / 6.9, ...){
      N <- n - 1
      x <- exp(-abs(seq.int(0, N) - N / 2) / p1)
      x / sum(x)
    }
)

# detrending time series
trend <- function(y, method = c("blockAVG", "linear", "linear_robust", "ma_360"), Hz_ts = 10) {
	n <- length(y)
	method <- method[1]
	switch(method, 
		"blockAVG" = {
			my <- .Internal(mean(y))
			fitted <- rep(my, n) 
			list(
				coefficients = c(intercept = my, slope = 0)
				, fitted = fitted
				, residuals = y - fitted
				) 
		}
		, "linear" = {
			my <- .Internal(mean(y))
			mx <- (n + 1) / 2
			xstr <- (x <- seq.int(n)) - mx
			b <- sum(xstr * (y - my)) / (n - 1) / (sum(xstr ^ 2) / (n - 1))
			a <- my - mx * b
			fitted <- a + b * x
			list(
				coefficients = c(intercept = a, slope = b)
				, fitted = fitted
				, residuals = y - fitted
				) 
		}
		, "linear_robust" = {
			mod <- MASS::rlm(y ~ seq.int(n))
			cfs <- mod$coefficients
			a <- cfs[1]
			b <- cfs[2]
			fitted <- mod$fitted
			if (!mod$converged) {
				stop("robust linear regression did not converge!")
			}
			list(
				coefficients = c(intercept = a, slope = b)
				, fitted = fitted
				, residuals = y - fitted
				) 
		}
		,{
            if (grepl('^ma_', method)) {
                ma <- round(as.numeric(sub("ma_", "", method)) * Hz_ts)
                fitted <- caTools::runmean(y, ma, "C", "mean")
                list(
                    coefficients = c(intercept = NA, slope = NA)
                    , fitted = fitted
                    , residuals = y - fitted
                ) 
            } else if (sub('_\\d+$', '', method) %in% names(filter_list)) {
                win <- round(as.numeric(sub(".*_", "", method)) * Hz_ts)
                filter_name <- sub('_\\d+$', '', method)
                # lowpass filter
                filt <- filter_list[[filter_name]](win)
                ly <- length(y)
                ihalf <- 1:(win / 2)
                yh1 <- mean(y[ihalf])
                yh2 <- mean(y[ly + 1 - ihalf])
                ym <- y - (yh2 - yh1)
                yp <- y - (yh1 - yh2)
                z <- c(ym[ly + 1 - ihalf], y, yp[ihalf])
                fitted <- stats::filter(z, filt, 'convolution', 2, circular = FALSE)[1:ly +
                    ihalf[length(ihalf)]]
                list(
                    coefficients = c(intercept = NA, slope = NA)
                    , fitted = fitted
                    , residuals = y - fitted
                )
            } else {
                stop(
                    'detrending method not valid! Should be one of "blockAVG",',
                    ' "linear", "linear_robust", "ma_xxx" or ',
                    'any valid filter name yyy_xxx from the "filter_list" where',
                    ' xxx represents the moving average or filter range in seconds'
                )
            }
		}
	)
}

# get WD from different coord. systems
csystem_wd <- function(c.system, um, vm) {
    if (missing(vm)) {
        theta <- um
    } else {
        theta <- atan2(vm, um)
    }
	switch(tolower(c.system)
        , "windmaster" = (180 - theta / pi * 180) %% 360 
        , "hs-wauwilermoos" = (150 - theta / pi * 180) %% 360
        , "art.ec1" = ((180 / pi) * -theta + 150 + 147) %% 360
	)
}

# two-axis rotation
rotate_twoaxis <- function(u, v, w, phi = NULL, c.system = "Windmaster") {
    # wind direction
	thetam <- atan2(.Internal(mean(v)), .Internal(mean(u)))
	# horizontal rotation
	u1 <- u * cos(thetam) + v * sin(thetam)
    # apply vertival rotation?
	if (is.null(phi)) {
		phi <- atan2(.Internal(mean(w)), .Internal(mean(u1)))
	}
    # output
	out <- list(
		wd = csystem_wd(c.system, thetam)
		, phi = phi
		, urot = u1 * cos(phi) + w * sin(phi)
		, vrot = -u * sin(thetam) + v * cos(thetam)
		, wrot = -u1 * sin(phi) + w * cos(phi)
	)
}

# planar fit
planar_fit <- function(x, coord_system, method = c('Wilczak2001', 'vanDik2004')[1], 
    avg_time = '30mins', wd_sectors = c(0, 360), u_thresh = 0, n_thresh = 10, 
    reg_fun = MASS::rlm, start_time, data_threshold, Hz, dev_north) {
    cat(paste0("~~~\nApply planar fit method (", method, ") based on"), 
        if (is.character(avg_time)) {
            avg_time
        } else {
            paste(parse_time_diff(avg_time), "secs")
        }
        , "intervals\n"
    )
    pf_avg_time <- parse_time_diff(avg_time)
    # PF periods:
    x[, pf_period := floor(
        as.numeric(Time - start_time[1], units = "secs") / 
            pf_avg_time)]
    ### tag 'complete' intervals
    x[, pf_keep := .N / (pf_avg_time * Hz) > data_threshold, by = pf_period]
    # get Uavg & provis. WD
    x[, c("um", "vm", "wm", "Uavg", "WDprov") := {
        um <- mean(u)
        vm <- mean(v)
        wm <- mean(w)
        .(
            um, vm, wm
            , sqrt(um ^ 2 + vm ^ 2 + wm ^ 2)
            , csystem_wd(coord_system, um, vm)
        )
    }, by = pf_period]
    x[, WDprov := (WDprov + dev_north[1]) %% 360]
    # fit to plane
    wd_sort <- unique(sort(wd_sectors))
    if (wd_sort[1] != 0 && wd_sort[length(wd_sort)] != 360) {
        wd_sec <- c(wd_sort[length(wd_sort)], wd_sort)
    } else {
        wd_sec <- wd_sort
    }
    no_sec <- length(wd_sec) - 1
    if (is.null(u_thresh)) {
        u_thresh <- 0
    }
    x[, c("alpha", "beta", "w_bias") := NA_real_]
    cat(sprintf("\tU threshold = %1.2f m/s (%1.1f%% of data set)\n", u_thresh, 
            x[Uavg > u_thresh, .N] / nrow(x) * 100))
    if (no_sec >= 1) {
        cat("\tFollowing wind sectors have been defined:\n")
        for (i in seq.int(no_sec)) {
            x_sub <- x[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1] & Uavg > u_thresh &
                pf_keep]
            n_ints <- x_sub[, uniqueN(pf_period)]
            n_data <- x_sub[, .N]
            cat(paste0("\tSector ", i, ": (", wd_sec[i], "\u00B0 to ", 
                    wd_sec[i + 1], "\u00B0] (", n_ints, " intervals", 
                    sprintf(" %1.1f%% of wind speed subset)", 
                        n_data / x[Uavg > u_thresh, .N] * 100)
                    , "\n"))
            if (n_ints < n_thresh) {
                if (n_data == 0) {
                    cat("\t\t|-> No data in this sector.\n")
                } else {
                    cat(paste0("\t\t|-> Not enough intervals (", 
                            n_ints, 
                            ") for planar fit with threshold of ", 
                            n_thresh, 
                            " (applying two axis rotation for this sector)\n"))
                }
            } else {
                # planar fitting
                PF <- x_sub[, .(um = um[1], vm = vm[1], wm = wm[1]), 
                    by = pf_period][, {
                    list2env(pf_fit(um, vm, wm, FUN = reg_fun, 
                            method = method)[
                        c("Pmat", "cw", "alpha", "beta")])
                }]
                # transform:
                x[WDprov > wd_sec[i] & WDprov <= wd_sec[i + 1], 
                    c("u", "v", "w", "alpha", "beta", "w_bias") := {
                    c(pf_transf(u, v, w, with(PF, Pmat), with(PF, cw)), 
                        list(with(PF, alpha), with(PF, beta), with(PF, cw)))
                }]
            }
        }
        # browser()
    } else {
        x_sub <- x[Uavg > u_thresh, ]
        n_ints <- x_sub[, uniqueN(pf_period)]
        n_data <- x_sub[, .N]
        cat("No wind sectors have been defined.\n")
        if (n_ints < n_thresh) {
            cat(paste0("\t-> Not enough intervals (", n_ints, 
                    ") for planar fit with threshold of ", n_thresh, 
                    " (applying two axis rotation)\n"))
        } else {
            cat( n_ints, " intervals available\n")
            # planar fitting
            PF <- x[Uavg > u_thresh, .(um = um[1], vm = vm[1], wm = wm[1]), 
                by = pf_period][, {
                list2env(pf_fit(um, vm, wm, FUN = reg_fun, 
                        method = method)[
                    c("Pmat", "cw", "alpha", "beta")])
            }]
            # transform:
            x[, c("u", "v", "w", "alpha", "beta", "w_bias") := .(
                pf_transf(u, v, w, with(PF, Pmat), with(PF, cw)), with(PF, alpha), 
                with(PF, beta), with(PF, cw))]
        }
    }
    ### subtract w_avg (is done by detrending) and remove pf_keep column:
    x[, c('pf_period', 'pf_keep', 'WDprov', 'Uavg') := NULL]
    x
}
pf_transf <- function(u, v, w, P, cw) {
	list(
		u = P["up", "u"] * u + P["up", "v"] * v + P["up", "w"] * (w - cw)
		, v = P["vp", "u"] * u + P["vp", "v"] * v + P["vp", "w"] * (w - cw)
		, w = P["wp", "u"] * u + P["wp", "v"] * v + P["wp", "w"] * (w - cw)
		)
}
pf_fit <- function(u, v, w, FUN = MASS::rlm, method = c("Wilczak2001", "vanDik2004"), ...) {
	if (method[1] %in% "Wilczak2001") {
		# Wilczak 2001 equation 39, wm = b0 + b1*um + b2*vm:	
		mod <- FUN(w ~ u + v, ...)
		b <- coef(mod)
		names(b) <- c("b0", "b1", "b2")
		cw <- b["b0"]
		b1 <- b["b1"]
		b2 <- b["b2"]
	} else {
		# van Dik 2004, wm = b1*um + b2*vm:	
		mod <- FUN(w ~ u + v - 1, ...)
		b <- coef(mod)
		names(b) <- c("b1", "b2")
		cw <- 0
		b1 <- b["b1"]
		b2 <- b["b2"]		
	}
	# equation 42, p3j:
	p31 <- -b1 / sqrt(b1 ^ 2 + b2 ^ 2 + 1)
	p32 <- -b2 / sqrt(b1 ^ 2 + b2 ^ 2 + 1)
	p33 <- 1 / sqrt(b1 ^ 2 + b2 ^ 2 + 1)
	# equation 44, angles:
	# tan_b <- -p32/p33
	sin_b <- -p32 / sqrt(p32 ^ 2 + p33 ^ 2)
	cos_b <- p33 / sqrt(p32 ^ 2 + p33 ^ 2)
	sin_a <- p31
	cos_a <- sqrt(p32 ^ 2 + p33 ^ 2)
	# equation 2, D + C:
	Cmat <- cbind(c(1, 0, 0), c(0, cos_b, sin_b), c(0, -sin_b, cos_b))
	Dmat <- cbind(c(cos_a, 0, -sin_a), c(0, 1, 0), c(sin_a, 0, cos_a))
	# equation 36, P = t(D)%*%t(C):
	Pmat <- t(Dmat) %*% t(Cmat)
	colnames(Pmat) <- c("u", "v", "w")
	rownames(Pmat) <- c("up", "vp", "wp")
	list(
        Pmat = Pmat, 
        cw = cw, 
        mod = mod, 
        alpha = asin(sin_a) / pi * 180, 
        beta = asin(sin_b) / pi * 180
    )
}

# find the dynamic lag time (search for max in cov fun)
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
	# find max relative to baseline (average):
    avg <- mean(x[m + seq(-500, 500)])
	maxis <- ind[which.max(abs(x[ind] - avg))]
	c(index = maxis, tau = maxis - m)
}

# reduce cospecta "resolution" for plotting
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

# estimate high-frequency damping by fitting a good-quality
# ogive to a pre-defined frequency range
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

# calculate wind statistics and MOST parameters
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

# calculate average wind speed from MOST profile
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

# plot time series including filtered trend
plot.tseries <- function(dat,wind,scal,selection,color,units){
	msg <- paste(c(format(dat[1,1],"%d.%m.%Y")," - time series"),collapse="")
	tsbeginning <- dat[1,1]
    # fix wind variables
	dat[, c("u", "v", "w", "T")] <- mapply("+", 
        wind[c("umrot", "vmrot", "wmrot", "Tmdet")], 
        wind[c("uprot", "vprot", "wprot", "Tdet")])
	### get trends:
	dat3 <- list2DF(wind[c("umrot","vmrot","wmrot","Tmdet")])
	names(dat3) <- c("u","v","w","T")
	if(!is.null(scal)){
        # fix scalars
        dat[, names(scal)] <- lapply(scal, \(x) {
            if (is.null(isna <- na.action(x$residuals))) {
                x$fitted + x$residuals
            } else {
                out <- rep(NA_real_, nrow(dat))
                out[-isna] <- x$fitted + x$residuals
                out
            }
        })
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
	dat2 <- reshape2::melt(dat[,c("st",selection)],id="st")
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
			# panel.grid(h=-1, v=-1, lty=3, col="gray80")
			panel.xyplot(x,y,...)
			y2 <- dat2[list(...)$subscripts,"trend"]
            if (!identical(y, y2)) {
                # panel.xyplot(x,y2,type="l",lwd=1.5, lty=3, col="gray30")
                # panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightblue")
                # panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightgrey")
                # panel.xyplot(x, y2, type = "l", lwd = 2, lty = 2, col = "#B37FDF")
                panel.xyplot(x, y2, type = "l", lwd = 3, col = "darkgrey")
            }
		}
	)  	
}

# plot covariance function
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

# plot cospectra and ogives
plot_cospec_ogive <- function(ogive, cospec, freq, ylab = NULL, xlim = NULL, cx = 1.5, 
    col = "lightblue", nred = floor(sqrt(sqrt(length(ogive))) * 1.5), model_par = NULL) {
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
	lines(cospec_reduced0$freq,y_cf,col="#E5E0E0")
    # modelled cospec/ogive
    if (!is.null(model_par) && !anyNA(model_par)) {
        # cospec: f * Co(f)
        cs <- freq * cospec_model(model_par['fx'], model_par['m'], model_par['mu'], 
            model_par['A0'], freq)
	    y_cs <- (cs - min(prCo)) / diff(range(prCo)) * diff(ylim) + ylim[1]
        lines(freq, y_cs, col = '#AE71EB99', lwd = 2)
        # ogive
        og <- ogive_model(model_par['fx'], model_par['m'], model_par['mu'], 
            model_par['A0'], freq)
        lines(freq, og, col = '#AE71EB99', lwd = 2)
    }
    # measured cospec/ogive
	lines(cospec_reduced$freq,y_crm,type="b",col="black",lwd=2)
	lines(freq,ogive,col=col,lwd=2)
}

# convert colors to hex
col2hex <- function(name, alpha) {
    if (!grepl('^#', name)) {
        # check if alpha provided
        if (!missing(alpha) && is.character(alpha)) {
                alpha <- as.integer(as.hexmode(alpha)) / 255
        }
        m <- col2rgb(name) / 255
        rgb(m[1, ], m[2, ], m[3, ], alpha)
    } else {
        # check if alpha provided
        if (missing(alpha)) {
            name
        } else {
            if (is.numeric(alpha)) {
                alpha <- as.hexmode(as.integer(alpha * 255))
            }
            paste0(sub('^(#\\d{6}).*', '\\1', name), alpha)
        }
    }
}

# plot damping from fitted ogives
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
		, rotation_method = c("two axis", "planar fit")
		, rotation_args = list(
            phi = NULL, 
            pf_avg_time = avg_period, 
            pf_FUN = MASS::rlm, 
            pf_method = c("Wilczak2001", "vanDik2004"), 
            pf_N_thresh = 10, 
            pf_wd_sectors = c(0, 360), 
            pf_U_thresh = 0
        )
        # detrending -> valid entries are blockAVG,linear,linear_robust,ma_xx (xx = time in seconds)
        , detrending = c(u = 'blockAVG', v = 'blockAVG', w = 'blockAVG', T = 'blockAVG', 
            nh3_ppb = 'blockAVG', nh3_ugm3 = 'blockAVG', h2o_mmolm3 = 'blockAVG', 
            co2_mmolm3 = 'blockAVG')
        , hard_limits = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, 
            nh3_ppb = TRUE, nh3_ugm3 = TRUE, h2o_mmolm3 = TRUE, 
            co2_mmolm3 = TRUE)
        , hard_limits_lower = c(u = -30, v = -30, w = -10, T = 243, 
            nh3_ppb = -100, nh3_ugm3 = -100, h2o_mmolm3 = -100, 
            co2_mmolm3 = -100)
        , hard_limits_upper = c(u = 30, v = 30, w = 10, T = 333, 
            nh3_ppb = 5000, nh3_ugm3 = 5000, h2o_mmolm3 = 5000, 
            co2_mmolm3 = 5000)
        , hard_limits_window = '5mins'
        , hard_limits_replace = FALSE
		, covariances = c('uxw', 'wxT', 'wxnh3_ugm3', 'wxh2o_mmolm3', 'wxco2_mmolm3')
        # fix lag in seconds
		, lag_fix = c(uxw = 0, wxT = 0, wxnh3_ppb = -0.4, wxnh3_ugm3 = -0.4, 
            wxh2o_mmolm3 = -0.2, wxco2_mmolm3 = -0.2)
        # dyn lag in seconds around lag_fix
		, lag_dyn = c(uxw = 0.2, wxT = 0.2, wxnh3_ppb = 1.5, wxnh3_ugm3 = 1.5, 
            wxh2o_mmolm3 = 1.5, wxco2_mmolm3 = 1.5)
        # re_rmse window
        , gamma_time_window = c(2.5, 5)
        # damping_reference: either a specific covariance as 'wxT' or best quality ogives
        #   such as either 'base_quality' (for best base_quality_fix/_dyn)
        #   or 'ogive_quality' (for best ogive_quality_fix/_dyn)
		, damping_reference = c(wxnh3_ppb = 'wxT', wxnh3_ugm3 = 'wxT', 
            wxh2o_mmolm3 = 'wxT', wxco2_mmolm3 = 'wxT')
        # lower & upper bounds of fitting ogives (in seconds)
		, damping_lower = c(wxnh3_ppb = 2, wxnh3_ugm3 = 2, wxh2o_mmolm3 = 2, wxco2_mmolm3 = 2)
		, damping_upper = c(wxnh3_ppb = 20, wxnh3_ugm3 = 20, wxh2o_mmolm3 = 20, wxco2_mmolm3 = 20)
        , low_cont_sec = 20
        , high_cont_sec = 2
        , cont_pts = 5
        # , subintervals = TRUE
        # , subint_n = 5
        # , subint_detrending = c(u = 'linear', v = 'linear', w = 'linear', T = 'linear', nh3_ppb = 'linear', nh3_ugm3 = 'linear', h2o_mmolm3 = 'linear', co2_mmolm3 = 'linear')
        , oss_threshold = 0
        , na_alarm_code = c(1:3, 5:8, 11, 13)
        , thresh_period = 0.75
		, create_graphs = TRUE
		, graphs_directory = NULL
		, add_name = ''
        , plot_timeseries = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, 
            ht_oss = TRUE, nh3_ppb = FALSE, nh3_ugm3 = TRUE, 
            li_co2ss = TRUE, h2o_mmolm3 = TRUE, co2_mmolm3 = TRUE)
        , plotting_var_units = c(u = 'm/s', v = 'm/s', w = 'm/s', T = 'K', 
            ht_oss = '-', nh3_ppb = 'ppb', nh3_ugm3 = 'ug/m3', 
            li_co2ss = '-', h2o_mmolm3 = 'mmol/m3', co2_mmolm3 = 'mmol/m3')
        , plotting_var_colors = c(u = 'gray20', v = 'gray20', w = 'gray20', T = 'orange', 
            ht_oss = 'grey', nh3_ppb = 'indianred', nh3_ugm3 = 'indianred', 
            li_co2ss = 'grey', h2o_mmolm3 = '#8FC1E6', co2_mmolm3 = 'seagreen4')
        , plotting_covar_units = c(uxw = 'm2/s2', wxT = 'K*m/s', wxnh3_ppb = 'ppb*m/s', 
            wxnh3_ugm3 = 'ug/m2/s', wxh2o_mmolm3 = 'mmol/m2/s', wxco2_mmolm3 = 'mmol/m2/s')
        , plotting_covar_colors = c(uxw = 'gray70', wxT = 'orange', wxnh3_ppb = 'indianred', 
            wxnh3_ugm3 = 'indianred', wxh2o_mmolm3 = '#8FC1E6', wxco2_mmolm3 = 'seagreen4')
		, ogives_out = FALSE
        , as_ibts = TRUE
	){

    # start processing
	script_start <- Sys.time()

    # check input
    if (create_graphs) {
        # TODO: eventually remove one of either arguments
        if (is.null(graphs_directory)) {
            stop('argument "create_graphs" is TRUE, but "graphs_directory" is not specified!')
        }
        library(reshape2)
        library(lattice)
    }
    if (create_graphs && (
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
    rotation_method <- match.arg(rotation_method, c('two axis', 'planar fit'))
    rot_args <- rotation_args
    rotation_args <- eval(formals(process_ec_fluxes)$rotation_args)
    rotation_args[names(rot_args)] <- rot_args
    detrending <- fix_defaults(detrending, variables)
    hard_limits <- fix_defaults(hard_limits, variables)
    hard_limits_lower <- fix_defaults(hard_limits_lower, variables)
    hard_limits_upper <- fix_defaults(hard_limits_upper, variables)
    lag_fix <- fix_defaults(lag_fix, covariances)
    lag_dyn <- fix_defaults(lag_dyn, covariances)
    damping_reference <- fix_defaults(damping_reference, covariances[scalar_covariances])
    damping_lower <- fix_defaults(damping_lower, covariances[scalar_covariances])
    damping_upper <- fix_defaults(damping_upper, covariances[scalar_covariances])
    # subint_detrending <- fix_defaults(subint_detrending, variables)
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
    # ts_vars <- union(names(plot_timeseries), ts_vars)
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
    # } else {
        # check if ht available
        # if (!ht_provided) {
        #     stop('neither ht8700 nor licor data or directory has been provided -> cannot process fluxes without concentration data!')
        # }
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
    }
    if (any(end_time - start_time <= 0)) {
        stop('argument "end_time" must be strictly greater than "start_time"!')
    }
    if (length(start_time) == 1) {
        start_time <- seq(start_time, end_time, by = avg_secs)
        end_time <- start_time + avg_secs
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

    # prepare interval times
    dates_utc <- unique(date(c(start_time, end_time)))
    dates_formatted <- gsub('-', '_', dates_utc, fixed = TRUE)

    # read raw data:
    # ------------------------------------------------------------------------------

    # TODO: fix missing intervals in sonic/ht/licor!!!
    #       fix old sonic files format
    #       fix old sonic time format
    #       fix gill bug in old sonic
			# Data[, w := w * ifelse(w < 0, 1.289, 1.166)]

    # get sonic data
    if (sonic_has_data) {
        cat('~~~\nSubsetting sonic data - ')
        # subset provided sonic data
        sonic <- sonic_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
    } else {
        cat('~~~\nReading sonic files\n')
        # select files
        sonic_selected <- sonic_files[
            sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', sonic_files) %in% dates_formatted
            ]
        # read sonic files
        if (length(sonic_selected)) {
            sonic <- rbindlist(lapply(
                    file.path(sonic_directory, sonic_selected), 
                    \(x) {
                        cat('\t')
                        read_sonic(x)
                    }
            ))
        } else {
            sonic <- NULL
        }
        # no UTC conversion needed, since this case is excluded up to now
        # # convert to UTC
        # sonic[, Time := with_tz(Time, 'UTC')]
    }
    cat('done\n')
    # check sonic data
    if (is.null(sonic) || nrow(sonic) == 0) {
        stop('No sonic data available within given time range! Aborting routine!\n')
    }

    # get ht8700 data
    if (ht_provided && ht_has_data) {
        cat('Subsetting HT8700 data - ')
        # subset
        ht <- ht_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
        cat('done\n')
    } else if (ht_provided && !ht_with_sonic) {
        cat('Reading HT8700 files\n')
        # select files
        ht_selected <- ht_files[
            sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', ht_files) %in% dates_formatted
            ]
        # read ht files
        if (length(ht_selected)) {
            ht <- rbindlist(lapply(
                    file.path(ht_directory, ht_selected), 
                    \(x) {
                        cat('\t')
                        read_ht8700(x)
                    }
            ))
        } else {
            ht <- NULL
        }
        # no UTC conversion needed, since this case is excluded up to now
        # # convert to UTC
        # ht[, Time := with_tz(Time, 'UTC')]
        cat('done\n')
    } else {
        # no ht data
        ht <- NULL
    }
    # check ht
    if (ht_provided && is.null(ht)) {
        cat('-> No HT8700 data within given time range.\n')
    }

    # get licor data
    if (licor_provided && licor_has_data) {
        cat('Subsetting LI-7500 data - ')
        # copy
        licor <- licor_directory[Time >= start_time & Time < tail(end_time, 1), ]
        cat('done\n')
    } else if (licor_provided && !licor_with_sonic) {
        cat('Reading LI-7500 files\n')
        # select files
        licor_selected <- licor_files[
            sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', licor_files) %in% dates_formatted
            ]
        if (length(licor_selected)) {
            # read new licor files
            licor <- rbindlist(lapply(
                    file.path(licor_directory, licor_selected), 
                    \(x) {
                        cat('\tFile:', x, '- ')
                        out <- read_licor(x)
                        cat('done\n')
                        out
                    }
            ))
        } else {
            licor <- NULL
        }
        # no UTC conversion needed, since this case is excluded up to now
        # # convert to UTC
        # licor[, Time := with_tz(Time, 'UTC')]
    } else {
        # no licor data
        licor <- NULL
    }
    # check licor
    if (licor_provided && is.null(licor)) {
        cat('-> No LI-7500 data within given time range.\n')
    }

    cat('Merging files - ')
    daily_data <- merge_data(sonic, ht, licor)
    cat('done\n~~~\n')

    # check variables and covars and subset by columns
    if (!all(cn_check <- variables %in% names(daily_data))) {
        stop(
            "cannot find variable(s): ", 
            paste(variables[!cn_check], collapse = ", "),
            "\nAvailable column names are: ", 
            paste(names(daily_data), collapse = ", ")
        )
    }
        
    # define bins & subset again & just make sure sonic has no missing data
    daily_data <- daily_data[, bin := getIntervals(Time, start_time, end_time)][
        bin != 0L & is.finite(u) & is.finite(v) & is.finite(w) & is.finite(T) &
            !is.na(Time)]

    # get Hz/frequency etc.
    rec_Hz <- daily_data[, {
        d_t <- diff(as.numeric(Time))
        round(1 / median(d_t, na.rm = TRUE), -1)
    }]
    if (rec_Hz < 10) {
        stop('Frequency is lower than 10 Hz! Aborting script!')
    }
    n_period <- avg_secs * rec_Hz
    n_threshold <- thresh_period * n_period
    freq <- rec_Hz * seq(floor(n_period / 2)) / floor(n_period / 2)

	# remove incomplete intervals
    cat('Check given "thresh_period"\n')
    np0 <- daily_data[, uniqueN(bin)]
	daily_data <- daily_data[, keep := .N >= n_threshold, by = bin][(keep)][,
        keep := NULL]
    np <- daily_data[, uniqueN(bin)]
    if (np0 - np > 1) {
        cat(np0 - np, 'intervals have less than', thresh_period * 100, '% of valid data\n')
    } else if (np0 - np > 0) {
        cat(np0 - np, 'interval has less than', thresh_period * 100, '% of valid data\n')
    }

    # check again if empty
    if (nrow(daily_data) == 0) {
        stop('No valid data available within given time range!')
    }

    # --------------------- read fix lags from lag_lookuptable ---------------------
    # TODO: only if necessary/wanted

    # fix and dyn lags:                                      
    # ------------------------------------------------------------------------------ 
    input_fix_lag <- round(lag_fix * rec_Hz)
    input_dyn_lag <- rbind(
        lower = round((lag_fix - lag_dyn) * rec_Hz)
        , upper = round((lag_fix + lag_dyn) * rec_Hz)
    )

    # be verbose
    cat("\n~~~\nCalculation will include",
        daily_data[, uniqueN(bin)]
        , "of", length(start_time), "intervals between", 
        format(start_time[1], format = "%Y-%m-%d %H:%M"), "and", 
        format(tail(end_time, 1), 
            format = "%Y-%m-%d %H:%M", usetz = TRUE), 
        "on a", avg_secs / 60, "minute basis\n~~~\n\n"
    )

    # raw data quality control I, i.e. hard limits = physical range
    # --------------------------------------------------------------------------
    cat("~~~\nchecking hard limits...\n")
    if (any(hard_limits)) {
        daily_data[, variables[hard_limits] := check_limits(
            mget(variables[hard_limits]), Time, rec_Hz, 
            lim_range[, variables[hard_limits], drop = FALSE], 
            hard_limits_window, hl_method
            )]
    } 

    # define coordinate system
    if (daily_data[, sonic[1] == 'HS']) {
        coord.system <-  'HS-Wauwilermoos'
    } else {
        coord.system <-  'Windmaster'
    }

    # declination
    current_declination <- daily_data[, mag_dec(Time[1]), by = bin][, V1]
    d_north <- dev_north + current_declination

    ## coordinate rotation
	cat("~~~\nrotating data...\n")
    # apply planar fit
    if (rotation_method == 'planar fit') {
        cat('applying planar-fit rotation...\n')
        daily_data <- planar_fit(daily_data, coord_system = coord.system, 
            method = rotation_args[['pf_method']],
            avg_time = rotation_args[['pf_avg_time']],
            wd_sectors = rotation_args[['pf_wd_sectors']],
            u_thresh = rotation_args[['pf_U_thresh']],
            n_thresh = rotation_args[['pf_N_thresh']],
            reg_fun = rotation_args[['pf_FUN']],
            start_time = start_time[1], data_threshold = thresh_period, Hz = rec_Hz, 
            dev_north = d_north
        )
    } else {
        cat('applying two-axis rotation...\n')
    }

    # apply two-axis rotation (if planar fit has not been applied successfully)
	daily_data[, c("WD", "phi", "urot", "vrot", "wrot") := rotate_twoaxis(u, v, w,
        phi = if (rotation_method[1] %in% "two axis") {
            # two-axis rotation
            rotation_args$phi 
        } else if (is.na(alpha[1])) {
            # planar fit failed
            NULL 
        } else {
            # planar fit successful
            0
        }, c.system = coord.system)
		, by = bin]

	### correct for sonic north deviation
	daily_data[, WD := (WD + d_north[.GRP]) %% 360, by = bin]

    # prepare ogive output
    if (ogives_out) {
        e_ogive <- new.env()
        e_ogive$Cospec_dyn_Out <- e_ogive$Cospec_fix_Out <- e_ogive$Covars_Out <- 
            e_ogive$Ogive_fix_Out <- e_ogive$Ogive_dyn_Out <- 
                setNames(vector("list", length(start_time)), format(start_time))
    }


    # loop over bins
    results <- daily_data[, {

        # be verbose
        cat("\n\n~~~~~~~~\n", format(start_time[.BY[[1]]]), " - interval ", .GRP, 
            " of ", .NGRP, "\n", sep = '')

        # check number of data
        if (.N > n_threshold) {

            # ugly copy because of inability to change .SD values
            SD <- copy(.SD)

            # check HT8700 quality: alarm codes & OSS
            if (ht_provided) {
                cat('~~~\nChecking HT-8700 OSS and alarm codes - ')
                # add alarm code
                if (!('ht_alarm_code' %in% names(SD))) {
                    SD[, ht_alarm_code := get_alarms(.SD)]
                }
                # create regex pattern
                na_alarm_pattern <- paste(paste0('\\b', na_alarm_code, '\\b'),
                    collapse = '|')
                # check alarms and set nh3 NA
                nh3_vars <- grep('nh3', names(SD), value = TRUE)
                na0 <- SD[, sum(is.na(get(nh3_vars[1])))]
                SD[grepl(na_alarm_pattern, ht_alarm_code),
                    (nh3_vars) := NA_real_]
                na1 <- SD[, sum(is.na(get(nh3_vars[1])))]
                # check oss
                SD[ht_oss < oss_threshold,
                    (nh3_vars) := NA_real_]
                cat('done\nBad alarms:', na1 - na0, '\nValues below OSS threshold:', 
                    SD[, sum(ht_oss < oss_threshold, na.rm = TRUE)], '\n')
            }

            # TODO: check licor ss!

            # get start of interval
            Int_Start <- start_time[.BY[[1]]]
            # TODO (maybe for later): include NA where d_t>2*mean, remove entries where d_t < 0.5*mean

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
            cat("~~~\ndeterending sonic data...\n")
            wind <- SD[, {
                ud <- trend(urot, detrending["u"], rec_Hz)
                vd <- trend(vrot, detrending["v"], rec_Hz)
                wd <- trend(wrot, detrending["w"], rec_Hz)
                Td <- trend(T, detrending["T"], rec_Hz)
                I(list(
                    uprot = ud$residuals
                    , vprot = vd$residuals
                    , wprot = wd$residuals
                    , Tdet = Td$residuals
                    , umrot = ud$fitted
                    , vmrot = vd$fitted
                    , wmrot = wd$fitted
                    , Tmdet = Td$fitted
                ))
            }]
            SD[, c("u", "v", "w", "T") := wind[c("uprot", "vprot", "wprot", "Tdet")]]

            # calculate some turbulence parameters and collect some wind parameters
            # -------------------------------------------------------------------------- 
            wind_stats <- wind_statistics(wind, z_canopy, z_ec)

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
                    detrending[scalars], MoreArgs = list(Hz_ts = rec_Hz), SIMPLIFY = FALSE)
            # calculate scalar averages and sd:
            # -------------------------------------------------------------------------- 
                scalar_means[scalars] <- sapply(scalar_list, mean)
                scalar_sd[scalars] <- sapply(scalar_list, sd)
            # assign to SD
                SD[, (scalars) := lapply(names(detrended_scalars), \(nms) {
                    if (!is.null(isna <- na.action(scalar_list[[nms]]))) {
                        out <- rep(NA_real_, .N)
                        x <- detrended_scalars[[nms]]$residuals
                        out[-isna] <- x
                        out
                    } else {
                        detrended_scalars[[nms]]$residuals
                    }
                    })]
            # } else {
                # cat('No scalars available -> skipping interval!\n')
                # next
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
            lo_range <- m - rev(gamma_time_window) * 60 * rec_Hz
            hi_range <- m + gamma_time_window * 60 * rec_Hz
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
                    xs <- fft(data.table::shift(x, lag, type = 'cyclic')) / N
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
                    xs <- fft(data.table::shift(x, lag, type = 'cyclic')) / N
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


            # calculate low frequency contribution
            i_hi <- which(1 / freq < high_cont_sec)[1]
            if (length(i_hi) != 1) stop('check argument "high_cont_sec"!')
            hi_cont_fix <- sapply(Ogive_fix, \(x) {
                mean(x[i_hi + seq(-cont_pts, cont_pts)])/ x[1]
            })
            hi_cont_dyn <- sapply(Ogive_dyn, \(x) {
                mean(x[i_hi + seq(-cont_pts, cont_pts)])/ x[1]
            })
            i_lo <- which(1 / freq <= low_cont_sec)[1]
            if (length(i_lo) != 1) stop('check argument "low_cont_sec"!')
            lo_cont_fix <- sapply(Ogive_fix, \(x) {
                (x[1] - mean(x[i_lo + seq(-cont_pts, cont_pts)])) / x[1]
            })
            lo_cont_dyn <- sapply(Ogive_dyn, \(x) {
                (x[1] - mean(x[i_lo + seq(-cont_pts, cont_pts)])) / x[1]
            })

            # fit theoretical ogive shape
            cat('\t- estimate flux quality\n')
            ogff <- function(pars) {
                ogive_model(pars['fx'], 3 / 4, pars['mu'], pars['A0'], freq)
            }
            fit_og <- function(pars, ogive, ind = NULL) {
                if (is.null(ind)) {
                    fit_ogive(pars, ogive, freq, 1, length(ogive))
                } else {
                    fit_ogive(pars, ogive, freq, i_lo, i_hi)
                }
            }

            # base range
            i_r <- i_lo:i_hi
            # set limit for valid mu values
            mu_lim <- 10

            # quality of ogive 
            ogive_quality_fix <- sapply(Ogive_fix, \(x) {
                # ini <- c(fx = 0.05, m = 3 / 4, mu = 1 / 6, A0 = x[i_lo] / 50)
                ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                opt <- optim(ini, fit_og, ogive = x, control = list(maxit = 5e3))
                if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                    return(rep(NA_real_, 6))
                }
                fog <- ogff(opt$par)
                d <- (fog - x) / fog[1]
                d0 <- diff(d)
                c(1 - mean(d0 ^ 6) ^ (1 / 6) * 6, fog[1], opt$par, 3 / 4)
            })
            ogive_quality_dyn <- sapply(Ogive_dyn, \(x) {
                ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                opt <- optim(ini, fit_og, ogive = x, control = list(maxit = 5e3))
                if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                    return(rep(NA_real_, 6))
                }
                fog <- ogff(opt$par)
                d <- (fog - x) / fog[1]
                d0 <- diff(d)
                c(1 - mean(d0 ^ 6) ^ (1 / 6) * 6, fog[1], opt$par, 3 / 4)
            })
            # quality of base part of ogive 
            base_quality_fix <- sapply(Ogive_fix, \(x) {
                ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                opt <- optim(ini, fit_og, ogive = x, ind = i_r, 
                    control = list(maxit = 5e3))
                if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                    return(NA_real_)
                }
                # d <- (ogff(opt$par)[i_r] - x[i_r]) / mean(x[i_r])
                # 1 - mean(d ^ 2) ^ (1 / 3)
                fog <- ogff(opt$par)
                d <- (fog - x) / fog[1]
                d0 <- diff(d)
                1 - mean(d0 ^ 6) ^ (1 / 6) * 6
            })
            base_quality_dyn <- sapply(Ogive_dyn, \(x) {
                ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                opt <- optim(ini, fit_og, ogive = x, ind = i_r, 
                    control = list(maxit = 5e3))
                if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                    return(NA_real_)
                }
                # d <- (ogff(opt$par)[i_r] - x[i_r]) / mean(x[i_r])
                # 1 - mean(d ^ 2) ^ (1 / 3)
                fog <- ogff(opt$par)
                d <- (fog - x) / fog[1]
                d0 <- diff(d)
                1 - mean(d0 ^ 6) ^ (1 / 6) * 6
            })
            
            # get Albrecht's ogive bias
            ogive_bias_fix <- sapply(Cospec_fix, \(x) {
                x[is.na(x)] <- 0
                ms <- median(sign(x))
                o <- sum(x)
                (ms * sum(abs(x)) - o) / o
            })
            ogive_bias_dyn <- sapply(Cospec_dyn, \(x) {
                x[is.na(x)] <- 0
                ms <- median(sign(x))
                o <- sum(x)
                (ms * sum(abs(x)) - o) / o
            })

            # fix matrix to vector for output
            rownames(ogive_quality_fix) <- rownames(ogive_quality_dyn) <- 
                c('q', 'f', 'fx', 'mu', 'A0', 'm')
            ogive_par_fix <- ogive_quality_fix[3:6, ]
            ogive_par_dyn <- ogive_quality_dyn[3:6, ]
            ogive_fitted_fix <- ogive_quality_fix[2, ]
            ogive_quality_fix <- ogive_quality_fix[1, ]
            ogive_fitted_dyn <- ogive_quality_dyn[2, ]
            ogive_quality_dyn <- ogive_quality_dyn[1, ]

            # should ogives be provided with output
            if (ogives_out) {
                int_start <- format(Int_Start)
                e_ogive$Cospec_fix_Out[[int_start]] <- c(
                    list(freq = freq, model_coef = ogive_par_fix), Cospec_fix)
                e_ogive$Cospec_dyn_Out[[int_start]] <- c(
                    list(freq = freq, model_coef = ogive_par_dyn), Cospec_dyn)
                e_ogive$Ogive_fix_Out[[int_start]] <- c(
                    list(freq = freq, model_coef = ogive_par_fix), Ogive_fix)
                e_ogive$Ogive_dyn_Out[[int_start]] <- c(
                    list(freq = freq, model_coef = ogive_par_dyn), Ogive_dyn)
                e_ogive$Covars_Out[[int_start]] <- c(
                    list(Hz = rec_Hz), Covars)
            }

            # empirical damping estimation, dyn and fix should have best reference (dyn/fix)...
            # ------------------------------------------------------------------------
            if (any(scalar_covariances)) {
                cat("\t- damping\n")
                damping_reference_fix <- damping_reference
                damping_reference_dyn <- damping_reference
                # check qualities
                if (any(bq <- grepl('base_quality', damping_reference))) {
                    best_base_fix <- names(base_quality_fix)[which.max(base_quality_fix)]
                    damping_reference_fix[bq] <- best_base_fix
                    best_base_dyn <- names(base_quality_dyn)[which.max(base_quality_dyn)]
                    damping_reference_dyn[bq] <- best_base_dyn
                }
                if (any(oq <- grepl('ogive_quality', damping_reference))) {
                    best_ogive_fix <- names(ogive_quality_fix)[which.max(ogive_quality_fix)]
                    damping_reference_fix[oq] <- best_ogive_fix
                    best_ogive_dyn <- names(ogive_quality_dyn)[which.max(ogive_quality_dyn)]
                    damping_reference_dyn[oq] <- best_ogive_dyn
                }
                Damping_fix <- mapply(damp_hac5, ogive = Ogive_fix[scalar_covariances], 
                    ogive_ref = Ogive_fix[damping_reference_fix], freq.limits = damp_region,
                    MoreArgs = list(freq = freq), SIMPLIFY = FALSE)
                Damping_dyn <- mapply(damp_hac5, ogive = Ogive_dyn[scalar_covariances], 
                    ogive_ref = Ogive_dyn[damping_reference_dyn], freq.limits = damp_region,
                    MoreArgs = list(freq = freq), SIMPLIFY = FALSE)
            } else {
                Damping_dyn <- Damping_fix <- NULL
            }


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
                    , et = end_time[.BY[[1]]]
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
                    , if (length(input_scalars) > 0) {
                        as.list(c(
                            setNames(scalar_means[input_scalars], 
                                paste0('avg_', input_scalars))
                            , setNames(scalar_sd[input_scalars], 
                                paste0('sd_', input_scalars))
                        ))
                    }
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
                        fix_lag_out[input_covariances] / rec_Hz, 
                        paste0('lag_fix_', input_covariances)
                    )
                    , setNames(
                        dyn_lag_out[input_covariances] / rec_Hz, 
                        paste0('lag_dyn_', input_covariances)
                    )
                    , setNames(
                        re_rmse[input_covariances],
                        paste0('re_rmse_', input_covariances)
                    )
                    , setNames(
                        hi_cont_fix[input_covariances],
                        paste0('hi_cont_fix_', input_covariances)
                    )
                    , setNames(
                        hi_cont_dyn[input_covariances],
                        paste0('hi_cont_dyn_', input_covariances)
                    )
                    , setNames(
                        lo_cont_fix[input_covariances],
                        paste0('lo_cont_fix_', input_covariances)
                    )
                    , setNames(
                        lo_cont_dyn[input_covariances],
                        paste0('lo_cont_dyn_', input_covariances)
                    )
                    , setNames(
                        ogive_fitted_fix[input_covariances],
                        paste0('ogive_fitted_fix_', input_covariances)
                    )
                    , setNames(
                        ogive_fitted_dyn[input_covariances],
                        paste0('ogive_fitted_dyn_', input_covariances)
                    )
                    , setNames(
                        ogive_quality_fix[input_covariances],
                        paste0('ogive_quality_fix_', input_covariances)
                    )
                    , setNames(
                        ogive_quality_dyn[input_covariances],
                        paste0('ogive_quality_dyn_', input_covariances)
                    )
                    , setNames(
                        base_quality_fix[input_covariances],
                        paste0('base_quality_fix_', input_covariances)
                    )
                    , setNames(
                        base_quality_dyn[input_covariances],
                        paste0('base_quality_dyn_', input_covariances)
                    )
                    , setNames(
                        ogive_bias_fix[input_covariances],
                        paste0('ogive_bias_fix_', input_covariances)
                    )
                    , setNames(
                        ogive_bias_dyn[input_covariances],
                        paste0('ogive_bias_dyn_', input_covariances)
                    )
                    , if (!is.null(Damping_fix) && length(scalar_covariances_only) > 0) {
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
                    , if (!is.null(Damping_dyn) && length(scalar_covariances_only) > 0) {
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

            ##### figures
            # -----------------------------------------------------------------------

            if (create_graphs) {
                if (!dir.exists(path_folder)) {
                    dir.create(path_folder, recursive = FALSE)
                }
                # get end of current interval time in UTC
                soi_utc <- start_time[.BY[[1]]]
                eoi_utc <- end_time[.BY[[1]]]
                # get date in correct format
                date_formatted <- format(soi_utc, '%Y%m%d')
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
                jpeg(file = paste0(path_folder, '/', plotname, ".jpg"), width = 600, 
                    height = (sum(plot_timeseries)) * 100, quality = 60)
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
                    plotname <- paste("plots", date_formatted, time2, 
                        covariances_plotnames[i], sep = "-")
                    # fix ylab
                    ylab <- sub('(.+)x(.+)', "<\\1'\\2'>", i)
                    jpeg(file = paste0(path_folder, '/', plotname, ".jpg"), 
                        width = 1350, height = 900, quality = 60)
                        par(mfrow=c(2,3))
                        # ----------------------- Covariance -----------------------
                        plot_covfunc(Covars[[i]], n_period / rec_Hz, dyn_lag_max[, i], 
                            fix_lag[i], ylab = ylab, xlim = c(-50, 50), cx = 1.5, 
                            cxmt = 1.25, cl = plotting_covar_colors[i])
                        # ---------------------- Co-Spec/Ogive fix lag -----------------------
                        plot_cospec_ogive(Ogive_fix[[i]], Cospec_fix[[i]], freq, 
                            ylab = paste0("ogive (fix lag) of ", ylab), cx = 1.5, 
                            col = plotting_covar_colors[i], 
                            model_par = ogive_par_fix[, i]
                        )
                        # ---------------------- Co-Spec/Ogive dyn lag -----------------------
                        plot_cospec_ogive(Ogive_dyn[[i]], Cospec_dyn[[i]], freq, 
                            ylab = paste0("ogive (dyn lag) of ", ylab), cx = 1.5, 
                            col = plotting_covar_colors[i], 
                            model_par = ogive_par_dyn[, i]
                        )
                        if (scalar_covariances[i]) {
                            # ---------------------- empirical damping -----------------------
                            plot_damping(Damping_fix[[i]], freq, ylab = 
                                paste0("ogive (fix lag) of ", i), cx = 1.5, 
                            col = plotting_covar_colors[i])
                            plot_damping(Damping_dyn[[i]], freq, ylab = 
                                paste0("ogive (dyn lag) of ", i), cx = 1.5, 
                            col = plotting_covar_colors[i])
                        }
                        title(paste0(ylab, " flux ", 
                            format(soi_utc, format = "(%H:%M:%S"), " - ", 
                            format(eoi_utc, format = "%H:%M:%S)"), 
                            if (scalar_covariances[i]) {
                                if (damping_reference[i] == 'ogive_quality') {
                                    reflab <- paste0(damping_reference[i], ': ',
                                        sub('(.+)x(.+)', "<\\1'\\2'>", best_ogive_fix),
                                        if (best_ogive_fix != best_ogive_dyn) {
                                            paste0(
                                                ' (fix), ',
                                                sub('(.+)x(.+)', "<\\1'\\2'>", best_ogive_dyn),
                                                ' (dyn) '
                                            )
                                        })
                                } else if (damping_reference[i] == 'base_quality') {
                                    reflab <- paste0(damping_reference[i], ': ',
                                        sub('(.+)x(.+)', "<\\1'\\2'>", best_base_fix),
                                        if (best_base_fix != best_base_dyn) {
                                            paste0(
                                                ' (fix), ',
                                                sub('(.+)x(.+)', "<\\1'\\2'>", best_base_dyn),
                                                ' (dyn) '
                                            )
                                        })
                                } else {
                                    reflab <- sub('(.+)x(.+)', "<\\1'\\2'>", 
                                        damping_reference[i])
                                }
                                paste0(" - damping reference flux: ", reflab)
                            }
                        ), outer = TRUE, line = -1)
                    dev.off()	
                }
            } # end plotting
            cat('~~~\nfinished interval.\n~~~~~~~~')
            # return out
            out
        } else {
            cat("less than ", round(thresh_period * 100),"% of data points in raw data - skipping interval.\n", sep = '')
            # return NULL
            NULL
        }
    }, by = bin]

	cat("\n************************************************************\n") 
	cat("operation finished @", format(Sys.time(), "%d.%m.%Y %H:%M:%S"), 
        "time elapsed: ", difftime(Sys.time(), script_start, unit = "mins"),
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

    # output incl. ogives
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

## function to fit theoretical ogive shape to measurement
# library(Rcpp)
cppFunction('
double fit_ogive(const NumericVector paras, const NumericVector ogive, const NumericVector f,
    const int ilo, const int ihi) {
    const double len = ogive.size();
    const double m = 3.0 / 4.0;
    const double fx = paras[0];
    const double mu = paras[1];
    const double A0 = paras[2];
    const double A0_fx = A0 / fx;
    const double m_mu_pow = (m + 1.0) / (2.0 * mu * m);
    double last_value = 0.0;
    double ss = 0.0;
    // loop in reverse over ogive and get cumsum
    for (int i = len - 1; i >= (ilo - 1); i--) {
        // get cospec value devided by f
        last_value = last_value + 
            A0_fx / (
                std::pow(1.0 + m * std::pow(f[i] / fx, (2.0 * mu)), m_mu_pow)
            );
        // get difference to ogive
        if (i < ihi) {
            ss += std::fabs(ogive[i] - last_value) * std::sqrt(1 / f[i]);
        }
    }
    return ss;
}
')
# xx <- fit_ogive(ini, og, freq)
# yy <- fit_og(ini, og)
# fcpp <- function(x) fit_ogive(x, og, freq)

# convenience functions for theoretical cospec/ogive models
cospec_model <- function(fx, m, mu, A0, f = freq) {
    A0 / fx / (
        (1 + m * (f / fx) ^ (2 * mu)) ^ ((m + 1) / (2 * mu * m))
    )
}
ogive_model <- function(fx, m, mu, A0, f = freq) {
    rev(cumsum(rev(cospec_model(fx, m, mu, A0, f))))
}

