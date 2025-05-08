
library(data.table)
library(ibts)
library(Rcpp)
library(sodium)
library(fields)

## helper functions ----------------------------------------

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
        whflgs <- hflgs > 0
		# since we're only intrested in specific time windows, find first NAs:
		isna <- lapply(dat[whflgs],function(x)which(is.na(x)))
		# replace if all NA? (Why could this happen?)
		l <- length(dat[[1]])
		replAll <- (l - lengths(isna)) < 2
		if(any(replAll)) dat[whflgs][replAll] <- rep(-99999,l)
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
			# dat[whflgs][[i]][isna[[i]]] <- sapply(ind, function(x) mean(dat[whflgs][[i]][x], na.rm = TRUE))
			dat[whflgs][[i]][isna[[i]]] <- sapply(ind, function(x) median(dat[whflgs][[i]][x], na.rm = TRUE))
		}	
		cat("number of replaced values\n*~~~~*\n", names(dat[whflgs]),"\n", lengths(isna),"\n*~~~~*\n")
	} else {
		cat(
            paste0(
                'number of values outside of hard limits (set to NA)\n*~~~~*\n',
                paste(
                    sprintf('%s: %i', names(dat[whflgs]), hflgs[whflgs > 0]), 
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
                    'detrending method "', method,'" is not valid!\n',
                    'Valid methods are one of "blockAVG",',
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
find_dynlag <- function(x, dyn) {
	n <- length(x)
	if (n %% 2) {
		# ungerade
		m <- (n + 1) / 2
	} else {
		# gerade
		m <- n / 2 + 1
	}
    # get searching window
	ind <- seq(dyn[1], dyn[2]) + m
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
            dat3 <- cbind(dat3, dat[, res, drop = FALSE])
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
plot_covfunc <- function(cov_func, avg_t, dynLag, fixLag, ylab = NULL, xlim = NULL, 
    cx = 1.5, cxmt = 1.25, cl = "black", re = NULL) {
	midP <- dynLag[1] - dynLag[2]
	n <- length(cov_func)
	fix_cov <- cov_func[midP + fixLag]
	dyn_cov <- cov_func[dynLag[1]]
	Hz <- n / avg_t
	if (is.null(xlim)) {
		xlim <- c(-avg_t / 2, avg_t / 2)
		x <- seq(-avg_t / 2, avg_t / 2, length.out = n)
	} else {
		xlim2 <- round(xlim * Hz)
		x_hi <- intersect(seq(xlim2[1], xlim2[2]), seq(n) - midP)
		x <- x_hi / Hz
		cov_func <- cov_func[midP + x_hi]
	}
    ylim <- range(cov_func, na.rm = TRUE)
    if (add_re <- !is.null(re)) {
        re2 <- c(-1.96, 1.96) * re
        ylim <- range(ylim, re2, na.rm = TRUE)
    }
	plot(0, type = "n", xlim = xlim, ylim = ylim, ylab = ylab, cex.axis = cx, 
        xlab = expression(paste(italic(tau), " (s)")), cex.lab = cx)
	abline(h = 0, lwd = 1.5, lty = 2, col = "gray60")
    abline(v = 0, lwd = 1.5, lty = 2, col = "gray60")
	abline(v = fixLag / Hz, lty = 3, col = "lightgrey")
	abline(v = dynLag[2] / Hz, lty = 4, col = "lightgrey")
    # add random error indication
    if (!is.null(re)) {
        abline(h = re2, col = cl, lty = 3)
    }
    # draw cov lines
	lines(x, cov_func, lwd = 2, col = cl)
	if (fixLag == dynLag[2]) {
		mtext(substitute(paste(italic(tau) == a, "s: ", y), 
                list(y = sprintf("%1.8f", dyn_cov), a = sprintf("%1.1f", fixLag / Hz))), 
            side = 3, cex = cxmt)
	} else {
		mtext(substitute(
                paste(italic(tau)[fix] == a, "s: ", x, " / ", 
                    italic(tau)[dyn] == b, "s: ", y), 
                list(x = sprintf("%1.8f", fix_cov), 
                    y = sprintf("%1.8f", dyn_cov), 
                    a = sprintf("%1.1f", fixLag / Hz), 
                    b = sprintf("%1.1f", dynLag[2] / Hz)
                    )), 
            side = 3, cex = cxmt)
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
        # extend/replace x with default
        if (nms[1] == '') {
            # replace all
            x <- setNames(rep(x, length(missing_vars)), missing_vars)
        } else {
            # get formals
            frms <- formals(sys.function(sys.parent(1L)))
            # get formal name
            fname <- as.character(substitute(x))
            # get default values
            default <- eval(frms[[fname]])
            # extend with defaults
            x <- c(x, default[missing_vars])
        }
    }
    # return
    x[vars]
}

## main flux processing function ----------------------------------------

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
        , gamma_time_window = c(5, 10)
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
        , co2ss_threshold = 0
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
        licor <- licor_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
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
    cat('checking interval data threshold...\n')
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

        # get start of interval
        interval_start <- start_time[.BY[[1]]]
        # TODO (maybe for later): include NA where d_t>2*mean, remove entries where d_t < 0.5*mean

        # be verbose
        cat("\n\n~~~~~~~~\n", format(interval_start), " - interval ", .GRP, 
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
                SD[grepl(na_alarm_pattern, ht_alarm_code), (nh3_vars) := NA_real_]
                na1 <- SD[, sum(is.na(get(nh3_vars[1])))]
                # check oss
                SD[ht_oss < oss_threshold, (nh3_vars) := NA_real_]
                cat('done\nBad alarms:', na1 - na0, '\nValues below OSS threshold:', 
                    SD[, sum(ht_oss < oss_threshold, na.rm = TRUE)], '\n')
            }

            # check licor ss
            if (licor_provided) {
                cat('~~~\nChecking LI-7500 signal strength - ')
                # get variables
                li_vars <- grep('^(co2|h2o)_', names(SD), value = TRUE)
                # check li_co2ss
                SD[li_co2ss < co2ss_threshold, (li_vars) := NA_real_]
                cat('done\nValues below "co2ss" threshold:', 
                    SD[, sum(li_co2ss < co2ss_threshold, na.rm = TRUE)], '\n')
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

            # # calculate auto-covariances
            # # -------------------------------------------------------------------------- 
            # Power <- lapply(names(FFTs), \(i) {
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
            # names(Power) <- flux_variables

            # # power-spectra for flux variables
            # # ------------------------------------------------------------------------ 			
            # cat("\t- power-spectra\n")
            # PowerSpec <- lapply(flux_variables, function(i) {
            #         # check scalars
            #         N <- .N
            #         if (i %in% scalars) {
            #             if (!is.null(na_ind <- na.action(scalar_list[[i]]))) {
            #                 # fix lengths
            #                 FFTs[[i]] <- FFTs[[i]][-na_ind]
            #                 # get N
            #                 N <- length(FFTs[[i]])
            #             }
            #         }
            #         re <- Re(Conj(FFTs[[i]]) * FFTs[[i]])[seq(N / 2) + 1] * N / 
            #             (N - 1) * 2
            #         # get missing
            #         n_missing <- n_period / 2 - length(re)
            #         if (sign(n_missing) >= 0) {
            #             c(re, rep(NA_real_, n_missing))
            #         } else {
            #             re[seq_len(n_period / 2)]
            #         }
            #     })
            # names(PowerSpec) <- names(Power)

            # calculate covariances with fix lag time:
            # -------------------------------------------------------------------------- 
            cat("\t- covariances\n")
            Covars <- SD[, I(lapply(covariances_variables, function(i) {
                # TODO -> get ffts here and add to output
                # get x
                x <- get(i[1])
                # get y
                y <- get(i[2])
                # get NA values
                isfinite <- is.finite(x) & is.finite(y)
                # get N
                N <- sum(isfinite)
                # get ffts
                xfft <- fft(x[isfinite] / N)
                yfft <- fft(y[isfinite] / N)
                # get Re
                re <- Re(fft(Conj(yfft) * xfft, inverse = TRUE))
                # subset
                if (N %% 2) {
                    out <- re[c(((N + 1) / 2 + 1):N, 1:((N + 1) / 2))] * N / (N - 1)
                } else {
                    out <- re[c((N / 2 + 1):N, 1:(N / 2))] * N / (N - 1)
                }
                # get missing
                n_missing <- n_period - N
                if (sign(n_missing) >= 0) {
                    if (n_missing %% 2) {
                        n1 <- rep(NA_real_, (n_missing + 1) / 2)
                        n2 <- rep(NA_real_, (n_missing - 1) / 2)
                    } else {
                        n1 <- n2 <- rep(NA_real_, n_missing / 2)
                    }
                    out <- c(n1, out, n2)
                } else {
                    if (-n_missing %% 2) {
                        ind <- ((1 - n_missing) / 2 + 1):(N - (-n_missing - 1) / 2)
                    } else {
                        ind <- (1 - n_missing / 2):(N + n_missing / 2)
                    }
                    out <- out[ind]
                }
                # attach ffts
                structure(out, ffts = list(x = xfft, y = yfft), isfinite = isfinite)
            }))]
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
            m <- ifelse(n_period %% 2, (n_period + 1) / 2, n_period / 2 + 1)
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

            # cospectra for fixed & dynamic lags
            # ------------------------------------------------------------------------ 			
            cat("\t- co-spectra\n")
            Cospec_fix <- mapply(function(var) {
                    # get covar
                    covars <- Covars[[paste(var, collapse = 'x')]]
                    # get ffts
                    xfft <- attr(covars, 'ffts')$x
                    # get length
                    N <- length(xfft)
                    # get y
                    yfft <- attr(covars, 'ffts')$y
                    # get cospec
                    re <- Re(Conj(yfft) * xfft)[seq(N / 2) + 1] * N / (N - 1) * 2
                    # get missing
                    n_missing <- n_period / 2 - length(re)
                    if (sign(n_missing) >= 0) {
                        c(re, rep(0, n_missing))
                    } else {
                        re[seq_len(n_period / 2)]
                    }
                }, var = covariances_variables, SIMPLIFY = FALSE)
            names(Cospec_fix) <- covariances
            Cospec_dyn <- mapply(function(var, lag) {
                    # get covar
                    covars <- Covars[[paste(var, collapse = 'x')]]
                    # get ffts
                    xfft <- attr(covars, 'ffts')$x
                    # get length
                    N <- length(xfft)
                    # get y
                    if (lag != 0) {
                        y <- SD[attr(covars, 'isfinite'), get(var[2])]
                        yfft <- fft(data.table::shift(y, lag, type = 'cyclic')) / N
                    } else {
                        yfft <- attr(covars, 'ffts')$y
                    }
                    # get cospec
                    re <- Re(Conj(yfft) * xfft)[seq(N / 2) + 1] * N / (N - 1) * 2
                    # get missing
                    n_missing <- n_period / 2 - length(re)
                    if (sign(n_missing) >= 0) {
                        c(re, rep(0, n_missing))
                    } else {
                        re[seq_len(n_period / 2)]
                    }
                }, var = covariances_variables, lag = dyn_lag_max[2, ], SIMPLIFY = FALSE)
            names(Cospec_dyn) <- covariances

            # ogives for fixed & dynamic lags 
            # ------------------------------------------------------------------------ 
            Ogive_fix <- lapply(Cospec_fix, function(x) {
                rev(cumsum(rev(x)))
            })
            names(Ogive_fix) <- covariances
            Ogive_dyn <- lapply(Cospec_dyn, function(x) {
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
                int_start <- format(interval_start)
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
                    st = interval_start
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
                soi_user <- with_tz(interval_start, tz_user)
                eoi_user <- with_tz(end_time[.BY[[1]]], tz_user)
                # get date in correct format
                date_formatted <- format(soi_user, '%Y%m%d')
                # plotting:
                # -------------------------------------------------------------------------- 
                # -------------------------------------------------------------------------- 
                cat("~~~\nplotting timeseries and fluxes\n")
                # time series:
                # -------------------------------------------------------------------------- 
                # plot and save (rotated) data time series with raw-data trends...
                # ------------------------------------------------------------------------
                ## TODO: -> fix tz!!! -> use UTC but indicate in name!!!
                time2 <- format(soi_user, format = "%H%M")
                plotname <- paste("timeseries", date_formatted, time2, sep="-") 
                ts_vars <- names(plot_timeseries)[plot_timeseries]
                jpeg(file = paste0(path_folder, '/', plotname, ".jpg"), width = 600, 
                    height = (sum(plot_timeseries)) * 100, quality = 60)
                    ts_plot <- plot.tseries(
                        cbind(st = Time, as.data.frame(SD)),
                        wind, detrended_scalars, ts_vars,
                        plotting_var_colors, plotting_var_units)
                    # fix time zone
                    attr(ts_plot$x.limits, 'tzone') <- tz_user
                    print(ts_plot)
                dev.off()

                # plot and save flux evaluation...
                # ------------------------------------------------------------------------
                for(i in covariances){
                    # i <- "w'TDL CH4'"
                    # i <- covariances[3]
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
                            cxmt = 1.25, cl = plotting_covar_colors[i], re = re_rmse[i])
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
                            format(soi_user, format = "(%H:%M:%S"), " - ", 
                            format(eoi_user, format = "%H:%M:%S)"), 
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

## wpl post-processing ----------------------------------------

# wpl function
# TODO: add units to output & improve output var names (self-explanatory)
wpl_correction <- function(ec_dat, fluxes = c('nh3_ugm3', 'co2_mmolm3'), 
    dynamic_lag = TRUE, return_extra_cols = FALSE, water = 'h2o_mmolm3', 
    temperature = 'li_temp_amb', pressure = 'li_press_amb') {
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
    # define mu (due to paper formulae
    mu <- M_dry / M_h2o 
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
        # calculate water vapor mole fraction
        out[, xv_avg := get(paste0('avg_', water)) * 8.314 * t_k / p_pa / 1000]
        # calculate Pe in kPa
        out[, alpha_v := 2.48]
        out[, pe_kPa := p_pa / 1000 * (1 + alpha_v * xv_avg)]
        # get kappa values
        out[, c('kappa_p', 'dkappa_p_at_t', 'dkappa_t_at_p') := 
            as.data.table(t(mapply(get_kappa, tval = t_k - 273.14, pval = p_pa / 1000)))]
        out[, c('kappa_pe', 'dkappa_pe_at_t', 'dkappa_t_at_pe') := 
            as.data.table(t(mapply(get_kappa, tval = t_k - 273.14, pval = pe_kPa)))]
        # calculate coefficients A, B & C
        out[, paste0('coef_', letters[1:3]) := {
            coef_a <- kappa_p
            coef_b <- 1 + (1 - (alpha_v + 1) * xv_avg) * alpha_v * (pe_kPa * 1000) * 
                dkappa_pe_at_t / coef_a
            coef_c <- 1 + (1 - xv_avg) * out[, t_k] * 
                dkappa_t_at_pe / coef_a + xv_avg * (coef_b - 1)
            .(coef_a, coef_b, coef_c)
        }]
        f_nh3 <- 1e-6
        # density
        out[, rho_nh3 := avg_nh3_ugm3 * f_nh3]
        # raw flux
        out[, flux_nh3_raw := get(paste0('flux_', dynfix, 'wxnh3_ugm3')) * f_nh3]
        # corrected flux
        out[, flux_nh3_gm2s := coef_a * (
            # first term
            flux_nh3_raw +
            # second term
            coef_b * rho_nh3 / rho_dry * flux_h2o +
            # third term
            coef_c * (1 + mu * rho_h2o / rho_dry) * 
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
        out[, flux_co2_gm2s :=
            # first term
            flux_co2_raw +
            # second term
            rho_co2 / rho_dry * flux_h2o +
            # third term
            (1 + mu * rho_h2o / rho_dry) * 
                rho_co2 / t_k * cov_wT
            ]
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
    # add dyn/fix to names
    new_cols <- setdiff(names(out), names(ec_dat))
    flux_cols <- grep('flux', new_cols, value = TRUE)
    setnames(out, new_cols, sub('flux_', paste0('flux_', dynfix), new_cols))
    # switch back to ibts
    if (is_ibts) {
        out <- as.ibts(out)
    }
    out
}

# get kappa value
get_kappa <- function(tval, pval) {
    tout <- tval + (-1:1) * 0.01
    pout <- pval + (-1:1) * 0.05
    t_deg <- attr(get_kappa, 't_deg')
    p_kPa <- attr(get_kappa, 'p_kPa')
    k_values <- attr(get_kappa, 'k_values')
    out <- fields::interp.surface.grid(list(x = t_deg, y = p_kPa, z = k_values),
        list(x = tout, y = pout))
    c(
        kappa = out$z[2, 2],
        dkappa_P_at_T = (out$z[2, 3] - out$z[2, 1]) / (out$y[3] - out$y[1]) / 1000,
        dkappa_T_at_P = (out$z[3, 2] - out$z[1, 2]) / (out$x[3] - out$x[1])
    )
}

## add kappa values for HT ----------------------------------------

### get data to script
dat2string <- function(x) {
    key <- keygen(as.raw(rep(1, 32)))
    nonce <- as.raw(rep(1, 24))
    xs <- qs2::qd_serialize(x)
    xr <- data_encrypt(xs, key, nonce = nonce)
    as.character(xr)
}
string2dat <- function(x) {
    key <- keygen(as.raw(rep(1, 32)))
    nonce <- as.raw(rep(1, 24))
    xr <- as.raw(paste0('0x', x))
    qs2::qd_deserialize(data_decrypt(xr, key, nonce = nonce))
}
pasteN <- function(x, n) {
    y <- rbind(x, rep(c(rep('","', n - 1), '",\n"'), length(x))[1:length(x)])
    paste0('c("', paste(y[-length(y)], collapse = ''), '")')
}
add_to_eof <- function(x, attr_name, header, f_name, ec_file) {
    con <- file(ec_file, open = 'a')
    xs <- dat2string(x)
    on.exit({close(con)})
    writeLines(paste0('# ', header), con)
    writeLines(paste0(
            'attr(', f_name, ', "', attr_name, '") <- string2dat(\n',
            pasteN(xs, 18), ')\n'), con)
}

if (FALSE) {
    # add kappa values to function
    library(xlsx)
    # k_tab <- read.xlsx('~/var/coef-table-ht.xls', 2)
    # k_values <- as.matrix(k_tab[, 2:6])
    # p_kPa <- as.numeric(sub('X(\\d+)(kPa)?', '\\1', names(k_tab)[2:6]))
    k_tab <- read.xlsx('~/var/kappa_table_20210621.xls', 1)
    k_values <- as.matrix(k_tab[, -1])
    p_kPa <- as.numeric(sub('X(\\d+)(kPa)?', '\\1', names(k_tab)[-1]))
    t_deg <- as.numeric(k_tab[, 1])
    # add values to script
    ec_file <- '~/repos/5_GitHub/gel-scripts/ec-evaluation.r'
    add_to_eof(k_values, 'k_values', '~~~ kappa table ~~~', 'get_kappa', ec_file)
    add_to_eof(p_kPa, 'p_kPa', '~~~ pressure (columns) ~~~', 'get_kappa', ec_file)
    add_to_eof(t_deg, 't_deg', '~~~ temperature (rows) ~~~', 'get_kappa', ec_file)
}

# ~~~ kappa table ~~~
attr(get_kappa, "k_values") <- string2dat(
c("ce","84","5b","28","a3","a3","1b","e5","bb","b4","66","8c","8f","14","f1","f1","f9","f6",
"5c","9c","ab","ef","39","07","dd","7f","3e","39","6a","07","fd","6f","de","be","09","07",
"72","ad","e2","f2","74","b4","ff","b1","55","3c","09","16","02","ea","b1","f4","2b","75",
"ea","69","a1","1f","0d","c5","90","cc","d5","21","fd","3a","39","27","0b","1a","b0","3a",
"6a","11","34","b2","53","e1","3b","a3","0b","e8","7f","11","96","0b","47","6a","a7","8d",
"7b","7c","87","ba","17","0e","83","f1","51","8a","1a","dd","7c","f0","f6","26","d9","67",
"53","09","53","44","63","2c","74","fa","32","e4","16","24","e1","88","c0","ee","90","7f",
"b7","0e","2e","89","ca","8e","38","11","2f","ce","78","64","ae","9f","95","b9","b7","23",
"a5","6b","64","04","b0","7c","1a","c8","6c","24","6c","11","df","1c","45","9b","0e","d6",
"cb","b9","73","1b","cf","43","cb","0c","83","63","a7","a2","72","6f","7f","5e","48","8b",
"b3","24","87","d3","5b","13","b9","ef","f8","fe","87","73","96","c0","93","44","17","7d",
"c6","82","f8","46","e7","44","54","f7","66","f6","a8","44","db","94","e2","2d","0a","12",
"f7","72","91","68","99","ab","88","69","bb","77","ae","a3","d0","0c","73","bb","c8","16",
"74","3c","0e","5a","8f","21","56","18","17","89","c7","8d","e5","72","5b","f2","bf","7a",
"c8","c0","79","81","cc","c5","48","a3","17","5a","a0","44","1f","9e","c8","5a","00","7c",
"54","e4","86","65","3f","ca","25","62","16","f8","95","c8","ea","f8","3d","9d","c1","7d",
"e1","7e","99","e8","53","44","44","d7","57","5a","46","26","a1","c9","6c","df","e5","01",
"33","49","a3","69","c8","50","7d","4b","fe","1d","ad","77","33","09","e9","0c","0b","b7",
"d6","3a","c8","da","d0","14","2f","17","56","35","fc","1b","3c","6b","b4","67","d4","5f",
"44","b6","89","f6","54","c3","e4","d8","ad","49","9e","39","0f","88","1c","11","4e","2c",
"5a","ee","82","44","97","93","bd","89","b3","fc","eb","e4","55","bc","43","31","18","b0",
"7f","bd","39","4c","d9","a4","87","80","78","41","ab","bf","39","4a","33","07","78","24",
"a1","33","8d","3b","64","90","4b","20","4d","be","56","0f","37","1b","b7","dd","09","22",
"f6","9f","71","d4","37","9b","01","92","3b","44","73","bc","91","11","8b","f7","ae","7b",
"f5","46","ed","79","a8","2e","98","bf","47","88","47","57","3f","6a","43","12","31","d3",
"73","70","45","b3","87","93","4e","de","1b","07","2a","ea","d9","a9","4a","37","28","5d",
"f9","96","c9","64","24","f4","e8","e3","a9","da","d7","f6","47","b7","8d","d1","e1","8b",
"2c","af","45","f7","68","8b","44","9b","6e","00","d5","a7","a2","d7","46","4f","70","47",
"19","38","06","84","73","d9","20","fa","26","45","44","0e","cf","7c","fe","e8","a5","4b",
"34","b1","79","90","0c","89","50","7b","35","2b","21","5a","18","6c","23","15","73","03",
"fe","e3","3e","01","cd","61","4e","19","a0","65","0c","54","14","36","21","47","83","7e",
"60","ab","dc","91","83","f1","f2","1c","0c","8e","fc","72","42","34","96","c3","0b","22",
"c5","05","92","96","10","82","c1","6a","46","a6","88","25","fd","29","a4","44","43","74",
"ae","03","f5","96","35","2e","e1","e6","0a","7a","89","c0","9c","ad","8f","38","47","0b",
"52","66","58","c8","c9","6f","40","da","91","54","25","4f","de","86","9b","9e","bf","63",
"5f","9e","ae","6f","4b","67","5c","53","5f","70","13","0d","3a","9a","8f","7a","a0","2e",
"98","82","81","f2","64","15","0d","d2","0c","db","ef","67","44","f7","74","10","49","e1",
"00","93","6f","3d","b0","37","18","be","64","d7","3a","59","fe","a3","1f","f0","c6","6c",
"95","9f","dd","fb","3b","ea","8f","46","c8","38","f0","3c","b4","34","87","f2","46","35",
"29","cc","25","99","d7","73","6f","eb","69","dd","63","cb","4b","75","be","f5","f3","62",
"94","b0","4a","ea","27","d3","6a","7f","34","5a","17","0b","d7","27","02","b7","4f","39",
"f9","e5","05","64","c7","3d","f4","72","7d","c9","e0","f0","2d","b4","ff","70","24","28",
"7a","db","a6","9f","1f","89","ed","23","7d","06","30","30","f4","74","8c","d4","c2","83",
"8e","57","fe","4e","61","f7","cd","9b","60","fd","27","be","4f","42","91","b9","04","cc",
"32","80","25","71","51","5a","c0","f5","6c","64","59","31","f1","ea","b4","c9","b7","08",
"91","7d","bb","74","9c","47","5a","0c","3f","61","8d","07","41","b3","ad","83","bb","46",
"d7","47","f8","53","84","94","73","08","3f","8f","01","10","93","00","49","ea","33","62",
"9d","0e","42","8d","b8","60","dd","24","15","08","d0","5f","23","e3","28","bf","a6","57",
"88","82","f7","2f","53","9b","6c","0e","9f","13","11","e5","5b","51","4c","e3","a6","b6",
"11","8f","c6","aa","2d","8c","96","8c","82","fe","1f","9e","94","cb","00","3a","83","fe",
"e8","7f","33","07","b1","1b","43","76","22","b0","ad","69","ea","00","6b","91","ec","b4",
"9a","3a","dd","d9","9d","bd","1d","bb","f4","a0","6c","5d","24","cd","c6","09","60","18",
"eb","7d","b3","ef","16","dc","f9","a0","05","70","3b","ba","b6","1a","5a","82","9e","6e",
"c6","57","75","24","d2","c9","a3","39","e4","ba","b6","04","19","3e","af","42","f9","19",
"0f","4e","0d","0b","d4","bb","92","f6","52","6d","18","10","13","3d","fb","63","69","6f",
"8e","3a","13","68","b2","a7","23","54","0a","5c","86","a2","c3","08","6e","e4","14","89",
"9d","cf","02","35","ff","ea","25","98","da","00","af","49","79","56","46","dc","a8","06",
"83","00","1a","02","8d","fa","f5","35","fa","91","5f","b0","23","39","23","02","ae","d0",
"a6","6f","4a","8d","02","7a","b3","e8","23","01","cc","d8","95","e3","58","9e","b4","d1",
"fd","91","9f","8e","38","01","eb","74","a1","d7","f1","5c","a5","c6","6c","ac","d6","38",
"c9","d3","3e","63","63","bc","96","79","8c","b8","cd","32","48","3f","e3","46","20","ba",
"f0","74","3d","67","f3","70","0d","f6","d2","70","1c","a0","16","19","8a","fa","10","1b",
"da","cd","07","34","14","d7","19","ca","42","d3","3c","b8","10","1a","f8","22","42","aa",
"b4","87","a1","12","7c","c3","36","ed","ed","a8","86","2d","e3","cf","fd","05","de","c4",
"6b","64","77","fa","9c","3a","f3","50","ad","02","86","69","eb","0f","3d","a4","a4","96",
"81","21","86","85","c5","5f","df","c2","72","4f","a9","78","7c","9b","50","45","d6","26",
"d1","ad","30","47","63","59","9d","78","e7","83","83","5c","2a","8e","6b","e2","ca","d0",
"22","42","56","dd","9b","8a","c8","00","6c","77","b6","f5","e4","08","9b","62","e0","73",
"5e","f6","07","00","14","d2","a8","6d","0b","59","34","3b","de","05","68","bd","4e","d5",
"67","9e","03","8a","1f","c1","ca","ce","26","54","de","2c","91","27","6a","6b","c1","4a",
"f3","5d","f3","e7","35","e2","da","e9","2d","78","b4","d1","a5","5f","99","fc","53","33",
"f3","cb","e9","10","09","4c","db","9d","e0","96","6f","0f","78","a5","5c","af","15","32",
"62","e8","7f","d4","b9","26","09","aa","0c","65","aa","c2","4e","c6","0e","2c","f3","e7",
"7d","5c","d8","df","97","ed","7d","25","69","81","45","cf","a0","6a","d9","46","65","f3",
"18","1f","2d","ac","ff","06","80","5b","f2","16","ff","29","da","64","76","4c","87","51",
"02","90","f8","11","dd","81","fa","b2","f5","92","f8","90","d4","95","97","68","b9","b4",
"80","15","65","85","1a","0e","e6","e0","dc","36","7c","6d","18","6c","3b","6e","08","3a",
"70","b1","43","03","b6","44","be","eb","80","ea","70","5f","b5","69","2c","92","c8","eb",
"f8","5d","0b","be","19","a6","58","2c","df","32","ab","4b","0e","2d","bb","01","34","c1",
"73","c3","e3","44","2a","db","bf","a3","86","3c","fe","07","ce","14","bb","e0","54","27",
"0c","3e","17","c6","56","f1","97","9a","08","e8","e4","c9","6b","f2","8e","6a","89","2d",
"2b","4f","fb","2c","2c","a0","98","7c","b8","7a","c9","d9","3e","db","65","8c","9b","db",
"95","65","87","56","ca","98","0c","55","c6","8e","43","1d","44","df","6f","fa","e4","4a",
"91","a2","76","37","3f","6d","17","70","c4","03","d8","24","6b","e1","82","c2","f6","97",
"32","02","37","fa","ac","2b","77","13","9f","94","ee","41","21","db","4a","62","dd","47",
"b1","84","94","30","2b","d8","5c","cd","00","db","20","d9","05","6c","a9","8b","3b","2a",
"a9","20","0f","b0","ea","3b","f9","04","4f","11","07","b2","69","b9","51","42","9f","99",
"32","cb","1f","a6","5d","1d","ab","c4","be","8c","20","61","cc","3d","69","31","00","d0",
"72","aa","6a","c9","8f","c2","38","8f","b1","ff","f0","13","79","d5","13","a0","ba","23",
"ba","ff","63","42","d0","1f","0b","af","6e","9e","77","28","f9","06","f7","dd","ce","6f",
"31","44","8e","72","a7","85","e2","0e","9d","59","41","dd","43","bb","d5","5d","c9","62",
"fd","79","64","e4","61","0c","c2","42","ca","66","08","e0","21","4f","80","85","f6","01",
"4d","b4","e5","06","48","26","97","17","8e","e1","c9","57","e5","86","b3","c1","b2","e7",
"00","47","cf","92","19","cd","ba","25","50","2b","e2","8c","c3","0f","b2","20","2d","99",
"c4","56","fb","16","b3","d3","e7","c2","cd","cf","a4","44","f7","03","92","4e","4e","7d",
"fa","65","3c","6d","5f","ff","cd","49","69","67","8a","a5","34","80","57","f5","eb","24",
"ca","fc","a1","e6","33","4b","e0","17","16","35","c7","cc","54","e5","0a","91","ce","7d",
"cf","fb","5f","cf","f4","7c","b5","41","3f","52","b7","6a","08","e8","be","5a","bd","e8",
"1d","f1","6a","73","cb","9c","df","7a","dd","f4","88","75","6f","76","2f","f3","c0","61",
"85","da","80","f6","c0","98","ca","fe","6c","2c","27","a4","6e","c9","87","0e","39","ee",
"e4","f4","59","aa","22","83","6d","1c","1a","52","44","84","7b","59","88","57","27","03",
"97","1b","a5","af","0c","f4","e0","37","fb","da","5f","a4","b6","42","3a","93","d2","c2",
"18","ba","2d","78","0f","60","65","f6","b5","81","27","d0","8d","94","3e","1d","e6","38",
"6d","71","78","36","c7","c3","37","54","86","88","a5","91","b6","ba","80","2d","71","9f",
"ec","82","1d","15","78","f4","57","66","4b","7f","73","67","c1","59","4d","1e","3d","13",
"1a","e4","db","ce","fd","75","be","cb","9b","b9","53","d0","c7","a7","48","5a","37","29",
"59","25","8c","32","09","de","6d","33","5a","8d","78","63","f4","54","0b","ee","c3","9e",
"ce","71","26","ec","2e","cc","0a","c2","d9","b9","3a","2f","1e","3c","5a","21","bc","07",
"f2","2d","6c","99","17","42","11","c1","09","ca","f7","ca","83","ff","4c","c7","7c","b9",
"9e","69","39","9f","37","92","f3","7f","3c","69","4d","86","1c","6c","95","51","84","1e",
"ef","95","ea","72","f4","cd","d5","43","66","82","0c","9f","a4","d2","01","66","c7","75",
"42","d8","a0","2a","92","78","7a","41","8f","81","74","a7","d8","f6","ad","fd","c4","7a",
"08","5c","2b","f7","dd","fc","8b","f3","a4","80","14","9a","b5","ae","4c","71","29","cc",
"4b","16","89","e7","a2","20","8c","60","7f","37","6d","11","db","74","9c","01","90","2f",
"9c","f0","e4","31","4b","21","8a","4a","32","0d","b1","cf","e8","00","8a","33","47","c4",
"a1","5a","cb","4b","2f","f6","ad","f8","b4","c1","b0","69","f8","e7","87","ce","dc","4b",
"6f","6d","35","be","8b","f2","02","a5","29","d5","81","b1","7c","25","62","33","dc","e6",
"8e","62","f2","96","52","63","da","cc","70","98","87","8b","f8","92","ea","4d","e0","0e",
"22","a1","4f","3f","4a","50","6b","83","74","64","c6","c0","1b","27","91","17","fb","b3",
"1a","22","36","f6","e8","a1","1d","36","d2","1a","82","67","01","ce","36","77","4a","0f",
"3f","a7","f4","af","bc","b7","5d","f2","cb","c1","c5","75","09","11","4d","33","66","10",
"d3","2d","f6","7a","09","6d","84","2a","56","d8","49","8c","5f","c9","db","ca","b9","44",
"c4","ed","2c","65","4a","8b","dc","a0","b6","97","ef","37","b0","89","d0","99","20","57",
"61","93","48","b9","c6","35","84","69","26","76","4a","ed","75","1c","68","e1","86","89",
"9a","6a","11","2d","b2","e9","e5","12","35","57","c6","15","b2","89","42","c9","e3","09",
"4a","87","b5","8c","14","13","ba","54","d6","5b","ef","8d","af","a7","26","34","5d","73",
"53","27","13","c7","2b","7e","82","e9","dd","dc","7f","f7","3b","e9","7f","ec","5e","eb",
"4e","0b","9c","16","8e","6e","78","6d","91","2e","23","4d","3b","42","a7","6b","20","86",
"a7","da","ac","39","f7","8f","37","0f","14","79","e0","e1","6f","c3","83","b9","52","2b",
"09","d6","05","cf","76","95","a6","b7","58","48","d6","29","8d","58","95","66","63","8f",
"27","39","00","68","52","cd","45","d7","b8","06","1a","fa","5f","98","4a","7d","05","d3",
"10","15","74","2b","b1","07","ab","7f","98","fa","1d","53","18","67","de","a8","f9","f4",
"70","c0","99","ca","c9","99","47","71","2c","8a","10","b2","41","83","59","32","67","9c",
"3c","2f","7c","c0","19","a3","0c","87","07","a3","71","e7","50","fa","69","42","fd","61",
"d4","62","79","5b","b7","2d","1b","ff","89","ca","32","94","ef","98","dc","67","55","36",
"46","c2","65","4c","f9","7b","9f","9b","88","53","99","f3","8e","6f","30","f7","13","8a",
"22","04","ad","1a","3d","5a","0c","04","21","f6","b6","e5","ef","73","00","58","4f","35",
"e6","90","1a","8b","24","f8","7d","a7","0d","76","74","ad","b1","94","6d","19","1d","fb",
"a4","15","7c","fd","87","b3","ab","8e","af","b2","d1","27","2e","e9","95","35","19","b3",
"3c","42","f0","ad","a0","28","2a","27","dd","16","5e","7b","d0","3f","cc","bf","46","c5",
"16","10","62","47","da","9b","3b","06","b5","20","a6","00","e8","45","1d","68","4d","01",
"9a","85","9d","5a","e0","e0","43","1b","de","d2","5c","fc","db","a5","bb","12","ef","de",
"12","b7","1b","e3","bb","8a","25","5d","43","9a","4a","d8","2b","d7","4f","de","ce","55",
"67","1d","50","a7","b0","63","d9","5e","bc","01","09","48","8d","42","5a","c0","6d","6f",
"ab","e9","a7","09","38","72","98","11","7b","d4","c1","8e","6b","2f","f2","4b","9f","de",
"59","dd","b2","f5","c8","dd","df","88","18","d8","a0","ee","68","d8","4e","14","3e","6d",
"b4","87","f5","81","07","f7","d5","a8","61","81","61","88","80","f2","53","f5","95","cb",
"68","10","e3","cd","55","f5","e2","e5","d8","7b","d4","02","47","f0","54","c8","9e","08",
"3c","a4","ef","5f","85","2e","fe","60","70","ec","51","e1","d8","2d","70","7e","42","b4",
"eb","6c","6f","c4","47","85","21","51","ef","ad","85","ae","28","e7","73","83","b0","ad",
"86","fa","a3","ee","7e","e0","c1","fb","0e","59","50","83","39","e9","26","46","40","07",
"ce","98","04","9a","13","b0","d8","bb","f7","9c","7e","0c","2d","e5","b2","e7","0a","59",
"c3","e4","36","8c","b7","f4","51","53","57","63","b9","16","aa","73","31","fc","09","33",
"bb","50","0d","fa","86","83","d3","20","c8","b4","46","f7","29","29","89","0c","24","b9",
"ea","cd","3b","8e","1d","81","24","b1","df","82","bd","54","33","8e","f4","ab","a6","72",
"84","d3","76","cc","97","d8","bf","94","27","89","ca","0c","82","f7","85","e5","d3","ec",
"f0","23","85","f8","4e","41","09","8b","d8","25","e1","1d","70","a4","b4","74","17","20",
"52","65","6c","ab","50","06","e8","00","a4","ee","b9","47","8f","ac","04","8d","ec","c1",
"dc","54","44","58","60","24","a2","68","ad","74","c1","b2","ef","73","ca","15","f1","86",
"de","18","1f","7a","e4","13","51","84","8f","db","91","51","76","29","b1","09","23","72",
"e5","39","a6","d2","e1","f1","e1","4b","16","df","bb","71","03","06","47","e1","49","12",
"71","21","e3","8b","90","ec","a9","01","24","b8","5b","5d","bb","78","b0","ce","7c","48",
"43","91","7e","06","f0","9c","88","b9","9f","d3","75","15","32","e0","d0","ab","95","53",
"7f","69","1b","c2","b5","4c","d8","a4","d7","bb","05","11","9b","bf","2b","b5","88","4e",
"e8","3b","f1","82","65","bf","19","73","b3","f1","80","60","68","00","cb","b1","90","39",
"64","1a","6a","b2","a7","7c","8a","11","0a","17","46","7a","b8","ca","8d","cb","46","93",
"61","ab","3c","9e","7d","9c","01","f1","9a","d9","d8","bc","c2","a2","5a","5c","03","9e",
"01","62","0f","c5","f1","b9","7a","ef","9e","35","01","e5","d3","47","9b","6a","56","88",
"b0","d2","5a","49","04","c2","36","b0","f7","b0","35","8c","8c","62","7f","f7","7c","61",
"3d","34","e6","70","42","07","02","55","fb","92","c6","b1","f7","55","8f","66","d8","02",
"87","4d","80","af","98","d6","44","40","83","9e","45","7a","70","5f","f8","c7","fd","fa",
"c9","5b","f3","10","56","fd","31","9b","34","ac","0f","41","4f","b9","36","99","e7","39",
"08","c2","9c","77","43","cd","3f","03","67","23","90","84","70","a6","02","ac","42","84",
"94","3f","a9","d2","6a","68","b2","b2","3e","23","28","b0","e4","13","2d","b6","81","b2",
"97","57","ec","bd","49","a7","cd","de","d9","39","7c","2e","5e","1c","65","36","f2","e0",
"ab","c7","d5","93","16","23","0f","ef","ec","2d","47","c9","75","82","50","cd","ac","1f",
"3d","8f","91","40","c1","42","a5","5e","87","20","f1","ad","24","52","fe","b7","48","70",
"c4","3f","f9","cb","39","fa","02","c3","f5","b9","ad","44","4d","58","07","b3","b8","05",
"f8","a6","c2","3f","4b","cd","24","b3","9b","4b","a9","8f","5d","ca","5c","f7","9a","a9",
"9b","84","8b","1b","82","ab","11","5a","1b","41","27","6e","3a","93","d1","c9","a3","2b",
"9f","73","72","89","f0","aa","39","80","c7","50","53","6b","55","73","e0","ab","c6","7d",
"ed","9b","38","db","ad","4f","8d","84","56","c0","a6","98","40","75","1b","1a","a9","09",
"af","cf","8e","eb","75","01","55","0e","67","c3","48","98","3a","44","92","15","23","f1",
"cc","d0","eb","b3","29","14","2c","ed","f1","f7","90","af","31","8a","b0","67","65","1f",
"a0","49","48","6d","fc","38","5c","cd","b2","3b","c5","b1","e1","dc","ee","50","0f","10",
"40","bc","89","6f","38","e3","f0","19","31","7c","2a","25","c9","6a","82","40","95","6f",
"ec","01","82","5c","25","76","d4","e2","ee","8d","1f","d8","60","5f","e5","fc","e8","9e",
"3b","8f","7c","a9","ed","26","54","ae","4b","e4","fe","b6","0d","46","2d","22","23","72",
"e2","4f","4f","64","e7","9e","88","e1","99","d7","8a","0d","0e","a3","e1","ae","3b","e2",
"39","57","53","a9","04","b5","2e","b8","3b","1a","f5","1b","7c","fd","c4","1f","e1","0b",
"10","6b","c7","b3","30","5c","04","9e","8b","96","ce","16","1c","56","1c","3c","78","12",
"88","36","05","ca","2e","10","a3","27","c3","71","6e","43","ae","8b","9b","8d","da","31",
"9d","95","03","95","60","2a","dc","f0","31","ce","b7","d8","dd","47","eb","1c","44","a7",
"3f","cf","ef","eb","e1","ad","e8","15","a3","eb","a5","20","e6","9e","a4","90","63","3d",
"f4","db","c5","e1","49","48","73","d0","be","a9","1d","96","de","3f","d5","0e","68","86",
"16","f9","a8","e9","af","88","62","05","e7","26","6f","b9","b8","a1","f3","63","d1","0c",
"4c","8a","0b","62","7b","60","3d","78","8e","44","b5","04","62","ee","35","5b","3e","5f",
"67","67","bd","19","2c","9c","d3","3b","45","d8","9e","62","a6","a9","75","d1","53","e2",
"ec","6b","6c","6d","ec","cf","20","82","ef","a9","14","88","89","ee","3e","32","52","cf",
"b1","7d","d3","cb","52","95","28","b4","75","89","3d","c6","89","dc","ae","40","b1","cb",
"95","dd","30","ca","51","bd","f8","c6","a9","3d","a7","4c","34","24","b7","7d","8e","4e",
"e7","79","68","1a","81","3e","06","b5","86","5b","5a","62","c3","ba","4d","02","f7","55",
"89","4f","2d","3c","e9","15","eb","8e","c5","fa","25","68","b7","6d","e1","ad","8e","f7",
"26","bd","b3","b9","7c","4a","03","41","f0","ea","35","b9","8b","a0","3c","b3","49","f8",
"b2","f0","0b","20","8a","12","2b","cd","57","11","f3","67","b3","3f","1d","9e","8c","33",
"cd","03","38","7a","44","eb","af","85","a3","ae","fd","3b","a9","96","f2","07","b6","6c",
"5a","67","12","ae","22","2a","81","2f","92","ef","88","d9","77","e8","32","d4","79","0b",
"23","20","6f","e1","24","e5","da","82","a9","b6","75","10","86","77","d2","8b","f4","9f",
"b2","9d","95","4b","34","0e","e8","d1","fb","76","cf","89","9f","7e","ae","f9","36","ae",
"5d","45","dc","93","2a","e3","0f","71","48","ee","55","47","63","19","19","eb","2d","3e",
"29","aa","46","f0","92","d2","4d","8f","88","9f","13","3b","2c","bf","3d","be","73","3a",
"c1","d0","de","09","d9","9f","a3","b0","93","c4","50","cf","90","5e","7a","c7","5e","d5",
"20","7d","bc","2e","7a","65","3c","e2","c4","d2","75","04","48","34","6e","e4","69","73",
"23","4f","f1","a9","22","64","76","ee","b5","da","f1","c2","bc","5f","db","20","39","12",
"72","68","40","97","0a","7d","fb","13","f6","39","6c","a0","4c","93","0e","ff","ee","93",
"73","3b","0c","41","25","be","5d","ad","4a","ce","0c","bd","36","5a","26","f3","e8","d6",
"8a","41","e2","53","11","01","bc","16","c1","6b","88","57","c9","99","b5","51","b2","ec",
"49","aa","ed","19","43","41","25","a0","91","5e","75","87","72","07","ac","fc","f5","e7",
"7d","2c","32","ca","f3","c2","ae","8f","65","92","93","23","ed","f2","fb","4f","69","4f",
"80","95","e1","1b","5e","a9","8b","f1","5c","1a","f0","80","3a","02","16","98","86","06",
"b4","7f","5b","ed","23","56","9a","63","97","75","4c","4f","f4","d7","02","87","00","58",
"ea","5a","66","c9","fd","41","35","ab","07","82","78","19","25","93","87","81","2a","52",
"b5","f0","82","34","7b","8d","50","22","30","54","0e","b0","4d","c1","7b","f1","04","96",
"ef","60","ee","1d","e1","ff","00","40","97","85","15","6d","9f","ab","a7","cd","de","38",
"7d","d2","9a","9d","f9","79","16","0a","ff","57","4d","25","3e","00","87","cd","09","56",
"62","09","fd","81","79","36","c9","fd","8c","a6","65","10","f6","b2","66","dd","fc","96",
"b1","46","2f","0a","7c","d0","64","81","aa","13","b7","f1","db","df","28","78","a1","69",
"cf","49","ad","33","35","20","c9","9b","0e","01","96","50","fe","af","4e","f0","66","1e",
"fa","a6","dc","7c","ba","0e","0f","48","c1","91","12","7f","0c","f9","3b","20","2f","b5",
"9d","ab","1f","c6","0b","57","87","72","ea","1c","fe","7a","9e","31","4d","31","8d","2b",
"0e","b7","dd","b5","98","6e","ba","62","e4","6a","d2","64","8a","e0","d2","36","c9","8b",
"35","c7","aa","21","f4","c0","f8","4b","d4","22","7d","21","7f","33","f2","ab","af","97",
"7e","97","ca","a0","ee","df","41","26","2f","e5","3f","4a","c7","f9","62","ac","2e","f5",
"d8","b9","9b","95","43","e2","94","63","01","9d","94","d4","36","22","9a","1d","78","92",
"56","58","bd","ad","13","f4","75","72","dd","a1","cb","ff","96","99","47","5b","88","83",
"c3","7b","61","1a","3a","d8","b3","ef","78","0c","e2","c5","61","23","90","4c","e3","86",
"8e","56","c7","e4","51","29","5a","eb","3f","6f","4e","a6","f9","ed","64","99","54","cc",
"00","c3","85","58","1d","6c","1b","72","f8","3b","ee","77","ef","76","de","49","fd","51",
"77","1a","c5","9a","27","16","13","6a","24","14","09","b8","15","51","20","24","38","f2",
"6a","37","64","38","b6","b5","9d","fb","ae","00","04","1c","b8","48","00","cd","12","12",
"de","97","40","b0","f9","0e","a5","82","ec","d3","1b","a1","24","df","78","fd","bc","dc",
"11","59","9b","e6","39","20","ea","55","a3","06","5f","c5","0f","69","58","95","fd","c2",
"c0","51","11","bb","eb","37","30","60","d2","87","4a","1b","b4","52","f6","7d","5d","75",
"a5","43","07","c0","60","25","97","cd","4d","62","24","8e","99","75","33","cc","1a","3d",
"5d","e4","43","c9","f4","eb","d0","bc","2c","3b","74","c3","fa","3f","d7","d9","5d","ea",
"6a","da","d4","40","5a","5d","9e","0a","4e","0f","e7","bd","70","cf","a2","9c","56","f8",
"9c","4c","c5","8e","7e","da","a2","5b","b0","c2","7e","17","ed","8c","ae","40","fa","d6",
"80","bd","bc","b0","e8","f8","d0","c9","8e","2b","46","64","18","19","ee","6d","2b","22",
"c0","ff","7d","70","4e","23","b7","fc","31","b8","de","2a","b5","0b","07","fe","c9","b1",
"84","9b","4f","c4","44","a3","ee","5e","56","8b","bb","6a","14","e5","68","cf","bc","d2",
"3b","5d","cc","ec","53","6f","4d","cb","6e","e9","1f","eb","61","0e","55","ee","09","8d",
"ed","b6","36","f8","39","66","55","70","b9","db","0d","30","96","b0","6f","ba","d9","b7",
"c7","b6","b2","fc","74","d8","50","c9","15","ed","4f","29","c7","3a","2f","a5","e5","72",
"31","3e","4d","c4","b1","47","b6","d0","68","f4","83","81","1c","7b","43","da","c3","47",
"44","6a","69","5e","b2","81","ef","5b","65","e8","00","6d","a9","e3","5f","3b","83","59",
"9a","23","6a","ce","ea","d2","75","19","45","9e","c3","31","69","36","68","6e","6e","7c",
"b2","d8","02","da","2d","fa","f6","94","c1","59","00","46","63","f0","f3","e9","d3","9c",
"8d","50","67","4e","46","d8","57","47","30","9c","69","d2","48","7e","e9","d3","53","77",
"66","81","2e","20","ca","8b","b4","df","96","e0","19","20","76","9b","40","9b","a5","e2",
"a5","2e","c8","50","9c","19","f6","b7","cf","93","11","9a","88","59","2d","39","29","30",
"75","6e","38","e4","d4","7e","d9","37","e2","f2","56","30","1d","77","d7","f0","07","ff",
"d1","9d","d5","5e","4b","27","4c","5a","09","98","d8","94","9d","a5","02","48","2b","67",
"de","ea","5e","09","0e","23","c1","1e","cd","d6","18","cf","0e","3b","b8","22","1b","ba",
"28","fc","1b","65","c2","81","33","d4","22","80","7f","99","1f","a4","a5","b6","21","8d",
"55","29","58","13","7f","b8","a6","25","f7","d1","51","4d","68","77","d7","ae","da","71",
"50","d9","b1","7f","8c","11","df","98","c7","4f","c3","7b","9f","71","29","32","ff","62",
"97","78","d7","31","ec","30","3b","f6","8a","d3","91","40","90","da","75","69","3b","74",
"35","97","12","20","d9","6d","73","a5","64","4c","7b","09","1f","8b","08","f2","02","b2",
"ac","f8","93","63","43","dd","eb","91","dd","3a","79","29","71","c0","66","1f","b2","a0",
"a4","d1","65","47","34","2b","d1","bd","37","e4","38","c7","c2","9b","29","bf","eb","97",
"07","87","ae","37","1f","d9","48","81","61","fe","7e","6e","b0","45","f7","c7","5c","98",
"21","13","25","59","1f","e5","4b","07","28","20","06","3d","b2","4b","4b","d9","c9","d2",
"20","61","71","0b","9c","9c","98","f1","f6","fe","d0","05","de","55","c1","35","c6","96",
"70","00","4c","0f","fa","96","1e","33","90","29","a1","03","eb","ea","6c","2c","00","43",
"c7","82","57","a6","c3","52","dc","6b","2d","23","84","59","72","36","be","be","c7","c7",
"c0","97","f4","05","4c","a3","45","96","a9","0a","be","98","ac","3a","0e","bd","af","86",
"73","b4","ab","dd","a1","5d","9c","fe","77","b3","22","1f","71","eb","03","4c","e1","7d",
"eb","91","4b","1e","34","d7","07","c7","e7","3e","61","9c","9b","6d","e1","af","cf","47",
"71","ec","3a","e0","ee","82","56","08","a8","56","63","b0","c7","f7","d1","e2","51","6c",
"00","85","5c","45","f4","bb","21","9b","b0","45","4b","d0","b4","53","32","20","ea","e8",
"54","3a","96","b6","ae","f2","3c","30","02","a7","d1","4c","3a","07","ae","32","06","ed",
"36","b6","a8","be","b1","aa","56","a5","f1","3d","f3","67","06","56","5c","4d","6b","18",
"bd","d4","0d","81","c9","39","29","af","d7","f0","81","0b","1b","09","5c","f4","45","ad",
"8f","89","4c","29","b8","78","95","d2","c8","b2","b5","91","b6","bc","86","45","c0","7c",
"dc","b5","7c","70","15","8f","20","7c","6e","73","51","98","56","bc","ce","e7","ab","0f",
"20","37","4a","9a","79","37","89","7c","85","38","3f","67","c6","aa","fb","9c","be","4b",
"43","84","57","28","41","19","e6","d8","d0","2b","c3","5f","3c","21","f8","c8","a2","f1",
"90","76","54","fc","cd","f1","bc","bd","9e","4c","cd","f8","50","4d","90","4e","5a","52",
"a8","b7","3a","40","2c","07","12","c4","4c","58","08","99","75","f9","d3","22","6e","71",
"85","73","e5","21","8a","3e","87","51","a3","0e","6c","aa","6c","1a","5d","17","f3","ae",
"50","06","b0","c2","50","18","2f","0b","1b","00","c9","47","52","73","c6","0d","95","4f",
"8d","4c","5d","0b","d8","c1","ff","6d","2b","2e","f4","16","2c","e0","2c","45","5f","ca",
"66","ae","5b","d4","49","ed","f9","22","78","88","26","6d","36","f5","b3","3f","3d","d2",
"84","56","9e","64","58","0a","b4","b0","e5","dd","62","55","9e","17","32","d7","fc","89",
"df","2c","f4","28","aa","8f","07","4f","56","8b","c3","0d","0d","c7","7d","41","c9","cb",
"5c","ef","76","36","6b","51","bb","e8","13","f6","14","07","fc","72","61","f0","f6","08",
"04","26","43","63","0f","67","e7","91","8e","bf","53","f1","93","33","5e","f0","21","71",
"f2","a6","af","1c","65","9f","d8","99","8f","26","82","f5","39","5d","d5","6e","3c","16",
"e9","27","89","09","5f","d3","97","96","9e","a5","20","fc","d4","e4","fe","3e","93","37",
"34","ad","32","e2","5a","9b","60","5c","f2","64","b3","8d","56","4c","23","9a","34","b6",
"2c","df","e9","ab","f2","ae","bd","57","88","90","56","7a","27","7f","c8","20","3e","82",
"57","45","30","1d","67","9c","fe","a7","08","7c","28","5b","d1","2c","83","36","08","ed",
"d0","ba","6b","aa","aa","b0","b2","9b","fd","76","03","5b","e4","da","38","ab","29","ae",
"bb","b6","99","d1","48","c6","69","e5","c5","a9","8b","f8","55","37","57","de","c6","d8",
"a9","2a","d3","c3","41","d3","20","83","9e","93","67","f0","8b","97","d3","02","e5","71",
"27","62","4c","ad","87","20","97","d3","53","3e","2f","49","ab","ec","eb","a2","44","60",
"24","bf","86","a2","38","31","19","8f","49","9d","91","13","c4","d4","cf","e3","46","2c",
"22","4b","41","ed","35","59","fa","33","48","87","f1","6e","22","4e","e1","94","f6","f0",
"22","13","6a","cd","06","9b","f6","eb","5c","bb","14","85","df","ff","20","5c","e7","42",
"40","a2","74","33","2c","f5","33","82","d7","3a","1c","6e","2e","bb","6a","db","ea","00",
"44","97","fc","23","3f","97","26","5c","91","36","13","66","25","51","c1","8b","6d","17",
"7f","06","9b","7e","a4","62","ab","5b","ba","24","29","15","99","56","ad","ba","95","ae",
"09","9d","c1","32","43","52","69","0a","07","dd","8f","d5","e5","6b","69","46","3f","2c",
"52","ce","15","38","9e","9f","80","df","ac","a0","35","74","3c","28","15","1d","4e","56",
"c8","ac","56","41","17","80","bd","38","75","98","3e","e6","5f","3f","4f","c7","46","75",
"30","d2","2b","70","bd","b3","f5","73","d2","f1","a4","22","db","a6","16","fb","ee","57",
"3c","9c","94","3d","ce","66","c6","f5","2e","1d","cd","31","fc","ac","a3","88","d4","60",
"34","b1","b4","5f","33","a9","f8","26","e9","7c","b6","73","02","14","e4","90","2d","4f",
"98","f1","04","03","a0","2f","ad","03","98","63","ed","ef","28","72","4f","70","be","b6",
"ba","9c","fd","08","07","3a","14","a6","f8","eb","08","e1","f9","05","90","71","c2","84",
"e0","3a","84","25","c6","af","7a","f5","db","60","da","63","87","22","f0","64","d4","ad",
"5b","f2","a1","55","9e","97","d5","36","72","8e","f8","9f","57","d7","09","09","ca","78",
"38","57","27","d7","fb","f6","a6","5a","ee","7c","54","81","b8","38","ae","95","27","a8",
"73","49","a3","60","74","30","b0","9f","39","ae","df","0a","87","5d","8c","ae","b5","9e",
"42","1a","d9","27","b3","b9","35","07","92","12","fe","8c","78","89","33","f3","4f","bc",
"71","fc","b0","9d","70","20","26","88","19","d4","35","69","1a","2c","41","40","36","62",
"d3","1a","f3","3d","57","33","6e","9c","9b","d5","30","fe","67","eb","4d","1b","e7","1b",
"e2","08","cf","d5","3b","58","fc","c2","eb","5e","ec","af","9f","b1","50","21","c3","82",
"c4","96","2d","cb","42","fa","10","04","2e","50","6f","e7","e1","39","fe","a7","f1","47",
"7d","13","7c","b2","5a","d7","ff","f7","4f","10","56","2f","90","69","98","6d","ea","18",
"77","e1","eb","b9","ec","e5","ca","f3","1d","f2","fb","ef","07","46","9f","cc","5c","9b",
"d4","df","c5","bb","29","90","a7","55","ce","89","fa","10","5f","34","74","5c","05","42",
"ad","c4","8c","76","2e","68","89","9e","41","25","f5","cd","79","75","4e","f4","dc","b0",
"71","7c","74","8b","85","2c","16","5e","13","33","94","ca","c8","93","f0","b5","a5","6c",
"88","76","52","d8","c2","73","0b","f6","67","a9","92","71","d7","d8","84","f4","21","f5",
"6b","24","6a","41","d8","96","e6","ce","40","f0","86","82","f1","30","27","e3","12","34",
"2b","87","71","be","0f","77","9b","01","63","ec","37","db","5f","52","53","6f","e3","f7",
"0c","d2","25","e3","98","d5","93","1a","f9","3a","d9","ed","12","21","98","c6","81","7d",
"dc","c4","15","72","18","19","82","e7","19","7c","2f","2b","ee","1c","58","5d","d3","0f",
"87","4c","af","78","f4","ec","b6","30","f7","c9","71","4a","61","47","6c","0e","07","02",
"0c","73","71","d6","52","13","c3","7b","7a","29","00","f2","54","96","7c","47","47","68",
"3b","52","6a","80","0f","96","bb","92","ca","db","f4","95","eb","ca","c6","b1","ef","f8",
"19","15","63","bc","3f","2d","f9","ba","be","c8","b4","20","87","de","f0","6d","38","30",
"00","0e","0f","ab","e4","3c","eb","5e","0b","4a","d3","64","ec","1f","1a","6a","63","44",
"e8","51","1b","4f","e7","da","ff","12","6b","9d","70","f7","2e","43","3a","09","2d","42",
"5f","0c","73","c2","3a","45","9a","b8","93","03","bd","a7","4a","54","c2","c1","50","19",
"41","5f","c2","cc","60","4d","e3","3f","d4","ba","68","42","9b","72","10","f8","f8","3d",
"69","c7","5b","e5","9a","3c","a8","4a","42","33","b0","85","b7","61","d1","ee","b0","cc",
"0f","01","1a","45","7a","a6","39","cc","46","fc","5e","43","05","fd","2e","c2","d9","01",
"a5","e9","90","8c","99","25","2c","c4","74","6d","8e","7c","01","15","cf","70","de","db",
"aa","75","5c","17","5e","db","aa","1d","64","7f","a3","b4","3f","93","5a","a9","9a","ea",
"66","2d","23","dc","1d","eb","49","18","28","c0","b9","1c","92","d8","d3","cc","b9","b1",
"31","3e","9c","eb","1a","bd","c7","71","dc","88","29","ea","0a","d9","22","b5","4c","ae",
"40","9c","00","5e","0b","c4","82","ef","0c","2f","2f","2c","50","19","0b","16","b3","fd",
"4e","3e","c7","9e","2b","5c","b1","88","a4","cf","de","bd","0b","90","ce","cd","43","4d",
"90","08","c5","e3","0f","74","d8","7e","68","a4","ce","fd","a5","7a","54","96","a6","29",
"94","31","c4","79","c1","ca","b9","d3","0a","3f","f3","87","0a","ba","65","bf","24","1e",
"33","2a","43","b4","6c","7f","7b","fa","b8","b6","95","60","9c","9d","6d","da","72","bb",
"99","e2","36","de","c2","cf","62","b9","94","2a","c3","a3","91","0c","84","a9","8e","80",
"50","ef","a1","2d","16","1c","fb","5c","62","3e","7d","cd","46","07","6b","6d","58","3f",
"b7","f7","36","46","3d","5d","f7","48","01","69","ea","b3","86","9a","e4","0d","9d","e8",
"03","ff","59","35","2b","65","41","b8","71","89","80","4d","21","ec","86","bb","07","ce",
"be","9e","15","a5","3a","a3","4c","90","10","54","49","b5","8f","4f","ce","8a","7c","17",
"ff","49","b5","fe","06","ca","ce","51","5f","e2","6e","49","1d","1b","b2","59","c2","d5",
"cb","4f","be","48","28","78","9a","1a","55","9a","c1","d5","67","a1","98","83","8a","fc",
"7b","a5","50","37","6a","bb","8b","8b","38","06","45","f5","25","7f","03","c2","04","2f",
"ba","24","5e","bc","41","72","17","38","f2","83","1e","77","10","0f","d7","82","59","87",
"40","10","b6","b6","1a","b9","e0","79","f8","24","d2","e3","17","eb","c3","2c","63","61",
"80","8b","82","d9","6c","44","9f","77","b6","e0","6b","9a","78","a5","7f","9a","32","df",
"48","be","c1","ff","6c","53","ec","47","55","60","99","b5","14","ed","e4","b0","d6","18",
"60","fc","79","7f","87","4c","0c","d1","b7","c6","f8","b1","df","03","42","76","68","8f",
"bc","b0","51","e0","f0","6f","17","7d","12","08","36","6a","ae","4e","f2","79","ea","de",
"b7","49","d4","7f","9a","1b","4d","16","51","9e","68","a8","97","6b","c9","02","9e","a8",
"bf","30","2c","be","62","4c","2f","5f","5d","86","a0","eb","e0","dc","34","ac","49","1d",
"fe","bc","9f","04","b1","88","f7","43","f8","9f","56","1e","8a","0b","b5","27","4a","72",
"8e","b4","b3","fd","85","03","32","35","0c","78","d1","bd","30","31","84","46","fd","3d",
"87","0f","3c","14","db","a1","24","9f","d5","b9","b3","f0","4f","ba","03","15","28","4a",
"67","8d","29","ef","15","3f","ae","4f","47","c1","9e","a0","1d","e8","eb","30","aa","e4",
"9f","ed","20","1f","f9","a2","e5","e6","21","7a","51","f2","61","c0","a3","ff","46","4f",
"79","2b","5d","e9","00","2a","ff","31","c0","71","9f","bc","50","16","1b","30","3d","c2",
"1e","bc","ad","9b","f2","43","70","56","8d","e0","3e","ab","83","64","a9","77","6a","b3",
"bb","51","8c","35","29","ca","25","f5","e5","42","74","44","0f","89","72","9c","ab","b4",
"8a","fa","91","d5","f8","c6","24","0e","94","2c","f5","f2","fc","d0","c9","15","b7","89",
"4f","32","ed","30","c9","49","59","02","db","63","4d","7b","6b","56","26","4b","b2","76",
"15","34","b5","46","52","11","0e","65","61","b4","5c","69","43","6c","d2","dc","6d","ab",
"d8","cc","f7","80","cb","1f","99","a2","a3","53","67","0b","de","2c","de","e0","c7","54",
"20","a4","25","f0","0d","e1","ba","be","3f","2d","fa","fe","95","aa","94","79","0b","46",
"89","41","a0","a5","0d","54","42","80","e1","df","7b","ce","01","f8","b3","e5","59","dc",
"aa","b6","95","f6","c6","cb","d3","1b","82","bd","55","1f","09","28","5e","e6","f6","1f",
"ba","1f","05","ce","3a","f0","4e","57","af","3a","b1","b1","87","4b","9f","7f","e6","0b",
"b9","5e","47","38","34","c3","c6","7f","81","35","92","a1","99","d9","ab","11","7d","a8",
"d4","e6","d4","6f","f3","a0","23","5e","39","54","14","aa","26","27","d7","28","31","7e",
"f4","66","d7","c3","a6","de","55","27","69","ad","b6","bd","a0","1e","42","10","f3","06",
"b1","e8","92","5e","e4","ef","36","27","02","5a","d9","a2","95","ad","de","66","f0","92",
"e5","b3","7a","1c","53","58","21","61","70","1a","57","ed","17","33","ea","17","ef","d1",
"85","40","ec","70","16","69","b9","e2","49","52","e6","84","2b","b7","0a","7e","37","c6",
"96","91","29","df","1b","4b","37","9f","49","61","9f","7f","70","cb","46","18","ad","9c",
"41","4d","cf","ae","d1","4c","9d","e3","c9","54","87","78","cc","c7","8a","be","d9","58",
"b3","49","04","bc","9f","e3","c8","35","1d","df","53","fb","f6","10","3a","53","72","70",
"c1","c1","17","70","1d","c7","63","ff","96","ed","8e","00","ce","8b","97","1a","5c","5d",
"69","e9","57","b6","bd","03","ac","21","d2","c0","9d","6f","26","d7","f2","b5","54","9a",
"26","31","0d","ba","48","c2","d6","2e","1c","83","60","8d","8a","0c","10","66","8a","09",
"c3","e5","8e","14","61","a8","df","8d","f4","e2","98","1c","d0","33","22","3b","ad","50",
"7f","a4","de","79","de","e2","67","d3","c2","c3","ea","1f","48","cb","d9","3a","62","7b",
"b4","01","06","fb","70","f9","51","41","0e","b7","b1","46","06","3a","73","86","99","d9",
"f0","11","3c","4d","13","91","3f","99","c2","6c","a2","42","0c","fd","2f","90","60","52",
"30","59","b8","8f","ff","20","9f","fe","3c","f5","a5","be","a2","9c","9f","07","ab","14",
"0d","a6","4b","e3","2c","09","d5","08","3c","01","2c","dc","d8","85","7f","e7","5b","f9",
"cb","fa","e1","17","16","9b","a6","c5","c7","6c","8e","b1","4b","3f","5a","59","bf","6c",
"25","72","84","34","4c","1c","7f","3d","d5","96","53","ff","a7","40","85","e6","ba","d7",
"7e","69","ba","57","e5","5a","47","e5","31","b9","23","4e","fc","6d","23","b7","f1","36",
"98","8c","2b","4e","1e","0c","f4","6d","62","e3","18","92","1d","16","35","10","0d","ab",
"5c","70","77","44","53","b7","1c","68","61","b3","0c","53","5a","52","ff","45","6f","89",
"5a","fd","cb","85","da","71","8b","8c","46","87","10","e9","d9","ee","d1","59","13","c3",
"ec","00","d2","d6","b5","6d","e9","0e","92","93","e3","3c","1f","f2","fd","87","b3","f4",
"2e","96","2f","fd","b3","b1","79","e6","35","54","54","3a","71","95","de","54","39","e1",
"35","8b","90","f8","77","8b","13","17","67","23","7c","e2","18","5a","49","d2","3a","70",
"11","4f","43","f3","3e","fe","8c","0a","85","0f","64","dd","77","13","08","d8","26","45",
"4a","6c","a2","3d","ee","01","c7","a2","58","d3","61","4d","ac","59","19","e2","45","ba",
"25","06","2e","0a","ed","57","5d","73","38","8a","b5","30","eb","49","77","45","cb","c0",
"45","f4","76","87","bf","c9","65","c4","f2","fd","e7","73","6d","2e","ae","51","dd","bc",
"c3","90","27","0b","3c","16","ca","83","90","b2","e8","77","aa","43","29","a1","a8","34",
"36","3d","b8","31","46","6a","05","39","86","20","1d","70","57","d7","a5","2c","be","08",
"80","de","78","b6","64","4f","f4","38","57","e1","b3","98","9d","ea","91","11","70","2b",
"7c","a5","97","3f","c2","c2","f0","91","ac","ac","14","62","81","7b","76","06","b9","11",
"60","71","d8","49","01","c1","35","2f","af","a3","01","43","9b","97","e0","5f","83","54",
"cd","39","5d","c5","6a","19","6f","aa","fb","27","db","c5","3f","8f","7a","4c","c1","8f",
"3a","b9","87","98","f0","7d","1c","56","4c","51","f6","98","31","1f","3b","05","3c","09",
"04","be","3c","34","48","26","be","09","ce","05","60","ed","b4","bb","68","b4","e7","d9",
"03","4e","05","93","fa","b2","00","2b","62","11","a5","35","17","9d","80","5d","d5","b5",
"55","16","61","89","dd","90","db","07","57","76","38","82","7f","1d","ea","68","44","d4",
"b0","08","13","48","71","8a","1f","5b","e5","64","43","6e","1c","dc","5e","85","50","21",
"9e","41","e8","f3","63","58","23","7b","f0","31","05","ec","79","a9","a4","01","cd","de",
"3e","28","bf","dd","99","48","30","36","17","a3","e9","70","ed","9e","02","7e","4f","c7",
"3d","d0","b5","58","1d","20","f2","b8","29","02","24","2f","d3","06","42","fa","bd","a1",
"43","68","50","24","d7","3e","3d","12","fd","25","8a","ec","2f","5b","eb","10","91","b8",
"a9","44","30","52","cb","4a","8b","87","0d","ba","83","da","4a","75","f7","5f","b8","e9",
"ab","29","32","79","45","92","71","79","86","9a","c1"))

# ~~~ pressure (columns) ~~~
attr(get_kappa, "p_kPa") <- string2dat(
c("7f","0f","67","dd","58","31","66","66","9f","c8","f6","4b","f3","07","01","59","f9","f6",
"5c","9c","ab","ef","39","07","dd","7f","3e","39","6a","07","fd","6f","2b","b2","d5","c8",
"33","cd","fd","7a","ad","aa","ff","b1","55","3c","09","16","42","59","28","88","df","9d",
"cc","d4","78","29","79","d5","90","f3","c3","04","b7","34","a3","36","c1","dd","ad","b0",
"b7","95","a9","15","c1","84","f9","c2","91","48","ea","c6","d0","d9","d2","de","3f","16",
"d7","88","88","95","45","1b","ec","59","4b","a4","63","4b","88","56","c6","78","45","35"))

# ~~~ temperature (rows) ~~~
attr(get_kappa, "t_deg") <- string2dat(
c("79","e5","c5","2b","68","a8","9d","fa","99","11","a4","30","fe","8d","80","f3","f9","f6",
"5c","9c","ab","ef","39","07","dd","7f","3e","39","6a","07","fd","6f","45","66","78","64",
"3b","c5","3e","be","44","aa","ff","b1","55","3c","09","16","02","1a","97","74","db","75",
"62","d3","24","59","a9","30","62","5d","fc","25","ce","06","85","c5","25","5d","a3","df",
"30","cc","11","a2","57","a1","82","1b","4e","dc","13","b2","63","ce","07","2b","2a","5c",
"2e","64","8e","fd","a3","86","b8","64","6f","41","ec","d3","d8","8c","64","af","92","f4",
"fa","f0","f4","53","30","19","e1","a9","75","b7","dc","61","82","35","e9","4b","26","99",
"80","20","85","f5","7e","15","bf","97","73","df","68","f9","2f","38","ef","59","e5","fc",
"27","b8","d4","2c","c3","7f","d7","fb","b1","e7","c9","9b","78","de","2d","56","b3","ee",
"bb","a6","36","fc","67","c5","d3","86","76","f7","7b","35","e3","2e","49","a1","a0","12",
"81","f5","17","49","a9","a8","e2","6f","c5","04","55","07","d1","ba","23","e4","eb","d6",
"b3","fb","ae","d0","31","ed","06","93","d0","5c","92","e3","05","12","a2"))

