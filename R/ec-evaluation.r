## documentation (parts) ec-processing function ----------------------------------------

# 1. despiking of time series
# default: scalar data are despiked, sonic data is not
# despike = c(u = FALSE, v = FALSE, w = FALSE, T = FALSE, 
#    nh3_ppb = TRUE, nh3_ugm3 = TRUE, h2o_mmolm3 = TRUE, 
#    co2_mmolm3 = TRUE)
#   => can be switched on/off for all, e.g. despike = FALSE
# despike_baseline_width = c(u = 10, v = 10, w = 10, T = 10, nh3_ppb = 10, 
#    nh3_ugm3 = 10, h2o_mmolm3 = 10, co2_mmolm3 = 10)
#   => Blackman-Nuttall hight-pass filter is applied on time series to get a baseline
#       "despike_baseline_width" is the filter window width in seconds
#       defaults to 10 secs
# despike_quantile = c(u = 0.95, v = 0.95, w = 0.95, T = 0.95, nh3_ppb = 0.95, 
#    nh3_ugm3 = 0.95, h2o_mmolm3 = 0.95, co2_mmolm3 = 0.95)
#   => statistics on the baseline, timeseries and their difference (d = timeseries - baseline) is gathered.
#       The statistics is calculated as sd(baseline) + mad(timeseries) + quantile(d, q),
#       where q is given by "despike_quantile". Defaults to 0.95
# despike_stats_width = c(u = 30, v = 30, w = 30, T = 30, nh3_ppb = 30, 
#    nh3_ugm3 = 30, h2o_mmolm3 = 30, co2_mmolm3 = 30)
#   => width (in seconds) of the statistics window, i.e. the window where the statistics is calculated within
#       defaults to 30 secs
# despike_stats_multiply = c(u = 4, v = 4, w = 4, T = 4, nh3_ppb = 4, 
#    nh3_ugm3 = 4, h2o_mmolm3 = 4, co2_mmolm3 = 4)
#   => multiplication of the statistics value s for the final filtering band.
#       the filtering criteria is given by abs(d) > (baseline + s * "despike_stats_multiply")

# 2. check hard limits and replace/interpolate NA values
# default: hard limits are applied to all time series and NA values are interpolated by a weighted (w = 1/r^2) mean
# na_limits = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, 
#    nh3_ppb = TRUE, nh3_ugm3 = TRUE, h2o_mmolm3 = TRUE, 
#    co2_mmolm3 = TRUE)
# limits_lower = c(u = -30, v = -30, w = -10, T = 243, 
#    nh3_ppb = -100, nh3_ugm3 = -100, h2o_mmolm3 = -100, 
#    co2_mmolm3 = -100)
# limits_upper = c(u = 30, v = 30, w = 10, T = 333, 
#    nh3_ppb = 22000, nh3_ugm3 = 15000, h2o_mmolm3 = 5000, 
#    co2_mmolm3 = 5000)
# na_limits_method = c('norepl', 'median', 'dist', 'squaredist')[4]
#   => norepl: no replacing, 
#       median: median(y_win, na.rm = TRUE), 
#       dist: mean(w * y_win, na.rm = TRUE) where w = 1/d_x,
#       squaredist: mean(w * y_win, na.rm = TRUE) where w = 1/d_x ^ 2 (closer values get higher weights)
# na_limits_window = c(pass = '10secs', replace = '5mins')
#   => windows -> pass: consecutive NA values longer than xx seconds are not replaced (default = 10 secs)
#              -> replace: window width on which interpolation of central NA value is based (default = 5mins)

# x. Sonic Rotation Method
# rotation_method = c("two axis", "planar fit")
# default: two-axis rotation
# two-axis rotation as usually done (rotat

## helper functions ----------------------------------------

# check (and optionally replace) 'hard' data limits
check_limits <- function(dat, limits, lim_window = 500, 
    lim_method = c('norepl', 'median', 'dist', 'squaredist')[4],
    d_t, bin_threshold = NULL){ 
    # check data.table
    if (!is.data.table(dat)) {
        cat('Fixme in "check_limits" - argument "dat" is not a data.table...\n')
        browser()
    }
    # get variable names
    hl_vars <- colnames(limits)
    # mark values outside of limits
    dat[, paste0(hl_vars, '_flag') := mapply(\(x, nm) {
        out <- as.integer(x < limits['lower', nm] | x > limits['upper', nm])
        out[is.na(out)] <- 2
        out
        }, x = .SD, nm = hl_vars, SIMPLIFY = FALSE), .SDcols = hl_vars]
    # check method
    if (!pmatch(lim_method, 'norepl', nomatch = 0)) {
        # re-flag windows which can't be replaced anyway
        # parse windows
        pass_window <- parse_time_diff(lim_window[['pass']]) * d_t
        # check max NA window lengths (pass) -> set flag to -1 if > window
        dat[, paste0(hl_vars, '_flag') := lapply(.SD, \(x) {
            y <- x > 0
            z <- y[-1L] != y[-.N]
            i <- which(z)
            i_to <- c(i, .N)
            i_from <- c(0, i)
            i_ind <- diff(c(0, i_to)) > pass_window & y[i_to]
            if (any(i_ind)) {
                ind <- unlist(mapply(':', i_from[i_ind], i_to[i_ind], SIMPLIFY = FALSE))
                x[ind] <- -1
            }
            x
        }), .SDcols = paste0(hl_vars, '_flag')]
        # check flagged
        hflgs <- dat[, sapply(.SD, \(x) sum(x > 0)), .SDcols = paste0(hl_vars, '_flag')]
        # get variable to replace
        var_hflgs <- hl_vars[hflgs > 0]
        if (any(hflgs > 0)) {
            # parse 'replace' window
            repl_window <- parse_time_diff(lim_window[['replace']])
            # get window limits
            win_half <- repl_window / 2
            # check method
            switch(lim_method
                , 'dist' = {
                    # get weighting (1 / r)
                    w1 <- 1 / seq_len(win_half * d_t)
                    w_1 <- c(rev(w1), 1, w1)
                }
                , 'squaredist' = {
                    # get weighting (1 / r ^ 2)
                    w1 <- 1 / seq_len(win_half * d_t) ^ 2
                    w_1 <- c(rev(w1), 1, w1)
                }
                , 'median' = {
                    w_1 <- NULL
                }
                , stop('replacement method not valid')
            )
            # loop over variables
            if (is.null(w_1)) {
                # be verbose
                cat("Replacing flagged values by window median...\n")
                # median
                dat[, (var_hflgs) := lapply(var_hflgs, \(var) {
                    cat('\t-', var, '- ')
                    x <- get(var)
                    x_flag <- get(paste0(var, '_flag'))
                    isna <- which(x_flag > 0)
                    x1 <- as.numeric(Time[isna]) - win_half
                    x2 <- as.numeric(Time[isna]) + win_half
                    # find window for each NA
                    cat('get windows ')
                    max_size <- win_half * 2 * Hz[1] * 1.2
                    ind <- find_window(as.numeric(Time), x1, x2, max_size)
                    # replace values
                    cat(' - replace values: ')
                    x[isna] <- sapply(seq_along(isna), \(i) {
                        verb <- paste(i, '/', length(isna))
                        cat(verb)
                        out <- median(x[ind[[i]]], na.rm = TRUE)
                        cat(paste(rep('\b', nchar(verb)), collapse = ''))
                        out
                    })
                    cat('\n')
                    x
                })]
            } else {
                # be verbose
                cat("Replacing flagged values by window weighted average\n")
                # weighted mean
                dat[, (var_hflgs) := lapply(var_hflgs, \(var) {
                    cat('\t-', var, '- ')
                    x <- get(var)
                    x_flag <- get(paste0(var, '_flag'))
                    isna <- which(x_flag > 0)
                    # dat[749718 + (-1500:1500),]
                    x1 <- as.numeric(Time[isna]) - win_half
                    x2 <- as.numeric(Time[isna]) + win_half
                    # find window for each NA
                    cat('get windows ')
                    max_size <- win_half * 2 * Hz[1] * 1.2
                    ind <- find_window(as.numeric(Time), x1, x2, max_size)
                    # replace values
                    cat(' - replace values: ')
                    x[isna] <- sapply(seq_along(isna), \(i) {
                        verb <- paste(i, '/', length(isna))
                        cat(verb)
                        j <- setdiff(ind[[i]], 0L)
                        if (length(j) != length(w_1)) {
                            w_i <- j - isna[i] + win_half * d_t + 2
                            out <- sum(x[j] * w_1[w_i], na.rm = TRUE) / 
                                sum(w_1[w_i][x_flag[j] == 0])
                        } else {
                            out <- sum(x[j] * w_1, na.rm = TRUE) / sum(w_1[x_flag[j] == 0])
                        }
                        cat(paste(rep('\b', nchar(verb)), collapse = ''))
                        out
                    })
                    cat('\n')
                    x
                })]
            }
            cat("\n*~~~~*\n")
        }
	} else {
        # set values outsid limits to NA
	    if (dat[, any(sapply(.SD, \(x) any(x == 1))), 
            .SDcols = paste0(hl_vars, '_flag')]) {
            cat('setting values outside limits to NA...\n')
            # set flagged values to NA
            dat[, (hl_vars) := lapply(hl_vars, \(var) {
                xf <- get(paste0(var, '_flag'))
                x <- get(var)
                x[xf == 1] <- NA_real_
                x
                })
            ]
        }
    }
    # check remaining NA values
    hflgs <- dat[, sapply(.SD, \(x) sum(is.na(x))), .SDcols = hl_vars]
    names(hflgs) <- hl_vars
    var_hflgs <- names(hflgs)[hflgs > 0]
    if (length(var_hflgs) > 0) {
        cat(
            paste0(
                '------\nremaining NA-values in data set:\n',
                paste(
                    sprintf('%s: %i', var_hflgs, hflgs[var_hflgs]), 
                    collapse = '\n'),
                "\n"
            )
        )
    }
    # check interval threshold
    if (!is.null(bin_threshold)) {
        dat[, paste0(hl_vars, '_flag') := lapply(hl_vars, \(var) {
            xv <- get(var)
            # check na
            isna <- is.na(xv)
            if (all(isna)) {
                # all NA values
                return(rep(-3, .N))
            }
            x <- get(paste0(var, '_flag'))
            # check run length again
            y <- !(isna | x < 0)
            z <- y[-1L] != y[-.N]
            i <- which(z)
            i_to <- c(i, .N)
            if (any(y[i_to])) {
                di <- diff(c(0, i_to))
                # get max.
                i1 <- which(y[i_to])
                i2 <- i1[which.max(di[i1])]
                if (di[i2] > bin_threshold) {
                    # get index
                    i_from <- c(0, i)
                    ind <- i_from[i2]:i_to[i2]
                    x[-ind] <- -1
                } else {
                    # invalid interval (too few continuous data)
                    x[] <- -2
                }
            } else {
                # invalid interval (no data)
                x[] <- -2
            }
            x
        }), by = bin]
    }
	invisible(dat)
}

# filter functions for higph-pass filtering time series
# (copied from minidoas scripts)
# https://en.wikipedia.org/wiki/Window_function
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
      x <- sin(pi * seq.int(0, N) / N) ^ 2 / (0.5 * N)
      x / sum(x)
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
      x / sum(x)
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
    "Exp" = function(n, p1 = min(n / 2, 30), ...){
      N <- (n - 1) / 2
      x <- exp(p1 * sqrt(1 - (seq.int(-N, N) / N) ^ 2)) / exp(p1)
      x <- na.omit(x)
      x / sum(x)
    },
    # Poisson:
    "Poisson" = function(n, p1 = n / 2 / 6.9, ...){
      N <- n - 1
      x <- exp(-abs(seq.int(0, N) - N / 2) / p1)
      x <- na.omit(x)
      x / sum(x)
    }
)

# detrending time series
trend <- function(y, method = c("blockAVG", "linear", "linear_robust", "ma_360"), Hz_ts = 10) {
    if (is.null(naa <- na.action(y))) {
        y <- na.omit(y)
        naa <- na.action(y)
    }
    yfin <- y
    n <- length(y)
    fitted <- rep(NA_real_, length(y))
	method <- method[1]
	switch(method, 
		"blockAVG" = {
			my <- mean(y, na.rm = TRUE)
			fitted <- rep(my, n) 
            resid <- y - fitted
            attr(resid, 'na.action') <- attr(fitted, 'na.action') <- naa
			list(
				coefficients = c(intercept = my, slope = 0)
				, fitted = fitted
				, residuals = resid
				) 
		}
		, "linear" = {
			my <- mean(yfin)
			mx <- (n + 1) / 2
			xstr <- (x <- seq.int(n)) - mx
			b <- sum(xstr * (yfin - my)) / (n - 1) / (sum(xstr ^ 2) / (n - 1))
			a <- my - mx * b
			fitted <- a + b * x
            resid <- y - fitted
            attr(resid, 'na.action') <- attr(fitted, 'na.action') <- naa
			list(
				coefficients = c(intercept = a, slope = b)
				, fitted = fitted
				, residuals = resid
				) 
		}
		, "linear_robust" = {
			mod <- MASS::rlm(yfin ~ seq.int(n))
			cfs <- mod$coefficients
			a <- cfs[1]
			b <- cfs[2]
			fitted <- mod$fitted
            resid <- y - fitted
            attr(resid, 'na.action') <- attr(fitted, 'na.action') <- naa
			if (!mod$converged) {
				stop("robust linear regression did not converge!")
			}
			list(
				coefficients = c(intercept = a, slope = b)
				, fitted = fitted
				, residuals = resid
				) 
		}
		,{
            if (grepl('^ma_', method)) {
                ma <- round(as.numeric(sub("ma_", "", method)) * Hz_ts)
                fitted <- caTools::runmean(yfin, ma, "C", "mean")
                resid <- y - fitted
                attr(resid, 'na.action') <- attr(fitted, 'na.action') <- naa
                list(
                    coefficients = c(intercept = NA, slope = NA)
                    , fitted = fitted
                    , residuals = resid
                ) 
            } else if (sub('_\\d+$', '', method) %in% names(filter_list)) {
                win <- round(as.numeric(sub(".*_", "", method)) * Hz_ts)
                if (win >= (length(y) + length(na.action(y)))) {
                    stop('(detrending) filter window ist larger than averaging time!')
                }
                if (win >= n) {
                    warning('(detrending) filter window ist larger than available data')
                    resid <- y - fitted
                    attr(resid, 'na.action') <- attr(fitted, 'na.action') <- naa
                    return(
                        list(
                            coefficients = c(intercept = NA, slope = NA)
                            , fitted = fitted
                            , residuals = resid
                        )
                    )
                }
                filter_name <- sub('_\\d+$', '', method)
                # lowpass filter
                filt <- filter_list[[filter_name]](win)
                ihalf <- 1:(win / 2)
                yh1 <- mean(yfin[ihalf])
                yh2 <- mean(yfin[n + 1 - ihalf])
                ym <- yfin - (yh2 - yh1)
                yp <- yfin - (yh1 - yh2)
                z <- c(ym[n + 1 - ihalf], yfin, yp[ihalf])
                fitted <- stats::filter(z, filt, 'convolution', 2, 
                    circular = FALSE)[1:n + ihalf[length(ihalf)]]
                resid <- y - fitted
                attr(resid, 'na.action') <- attr(fitted, 'na.action') <- naa
                list(
                    coefficients = c(intercept = NA, slope = NA)
                    , fitted = fitted
                    , residuals = resid
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
	thetam <- atan2(mean(v, na.rm = TRUE), mean(u, na.rm = TRUE))
	# horizontal rotation
	u1 <- u * cos(thetam) + v * sin(thetam)
    # apply vertival rotation?
	if (is.null(phi)) {
		phi <- atan2(mean(w, na.rm = TRUE), mean(u1, na.rm = TRUE))
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
    mi <- seq(max(m - 500, 1), min(m + 500, length(x)))
    avg <- mean(x[mi], na.rm = TRUE)
	maxis <- ind[which.max(abs(x[ind] - avg))]
	c(index = maxis, tau = maxis - m)
}
fixlag_index <- function(x, fix) {
	n <- length(x)
	if (n %% 2) {
		# ungerade
		m <- (n + 1) / 2
	} else {
		# gerade
		m <- n / 2 + 1
	}
    m + fix
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

# helper function for ec_main
detrend_sonic_data <- function(x, detr, rhz) {
    wind <- x[, {
        ud <- trend(urot, detr["u"], rhz)
        vd <- trend(vrot, detr["v"], rhz)
        wd <- trend(wrot, detr["w"], rhz)
        Td <- trend(T, detr["T"], rhz)
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
    # keep T for plotting
    x[, Trot := T]
    # replace detrended variables
    ind <- x[, seq_len(.N)]
    vkey <- c("u" = 'uprot', "v" = 'vprot', "w" = 'wprot', "T" = 'Tdet')
    for (w in names(vkey)) {
        vcol <- vkey[[w]]
        naa <- na.action(wind[[vcol]])
        if (is.null(naa)) {
            set(x, , w, wind[[vcol]])
        } else {
            set(x, i = ind[-naa], w, wind[[vcol]])
        }
    }
    # return wind
    wind
}

# calculate wind statistics and MOST parameters
wind_statistics <- function(wind, z_canopy, z_sonic, 
    ustar_method = c('neg_sqrt', 'double_sqrt', 'fallback')[1]) {
    ustar_method <- match.arg(ustar_method, 
        c('neg_sqrt', 'double_sqrt', 'fallback'))
	Cov_sonic <- cov(list2DF(wind[1:4]), use = 'na.or.complete')
	Var_sonic <- diag(Cov_sonic)
	names(Var_sonic) <- c('var_u', 'var_v', 'var_w', 'var_T')
	Cov_sonic <- Cov_sonic[cbind(c("uprot","uprot","uprot","vprot","vprot","wprot"),c("vprot","wprot","Tdet","wprot","Tdet","Tdet"))]
	names(Cov_sonic) <- c('cov_uv', 'cov_uw', 'cov_uT', 'cov_vw', 'cov_vT', 'cov_wT')
    if (ustar_method %in% c('neg_sqrt', 'fallback')) {
        suppressWarnings(Ustar <- c(sqrt(-Cov_sonic["cov_uw"]),use.names = FALSE))
    }
    if (
        ustar_method == 'double_sqrt' || 
        (ustar_method == 'fallback' && is.na(Ustar))
    ) {
        Ustar <- c(sqrt(sqrt(Cov_sonic['cov_uw'] ^ 2 + Cov_sonic['cov_vw'] ^ 2)),
            use.names = FALSE)
    }
	T_K <- mean(wind$Tmdet + wind$Tdet, na.rm = TRUE)
	U <- mean(wind$umrot + wind$uprot, na.rm = TRUE)
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

# aerodynamic resistance of heat & trace gases
Ra <- function(z, ustar, L, z0, d = 0, kv = 0.4) {
    z_d <- z - d
    1 / (kv * ustar) * (log(z_d / z0) - psi_h(z_d, L) + psi_h(z0, L))
}
psi_h <- function(z_d, L) {
    zL <- z_d / L
    suppressWarnings(ifelse(zL >= 0, -5 * zL, 2 * log((1 + sqrt(1 - 16 * zL)) / 2)))
}
# pseudo-laminar layer resistance Rb of NH3
Rb_nh3 <- function(ustar, Tc, Zo, p_hPA = 960) {
    Diff.Ammon <- function(Tc, p_hPA) {
        # Temp <- c(293, 296, 298)
        # Dm <- c(173, 180, 173) 
        # DopoTo <- mean(Dm / (760 / 1013.25) / Temp ^ 1.5)
        out <- 0.04598255 / p_hPA * (Tc + 273.15) ^ 1.5 
        out * 10 ^ -4 
    }
    kin.visc <- function(Tc, p_hPA) {
        Tk <- Tc + 273.15
        b <- 1.458 * 10 ^ -6
        S <- 110.4
        p <- p_hPA * 100
        R <- 8.314
        M <- 29
        (b * Tk ^ (5 / 2) * R * 1000) / (p * M * (Tk + S))
    }
    1.45 * (Zo * ustar / kin.visc(Tc, p_hPA)) ^ 0.24 * 
        (kin.visc(Tc, p_hPA) / Diff.Ammon(Tc, p_hPA)) ^ 0.8 / ustar
}

# ec_main helpers
get_covariance <- function(dt, cvars, cvars_name, nperiod) {
    cov_out <- dt[, I(lapply(cvars, function(i) {
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
        n_missing <- nperiod - N
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
        structure(out, ffts = list(x = xfft, y = yfft), 
            isfinite = isfinite)
    }))]
    # add names
    names(cov_out) <- cvars_name
    cov_out
}
get_cospectra <- function(dt, cvrs, cvars, cvars_name, lag_tau, nperiod) {
    out <- mapply(function(var, lag) {
            # get covar
            covars <- cvrs[[paste(var, collapse = 'x')]]
            # get ffts
            xfft <- attr(covars, 'ffts')$x
            # get length
            N <- length(xfft)
            # get y
            if (lag != 0) {
                y <- dt[attr(covars, 'isfinite'), get(var[2])]
                yfft <- fft(data.table::shift(y, lag, type = 'cyclic')) / N
            } else {
                yfft <- attr(covars, 'ffts')$y
            }
            # get cospec
            re <- Re(Conj(yfft) * xfft)[seq(N / 2) + 1] * N / (N - 1) * 2
            # get missing
            n_missing <- nperiod / 2 - length(re)
            if (sign(n_missing) >= 0) {
                c(re, rep(0, n_missing))
            } else {
                re[seq_len(nperiod / 2)]
            }
        }, var = cvars, lag = lag_tau, SIMPLIFY = FALSE)
    # add names
    names(out) <- cvars_name
    out
}


# plot time series including filtered trend
plot.tseries <- function(dat,wind,scal,selection,color,units,st_sub=NULL){
	msg <- paste(c(format(dat[1,1],"%d.%m.%Y")," - time series"),collapse="")
    # fix wind variables
	dat[, c("u", "v", "w", "T")] <- dat[, paste0(c('u', 'v', 'w', 'T'), 'rot')]
	# dat[, c("u", "v", "w", "T")] <- mapply("+", 
        # wind[c("umrot", "vmrot", "wmrot", "Tmdet")], 
        # wind[c("uprot", "vprot", "wprot", "Tdet")])
	### get trends:
	dat3 <- list2DF(wind[c("umrot","vmrot","wmrot","Tmdet")])
	names(dat3) <- c("u","v","w","T")
	if(!is.null(scal)){
        # fix scalars
        dat[, names(scal)] <- lapply(scal, \(x) {
            if (is.null(isna <- na.action(x$fitted))) {
                x$fitted + x$residuals
            } else {
                out <- rep(NA_real_, nrow(dat))
                out[-isna] <- x$fitted + x$residuals
                out
            }
        })
        # add scalars
		dat3 <- cbind(dat3, lapply(scal, \(x) {
                if (is.null(isna <- na.action(x$fitted))) {
                    x$fitted
                } else {
                    out <- rep(NA_real_, nrow(dat3))
                    out[-isna] <- x$fitted
                    out
                }
            })
        )
    }
    # add residuals to dat3
    res <- setdiff(selection, names(dat3))
    if (length(res) > 0) {
        dat3 <- cbind(dat3, dat[, res, drop = FALSE])
    }
	### melt and add trends:
	dat4 <- reshape2::melt(dat3[, selection], id = NULL, value.name = "trend")
	dat2 <- reshape2::melt(dat[,c("st",selection)],id="st")
	dat2[, "trend"] <- dat4[, "trend"]
	myxscale.component <- function(...) {
		ans <- lattice::xscale.components.default(...)
		ans$top <- ans$bottom
		ans
	}
	myyscale.component <- function(...) {
		ans <- lattice::yscale.components.default(...)
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
            lattice::panel.xyplot(x,y,...)
			y2 <- dat2[list(...)$subscripts,"trend"]
            if (!identical(y, y2)) {
                # panel.xyplot(x,y2,type="l",lwd=1.5, lty=3, col="gray30")
                # panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightblue")
                # panel.xyplot(x,y2,type="l",lwd=2, lty=2, col="lightgrey")
                # panel.xyplot(x, y2, type = "l", lwd = 2, lty = 2, col = "#B37FDF")
                lattice::panel.xyplot(x, y2, type = "l", lwd = 3, col = "darkgrey")
            }
            if (!is.null(st_sub)) lattice::panel.abline(v = st_sub, lty = 2, lwd = 2, col = 'grey')
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
    col = "lightblue", nred = floor(sqrt(sqrt(length(ogive))) * 1.5), model_par = NULL,
    model_cols = '#AE71EB99') {
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
	plot(1, xlim = xlim, ylim = ylim, cex.axis = cx, cex.lab = cx, type = "n", 
        log = "x", xaxt = "n", yaxt = "n", xlab = "frequency [Hz]", ylab = "", 
        panel.first = abline(h = 0, col = col, lty = 2))
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
        if (length(model_cols) == 1) model_cols <- rep(model_cols, 2)
        # cospec: f * Co(f)
        cs <- freq * cospec_model(model_par['fx'], model_par['m'], model_par['mu'], 
            model_par['A0'], freq)
	    y_cs <- (cs - min(prCo)) / diff(range(prCo)) * diff(ylim) + ylim[1]
        lines(freq, y_cs, col = model_cols[1], lwd = 2)
        # ogive
        og <- ogive_model(model_par['fx'], model_par['m'], model_par['mu'], 
            model_par['A0'], freq)
        lines(freq, og, col = model_cols[2], lwd = 2)
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

## convenience wrapper for turbulence only
process_turbulence <- function(
    sonic_directory
    , start_time = NULL
    , end_time = NULL
    , z_sonic = NULL
    , dev_north = NULL
    , declination = NULL
    , z_canopy = NULL
	, avg_period = '30mins'
    , subintervals = FALSE
    , detrending = c(u = 'blockAVG', v = 'blockAVG', w = 'blockAVG', T = 'blockAVG')
    , ustar_method = c('neg_sqrt', 'double_sqrt', 'fallback')[1]
    , na_limits = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE)
    , limits_lower = c(u = -30, v = -30, w = -10, T = 243)
    , limits_upper = c(u = 30, v = 30, w = 10, T = 333)
    , na_limits_method = c('norepl', 'median', 'dist', 'squaredist')[4]
    , na_limits_window = c(pass = '10secs', replace = '5mins')
    , covariances = c('uxw', 'wxT')
    , create_graphs = FALSE
    , ogive_out = FALSE
    , ...
    ) {
    # remove covar/cospec/ogive calculation if graphs are not created and ogives are not
    #   included in output
    if (!create_graphs && !ogive_out) {
        covariances <- NULL
    }
    process_ec_fluxes(
        sonic_directory = sonic_directory
        , start_time = start_time
        , end_time = end_time
        , z_ec = z_sonic
        , dev_north = dev_north
        , declination = declination
        , z_canopy = z_canopy
		, avg_period = avg_period
        , subintervals = subintervals
        , detrending = detrending
        , na_limits = na_limits
        , limits_lower = limits_lower
        , limits_upper = limits_upper
        , na_limits_window = na_limits_window
        , na_limits_method = na_limits_method
        , covariances = covariances
        , create_graphs = create_graphs
        , ogive_out = ogive_out
        , ustar_method = ustar_method
        , ...
    )
}

## main flux processing function
process_ec_fluxes <- function(
		sonic_directory
        , ht_directory = NULL
        , licor_directory = NULL
        , miro_directory = NULL
		, start_time = NULL
		, end_time = NULL
		, avg_period = '30mins'
		, tz_user = 'UTC'
		, dev_north = NULL
        , declination = NULL
		, z_ec = NULL
		, z_canopy = NULL
        , despike = c(u = FALSE, v = FALSE, w = FALSE, T = FALSE, 
            nh3_ppb = TRUE, nh3_ugm3 = TRUE, h2o_mmolm3 = TRUE, 
            co2_mmolm3 = TRUE, h2o_molfrac = TRUE, ch4_molfrac = TRUE,
            n2o_molfrac = TRUE
        )
        , despike_baseline_width = c(u = 10, v = 10, w = 10, T = 10, 
            nh3_ppb = 10, nh3_ugm3 = 10, h2o_mmolm3 = 10, co2_mmolm3 = 10,
            h2o_molfrac = 10, ch4_molfrac = 10, n2o_molfrac = 10
        )
        , despike_quantile = c(u = 0.95, v = 0.95, w = 0.95, T = 0.95, 
            nh3_ppb = 0.95, nh3_ugm3 = 0.95, h2o_mmolm3 = 0.95, 
            co2_mmolm3 = 0.95, h2o_molfrac = 0.95, ch4_molfrac = 0.95,
            n2o_molfrac = 0.95
        )
        , despike_stats_width = c(u = 30, v = 30, w = 30, T = 30, 
            nh3_ppb = 30, nh3_ugm3 = 30, h2o_mmolm3 = 30, co2_mmolm3 = 30,
            h2o_molfrac = 30, ch4_molfrac = 30, n2o_molfrac = 30
        )
        , despike_stats_multiply = c(u = 4, v = 4, w = 4, T = 4, 
            nh3_ppb = 4, nh3_ugm3 = 4, h2o_mmolm3 = 4, co2_mmolm3 = 4,
            h2o_molfrac = 4, ch4_molfrac = 4, n2o_molfrac = 4
        )
        , na_limits = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, 
            nh3_ppb = TRUE, nh3_ugm3 = TRUE, h2o_mmolm3 = TRUE, 
            co2_mmolm3 = TRUE, h2o_molfrac = TRUE, ch4_molfrac = TRUE,
            n2o_molfrac = TRUE
        )
        , limits_lower = c(u = -30, v = -30, w = -10, T = 243, 
            nh3_ppb = -100, nh3_ugm3 = -100, h2o_mmolm3 = -100, 
            co2_mmolm3 = -100, h2o_molfrac = -0.1, ch4_molfrac = -0.1,
            n2o_molfrac = -0.1
        )
        , limits_upper = c(u = 30, v = 30, w = 10, T = 333, 
            nh3_ppb = 2.2e4, nh3_ugm3 = 1.5e4, h2o_mmolm3 = 5e3, 
            co2_mmolm3 = 5e3, h2o_molfrac = 0.10, ch4_molfrac = 200e-6,
            n2o_molfrac = 20e-6
        )
        , na_limits_window = c(pass = '10secs', replace = '5mins')
        , na_limits_method = c('norepl', 'median', 'dist', 'squaredist')[4]
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
        , detrending = c(u = 'blockAVG', v = 'blockAVG', w = 'blockAVG', 
            T = 'blockAVG', nh3_ppb = 'blockAVG', nh3_ugm3 = 'blockAVG', 
            h2o_mmolm3 = 'blockAVG', co2_mmolm3 = 'blockAVG',
            h2o_molfrac = 'blockAVG', ch4_molfrac = 'blockAVG',
            n2o_molfrac = 'blockAVG'
        )
		, covariances = c('uxw', 'wxT', 'wxnh3_ugm3', 'wxh2o_mmolm3', 
            'wxco2_mmolm3', 'wxh2o_molfrac', 'wxch4_molfrac', 
            'wxn2o_molfrac'
        )
        # fix lag in seconds
		, lag_fix = c(uxw = 0, wxT = 0, wxnh3_ppb = -0.4, 
            wxnh3_ugm3 = -0.4, wxh2o_mmolm3 = -0.2, wxco2_mmolm3 = -0.2,
            wxh2o_molfrac = -0.7, wxch4_molfrac = -0.7, wxn2o_molfrac = -0.7
        )
        # dyn lag in seconds around lag_fix
		, lag_dyn = c(uxw = 0.5, wxT = 0.5, wxnh3_ppb = 1.5, 
            wxnh3_ugm3 = 1.5, wxh2o_mmolm3 = 1.5, wxco2_mmolm3 = 1.5,
            wxh2o_molfrac = 1.5, wxch4_molfrac = 1.5, wxn2o_molfrac = 1.5
        )
        # which dyn lag approach should be taken?
		, lag_dyn_calc_pw = c(uxw = FALSE, wxT = FALSE, wxnh3_ugm3 = FALSE, 
            wxh2o_mmolm3 = FALSE, wxco2_mmolm3 = FALSE, 
            wxh2o_molfrac = FALSE, wxch4_molfrac = FALSE, 
            wxn2o_molfrac = FALSE
        )
        , lag_dyn_method = c('raw-cov', 'simple-pw', 'boot-pw')
        # , lag_dyn_wdt = 5 # suggested in RFlux
        , lag_dyn_wdt = 7 
        , lag_dyn_lagmax = 10 # dito
        , lag_dyn_model = c('ar', 'arima')[1]
        # , lag_dyn_smooth = 'mean' # original Vitale smoothing
        , lag_dyn_smooth = function(x, wdt) {
            flt <-  getOption('md.filter.function.list')$BmNuttall(wdt)
            sum(x * flt, na.rm = TRUE) / sum(flt, na.rm = TRUE)
        }
        # re_rmse window
        , gamma_time_window = c(5, 10)
        # damping_reference: either a specific covariance as 'wxT' or best quality ogives
        #   such as either 'base_quality' (for best base_quality_fix/_dyn)
        #   or 'ogive_quality' (for best ogive_quality_fix/_dyn)
		, damping_reference = c(wxnh3_ppb = 'wxT', wxnh3_ugm3 = 'wxT', 
            wxh2o_mmolm3 = 'wxT', wxco2_mmolm3 = 'wxT', 
            wxh2o_molfrac = 'wxT', wxch4_molfrac = 'wxT', 
            wxn2o_molfrac = 'wxT'
        )
        # lower & upper bounds of fitting ogives (in seconds)
		, damping_lower = c(wxnh3_ppb = 2, wxnh3_ugm3 = 2, 
            wxh2o_mmolm3 = 2, wxco2_mmolm3 = 2, wxh2o_molfrac = 2, 
            wxch4_molfrac = 2, wxn2o_molfrac = 2
        )
		, damping_upper = c(wxnh3_ppb = 20, wxnh3_ugm3 = 20, 
            wxh2o_mmolm3 = 20, wxco2_mmolm3 = 20, wxh2o_molfrac = 20,
            wxch4_molfrac = 20, wxn2o_molfrac = 20
        )
        , low_cont_sec = 20
        , high_cont_sec = 2
        , cont_pts = 5
        , subintervals = TRUE
        , subint_prefix = 'subint_'
        , subint_n = 5
        , subint_detrending = c(u = 'linear', v = 'linear', w = 'linear', 
            T = 'linear', nh3_ppb = 'linear', nh3_ugm3 = 'linear', 
            h2o_mmolm3 = 'linear', co2_mmolm3 = 'linear', 
            h2o_molfrac = 'linear', ch4_molfrac = 'linear', 
            n2o_molfrac = 'linear'
        )
        , subint_return = FALSE
        , oss_threshold = 0
        , co2ss_threshold = 0
        , na_alarm_code = c(1:3, 5:8, 11, 13)
        , thresh_period = 0.75
        , ustar_method = c('neg_sqrt', 'double_sqrt', 'fallback')[1]
		, create_graphs = TRUE
		, create_dailygraphs = TRUE
		, graphs_directory = NULL
		, add_name = ''
        , plot_timeseries = c(u = TRUE, v = TRUE, w = TRUE, T = TRUE, 
            ht_oss = TRUE, nh3_ppb = FALSE, nh3_ugm3 = TRUE, 
            li_co2ss = TRUE, h2o_mmolm3 = TRUE, co2_mmolm3 = TRUE)
        , plotting_var_units = c(u = 'm/s', v = 'm/s', w = 'm/s', 
            T = 'K', ht_oss = '-', nh3_ppb = 'ppb', nh3_ugm3 = 'ug/m3', 
            li_co2ss = '-', h2o_mmolm3 = 'mmol/m3', co2_mmolm3 = 'mmol/m3',
            h2o_molfrac = 'mol/mol', ch4_molfrac = 'mol/mol', 
            n2o_molfrac = 'mol/mol'
        )
        , plotting_var_colors = c(u = 'gray20', v = 'gray20', w = 'gray20', 
            T = 'orange', ht_oss = 'grey', nh3_ppb = 'indianred', 
            nh3_ugm3 = 'indianred', li_co2ss = 'grey', 
            h2o_mmolm3 = '#8FC1E6', co2_mmolm3 = 'seagreen4',
            h2o_molfrac = '#3A9BBB', ch4_molfrac = '#EBD400', 
            n2o_molfrac = 'purple'
        )
        , plotting_covar_units = c(uxw = 'm2/s2', wxT = 'K*m/s', 
            wxnh3_ppb = 'ppb*m/s', wxnh3_ugm3 = 'ug/m2/s', 
            wxh2o_mmolm3 = 'mmol/m2/s', wxco2_mmolm3 = 'mmol/m2/s',
            wxh2o_molfrac = 'mol/mol*m/s', wxch4_molfrac = 'mol/mol*m/s',
            wxn2o_molfrac = 'mol/mol*m/s'
        )
        , plotting_covar_colors = c(uxw = 'gray70', wxT = 'orange', 
            wxnh3_ppb = 'indianred', wxnh3_ugm3 = 'indianred', 
            wxh2o_mmolm3 = '#8FC1E6', wxco2_mmolm3 = 'seagreen4',
            wxh2o_molfrac = '#3A9BBB', wxch4_molfrac = '#EBD400',
            wxn2o_molfrac = 'purple'
        )
        , model_colors = c(cospec = '#F02E42', ogive = '#9A33DA')
		, ogives_return = FALSE
        , as_ibts = TRUE
        , ncores = 1
        , psock_args = NULL
        , debug_parallel = FALSE
        , parallel_mem_limit = NULL
        , processing_strategy = c('sequential', 'all-in-one')[1]
        , minimal_output = FALSE
        , parallel_tmpdir = tempdir()
        , ...
	){
    # ARGUMENTS

    # start processing
	script_start <- Sys.time()

    # check parallelism startegy
    processing_strategy <- match.arg(processing_strategy, 
        choices = c('sequential', 'all-in-one', 'recursive'))

    # get current environment
    current_env <- environment()

    # if/else RECURSIVE
    if (processing_strategy != 'recursive') {

        # check input
        if (!missing(create_graphs) && create_graphs) {
            # TODO: eventually remove one of either arguments
            if (is.null(graphs_directory)) {
                stop('argument "create_graphs" is TRUE, but "graphs_directory" is not specified!')
            }
        }
        if (!missing(create_dailygraphs) && create_dailygraphs) {
            # TODO: eventually remove one of either arguments
            if (is.null(graphs_directory)) {
                stop('argument "create_dailygraphs" is TRUE, but "graphs_directory" is not specified!')
            }
        }
        # don't create figures if graphs_directory is missing
        if (missing(graphs_directory)) {
            create_dailygraphs <- create_graphs <- FALSE
        }
        if ((create_graphs || create_dailygraphs) && (
                !is.character(graphs_directory) || !dir.exists(graphs_directory)
                )) {
            stop('argument "graphs_directory": directory "', graphs_directory, 
                '" does not exist!')
        }

        # check data for nh3, co2 & h2o
        if (is.data.table(sonic_directory)) {
            if (any(grepl('nh3', names(sonic_directory)))) {
                ht_directory <- 'data provided with sonic'
            }
            if (any(grepl('co2', names(sonic_directory)))) {
                licor_directory <- 'data provided with sonic'
            }
            if (any(grepl('n2o', names(sonic_directory)))) {
                miro_directory <- 'data provided with sonic'
            }
        }

        # get scalars and fix missing instruments
        scalars <- sub('wx', '', grep('wx[^T].+', covariances, value = TRUE))
        # fix missing ht8700
        if (ht_null <- is.null(ht_directory)) {
            scalars <- grep('nh3', scalars, value = TRUE, invert = TRUE)
            covariances <- grep('nh3', covariances, value = TRUE, invert = TRUE)
        }
        # fix missing licor
        if (licor_null <- is.null(licor_directory)) {
            scalars <- grep('mmolm3', scalars, value = TRUE, invert = TRUE)
            covariances <- grep('mmolm3', covariances, value = TRUE, invert = TRUE)
        }
        # fix missing miro
        if (miro_null <- is.null(miro_directory)) {
            scalars <- grep('molfrac', scalars, value = TRUE, invert = TRUE)
            covariances <- grep('molfrac', covariances, value = TRUE, invert = TRUE)
        }

        # prepare variables etc.
        variables <- c('u', 'v', 'w', 'T', scalars)
        scalar_covariances <- grepl('wx[^T].+', covariances)
        names(scalar_covariances) <- covariances
        covariances_plotnames <- make.names(sub("x", "", covariances))
        names(covariances_plotnames) <- covariances
        covariances_variables <- strsplit(covariances, "x")
        if (is.null(z_ec) || !is.numeric(z_ec)) {
            stop('argument "z_ec" must be provided as numeric value (height in m a.g.l)!')
        }
        if (length(z_ec) > 1) {
            stop('argument "z_ec" has length > 1!')
        }
        if (is.null(z_canopy) || !is.numeric(z_canopy)) {
            stop('argument "z_canopy" must be provided as numeric value ',
                '(height of canopy in meters)!')
        }
        if (length(z_canopy) > 1) {
            stop('argument "z_canopy" has length > 1 which is not yet accepted!')
        }
        if (is.null(dev_north) || !is.numeric(dev_north)) {
            stop('argument "dev_north" must be provided as numeric value!')
        }
        if (is.null(declination)) {
            stop('argument "declination" must be provided!',
            ' -> it is also possible to provide a list with lon/lat entries...')
        } else if (is.list(declination) || length(declination) == 2) {
            if (!requireNamespace('oce')) {
                stop('R package "oce" must be installed when lat/lon is provided')
            }
            mag_dec <- \(x) oce::magneticField(declination[['lon']], 
                declination[['lat']], x)$declination
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
        rm(rot_args)
        detrending <- fix_defaults(detrending, variables)
        na_limits <- fix_defaults(na_limits, variables)
        limits_lower <- fix_defaults(limits_lower, variables)
        limits_upper <- fix_defaults(limits_upper, variables)
        despike <- fix_defaults(despike, variables)
        despike_baseline_width <- fix_defaults(despike_baseline_width, variables)
        despike_quantile <- fix_defaults(despike_quantile, variables)
        despike_stats_width <- fix_defaults(despike_stats_width, variables)
        despike_stats_multiply <- fix_defaults(despike_stats_multiply, variables)
        lag_fix <- fix_defaults(lag_fix, covariances)
        lag_dyn <- fix_defaults(lag_dyn, covariances)
        lag_dyn_method <- match.arg(lag_dyn_method)
        lag_dyn_calc_pw <- fix_defaults(lag_dyn_calc_pw, covariances)
        # check lag_dyn_method vs lag_dyn_calc_pw
        if (lag_dyn_method != 'raw-cov' && all(!lag_dyn_calc_pw)) {
            stop('argument "lag_dyn_method" != "raw-cov" but pre-whitening',
                ' method is switched off for all covariances!')
        }
        damping_reference <- fix_defaults(damping_reference, covariances[scalar_covariances])
        damping_lower <- fix_defaults(damping_lower, covariances[scalar_covariances])
        damping_upper <- fix_defaults(damping_upper, covariances[scalar_covariances])
        # sub intervals (4 to 8 are valid)
        if (!(is.numeric(subint_n) && (subint_n %% 1 == 0) && subint_n >= 4 &&
            subint_n <= 8)) {
            stop('argument "subint_n" must be an integer number between 3 and 10')
        }
        subintervals <- subintervals || subint_return
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
        # ts_vars <- union(names(plot_timeseries), ts_vars)
        plot_timeseries <- fix_defaults(plot_timeseries, ts_vars)
        # only take TRUE
        plot_timeseries <- plot_timeseries[plot_timeseries]
        plotting_var_colors <- fix_defaults(plotting_var_colors, ts_vars)[names(plot_timeseries)]
        plotting_var_units <- fix_defaults(plotting_var_units, ts_vars)[names(plot_timeseries)]
        plotting_covar_units <- fix_defaults(plotting_covar_units, covariances)
        plotting_covar_colors <- fix_defaults(plotting_covar_colors, covariances)
        model_color_names <- names(formals(process_ec_fluxes)$model_colors)
        if (is.null(names(model_colors)) && 
            length(model_colors) == length(model_color_names)) {
            names(model_colors) <- model_color_names
        } else {
            model_colors <- fix_defaults(model_colors, model_color_names)
        }

        lim_range <- rbind(lower = limits_lower, upper = limits_upper)
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
            # remove duplicated
            if (sonic_old_format) {
                # remove duplicated old (no ending)
                sonic_files <- c(
                    setdiff(sonic_files[!i_gz], sub('\\.gz$', '', sonic_files[i_gz])),
                    sonic_files[i_gz]
                )
            } else {
                # remove duplicated non-py
                sonic_files <- c(
                    setdiff(sonic_files[!i_gz], sub('\\.gz$', '.csv', sonic_files[i_gz])),
                    sonic_files[i_gz]
                )
            }
            # get date
            sonic_dates <- sub('^((py_)?fnf_0\\d_sonic_|data_sonic-._)', '', sonic_files)
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
                ht_dates <- sub('^(py_)?fnf_0\\d_ht8700_', '', ht_files)
                # sort by date
                ht_files <- ht_files[order(ht_dates)]
                # get start
                if (ht_old_format) {
                    ht_pattern <- c('.*_(\\d{8})_(\\d{6}).*', '\\1\\2', 'Etc/GMT-1', '\\1000000')
                } else {
                    ht_pattern <- c('.*_(\\d{4})_(\\d{2})_(\\d{2}).csv', '\\1\\2\\3000000', 'UTC')
                    ht_pattern[4] <- ht_pattern[2]
                }
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
                licor_dates <- sub('^(py_)?fnf_0\\d_licor_', '', licor_files)
                # sort by date
                licor_files <- licor_files[order(licor_dates)]
                # get start & end
                licor_pattern <- c('.*_(\\d{4})_(\\d{2})_(\\d{2}).csv', '\\1\\2\\3000000',
                    'UTC')
            } else if (licor_with_sonic) {
                cat('LI-7500: raw data provided with sonic input...\n')
            }
        # } else {
            # check if ht available
            # if (!ht_provided) {
            #     stop('neither ht8700 nor licor data or directory has been provided -> cannot process fluxes without concentration data!')
            # }
        }

        # optional miro data
        miro_with_sonic <- FALSE
        if (miro_provided <- !is.null(miro_directory)) {
            if (miro_has_data <- inherits(miro_directory, 'data.table')) {
                cat('MIRO: raw data provided in input...\n')
                miro_files <- NULL
                setDT(miro_directory)
                # set time zone to UTC
                miro_directory[, Time := with_tz(Time, 'UTC')]
            } else if (!(miro_with_sonic <- miro_directory == 'data provided with sonic')) {
                cat('MIRO: data provided...\n')
                miro_files <- dir(miro_directory, pattern = '^(py_)?fnf_.*_miro_.*')
                if (length(miro_files) == 0) {
                    stop('No miro data available in "', miro_directory, '"!')
                }
                # prioritize py_ files
                # get py_
                i_py <- grepl('^py_', miro_files)
                # remove duplicated non-py
                miro_files <- c(
                    setdiff(miro_files[!i_py], sub('^py_', '', miro_files[i_py])),
                    miro_files[i_py]
                )
                # prioritize gz files
                i_gz <- grepl('\\.gz$', miro_files)
                # remove duplicated non-py
                miro_files <- c(
                    setdiff(miro_files[!i_gz], sub('\\.gz$', '.csv', miro_files[i_gz])),
                    miro_files[i_gz]
                )
                # get date
                miro_dates <- sub('^(py_)?fnf_0\\d_miro_', '', miro_files)
                # sort by date
                miro_files <- miro_files[order(miro_dates)]
                # get start & end
                miro_pattern <- c('.*_(\\d{4})_(\\d{2})_(\\d{2}).csv', '\\1\\2\\3000000',
                    'UTC')
            } else if (miro_with_sonic) {
                cat('MIRO: raw data provided with sonic input...\n')
            }
        # } else {
            # check if ht available
            # if (!ht_provided) {
            #     stop('neither ht8700 nor miro data or directory has been provided -> cannot process fluxes without concentration data!')
            # }
        }
        cat("************************************************************\n")

        # parse time diff
        avg_secs <- parse_time_diff(avg_period)

        # check gamma_time_window
        if (max(gamma_time_window) * 60 >= avg_secs) {
            warning('"gamma_time_window" is too large compared to "avg_period"',
                ' -> changing to smaller values!')
            gamma_time_window <- c(5, 10) / 30 * avg_secs / 60
        }

        # check times
        if (is.null(start_time) && is.null(end_time)) {
            start_time <- 'first'
            end_time <- 'last'
        } 

        # fix start_time (only use sonic)
        if (length(start_time) == 1 && is.character(start_time) && 
            start_time == 'first') {
            start_time <- with_tz(start_sonic, tz_user)
        } else if (!is.null(start_time)) {
            start_time <- parse_date_time3(start_time, tz = tz_user)
        }

        # fix end_time (only use sonic)
        if (is.null(end_time)) {
            end_time <- start_time + avg_secs
        } else if (length(end_time) == 1 && is.character(end_time) && 
            end_time == 'last') {
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
        # get full vector if not yet provided
        if (length(start_time) == 1) {
            # "fix" end_time not included
            start_time <- seq(start_time, end_time - 1e-4, by = avg_secs)
            end_time <- start_time + avg_secs
        }

        # prepare dates
        dates_utc <- unique(lubridate::date(c(start_time, end_time - 1e-4)))
        dates_formatted <- gsub('-', '_', dates_utc, fixed = TRUE)

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
            create_dailygraphs <- create_graphs <- FALSE
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

        # setup parallelism
        if (run_parallel <- (!is.numeric(ncores) || isFALSE(ncores == 1))) {
            if (inherits(ncores, 'cluster')) {
                # copy cluster
                cl <- ncores
            } else {
                cat('\nSetting up parallelism...\n')
                # start cluster
                cl <- do.call(.makePSOCKcluster,
                    c(list(names = ncores, memory_limit = parallel_mem_limit),
                        psock_args))
                # stop cluster on exit
                on.exit(parallel::stopCluster(cl))
                # set data.table threads to 1 on slaves
                parallel::clusterEvalQ(cl, data.table::setDTthreads(1L))
                # cat('-> exporting R objects...')
                # # export ec functions + all reading functions + cpp codes
                # parallel::clusterExport(cl, c(.ec_evaluation, .read_sonic,
                #         .ht8700_functions, .licor_functions))
                # cat(' done.\n')
            }
            cat('\n=> Processing in parallel using', length(cl), 'cores.\n\n')
        }

        # get temporary file name base
        tf <- tempfile(tmpdir = parallel_tmpdir)
    # if/else RECURSIVE else
    } else {
        # assign dots to current env
        nms <- ...names()
        for (i in seq_len(...length())) {
            assign(nms[i], ...elt(i))
        }
        # read cobj dots
        dots <- qs2::qs_read(tf_cobj)
        # assign dots to current env
        for (what in names(dots)) {
            assign(what, dots[[what]])
            dots[[what]] <- NULL
        }
        rm(dots)
        mag_dec <- function(x) {}
        body(mag_dec) <- parse(text = mag_dec_body)
        # read sonic_directory
        sonic_directory <- qs2::qd_read(tf_sonic)
        ht_directory <- qs2::qd_read(tf_ht)
        licor_directory <- qs2::qd_read(tf_licor)
        # check times (for parallel calculation with less intervals than cores)
        if (length(start_time) == 0 || length(end_time) == 0) {
            # return NULL if length is zero
            return(NULL)
        }

    # if/else RECURSIVE end
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # check all-in-one vs. sequential
    if (length(start_time) > 1 && processing_strategy == 'sequential') {
        # SEQUENTIAL
        # get current objects as string
        #   -> strip off objects which need changes
        cobj <- setdiff(ls(envir = current_env), 
            c('start_time', 'end_time', 'processing_strategy', 'dates_utc', 
                'dates_formatted', 'as_ibts', 'current_env', 'mag_dec',
                'ncores', 'run_parallel', 'cl', 
                'sonic_directory', 'ht_directory', 'licor_directory'))
        # get start/end time dates
        st_dates <- as.Date(start_time)
        et_dates <- as.Date(end_time - 1e-4)
        # create tempfile paths
        tf_cobj <- paste0(tf, '-cobj.qs2')
        tf_resid <- paste0(tf, '-resid.qdata')
        tf_sonic <- paste0(tf, '-sonic.qdata')
        tf_ht <- paste0(tf, '-ht.qdata')
        tf_licor <- paste0(tf, '-licor.qdata')
        # save cobj exports as qs2 & qdata
        dots <- mget(cobj, envir = current_env)
        # fix saving mag_dec function
        mag_dec_body <- deparse1(body(mag_dec))
        #   -> fix parallel settings & mag_dec
        # set ncores to 1 & run_parallel to FALSE
        dots <- c(dots, list(mag_dec_body = mag_dec_body, ncores = 1, run_parallel = FALSE))
        qs2::qs_save(dots, tf_cobj, compress_level = -7)
        rm(dots)
        qs2::qd_save(sonic_directory, tf_sonic, warn_unsupported_types = FALSE, 
            compress_level = -7)
        qs2::qd_save(ht_directory, tf_ht, warn_unsupported_types = FALSE,
            compress_level = -7)
        qs2::qd_save(licor_directory, tf_licor, warn_unsupported_types = FALSE,
            compress_level = -7)
        # TODO: clean up tempfiles later!!!
        rm(sonic_directory, ht_directory, licor_directory)
        for (i in 1:10) gc()
        # copy files to # of workers ??
        # save on ssd (cruncher)
        # ..!
        # call main function recursively
        if (run_parallel) {
            # distribute intervals
            # get number of cluster
            nc <- length(cl)
            # get unique dates
            st_udates <- unique(st_dates)
            # get runlength
            nrle <- rle(as.numeric(st_dates))
            # number of days > 20h (40 intervals)
            full_day <- which(nrle$lengths > 40)
            nd <- length(full_day)
            incomplete_day <- which(nrle$lengths <= 40)
            # split full days
            if (nd >= nc) {
                # get number of full days to split
                ifd <- nc * floor(nd / nc)
                # check before first full day
                if (ifd != nd && length(incomplete_day) > 0 && 
                    any(incomplete_day < full_day[1])) {
                    ind_fd <- full_day[length(full_day) + 1 - rev(seq_len(ifd))]
                } else {
                    ind_fd <- full_day[seq_len(ifd)]
                }
                # extend incomplete days
                incomplete_day <- sort(c(incomplete_day, setdiff(full_day, ind_fd)))
                # split days
                ind_split <- lapply(ind_fd, \(x) which(st_dates == st_udates[x]))
            } else {
                # prepare index list
                ind_split <- list()
                # add full_day to incomplete_day
                incomplete_day <- sort(c(incomplete_day, full_day))
            }
            # split incomplete days
            if (length(incomplete_day)) {
                # get intervals
                ints <- unlist(lapply(incomplete_day, \(x) which(st_dates == st_udates[x])))
                # split incomplete days
                ind_split <- c(
                    ind_split,
                    parallel::clusterSplit(cl, ints)
                )
            }
            # save residual exports
            qd_save(list(
                    st_dates = st_dates, 
                    et_dates = et_dates, 
                    start_time = start_time, 
                    end_time = end_time
                    ), tf_resid, warn_unsupported_types = FALSE)
            if (isTRUE(debug_parallel)) {
                # be verbose
                cat('~~~\nDebugging main function in parallel mode...\n')
                # debug parallelism
                out_list_paths <- list()
                for (ispl in ind_split) {
                    out_list_paths <- c(out_list_paths,
                        .pef_wrapper(ispl, tf_cobj = tf_cobj, 
                            tf_resid = tf_resid, tf_sonic = tf_sonic, 
                            tf_ht = tf_ht, tf_licor = tf_licor
                        )
                    )
                }
            } else {
                # be verbose
                cat('~~~\nRunning main function in parallel...\n')
                # # call main function in parallel
                out_list_paths <- .clusterApplyLB(cl, ind_split, 
                    .pef_wrapper, tf_cobj = tf_cobj, tf_resid = tf_resid, 
                    tf_sonic = tf_sonic, tf_ht = tf_ht, tf_licor = tf_licor
                )
            }
            # remove files
            unlink(tf_resid)
            # loop over list
            out_list <- vector('list', length(out_list_paths))
            for (i in seq_along(out_list_paths)) {
                # check on try errors
                if (inherits(out_list_paths[[i]], 'try-error')) {
                    # add times
                    cat('try errors in parallel calls!',
                        'Returning worker results list for checking..\n')
                    return(out_list_paths)
                }
                # read file
                if (!is.null(out_list_paths[[i]])) {
                    out_list[[i]] <- qs2::qd_read(out_list_paths[[i]])
                    # delete file
                    unlink(out_list_paths[[i]])
                }
            }
        } else {
            # loop over dates
            out_list <- lapply(
                unique(st_dates), 
                \(udate) {
                    # get indices for date
                    ind <- st_dates == udate
                    # get utc dates
                    utc_dates <- unique(c(udate, et_dates[ind]))
                    formatted_dates <- gsub('-', '_', utc_dates, fixed = TRUE)
                    process_ec_fluxes(
                        dates_utc = utc_dates,
                        dates_formatted = formatted_dates,
                        start_time = start_time[ind],
                        end_time = end_time[ind],
                        as_ibts = FALSE,
                        processing_strategy = 'recursive',
                        tf_cobj = tf_cobj,
                        tf_sonic = tf_sonic,
                        tf_ht = tf_ht,
                        tf_licor = tf_licor
                    )
                }
            )
        }
        # remove files
        unlink(tf_cobj)
        unlink(tf_ht)
        unlink(tf_sonic)
        unlink(tf_licor)
        # fix subints
        if (subint_return) {
            subints <- rbindlist(lapply(out_list, attr, 'subintervals'),
                fill = TRUE, use.names = TRUE)
        }
        # rbind output list
        results <- rbindlist(out_list, fill = TRUE)
        if (subint_return) {
            setattr(results, 'subintervals', subints)
        }
        if (!(create_dailygraphs || ogives_return)) {
            rm(out_list)
            for (i in 1:10) gc()
        }
        # SEQUENTIAL END
    } else {
        # ALL-IN-ONE/RECURSIVE

        # read raw data:
        # ------------------------------------------------------------------------------

        # get sonic data
        if (sonic_has_data) {
            cat('~~~\nSubsetting sonic data - ')
            # subset provided sonic data
            sonic_raw <- sonic_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
            rm(sonic_directory)
            for (i in 1:10) gc()
        } else {
            cat('~~~\nReading sonic files\n')
            # select files
            if (sonic_old_format) {
                sonic_selected <- sonic_files[
                    sub('data_sonic-._(\\d{4})(\\d{2})(\\d{2})_.*', '\\1_\\2_\\3', sonic_files) %in% dates_formatted
                    ]
            } else {
                sonic_selected <- sonic_files[
                    sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', sonic_files) %in% dates_formatted
                    ]
            }
            # read sonic files
            if (length(sonic_selected)) {
                if (run_parallel) {
                    sonic_raw <- rbindlist(
                        .clusterApplyLB(
                            cl,
                            file.path(sonic_directory, sonic_selected), 
                            read_sonic
                        )
                    )
                } else {
                    sonic_raw <- rbindlist(lapply(
                            file.path(sonic_directory, sonic_selected), 
                            \(x) {
                                cat('\t')
                                read_sonic(x)
                            }
                    ))
                }
            } else {
                sonic_raw <- NULL
            }
            # check sonic data
            if (is.null(sonic_raw) || nrow(sonic_raw) == 0) {
                cat('No sonic data available within given time range...\n')
                return(NULL)
            }
            # get time zone
            sonic_tz <- sonic_raw[, tz(Time)]
            if (sonic_tz != 'UTC') {
                # convert to UTC
                sonic_raw[, Time := with_tz(Time, 'UTC')]
                # fix start_time & end_time
                start_time <- with_tz(force_tz(start_time, sonic_tz), 'UTC')
                end_time <- with_tz(force_tz(end_time, sonic_tz), 'UTC')
            }
        }
        cat('done\n')

        # get ht8700 data
        if (ht_provided && ht_has_data) {
            cat('Subsetting HT8700 data - ')
            # subset
            ht <- ht_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
            rm(ht_directory)
            for (i in 1:10) gc()
            cat('done\n')
        } else if (ht_provided && !ht_with_sonic) {
            cat('Reading HT8700 files\n')
            # select files
            ht_selected <- ht_files[
                sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', ht_files) %in% dates_formatted
                ]
            # read ht files
            if (length(ht_selected)) {
                if (run_parallel) {
                    ht <- rbindlist(
                        .clusterApplyLB(
                            cl,
                            file.path(ht_directory, ht_selected), 
                            read_ht8700
                        )
                    )
                } else {
                    ht <- rbindlist(lapply(
                            file.path(ht_directory, ht_selected), 
                            \(x) {
                                cat('\t')
                                read_ht8700(x)
                            }
                    ))
                }
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
        if (!is.null(ht) && nrow(ht) == 0) ht <- NULL
        if (ht_provided && is.null(ht)) {
            cat('-> No HT8700 data within time range.\n')
        }

        # get licor data
        if (licor_provided && licor_has_data) {
            cat('Subsetting LI-7500 data - ')
            # copy
            licor <- licor_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
            rm(licor_directory)
            for (i in 1:10) gc()
            cat('done\n')
        } else if (licor_provided && !licor_with_sonic) {
            cat('Reading LI-7500 files\n')
            # select files
            licor_selected <- licor_files[
                sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', licor_files) %in% dates_formatted
                ]
            if (length(licor_selected)) {
                # read new licor files
                if (run_parallel) {
                    licor <- rbindlist(
                        .clusterApplyLB(
                            cl,
                            file.path(licor_directory, licor_selected), 
                            read_licor
                        )
                    )
                } else {
                    licor <- rbindlist(lapply(
                            file.path(licor_directory, licor_selected), 
                            \(x) {
                                cat('\tFile:', x, '- ')
                                out <- read_licor(x)
                                cat('done\n')
                                out
                            }
                    ))
                }
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
        if (!is.null(licor) && nrow(licor) == 0) licor <- NULL
        if (licor_provided && is.null(licor)) {
            cat('-> No LI-7500 data within time range.\n')
        }

        # get miro data
        if (miro_provided && miro_has_data) {
            cat('Subsetting MIRO data - ')
            # copy
            miro <- miro_directory[Time >= start_time[1] & Time < tail(end_time, 1), ]
            rm(miro_directory)
            for (i in 1:10) gc()
            cat('done\n')
        } else if (miro_provided && !miro_with_sonic) {
            cat('Reading MIRO files\n')
            # select files
            miro_selected <- miro_files[
                sub('.*_(\\d{4}_\\d{2}_\\d{2})\\..*', '\\1', miro_files) %in% dates_formatted
                ]
            if (length(miro_selected)) {
                # read new miro files
                if (run_parallel) {
                    miro <- rbindlist(
                        .clusterApplyLB(
                            cl,
                            file.path(miro_directory, miro_selected), 
                            read_miro
                        )
                    )
                } else {
                    miro <- rbindlist(lapply(
                            file.path(miro_directory, miro_selected), 
                            \(x) {
                                # cat('\tFile:', x, '- ')
                                out <- read_miro(x)
                                # cat('done\n')
                                out
                            }
                    ))
                }
            } else {
                miro <- NULL
            }
            # no UTC conversion needed, since this case is excluded up to now
            # # convert to UTC
            # miro[, Time := with_tz(Time, 'UTC')]
        } else {
            # no miro data
            miro <- NULL
        }
        # check miro
        if (!is.null(miro) && nrow(miro) == 0) miro <- NULL
        if (miro_provided && is.null(miro)) {
            cat('-> No MIRO data within time range.\n')
        }

        # get Hz/frequency etc.
        rec_Hz <- sonic_raw[, {
            d_t <- diff(as.numeric(Time))
            round(1 / median(d_t, na.rm = TRUE), -1)
        }]
        if (rec_Hz < 10) {
            stop('Frequency is lower than 10 Hz! Aborting script!')
        }

        # fix time to exactly xx Hz
        t_basis <- unlist(mapply(\(x_st, x_et) {
                seq(x_st, x_et - 1 / rec_Hz, by = 1 / rec_Hz)
            }, x_st = start_time, x_et = end_time,
            SIMPLIFY = FALSE
        ))
        t0 <- t_basis[1]
        # get matching indices
        t_indices <- sonic_raw[, {
            I(match_times(t_basis - t0, as.numeric(Time, units = 'secs') - t0, 
                rec_Hz * 1.1
            ))
        }]
        rm(t0)
        # sonic including all times
        full_sonic <- data.table(
            Time = .POSIXct(t_basis, tz = 'UTC'), Hz = rec_Hz,
            u = NA_real_, v = NA_real_, w = NA_real_, T = NA_real_,
            sonic = sonic_raw[, sonic[1]]
        )
        rm(t_basis)
        # fill sonic data
        full_sonic[t_indices[[1]], c('u', 'v', 'w', 'T') := 
            sonic_raw[t_indices[[2]], .SD, .SDcols = c('u', 'v', 'w', 'T')]]
        rm(sonic_raw, t_indices)
        for (i in 1:10) gc()

        cat('Merging files - ')
        daily_data <- merge_data(full_sonic, ht, licor, miro)
        rm(full_sonic, ht, licor, miro)
        for (i in 1:10) gc()
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
        daily_data[, bin := getIntervals(Time, start_time, end_time)]

        # TODO: add measures on interval basis
        n_period <- avg_secs * rec_Hz
        n_threshold <- thresh_period * n_period
        freq <- rec_Hz * seq(floor(n_period / 2)) / floor(n_period / 2)

        # --------------------- read fix lags from lag_lookuptable ---------------------
        # TODO: only if necessary/wanted

        # fix and dyn lags:                                      
        # ------------------------------------------------------------------------------ 
        # check fix lag (convert to list of functions)
        # possible ways of providing fix lag:
        # list(uxw = 0, wxnh3_ugm3 = c('[-45, 135)' = -0.2, 0.4),
            # wxh2o_mmolm3 = list('[-45, 60)' = 0.2, '[60, 180)' = 0.3), # default from argument lag_fix
            # wxh2o_mmolm3 = list('[-45, 60)' = 0.2, '[50, 180)' = 0.3, 0), # should fail
            # if missing => defaults
            # wxco2_mmolm3 = list('[-45, 60)' = -0.2, '[60, 180)' = 0.3, '[180, -45)' = 0)
        # also allowed: 60 - 90 or 60/90 or 60//90
        input_lag_functions <- sapply(names(lag_fix), \(nx) {
            x <- lag_fix[[nx]]
            if (length(x) == 1) {
                fun_out <- function(wd, dyn = FALSE) {
                    if (dyn) {
                        c(
                            lower = lag_value - dlag,
                            upper = lag_value + dlag
                        )
                    } else {
                        lag_value
                    }
                }
                efun <- environment(fun_out)
                efun$lag_value <- round(x * rec_Hz)
            } else {
                # get ranges
                rgs <- names(x)
                # is there a residual lag time?
                is_resid <- rgs == ''
                if (sum(is_resid) > 1) {
                    stop('fix lag by wind sector has more than one default (unnamed) lag time!')
                } else if (any(is_resid)) {
                    # convert to numeric
                    resid_value <- as.numeric(x[is_resid])
                } else {
                    # get default from arguments
                    resid_value <- eval(formals(process_ec_fluxes)$lag_fix)[nx]
                }
                # subset
                rgs_sub <- rgs[!is_resid]
                x_sub <- as.numeric(x[!is_resid])
                # check ranges validity
                pdig <- '(-|+)?\\d+[.]?\\d*'
                pattern <- paste0(
                    '^\\s*([[]|[(])?\\s*',
                    pdig,
                    '\\s*([,]|[/][/]?|\\s[-]\\s)\\s*',
                    pdig,
                    '\\s*([)]|[]])?\\s*$'
                )
                is_valid <- grepl(pattern, rgs_sub)
                if (any(!is_valid)) {
                    stop('fix lag named entries: ', paste(paste0('"', rgs_sub[!is_valid], '"'), collapse = '; '), ' are not valid!')
                }
                # check "closed"
                left_closed <- grepl('[', rgs_sub, fixed = TRUE)
                right_closed <- grepl(']', rgs_sub, fixed = TRUE)
                # fix missing closed
                left_closed[!left_closed & !right_closed] <- TRUE
                # split
                xs <- strsplit(trimws(gsub('[[]|[(]|[)]|[]]|,|/|(\\s[-]\\s)', ' ', rgs_sub)), split = '\\s+')
                # convert to numeric
                xs <- lapply(xs, as.numeric)
                # get lower & upper
                lower <- sapply(xs, '[[', 1) %% 360
                upper <- sapply(xs, '[[', 2) %% 360
                value <- as.numeric(x[!is_resid]) * rec_Hz
                # fill up
                out <- data.frame(lower = numeric(0), upper = numeric(0), 
                    value = numeric(0), right_closed = logical(0), 
                    left_closed = logical(0), pass_360 = logical(0))
                for (current in seq_along(lower)) {
                    prev <- if (current == 1) length(lower) else current - 1
                    # check overlapping & ranges coverage
                    d <- lower[current] - upper[prev]
                    if (d != 0) {
                        # fill with default
                        out <- rbind(out,
                            data.frame(
                                lower = upper[prev],
                                upper = lower[current],
                                value = resid_value,
                                right_closed = !left_closed[current],
                                left_closed = !right_closed[prev],
                                pass_360 = lower[current] < upper[prev]
                            )
                        )
                    }
                    # add current
                    out <- rbind(out,
                        data.frame(
                            lower = lower[current],
                            upper = upper[current],
                            value = x_sub[current],
                            right_closed = right_closed[current],
                            left_closed = left_closed[current],
                            pass_360 = upper[current] < lower[current]
                        )
                    )
                }
                # check overlapping sectors
                degs <- sum((out$upper - out$lower) %% 360)
                if (degs > 360) {
                    stop('overlapping wind sectors defined in fix lag!')
                }
                # change secs -> steps
                out$value <- out$value * rec_Hz
                # build function
                fun_out <- function(wd, dyn = FALSE) {
                    mwd <- wd %% 360
                    # check left closed?
                    is_left <- mwd == tab$lower & tab$left_closed
                    is_right <- mwd == tab$upper & tab$right_closed
                    is_normal <- mwd > tab$lower & mwd < tab$upper & !tab$pass_360
                    is_p360 <- (mwd > tab$lower | mwd < tab$upper) & tab$pass_360
                    out <- tab$value[is_left | is_right | is_normal | is_p360]
                    if (dyn) {
                        c(
                            lower = out - dlag,
                            upper = out + dlag
                        )
                    } else {
                        out
                    }
                }
                # assign table to function environment
                efun <- environment(fun_out)
                efun$tab <- out
            }
            # add dlag
            efun$dlag <- lag_dyn[[nx]] * rec_Hz
            # return function
            fun_out
        }, simplify = FALSE)

        # input_fix_lag <- round(lag_fix * rec_Hz)
        # input_dyn_lag <- rbind(
        #     lower = round((lag_fix - lag_dyn) * rec_Hz)
        #     , upper = round((lag_fix + lag_dyn) * rec_Hz)
        # )

        # be verbose
        cat("\n~~~\nCalculation will include",
            daily_data[, uniqueN(bin)] , "intervals between", 
            format(start_time[1], format = "%Y-%m-%d %H:%M"), "and", 
            format(tail(end_time, 1), 
                format = "%Y-%m-%d %H:%M", usetz = TRUE), 
            "on a", avg_secs / 60, "minute basis\n~~~\n\n"
        )

        # check HT8700 quality: alarm codes & OSS
        if (ht_provided) {
            cat('~~~\nChecking HT-8700 OSS and alarm codes - ')
            # add alarm code
            if (!('ht_alarm_code' %in% names(daily_data))) {
                daily_data[, ht_alarm_code := get_alarms(.SD)]
            }
            # create regex pattern
            na_alarm_pattern <- paste(paste0('\\b', na_alarm_code, '\\b'),
                collapse = '|')
            # check alarms and set nh3 NA
            nh3_vars <- grep('nh3', names(daily_data), value = TRUE)
            na0 <- daily_data[, sum(is.na(get(nh3_vars[1])))]
            daily_data[grepl(na_alarm_pattern, ht_alarm_code), (nh3_vars) := NA_real_]
            na1 <- daily_data[, sum(is.na(get(nh3_vars[1])))]
            # check oss
            daily_data[ht_oss < oss_threshold, (nh3_vars) := NA_real_]
            cat('done\nBad alarms:', na1 - na0, '\nValues below OSS threshold:', 
                daily_data[, sum(ht_oss < oss_threshold, na.rm = TRUE)], '\n')
        }

        # check licor ss
        if (licor_provided) {
            cat('~~~\nChecking LI-7500 signal strength - ')
            # get variables
            li_vars <- grep('^(co2|h2o)_', names(daily_data), value = TRUE)
            # check li_co2ss
            daily_data[li_co2ss < co2ss_threshold, (li_vars) := NA_real_]
            cat('done\nValues below "co2ss" threshold:', 
                daily_data[, sum(li_co2ss < co2ss_threshold, na.rm = TRUE)], '\n')
        }

        # raw data quality control I, despiking
        # --------------------------------------------------------------------------

        # despiking procedure
        if (any(despike)) {

            # backup original data for plotting
            daily_data[, paste0(names(despike), '_original') := copy(.SD),
                .SDcols = names(despike)]

            cat("~~~\nDespiking raw time series...\n")
            # routine
            for (s in names(despike)[despike]) {
                # s <- 'co2_mmolm3'
                despike_timeseries(daily_data, scalar = s,
                    filter_width = despike_baseline_width[[s]],
                    qval = despike_quantile[[s]],
                    qwidth = despike_stats_width[[s]],
                    qmult = despike_stats_multiply[[s]]
                )
            }
            # plottin is done in intervals

            # backup despiked data for plotting
            daily_data[, paste0(names(despike), '_despiked') := copy(.SD),
                .SDcols = names(despike)]

        }

        # raw data quality control II, i.e. hard limits = physical range
        # --------------------------------------------------------------------------
        cat("~~~\nchecking NA-values and hard limits...\n")
        hl_vars <- names(na_limits)[na_limits]
        if (any(!(hl_vars %in% colnames(lim_range)))) {
            hl_missing <- !(hl_vars %in% colnames(lim_range))
            stop(
                'hard limits are missing for variables: ',
                paste(names(na_limits)[hl_missing], collapse = ', ')
            )
        }
        check_limits(daily_data, lim_range[, hl_vars], na_limits_window, 
            na_limits_method, rec_Hz, n_threshold)


        # define coordinate system
        if (daily_data[, sonic[1] == 'HS']) {
            coord.system <-  'HS-Wauwilermoos'
        } else {
            coord.system <-  'Windmaster'
        }

        # declination
        current_declination <- daily_data[, mag_dec(Time[1]), by = bin][, 
            setNames(V1, bin)]
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
        daily_data[, WD := (WD + d_north[.BY[[1]]]) %% 360, 
            by = as.character(bin)]

        # get subinterval bins
        if (subintervals) {
            # add subint bins
            avg_sub <- avg_secs / subint_n
            i_sub <- 1:subint_n
            daily_data[, subint := {
                # get st & et
                st_int <- start_time[.BY[[1]]]
                et_int <- end_time[.BY[[1]]]
                # get sub times
                st_sub <- seq(st_int, et_int, by = avg_sub)[i_sub]
                et_sub <- st_sub + avg_sub
                # split into sub-intervals
                sub_ints <- getIntervals(Time, st_sub, et_sub)
                paste(.BY[[1]], sub_ints, sep = '-')
            }, by = bin]
            cat('~~~\nrotate subinterval data...\n')
            # rotate raw data
            daily_data[, paste0('subint_', 
                c("WD", "phi", "urot", "vrot", "wrot")) := 
                rotate_twoaxis(u, v, w,
                    phi = if (rotation_method[1] %in% "two axis") {
                        # two-axis rotation
                        rotation_args$phi 
                    } else if (is.na(alpha[1])) {
                        # planar fit failed
                        NULL 
                    } else {
                        # planar fit successful
                        0
                    }, c.system = coord.system
                ), by = subint]
            # correct for sonic north deviation
            daily_data[, subint_WD := (subint_WD + d_north[.BY[[1]]]) %% 360, by = as.character(bin)]
            # backup original T_sonic + scalars
            daily_data[, paste0('subint_', c('T', scalars)) := .SD,
                .SDcols = c('T', scalars)]
        }

        # get relevant environment objects
        env_obj <- setdiff(ls(envir = current_env), c(
            'cl', 'ncores', 'run_parallel', 'daily_data', 
            'current_env'#, 'rotation_args', 'mag_dec'
        ))
        eobj <- mget(env_obj, envir = current_env)

        # loop over intervals: call MAIN function
        if (run_parallel) {
            # save env objects
            tf_env <- paste0(tf, '-env.qs2')
            qs2::qs_save(eobj, tf_env)
            # save daily_data
            dd <- daily_data[, {
                tf_sd <- paste0(tf, '-SD-', .BY[[1]], '.qdata')
                qs2::qd_save(cbind(.SD, bin = .BY[[1]]), tf_sd, 
                    warn_unsupported_types = FALSE)
                .(tf_sd)
            }, by = bin]
            if (isTRUE(debug_parallel)) {
                cat('~~~\ndebugging fluxes in parallel mode...\n')
                # debug parallelism
                out_list <- list()
                for (tf_sd in dd[, tf_sd]) {
                    out_list <- c(out_list, .wrapper_main(tf_sd, tf_env))
                }
            } else {
                cat('~~~\nprocessing fluxes in parallel...\n')
                # run main function
                out_list <- .clusterApplyLB(cl, dd[, tf_sd], .wrapper_main, 
                    tf_env = tf_env)
            }
            results <- rbindlist(lapply(out_list, \(x) {
                # read results
                out <- qs2::qd_read(x)
                # delete temporary file
                unlink(x)
                # return
                out
            }), fill = TRUE)
            # delete temporary files
            unlink(tf_env)
        } else {
            cat('~~~\nprocessing fluxes sequentially...')
            results <- .ec_main(daily_data, eobj)
        }
        # ALL-IN-ONE/RECURSIVE end
    }

    # check if no results
    if (is.null(results) || nrow(results) == 0) {
        results <- NULL
    } else {
        # fix ogives out
        if (create_dailygraphs || ogives_return) {
            if (processing_strategy == 'sequential') {
                # get ogives etc.
                Covars_Out <- unlist(lapply(out_list, attr, 'covars'), 
                    recursive = FALSE)
                Cospec_fix_Out <- unlist(lapply(out_list, attr, 'cospec_fix'), 
                    recursive = FALSE)
                Cospec_dyn_Out <- unlist(lapply(out_list, attr, 'cospec_dyn'), 
                    recursive = FALSE)
                Ogive_fix_Out <- unlist(lapply(out_list, attr, 'ogv_fix'), 
                    recursive = FALSE)
                Ogive_dyn_Out <- unlist(lapply(out_list, attr, 'ogv_dyn'), 
                    recursive = FALSE)
                if (ogives_return) {
                    # assign to output
                    results <- structure(
                        results, 
                        covars = Covars_Out,
                        cospec_fix = Cospec_fix_Out, 
                        cospec_dyn = Cospec_dyn_Out,
                        ogv_fix = Ogive_fix_Out, 
                        ogv_dyn = Ogive_dyn_Out
                    )
                }
            } else if (processing_strategy == 'all-in-one' && create_dailygraphs) {
                # get ogives etc. for plotting
                Covars_Out <- attr(results, 'covars')
                Cospec_fix_Out <- attr(results, 'cospec_fix')
                Cospec_dyn_Out <- attr(results, 'cospec_dyn')
                Ogive_fix_Out <- attr(results, 'ogv_fix')
                Ogive_dyn_Out <- attr(results, 'ogv_dyn')
            }
        }
        # create daily graphs
        if (processing_strategy != 'recursive' && create_dailygraphs) {
            if (avg_period != '30mins') {
                stop('daily figures not yet implemented for avg_period != "30mins"')
            }
            cat('~~~\ncreating daily graphs...')
            if (!dir.exists(path_folder)) {
                dir.create(path_folder, recursive = FALSE)
            }
            # get indices split by date
            idate <- results[, I(split(format(st, '%Y-%m-%d %H:%M'), date(st)))]
            # get freq (this might fail!)
            freq <- Cospec_dyn_Out[[1]]$freq
            # get covariance range
            sec_hr <- 50
            rec_Hz <- unique(sapply(Covars_Out, \(x) x$Hz))
            n_period <- unique(sapply(Covars_Out, \(x) x$N))
            x_cov <- seq(-sec_hr * rec_Hz, sec_hr * rec_Hz) / rec_Hz
            cov_sub <- n_period / 2 + seq(-sec_hr * rec_Hz, sec_hr * rec_Hz)
            # set resolution
            res <- 200
            # loop over covariances
            for (flux in covariances) {
                # flux <- covariances[1]
                # loop over dates
                for (day in names(idate)) {
                    # day <- names(idate)[1]
                    # create covariances figure
                    fn_covars <- paste0(path_folder, '/dailygraphs-covariances-', 
                        flux, '-', day, '.png')
                    png(fn_covars, width = 14, height = 14, units = 'in', res = res, ...)
                    # x11(width = 14, height = 14)
                    par(mfrow = c(8, 6), mar = c(2, 4, 2, 2))
                    # loop over intervals
                    for (i in idate[[day]]) {
                        if (!is.null(Covars_Out[[i]][[flux]])) {
                            itime <- sub('[^ ]+ ', '', i)
                            plot(x_cov, Covars_Out[[i]][[flux]][cov_sub], type = 'l', 
                                main = itime,
                                xlab = '', ylab = flux, 
                                col = plotting_covar_colors[flux],
                                panel.first = {grid(); abline(h = 0, col = 'lightgrey')})
                        }
                    }
                    dev.off()
                    # cospec & ogive (dyn)
                    fn_dynco <- paste0(path_folder, '/dailygraphs-ogives-dyn-', 
                        flux, '-', day, '.png')
                    png(fn_dynco, width = 14, height = 14, units = 'in', res = res, ...)
                    par(mfrow = c(8, 6), mar = c(2, 4, 2, 2))
                    # loop over intervals
                    for (i in idate[[day]]) {
                        # i <- idate[[day]][1]
                        if (!is.null(Covars_Out[[i]][[flux]])) {
                            itime <- sub('[^ ]+ ', '', i)
                            plot_cospec_ogive(Ogive_dyn_Out[[i]][[flux]], 
                                Cospec_dyn_Out[[i]][[flux]], freq,
                                ylab = flux, col = plotting_covar_colors[flux], 
                                model_cols = c('#F02E42', '#9A33DA'),
                                model_par = Cospec_dyn_Out[[i]]$model_coef[, flux]
                            )
                            om <- Ogive_dyn_Out[[i]][[flux]][1]
                            y <- par('usr')[3:4]
                            text(freq[1], y[which.max(abs(y - om))] + diff(y) * sign(om) * 0.05,
                                labels = itime, pos = 2)
                            # add tau
                            text(tail(freq, 1), y[which.min(abs(y - om))] - sign(om) * 
                                diff(y) * 0.15, pos = 4, labels = sprintf(
                                    'lag time =\n%1.1f secs', 
                                results[format(st, '%Y-%m-%d %H:%M') == i, .SD[[1]], 
                                    .SDcols = paste0('lag_dyn_', flux)]
                                ), cex = 1
                            )
                        }
                    }
                    dev.off()
                    # cospec & ogive (fix)
                    fn_fixco <- paste0(path_folder, '/dailygraphs-ogives-fix-', 
                        flux, '-', day, '.png')
                    png(fn_fixco, width = 14, height = 14, units = 'in', res = res, ...)
                    par(mfrow = c(8, 6), mar = c(2, 4, 2, 2))
                    # loop over intervals
                    for (i in idate[[day]]) {
                        if (!is.null(Covars_Out[[i]][[flux]])) {
                            itime <- sub('[^ ]+ ', '', i)
                            plot_cospec_ogive(Ogive_fix_Out[[i]][[flux]], 
                                Cospec_fix_Out[[i]][[flux]], freq,
                                ylab = flux, col = plotting_covar_colors[flux], 
                                model_cols = model_colors,
                                model_par = Cospec_fix_Out[[i]]$model_coef[, flux]
                            )
                            om <- Ogive_fix_Out[[i]][[flux]][1]
                            y <- par('usr')[3:4]
                            text(freq[1], y[which.max(abs(y - om))] + diff(y) * sign(om) * 0.05,
                                labels = itime, pos = 2)
                            # add tau
                            text(tail(freq, 1), y[which.min(abs(y - om))] - sign(om) * 
                                diff(y) * 0.15, pos = 4, labels = sprintf(
                                    'lag time =\n%1.1f secs', 
                                results[format(st, '%Y-%m-%d %H:%M') == i, .SD[[1]], 
                                    .SDcols = paste0('lag_fix_', flux)]
                                ), cex = 1
                            )
                        }
                    }
                    dev.off()
                }
            }
            cat(' done.\n')
        }
        # fix time zone
        results[, ':='(
            st = with_tz(st, tz_user),
            et = with_tz(et, tz_user)
        )]
        # transfer attributes
        if (processing_strategy != 'recursive' && (subint_return || ogives_return)) {
            anms <- names(attributes(results))
            tatts <- sapply(grep('covars|cospec|ogv|subintervals', anms, value = TRUE),
                \(a) attr(results, a), simplify = FALSE)
        }
        # check if minimal output is requested
        if (minimal_output) {
            results <- results[, .SD, .SDcols = c('st', 'et', 'n_values', 'Hz',
                'WD', 'Ustar', 'L', 'Zo', 'sUu', 'sVu', 'sWu', 'd', 'z_sonic', 
                'U_sonic', 'T_sonic', 'Ra', 'Rb_nh3',
                grep('^(phi|alpha|beta|w_bias)$', names(results), value = TRUE),
                grep('^(avg_|ht_|li_)', names(results), value = TRUE),
                grep('^flux_', names(results), value = TRUE)
                )
            ]
        }
        # convert to ibts
        if (processing_strategy != 'recursive' && as_ibts) {
            results <- as.ibts(results)
        }
        # assign attributes back
        if (processing_strategy != 'recursive' && (subint_return || ogives_return)) {
            for (a in names(tatts)) {
                setattr(results, a, tatts[[a]])
            }
        }
    }

    if (processing_strategy != 'recursive') {
        cat("\n************************************************************\n") 
        cat("operation finished @", format(Sys.time(), "%d.%m.%Y %H:%M:%S"), 
            "time elapsed: ", sprintf('%1.1f', 
                d <- difftime(Sys.time(), script_start)), 
            attr(d, 'units'), "\n")
        cat("************************************************************\n")  
    }

    return(results)

}

# convenience functions for theoretical cospec/ogive models
cospec_model <- function(fx, m, mu, A0, f = freq) {
    A0 / fx / (
        (1 + m * (f / fx) ^ (2 * mu)) ^ ((m + 1) / (2 * mu * m))
    )
}
ogive_model <- function(fx, m, mu, A0, f = freq) {
    rev(cumsum(rev(cospec_model(fx, m, mu, A0, f))))
}

# main wrapper function for parallel call
.wrapper_main <- function(tf_sd, tf_env) {
    # call main function
    out <- .ec_main(qs2::qd_read(tf_sd), qs2::qs_read(tf_env))
    # save to qdata
    qs2::qd_save(out, tf_sd)
    # return path
    tf_sd
}

# EC MAIN function
.ec_main <- function(dat, env_list) {
    # get variables
    for (what in names(env_list)) {
        assign(what, env_list[[what]])
    }
    # rm(env_list)
    # prepare ogive output
    if (create_dailygraphs || ogives_return) {
        e_ogive <- new.env()
        e_ogive$Cospec_dyn_Out <- e_ogive$Cospec_fix_Out <- e_ogive$Covars_Out <- 
            e_ogive$Ogive_fix_Out <- e_ogive$Ogive_dyn_Out <- list()
    }
    # call main function
    out <- dat[, {

        # get start of interval
        interval_start <- start_time[.BY[[1]]]
        # TODO (maybe for later): include NA where d_t>2*mean, remove entries where d_t < 0.5*mean

        # be verbose
        if ('parent_interval' %in% names(env_list)) {
            wint <- ' - subinterval '
        } else {
            wint <- ' - interval '
        }
        cat("\n\n~~~~~~~~\n", format(interval_start), wint, .GRP, 
            " of ", .NGRP, "\n", sep = '')

        # ugly copy because of inability to change .SD values
        SD <- copy(.SD)

        # check if any column contains finite data
        if (
            .N >= n_threshold && 
            SD[1, any(sapply(.SD, \(x) x > -3)), .SDcols = paste0(hl_vars, '_flag')]
        ) {

            # check NA values in scalars
            scalars <- input_scalars
            covariances <- input_covariances
            covariances_variables <- input_covariances_variables
            covariances_plotnames <- input_covariances_plotnames
            scalar_covariances <- input_scalar_covariances
            flux_variables <- input_flux_variables
            plot_timeseries <- input_plot_timeseries
            # fix_lag <- input_fix_lag
            # dyn_lag <- input_dyn_lag
            lag_functions <- input_lag_functions
            damping_reference <- input_damping_reference
            damp_region <- input_damp_region
            if (
                !is.null(flux_variables) && 
                SD[, anyNA(.SD), .SDcols = flux_variables]
            ) {
                # loop over flux_variables & check
                for (fv in unique(c(flux_variables, 'u', 'v', 'w', 'T'))) {
                    x <- SD[, v, env = list(v = fv)]
                    if (all(is.na(x))) {
                        cat('=> ', fv, ': measurement contains only NA values',
                            ' -> exclude from current interval\n', sep = '')
                        scalars <- scalars[!(scalars %in% fv)]
                        flux_variables <- flux_variables[!(flux_variables %in% fv)]
                        plot_timeseries <- plot_timeseries[!(names(plot_timeseries) %in% fv)]
                        sind <- grep(fv, covariances)
                        if (length(sind)) {
                            damping_reference <- damping_reference[
                                !(names(damping_reference) %in% covariances[sind])]
                            damp_region <- damp_region[
                                !(names(damp_region) %in% covariances[sind])]
                            covariances <- covariances[-sind]
                            # fix_lag <- fix_lag[-sind]
                            # dyn_lag <- dyn_lag[, -sind, drop = FALSE]
                            covariances_variables <- covariances_variables[-sind]
                            covariances_plotnames <- covariances_plotnames[-sind]
                            scalar_covariances <- scalar_covariances[-sind]
                        }
                    }
                }
            }

            # detrend sonic data (T, u, v + w)
            # ---------------------------------------------------------------------- 
            cat("~~~\ndeterending sonic data...\n")
            wind <- detrend_sonic_data(SD, detrending, rec_Hz)

            # calculate wind statistics: turbulence values and MOST parameters
            # ---------------------------------------------------------------------- 
            wind_stats <- wind_statistics(wind, z_canopy[[1]], z_ec[[1]], 
                ustar_method = ustar_method)

            # get dyn & fix lag
            dyn_lag <- sapply(lag_functions, \(x) x(SD[, WD[1]], dyn = TRUE))
            fix_lag <- sapply(lag_functions, \(x) x(SD[, WD[1]]))

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

                # despiking figures
                if (any(despike) && create_graphs) {
                    if (!dir.exists(path_folder)) {
                        dir.create(path_folder, recursive = FALSE)
                    }
                    # get end of current interval time in UTC
                    soi_user <- with_tz(interval_start, tz_user)
                    # eoi_user <- with_tz(end_time[.BY[[1]]], tz_user)
                    # get date in correct format
                    date_formatted <- format(soi_user, '%Y%m%d')
                    time2 <- format(soi_user, format = "%H%M")
                    plotname <- paste("despiked-timeseries", date_formatted, time2, 
                        sep="-") 
                    despike_vars <- names(despike)[despike]
                    jpeg(filename = paste0(path_folder, '/', plotname, ".jpg"), 
                        width = 900, height = sum(despike) * 150, 
                        quality = 80)
                        par(mfrow = c(length(despike_vars), 1), 
                            mar = c(2, 4, 2, 2), oma = c(2, 0, 0, 0))
                        for (d in despike_vars) {
                            SD[, {
                                n_spikes <- sum(is.na(dspk)) - sum(is.na(orig))
                                plot(Time, orig, type = 'l', col = 'indianred',
                                    ylab = d, xlab = '')
                                lines(Time, dspk)
                                mtext(text = sprintf('%i spikes removed', n_spikes),
                                    adj = 0
                                )
                            }, env = list(
                                orig = paste0(d, '_original'),
                                dspk = paste0(d, '_despiked')
                            )]
                        }
                    dev.off()
                }

                # detrend scalars
                # ------------------------------------------------------------------ 
                cat("~~~\ndetrending scalars...\n")
                detrended_scalars <- mapply(trend, y = scalar_list, method = 
                    detrending[scalars], MoreArgs = list(Hz_ts = rec_Hz), 
                    SIMPLIFY = FALSE
                )

                # calculate scalar averages and sd:
                # ------------------------------------------------------------------ 
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
            }

            # fix detrending resulting in NA
            for (fv in unique(c(scalars, 'u', 'v', 'w', 'T'))) {
                if (SD[, sum(is.finite(var)) < n_threshold, env = list(var = fv)]) {
                    cat('=> ', fv, ': measurement contains too many NA values',
                        ' -> exclude from flux processing\n', sep = '')
                    flux_variables <- flux_variables[!(flux_variables %in% fv)]
                    sind <- grep(fv, covariances)
                    if (length(sind)) {
                        damping_reference <- damping_reference[
                            !(names(damping_reference) %in% covariances[sind])]
                        damp_region <- damp_region[
                            !(names(damp_region) %in% covariances[sind])]
                        covariances <- covariances[-sind]
                        covariances_variables <- covariances_variables[-sind]
                        covariances_plotnames <- covariances_plotnames[-sind]
                        scalar_covariances <- scalar_covariances[-sind]
                    }
                }
            }

            # check if any fluxes can be derived
            if (has_flux <- length(covariances_variables) > 0) {
                # start of flux relevant data manipulation
                # ------------------------------------------------------------------ 
                # ------------------------------------------------------------------ 
                cat("~~~\nstarting flux evaluation...\n")

                # calculate covariances with fix lag time:
                # ------------------------------------------------------------------ 
                cat("\t- covariances\n")
                Covars <- get_covariance(SD, covariances_variables, 
                    covariances, n_period)

                # find maximum in dynamic lag time range:
                # ------------------------------------------------------------------ 
                cat("\t- dyn lag\n")
                dlag_max <- sapply(covariances, function(i, x, lag) {
                    find_dynlag(x[[i]], lag[, i])
                }, x = Covars, lag = dyn_lag)

                # find dynlag using pre-whitening
                dyn_lag_max <- dlag_max
                if (any(lag_dyn_calc_pw[covariances])) {
                    pw_covs <- intersect(covariances, 
                        names(lag_dyn_calc_pw)[lag_dyn_calc_pw])
                    pw_vars <- strsplit(pw_covs, split = 'x', fixed = TRUE)
                    names(pw_vars) <- pw_covs
                    cat('\t  ** pre-whitened, bootstrapped dyn lag\n')
                    if (create_graphs) {
                        if (!dir.exists(path_folder)) {
                            dir.create(path_folder, recursive = FALSE)
                        }
                        # get end of current interval time in UTC
                        soi_user <- with_tz(interval_start, tz_user)
                        # eoi_user <- with_tz(end_time[.BY[[1]]], tz_user)
                        # get date in correct format
                        date_formatted <- format(soi_user, '%Y%m%d')
                        time2 <- format(soi_user, format = "%H%M")
                        plotname <- paste("pwb-timelag", date_formatted, time2, 
                            sep="-") 
                    }
                    tlag_pw <- sapply(pw_covs,
                        \(pv) {
                            v <- pw_vars[[pv]]
                            cat('\t    - ', pv, '\n')
                            # lws, uws & LAG.MAX in secs
                            lws <- dyn_lag['lower', pv] / rec_Hz
                            uws <- dyn_lag['upper', pv] / rec_Hz
                            tlag_detection(
                                SD[, .SD, .SDcols = v],
                                mfreq = rec_Hz, Rboot = 100, 
                                lws = lws, uws = uws, 
                                LAG.MAX = lag_dyn_lagmax,
                                model = lag_dyn_model,
                                wdt = lag_dyn_wdt,
                                plot.it = create_graphs,
                                plot.dir = path_folder,
                                plot.name = paste0(plotname, '-', pv, '.jpg'),
                                smooth_fun = lag_dyn_smooth,
                                dyn_lag_ext = dlag_max['tau', pv]
                            )
                        }, simplify = FALSE
                    )
                    names(tlag_pw) <- pw_covs
                    # prepare output
                    pw_out <- unlist(tlag_pw)
                    names(pw_out) <- paste0('dyn_lag_pw_', names(pw_out))
                    # which dyn lag approach should be taken?
                    switch(lag_dyn_method
                        , 'simple-pw' = {
                            out <- sapply(tlag_pw, '[[', 'tl_pw')
                            out <- rbind(
                                dlag_max[1, pw_covs] - dlag_max[2, pw_covs] + out,
                                out
                            )
                            dyn_lag_max[, pw_covs] <- out
                        }
                        , 'boot-pw' = {
                            out <- sapply(tlag_pw, '[[', 'tl_pwb')
                            out <- rbind(
                                dlag_max[1, pw_covs] - dlag_max[2, pw_covs] + out,
                                out
                            )
                            dyn_lag_max[, pw_covs] <- out
                        }
                    )
                    # fix NA values (fallback to previous max cov method)
                    dyn_lag_max[, is.na(dyn_lag_max[2, ])] <- dlag_max[, is.na(dyn_lag_max[2, ])]
                } else if (lag_dyn_method != 'raw-cov') {
                    warning('fallback to "raw-cov" dyn lag')
                }
                # covariance function's standard deviation and mean values left and right of fix lag
                # ----------------------------------------------------------------
                # RE_RMSE (Eq. 9 in Langford et al. 2015)
                # -> ranges lo/hi (-/+180 to -/+150 secs? => define range as argument)
                m <- ifelse(n_period %% 2, (n_period + 1) / 2, 
                    n_period / 2 + 1)
                lo_range <- m - rev(gamma_time_window) * 60 * rec_Hz
                lo_ind <- seq(lo_range[1], lo_range[2])
                lo_ind <- lo_ind[lo_ind > 1]
                hi_range <- m + gamma_time_window * 60 * rec_Hz
                hi_ind <- seq(hi_range[1], hi_range[2])
                hi_ind <- hi_ind[hi_ind <= n_period]
                # -> sd_cov_low
                sd_cov_lo <- sapply(Covars, \(x) sd(x[lo_ind], na.rm = TRUE))
                # -> avg_cov_low
                avg_cov_lo <- sapply(Covars, \(x) mean(x[lo_ind], na.rm = TRUE))
                # -> sd_cov_hi
                sd_cov_hi <- sapply(Covars, \(x) sd(x[hi_ind], na.rm = TRUE))
                # -> avg_cov_hi
                avg_cov_hi <- sapply(Covars, \(x) mean(x[hi_ind], na.rm = TRUE))
                re_rmse <- sqrt(0.5 * (sd_cov_lo ^ 2 + avg_cov_lo ^ 2 +
                        sd_cov_hi ^ 2 + avg_cov_hi ^ 2))

                # cospectra for fixed & dynamic lags
                # ------------------------------------------------------------------
                cat("\t- co-spectra\n")
                # fix lag
                Cospec_fix <- get_cospectra(SD, Covars, covariances_variables, 
                    covariances, fix_lag[covariances], n_period)
                # dyn lag
                Cospec_dyn <- get_cospectra(SD, Covars, covariances_variables, 
                    covariances, dyn_lag_max[2, covariances], n_period)

                # ogives for fixed & dynamic lags 
                # ------------------------------------------------------------------
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
                    hi_ind <- seq(max(1, i_hi - cont_pts), min(length(x), i_hi + cont_pts))
                    (x[1] - mean(x[hi_ind]))/ x[1]
                })
                hi_cont_dyn <- sapply(Ogive_dyn, \(x) {
                    hi_ind <- seq(max(1, i_hi - cont_pts), min(length(x), i_hi + cont_pts))
                    (x[1] - mean(x[hi_ind]))/ x[1]
                })
                i_lo <- which(1 / freq <= low_cont_sec)[1]
                if (length(i_lo) != 1) stop('check argument "low_cont_sec"!')
                lo_cont_fix <- sapply(Ogive_fix, \(x) {
                    lo_ind <- seq(max(1, i_lo - cont_pts), min(length(x), i_lo + cont_pts))
                    (x[1] - mean(x[lo_ind])) / x[1]
                })
                lo_cont_dyn <- sapply(Ogive_dyn, \(x) {
                    lo_ind <- seq(max(1, i_lo - cont_pts), min(length(x), i_lo + cont_pts))
                    (x[1] - mean(x[lo_ind])) / x[1]
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
                # get bins
                q_bounds <- exp(seq(log(freq[1]), log(tail(freq, 1)), length.out = 4))
                q_bins <- findInterval(freq, q_bounds)

                # quality of ogive 
                ogive_quality_fix <- sapply(Ogive_fix, \(x) {
                    # ini <- c(fx = 0.05, m = 3 / 4, mu = 1 / 6, A0 = x[i_lo] / 50)
                    ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                    opt <- optim(ini, fit_og, ogive = x, control = list(maxit = 5e3))
                    if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                        return(rep(NA_real_, 6))
                    }
                    fog <- ogff(opt$par)
                    d <- fog - x
                    md <- tapply(d, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    mm <- tapply(fog, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    qv <- mean(md) / mean(mm)
                    c(1 - qv, fog[1], opt$par, 3 / 4)
                })
                ogive_quality_dyn <- sapply(Ogive_dyn, \(x) {
                    ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                    opt <- optim(ini, fit_og, ogive = x, control = list(maxit = 5e3))
                    if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                        return(rep(NA_real_, 6))
                    }
                    fog <- ogff(opt$par)
                    d <- fog - x
                    md <- tapply(d, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    mm <- tapply(fog, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    qv <- mean(md) / mean(mm)
                    c(1 - qv, fog[1], opt$par, 3 / 4)
                })
                # quality of base part of ogive 
                base_quality_fix <- sapply(Ogive_fix, \(x) {
                    ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                    opt <- optim(ini, fit_og, ogive = x, ind = i_r, 
                        control = list(maxit = 5e3))
                    if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                        return(NA_real_)
                    }
                    fog <- ogff(opt$par)
                    d <- fog - x
                    md <- tapply(d, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    mm <- tapply(fog, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    qv <- mean(md) / mean(mm)
                    1 - qv
                })
                base_quality_dyn <- sapply(Ogive_dyn, \(x) {
                    ini <- c(fx = 0.05, mu = 1 / 6, A0 = x[i_lo] / 50)
                    opt <- optim(ini, fit_og, ogive = x, ind = i_r, 
                        control = list(maxit = 5e3))
                    if (opt$convergence != 0 || opt$par['mu'] > mu_lim) {
                        return(NA_real_)
                    }
                    fog <- ogff(opt$par)
                    d <- fog - x
                    md <- tapply(d, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    mm <- tapply(fog, q_bins, \(x) mean(x ^ 6) ^ (1/6))
                    qv <- mean(md) / mean(mm)
                    1 - qv
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
                flux_modelled_fix <- ogive_quality_fix[2, ]
                ogive_quality_fix <- ogive_quality_fix[1, ]
                flux_modelled_dyn <- ogive_quality_dyn[2, ]
                ogive_quality_dyn <- ogive_quality_dyn[1, ]

                # should ogives be provided with output
                if (ogives_return || create_dailygraphs) {
                    int_start <- format(interval_start, '%Y-%m-%d %H:%M')
                    e_ogive$Cospec_fix_Out[[int_start]] <- c(
                        list(freq = freq, model_coef = ogive_par_fix), Cospec_fix)
                    e_ogive$Cospec_dyn_Out[[int_start]] <- c(
                        list(freq = freq, model_coef = ogive_par_dyn), Cospec_dyn)
                    e_ogive$Ogive_fix_Out[[int_start]] <- c(
                        list(freq = freq, model_coef = ogive_par_fix), Ogive_fix)
                    e_ogive$Ogive_dyn_Out[[int_start]] <- c(
                        list(freq = freq, model_coef = ogive_par_dyn), Ogive_dyn)
                    e_ogive$Covars_Out[[int_start]] <- c(
                        list(Hz = rec_Hz, N = n_period), Covars)
                }

                # empirical damping estimation, dyn and fix should have best reference (dyn/fix)...
                # ------------------------------------------------------------------
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

                if (length(scalar_covariances_only)) {
                    if (!is.null(Damping_fix)) {
                        fix_damping_pbreg <- sapply(Damping_fix, "[[", 1)
                        fix_damping_deming <- sapply(Damping_fix, "[[", 2)
                        names(fix_damping_pbreg) <- names(Damping_fix)
                        names(fix_damping_deming) <- names(Damping_fix)
                    } else {
                        fix_damping_pbreg <- fix_damping_deming <- setNames(
                            rep(NA_real_, length(scalar_covariances_only)),
                            scalar_covariances_only
                            )
                    }
                    if (!is.null(Damping_dyn)) {
                        dyn_damping_pbreg <- sapply(Damping_dyn, "[[", 1)
                        dyn_damping_deming <- sapply(Damping_dyn, "[[", 2)
                        names(dyn_damping_pbreg) <- names(Damping_dyn)
                        names(dyn_damping_deming) <- names(Damping_dyn)
                    } else {
                        dyn_damping_pbreg <- dyn_damping_deming <- setNames(
                            rep(NA_real_, length(scalar_covariances_only)),
                            scalar_covariances_only
                            )
                    }
                }
            }

            # get statistics on NA values and hard limits
            v_stats <- unlist(sapply(variables, \(x) {
                v <- get(x)
                if (na_limits[x]) {
                    # get -1, 1, 2 flags
                    vf <- get(paste0(x, '_flag'))
                    c(
                        below_thresh = vf[1] == -1,
                        na_before = sum(vf == 1),
                        hl_values = sum(vf == 2),
                        na_remaining = sum(is.na(v))
                    )
                } else {
                    # get # NA values only
                    c(
                        below_thresh = NA_integer_,
                        na_before = NA_integer_,
                        hl_values = NA_integer_,
                        na_remaining = sum(is.na(v))
                    )
                }
            }, simplify = FALSE))

            # write results:
            # -----------------------------------------------------------------------
            out <- c(
                list(
                    st = interval_start
                    , et = end_time[.BY[[1]]]
                    , Hz = rec_Hz
                    , n_values = .N
                    , n_subint =  if (subintervals) subint_n else NA_integer_
                )
                , as.list(v_stats)
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
                    # wind statistics
                    wind_stats
                    # bLSmodelR statistics
                    , setNames(
                        sqrt(wind_stats[c('var_u', 'var_v', 'var_w')]) / 
                            wind_stats['Ustar'],
                        c('sUu', 'sVu', 'sWu')
                    )
                    # wind direction
                    , WD = SD[, WD[1]]
                    # rotation parameters
                    , if (rotation_method == 'planar fit') {
                        c(
                            alpha = SD[, alpha[1]],
                            beta = SD[, beta[1]],
                            w_bias = SD[, w_bias[1]]
                        )
                    } else {
                        c (
                            phi = SD[, phi[1]]
                        )
                    }
                    # Ra(z_ec - d) & Rb
                    , Ra = Ra(wind_stats[['z_sonic']], wind_stats[['Ustar']], wind_stats[['L']], 
                        wind_stats[['Zo']], d = wind_stats[['d']], kv = 0.4)
                    , Rb_nh3 = Rb_nh3(wind_stats[['Ustar']], wind_stats[['T_sonic']] - 273.14, 
                        wind_stats[['Zo']], p_hPA = 960)
                    # scalar avg & sd
                    , if (length(input_scalars) > 0) {
                        as.list(c(
                            setNames(scalar_means[input_scalars], 
                                paste0('avg_', input_scalars))
                            , setNames(scalar_sd[input_scalars], 
                                paste0('sd_', input_scalars))
                        ))
                    }
                    , if (has_flux) {
                        as.list(c(
                            # fluxes
                            setNames(
                                flux_fix_lag[input_covariances],
                                paste0('flux_fix_', input_covariances)
                            )
                            , setNames(
                                flux_dyn_lag[input_covariances],
                                paste0('flux_dyn_', input_covariances)
                            )
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
                                flux_modelled_fix[input_covariances],
                                paste0('flux_modelled_fix_', input_covariances)
                            )
                            , setNames(
                                flux_modelled_dyn[input_covariances],
                                paste0('flux_modelled_dyn_', input_covariances)
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
                            , if (length(scalar_covariances_only) > 0) {
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
                            , if (length(scalar_covariances_only) > 0) {
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
                            , if (any(lag_dyn_calc_pw)) {
                                # print all lag times again
                                c(
                                    # raw-cov
                                    setNames(
                                        dlag_max['tau', pw_covs],
                                        paste0('dyn_lag_raw-cov_', 
                                            pw_covs)
                                    ),
                                    # pre-whitening
                                    pw_out
                                )
                            }
                        ))
                    }
                    ))
            )

            # sub-intervals:
            # -----------------------------------------------------------------------
            if (subintervals) {

                # be verbose
                cat('\n\n~~~ START SUBINTERVALS\n')

                cat(paste0('\n~~~\nprocessing sub-intervals (', subint_n, 
                        ' intervals @ ', round(avg_secs / subint_n / 60, 1), 
                        ' mins)\n'))

                # copy original .SD for sub-interval processing
                SDsub <- copy(SD)

                # copy environment
                env_sub <- copy(env_list)

                # fix entries to match subint
                env_sub$avg_secs <- env_list$avg_secs / subint_n
                env_sub$start_time <- seq(start_time[.BY[[1]]], end_time[.BY[[1]]], by = env_sub$avg_secs)[1:subint_n]
                env_sub$end_time <- env_sub$start_time + env_sub$avg_secs
                env_sub$n_period <- env_sub$avg_secs * rec_Hz
                env_sub$n_threshold <- env_sub$n_period * env_sub$thresh_period
                env_sub$freq <- rec_Hz * seq(floor(env_sub$n_period / 2)) / floor(env_sub$n_period / 2)
                env_sub$subintervals <- FALSE
                env_sub$subint_return <- FALSE
                env_sub$detrending <- env_list$subint_detrending
                # despiking has been done
                env_sub$despike[] <- FALSE

                # fix gamma_time_window
                if (max(env_sub$gamma_time_window) * 60 >= env_sub$avg_secs) {
                    env_sub$gamma_time_window <- c(5, 10) / 30 * env_sub$avg_secs / 60
                }

                # add extra figures subint sub-directory
                if (create_graphs) {
                    env_sub$path_folder <- paste0(env_list$path_folder, '/subintervals')
                    if (!dir.exists(env_sub$path_folder)) {
                        dir.create(env_sub$path_folder, recursive = TRUE)
                    }
                }

                # add bins
                SDsub[, bin := rleid(subint)]

                # add parent interval
                env_sub$parent_interval <- .BY[[1]]
                # env_sub$parent_interval <- list(
                #     st = start_time[.BY[[1]]],
                #     et = end_time[.BY[[1]]],
                #     bin = .BY[[1]]
                # )

                # use already rotated raw data
                SDsub[, c("WD", "phi", "urot", "vrot", "wrot") := .SD,
                    .SDcols = paste0('subint_', c("WD", "phi", "urot",
                            "vrot", "wrot"))]

                # use original T_sonic & scalars
                SDsub[, c('T', scalars) := .SD, 
                    .SDcols = paste0('subint_', c('T', scalars))]

                # run subintervals
                res_sub <- .ec_main(SDsub, env_sub)

                # get averages
                if (has_flux) {
                    avg_subint <- res_sub[, lapply(.SD, mean), 
                        .SDcols = c(names(wind_stats), 
                            paste0('flux_fix_', names(flux_fix_lag)), 
                            paste0('flux_dyn_', names(flux_dyn_lag)))]
                } else {
                    avg_subint <- res_sub[, lapply(.SD, mean), 
                        .SDcols = names(wind_stats)]
                }
                # unselect z_sonic etc.
                avg_sel <- grep('^(z_sonic|z_canopy|d)$', names(avg_subint), 
                    invert = TRUE, value = TRUE)
                # recalc wind_statistics
                wstats_subint <- avg_subint[, {
                    # Ustar
                    if (ustar_method[1] %in% c('neg_sqrt', 'fallback')) {
                        suppressWarnings(
                            Ustar <- sqrt(-cov_uw)
                        )
                    }
                    if (ustar_method[1] == 'double_sqrt' || (is.na(Ustar) &&
                            ustar_method[1] == 'fallback')) {
                        Ustar <- sqrt(sqrt(cov_uw ^ 2 + cov_vw ^ 2))
                    }
                    # L
                    L <- -Ustar ^ 3 * T_sonic / (0.4 * 9.8062 * cov_wT)
                    # Zo
                    suppressWarnings(
                        z0 <- optimize(
                            function(x, ustar, L, z, d, U) {
                                abs(U - calcU(ustar, x, L, z - d))
                            }, c(0, z_sonic * 1.1), ustar = Ustar, L = L, 
                            U = U_sonic, z = z_sonic, d = d)$minimum
                    )
                    if (z0 >= z_sonic * 1.09) z0 <- NA_real_
                    .(Ustar, L, z0)
                }]
                # get var
                var_subint <- res_sub[, lapply(.SD, var), .SDcols = avg_sel]
                # subset avg
                avg_subint <- avg_subint[, .SD, .SDcols = avg_sel]
                # fix names
                setnames(avg_subint, paste0(subint_prefix, 'avg_', names(avg_subint)))
                setnames(var_subint, paste0(subint_prefix, 'var_', names(var_subint)))
                setnames(wstats_subint, paste0(subint_prefix, 'recalc_', 
                        names(wstats_subint)))
                # bind together and append to output
                out <- c(out, as.list(cbind(avg_subint, wstats_subint, var_subint)))
                # attach subinterval results
                if (subint_return) {
                    # add parent interval
                    res_sub[, parent_interval := .BY[[1]]]
                    # attach
                    setattr(out, 'subintervals', res_sub)
                }

                cat('\n~~~ END SUBINTERVALS\n\n')

            } # END subintervals


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
                ts_vars <- names(plot_timeseries)[plot_timeseries]
                if ('parent_interval' %in% names(env_list)) {
                    plotname <- paste("timeseries", date_formatted, 
                        time2, 'parent', parent_interval, sep="-") 
                } else {
                    plotname <- paste("timeseries", date_formatted, time2, sep="-") 
                }
                jpeg(filename = paste0(path_folder, '/', plotname, ".jpg"), width = 600, 
                    height = (sum(plot_timeseries)) * 100, quality = 80)
                    ts_plot <- plot.tseries(
                        cbind(st = Time, as.data.frame(SD)),
                        wind, detrended_scalars, ts_vars,
                        plotting_var_colors, plotting_var_units,
                        if (subintervals) c(env_sub$start_time[1],
                            env_sub$end_time) else NULL
                    )
                    # fix time zone
                    attr(ts_plot$x.limits, 'tzone') <- tz_user
                    print(ts_plot)
                dev.off()

                # plot and save flux evaluation...
                # ------------------------------------------------------------------------
                for(i in covariances){
                    # i <- "w'TDL CH4'"
                    # i <- covariances[3]
                    if ('parent_interval' %in% names(env_list)) {
                        plotname <- paste("plots", date_formatted, 
                            time2, covariances_plotnames[i], 
                            'parent', parent_interval,
                            sep = "-")
                    } else {
                        plotname <- paste("plots", date_formatted, time2, 
                            covariances_plotnames[i], sep = "-")
                    }
                    # fix ylab
                    ylab <- sub('(.+)x(.+)', "<\\1'\\2'>", i)
                    jpeg(filename = paste0(path_folder, '/', plotname, ".jpg"), 
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
            list(list(out))
        } else {
            cat("less than ", round(thresh_period * 100),"% of data points in raw data - skipping interval.\n", sep = '')
            # return NULL
            NULL
        }
    }, by = bin]

    # check rows
    if (nrow(out) == 0) {
        return(NULL)
    }

    # get subintervals
    if (subint_return) {
        subints <- out[, rbindlist(lapply(V1, attr, 'subintervals'),
            fill = TRUE, use.names = TRUE)]
    }

    # bind lists to one data.table
    out <- out[, rbindlist(V1, fill = TRUE, use.names = TRUE)]

    if (subint_return) {
        setattr(out, 'subintervals', subints)
    }

    # output incl. ogives
    if (create_dailygraphs || ogives_return) {
        nogv <- names(e_ogive)
        natt <- tolower(sub('_Out', '', nogv))
        natt <- sub('ogive', 'ogv', natt)
        names(natt) <- nogv
        for (nm in nogv) setattr(out, natt[nm], e_ogive[[nm]])
    }

    # return
    out
}
# EC MAIN end


# parallel helper
.pef_wrapper <- function(ind, tf_cobj, tf_resid, tf_sonic, tf_ht,
    tf_licor) {
    if (length(ind) == 0) {
        return(NULL)
    }
    resid_list <- qs2::qd_read(tf_resid)
    utc_dates <- unique(c(resid_list$st_dates[ind], 
            resid_list$et_dates[ind]))
    formatted_dates <- gsub('-', '_', utc_dates, fixed = TRUE)
    out <- try(process_ec_fluxes(
        dates_utc = utc_dates,
        dates_formatted = formatted_dates,
        start_time = resid_list$start_time[ind],
        end_time = resid_list$end_time[ind],
        as_ibts = FALSE,
        processing_strategy = 'recursive',
        tf_cobj = tf_cobj,
        tf_sonic = tf_sonic,
        tf_ht = tf_ht,
        tf_licor = tf_licor
    ))
    # save to tmpfile on success
    if (inherits(out, 'try-error')) {
        # add intervals
        attr(out, 'start_time') <- resid_list$start_time[ind]
        attr(out, 'end_time') <- resid_list$end_time[ind]
        # return error
        return(out)
    }
    # get path
    tf_out <- sub('resid', digest::digest(ind), tf_resid)
    # save file
    qs2::qd_save(out, tf_out, warn_unsupported_types = FALSE)
    # return path
    tf_out
}

## quality criteria (instationarity) post-processing ----------------------------------------

# function to classify fluxes by Foken & Wichura 1996 criterion based on Foken et al., 2004 classification
# an altered, more robust classification scheme is provided as well
foken_wichura <- function(x, altered = FALSE, subint_pattern = '.*_avg_(flux_dyn|cov)_.*') {
    # get names
    xnms <- names(x)
    # check if subintervals exist
    vars <- grep(subint_pattern, xnms, value = TRUE)
    if (length(vars) > 0) {
        # convert to data.table
        if (is_ibts <- is.ibts(x)) {
            x <- as.data.table(x)
        }
        # get functions & table
        # Foken et al., 2004 (Handbook of Micromet.)
        breaks <- c(0, 15, 30, 50, 75, 100, 250, 500, 1000, Inf) / 100
        if (altered) {
            # robust version of Foken & Wichura 1996
            fun <- function(sub, all) abs((sub - all) / (sub + all))
            # adapted table
            breaks <- breaks / (2 + breaks)
            breaks[length(breaks)] <- Inf
        } else {
            # Foken & Wichura 1996
            fun <- function(sub, all) abs((sub - all) / all)
        }
        # get variable "reference"
        vref <- sub('.*_avg_', '', vars)
        # check
        if (!all(vr <- vref %in% xnms)) {
            stop('Cannot find refrence columns of: ', 
                paste0(vars[!vr], collapse = ', '))
        }
        # loop over variables
        for (i in seq_along(vars)) {
            vr <- vref[i]
            vs <- vars[i]
            vp <- sub(vr, '', vs)
            vd <- sub('_avg_', '_dev_', vp)
            vc <- sub('_avg_', '_fwclass_', vp)
            # add deviation
            x[, paste0(vd, vr) := fun(s, r), env = list(s = vs, r = vr)]
            # add criterion
            x[, paste0(vc, vr) := findInterval(dev, breaks), env = 
                list(dev = paste0(vd, vr))]
        }
        # convert back to ibts
        if (is_ibts) {
            x <- as.ibts(x)
        }
    } else {
        stop('no column names match argument "subint_pattern"')
    }
    # return
    x
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
    cat('wpl flux correction using', ifelse(dynamic_lag, 'dynamic', 'fixed'), 'lag fluxes\n')
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
        # ug/m3 -> g/m3
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
        # correct nh3_ugm3 (McDermitt et al. 2011, eq. 2)
        cat('correcting nh3 concentration due to h2o\n')
        out[, paste0(c('avg_nh3_ugm3', 'sd_nh3_ugm3'), '_raw') := 
            .(avg_nh3_ugm3, sd_nh3_ugm3)]
        out[, c('avg_nh3_ugm3', 'sd_nh3_ugm3') := .(
            avg_nh3_ugm3 * kappa_pe / kappa_p,
            sd_nh3_ugm3 * kappa_pe / kappa_p
        )]
    }
    # get co2 density/flux in g/m3
    if (co2_mmolm3 <- 'co2_mmolm3' %in% fluxes) {
        M_co2 <- 44.01
        # mmol/m3 -> g/m3
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
    out <- fields::interp.surface.grid(list(x = gel::T_deg, y = gel::p_kPa, 
            z = gel::kappa), list(x = tout, y = pout))
    c(
        kappa = out$z[2, 2],
        dkappa_P_at_T = (out$z[2, 3] - out$z[2, 1]) / (out$y[3] - out$y[1]) / 1000,
        dkappa_T_at_P = (out$z[3, 2] - out$z[1, 2]) / (out$x[3] - out$x[1])
    )
}

#' Merge Sonic and Other Instrument Data Based on Time
#'
#' This function merges sonic and HT8700 data based on time. The output contains the same times as the `basis` input. Values from `draw` will be repeated or dropped to match `basis` times. Licor data is optional and must be provided by `licor`.
#'
#' @param sonic A data.table containing the sonic data to be used as the basis for merging.
#' @param ht A data.table containing the HT8700 data to be merged. Default is NULL.
#' @param licor A data.table containing the Licor data to be merged. Default is NULL.
#' @param miro A data.table containing the Licor data to be merged. Default is NULL.
#' @return A data.table containing the merged data with the same times as `basis`.
#'
# merge sonic & other instrument data based on time
#   -> output contains the same times as 'basis'
#   -> values from 'draw' will be repeated or dropped to match 'basis' times
# merge_data <- function(sonic, ht = NULL, licor = NULL, miro = NULL, 
#     basis = 'sonic') {
merge_data <- function(basis, ..., ec_subset = TRUE) {
    # initialize all
    miro <- licor <- ht <- snc <- NULL
    # capture dots
    for (i in seq_len(...length())) {
        if (!is.null(...elt(i))) {
            switch(names(...elt(i))[3]
                , 'u' = snc <- ...elt(i)
                , 'nh3_ppb' = ht <- ...elt(i)
                , 'CO2D' = licor <- ...elt(i)
                , 'ch4_wet' = miro <- ...elt(i)
            )
        }
    }
    # prepare output
    n_out <- nrow(basis)
    out <- data.table(
        Time = basis[, Time],
        Hz = rep(NA_real_, n_out),
        u = rep(NA_real_, n_out),
        v = rep(NA_real_, n_out),
        w = rep(NA_real_, n_out),
        T = rep(NA_real_, n_out),
        sonic = rep(NA_character_, n_out),
        nh3_ppb = rep(NA_real_, n_out),
        nh3_ugm3 = rep(NA_real_, n_out), 
        ht_temp_amb = rep(NA_real_, n_out), 
        ht_press_amb = rep(NA_real_, n_out), 
        ht_oss = rep(NA_real_, n_out), 
        ht_peak_pos = rep(NA_real_, n_out),
        ht_alarm_code = rep(NA_character_, n_out),
        h2o_mmolm3 = rep(NA_real_, n_out),
        co2_mmolm3 = rep(NA_real_, n_out),
        li_temp_amb = rep(NA_real_, n_out),
        li_press_amb = rep(NA_real_, n_out),
        li_co2ss = rep(NA_real_, n_out),
        h2o_molfrac = rep(NA_real_, n_out),
        ch4_molfrac = rep(NA_real_, n_out),
        n2o_molfrac = rep(NA_real_, n_out)
    )
    if (n_out == 0) {
        return(out)
    }
    # fill sonic
    sonic_vars <- names(out)[2:7]
    sonic_orig <- sonic_vars
    if (!is.null(snc) && nrow(snc) > 0) {
        # times
        t_basis <- basis[, as.numeric(Time)]
        t_draw <- snc[, as.numeric(Time)]
        t0 <- t_basis[1]
        # ~ 10 Hz
        d_t <- median(diff(t_basis)) * 1.1
        # get matching indices
        indices <- match_times(t_basis - t0, t_draw - t0, d_t)
        # fill values
        out[indices[[1]], (sonic_vars) := snc[indices[[2]], sonic_orig,
            with = FALSE]]
    } else {
        # check if all columns ok
        if (all(sonic_vars %in% names(basis))) {
            # re-add sonic data
            add <- basis[, sonic_vars, with = FALSE]
            out[, (sonic_vars) := copy(add)]
        } else {
            # warn about missing sonic data
            warning('no sonic data provided!')
        }
    }
    # fill ht
    ht_vars <- names(out)[8:14]
    ht_orig <- c('nh3_ppb', 'nh3_ugm3', 'temp_amb', 'press_amb', 'oss', 'peak_pos', 
        'alarm_code')
    if (!is.null(ht) && nrow(ht) > 0) {
        # check alarm codes
        if (!('alarm_code' %in% names(ht))) {
            ht <- copy(ht)
            ht[, alarm_code := get_alarms(.SD)]
        }
        # times
        t_basis <- basis[, as.numeric(Time)]
        t_draw <- ht[, as.numeric(Time)]
        t0 <- t_basis[1]
        # ~ 10 Hz
        d_t <- median(diff(t_basis)) * 1.1
        # get matching indices
        indices <- match_times(t_basis - t0, t_draw - t0, d_t)
        # fill values
        out[indices[[1]], (ht_vars) := ht[indices[[2]], ht_orig, with = FALSE]]
    } else {
        # check if original names
        if ('oss' %in% names(basis)) {
            # check alarm codes
            if (!('alarm_code' %in% names(basis))) {
                basis <- copy(basis)
                basis[, alarm_code := get_alarms(.SD)]
            }
            # re-add ht data
            add <- basis[, ht_orig, with = FALSE]
            setnames(add, ht_orig, ht_vars)
            out[, (ht_vars) := copy(add)]
        } else if ('ht_oss' %in% names(basis)) {
            # re-add ht data
            add <- basis[, ht_vars, with = FALSE]
            out[, (ht_vars) := copy(add)]
        }
    }
    # fill licor
    licor_vars <- names(out)[15:19]
    licor_orig <- c('H2OD', 'CO2D', 'Temp', 'Pres', 'CO2SS')
    if (!is.null(licor) && nrow(licor) > 0) {
        t_basis <- basis[, as.numeric(Time)]
        t_licor <- licor[, as.numeric(Time)]
        t0 <- t_basis[1]
        # ~ 10 Hz
        d_t <- median(diff(t_basis)) * 1.1
        # get matching indices
        indices <- match_times(t_basis - t0, t_licor - t0, d_t)
        # fill values
        out[indices[[1]], (licor_vars) := licor[indices[[2]], licor_orig, 
                with = FALSE]]
    } else {
        # check if original names
        if ('CO2D' %in% names(basis)) {
            # re-add licor data
            add <- basis[, licor_orig, with = FALSE]
            setnames(add, licor_orig, licor_vars)
            out[, (licor_vars) := copy(add)]
        } else if ('co2_mmolm3' %in% names(basis)) {
            # re-add licor data
            add <- basis[, licor_vars, with = FALSE]
            out[, (licor_vars) := copy(add)]
        }
    }
    # fill miro
    miro_vars <- names(out)[20:22]
    miro_orig <- c('h2o', 'ch4_dry', 'n2o_dry')
    if (!is.null(miro) && nrow(miro) > 0) {
        t_basis <- basis[, as.numeric(Time)]
        t_miro <- miro[, as.numeric(Time)]
        t0 <- t_basis[1]
        # ~ 10 Hz
        d_t <- median(diff(t_basis)) * 1.1
        # get matching indices
        indices <- match_times(t_basis - t0, t_miro - t0, d_t)
        # fill values
        out[indices[[1]], (miro_vars) := miro[indices[[2]], miro_orig, 
                with = FALSE]]
    } else {
        # check if original names
        if ('ch4_dry' %in% names(basis)) {
            # re-add miro data
            add <- basis[, miro_orig, with = FALSE]
            setnames(add, miro_orig, miro_vars)
            out[, (miro_vars) := copy(add)]
        } else if ('ch4_molfrac' %in% names(basis)) {
            # re-add miro data
            add <- basis[, miro_vars, with = FALSE]
            out[, (miro_vars) := copy(add)]
        }
    }
    # return
    out
}

# # add sonic data to ht8700 data
# add_sonic <- function(x, path) {
#     # get from/to from ht-data
#     tr <- x[, range(Time)]
#     ft <- as.integer(format(tr, format = '%Y%m%d'))
#     # get sonic file paths
#     files <- dir(path, pattern = 'sonic')
#     int_files <- as.integer(sub('.*_(\\d{8})_\\d{6}([.]gz)?$', '\\1', files))
#     ind_files <- int_files >= ft[1] & int_files <= ft[2]
#     # read sonic data
#     dat <- rbindlist(lapply(file.path(path, files[ind_files]), read_windmaster_ascii))
#     # merge data
#     merge_data(x, dat)
# }

## pre-whitening dyn lag ----------------------------------------

# original source: https://github.com/domvit81/RFlux
tlag_detection <- function (dat, mfreq = 10, wdt = 5, 
    model = "ar", LAG.MAX = 10, lws = -LAG.MAX, uws = LAG.MAX, 
    Rboot = 100, plot.it = FALSE, boot_parallel = 'snow', 
    boot_ncpus = getOption('boot.ncpus', 1L), boot_cl = NULL,
    plot.dir = NULL, plot.name = NULL, smooth_fun = 'mean',
    dyn_lag_ext = NULL
) 
{

    require(zoo)
    require(boot)
    require(HDInterval)
    require(bayestestR)

    # get names
    nms <- names(dat[, 1:2])


    # replace NA values
    set <- na.omit(
        cbind(
            zoo::na.approx(dat[, v, env = list(v = nms[1])], na.rm = FALSE), 
            zoo::na.approx(dat[, v, env = list(v = nms[2])] , na.rm = FALSE)
        )
    )

    # convert & fix LAG.MAX, lws & uws
    lws <- max(-LAG.MAX, lws) * mfreq
    uws <- min(LAG.MAX, uws) * mfreq
    LAG.MAX <- LAG.MAX * mfreq

    # define lag window 
    lag_win <- LAG.MAX + (lws + 1):(uws + 1)

    # Unit root test based upon Breitung's variance ratio -> stationary ts???
    stat <- egcm_bvr.test(set[, 1])$p.val < 0.01 && egcm_bvr.test(set[, 2])$p.val < 0.01
    if (is.na(stat)) {
        # return all NA
        return(
            list(
                "tl_mc" = NA_integer_,
                "tl_pw" = NA_integer_,
                "tl_pwb" = NA_integer_,
                "tl_pwb_lci"= NA_integer_,
                "tl_pwb_uci"= NA_integer_,
                "cor_pw" = NA_real_,
                "cor_pwb" = NA_real_,
                "cov_mc" = NA_real_, 
                "cov_pw" = NA_real_, 
                "cov_pwb" = NA_real_  
            )
        )
    } else if (stat) {
        x <- set[, 1]
        y <- set[, 2]
    } else {
        x <- diff(set[, 1])
        y <- diff(set[, 2])
    }

    if (model == "ar") {
        o.max <- floor(10 ^ 2 * log10(nrow(set)))
        ar.resx <- ar(x, aic = TRUE, order.max = o.max)
        ar.resy <- ar(y, aic = TRUE, order.max = o.max)
        x1 <- stats::filter(x, filter = c(1, -ar.resx$ar), method = "convolution", 
            sides = 1)
        y1 <- stats::filter(y, filter = c(1, -ar.resx$ar), method = "convolution", 
            sides = 1)
        x2 <- stats::filter(x, filter = c(1, -ar.resy$ar), method = "convolution", 
            sides = 1)
        y2 <- stats::filter(y, filter = c(1, -ar.resy$ar), method = "convolution",
            sides = 1)
    }


    if (model == "arima") {
        require(forecast)
        filter.mod <- function(x, model) {
            x <- x - mean(x, na.rm = TRUE)
            if (length(model$Delta) >= 1) 
                x <- stats::filter(x, filter = c(1, -model$Delta), 
                    method = "convolution", sides = 1)
            if (length(model$theta) >= 1 && any(model$theta != 0)) 
                x <- stats::filter(na.omit(x), filter = -model$theta,
                    method = "recursive", sides = 1)
            if (length(model$phi) >= 1 && any(model$phi != 0)) 
                x <- stats::filter(x, filter = c(1, -model$phi), 
                    method = "convolution", sides = 1)
            x
        }
        mod1 <- forecast::auto.arima(x, d = 0)
        mod2 <- forecast::auto.arima(y, d = 0)
        x1 <- filter.mod(x, model = mod1$model)
        y1 <- filter.mod(y, model = mod1$model)
        x2 <- filter.mod(x, model = mod2$model)
        y2 <- filter.mod(y, model = mod2$model)
    }


    ## Maximum Covariance Procedure    
    ccf_mc <- ccf(x = as.vector(egcm_detrend(set[, 1])), 
        y = as.vector(egcm_detrend(set[, 2])), lag.max = LAG.MAX, plot = FALSE, 
        type = "covariance", na.action = na.pass
    )$acf

    # get time lag + window
    tl_mcw <- which.max(
        abs(ccf_mc)[lag_win]
    ) + LAG.MAX + lws

    # get covariance at time lag
    cov_mcw <- ccf_mc[tl_mcw]

    ## PreWhitening - Standard procedure
    ccf_pw <- ccf(x1, y1, na.action = na.pass, plot = FALSE, lag.max = LAG.MAX)$acf

    # get time lag + window
    tl_pw <- which.max(
        abs(ccf_pw)[lag_win]
    ) + LAG.MAX + lws

    # get covariance at time lag
    cor_pw <- ccf_pw[tl_pw]
    cov_pw <- ccf_mc[tl_pw]


    ## PreWhitening + BOOTSTRAPPING + Smoothing
    boot_fun <- function(x, y) {
        boot::tsboot(cbind(x, y), function(x) ccf(x[, 1], x[, 2], na.action = na.pass, 
            plot = FALSE, lag.max = LAG.MAX)$acf, R = Rboot, sim = 'fixed', 
            l = LAG.MAX * 2, parallel = boot_parallel, ncpus = boot_ncpus, cl = boot_cl)
    }

    # bootstrap
    bootccf_cs <- boot_fun(x1, y1)
    ccf_cs <- colMeans(bootccf_cs$t, na.rm = TRUE)
    bootccf_sc <- boot_fun(x2, y2)
    ccf_sc <- colMeans(bootccf_sc$t, na.rm = TRUE)

    # fix wdt usage
    fnm <- names(formals(smooth_fun))
    if ('wdt' %in% fnm) {
        old_fun <- smooth_fun
        smooth_fun <- \(x) old_fun(x, wdt = wdt)
    }

    # smooth
    ccfs_cs <- zoo::na.locf(
        zoo::na.locf(
            data.table::frollapply(ccf_cs, n = wdt, FUN = smooth_fun, 
                align = 'center', fill = NA), na.rm = FALSE
        ), fromLast = TRUE
    )
    ccfs_sc <- zoo::na.locf(
        zoo::na.locf(
            data.table::frollapply(ccf_sc, n = wdt, FUN = smooth_fun, 
                align = 'center', fill = NA), na.rm = FALSE
        ), fromLast = TRUE
    )
    ccfs_csb <- apply(bootccf_cs$t, MARGIN = 1, function(x) {
        which.max(
            abs(zoo::na.locf(
            zoo::na.locf(
                data.table::frollapply(x, n = wdt, FUN = smooth_fun, 
                align = 'center', fill = NA),
            na.rm = FALSE), fromLast = TRUE)[lag_win])
        ) + LAG.MAX + lws
    })
    ccfs_scb <- apply(bootccf_sc$t, MARGIN = 1, function(x) {
        which.max(
            abs(zoo::na.locf(
            zoo::na.locf(
                data.table::frollapply(x, n = wdt, FUN = smooth_fun, 
                align = 'center', fill = NA),
            na.rm = FALSE), fromLast = TRUE)[lag_win])
        ) + LAG.MAX + lws
    })


    # get ci interval
    hdis_cs <- HDInterval::hdi(ccfs_csb, credMass = .95);
    hdis_sc <- HDInterval::hdi(ccfs_scb, credMass = .95);

    # ??
    maps <- round(
        c(
            bayestestR::map_estimate(abs(ccfs_csb + rnorm(length(ccfs_csb), 0, 0.0001)))$MAP_Estimate,
            bayestestR::map_estimate(abs(ccfs_scb + rnorm(length(ccfs_scb), 0, 0.0001)))$MAP_Estimate
        ), 0
    )

    corr_est_s <- c(ccfs_cs[maps[1]], ccfs_sc[maps[2]])
    corr_ind <- which.max(abs(corr_est_s))
    corr_max <- corr_est_s[corr_ind]
    peak_ref <- maps[corr_ind]
    switch(corr_ind
        , hdis <- as.vector(hdis_cs)
        , hdis <- as.vector(hdis_sc)
    )
    # get boot time lag cov
    cov_pwb <- ccf_mc[peak_ref]

    ## PLOT
    if (plot.it) {
        
        jpeg(filename = file.path(plot.dir, plot.name), width = 700, 
                    height = 700, quality = 100)

            n_xax <- 9
            ylab_paren <- sprintf('(%s,%s)', nms[1], nms[2])
            ylab_paren_inv <- sprintf('(%s,%s)', nms[2], nms[1])
            par(mfrow = c(3, 1), mar = c(5, 4, 2, 1), oma = c(1, 1, 5, 0.5), las = 0, 
                cex.axis = 1.3, cex.lab = 1.3)
            # cross-cov
            plot((-LAG.MAX:LAG.MAX), ccf_mc, ylab = paste("cross-cov", ylab_paren), 
                xlab = "Lag (sec)", type = "h", col = "grey68",  #xlim = c(lws, uws),
                ylim = c(min(ccf_mc * 1.05, 0), max(ccf_mc * 1.05, 0)), xaxt = "n")
            axis(side = 1, at = seq(-LAG.MAX, LAG.MAX, length.out = n_xax), 
                labels = seq(-LAG.MAX, LAG.MAX, length.out = n_xax) / mfreq)  
            # axis(side = 1, at = seq(lws, uws, length.out = n_xax), 
            #     labels = seq(lws, uws, length.out = n_xax) / mfreq)  
            polygon(x = c(hdis[1]: hdis[2], hdis[2]: hdis[1]) - LAG.MAX - 1, 
                y = c(ccf_mc[hdis[1]: hdis[2]], rep(0, hdis[2] - hdis[1] + 1)), 
                col = "lightblue", border = "lightblue")
            # boot time lag cov
            segments(x0 = peak_ref - LAG.MAX - 1, y0 = 0, y1 = cov_pwb, col = 2, 
                lwd = 2)
            # max cov
            segments(x0 = tl_mcw - LAG.MAX - 1, y0 = 0, y1 = cov_mcw, col = 1, 
                lwd = 1)
            # cov as line
            lines((-LAG.MAX:LAG.MAX), ccf_mc, col = 1, lwd = 2)
            mtext(side = 3, line = .5, adj = 0, 
                paste0("Peak at ", (tl_mcw - LAG.MAX - 1) / mfreq, " sec"), cex = 1.1) 
            # mtext(side = 3, line = .5, adj = 1, "a", cex = 1.5, font = 2) 
            box(lwd = 1.5)
            mtext(side = 3, line = 3, cex = 1.25, col = 2,
                paste0("Time lag at: ", (peak_ref - LAG.MAX - 1) / mfreq, " sec"))
            # add dyn lag window
            abline(v = c(lws, uws), lwd = 2 , col = 1, lty = 3)
            if (!is.null(dyn_lag_ext)) {
                tl_ext <- dyn_lag_ext + LAG.MAX + 1
                cov_ext <- ccf_mc[tl_ext]
                # external dyn lag
                segments(x0 = dyn_lag_ext, y0 = 0, y1 = cov_ext, 
                    col = 'dodgerblue3', lwd = 1)
                mtext(side = 3, line = .5, adj = 1, 
                    paste0("ext. lag at ", dyn_lag_ext / mfreq, " sec"), cex = 1.1) 
            }
            # cross-cor 1
            ylim <- c(
                min(ccf_sc, ccf_cs, -4 / sqrt(length(x))), 
                max(ccf_sc, ccf_cs, 4 / sqrt(length(x)))
            )
            plot((-LAG.MAX:LAG.MAX), ccf_cs, ylab = paste("pwb cross-cor", ylab_paren), 
                xlab = "Lag (sec)", type = "h", col = "grey68", #xlim = c(lws, uws),
                ylim = ylim, xaxt = "n")    
            lines((-LAG.MAX:LAG.MAX), ccfs_cs, col = 1, lwd = 2)
            axis(side = 1, at = seq(-LAG.MAX, LAG.MAX, length.out = n_xax), 
                labels = seq(-LAG.MAX, LAG.MAX, length.out = n_xax) / mfreq)  
            abline(h = c(-3.291, 3.291) / sqrt(length(x) * 13), col = 4, lty = 2, lwd = 2)
            points(x = maps[1] - LAG.MAX - 1, y = ylim[1], 
                pch = 24, col = 1, cex = 1.25, bg = "red")
            # mtext(side = 3, line = 1, adj = 1, "b", cex = 1.5, font = 2) 
            mtext(side = 3, line = .5, adj = 0, 
                paste0("Peak at ", (maps[1] - LAG.MAX - 1) / mfreq , " sec"), cex = 1.1) 
            box(lwd = 1.5)
            # add dyn lag window
            abline(v = c(lws, uws), lwd = 2 , col = 1, lty = 3)
            # cross-cor 2
            plot((-LAG.MAX:LAG.MAX), ccf_sc, ylab = paste("pwb cross-cor", ylab_paren_inv), 
                xlab = "Lag (sec)", type = "h", col = "grey68", #xlim = c(lws, uws),
                ylim = ylim, xaxt = "n")    
            lines((-LAG.MAX:LAG.MAX), ccfs_sc, col = 1, lwd = 2)
            axis(side = 1, at = seq(-LAG.MAX, LAG.MAX, length.out = n_xax), 
                labels = seq(-LAG.MAX, LAG.MAX, length.out = n_xax) / mfreq)  
            abline(h = c(-3.291, 3.291) / sqrt(length(x) * 13), col = 4, lty = 2, lwd = 2)
            points(x = maps[2] - LAG.MAX - 1, y = ylim[1], 
                pch = 24, col = 1, cex = 1.25, bg = "red")
            # mtext(side = 3, line = 1, adj = 1, "c", cex = 1.5, font = 2) 
            mtext(side = 3, line = .5, adj = 0, 
                paste0("Peak at ", (maps[2] - LAG.MAX - 1) / mfreq , " sec"), cex = 1.1) 
            box(lwd = 1.5)
            # add dyn lag window
            abline(v = c(lws, uws), lwd = 2 , col = 1, lty = 3)

        dev.off()

    }

    list(
        "tl_mc" = tl_mcw - LAG.MAX - 1, 
        "tl_pw" = tl_pw - LAG.MAX - 1, 
        "tl_pwb" = peak_ref - LAG.MAX - 1, 
        "tl_pwb_lci"= hdis[1] - LAG.MAX - 1, 
        "tl_pwb_uci"= hdis[2] - LAG.MAX - 1, 
        "cor_pw" = cor_pw, 
        "cor_pwb" = corr_max, 
        "cov_mc" = cov_mcw, 
        "cov_pw" = cov_pw, 
        "cov_pwb" = cov_pwb  
    )

}

## despiking of timeseries ----------------------------------------

# move outside
despike_filter1 <- function(x, flt = filt) {
    m <- median(x, na.rm = TRUE)
    w <- 1 / (1 + abs(x - m) ^ 2)
    sum(x * w * flt, na.rm = TRUE) / sum(w * flt, na.rm = TRUE)
}
despike_filter2 <- function(i, ma, c_orig, n = n_filt2, 
    quant = 0.95) {
    sd(ma[i], na.rm = TRUE) + 
        quantile(abs(ma[i] - c_orig[i]), quant, na.rm = TRUE) +
        # sd(c_orig[i], na.rm = TRUE)
        mad(c_orig[i], na.rm = TRUE)
}
despike_timeseries <- function(dat, scalar, filter_width = 10, 
    qval = 0.95, qwidth = 30, qmult = 4, 
    filter1 = despike_filt1, filter2 = despike_filter2) {
    cat('despiking', scalar, '\n')
    Hz <- dat[, Hz[1]]
    md_filters <- getOption('md.filter.function.list')
    filt <- md_filters$BmNuttall(filter_width * Hz)
    n_filt2 <- qwidth * Hz
    sna_before <- dat[, sum(is.na(scal)), env = list(scal = scalar)]
    # loop over bins
    for (b in dat[, unique(bin)]) {
        cat('\r\r', b, '/', dat[, max(bin)])
        desp_conc <- dat[, {
            # get index
            ext <- ind <- which(bin == b)
            # # debugging
            # orig <- scal[ind]
            # get length of filter
            n <- length(filt)
            n2 <- 2 * n
            ind2 <- seq_along(ind) + n2
            # extended index
            if (ext[1] <= n2) {
                ext <- c(rep(ext[1], n2), ext)
            } else {
                ext <- c(ext[1] - (n2:1), ext)
            }
            if (ext[length(ext)] >= .N - n2 + 1) {
                ext <- c(ext, rep(ext[length(ext)], n2))
            } else {
                ext <- c(ext, ext[length(ext)] + 1:n2)
            }
            # filter 1
            # cat('-> apply filter 1\n')
            ma0 <- frollapply(c_ext <- scal[ext], n = n, 
                FUN = despike_filter1, align = 'center', flt = filt
            )
            # filter 2
            # cat('-> apply filter 2\n')
            mq <- frollapply(seq_along(ma0), n = n_filt2, FUN = despike_filter2, 
                align = 'center', ma = ma0, c_orig = c_ext, quant = 0.75)
            c1 <- scal[ind]
            ma <- ma0[ind2]
            d <- ma - scal[ind]
            qd <- mq[ind2]
            # qd <- mq[ind2] + sd(ma, na.rm = TRUE)
            qthresh <- qd * qmult
            flag <- abs(d) > qthresh
            # c1 <- orig
            c1[flag] <- NA
            # plot(Time[ind], orig, type = 'l')
            # lines(Time[ind], ma, col = 'red', lwd = 2)
            # lines(Time[ind], ma + qthresh, col = 'red')
            # lines(Time[ind], ma - qthresh, col = 'red')
            # points(Time[ind][flag], orig[flag], pch = 20, col = 'blue')
            # browser()
            c1
        }, env = list(scal = scalar)]
        # assign
        dat[b == bin, (scalar) := desp_conc]
    }
    sna_after <- dat[, sum(is.na(scal)), env = list(scal = scalar)]
    cat(' ->', sna_after - sna_before, 'spikes removed\n')
    invisible(dat)
}


## egcm package related functions for pwb lag time ----------------------------------------

# copy-paste from depricated egcm package
# => egcm package is no longer required

# first function
egcm_bvr.test <- function(Y, detrend=FALSE) {
	# Tests for a unit root of an AR(1) process using the variance ratio
	# test of Breitung.

	if (!exists("egcm_bvr_qtab")) {
		stop("Could not find quantile table egcm_bvr_qtab")
	}
	
    DNAME <- deparse(substitute(Y))
	STAT <- egcm_bvr_rho (Y, detrend=detrend)
	
	if (detrend) {
		PVAL <- egcm_quantile_table_interpolate(egcm_bvr_qtab_detrended, length(Y), STAT)
	} else {
		PVAL <- egcm_quantile_table_interpolate(egcm_bvr_qtab, length(Y), STAT)
	}
    METHOD <- "Breitung Variance Ratio Test for a Unit Root"
    names(STAT) <- "rho"
    alternative <- "stationary"
    structure(list(statistic = STAT, alternative = alternative, 
        p.value = PVAL, method = METHOD, data.name = DNAME), 
        class = "htest")
    
}

# second function
egcm_detrend <- function (Y) {
	# Removes a linear trend from Y
	if (is(Y, "zoo")) {
		X <- as.numeric(index(Y))
		X <- X - X[1]
		Y <- coredata(Y)
	} else {
		X <- 1:length(Y)
	}
	beta = cov(X,Y)/var(X)
	alpha = mean(Y) - beta * mean(X)
	eps = Y - alpha - beta * X
	eps
}

##  • helper functions ====================

egcm_bvr_rho <- function (Y, detrend=FALSE) {
	# Calculates the variance ratio statistic rho described in equation (5) of
	#   Breitung, Jorg (2001).  Nonparametric tests for unit roots and cointegration,
	#   Journal of Econometrics, 108, 343-363.

	if (detrend) {
		y <- detrend(Y)
	} else {
		y <- coredata(Y) - mean(coredata(Y))
	}

	ys <- cumsum(y)
	n <- length(y)
	rho_num <- (1 / n^2) * sum(ys^2)
	rho_den <- (1 / n) * sum(y^2)
	# Note factor of 1/n has been added that was not in the original paper.
	rho <- (rho_num / rho_den) / n
	rho
}

egcm_bvr_quantiles <- function(sample_size=100, nrep=40000, 
	q=c(0.001, 0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999),
	sd=1, detrend=FALSE) {
	# Calculates quantiles of the egcm_bvr_rho function under the assumptions
	# a_0 = 0 and a_1 = 1.
	qvals <- replicate(nrep, egcm_bvr_rho(egcm_rar1(sample_size, sd=sd), detrend=detrend))
	quantile(qvals, q)
}

egcm_bvr_quantile_table <- function(nrep=40000,
	q=c(0.001, 0.01, 0.025, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.999),
	n=c(25, 50, 100, 250, 500, 750, 1000, 1250), sd=1, detrend=FALSE) {
	# Calculates a table of quantile values by sample size of the egcm_bvr_rho function
	# under the assumption a_0=0 and a_1=1.
	df <- do.call("cbind", lapply(n, function(nv) c(nv, egcm_bvr_quantiles(nv, nrep, q, sd, detrend=detrend))))
	df <- as.data.frame(df)
	colnames(df) <- n
	df <- cbind(data.frame(quantile=c(NA,q)), df) 
	df
}

# The following table was generated using the call
#   egcm_bvr_qtab <- egcm_bvr_quantile_table()
egcm_bvr_qtab <- structure(list(quantile = c(NA, 0.001, 0.01, 0.025, 0.05, 0.1, 
0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999), `25` = c(25, 0.00332637455591009, 
0.00585745961329492, 0.00798960027086461, 0.0105373080553876, 
0.0150004420791956, 0.0220490768870461, 0.0517838845245987, 0.0798377158634759, 
0.0879813700851962, 0.0921303890798318, 0.0946062893660788, 0.0965695782064616, 
0.0988687421743598), `50` = c(50, 0.00312981685753421, 0.00580422085442258, 
0.00784802426252351, 0.0103672783549273, 0.0144821609570802, 
0.0215786933909018, 0.0514264840801606, 0.0796085968951424, 0.0874329498630163, 
0.0917242451170734, 0.0941587848720424, 0.0961805629259171, 0.0985460669231598
), `100` = c(100, 0.0030302822493761, 0.00558912528222194, 0.00772184312229409, 
0.0104166077319071, 0.0146476080970805, 0.0216939692835852, 0.0513220929469567, 
0.0797959488410655, 0.087631621168659, 0.0919207577499628, 0.094339872897165, 
0.0962604330637953, 0.0986130828033007), `250` = c(250, 0.00289850174980398, 
0.00553139542288896, 0.00762950737159783, 0.0101319765011695, 
0.0142733424169059, 0.0214928641627112, 0.0514630620014147, 0.0793934536098819, 
0.0873851099727652, 0.0915855661457587, 0.0941517184573074, 0.0962902144155935, 
0.098677267337584), `500` = c(500, 0.00305710661105332, 0.00544600639541136, 
0.00754559550481492, 0.0101011442885638, 0.0142472134777352, 
0.0212781070264034, 0.0511920886616309, 0.0793889207266449, 0.0873707165385677, 
0.0917426949005179, 0.0941981056662748, 0.0960242644480482, 0.0985074830761626
), `750` = c(750, 0.00301007161402827, 0.00548558465584746, 0.00767433432579238, 
0.0101582135161997, 0.0143894210957233, 0.0215221766847057, 0.0509050619344435, 
0.0795339082392001, 0.0876363386055488, 0.0918635283740264, 0.0942589928805052, 
0.0962438714606306, 0.0986682091672022), `1000` = c(1000, 0.00298323901478343, 
0.00551493610235411, 0.0076654895410273, 0.010276288911672, 0.0144242639907012, 
0.0213628350811129, 0.051228078374354, 0.0794048282306855, 0.087453070396734, 
0.0916568658959338, 0.0941431922060523, 0.0961921269254635, 0.0985127878537571
), `1250` = c(1250, 0.00306253528257504, 0.00537748023783321, 
0.0074666813665504, 0.0100679219782415, 0.014293190071986, 0.0215756981612748, 
0.0515638225099148, 0.0796758586257054, 0.0877612789797003, 0.0918409286997245, 
0.0942173357670156, 0.0962726790088758, 0.0985953726395397)), .Names = c("quantile", 
"25", "50", "100", "250", "500", "750", "1000", "1250"), row.names = c("", 
"0.1%", "1%", "2.5%", "5%", "10%", "20%", "50%", "80%", "90%", 
"95%", "97.5%", "99%", "99.9%"), class = "data.frame")

# The following table was generated using the call
#   egcm_bvr_qtab <- egcm_bvr_quantile_table(detrend=TRUE)
egcm_bvr_qtab_detrended <- structure(list(quantile = c(NA, 0.001, 0.01, 0.025, 0.05, 0.1, 
0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.999), `25` = c(25, 0.00168712405179854, 
0.00253217722042582, 0.00312181424925527, 0.00376926655109628, 
0.00474287543098176, 0.00633222666543711, 0.0106271227725071, 
0.0162771978932732, 0.0189023637742038, 0.0205143658341207, 0.0215287042894976, 
0.0225526147575758, 0.023875900170275), `50` = c(50, 0.00140880842231724, 
0.00229920683777919, 0.00283678960906636, 0.00349924284771183, 
0.00443242586703879, 0.00601680334452872, 0.0102746842717097, 
0.0159982270657943, 0.0185884185297237, 0.020253402175764, 0.0213145375046391, 
0.022280873629319, 0.0236687999417853), `100` = c(100, 0.00133476253824772, 
0.00220190519602716, 0.00281996248275345, 0.0034627339157511, 
0.00441240940690256, 0.00597919462373363, 0.0102019645534696, 
0.0159436827883071, 0.0185459998901182, 0.0202027655899719, 0.0212572406658586, 
0.0222141349110763, 0.0234491859902589), `250` = c(250, 0.00139751874609496, 
0.00220312428847516, 0.00277934566074049, 0.00340697821574804, 
0.00442693527992565, 0.00600849012813196, 0.0102029787959123, 
0.0159475349576176, 0.0184931201760397, 0.0201428921249319, 0.0212207606044985, 
0.0221897505899013, 0.0235601948472796), `500` = c(500, 0.00138222144173778, 
0.00219573407632577, 0.00278365119795981, 0.00340897111860022, 
0.00438328737782025, 0.00591410492263717, 0.010134997154612, 
0.0159497730503104, 0.0185564871388565, 0.0201649557104934, 0.021305788002882, 
0.0222391953355076, 0.0234982443840029), `750` = c(750, 0.00136926617120187, 
0.00220627938658952, 0.00278207347911421, 0.00343723108464136, 
0.00441177094147542, 0.00593155938053703, 0.0101212384027732, 
0.0158329539639887, 0.0184505021492885, 0.0200971092094903, 0.0212134064888182, 
0.0221790396015763, 0.0234741784797023), `1000` = c(1000, 0.00137427108101005, 
0.00220531383360885, 0.00280374191430434, 0.00345851096304081, 
0.00439460350000821, 0.00595719408311324, 0.0101981868726643, 
0.0159558007598192, 0.018541579368193, 0.0201872908757718, 0.0213201944053304, 
0.022205002688817, 0.0234768963658856), `1250` = c(1250, 0.00135100041059188, 
0.00219731203488951, 0.00278827000837087, 0.00342866371049778, 
0.00437121227906767, 0.00594994735644952, 0.0101282778302252, 
0.0157930542963664, 0.0183884536165205, 0.0200778640062674, 0.0211651097698149, 
0.0221193040226837, 0.0233513907649339)), .Names = c("quantile", 
"25", "50", "100", "250", "500", "750", "1000", "1250"), row.names = c("", 
"0.1%", "1%", "2.5%", "5%", "10%", "20%", "50%", "80%", "90%", 
"95%", "97.5%", "99%", "99.9%"), class = "data.frame")

egcm_rar1 <- function (n, a0=0, a1=1, trend=0, sd=1, x0=0) {
	# Generates a vector of length n representing a simulation of an AR(1)
	# process   X[k] = a0 +  a1 * X[k-1] + eps
	# where eps is an i.i.d. normal random variate with mean 0 and standard
	# deviation sd.  
    #
    # If trend is non-zero, then returns a realization of a trend-stationary 
    # AR(1) process.  E.g., the process is defined by the relations:
    #    R[k] = a0 + a1 * R[k-1] + eps
    #    X[k] = k * trend + R[k]
	eps <- rnorm(n, 0, sd)
	X <- numeric(n)
	xp <- x0
	for (k in 1:n) {
		X[k] <- xp <- a0 + a1 * xp + eps[k]
	}
	X + trend * (1:n)
}

egcm_quantile_table_interpolate <- function (qtab, sample_size, stat, stop.on.na=FALSE) {
	# On input, qtab is a dataframe of quantiles.  Each column corresponds to
	# a sample size, and each row corresponds to a quantile value.  The sample
	# sizes are given in the first row, and the quantiles are given in the
	# first column.  
	n <- nrow(qtab)
	i <- findInterval(sample_size, qtab[1,2:ncol(qtab)])+1
	if (i == 1) {
		parent_name <- as.character(sys.call(-1)[[1]])
		if (stop.on.na) {
			stop (parent_name," requires a minimum of ", qtab[1,2], " observations.")
		} else {
			warning (parent_name, " requires a minimum of ", qtab[1,2], " observations.")
			return(NA)
		}
	}
	y1 <- approx(qtab[2:n, i], qtab[2:n, 1], stat, rule=2)$y
	if (i < ncol(qtab)) {
		y2 <- approx(qtab[2:n, i+1], qtab[2:n, 1], stat, rule=2)$y
		n1 <- qtab[1,i]
		n2 <- qtab[1,i+1]
		y <- y1 * (n2 - sample_size) / (n2 - n1) + y2 * (sample_size - n1)/(n2 - n1)
	} else {
		y <- y1
	}
	y
}
