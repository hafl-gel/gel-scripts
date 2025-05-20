
# source('~/repos/3_Scripts/4_MiniDOASAuswertung/RScripts/miniDOAS_functions_hac5.R')


#### read doas data
read_data <- function(folder, from, to = NULL, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), ncores = 1,
    Serial = NULL, force.write.daily = FALSE, rawdataOnly = missing(path_dailyfiles),
    path_dailyfiles = folder){
    if(length(from) > 1){
        to <- from[2]
        from <- from[1]
    } else if(is.null(to)){
        to <- parse_timerange(from, tz = tz)
        from <- to[1]
        to <- to[2]
    }
    to <- parse_date_time3(to, tz = tz)
    from <- parse_date_time3(from, tz = tz)
    stopifnot(to >= from)
    di <- getDOASinfo(doas, timerange = c(from, to), tzone = tz, Serial = Serial)
    structure(
        readDOASdata(di, folder, rawdataOnly = rawdataOnly, 
            force.write.daily = force.write.daily, ncores = ncores,
            path_dailyfiles = path_dailyfiles
        )
        , class = 'rawdat'
        )
}

#### print specdat method
print.rawdat <- function(x, ...){
    nc <- length(x$RawData)
    wl <- get_wl(x)
    cat('~~~~\n')
    cat('\t', x$DOASinfo$DOASmodel, '/', x$DOASinfo$Spectrometer$Serial, '- raw data\n')
    cat('\t', nc, 'entries\n')
    cat('\t', min(wl), 'to', max(wl), 'nm', sprintf('(%s pixel)\n', length(wl)))
    cat('\t recorded between', format(x$Header$st[1]),'and', format(x$Header$et[nc]),'\n')
    cat('~~~~\n')
}

#### plot raw spectra
plot.rawdat <- function(x, y, sweep_fun = NULL, sweep_stats = median, log = 'y', xlim = c(190, 230), ylim = NULL, 
    col = 'black', lwd = 1, lty = 1, ...) {
    # get wl
    wl <- get_wl(x)
    # get specs
    specs <- x$RawData
    # get subset xlim
    if (is.null(xlim)) {
        xlim <- range(wl)
    } else {
        ind <- which(wl >= xlim[1] & wl <= xlim[2])
        wl <- wl[ind]
        specs <- lapply(specs, function(y) y[ind])
    }
    # calculate difference to median
    if (!is.null(sweep_fun)) {
        log <- ''
        stats <- apply(data.frame(specs), 1, sweep_stats)
        specs <- lapply(specs, sweep_fun, stats)
    }
    # get ylim
    if (is.null(ylim)) ylim <- range(specs)
    # plot empty + lines
    plot(xlim, ylim, type = 'n', log = log, ...)
    for (spec in specs) {
        lines(wl, spec, col = col, lwd = lwd, lty = lty)
    }
}
lines.rawdat <- function(x, y, sweep_fun = NULL, sweep_stats = median, col = 'black', lwd = 1, lty = 1, ...) {
    # get wl
    wl <- get_wl(x)
    # get specs
    specs <- x$RawData
    # calculate difference to median
    if (!is.null(sweep_fun)) {
        log <- ''
        stats <- apply(data.frame(specs), 1, sweep_stats)
        specs <- lapply(specs, sweep_fun, stats)
    }
    for (spec in specs) {
        lines(wl, spec, col = col, lwd = lwd, lty = lty, ...)
    }
}

# split raw data depending on difference in light level
split_raw <- function(rawdat, max_dist = 5e-3, avg_dist = 3e-3, sd_dist = 3e-4) {
    # get specs
    specs <- rawdat$RawData
    # get windows
    wins <- get_wins(rawdat)
    # calculate distances between all
    dist <- outer(seq_along(specs), seq_along(specs), function(i, j) {
        # change to max between 190 and 230
        mapply(function(x, y) {
            # subsets
            # xs <- x[wins$pixel_filter]
            # ys <- y[wins$pixel_filter]
            xs <- x[wins$pixel_fit]
            ys <- y[wins$pixel_fit]
            ds <- abs(xs - ys) / (xs + ys)
            out <- c(max = max(ds), avg = mean(ds), std = sd(ds))
            # set dist > max_dist etc. to NA
            if (any(out > c(max_dist, avg_dist, sd_dist))) {
                NA_real_
            } else {
                out[1]
            }
            }, x = specs[i], y = specs[j], SIMPLIFY = TRUE)
    })
    # dist <- lapply(c(max, mean, sd), \(fun) outer(seq_along(specs), seq_along(specs), function(i, j) {
    #     # change to max between 190 and 230
    #     mapply(function(x, y) {
    #         # subsets
    #         xs <- x[wins$pixel_filter]
    #         ys <- y[wins$pixel_filter]
    #         ds <- abs(xs - ys) / (xs + ys)
    #         fun(ds)
    #         }, x = specs[i], y = specs[j], SIMPLIFY = TRUE)
    # }))
    diag(dist) <- NA
    out <- list()
    while (!all(is.na(dist))) {
        # print(dist)
        # find smallest 
        i_min <- which.min(dist)
        i_row <- ceiling(i_min / nrow(dist))
        i_col <- (i_min - 1) %% nrow(dist) + 1
        # print(i_row)
        # print(i_col)
        # check if can be connected
        row_exists <- i_row %in% unlist(out)
        col_exists <- i_col %in% unlist(out)
        if (row_exists && !col_exists) {
            # which set
            i_set <- which(sapply(out, function(x) i_row %in% x))
            # check i_col
            check_me <- out[[i_set]] 
            # remove i_row
            check_me <- check_me[!(check_me %in% i_row)]
            # don't add if any NA
            if (!anyNA(dist[i_col, check_me])) {
                out[[i_set]] <- c(out[[i_set]], i_col)
            }
            # set all to NA
            dist[i_col, check_me] <- NA
            dist[check_me, i_col] <- NA
        } else if (col_exists && !row_exists) {
            # which set
            i_set <- which(sapply(out, function(x) i_col %in% x))
            # check i_row
            check_me <- out[[i_set]] 
            # remove i_col
            check_me <- check_me[!(check_me %in% i_col)]
            # don't add if any NA
            if (!anyNA(dist[i_row, check_me])) {
                out[[i_set]] <- c(out[[i_set]], i_row)
            }
            # set all to NA
            dist[i_row, check_me] <- NA
            dist[check_me, i_row] <- NA
        } else if (!col_exists && !row_exists) {
            # add new set
            out <- c(out, list(c(i_row, i_col)))
        } else {
            # check if sets can be joined together
            # which sets
            col_set <- which(sapply(out, function(x) i_col %in% x))
            row_set <- which(sapply(out, function(x) i_row %in% x))
            # get set entries
            col_val <- out[[col_set]]
            row_val <- out[[row_set]]
            # anyNA?
            if (!anyNA(dist[col_val, row_val])) {
                # join
                out[[col_set]] <- c(out[[col_set]], out[[row_set]])
                # remove
                out[[row_set]] <- NULL
            }
            # set NA
            dist[col_val, row_val] <- NA
            dist[row_val, col_val] <- NA
        }
        # set to NA
        dist[c(i_row, i_col), c(i_col, i_row)] <- NA
        # print(out)
    }
    # return list
    lapply(out, function(x) {
        filter_index(rawdat, sort(x))
    })
}


#### print avgdat method
print.avgdat <- function(x, ...) {
    has_dc <- attr(x, 'dark.corrected')
    has_sc <- attr(x, 'straylight.corrected')
    has_lc <- attr(x, 'linearity.corrected')
    x <- attr(x, 'RawData')
    nc <- length(x$RawData)
    wl <- get_wl(x)
    cat('~~~~\n')
    cat('\t', x$DOASinfo$DOASmodel, '/', x$DOASinfo$Spectrometer$Serial, '- average spectrum\n')
    cat('\t', nc, 'spectra averaged\n')
    cat('\t', min(wl), 'to', max(wl), 'nm', sprintf('(%s pixel)\n', length(wl)))
    cat('\t recorded between', format(x$Header$st[1]),'and', format(x$Header$et[nc]),'\n')
    cat('\t dark corrected:', has_dc, '\n')
    cat('\t straylight corrected:', has_sc, '\n')
    cat('\t linearity corrected:', has_lc, '\n')
    cat('~~~~\n')
}

#### extract raw spectra
get_specs <- function(folder, from, to = NULL, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.straylight = FALSE, correct.linearity = FALSE) {
    # get rawdata or subset
    if(!inherits(folder, 'rawdat')){
        folder <- read_data(folder, from, to, tz, doas, Serial)
    } else {
        if (is.null(to) && is.character(from) && length(from) == 1) {
            to <- parse_timerange(from, tz = tz)
            from <- to[1]
            to <- to[2]
        }
        # convert from/to to POSIXct
        from <- parse_date_time3(from, tz = tz)
        to <- parse_date_time3(to, tz = tz)
        # get indices
        ind <- which(folder$Header[['et']] > from & folder$Header[['st']] < to)
        if (length(ind)) {
            folder$RawData <- folder$RawData[ind]
            folder$Header <- folder$Header[ind, , drop = FALSE]
        } else {
            stop('No data within specified timerange available')
        }
    }
    # correct for dark current
    if (correct.straylight) {
        folder$RawData <- correct_straylight(folder)
    }
    # correct linearity
    if (correct.linearity) {
        folder$RawData <- correct_linearity(folder)
    }
    structure(
        folder$RawData
        , RawData = folder
        , class = 'specs'
        , straylight.corrected = correct.straylight
        , linearity.corrected = correct.linearity
        )
}

#### extract single spectrum
single_specs <- function(folder, at, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.straylight = TRUE, correct.linearity = TRUE, dark = NULL) {
    # convert at to POSIXct
    at <- parse_date_time3(at, tz = tz)
    # get indices
    ind <- which(folder$Header[['et']] >= at & folder$Header[['st']] <= at)
    if (length(ind)) {
        folder$RawData <- folder$RawData[ind]
        folder$Header <- folder$Header[ind, , drop = FALSE]
        # update timerange
        folder[['DOASinfo']][['timerange']] <- c(folder[['Header']][1, 'st'], folder[['Header']][1, 'et'])
    } else {
        stop('No data at specified time available')
    }
    # correct dark current
    if (is.null(dark)) {
        warning('no dark spectrum provided')
    } else {
        folder$RawData <- correct_dark(folder, dark)
    }
    # correct for dark current
    if (correct.straylight) {
        folder$RawData <- correct_straylight(folder)
    }
    # correct linearity
    if (correct.linearity) {
        folder$RawData <- correct_linearity(folder)
    }
    structure(
        folder$RawData[[1]]
        , RawData = folder
        , class = 'single_spec'
        , straylight.corrected = correct.straylight
        , linearity.corrected = correct.linearity
        , dark.corrected = !is.null(dark)
        )
}

### get spec within wavelength range

### methods for single_spec class
print.single_spec <- function(x, lo = 200, hi = 230, ...) {
    rd <- attr(x, 'RawData')
    xi <- cut_wl(x, lo, hi)
    cat('***\nSingle spectrum:\n')
    cat('   recorded between', format(rd$Header[['st']]), 'and', format(rd$Header[['et']]), '\n')
    cat('\t dark corrected:', attr(x, 'dark.corrected'), '\n')
    cat('\t straylight corrected:', attr(x, 'straylight.corrected'), '\n')
    cat('\t linearity corrected:', attr(x, 'linearity.corrected'), '\n')
    cat(sprintf('   I (%i to %i nm) min/avg/max: %1.0f/%1.0f/%1.0f\n***\n', lo, hi, min(xi, na.rm = TRUE), mean(xi, na.rm = TRUE), max(xi, na.rm = TRUE)))
}
plot.single_spec <- function(x, y, ylab = 'counts', xlab = 'nm', xlim = c(190, 230), ...) {
    if(is.null(xlim)) xlim <- range(get_wl(x))
    x <- cut_wl(x, xlim[1], xlim[2])
    y <- x
    class(y) <- 'numeric'
    x <- get_wl(x)
    plot(x, y, xlab = xlab, ylab = ylab, ...)
}
lines.single_spec <- function(x, ...) {
    y <- x
    class(y) <- 'numeric'
    x <- get_wl(x)
    lines(x, y, ...)
}
points.single_spec <- function(x, ...) {
    y <- x
    class(y) <- 'numeric'
    x <- get_wl(x)
    points(x, y, ...)
}

# helper functions to process calibration data
save_calref <- function(x, qs_preset = c('high', 'archive')[2], path_rsaves = 'rsaves') {
    # derive name
    tr <- as.POSIXct(range(lapply(x, function(y) get_timerange(y[['dc']]))), origin = '1970-01-01', tz = 'Etc/GMT-1')
    name <- paste0('miniDOAS-ref_cal_spec-',
        paste(format(tr[1], '%y%m%d'), paste(format(tr, '%H%M'), collapse = '_'), sep = '-')
        , '-', attr(x[['nh3']][['dc']], 'meas')$DOASinfo$DOASmodel, '.qs')
    # name <- deparse(substitute(x))
    qsave(x, file.path(path_rsaves, name), preset = qs_preset[1])
}
create_calref <- function(spec_set) {
    out <- list()
    out[['nh3']][['cal_spec']] <- read_cal(spec_set, 'dat.NH3', dont_warn_dark = TRUE)
    out[['nh3']][['ref_spec']] <- read_cal(spec_set, 'dat.N2.NH3', dont_warn_dark = TRUE)
    out[['no']][['cal_spec']] <- read_cal(spec_set, 'dat.NO', dont_warn_dark = TRUE)
    out[['no']][['ref_spec']] <- read_cal(spec_set, 'dat.N2.NO', dont_warn_dark = TRUE)
    out[['so2']][['cal_spec']] <- read_cal(spec_set, 'dat.SO2', dont_warn_dark = TRUE)
    out[['so2']][['ref_spec']] <- read_cal(spec_set, 'dat.N2.SO2', dont_warn_dark = TRUE)
    out[['nh3']][['dc']] <- calc_dc(out[['nh3']][['cal_spec']], out[['nh3']][['ref_spec']])
    out[['no']][['dc']] <- calc_dc(out[['no']][['cal_spec']], out[['no']][['ref_spec']])
    out[['so2']][['dc']] <- calc_dc(out[['so2']][['cal_spec']], out[['so2']][['ref_spec']])
    structure(out, class = 'calref')
}
read_calref <- function(...) {
    create_calref(getSpecSet(...))
}

process_callist <- function(callist, all = 1, nh3 = all, no = all, so2 = all, 
    n2 = list(nh3 = all, no = all, so2 = all), cuvette.length = 0.075, 
    cuvette.mgm3 = list(nh3 = 193.4095, no = 593.9938, so2 = 76.29128),
    dark_spec = NULL) {
    # capture n2 arg
    n2_args <- list(nh3 = all, no = all, so2 = all)
    if (is.list(n2)) {
        n2 <- n2[names(n2) %in% names(n2_args)]
        if (length(n2)) {
            for (nm in names(n2)) {
                n2_args[[nm]] <- n2[[nm]]
            }
        }
    } else {
        n2_args <- setNames(rep(list(n2), 3), names(n2_args))
    }
    # capture args in list
    args <- list(nh3 = nh3, no = no, so2 = so2, n2 = n2_args)
    # use lowercase names for callist
    nms <- names(callist) <- tolower(names(callist))
    # loop over cals
    cal_read <- setNames(lapply(nms, function(x) {
        # read cal
        if (x == 'n2') {
            # get n2 args and names
            n2_args <- args[[x]]
            n2_nms <- names(n2_args)
            # only select references for provided gases
            n2_nms <- n2_nms[n2_nms %in% nms]
            # loop over n2 entries
            setNames(lapply(n2_nms, function(y) {
                # average spectra?
                if (length(n2_args[[y]]) != 1) {
                    read_cal(.avg_avgs(callist[[x]][['avgs']][n2_args[[y]]]), dark = dark_spec)
                } else {
                    read_cal(callist[[x]][['avgs']][[n2_args[[y]]]], dark = dark_spec)
                }
            }), n2_nms)
        } else {
            # average spectra?
            if (length(args[[x]]) != 1) {
                read_cal(.avg_avgs(callist[[x]][['avgs']][args[[x]]]), dark = dark_spec)
            } else {
                read_cal(callist[[x]][['avgs']][[args[[x]]]], dark = dark_spec)
            }
        }
    }), nms)  
    # correct cuvette concentration
    if (!missing(cuvette.mgm3)) {
        if (isTRUE(cuvette.mgm3) || 
            (isTRUE(is.character(cuvette.mgm3)) && cuvette.mgm3 == 'revolver')
        ) {
            cuvette.mgm3 <- eval(formals(process_callist)$cuvette.mgm3)
        }
        for (nm in names(cuvette.mgm3)) {
            # fix conc
            cal_read[[nm]]$Calinfo$cuvette.conc <- cuvette.mgm3[[nm]]
            ind <- grep('cuvette conc', cal_read[[nm]]$calref.info)
            cal_read[[nm]]$calref.info[ind] <- sub('([0-9.]*)$', cuvette.mgm3[[nm]], cal_read[[nm]]$calref.info[ind])
            # fix gas
            cal_read[[nm]]$Calinfo$cuvette.gas <- toupper(nm)
            ind <- grep('cuvette gas', cal_read[[nm]]$calref.info)
            cal_read[[nm]]$calref.info[ind] <- sub(':.*$', paste(':', toupper(nm)), cal_read[[nm]]$calref.info[ind])
            # fix cuvette length
            cal_read[[nm]]$Calinfo$cuvette.path <- cuvette.length
        }
    }
    # assamble output
    out <- list()
    if ('nh3' %in% nms) {
        out[['nh3']][['cal_spec']] <- cal_read[['nh3']]
        out[['nh3']][['ref_spec']] <- cal_read[['n2']][['nh3']]
        out[['nh3']][['dc']] <- calc_dc(out[['nh3']][['cal_spec']], out[['nh3']][['ref_spec']])
    }
    if ('no' %in% nms) {
        out[['no']][['cal_spec']] <- cal_read[['no']]
        out[['no']][['ref_spec']] <- cal_read[['n2']][['no']]
        out[['no']][['dc']] <- calc_dc(out[['no']][['cal_spec']], out[['no']][['ref_spec']])
    }
    if ('so2' %in% nms) {
        out[['so2']][['cal_spec']] <- cal_read[['so2']]
        out[['so2']][['ref_spec']] <- cal_read[['n2']][['so2']]
        out[['so2']][['dc']] <- calc_dc(out[['so2']][['cal_spec']], out[['so2']][['ref_spec']])
    }
    structure(out, class = 'calref')
}
.avg_avgs <- function(...) {
    dots <- list(...)
    if (length(dots) == 1 && is.list(dots) && is.null(names(dots))) {
        dots <- dots[[1]]
    }
    rawdat_list <- lapply(dots, attr, 'RawData')
    if (uniqueN(sapply(rawdat_list, \(x) x$DOASinfo$DOASmodel)) != 1) {
        stop('data is not from one single minidoas!')
    }
    if (sum(sapply(dots, \(x) attr(x, 'straylight.corrected'))) %in% c(1, 2)) {
        stop('straylight correction is inconsistent!')
    }
    if (sum(sapply(dots, \(x) attr(x, 'linearity.corrected'))) %in% c(1, 2)) {
        stop('linearity correction is inconsistent!')
    }
    if (sum(sapply(dots, \(x) attr(x, 'dark.corrected'))) %in% c(1, 2)) {
        stop('dark correction is inconsistent!')
    }
    rd <- do.call('c', lapply(rawdat_list, '[[', 'RawData'))
    hd <- do.call('rbind', lapply(rawdat_list, '[[', 'Header'))
    di <- rawdat_list[[1]]$DOASinfo
    di$timerange <- range(hd$st, hd$et)
    RawData <- structure(list(
        RawData = rd,
        Header = hd,
        DOASinfo = di
    ), class = 'rawdat')
    structure(
        suppressWarnings(avg_spec(RawData)),
        RawData = RawData,
        straylight.corrected = FALSE,
        linearity.corrected = FALSE,
        dark.corrected = FALSE
    )
}

plot.calref <- function(x, add_cheng = TRUE, per_molecule = TRUE, log = '', save.path = NULL, 
    scale_cheng = 1, ylim = c('fix', 'free')[1], dc_grid = TRUE, ...) {
    # save figure?
    if (!is.null(save.path)) {
        # derive figure name
        tr <- as.POSIXct(range(lapply(x, function(y) get_timerange(y[['dc']]))), origin = '1970-01-01', tz = 'Etc/GMT-1')
        name <- paste0('miniDOAS-', attr(x[['nh3']][['dc']], 'meas')$DOASinfo$DOASmodel,
            '-ref_cal_spec-', paste(
                format(tr[1], '%y%m%d'), 
                paste(format(tr, '%H%M'), collapse = '_')
                , sep = '-'
            ) , '.jpg'
        )
        # be verbose
        cat('saving figure to', file.path(save.path, name), '\n')
        height <- 480 * 1.5
        jpeg(file.path(save.path, name), width = height * 1.4, height = height)
            plot(x, add_cheng = add_cheng, per_molecule = per_molecule, log = log, save.path = NULL, ...)
        dev.off()
    }
    if (dc_grid) {
        panel_dc <- function() {
            grid()
            abline(h = 0, col = 'lightgrey')
        }
    } else {
        panel_dc <- function() NULL
    }
    ylim_arg <- ylim
    par(mfrow = c(3, 3))
    for (i in seq_along(x)) {
        # plot cal spec
        plot(x[[i]][['cal_spec']], log = log, ...)
        legend('topleft', bty = 'n', legend = paste0(
                names(x)[i], ' cal spec, measured between:\n',
                deparse_timerange(get_timerange(x[[i]][['cal_spec']]), sep = ' and ')) )
        # plot ref spec
        plot(x[[i]][['ref_spec']], type = 'n', log = log, ...)
        lines(x[[i]][['cal_spec']], col = 'darkgrey')
        lines(x[[i]][['ref_spec']], col = 'black')
        legend('topleft', bty = 'n', legend = paste0(
                names(x)[i], ' ref spec, measured between:\n',
                deparse_timerange(get_timerange(x[[i]][['ref_spec']]), sep = ' and ')) )
        # plot dc
        if (add_cheng && names(x)[i] == 'nh3') {
            if (is.na(x[['nh3']][['cal_spec']]$Calinfo$cuvette.conc)) {
                stop('NH3 calibration concentration is not specified!\n',
                    'Most likely, the calibration was done without the cuvette revolver being installed!')
            }
            cheng <- find_cheng(x[['nh3']][['dc']], show = FALSE, return.cheng.dc = TRUE, mgm3 = x[['nh3']][['cal_spec']]$Calinfo$cuvette.conc)
            s_cheng <- dc2sigma(cheng$cheng, copy = TRUE)
            # scale cheng
            s_cheng$cnt <- s_cheng$cnt * scale_cheng
            s_dc <- dc2sigma(x[['nh3']][['dc']], copy = TRUE)
            if (isTRUE(is.character(ylim[1]))) {
                ylim <- switch(ylim_arg[1]
                    , free = 
                    , range = {
                        range(c(s_cheng$cnt, s_dc$cnt), na.rm = TRUE)
                    }
                    , fix = 
                    , fixed = 
                    , given = {
                        c(-3.5, 2) * 1e-18
                    }
                    , stop('argument ylim is not valid!')
                )
            } else if (is.null(ylim_arg)) {
                ylim <- range(c(s_cheng$cnt, s_dc$cnt), na.rm = TRUE)
            }
            plot(x[[i]][['dc']], per_molecule = per_molecule, type = 'n', ylim = ylim, 
                panel.first = panel_dc(), ...)
            lines(cheng$cheng, fctr = scale_cheng, col = 'indianred', per_molecule = per_molecule)
            lines(x[[i]][['dc']], per_molecule = per_molecule, col = 'black')
            legend('bottomright', legend = sprintf('span = %1.3f (+/- %1.3f)\nshift = %1.2f', 
                    cheng$coefs[2], cheng$se[2], cheng$shift), bty = 'n', inset = 0.05)
        } else {
            if (isTRUE(is.character(ylim_arg[1]))) {
                ylim <- switch(ylim_arg[1]
                    , free = 
                    , range = {
                        range(s_dc$cnt, na.rm = TRUE)
                    }
                    , fix = 
                    , fixed = 
                    , given = {
                        switch(names(x)[i]
                            , no = c(-1.8, 0.6) * 1e-18
                            , so2 = c(-4.8, 3) * 1e-18
                        )
                    }
                    , stop('argument ylim is not valid!')
                )
            } else if (is.null(ylim_arg)) {
                ylim <- range(s_dc$cnt, na.rm = TRUE)
            }
            plot(x[[i]][['dc']], per_molecule = per_molecule, ylim = ylim, 
                panel.first = panel_dc(), ...)
        }
    }
}

ppm_mgm3 <- function(ppm, T_deg, p_hPa, molar_mass = 17) {
    ppm * molar_mass * p_hPa / 10 / (8.3144598 * (T_deg + 273.15))
}
mgm3_ppm <- function(mgm3, T_deg, p_hPa, molar_mass = 17) {
    mgm3 * (8.3144598 * (T_deg + 273.15)) * 10 / p_hPa / molar_mass
}


#### process calref rawdata
read_gas <- function(gas, path_data, from, show = TRUE, max.dist = 5e-3, 
    avg.dist = 3e-3, sd.dist = 3e-4, min.num = 2) {
    if (inherits(path_data, 'rawdat')) {
        # pass on
        rawdata <- path_data
    } else {
        # read data from timerange
        rawdata <- read_data(path_data, from = from)
    }
    if (is.character(gas) && tolower(gas[1]) %in% c('nh3', 'so2', 'no', 'n2')) {
        # filter for revolver position
        raw <- filter_position(rawdata, gas)
        time_filtered <- FALSE
    } else if (is.character(gas) && tolower(gas[1]) %in% 'closed') {
        # filter for shutter closed
        raw <- filter_closed(rawdata)
        time_filtered <- FALSE
    } else {
        day <- format(get_timerange(rawdata)[1], '%d.%m.%Y ')
        # split times
        split_times <- strsplit(gas, split = paste(getOption('time.separators'), collapse = '|'))
        # filter by time
        times <- lapply(split_times, function(x) {
            # trim whitespace
            x <- trimws(x)
            # get days
            days <- sub('(.*)[0-9]{2}:[0-9]{2}(:[0-9]{2})?$', '\\1', x)
            days[days == ''] <- day
            # get times
            times <- sub('(.*)([0-9]{2}:[0-9]{2}(:[0-9]{2})?)$', '\\2', x)
            # paste together again
            x <- paste0(days, times)
            format(ibts::parse_timerange(x, tz = 'Etc/GMT-1'))
        })
        raw <- filter_time(rawdata, sapply(times, '[', 1), sapply(times, '[', 2))
        gas <- paste(gas, collapse = '\n')
        time_filtered <- TRUE
    }
    # split
    sets <- split_raw(raw, max.dist, avg.dist, sd.dist)
    # remove sets with less than min.num
    sets <- sets[sapply(sets, function(x) length(x$RawData) >= min.num)]
    if (show) {
        # get color palette
        require(colorspace)
        cpal <- qualitative_hcl(length(sets))
        # sweep before
        med <- apply(data.frame(raw$RawData), 1, median)
        swept <- raw
        swept$RawData <- lapply(raw$RawData, '/', med)
        # raw
        x11()
        par(mfrow = c(2, 1))
        plot(raw, main = paste(gas, '- unfiltered'))
        # use different colors for different sets (use color palette also in set figures)
        plot(swept, ylim = 1 + c(-1, 1) * max(abs(quantile(unlist(swept$RawData), c(0.025, 0.975)) - 1)) * 2, 
            col = 'lightgrey')
        # loop over sets
        for (s in seq_along(sets)) {
            lines(filter_index(swept, which(names(swept$RawData) %in% names(sets[[s]]$RawData))), col = cpal[s])
        }
        # sets
        lapply(seq_along(sets), function(i) {
            x11()
            par(mfrow = c(2, 1))
            plot(sets[[i]], main = paste(
                    if (time_filtered) paste(sets[[i]]$Header$st[1], 
                        sets[[i]]$Header$et[length(sets[[i]]$RawData)], sep = ' to ') else gas, 
                    '- set', i, ':', length(sets[[i]]$RawData), 'specs - ', 
                    sprintf('%1.0f', max(sapply(sets[[i]]$RawData, max, na.rm = TRUE)))), 
                col = cpal[i])
            plot(sets[[i]], sweep_fun = '/', col = '#00000022', ylim = 1 + c(-1, 1) * 
            max.dist)
            })
    }
    # average spectra
    suppressWarnings(avgs <- lapply(sets, avg_spec, correct.straylight = FALSE,
            correct.linearity = FALSE))
    # get i_max
    i_max <- sapply(avgs, max)
    # return list
    list(
        raw = raw,
        sets = sets,
        avgs = avgs,
        i_max = i_max
        )
}

read_all_gases <- function(path_data, timerange, show = TRUE, 
    gases = c('N2', 'NH3', 'NO', 'SO2'), max.dist = 5e-3, avg.dist = 3e-3,
    sd.dist = 3e-4, min.num = 2) {
    if (inherits(path_data, 'rawdat')) {
        # pass on
        rawdata <- path_data
    } else {
        # read data from timerange
        rawdata <- read_data(path_data, from = timerange)
    }
    # check max.dist
    md <- list(nh3 = 5e-3, no = 5e-3, so2 = 5e-3, n2 = 5e-3, closed = 1e-3)
    if (is.list(max.dist)) {
        max.dist <- max.dist[names(max.dist) %in% names(md)]
        if (length(max.dist)) {
            for (nm in names(max.dist)) {
                md[[nm]] <- max.dist[[nm]]
            }
        }
    } else if (!missing(max.dist) && is.numeric(max.dist)) {
        md <- setNames(rep(list(max.dist), length(md)), names(md))
    }
    ad <- list(nh3 = 3e-3, no = 3e-3, so2 = 3e-3, n2 = 3e-3, closed = 3e-3/5)
    if (is.list(avg.dist)) {
        avg.dist <- avg.dist[names(avg.dist) %in% names(ad)]
        if (length(avg.dist)) {
            for (nm in names(avg.dist)) {
                ad[[nm]] <- avg.dist[[nm]]
            }
        }
    } else if (!missing(avg.dist) && is.numeric(avg.dist)) {
        ad <- setNames(rep(list(avg.dist), length(ad)), names(ad))
    }
    sd <- list(nh3 = 3e-4, no = 3e-4, so2 = 3e-4, n2 = 3e-4, closed = 3e-4/5)
    if (is.list(sd.dist)) {
        sd.dist <- sd.dist[names(sd.dist) %in% names(sd)]
        if (length(sd.dist)) {
            for (nm in names(sd.dist)) {
                sd[[nm]] <- sd.dist[[nm]]
            }
        }
    } else if (!missing(sd.dist) && is.numeric(sd.dist)) {
        sd <- setNames(rep(list(sd.dist), length(sd)), names(sd))
    }
    # loop over gases
    if (is.list(gases)) {
        setNames(mapply(read_gas, gases, max.dist = md[tolower(names(gases))], 
                avg.dist = ad[tolower(names(gases))],
                sd.dist = sd[tolower(names(gases))],
                MoreArgs = list(path_data = rawdata, show = show, min.num = min.num), 
                SIMPLIFY = FALSE), names(gases))
    } else {
        setNames(mapply(read_gas, gases, max.dist = md[tolower(gases)], 
                avg.dist = ad[tolower(gases)],
                sd.dist = sd[tolower(gases)],
                MoreArgs = list(path_data = rawdata, show = show, min.num = min.num), 
                SIMPLIFY = FALSE), gases)
    }
}
read_shutter_closed <- function(path_data, timerange, gases = 'closed', ...) {
    read_all_gases(path_data, timerange, gases = gases, ...)[['closed']]
}



    #### average raw data
avg_spec <- function(folder, from = NULL, to = NULL, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.straylight = TRUE, correct.linearity = TRUE, dark = NULL) {
    if(inherits(folder, 'rawdat')){
        if (!is.null(from)) {
            if (is.null(to) && is.character(from) && length(from) == 1) {
                to <- parse_timerange(from, tz = tz)
                from <- to[1]
                to <- to[2]
            }
            # convert from/to to POSIXct
            from <- parse_date_time3(from, tz = tz)
            to <- parse_date_time3(to, tz = tz)
            # get indices
            ind <- which(folder$Header[['et']] > from & folder$Header[['st']] < to)
        } else {
            ind <- seq_along(folder$Header$et)
        }
        if (length(ind)) {
            folder$RawData <- folder$RawData[ind]
            folder$Header <- folder$Header[ind, , drop = FALSE]
            # update timerange
            folder[['DOASinfo']][['timerange']] <- c(folder[['Header']][1, 'st'], folder[['Header']][length(ind), 'et'])
        } else {
            stop('No data within specified timerange available')
        }
    } else if (is.list(folder) && 'SpecAvg' %in% names(folder)) {
        folder <- structure(
            list(
                RawData = folder[['Specs']],
                Header = folder[['Info']],
                DOASinfo = folder[['DOASinfo']]
            ), class = 'rawdat'
        )
    } else {
        folder <- read_data(folder, from, to, tz, doas, Serial)
    }
    # correct for dark current
    if (is.null(dark)) {
        warning('no dark spectrum provided')
    } else {
        folder$RawData <- correct_dark(folder, dark)
    }
    # correct for straylight
    if (correct.straylight) {
        folder$RawData <- correct_straylight(folder)
    }
    # correct linearity
    if (correct.linearity) {
        folder$RawData <- correct_linearity(folder)
    }
    structure(
        rowMeans(as.data.frame(folder$RawData))
        , RawData = folder
        , class = 'avgdat'
        , straylight.corrected = correct.straylight
        , linearity.corrected = correct.linearity
        , dark.corrected = !is.null(dark)
        )
}

#### helper function to cut data to correct wvl
cut_wl <- function(x, lo = 190, hi = 230) {
    wl <- get_wl(x)
    ind <- which(wl >= lo & wl <= hi)
    if (inherits(x, 'rawdat')) {
        out <- x
        out$RawData <- lapply(x$RawData, function(y) y[ind])
        out$DOASinfo$Spectrometer$pixel <- x$DOASinfo$Spectrometer$pixel[ind]
        out$DOASinfo$Spectrometer$wavelength <- wl[ind]
    } else if (inherits(x, 'avgdat')) {
        out <- x[ind]
        attr(out, 'RawData') <- cut_wl(attr(x, 'RawData'), lo, hi)
        class(out) <- 'avgdat'
    } else if (inherits(x, 'single_spec')) {
        out <- x[ind]
        attr(out, 'RawData') <- cut_wl(attr(x, 'RawData'), lo, hi)
        class(out) <- 'single_spec'
    } else if (inherits(x, 'caldat') || inherits(x, 'chen')) {
        out <- x
        out$data <- x$data[ind, ]
        stop('fix caldat method')
    }
    out
}

#### helper function to get wl
get_wl <- function(x) {
    if (inherits(x, 'rawdat')) {
        x$DOASinfo$Spectrometer$wavelength
    } else if (inherits(x, 'avgdat') || inherits(x, 'single_spec')) {
        attr(x, 'RawData')$DOASinfo$Spectrometer$wavelength
    } else if (inherits(x, 'caldat') || inherits(x, 'chen')) {
        x$data[, wl]
    } else if (inherits(x, 'dc')) {
        attr(x, 'meas')$DOASinfo$Spectrometer$wavelength
    }
}

#### plot method avgdat (use y as ref spec, types, xlim)
plot.avgdat <- function(x, y = NULL, what = c('average', 'residuals', 'originals'), xlim = c(190, 230), 
    ylim = NULL, type = 'l', main = NULL, ylab = NULL, log = '', ...) {
    # cut to xlim
    lo <- xlim[1] - diff(xlim) * 0.04
    hi <- xlim[2] + diff(xlim) * 0.04
    x <- cut_wl(x, lo = lo, hi = hi)
    # switch plot kind
    switch(pmatch(what[1], c('average', 'residuals', 'originals'), nomatch = 4)
        # average
        , yp <- x
        # residuals
        , yp <- sweep(attr(x, 'RawData')$RawData, 1, x)
        # originals
        , yp <- attr(x, 'RawData')[['RawData']]
        # nonsense
        , stop('argument "what" is not valid')
        )
    # prep
    if (is.null(ylab)) ylab <- 'counts'
    if (is.null(main)) main <- paste(attr(x, 'RawData')$DOASinfo$DOASmodel, 
        attr(x, 'RawData')$DOASinfo$Spectrometer$Serial, sep = ' / ')
    xp <- get_wl(x)
    if (is.null(ylim)) ylim <- range(yp, na.rm = TRUE)
    if (NCOL(yp) > 1) {
        plot(xp, x, type = 'n', xlim = xlim, ylim = ylim, xlab = 'nm', ylab = ylab, main = main, log = log, ...)
        for (i in seq_len(ncol(yp))) {
            lines(xp, yp[, i], type = type, ...)
        }
    } else {
        plot(xp, yp, xlim = xlim, ylim = ylim, type = type, xlab = 'nm', ylab = ylab, main = main, log = log, ...)
    }
}

#### lines method avgdat
lines.avgdat <- function(x, y = NULL, what = c('average', 'residuals', 'originals'), ...) {
    # cut to xlim
    usr <- par('usr')[1:2]
    x <- cut_wl(x, lo = usr[1], hi = usr[2])
    # switch plot kind
    switch(pmatch(what[1], c('average', 'residuals', 'originals'), nomatch = 4)
        # average
        , yp <- x
        # residuals
        , yp <- sweep(attr(x, 'RawData')$RawData, 1, x)
        # originals
        , yp <- attr(x, 'RawData')[['RawData']]
        # nonsense
        , stop('argument "what" is not valid')
        )
    # prep
    xp <- get_wl(x)
    if (NCOL(yp) > 1) {
        for (i in seq_len(ncol(yp))) {
            lines(xp, yp[, i], ...)
        }
    } else {
        lines(xp, yp, ...)
    }
}

#### helper function to correct for dark current
correct_dark <- function(x, y) {
    lapply(x$RawData, '-', y$data[, cnt])
}

correct_straylight <- function(x) {
    win <- getWindows(x$DOASinfo)
    # if (is.data.frame(x$RawData)) {
    #     sweep(x$RawData, 2, colMeans(x$RawData[win$pixel_straylight, , drop = FALSE]))
    # } else {
        lapply(x$RawData, '-', sum(sapply(x$RawData, function(y) mean(y[win$pixel_straylight]))) / length(x$RawData))
    # }
}

#### helper function to correct for non-linearity
correct_linearity <- function(x) {
    lin.coef <- x$DOASinfo$Spectrometer$"Linearity Coefficients"
    # x$RawData / linearity.func(x$RawData, lin.coef)
    lapply(x$RawData, function(y) y / linearity.func(y, lin.coef))
}


#### get Cheng 2006 cross-section
get_Cheng <- function() {
    structure(
        list(
            data = attr(get_Cheng, 'cheng_data')
        ),
        class = 'chen'
    )
}

#### read calibration spectrum
read_cal <- function(file, spec = NULL, tz = 'Etc/GMT-1', Serial = NULL, is_dark = FALSE,
    correct.straylight = !is_dark, correct.linearity = !is_dark, lin_before_dark = FALSE, 
    dark = if (is_dark) is_dark else NULL, dont_warn_dark = FALSE) {
    if (is.character(file)) {
        cal <- fread(file, sep = '\n', header = FALSE)
        if (nrow(cal) == 1069) {
            info <- cal[, {
                x <- strsplit(grep(': ', V1, value = TRUE) , split = ': ')
                list(
                    var = sapply(x, '[', 1),
                    val = sapply(x, function(y) paste(y[-1], collapse = ': '))
                    )
            }]
            tracer <- info[7, val]
            out <- list(
                data = cal[, fread(text = V1[26:1069], col.names = c('wl', 'cnt'))],
                Calinfo = list(
                    info = info,
                    spec.name = NULL,
                    cuvette.gas = tracer,
                    cuvette.conc = if (tracer == 'N2') '-' else info[8, as.numeric(val)],
                    cuvette.path = info[9, as.numeric(val)]
                    ),
                DOASinfo = getDOASinfo(info[2, val], timerange = info[1, val], tzone = tz, Serial = Serial),
                calref.info = cal[1:10, V1]
            )
        } else if (nrow(cal) == 1054) {
            info <- cal[, {
                x <- strsplit(grep(': ', V1, value = TRUE) , split = ': ')
                list(
                    var = sapply(x, '[', 1),
                    val = sapply(x, function(y) paste(y[-1], collapse = ': '))
                    )
            }]
            DOASinfo <- getDOASinfo(sub('.*(S[1-6]).*', '\\1', file), timerange = info[1, fast_strptime(val, '%Y%m%d%H%M', tz = tz, lt = FALSE)], tzone = tz, Serial = Serial)
            tracer <- info[3, sub('.*(NH3|SO2|NO|N2).*', '\\1', val)]
            out <- list(
                data = cal[11:1054, .(wl = DOASinfo$Spectrometer$wavelength, cnt = as.numeric(V1))],
                Calinfo = list(
                    info = info,
                    spec.name = NULL,
                    cuvette.gas = tracer,
                    cuvette.conc = if (tracer == 'N2') '-' else info[3, convert2mg(val, tracer)],
                    cuvette.path = info[4, as.numeric(val)]
                    ),
                DOASinfo = DOASinfo,
                calref.info = cal[1:10, V1]
            )
        }
    } else if (inherits(file, 'avgdat')) {
        DOASinfo <- attr(file, 'RawData')$DOASinfo
        Header <- attr(file, 'RawData')$Header
        gas <- unique(Header[, 'RevPos'])
        if (length(gas) > 1) stop('data has more than one unique entry in revolver position!')
        cuvetteConc_mg <- switch(gas, NH3 = 193.4095, NO = 593.9938, SO2 = 76.29128, NA)
        cuvetteLength <- 0.075
        i_max <- sapply(attr(file, 'RawData')$RawData, max, na.rm = TRUE)
        Imax <- list(range = range(i_max), mean = mean(i_max))
        txt <- c(
            "averaged raw spectra file"
            ,paste0("calculation performed on: ",format(Now <- Sys.time()))
            ,paste0("miniDOAS: ",DOASinfo$DOASmodel)
            ,paste0("spectrometer: ",DOASinfo$Spectrometer$"Spectrometer Name")
            ,paste0("revolver position: ",paste(unique(Header[,"RevPos"]),sep=","))
            ,paste0("shutter position: ",paste(unique(Header[,"ShuPos"]),sep=","))
            ,"~~~ cuvette measurement ~~~"
            ,paste0("cuvette installed: ",ifelse((gas=="none"|is.na(gas)),"no","yes"))
            ,paste0("cuvette gas: ",gas)
            ,paste0("cuvette concentration (mg/m3): ",cuvetteConc_mg)
            ,paste0("cuvette length: ",cuvetteLength)
            ,"~~~ path measurement ~~~"
            ,paste0("path length: ",NA)
            ,paste0("ambient NH3 concentration (ug/m3): ",NA)
            ,paste0("ambient SO2 concentration (ug/m3): ",NA)
            ,paste0("ambient NO concentration (ug/m3): ",NA)
            ,"~~~ "
            ,paste0("averaging period: ",format(Header[1,"st"])," to ",format(rev(Header[,"et"])[1]))
            ,paste0("averaged ",length(Header[,"et"])," spectra")
            ,paste0("number of accumulations: ",paste(range(Header[, 'AccNum']),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(Header[, 'AccNum'])),")")
            ,paste0("exposure time (ms): ",paste(range(Header[, 'Expos']),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(Header[, 'Expos'])),")")
            ,paste0("TEC temperature (deg C): ",paste(range(as.numeric(Header[, 'TECTemp'])),collapse=" to ")," (avg: ",sprintf("%1.1f",mean(as.numeric(Header[, 'TECTemp']))),")")
            ,paste0("ambient temperature (deg C): ",paste(NA,collapse=" to ")," (avg: ",sprintf("%1.1f",NA),")")
            ,paste0("Imax: ",sprintf("%1.1f to %1.1f",Imax$range[1],Imax$range[2])," (avg: ",sprintf("%1.1f",Imax$mean),")")
            ,"-----")
        out <- list(
            data = data.table(wl = get_wl(file), cnt = as.numeric(file)),
            Calinfo = list(
                info = data.table(var = 'timerange', val = paste0(DOASinfo$timerange, collapse = ' and ')),
                spec.name = spec,
                cuvette.gas = gas,
                cuvette.conc = cuvetteConc_mg,
                cuvette.path = cuvetteLength
                ),
            DOASinfo = DOASinfo,
            calref.info = txt
            )
        # check corrections:
        if (any(attr(file, 'straylight.corrected'), attr(file, 'linearity.corrected'),
                attr(file, 'dark.corrected'))) {
            # TODO: set corrections to FALSE
            # ...
            # pass on starylight/dark values
            stop('Averaged spectra have been corrected before. Fix passing on of straylight and dark values!')
        }
    } else {
        # if file is result from getSpecSet
        if (is.null(spec)) {
            stop('argument spec is unset. spec defines which SpecSet list entry will be selected')
            # valid names:
            # "ref.spec"
            # "ref.dark.spec"
            # "dark.spec"
            # "NH3.cal.spec"
            # "SO2.cal.spec"
            # "NO.cal.spec"
            # "N2.NH3.cal.spec"
            # "N2.SO2.cal.spec"
            # "N2.NO.cal.spec"
            # "N2.dark.cal.spec"
        }
        if (!(spec %in% names(file)[-1])) {
            stop('Argument spec - recognized names are:', 
                sprintf('    %s\n', names(file)[-1])
                )
        }
        if (spec %in% c('dat.ref.dark', 'dat.dark', 'dat.N2.dark')) is_dark <- TRUE
        # get entry
        spec_out <- file[[spec]]
        # get doas info
        DOASinfo <- getDOASinfo(file[[1]], timerange = spec_out$timerange)
        out <- list(
            data = data.table(wl = DOASinfo$Spectrometer$wavelength, cnt = as.numeric(spec_out[['dat.spec']])),
            Calinfo = list(
                info = data.table(var = 'timerange', val = paste0(spec_out$timerange, collapse = ' and ')),
                spec.name = spec,
                cuvette.gas = spec_out[['cuvette']][['cuvetteGas']],
                cuvette.conc = spec_out[['cuvette']][['cuvetteConc_mg']],
                cuvette.path = spec_out[['cuvette']][['cuvetteLength']]
                ),
            DOASinfo = DOASinfo,
            calref.info = spec_out$info.spec
            )
    }
    # save original
    out$Calinfo$raw.spec <- out$data[, cnt]
    # correct dark
    out$Calinfo$dark.corrected = FALSE
    if (!is_dark && is.null(dark)) {
        if (!dont_warn_dark) warning('no dark spectrum provided')
    } else if(!is_dark) {
        if (all(c('type', 'tracer', 'dir', 'timerange') %in% names(dark))) {
            dark <- getSpec(dark, DOASmodel = out$DOASinfo$DOASmodel, SpecName = 'dark.spec')
        }
        if (all(c('SpecAvg', 'Info', 'Specs', 'DOASinfo', 'txt') %in% names(dark))) {
            cdark <- dark$SpecAvg
        } else if (all(c('dat.spec', 'timerange', 'info.spec') %in% names(dark))) {
            cdark <- dark$dat.spec
        } else {
            cdark <- dark$data[, cnt]
        }
        out$data[, cnt := cnt - cdark]
        out$Calinfo$dark.corrected <- TRUE
        out$Calinfo$dark.spec <- cdark
    }
    # correct straylight
    out$Calinfo$straylight.corrected = FALSE
    if (correct.straylight) {
        win <- getWindows(out$DOASinfo)
        dark <- out$data[, mean(cnt[win$pixel_straylight])]
        out$Calinfo$straylight.value <- dark
        if (!lin_before_dark) {
            out$data[, cnt := cnt - dark]
            out$Calinfo$straylight.corrected <- TRUE
        }
    }
    # correct linearity
    out$Calinfo$linearity.corrected = FALSE
    if (correct.linearity) {
        lin.coef <- out$DOASinfo$Spectrometer$'Linearity Coefficients'
        out$data[, cnt := cnt / linearity.func(cnt, lin.coef)]
        out$Calinfo$linearity.corrected <- TRUE
    }
    # correct.straylight afterwards
    if (correct.straylight && lin_before_dark) {
        out$data[, cnt := cnt - dark]
        out$Calinfo$straylight.corrected <- TRUE
    }
    structure(
        out,
        class = 'caldat'
        , straylight.corrected = correct.straylight
        , linearity.corrected = correct.linearity
        )
}

convert_calref <- function(obj, ref.spec = NULL, dark.spec = NULL) {
    if (is.character(obj)) obj <- qread(obj)
    if (!inherits(obj, 'calref')) stop('conversion only from "calref" object to spec set')
    if (is.null(dark.spec)) {
        dark.spec <- .conv_calref(obj, 'nh3', dark = TRUE)
    }
    list(
        DOAS.model = obj[[1]]$cal_spec$DOASinfo$DOASmodel
        , dat.ref = ref.spec
        , dat.ref.dark = dark.spec
        , dat.dark = dark.spec
        , dat.NH3 = .conv_calref(obj, 'nh3')
        , dat.N2.NH3 = .conv_calref(obj, 'nh3', TRUE)
        , dat.SO2 = .conv_calref(obj, 'so2')
        , dat.N2.SO2 = .conv_calref(obj, 'so2', TRUE)
        , dat.NO = .conv_calref(obj, 'no')
        , dat.N2.NO = .conv_calref(obj, 'no', TRUE)
        , dat.N2.dark = dark.spec
    )
}
.conv_calref <- function(obj, gas, n2 = FALSE, dark = FALSE) {
    # check for calref in main function
    what <- if (n2) 'ref_spec' else 'cal_spec'
    gas <- tolower(gas)
    if (dark) {
        list(
            dat.spec = obj[[gas]][[what]]$Calinfo$dark.spec
            , ambient = list(
                pathLength = NA,
                NH3ambient_ug = NA,
                SO2ambient_ug = NA,
                NOambient_ug = NA 
            )
            , cuvette = list(
                cuvetteLength = NA,
                cuvetteConc_mg = NA,
                cuvetteGas = NA
            )
            , timerange = NA
            , info.spec = 'no info available'
        )
    } else {
        list(
            dat.spec = obj[[gas]][[what]]$Calinfo$raw.spec
            , ambient = list(
                pathLength = NA,
                NH3ambient_ug = NA,
                SO2ambient_ug = NA,
                NOambient_ug = NA 
            )
            , cuvette = list(
                cuvetteLength = obj[[gas]][[what]]$Calinfo$cuvette.path,
                cuvetteConc_mg = obj[[gas]][[what]]$Calinfo$cuvette.conc,
                cuvetteGas = obj[[gas]][[what]]$Calinfo$cuvette.gas
            )
            , timerange = obj[[gas]][[what]]$DOASinfo$timerange
            , info.spec = obj[[gas]][[what]]$calref.info
        )
    }
}

# helper function: find best reference period
find_refperiod <- function(data, ref_duration = 20, n = 1, dn = ref_duration - 1,
    limits = c(i_min = 5e4, d_imax = 5e3, d_nh3 = 0.5, d_so2 = 1, d_no = 1),
    wts = list(
        abs = c(nh3 = 10, so2 = 1, no = 1, i_max = 1e-4),
        diff = c(nh3 = 10, so2 = 5, no = 5, i_max = 2e-3),
        mad = c(nh3 = 20, so2 = 10, no = 10, i_max = 4e-3),
        i_max = 6.6e4
    ), cols = c('nh3', 'so2', 'no', 'i_max')) {
    # subset data
    data <- data[, cols]
    # get i_max subset
    ind_imax <- data$i_max > limits[['i_min']]
    # get individual series
    n_nh3 <- cont_within_range(data$nh3, ind_imax, limits[['d_nh3']], ref_duration)
    n_so2 <- cont_within_range(data$so2, ind_imax, limits[['d_so2']], ref_duration)
    n_no <- cont_within_range(data$no, ind_imax, limits[['d_no']], ref_duration)
    n_imax <- cont_within_range(data$i_max, ind_imax, limits[['d_imax']], ref_duration)
    # find joint series
    n_all <- pmin(n_nh3, n_so2, n_no, n_imax)
    # which are above ref_duration?
    ind_all <- which(n_all == ref_duration - 1)
    # get minimum of nh3, so2, no
    mins <- apply(data[ind_imax, ], 2, min, na.rm = TRUE)
    # calculate penalties
    penalties <- calc_penalty(data, wts, ind_all, mins, ref_duration)
    # return best index
    ord <- ind_all[order(penalties)]
    out <- character(n)
    for (i in seq_len(n)) {
        out[i] <- deparse_timerange(x = st(data)[ord[i]], 
            y = et(data)[ord[i] + ref_duration - 1], sep = ' to ')
        ord <- ord[abs(ord - ord[i]) > dn]
        if (length(ord) == 0) break
    }
    out
}
# get limits from ibts object
get_limits <- function(data, expansion = 1, cols = c('i_max', 'nh3', 'so2', 'no')) {
    sapply(cols, \(x) {
        r <- range(data[, x], na.rm = TRUE)
        if (expansion != 1) {
            r <- mean(r) + diff(r) * c(-1, 1) * expansion
        }
        r
    }, simplify = FALSE)
}
# get statistics
get_stats <- function(data, cols = c('i_max', 'nh3', 'so2', 'no')) {
    out <- sapply(cols, \(x) {
        ind <- is.finite(data[[x]])
        dat <- data[ind,][[x]]
        r <- range(dat)
        inner_out <- c(
            n = sum(ind),
            avg = mean(dat),
            min = r[1],
            # median = median(dat),
            max = r[2],
            maxd = r[2] - r[1],
            sd = sd(dat)
        )
        if (x == 'i_max') {
            inner_out[-1] <- round(inner_out[-1], 0)
        } else {
            inner_out[-1] <- round(inner_out[-1], 2)
        }
    }, simplify = FALSE)
    do.call(rbind, out)
}
# check these reference times visually
check_reftimes <- function(data, ref_times = NULL, 
    x_expansion = '10hours', y_expansion = 1.2, ...) {
    if (is.null(ref_times)) {
        ref_times <- find_refperiod(data, ...)
    }
    ylims <- get_limits(data[ref_times], expansion = y_expansion)
    stats <- list()
    for (i in seq_along(ref_times)) {
        ind2 <- ref_times[i]
        ind <- ind2 %+% paste0(c('-', ''), x_expansion)
        x11(width = 14, height = 10)
        par(mfrow = c(4, 1), mar = c(3, 4, 2, 2))
        # imax
        plot(data[ind, "i_max"], col = "#D9D9D9", ylim = ylims[['i_max']])
        lines(data[ind2, "i_max"], col = "darkgrey")
        # nh3
        plot(data[ind, "nh3"], col = "#DCA0A0", ylim = ylims[['nh3']])
        lines(data[ind2, "nh3"], col = "indianred")
        # so2
        plot(data[ind, "so2"], col = "#A2C4E5", ylim = ylims[['so2']])
        lines(data[ind2, "so2"], col = "dodgerblue3")
        # no
        plot(data[ind, "no"], col = "#81BC9B", ylim = ylims[['no']])
        lines(data[ind2, "no"], col = "seagreen")
        # get statistics
        stats[[ref_times[i]]] <- get_stats(data[ind2])
        cat(ref_times[i], '\n')
        print(stats[[ref_times[i]]])
    }
    structure(ref_times, stats = stats)
}

require(Rcpp)
sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <map>
#include <string>
#include <functional>

// even
double med_even(std::vector<double> vIn, int nIn) {
    return (vIn[nIn / 2 - 1] + vIn[nIn / 2]) / 2;
}

// uneven
double med_uneven(std::vector<double> vIn, int nIn) {
    return vIn[(nIn - 1) / 2];
}

// TODO:
// compare both?: 
//   1) sort first and get ranks (efficient for large values of slices/total)
//   2) get slice and sort afterwards (efficient for few slices in many data points)

// [[Rcpp::export]]
std::vector<double> calc_penalty(
        Rcpp::DataFrame mddataIn, Rcpp::List weightsIn,
        std::vector<int> startsIn, NumericVector minimaIn,
        int reftimeIn) {

    // initialize
    int nrow = mddataIn.nrows();
    int nstarts = startsIn.size();
    NumericVector vcol(nrow);
    NumericVector wabs = weightsIn["abs"];
    NumericVector wdiff = weightsIn["diff"];
    NumericVector wmad = weightsIn["mad"];
    double i_max = weightsIn["i_max"];
    CharacterVector colnames = mddataIn.names();
    std::vector<double> out(nstarts, 0.0);
    std::vector<double> sub(reftimeIn);
    std::vector<double> absdiff(reftimeIn);
    std::string mcall, cname;
    std::map<std::string, std::function<double(std::vector<double>, int)>> median =
        {
            {"even", med_even},
            {"uneven", med_uneven}
        };
    double med, subsum;

    // which median function to call
    if (reftimeIn % 2 == 0) {
        mcall =  "even";
    } else {
        mcall =  "uneven";
    }

    // loop over columns
    for (int col = 0; col < mddataIn.size(); col++) {

        // get column
        cname = colnames[col];
        vcol = mddataIn[cname];

        // loop over startsIn
        for (int i = 0; i < nstarts; i++) {

            // get slice (from/to inclusive)
            sub = std::vector<double>(vcol.begin() + startsIn[i] - 1, vcol.begin() + startsIn[i] + reftimeIn - 1);

            // sort vector subset
            std::sort(sub.begin(), sub.end());

            // 2) add penalty for max diff (second list entry)
            out[i] = out[i] + (sub[reftimeIn - 1] - sub[0]) * wdiff[cname];

            // get median
            med = median[mcall](sub, reftimeIn);

            // set sum-up value to zero
            subsum = 0.0;

            // loop over sorted subset
            for (int s = 0; s < reftimeIn; s++) {

                // sum up
                subsum = subsum + sub[s];

                // get absolut difference
                absdiff[s] = std::fabs(sub[s] - med);

            }

            // sort absdiff
            std::sort(absdiff.begin(), absdiff.end());

            // 3) add penalty for mad
            out[i] = out[i] + median[mcall](absdiff, reftimeIn) * wmad[cname] * 1.4826;

            // 1) add penalty for absolut level
            if (cname == "i_max") {
                // i_max
                out[i] = out[i] + (i_max - subsum / reftimeIn) * wabs[cname];
            } else {
                out[i] = out[i] + (subsum / reftimeIn - minimaIn[cname]) * wabs[cname];
            }

        }

    }

    // return penalties
    return out;
}

// [[Rcpp::export]]
// find continuous series within given range
std::vector<int> cont_within_range (
        std::vector<double> xIn, std::vector<bool> bIn, 
        const double dxIn, const int refLengthIn
        ) {

    // vars
    const int n = xIn.size();
    std::vector<int> out(n);
    int count, j;

    // loop over vector
    for (int i = 0; i < n - 1; i++) {

        // check NA
        if (bIn[i]) {

            // set count to zero
            count = 0;

            // scan forward
            j = i + 1;
            while (j < std::min(n, i + refLengthIn) && bIn[j] && std::fabs(xIn[j] - xIn[i]) <= dxIn) {
                j++;
                count++;
            }

            // assign counts
            out[i] = count;

        } else {
            out[i] = 0;
        }
    }

    // add last element
    out[n - 1] = 0;

    // return result
    return out;
}
')


#### print caldat method
print.caldat <- function(x, ...){
    wl <- get_wl(x)
    cat('~~~~\n')
    cat('\t', x$DOASinfo$DOASmodel, '/', x$DOASinfo$Spectrometer$Serial, '- calibration spectrum:', x$Calinfo$spec.name, '\n\n')
    switch(as.character(nrow(x$Calinfo$info))
        , '1' = cat('\t recorded between', x$Calinfo$info[1, val], '\n')
        , '9' = cat('\t recorded between', x$Calinfo$info[14, val], '\n')
        , cat('\t recorded between', x$Calinfo$info[6, val], '\n')
    )
    cat('\t', min(wl), 'to', max(wl), 'nm', sprintf('(%s pixel)\n\n', length(wl)))
    cat('\t cuvette gas:', x$Calinfo$cuvette.gas, '\n')
    cat('\t cuvette concentration (mg/m3):', x$Calinfo$cuvette.conc, '\n')
    cat('\t cuvette length (m):', x$Calinfo$cuvette.path, '\n\n')
    cat('\t dark corrected:', x$Calinfo$dark.corrected, '\n')
    cat('\t straylight corrected:', x$Calinfo$straylight.corrected, '\n')
    cat('\t linearity corrected:', x$Calinfo$linearity.corrected, '\n')
    cat('~~~~\n')
}

### convert ppm to mg
convert2mg <- function(s, tracer) {
    M <- switch(tolower(tracer), nh3 = 17, no = 30, so2 = 64)
    sp <- as.numeric(strsplit(strsplit(s, split = ': ')[[1]][2], split = ',')[[1]])
    sp[1] * sp[3] / 8.314 / sp[2] / 10 * M
}

#### helper function to get counts
get_cnt <- function(x) {
    if (inherits(x, 'avgdat')) {
        as.numeric(x)
    } else if (inherits(x, 'caldat') || inherits(x, 'chen')) {
        x$data[, cnt]
    } else if (inherits(x, 'rawdat')){
        if (NCOL(x) > 1) {
            x$RawData
        } else {
            x$RawData[[1]]
        }
    }
}

get_integration_time <- function(x) {
    h <- get_header(x)
    int <- unlist(strsplit(grep('(integration|exposure)', h, value = TRUE), split = ': '))
    if (length(int) == 3) {
        as.numeric(sub(')', '', int[3]))
    } else {
        as.numeric(int[2])
    }
}
get_header <- function(rawdat) {
    if (inherits(rawdat, 'character')) {
        stop('reading in raw data from file not yet implemented')
    }
    if (inherits(rawdat, 'dc')) {
        return(get_header(attr(rawdat, 'meas')))
    }
    rawdat$calref.info
}

get_refcal_spec <- function(rawdat, timerange, SpecName, tz = 'Etc/GMT-1', lite = FALSE) {
    tr <- parse_timerange(timerange, tz = tz)
    rd <- list(
        type = '',
        tracer = '',
        rawdat = filter_time(rawdat, tr[1], tr[2])
    )
    getSpec(rd, rd$rawdat$DOASinfo$DOASmodel, SpecName = SpecName, lite = lite)
}



#### doascurve
calc_dc <- function(meas, ref, ftype = NULL, fstrength = NULL, fwin = NULL,
    fitwin = NULL, shift = NULL, correct.straylight = TRUE, correct.linearity = TRUE,
    lin_before_dark = FALSE, do_lowpass_filtering = FALSE, lp.type = 'Rect', lp.strength = 5) {
    if (is.character(meas)) {
        # file path or chen
        if (tolower(meas) == 'chen') {
            # call cheng2dc
            cat('Fix arguments!!!\n')
            browser()
            return(cheng2dc(meas, ref, ftype, fstrength, fwin, fitwin, shift))
        } else {
            meas <- read_cal(meas, correct.straylight = correct.straylight, correct.linearity = correct.linearity, lin_before_dark = lin_before_dark)
        }
    } else if (inherits(meas, 'single_spec') || inherits(meas, 'avgdat')) {
        meas <- attr(meas, 'RawData')
    }
    if (is.character(ref)) ref <- read_cal(ref, correct.straylight = correct.straylight, correct.linearity = correct.linearity, lin_before_dark = lin_before_dark)
    # get counts
    m <- get_cnt(meas)
    r <- get_cnt(ref)
    # get windows
    win <- getWindows(meas$DOASinfo, filter.type = ftype, 
        filter.strength = fstrength, filter.window = fwin, fit.window = fitwin,
        tau.shift = shift)
    # lowpass filtering
    if (do_lowpass_filtering) {
        m <- lowpass.filter(m, lp.type, lp.strength)
        r <- lowpass.filter(r, lp.type, lp.strength)
    }
    # calc diffspec
    ds <- suppressWarnings(log(m / r))
    ds[!is.finite(ds)] <- NA_real_
    # calc dc
    dc <- highpass.filter(ds, win)
    structure(
        list(
            wl = get_wl(meas)[win$pixel_filter],
            cnt = as.numeric(dc)
            ),
        class = 'dc',
        win = win,
        meas = meas,
        ref = ref,
        ds = ds
        )
}

# get local minima values for dc
local_minima <- function(dc, zero_value = -0.2e-20, wl_range = c(203, 219), show = FALSE) {
    # convert to sigma
    dc <- dc2sigma(dc, mgm3 = 193.4095, molar_mass = 17, copy = TRUE)
    if (show) dc_orig <- dc
    # select wl range
    dc$cnt[dc$wl <= wl_range[1] | dc$wl >= wl_range[2]] <- NA
    # dc pixels + wavelength
    dc_px <- attr(dc, 'win')[['pixel_filter']]
    dc_wl <- get_wl(dc)[dc_px]
    # counts
    cnt <- dc$cnt
    # get exact zeros:
    i_exact <- which(abs(diff(sign(cnt))) == 1L)
    i_lwr <- which(abs(diff(sign(cnt))) > 1L)
    i_upr <- i_lwr + 1
    abs_cnt <- abs(cnt)
    ind_zero <- c(i_exact, (i_lwr * abs_cnt[i_lwr] + i_upr * abs_cnt[i_upr]) / (abs_cnt[i_lwr] + abs_cnt[i_upr]))
    # NOTE: -> round: nearer, floor: lower, ceiling: upper indices...
    # local minima
    i_neg <- is.finite(cnt) & cnt < zero_value
    # build groups and find their minima or local minima
    grps <- data.table(
        cnt = cnt, wl = dc$wl, pix = seq.int(cnt),
        grp = cumsum(c(0,diff(as.integer(i_neg))) > 0) * as.integer(i_neg)
    )
    grps[, minima := min(cnt) == cnt, by = grp]
    # find exact minima
    min_exact <- grps[grp != 0, {
        m <- which(minima)
        if (m == 1 || m == .N) {
            px <- as.numeric(pix[m])
        } else if (cnt[m - 1] > cnt[m + 1]) {
            dl <- cnt[m] - cnt[m - 1]
            dr <- cnt[m + 1] - cnt[m + 2]
            px <- (dl * pix[m] + dr * pix[m + 1]) / (dl + dr)
        } else {
            dl <- cnt[m - 1] - cnt[m - 2]
            dr <- cnt[m] - cnt[m + 1]
            px <- (dl * pix[m - 1] + dr * pix[m]) / (dl + dr)
        }
        .(
            px = px,
            wl = dc_wl[floor(px)] + diff(dc_wl[floor(px) + (0:1)]) * (px %% 1)
            )
    }, by = grp]
    # return values
    out <- list(
        pixel_zero = sort(ind_zero) + dc_px[1] - 1,
        pixel_exact = min_exact[, px + dc_px[1] - 1],
        wl_exact = min_exact[, wl],
        pixel_minima = grps[, which(minima) + dc_px[1] - 1],
        sigmas_minima = grps[, cnt[which(minima)]]
        )
    # plot?
    if (show) {
        dc_neg <- dc_orig
        dc_neg$cnt[!i_neg] <- NA
        plot(dc_orig, col = 'lightgrey')
        abline(v = wl_range, lwd = 2)
        abline(h = zero_value, lwd = 2)
        lines(dc_neg, col = 'black', lwd = 2)
        points(dc_wl[grps[, which(minima)]], out$sigmas_minima, col = 'indianred', pch = 20, cex = 1.5)
    }
    out
}

#### helper function to calculate wavelength (correct incl blinds)
calc_wl <- function(x, pixel = NULL) {
    cal.coef <- x$DOASinfo$Spectrometer$'Calibration Coefficients'
    if (is.null(pixel)) pixel <- get_pixel(x)
    cal.coef[1] + cal.coef[2]*pixel + cal.coef[3]*pixel^2 + cal.coef[4]*pixel^3
}

#### helper function to get pixels
get_pixel <- function(x) {
    x$DOASinfo$Spectrometer$pixel
}

# new version of old_cheng2dc (argh!!!)
cheng2dc <- function(cheng, ref, shift = 0, filter = TRUE, fstrength = 5, ftype = 'Rect', ...) {
    # copy cheng data.table
    dta <- copy(cheng$data)
    # create synthetic cal spec
    # -> log(ref) - 193 * cheng2mg_m3 (T = exp(-sum(t)), t = sigma_i * c_i)
    # get/add 'straylight'
    if (inherits(ref, 'avgdat')) {
        old <- ref
        ref <- attr(old, 'RawData')
        # add Calinfo etc.
        ref <- structure(list(
            data = data.table(wl = ref$DOASinfo$Spectrometer$wavelength, cnt = as.numeric(old)),
            Calinfo = list(
                info = data.table(val = 'todo'),
                spec.name = 'dat.NH3.N2',
                cuvette.gas = 'todo',
                cuvette.conc = '-',
                cuvette.path = 'todo'
                ),
            DOASinfo = ref$DOASinfo
        ), class = 'caldat')
    } else {
        ref$data <- copy(ref$data)
    }
    win <- getWindows(ref$DOASinfo$DOASmodel)
    stray <- ref$data[win$pixel_straylight, mean(cnt)]
    if (stray > 1500) {
        # stray <- 2400
        ref$data[, cnt := cnt - stray]
    }
    # cheng (sigma in cm2 / molecule) to m2 / ug
    dta[, cal := sigma2dc(cnt)]
    # add cheng to ref spec
    pseudo_meas <- ref
    pseudo_meas[['Calinfo']][['info']] <- data.table(val = 'Cheng 2006')
    pseudo_meas[['Calinfo']][['cuvette.gas']] <- 'NH3'
    pseudo_meas[['Calinfo']][['cuvette.conc']] <- 193.4095
    pseudo_meas[['Calinfo']][['dark.corrected']] <- TRUE
    # map cal to ref wl
    cal_ref <- cnt2wl(dta[, wl], ref$data[, wl], dta[, cal], shift = shift)
    pseudo_meas$data <- copy(ref$data)
    pseudo_cal <- ref$data[, exp(log(cnt) - cal_ref)]
    pseudo_meas$data[, cnt := pseudo_cal]
    # add filter
    if (filter) {
        pseudo_meas$data[, cnt := lowpass.filter(cnt, ftype, fstrength, ...)]
        ref$data[, cnt := lowpass.filter(cnt, ftype, fstrength, ...)]
    }
    # calc dc
    calc_dc(pseudo_meas, ref)
}
# cheng (sigma in cm2 / molecule) to m2 / ug -> dann mit cuvetten pfad integr. conc mult
# dc werden jeweils durch cuvetten pfad conc in ug / m3 * m geteilt!
sigma2dc <- function(sigma) {
    # cm2 -> m2
    sigma * 1e-4 *
    # 1 / molecule -> 1 / mol
    6.02214076e23 *
    # 1 / mol -> mult mit cuvetten conc in mol / m3 * m
    # cuvetten conc (pfad int.: mol / m3 * m)
    193.4095 / 17e3 * 0.075
}
# map cnt to different wl
cnt2wl <- function(wl_from, wl_to, cnt, shift_nm = 0) {
    dfrom <- diff(wl_from) / 2
    from1 <- wl_from - c(dfrom[1], dfrom) + shift_nm
    from2 <- wl_from + c(dfrom, dfrom[length(dfrom)]) + shift_nm
    dto <- diff(wl_to) / 2
    to1 <- wl_to - c(dto[1], dto)
    to2 <- wl_to + c(dto, dto[length(dto)])
    bins <- cutIntervals(from1, from2, to1, to2)
    unlist(
        lapply(bins, function(b) {
            if (is.null(b)) {
                NA_real_
            } else {
                sum(cnt[b[, 1]] * b[, 2]) / sum(b[, 2])
            }
        })
    )
}

# get Cheng factor
cheng_factor <- function(dc, shift_cheng = 0, show = FALSE, mgm3 = 193.4095) {
    # S5 cheng dc
    cheng <- suppressWarnings(cheng2dc(get_Cheng(), attr(dc, 'ref'), shift = shift_cheng))
    # get sigmas
    dc2sigma(cheng)
    dc <- dc2sigma(dc, copy = TRUE, mgm3 = mgm3, molar_mass = 17)
    # fit
    mod <- lm(dc$cnt ~ cheng$cnt)
    # show?
    if (show) {
        par(mfrow = c(2, 1))
        plot(dc, lwd = 2)
        lines(cheng, col = 'indianred')
        legend('bottomright', bty = 'n', legend = c('meas', 'Cheng2006'), 
            lwd = c(2, 1), col = c('black', 'indianred'))
        plot(cheng$cnt, dc$cnt, xlab = 'Cheng2006', ylab = 'meas')
        legend('bottomright', bty = 'n', legend = sprintf('span = %1.2f', 
                coef(mod)[2]))
        abline(0, 1)
        abline(mod, col = 'lightblue')
    }
    # return
    list(
        coefs = setNames(coef(mod), c('offset', 'span')),
        se = setNames(summary(mod)$coefficients[, 2], c('offset', 'span')),
        rmse = summary(mod)$sigma
        )
}

find_cheng <- function(dc, show = FALSE, interval = c(-1, 1), 
    return.cheng.dc = FALSE, mgm3 = 193.4095) {
    # S5 cheng dc
    cheng <- suppressWarnings(cheng2dc(get_Cheng(), attr(dc, 'ref')))
    lmc <- local_minima(cheng)$wl_exact
    lmm <- local_minima(dc)$wl_exact
    dl <- length(lmm) - length(lmc)
    if (dl != 0) {
        adl <- abs(dl)
        sdl <- sign(dl)
        if (sdl > 0) {
            ref <- lmc
            too_long <- lmm
        } else {
            ref <- lmm
            too_long <- lmc
        }
        # get differences from leave dl out and find minimum
        l <- length(too_long)
        cn <- combn(l, adl)
        combs <- apply(cn, 2, \(x) {
            d <- ref - too_long[-x]
            max(abs(d - median(d)))
        })
        fix <- too_long[-cn[, which.min(combs)]]
        if (sdl > 0) {
            lmm <- fix
        } else {
            lmc <- fix
        }
    }
    ms <- median(lmm - lmc)
    # optimize
    par <- optimize(function(x) {
        cheng_factor(dc, x, mgm3 = mgm3)$rmse
        }, interval = interval + ms)
    c(
        cheng_factor(dc, par$minimum, show = show, mgm3 = mgm3),
        shift = par$minimum,
        # return shifted cheng spectrum
        cheng = if (return.cheng.dc) list(cheng2dc(get_Cheng(), attr(dc, 'ref'), shift = par$minimum))
        )
}

get_cal_infos <- function(dc, compact = TRUE, show = !compact) {
    # local minima
    loc_min <- local_minima(dc[['nh3']][['dc']])#, show = TRUE)
    # -> Imax? Iavg(fit) Imin(fit) Imax(fit)
    cinfo <- attr(dc[['nh3']][['dc']], 'meas')$Calinfo
    # straylight + dark
    dark <- 0
    if (cinfo$dark.corrected) dark <- dark + mean(cinfo$dark.spec, na.rm = TRUE)
    if (cinfo$straylight.corrected) dark <- dark + cinfo$straylight.value
    # cheng factor
    cheng <- find_cheng(dc[['nh3']][['dc']], show)
    out <- c(
        loc_min,
        list(
            integ_ms = get_integration_time(dc[['nh3']][['dc']]),
            light_max = max(cinfo$raw.spec),
            wl_max =  get_wl(dc[['nh3']][['dc']])[which.max(cinfo$raw.spec)],
            light_peaks = cinfo$raw.spec[loc_min$pixel_minima],
            light_dark = dark
        ),
        cheng
        )
    if (compact) {
        list(
            cheng = unlist(cheng)[c(2, 4, 5, 6)],
            light = unlist(out[c(6, 7, 9, 10)]),
            spec = unlist(out[c(2, 5)])
                )
    } else {
        out
    }
}



#### average cheng to reference doas
cheng2doas <- function(cheng, ref, wvlim = NULL) {
    stopifnot(!is.null(ref$DOASinfo))
    if (!is.null(wvlim)) ref <- cut_wl(ref, wvlim[1], wvlim[2])
    # pixel starten erst bei 5...
    pix <- get_pixel(ref) - 4
    wlo <- calc_wl(ref, c(pix[1] - 1, pix, pix[length(pix)] + 1))
    wd <- diff(wlo)
    brks <- wlo[-1] - wd / 2
    wlc <- get_wl(cheng)
    cnt <- get_cnt(cheng)
    # hier umrechnung auf m2/ug 
    cnt <- convertCS(cnt)
    as.numeric(tapply(cnt, cut(wlc, brks), mean))
}

#### convert cm2 / molecule to m2 / ug
# cm2 to m2: factor 1e-4
# molec-1 to ug-1: factor 1/M x Avg x 1e-6
convertCS <- function(x) {
    x / 17 * 6.022 * 10 ^ (23 - 10)
}


#### convert cheng to dc
old_cheng2dc <- function(dc, cheng = NULL, shift = FALSE) {
    stopifnot(!is.null(attr(dc, 'meas')$DOASinfo))
    if (is.null(cheng)) cheng <- get_Cheng()
    if (shift) {
        # recursive call optim
        # return result
    }
    # average cheng to ref
    cdoas <- cheng2doas(cheng, attr(dc, 'meas'))
    # calc dc
    highpass.filter(-cdoas, attr(dc, 'win'))
}

#### plot method for caldat
plot.caldat <- function(x, tau = 0, type = 'l', xlim = c(190, 230), ...) {
    if (tau != 0) {
        # copy data.table!
        x$data <- copy(x$data)
        # fix tau shift
        x$data[, wl := {
            ind <- match(wl, x[['DOASinfo']][['Spectrometer']][['wavelength']])
            # shift indices
            ind_shifted <- ind - tau
            # < 1 should give NA
            ind_shifted[ind_shifted < 1] <- length(ind) + 1
            # get new wavelengths
            x$wl <- x[['DOASinfo']][['Spectrometer']][['wavelength']][ind_shifted]
        }]
    }
    # fix xlim = NULL
    if (is.null(xlim)) xlim <- x$data[, range(wl, na.rm = TRUE)]
    x$data[wl >= xlim [1] & wl <= xlim[2], plot(wl, cnt, type = type, xlim = xlim, ...)]
    invisible()
}
lines.caldat <- function(x, tau = 0, ...) {
    if (tau != 0) {
        # copy data.table!
        x$data <- copy(x$data)
        # fix tau shift
        x$data[, wl := {
            ind <- match(wl, x[['DOASinfo']][['Spectrometer']][['wavelength']])
            # shift indices
            ind_shifted <- ind - tau
            # < 1 should give NA
            ind_shifted[ind_shifted < 1] <- length(ind) + 1
            # get new wavelengths
            x$wl <- x[['DOASinfo']][['Spectrometer']][['wavelength']][ind_shifted]
        }]
    }
    x$data[, lines(wl, cnt, ...)]
    invisible()
}

#### plot method for dc
plot.dc <- function(x, fctr = 1, tau = 0, per_molecule = FALSE, type = 'l', xlab = 'nm', ylab = 'doascurve', ...) {
    x$cnt <- x$cnt * fctr
    if (tau != 0) {
        # fix tau shift
        ind <- match(x$wl, attr(x, 'meas')[['DOASinfo']][['Spectrometer']][['wavelength']])
        # shift indices
        ind_shifted <- ind - tau
        # < 1 should give NA
        ind_shifted[ind_shifted < 1] <- length(ind) + 1
        # get new wavelengths
        x$wl <- attr(x, 'meas')[['DOASinfo']][['Spectrometer']][['wavelength']][ind_shifted]
    }
    if (per_molecule) {
        # sigma2dc inverse
        dc2sigma(x) 
    }
    plot(x$wl, x$cnt, type = type, xlab = xlab, ylab = ylab, ...)
}
lines.dc <- function(x, fctr = 1, tau = 0, per_molecule = FALSE, ...) {
    x$cnt <- x$cnt * fctr
    if (tau != 0) {
        # fix tau shift
        ind <- match(x$wl, attr(x, 'meas')[['DOASinfo']][['Spectrometer']][['wavelength']])
        # shift indices
        ind_shifted <- ind - tau
        # < 1 should give NA
        ind_shifted[ind_shifted < 1] <- length(ind) + 1
        # get new wavelengths
        x$wl <- attr(x, 'meas')[['DOASinfo']][['Spectrometer']][['wavelength']][ind_shifted]
    }
    if (per_molecule) {
        # sigma2dc inverse
        dc2sigma(x) 
    }
    lines(x$wl, x$cnt, ...)
}
points.dc <- function(x, fctr = 1, tau = 0, per_molecule = FALSE, type = 'p', ...) {
    lines(x, fctr = 1, tau = 0, per_molecule = FALSE, type = type, ...)
}
# sigma2dc inverse
dc2sigma <- function(dc, mgm3 = NULL, molar_mass = NULL, copy = FALSE) {
    dc_name <- deparse(substitute(dc))
    if (is.null(mgm3)) {
        meas <- attr(dc, 'meas')
        mgm3 <- switch(class(meas)[1]
            , rawdat = meas[['Header']][['cuvetteConc_mg']]
            , caldat = attr(dc, 'meas')[['Calinfo']][['cuvette.conc']]
            , {
                cat('dc2sigma not yet implemented for class', class(meas)[1], '!\n')
                browser()
            }
            )
    }
    if (is.null(molar_mass)) {
        meas <- attr(dc, 'meas')
        gas <- switch(class(meas)[1]
            , rawdat = meas[['Header']][['cuvetteGas']]
            , caldat = attr(dc, 'meas')[['Calinfo']][['cuvette.gas']]
            , {
                cat('dc2sigma not yet implemented for class', class(meas)[1], '!\n')
                browser()
            }
            )
        molar_mass <- switch(tolower(gas),
            nh3 = 17,
            so2 = 64,
            no = 30
            )
    }
    dc$cnt <- dc$cnt / 1e-4 / 6.02214076e23 / mgm3 * molar_mass * 1e3 / 0.075
    if (!copy) {
        assign(dc_name, dc, sys.frame(-1))
    }
    invisible(dc)
}

#### plot chen
plot.chen <- function(x, xlim = c(190, 230), type = 'l', ...) {
    x$data[, plot(wl, cnt, type = type, xlim = xlim, ...)]
}


#### get timerange
get_timerange <- function(obj) {
    switch(class(obj)[1]
        , 'avgdat' = 
        , 'single_spec' = attr(obj, 'RawData')[['DOASinfo']][['timerange']]
        , 'dc' = attr(obj, 'meas')[['DOASinfo']][['timerange']]
        , 'caldat' = 
        , 'rawdat' = obj[['DOASinfo']][['timerange']]
        , stop('class not yet implemented')
    )
}


#### further plotting helpers
ann_time <- function(obj, x = NULL, y = NULL, x_adj = 0, y_adj = 0, cex = 0.6, ...) {
    # get timerange
    timerange <- get_timerange(obj)
    if (is.character(x)) {
        text(x, labels = ibts::deparse_timerange(timerange, sep = ' to '), cex = cex, ...)
    } else {
        # get usr coordinates
        usr <- par('usr')
        # default position
        if (is.null(x)) {
            # 3/4 of range
            x <- usr[1] + (usr[2] - usr[1]) * 3 / 4
        }
        if (is.null(y)) {
            # 3/4 of usr range
            if (par('ylog')) {
                y <- 10 ^ (usr[3] + (usr[4] - usr[3]) * 3 / 4)
            } else {
                y <- usr[3] + (usr[4] - usr[3]) * 3 / 4
            }
        }
        # adjustments
        x_add <- (usr[2] - usr[1]) * x_adj
        if (par('ylog')) {
            y_add <- 10 ^ (log(y, 10) + (usr[4] - usr[3]) * y_adj) - y
        } else {
            y_add <- (usr[4] - usr[3]) * y_adj
        }
        # add text
        text(x + x_add, y + y_add, labels = ibts::deparse_timerange(timerange, sep = ' to '), cex = cex, ...)
    }
}

annotate_revolver <- function(obj, side = 'top', lcol = 'lightblue',
    tcol = 'black', add_vlines = FALSE, add_jitter = TRUE, 
    jitter_amount = c(n2 = 0.05, nh3 = 0, no = -0.05, so2 = -0.1)
    ) {
    st_ind <- st(obj)
    et_ind <- et(obj)
    rl <- rle(obj$revolver)
    sts <- cumsum(c(1, rl[[1]]))
    ets <- sts[-1] - 1
    usr <- par('usr')
    if (is.character(side)) {
        side <- c(top = 0.9, bottom = 0.1)[side]
    }
    if (missing(jitter_amount) && length(side) == 1 && side < 0.5) {
        jitter_amount <- -jitter_amount
    }
    if (add_jitter) {
        side <- side + jitter_amount[tolower(rl[[2]])]
    } else {
        side <- rep(side, length(ets))[seq_along(ets)]
    }
    y <- usr[3] + (usr[4] - usr[3]) * side
    for (i in seq_along(ets)) {
        x1 <- st_ind[sts[i]]
        x2 <- et_ind[ets[i]]
        at <- x1 + (x2 - x1) / 2
        lines(c(x1, x2), c(y[i], y[i]), col = lcol)
        text(at, y[i], labels = rl[[2]][i], col = tcol)
        if (add_vlines) {
            abline(v = c(x1, x2), col = lcol)
        }
    }
}

gather_callists <- function(nh3_callist, no_callist, so2_callist) {
    list(
        nh3 = nh3_callist$nh3, 
        no = no_callist$no, 
        so2 = so2_callist$so2,
        n2 = list(
            raw = list(nh3_callist$n2$raw, no_callist$n2$raw, so2_callist$n2$raw),
            sets = c(nh3_callist$n2$sets, no_callist$n2$sets, so2_callist$n2$sets),
            avgs = c(nh3_callist$n2$avgs, no_callist$n2$avgs, so2_callist$n2$avgs),
            i_max = c(nh3_callist$n2$i_max, no_callist$n2$i_max, so2_callist$n2$i_max)
        )
    )
}


#### save calibration to files
save_calspec <- function(x, save_cal = FALSE, save_online = FALSE) {
    md <- sub('^(S[1-6])_.*', '\\1', deparse(substitute(x)))
    folder <- file.path(calsave_dir, 'Calibration-Aula-Feb2021', md)
    dir.create(folder, showWarnings = FALSE, recursive = TRUE)
    rawdat <- readDOASdata(md, x$dir, TRUE, timerange = x$timerange, timezone = x$tz)
    out <- avgSpec(
        rawdat,
        type = x$type,
        tracer = x$tracer,
        cuvetteConc_mg = x$cuvetteConc_mg,
        Dirname = folder,
        saveToFile = save_cal && !save_online
    )
    if (save_cal) {
        if (save_online) {
            folder <- file.path(refDir, 'OnlineSpectra', md)
            dir.create(folder, showWarnings = FALSE, recursive = TRUE)
            fwrite(rbind(data.frame(a = rep('', 7)),
                    'Number of accumulations=1', data.frame(a = out$SpecAvg)), file.path(refDir, 'OnlineSpectra', md, 
                    paste0('calibration-aula-feb2021-', x$tracer, '-', x$type)), 
                col.names = FALSE)
        }
        invisible(NULL)
    } else {
        out
    }
}

#### get windows
get_wins <- function(x) {
    switch(class(x)[1]
        , 'character' = {
            doas_info <- getDOASinfo(x)
            getWindows(doas_info)
        }
        , 'single_spec' = {
            getWindows(attr(x, 'RawData')[['DOASinfo']])
        }
        , 'dc' = {
            getWindows(attr(x, 'meas')[['DOASinfo']])
        }
        , 'rawdat' = {
            getWindows(x[['DOASinfo']])
        }
        , stop('no method definded for current class')
        )
}

#### fit calibration spectra to measured single spectrum
fit_dc <- function(x, calrefspecs = NULL, tau.shift = 0, path.length = 1,
    dcnh3 = NULL, dcno = NULL, dcso2 = NULL, corNH3 = 1.16, robust = TRUE) {
    wins <- get_wins(x)
    if (!is.null(calrefspecs)) {
        # get doascurves from calibration spectra
        # nh3
        if (is.null(dcnh3)) {
            snh3 <- read_cal(calrefspecs, spec = 'dat.NH3')
            snh3_ref <- read_cal(calrefspecs, spec = 'dat.N2.NH3')
            dcnh3 <- calc_dc(snh3, snh3_ref)
        }
        # no
        if (is.null(dcno)) {
            sno <- read_cal(calrefspecs, spec = 'dat.NO')
            sno_ref <- read_cal(calrefspecs, spec = 'dat.N2.NO')
            dcno <- calc_dc(sno, sno_ref)
        }
        # so2
        if (is.null(dcso2)) {
            sso2 <- read_cal(calrefspecs, spec = 'dat.SO2')
            sso2_ref <- read_cal(calrefspecs, spec = 'dat.N2.SO2')
            dcso2 <- calc_dc(sso2, sso2_ref)
        }
    }
    xreg <- cbind(nh3 = dcnh3[['cnt']], no = dcno[['cnt']], so2 = dcso2[['cnt']])[wins[['pixel_dc']], ]
    if (robust) {
        out <- fit.curves.rob(x[['cnt']], wins[['pixel_dc']], xreg, tau.shift, path.length, all_coefs = TRUE)
    } else {
        out <- fit.curves(x[['cnt']], wins[['pixel_dc']], xreg, tau.shift, path.length, all_coefs = TRUE)
    }
    names(out) <- c(colnames(xreg), paste0(colnames(xreg), '_se'), 'tau.best', 'intercept')
    # return result incl calibration
    structure(out,
        nh3_cal = attr(dcnh3, 'meas')[['Calinfo']][['cuvette.conc']] * attr(dcnh3, 'meas')[['Calinfo']][['cuvette.path']] * 1e3 * corNH3,
        no_cal = attr(dcno, 'meas')[['Calinfo']][['cuvette.conc']] * attr(dcno, 'meas')[['Calinfo']][['cuvette.path']] * 1e3,
        so2_cal = attr(dcso2, 'meas')[['Calinfo']][['cuvette.conc']] * attr(dcso2, 'meas')[['Calinfo']][['cuvette.path']] * 1e3,
        tau = tau.shift,
        path = path.length,
        dc_nh3 = dcnh3,
        dc_no = dcno,
        dc_so2 = dcso2,
        dc_meas = x,
        cor_nh3 = corNH3,
        class = 'doas.fit'
    )
}

#### convert fit result to ug/m3
fit2ug <- function(fit, path.length = 1) {
    out <- fit
    ind <- grep('nh3|no|so2', names(fit))
    for (i in ind) {
        cal <- sub('(nh3|no|so2).*', '\\1_cal', names(fit)[i])
        out[[i]] <- out[[i]] * attr(fit, cal) * attr(fit, 'path') / path.length
    }
    out
}

#### get residues
resid_dc <- function(x, calrefspecs = NULL, tau.shift = 0, path.length = 1,
    dcnh3 = NULL, dcno = NULL, dcso2 = NULL, corNH3 = 1.16) {
    if (inherits(x, 'doas.fit')) {
        tau.shift <- attr(x, 'tau')
        path.length <- attr(x, 'path')
        fitted <- x
        x <- attr(x, 'dc_meas')
    } else {
        fitted <- fit_dc(x, calrefspecs, tau.shift, path.length,
            dcnh3, dcno, dcso2, corNH3)
    }
    structure(
        list(
            wl = x$wl,
            cnt = shift_tau(x[['cnt']], tau.shift) - 
                fitted[['nh3']] * attr(fitted, 'dc_nh3')[['cnt']] -
                fitted[['no']] * attr(fitted, 'dc_no')[['cnt']] -
                fitted[['so2']] * attr(fitted, 'dc_so2')[['cnt']] - 
                fitted[['intercept']]
            ),
        fit = fitted,
        class = 'dc_resid')
}
#### helper to shift by tau
shift_tau <- function(x, tau) {
    ind <- seq_along(x) + tau
    ind[ind < 1] <- NA
    x[ind]
}

#### dc operators
Ops.dc <- function(e1, e2) {
    # get classes
    cl1 <- class(e1)
    cl2 <- class(e2)
    # dc
    if ('dc' %in% cl1 && 'dc' %in% cl2) {
        gfun <- get(.Generic)
        e <- e1
        e[['cnt']] <- gfun(e1[['cnt']], e2[['cnt']])
        return(e)
    }
    # dc_resid
    if (any(which_resid <- c(cl1, cl2) %in% 'dc_resid')) {
        if (which_resid[1]) {
            e <- e2
            # shift resid
            t1 <- -attr(attr(e1, 'fit'), 'tau')
            t2 <- 0
        } else {
            e <- e1
            # shift resid
            t1 <- 0
            t2 <- -attr(attr(e2, 'fit'), 'tau')
        }
        fun <- get(.Generic)
        e[['cnt']] <- fun(shift_tau(e1[['cnt']], t1), shift_tau(e2[['cnt']], t2))
        return(e)
    }
    # fix_pattern
    if (any(which_fp <- c(cl1, cl2) %in% 'fix_pattern')) {
        if (which_fp[1]) {
            stop('add/substract fix patterns from doascurves, not vice versa')
        } else {
            # check dimension of e2
            if (NCOL(e2) > 1) stop('Only one fix pattern can be substracted')
            e <- e1
            # get fp tau shift
            t1 <- 0
            t2 <- -attr(e2, 'DOASwin')[['tau.shift']]
            # get fp indices
            ind_fp <- attr(e2, 'DOASwin')[['pixel_dc']]
            f <- e1
            f[['cnt']][] <- 0
            f[['cnt']][ind_fp] <- as.numeric(e2)
            fun <- get(.Generic)
            e[['cnt']] <- fun(shift_tau(e1[['cnt']], t1), shift_tau(f[['cnt']], t2))
            return(e)
        }
    }
    # numeric
    if (any(which_num <- c(cl1, cl2) %in% c('numeric', 'integer'))) {
        gfun <- get(.Generic)
        if (which_num[1]) {
            e <- e2
            e[['cnt']] <- gfun(e1, e2[['cnt']])
        } else {
            e <- e1
            e[['cnt']] <- gfun(e1[['cnt']], e2)
        }
        return(e)
    }
    stop('Ops.dc method not yet implemented')
}


#### plot method for dc_resid
plot.dc_resid <- function(x, type = 'l', xlab = 'nm', ylab = 'doascurve', ...) {
    plot(x$wl, x$cnt, type = type, xlab = xlab, ylab = ylab, ...)
}
lines.dc_resid <- function(x, ...) {
    lines(x$wl, x$cnt, ...)
}

print.doas.fit <- function(x, ...) {
    for (nms in names(x)) {
        cat(nms, ':', x[[nms]], '\n')
    }
    invisible(NULL)
}

get_fixpattern <- function(x, granularity = '1days', ref.time = colnames(x)[1], fun = NULL, ...) {
    # parse ref time
    t_zero <- parse_date_time3(ref.time)
    # get times (st)
    nms <- colnames(x)
    abs_times <- parse_date_time3(nms)
    rel_times <- abs_times - t_zero
    # binning intervals
    gr_time <- parse_time_diff(granularity)
    st_bins <- floor(rel_times / gr_time)
    # get times
    d_bins <- which(diff(st_bins) > 0)
    u_times <- paste0(abs_times[c(1, d_bins + 1)], ' to ', abs_times[c(d_bins, length(abs_times))])
    # loop over bins
    u_bins <- unique(st_bins)
    out <- matrix(NA_real_, nrow = nrow(x), ncol = length(u_bins), dimnames = list(NULL, u_times))
    for (b in u_bins) {
        if (is.null(fun)) {
            out[, u_bins == b] <- rowMeans(x[, st_bins == b, drop = FALSE])
        } else {
            out[, u_bins == b] <- apply(x[, st_bins == b, drop = FALSE], 1, fun, ...)
        }
    }
    structure(out
        , DOASinfo = attr(x, 'DOASinfo')
        , DOASwin = attr(x, 'DOASwin')
        , class = 'fix_pattern'
        )
}

plot.fix_pattern <- function(x, max.plots = 15, scale_y = TRUE, ylim = NULL, ylab = '', tau = NULL, ...) {
    nc <- NCOL(x)
    # do we have to many plots to display at once?
    if (nc > max.plots) {
        devAskNewPage(TRUE)
        cols <- rows <- 2
    } else {
        # get number of cols & rows
        cols <- floor(sqrt(nc))
        rows <- ceiling(nc / cols)
    }
    # start plotting
    if (nc > 1) {
        opar <- par()
        par(mfrow = c(rows, cols))
    }
    if (scale_y && is.null(ylim)) ylim <- range(x, na.rm = TRUE)
    if (is.null(tau)) tau <- 0
    ind <- attr(x, 'DOASwin')[['pixel_fit']] + tau
    wl <- attr(x, 'DOASinfo')[['Spectrometer']][['wavelength']][ind]
    cnms <- colnames(x)
    class(x) <- 'matrix'
    for (i in seq.int(nc)) {
        plot(wl, x[, i], type = 'l', ylim = ylim, ylab = ylab, main = cnms[i], ...)
    }
    if (nc > 1) {
        suppressWarnings(par(opar))
    }
    invisible(NULL)
}
lines.fix_pattern <- function(x, tau = NULL, ...) {
    nc <- NCOL(x)
    # start plotting
    if (is.null(tau)) tau <- 0
    ind <- attr(x, 'DOASwin')[['pixel_fit']] + tau
    wl <- attr(x, 'DOASinfo')[['Spectrometer']][['wavelength']][ind]
    cnms <- colnames(x)
    class(x) <- 'matrix'
    for (i in seq.int(nc)) {
        lines(wl, x[, i], ...)
    }
    invisible(NULL)
}

print.dc <- function(x, ...) {
    meas <- attr(x, 'meas')
    ref <- attr(attr(x, 'ref'), 'RawData')
    if (is.null(meas$Header[['st']])) {
        meas_txt <- meas$Calinfo$info[, val]
    } else {
        meas_txt <- paste(format(meas$Header[['st']][1]), 'and', format(meas$Header[['et']][length(meas$Header[['et']])]))
    }
    if (is.null(ref$Header[['st']])) {
        ref_txt <- attr(x, 'ref')$Calinfo$info[, val]
    } else {
        ref_txt <- paste(format(ref$Header[['st']][1]), 'and', format(ref$Header[['et']][length(ref$Header[['et']])]))
    }
    cat('***\ndoascurve:\n')
    cat('   measurement spectrum recorded between', meas_txt, '\n')
    cat('   reference spectrum recorded between', ref_txt, '\n')
    cat('***\n')
}

print.fix_pattern <- function(x, ...) {
    cat('***\nfix patterns calculated for:\n')
    cnms <- colnames(x)
    for (i in seq_along(cnms)) {
        cat('   ', sprintf('%2i', i), ': ', cnms[i], '\n')
    }
    cat('***\n')
}

#' Filter miniDOAS Raw Data By Time
#'
#' function to select a subset (filter) from a miniDOAS raw data list
#' by providing either a time range or individual times
#'
#' @param rawdat     A \code{list} with entries \code{'RawData'}, \code{'Header'}, 
#'                   \code{'DOASinfo'}, obtained from a call to \code{readDOASdata} 
#'                   or \code{read_data}.
#' @param from       A vector of time indices used to subset
#' @param to         An optional time vector providing the end of a time range. 
#'                   Argument \code{from} has to be of length 1 to provide
#'                   the start of the time range. Defaults to \code{NULL}.
#' @param including  logical. Only used if time range is provided. Should
#'                   data at the range edges be included?
#' @return A subset of the initial \code{rawdat} object
filter_time <- function(rawdat, from, to = NULL, including = TRUE) {
    rd_st <- rawdat[['Header']][['st']]
    rd_et <- rawdat[['Header']][['et']]
    rd_tz <- tzone(rd_st)
    # check to
    if (!is.null(to)) {
        # parse from/to
        if (is.character(from)) {
            from <- parse_date_time3(from, tz = rd_tz)
        }
        if (is.character(to)) {
            to <- parse_date_time3(to, tz = rd_tz)
        }
        ind <- integer(0)
        for (i in seq_along(to)) {
            # from / to
            if (including) {
               ind <- c(ind, which(rd_et > from[i] & rd_st < to[i]))
            } else {
               ind <- c(ind, which(rd_st >= from[i] & rd_et <= to[i]))
            }
        }
        # sort index
        ind <- sort(ind)
    } else {
        # from only
        # check if timerange
        if (is.character(from)) {
            seps <- getOption("time.separators")                                                                                                                
            split_char <- seps[sapply(seps, grepl, from)]
            if (length(split_char) != 0) {
                fromto <- lapply(from, parse_timerange, tz = rd_tz)
                from_new <- sapply(fromto, '[[', 1)
                to_new <- sapply(fromto, '[[', 2)
                if (is.na(from_new)) from_new <- rd_st[1]
                if (is.na(to_new)) to_new <- rd_et[length(rd_et)]
                return(
                    filter_time(rawdat, from_new, to_new, including = including)
                )
            }
        }
        # parse from
        if (is.character(from)) {
            from <- parse_date_time3(from, tz = rd_tz)
        }
        # find intervals
        ind <- ibts::findI_st(as.numeric(from), as.numeric(rd_st), as.numeric(rd_et))
    }
    # return rawdat class
    structure(
        list(
            RawData = rawdat[['RawData']][ind],
            Header = rawdat[['Header']][ind, ],
            DOASinfo = rawdat[['DOASinfo']]
            ),
        class = 'rawdat'
        )
}

# filter by simple index
filter_index <- function(rawdat, index) {
    # return rawdat class
    structure(
        list(
            RawData = rawdat[['RawData']][index],
            Header = rawdat[['Header']][index, ],
            DOASinfo = rawdat[['DOASinfo']]
            ),
        class = 'rawdat'
        )
}

# filter raw data by revolver position
filter_closed <- function(rawdat) {
    # get indices
    ind <- which(rawdat[['Header']][['ShuPos']] %in% c('close', 'closed'))
    # throw error if position doesn't exist
    if (length(ind) == 0) {
        stop(
            'Shutter positions "close" or "closed" do not exist!\n',
            '  Existing positions are:\n    ', 
            paste(unique(rawdat[['Header']][['ShuPos']]), collapse = ', '),
            '\n Fix in function "filter_closed" if necessary'
            )
    }
    # check any unlikely errors
    stopifnot(
        length(unique(rawdat[['Header']][ind, 'Spectrometer'])) == 1 || 
        length(unique(rawdat[['Header']][ind, 'DOASmodel'])) == 1
    )
    # return rawdat class
    structure(
        list(
            RawData = rawdat[['RawData']][ind],
            Header = rawdat[['Header']][ind, ],
            DOASinfo = rawdat[['DOASinfo']]
            ),
        class = 'rawdat'
        )
}

# filter raw data by revolver position
filter_position <- function(rawdat, position) {
    # get indices
    ind <- which(tolower(rawdat[['Header']][['RevPos']]) %in% tolower(position))
    # throw error if position doesn't exist
    if (length(ind) == 0) {
        stop(
            'Revovler position ', position, ' does not exist!\n',
            '  Existing positions are:\n    ', 
            paste(unique(rawdat[['Header']][['RevPos']]), collapse = ', ')
            )
    }
    # check any unlikely errors
    stopifnot(
        length(unique(rawdat[['Header']][ind, 'Spectrometer'])) == 1 || 
        length(unique(rawdat[['Header']][ind, 'DOASmodel'])) == 1
    )
    # return rawdat class
    structure(
        list(
            RawData = rawdat[['RawData']][ind],
            Header = rawdat[['Header']][ind, ],
            DOASinfo = rawdat[['DOASinfo']]
            ),
        class = 'rawdat'
        )
}
# evaluate rawdata
eval_raw <- function(rawdat, calspecs, path.length, tau.shift = 0, lite = TRUE, ...) {
    evalOffline(
        path.length = path.length,
        tau.shift = tau.shift,
        RawData = rawdat,
        CalRefSpecs = calspecs,
        lite = lite,
        ...
        )
}

# function to clean calibration blocks
clean_cal <- function(x_dt, max_time_gap = 5 * 60, show = TRUE,
    dev_med = function(x) 3 * mad(x, na.rm = TRUE), 
    exclude_times = NULL, subset = NULL) {
    # check x_dt
    if (!is.data.table(x_dt)) {
        stop('Please provide a data.table containing minidoas evaluation results')
    } else {
        x_dt <- copy(x_dt)
    }
    # check dev_med
    if (length(dev_med) != 4) {
        dev_med <- setNames(lapply(1:4, function(i)
            dev_med), c('nh3', 'no', 'so2', 'i_max'))
    }
    # exclude times
    if (!is.null(exclude_times)) {
        cat('Fix argument "exclude_times"!\n')
        browser()
    }
    if (!missing(subset)) {
        subs <- substitute(subset)
        x_dt[, excl_subset := !eval(subs)]
    } else {
        x_dt[, excl_subset := FALSE]
    }
    # find blocks of cal (mt diff > 5? min)
    x_dt[, cal_block := {
        xs <- as.numeric(st)
        xe <- as.numeric(et)
        paste0(revolver, '_', cumsum(c(0, as.integer((xe[-1] - xs[-.N]) > max_time_gap))))
    }]
    # find median values by block & highlight deviation > med +/- ?
    nms <- names(dev_med)
    x_dt[, paste0('dev_', nms) := {
        lapply(nms, function(i) {
            iv <- get(i)
            abs(iv - median(iv, na.rm = TRUE)) > dev_med[[i]](iv)
        })
    }, by = cal_block]
    # time excludes
    # plot
    if (show) {
        x <- as.ibts(x_dt)
        x_dev <- x[, c(nms, paste0('dev_', nms))]
        x_all <- x[, c(nms, paste0('dev_', nms))]
        x_dev[!x[['dev_nh3']], 'nh3'] <- NA_real_
        x_dev[!x[['dev_no']], 'no'] <- NA_real_
        x_dev[!x[['dev_so2']], 'so2'] <- NA_real_
        x_dev[!x[['dev_i_max']], 'i_max'] <- NA_real_
        x_all[x[['dev_nh3']] | x[['dev_no']] | x[['dev_so2']] | x[['dev_i_max']], ] <- NA_real_
        # add_blocks <- function(){
        #     usr <- par('usr')
        #     md <- (usr[4] - usr[3]) * 0.05
        #     m <- usr[4] + md
        #     clip(usr[1], usr[2], usr[4], m + md * 2)
        #     x_dt[, {
        #         lines(c(st[1], et[.N]), rep(m, 2))
        #         lines(rep(st[1], 2), m + c(-1 , 1) * md / 3)
        #         lines(rep(et[.N], 2), m + c(-1 , 1) * md / 3)
        #         text(mean(c(st[1], et[.N])), y = m + md, labels = .BY[[1]])
        #     }, by = cal_block]
        # }
        par(mfrow = c(4, 1))
        plot(x[, 'nh3'], gap.size.max = '1hours', col = 'lightgrey', ylim = range(x_all[, 'nh3'], na.rm = TRUE))
        lines(x_dev[, 'nh3'], gap.size.max = '1hours', col = 'red')
        lines(x_all[, 'nh3'], gap.size.max = '1hours')
        # add_blocks()
        plot(x[, 'no'], gap.size.max = '1hours', col = 'lightgrey', ylim = range(x_all[, 'no'], na.rm = TRUE))
        lines(x_dev[, 'no'], gap.size.max = '1hours', col = 'red')
        lines(x_all[, 'no'], gap.size.max = '1hours')
        plot(x[, 'so2'], gap.size.max = '1hours', col = 'lightgrey', ylim = range(x_all[, 'so2'], na.rm = TRUE))
        lines(x_dev[, 'so2'], gap.size.max = '1hours', col = 'red')
        lines(x_all[, 'so2'], gap.size.max = '1hours')
        plot(x[, 'i_max'], gap.size.max = '1hours', col = 'lightgrey', ylim = range(x_all[, 'i_max'], na.rm = TRUE))
        lines(x_dev[, 'i_max'], gap.size.max = '1hours', col = 'red')
        lines(x_all[, 'i_max'], gap.size.max = '1hours')
        title(x[['revolver']][1], outer = TRUE, line = -1.5, cex.main = 1.5)
    }
    # add overall index
    if ('excl_subset' %in% names(x_dt)) {
        x_dt[, ok := !(dev_nh3 | dev_no | dev_so2 | dev_i_max | excl_subset)]
    } else {
        x_dt[, ok := !(dev_nh3 | dev_no | dev_so2 | dev_i_max)]
    }
    # add column block
    x_dt[, block := NA_integer_][(ok), block := frank(as.integer(sub('.*_', '', cal_block)), ties.method = 'dense')]
    # return ibts
    as.ibts(x_dt[(ok)])
}

# plot all important cal values
plot_cal <- function(x, rev = if (length(unique(x[['revolver']])) == 1) x[['revolver']][1] else 'N2',
    fac = c(nh3 = 1, so2 = 1, no = 1), full.range = FALSE, ...) {
    # scale before plotting
    for (f in names(fac)) {
        x[, f] <- x[, f] * fac[f]
    }
    if (full.range) rev <- unique(x[['revolver']])
    par(mfrow = c(4, 1))
    plot(x[, 'nh3'], gap.size.max = '1hours', ylim = range(x[x$revolver %in% rev, 'nh3']), ...)
    plot(x[, 'no'], gap.size.max = '1hours', ylim = range(x[x$revolver %in% rev, 'no']), ...)
    plot(x[, 'so2'], gap.size.max = '1hours', ylim = range(x[x$revolver %in% rev, 'so2']), ...)
    plot(x[, 'i_max'], gap.size.max = '1hours', ylim = range(x[x$revolver %in% rev, 'i_max']), ...)
}


if (FALSE) {
    FileAulaNH3 <- '~/repos/3_Scripts/4_MiniDOASAuswertung/ReferenceSpectras/S5/miniDOAS_S5_hafl_aula_spectra_2104081407/miniDOAS_S5_NH3_cal_spec_210222592007-210222172106_202104081407.txt'
    FileAulaN2 <- '~/repos/3_Scripts/4_MiniDOASAuswertung/ReferenceSpectras/S5/miniDOAS_S5_hafl_aula_spectra_2104081407/miniDOAS_S5_NH3_ref_spec_210222052205-210222212205_202104081407.txt'
    FileOldNH3 <- '~/repos/3_Scripts/4_MiniDOASAuswertung/ReferenceSpectras/miniDOAS-ref_cal_spec-Eval201505111640_S5/miniDOAS_NH3calSpec_201410281529-201410281533-Eval201505111640.txt'
    FileOldN2 <- '~/repos/3_Scripts/4_MiniDOASAuswertung/ReferenceSpectras/miniDOAS-ref_cal_spec-Eval201505111640_S5/miniDOAS_N2calSpec_201410281537-201410281541-Eval201505111640.txt'
    lin_first <- FALSE
    c_dark <- TRUE
    c_lin <- TRUE
    # c_dark <- FALSE
    # c_lin <- FALSE
    # lin_first <- TRUE
    dcAula <- calc_dc(FileAulaNH3, FileAulaN2, correct.straylight = c_dark, correct.linearity = c_lin, lin_before_dark = lin_first)
    dcOld <- calc_dc(FileOldNH3, FileOldN2, correct.straylight = c_dark, correct.linearity = c_lin, lin_before_dark = lin_first)
    dsAula <- attr(dcAula, 'ds')[attr(dcAula, 'win')$pixel_filter]
    dsOld <- attr(dcOld, 'ds')[attr(dcOld, 'win')$pixel_filter]
    chAula <- old_cheng2dc(dcAula)
    chOld <- old_cheng2dc(dcOld)

    type <- 'l'
    plot(dcAula$wl, dcAula$cnt, type = type)
    # lines(dcAula$wl, chAula * (193e3 * 0.075), col = 'red', type = type)
    lines(dcOld$wl, dcOld$cnt, col = 'orange', type = type)

    # plot(dcAula$wl, dsAula, type = 'l', ylim = c(-0.8, 0))
    plot(calc_wl(attr(dcAula, 'meas'), attr(dcAula, 'win')$pixel_filter - 10), dsAula, type = 'l', ylim = c(-0.8, 0))
    lines(get_wl(nh3), cheng * 193e3 * 0.075, col = 'orange')

    dc1 <- calc_dc(FileAulaNH3, FileAulaN2, correct.straylight = c_dark, correct.linearity = c_lin, lin_before_dark = FALSE)
    dc2 <- calc_dc(FileAulaNH3, FileAulaN2, correct.straylight = c_dark, correct.linearity = c_lin, lin_before_dark = TRUE)
    type <- 'l'
    plot(dc1$wl, dc1$cnt, type = type)
    lines(dc2$wl, dc2$cnt, col = 'orange', type = type)

    Anh3 <- read_cal(FileAulaNH3)
    Onh3 <- read_cal(FileOldNH3)
    An2 <- read_cal(FileAulaN2)
    On2 <- read_cal(FileOldN2)
    # plot(Anh3, xlim = c(137, 200), ylim = c(0, 400))
    # unterer Bereich des Fits ist in einem sehr niedrig count Bereich


    # ymax <- 100; xlim <- c(137, 200)
    # ymax <- 500; xlim <- c(180, 200)
    # ymax <- 1500; xlim <- c(180, 210)
    # ymax <- 5000; xlim <- c(180, 210)
    # ymax <- 15000; xlim <- c(180, 215)
    ymax <- 25000; xlim <- c(180, 220)
    par(mfrow = c(2, 1))
    plot(Anh3, xlim = xlim, ylim = c(0, ymax))
    lines(An2, col = 'red')
    plot(Onh3, xlim = xlim, ylim = c(0, ymax))
    lines(On2, col = 'red')

    # -> rohe, unkorrigierte Spektren noch vergleichen!!!

    Anh3 <- read_cal(FileAulaNH3, correct.straylight = FALSE, correct.linearity = FALSE)
    Onh3 <- read_cal(FileOldNH3, correct.straylight = FALSE, correct.linearity = FALSE)
    An2 <- read_cal(FileAulaN2, correct.straylight = FALSE, correct.linearity = FALSE)
    On2 <- read_cal(FileOldN2, correct.straylight = FALSE, correct.linearity = FALSE)
    # plot(Anh3, xlim = c(137, 200), ylim = c(0, 400))

    ylim <- c(2000, 2700); xlim <- c(137, 190)
    # ylim <- 500; xlim <- c(180, 200)
    # ylim <- 1500; xlim <- c(180, 210)
    # ylim <- 5000; xlim <- c(180, 210)
    # ylim <- 15000; xlim <- c(180, 215)
    # ylim <- 25000; xlim <- c(180, 220)
    par(mfrow = c(2, 1))
    plot(Anh3, xlim = xlim, ylim = ylim)
    lines(An2, col = 'red')
    plot(Onh3, xlim = xlim, ylim = ylim)
    lines(On2, col = 'red')
    # unterer Bereich des Fits ist in einem extrem niedrigen count Bereich

    lin.coef <- Anh3$DOASinfo$Spectrometer$'Linearity Coefficients'
    linearity.func(c(100, 1000, 2000, 5000, 10000), lin.coef)


    # hier bin ich!!!!!
    # -> to do: check Abzug blinds anstatt pixel_straylight, ev mit linearity vorher (oder gar ohne lin?)...


    cheng <- -cheng2doas(get_Cheng(), nh3)

    ref <- nh3
    cheng <- get_Cheng()
    ftype = 'BmHarris'
    fstrength = 25
    fwin = c(190, 230)
    fitwin = NULL
    shift = FALSE

    n2 <- read_cal(FileN2)


    chen <- get_Cheng()
    cut_wl(chen)



    # dat <- read_data('~/data/projects/01_MD-calibration/data/S5', '22.02.2021 20:00 to 22.02.2021 22:00')
    dat <- read_data('~/LFE/02_Daten/1_MiniDOAS/Calibration-Aula2021/S5', '22.02.2021 20:00 to 22.02.2021 22:00')
    dat2 <- read_data('~/LFE/02_Daten/1_MiniDOAS/Calibration-Aula2021/S2', '22.02.2021 20:00 to 22.02.2021 22:00')
    avg <- avg_spec(dat)
    avg2 <- avg_spec(dat2)
    # sub <- cut_wl(avg)
    # get_wl(sub)

    # plot(avg)
    plot(avg, log = 'y')
    lines(avg2, col = 'blue')

    plot(avg, what = 'r', col = '#00000022')
}

### get data to script
library(sodium)
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
            pasteN(xs, 30), ')\n'), con)
}

if (FALSE) {
    cheng_data <- data.table::fread(
        '~/repos/3_Scripts/4_MiniDOASAuswertung/ReferenceSpectras/Literature_absorption_crossections/NH3_Cheng(2006)_298K_140-230nm(0.02nm).txt',
        col.names = c('wl', 'cnt'))
    add_to_eof(cheng_data, 'cheng_data', '~~~ Cheng absorption cross-section ~~~',
        'get_Cheng', '~/repos/5_GitHub/gel-scripts/minidoas/miniDOAS-functions-tools.r')
}

# ~~~ Cheng absorption cross-section ~~~
attr(get_Cheng, "cheng_data") <- string2dat(
c("aa","e2","9b","83","f0","2f","cb","aa","ee","02","e2","cf","63","a4","8e","ff","f9","f6","5c","9c","ab","ef","39","07","dd","7f","3e","39","6a","07",
"fd","6f","de","c3","b6","cd","e5","3b","75","86","94","f5","ff","b1","55","3c","09","16","c2","ae","8c","88","df","18","7b","d6","b4","db","3a","36",
"32","a7","6d","37","da","05","90","c1","48","f9","02","90","99","ae","ec","2c","2d","be","4a","ed","33","10","3f","db","d6","2c","a3","ec","78","62",
"cc","e7","ea","e7","49","7c","0e","ac","eb","4a","77","ca","7f","fb","83","56","37","e4","0c","03","95","c4","65","48","72","3d","55","78","7e","6c",
"03","70","6a","4c","16","c9","e6","d9","16","99","58","36","ee","d9","05","36","60","da","b5","a7","a6","8b","00","a3","bd","ca","77","7d","94","91",
"17","f6","72","64","71","65","f5","fa","a6","9d","57","2b","72","cd","51","fd","54","6b","96","03","a4","53","ea","4c","86","12","08","09","5c","6f",
"e1","83","f6","61","f5","27","20","a2","04","ee","67","04","2e","fb","04","24","28","3c","ce","0a","50","ec","0a","7e","67","d3","e8","89","e0","99",
"a9","d8","c5","27","91","1d","bb","e8","ef","72","e8","68","a3","8f","1c","69","60","4d","5a","16","4d","f4","b4","51","40","33","10","01","f1","41",
"d5","d8","1a","e0","f9","e2","bf","70","fb","c5","ac","0d","d2","77","46","25","7d","b8","0b","97","fb","75","ec","21","01","7f","1d","4d","86","37",
"7b","8f","4e","b6","f4","5d","7c","f7","ef","36","d1","8f","fa","d8","51","8b","9f","d4","70","18","f1","a1","b5","10","da","c2","4e","6b","d3","0d",
"8d","d3","d7","cd","0d","35","f1","c6","e4","48","0c","6b","29","19","42","1a","40","65","c2","76","42","f7","a4","7b","ca","2e","16","94","88","14",
"1b","1b","bc","f3","1b","27","8a","00","b4","20","9f","5f","12","e7","be","8b","2d","80","df","51","ab","21","11","5a","2f","97","36","57","35","e9",
"65","70","50","09","1d","19","b2","87","21","2a","0f","18","7c","dd","50","cd","1c","5a","dc","03","58","8b","12","13","a7","90","bd","50","64","d2",
"d0","cf","46","41","bf","d5","90","9e","21","e4","88","de","b2","7b","df","47","58","c4","54","91","ee","ca","8f","ca","c5","e6","db","db","48","19",
"22","49","96","dd","f6","30","f6","e1","56","da","d2","4c","58","ec","d9","98","1a","d1","27","01","f0","dc","49","0a","6f","40","6c","d4","7e","9c",
"03","c1","fa","88","dc","39","32","71","10","63","40","2d","9f","0e","e6","d9","9a","0b","95","15","14","8e","f9","60","c6","28","de","7b","9b","22",
"11","ce","08","ed","a9","51","2d","c7","2d","4d","6a","33","53","25","39","05","75","e9","25","05","d6","6a","6a","ae","9e","85","b1","f9","ff","0d",
"a4","9f","44","3f","4d","6e","cc","e7","b3","f6","c3","9c","3c","07","c3","68","33","e4","4a","9e","3e","0d","da","71","74","ce","84","a2","94","28",
"78","91","29","ed","f3","77","13","14","6e","57","d3","11","8f","b5","e9","a3","cc","fb","87","c1","3c","24","12","de","ae","00","51","97","23","a2",
"e6","9e","b0","a5","25","ce","b0","80","d8","f4","ea","2b","27","c2","d6","20","21","ec","8d","7c","f3","92","20","fd","ad","1d","0e","20","40","68",
"fd","40","b1","55","a0","64","42","44","77","21","2f","5d","42","4b","a5","1d","13","21","28","0b","75","f6","35","3e","95","01","cc","56","44","d3",
"69","9e","94","4b","88","74","f2","97","9a","c9","98","4c","a9","e5","3b","4c","d4","12","52","e1","ee","ac","e0","cf","75","c2","7b","0a","bc","cb",
"87","a6","03","b6","ce","5a","55","c2","23","e6","4e","70","42","c0","84","a2","49","3f","b1","1d","fd","f8","aa","45","b3","f1","b1","07","29","29",
"8a","15","38","7f","b3","23","f8","55","6a","64","e6","43","54","a3","01","07","1f","30","5a","62","f3","bf","de","e8","06","f4","bc","9b","62","b5",
"f6","c5","21","00","4e","c0","24","5e","95","19","29","67","69","20","99","3e","3c","71","46","e3","90","86","2b","90","ce","f9","00","1c","be","c0",
"71","59","c8","6f","0a","e1","36","41","f5","6b","80","55","75","eb","f4","be","cf","5f","74","92","1e","22","0e","c8","7c","4d","53","ed","74","15",
"c6","6a","25","25","e2","0a","b0","6b","d1","be","c2","dc","c6","0d","ea","47","bb","a6","bc","82","ff","98","63","4e","26","cb","29","a1","04","29",
"bb","d6","5b","84","9d","04","f8","26","54","19","47","0d","74","82","f4","14","6a","26","7b","4e","28","08","af","00","6f","19","47","77","f3","9e",
"fe","5e","dd","37","e0","e1","8a","4f","f5","e4","69","0e","0a","af","7b","36","56","9b","22","22","61","9d","49","7f","09","96","39","1f","01","de",
"a4","63","7e","3b","d0","7d","df","ab","86","0d","5b","32","c9","af","cc","7a","1f","36","c9","97","b4","b4","8e","f1","22","55","4e","7a","08","f5",
"0b","71","bc","f6","fc","dc","49","c2","dd","61","b4","a8","b9","ad","8f","c6","7e","8e","ce","ed","38","66","2d","c8","73","42","63","fb","86","10",
"2b","bf","a6","08","96","c7","d8","c8","ff","7d","d4","af","00","9e","32","33","72","a2","b4","ec","ea","ea","8a","58","89","7d","92","92","66","d2",
"19","71","39","9f","e9","f9","17","8a","de","d2","60","bd","8b","2d","52","46","17","41","8b","ef","a5","ae","3e","ac","0a","69","fd","93","a6","00",
"41","8d","e4","6e","86","60","da","bd","05","06","23","c7","69","dc","40","68","c5","30","b7","18","fc","c4","c7","24","b2","88","1d","85","da","17",
"bb","4e","d5","e4","9f","73","e2","25","3e","5b","43","cb","f9","aa","cb","10","58","c2","ea","e2","36","ce","e6","e9","97","6a","1f","f4","81","d8",
"c5","58","68","8c","85","bd","83","4a","99","8b","65","90","c1","4f","71","24","c3","42","2c","b6","e7","8b","ae","3a","47","5d","f1","4f","44","fa",
"9a","91","85","1d","c4","8d","4f","b4","ee","1c","96","4f","33","87","95","b1","01","d6","b6","d5","13","0a","9c","6f","2c","d2","b3","ee","a3","70",
"54","76","b6","58","97","58","cb","4d","9b","e2","02","0a","49","bb","5f","01","11","e3","b7","32","46","e3","65","7b","51","c7","81","37","f9","5e",
"62","c5","1c","43","0b","e7","1a","23","c9","07","ff","28","08","dd","af","8d","fb","a5","e6","88","06","bd","ce","ee","25","39","d6","05","fe","17",
"3c","fa","5a","85","e1","bc","e2","e9","a5","3a","f9","33","49","9b","d4","49","d2","ff","bf","b3","49","62","ed","69","b6","5f","65","8b","71","16",
"36","e9","20","99","95","56","97","b8","f9","28","61","a7","27","7c","5c","27","8e","d5","1e","5c","86","80","e5","f3","76","ec","bb","40","45","cf",
"8c","9f","17","c7","19","54","2c","e5","28","36","a9","dd","70","f0","40","11","ab","7a","77","38","b0","f3","96","05","88","af","0a","28","00","1c",
"1d","22","44","cf","5e","e2","51","6b","19","f8","50","9a","4d","c4","77","ac","de","0a","88","0f","d4","8b","5d","63","51","51","14","b8","83","bb",
"12","e9","22","cc","41","74","9d","ea","79","6a","7b","4a","d0","c5","42","b3","8f","81","83","bf","1a","90","11","be","52","b5","cd","ee","60","9d",
"db","33","c8","4c","87","3b","d2","10","74","0d","53","97","c0","61","93","ca","ea","2d","bb","4c","f8","86","85","0a","d0","4a","16","22","a1","23",
"e4","2e","69","de","ae","85","99","85","1d","2c","49","52","39","f8","d8","c9","ee","c9","0f","2a","07","d8","ef","ab","32","b7","f7","30","10","be",
"e1","6b","a5","b4","ad","bf","00","1b","b8","4e","02","39","6b","4c","85","9d","19","0d","56","30","ae","2b","12","0b","8f","51","06","d1","6e","51",
"28","f1","84","96","a4","6b","38","d9","4d","8b","7f","3c","3e","c7","c1","04","bf","23","e2","39","4a","3d","b8","15","b7","65","16","e8","62","b3",
"15","83","a6","39","c5","59","f5","38","2d","9b","3d","2b","16","04","d1","23","0c","4a","46","e7","d9","86","b5","c1","45","20","83","7a","a4","a7",
"0e","6c","52","c5","d0","8a","8a","e8","ee","50","f4","e5","22","1d","fd","32","b0","5b","01","a7","ee","dd","3f","0a","f7","33","40","7c","73","d2",
"16","1c","a6","fe","dc","f9","3b","87","2a","8e","d0","ec","f3","b8","c7","fa","26","60","9b","57","a8","7d","b5","02","a9","40","36","95","9c","db",
"71","6c","c7","1a","d8","e0","1d","c4","67","02","18","2c","38","22","ae","32","fe","d0","d8","62","01","fb","0f","b5","36","4a","30","66","61","2d",
"5a","6e","92","8e","ec","b2","4b","d1","ff","e0","7e","51","27","ee","c0","a9","dc","aa","de","bb","32","ec","88","89","b3","23","73","04","7d","b4",
"40","c6","4e","9d","a1","4d","76","87","f8","ef","0d","43","b4","56","09","1f","10","2a","e2","5d","96","53","9b","58","47","41","e3","75","4f","6e",
"7f","52","5f","3d","54","3e","0b","06","3e","e3","c7","0d","06","43","1c","67","e1","da","2c","14","f5","86","91","f2","72","9e","8c","fc","15","cd",
"af","9d","63","63","95","eb","03","bc","0f","e5","58","95","59","8e","39","9b","dc","c6","bd","78","d7","8a","14","ee","45","25","76","87","60","65",
"2b","82","18","2a","58","2e","d1","de","75","73","48","39","c5","d3","08","59","07","21","f3","a3","f2","17","e9","cc","61","ac","91","45","90","a0",
"68","8e","64","7f","e0","28","20","cf","db","95","d6","02","96","11","4f","1c","04","2b","c0","1c","05","81","37","49","08","6d","0c","28","71","74",
"63","cf","f1","f0","b3","38","27","a5","0e","84","12","2e","e6","09","dc","07","5a","f3","c9","87","80","2a","5b","51","fe","4c","19","d2","74","34",
"cc","3f","fd","7f","e1","d6","65","43","8e","46","ae","42","35","d2","9e","08","90","5f","36","33","dd","94","8a","e8","44","45","85","af","d8","a1",
"de","37","ad","f9","7a","21","8c","02","d6","aa","d6","d9","fc","de","a2","10","8b","fe","6c","29","9c","0e","c7","4d","7e","51","7a","e1","77","dd",
"98","dc","d8","00","98","02","af","bf","e7","94","e9","46","eb","8b","b6","5f","d3","a9","89","29","da","4d","ae","54","89","98","3f","ed","2c","6c",
"5d","4a","79","38","71","40","1a","e9","72","79","90","aa","aa","aa","49","ce","00","bf","81","d8","5c","56","c6","d5","96","b7","80","e4","4a","1b",
"d6","ad","87","8c","16","6f","52","e9","15","15","70","7e","4e","59","61","76","26","d1","43","77","81","89","3e","f7","92","fd","8e","18","19","76",
"de","d4","4a","2f","99","0c","50","69","22","cb","1a","34","ee","3e","4a","a0","f7","d8","14","73","2e","89","9a","97","39","17","50","19","25","05",
"a4","30","42","be","23","da","e5","31","1e","ae","38","4e","03","c3","4d","84","8c","be","cf","8a","b2","33","3d","d4","71","db","9b","1d","40","b7",
"69","1e","f8","f8","0a","d6","c6","25","50","39","0f","d2","98","07","86","1a","36","49","bf","f4","d5","d6","da","44","84","cd","bd","ca","94","41",
"fb","41","71","b3","e1","9a","d6","06","fb","42","69","c6","9a","f5","d6","a0","c5","5f","16","56","55","8c","50","20","6d","8e","03","da","38","06",
"e7","8f","cf","8b","42","50","36","22","87","87","50","44","a8","80","3e","f3","81","0c","57","24","95","19","de","0c","e2","cb","ae","ef","7b","12",
"6e","dd","80","1d","c5","98","97","69","6e","2d","e3","b6","db","49","4b","7c","5e","f7","99","47","2b","46","80","e4","bc","d1","05","4f","8a","58",
"85","fb","d1","e8","69","1d","f1","64","28","55","e6","d6","e4","a7","9b","0a","91","1f","2d","e2","4b","21","74","f1","16","ea","83","39","4a","f8",
"ee","1f","69","be","36","d7","71","77","f8","95","3b","e3","0e","99","cd","01","59","fe","f4","13","4a","1f","c9","59","e0","3f","7a","16","df","d7",
"74","a8","ff","0a","d8","25","4d","4a","c3","63","f7","a3","67","a6","2b","2d","23","d1","1c","a4","9e","44","7f","6d","80","f7","23","46","14","0e",
"df","3c","5d","6b","1f","f9","9f","c4","1e","6f","6a","ff","8a","fd","60","45","5d","dd","7c","43","77","1f","78","99","63","a3","cd","57","40","0b",
"a2","1e","32","44","43","e6","90","53","8c","aa","0f","4a","7f","a2","59","87","62","f9","c6","e4","c7","34","2a","ea","e5","88","47","85","6b","2d",
"5f","7e","27","98","30","58","dd","76","64","23","d2","20","75","8c","f4","a7","a6","94","ea","0c","49","50","a8","d7","d6","38","a3","df","a1","32",
"7c","bf","eb","c3","7c","be","81","aa","e9","69","ce","47","60","4b","ef","00","45","88","77","77","8f","40","e8","71","e7","7b","67","3c","e4","d5",
"23","84","d9","31","54","16","b4","bb","01","e8","78","d8","00","4d","fc","f6","98","55","48","79","c7","f0","19","1d","01","9a","12","10","2e","4d",
"d5","eb","55","f5","6c","80","6b","9b","9e","4a","6a","50","ff","2c","db","9a","99","67","b8","b0","de","76","3e","06","20","ea","a1","13","68","bc",
"49","50","8f","7c","5f","e3","b3","61","7a","b7","ca","e7","3c","67","c2","3c","24","00","5d","d6","0f","ff","04","aa","14","37","80","d2","d8","b0",
"37","cb","12","37","b5","26","49","45","7a","0c","16","56","a0","3b","32","2a","ac","c1","68","41","87","db","ce","c7","cd","e0","2d","6e","4b","d9",
"7b","7e","69","0e","c3","2b","e1","32","0f","4a","1d","83","42","3a","b4","34","32","55","c2","fd","8b","fc","b9","5a","e9","bb","b9","29","bd","bd",
"99","ca","e1","bb","aa","03","7f","f3","d5","b5","0c","72","e7","f5","97","18","33","f5","f6","b9","e5","f8","75","6d","52","56","20","9e","a5","2e",
"08","04","29","c8","3e","4c","a4","c5","a8","4e","f3","8d","f7","fd","a3","c8","23","d5","bc","fe","30","9c","00","cf","b0","9a","bb","94","09","f2",
"89","c5","6a","11","53","31","cb","d4","5b","0e","b3","21","2f","5a","b5","53","9d","f1","0d","7b","a7","f1","02","00","78","ce","18","d3","1e","4f",
"a9","75","ce","8c","f1","b8","79","e0","8e","10","90","e7","b5","91","f4","66","a9","0b","a5","62","07","0c","48","f0","be","91","d6","a6","01","3b",
"72","96","ec","f8","e4","16","e8","40","1a","a3","b9","8c","ad","46","bf","ec","ee","d6","cc","e1","a4","81","b6","29","44","81","76","30","85","fb",
"f2","a3","65","18","4f","e7","75","18","b7","3d","ad","56","51","62","82","bd","9a","cb","48","0d","65","f4","0d","10","51","8b","d7","63","38","17",
"e4","9a","53","72","bc","a4","a5","a8","cd","a3","74","61","9d","e6","12","ab","82","95","64","0d","6a","e2","f2","92","68","f8","b0","56","b0","ff",
"5f","e8","75","88","9d","98","01","ca","bd","8b","82","0a","f7","28","3a","04","6f","9f","f7","54","26","52","13","50","03","4c","5b","4e","c4","39",
"94","cc","9c","a9","5e","f4","3b","1f","3b","c2","4c","9d","fa","61","6c","91","04","b4","97","06","07","0c","6a","67","6e","22","79","14","7c","a9",
"ad","f3","fb","eb","66","bf","95","4a","39","ef","d1","35","1e","35","ea","4b","69","05","20","37","c1","95","5a","8f","d1","c1","1b","d1","3c","ca",
"b0","f2","05","9e","27","04","c2","6a","32","07","4c","c0","be","32","f9","be","a9","86","fd","67","04","59","e1","2b","9b","d4","2e","5b","25","77",
"82","0a","e5","e5","7c","b9","b9","ff","77","3f","f8","41","5b","c8","27","95","23","f9","75","7e","17","a0","5c","89","04","47","fb","f7","6b","f2",
"df","2a","b4","63","fc","e9","62","9c","8d","a8","33","38","a9","d8","55","8a","86","be","a2","ba","9b","e6","69","7d","21","d2","1a","44","63","99",
"db","23","43","c7","42","a8","b3","1b","a4","5e","7a","93","39","9b","45","6a","96","bb","d8","31","df","ca","22","b5","79","bf","a8","eb","60","28",
"ba","f5","3c","47","57","8d","70","96","8f","b3","06","d3","30","14","0a","2e","ce","be","de","b1","c9","5d","a4","d9","b5","4b","3f","66","02","8a",
"91","f6","ca","8d","56","33","12","d5","bf","2e","71","a3","a5","2b","15","cd","c4","c6","40","f0","4d","14","a9","4b","67","4a","6f","82","9c","a5",
"38","1f","a9","28","ba","d3","bf","e9","05","08","c3","b8","1a","47","7b","88","7f","20","22","36","91","ba","6a","bd","44","8d","ef","45","75","42",
"b1","65","5d","df","6c","63","e8","d6","66","47","c1","46","66","05","a0","ce","95","81","56","72","df","ed","4d","80","34","67","96","16","aa","72",
"10","58","44","6b","15","cf","9e","a7","c8","3d","0e","27","16","51","12","d3","f0","f2","9e","d6","aa","53","0e","36","c3","59","84","7e","1b","c7",
"1f","0f","92","ea","93","a6","7e","3d","b4","2e","a0","a0","17","f2","7c","a4","86","eb","72","46","33","01","11","1a","df","21","6c","bc","d0","f5",
"c4","48","c3","8d","0e","df","08","b2","de","bb","a5","ab","02","c7","66","1c","b3","e0","49","03","bb","a7","41","67","dd","b8","d3","2d","6c","aa",
"b7","b1","3f","d9","57","55","55","ab","61","7d","fa","f5","c3","1a","45","29","25","df","b0","2b","34","cb","a4","83","8b","61","82","21","04","c3",
"11","a7","55","e1","63","09","b7","d5","d3","87","3e","9a","24","30","04","0b","20","3c","1c","cd","ea","e7","90","0d","2d","4b","64","6c","52","e0",
"a4","bb","6c","8e","47","b5","4e","5f","67","53","e3","30","f0","42","c1","3f","9d","de","2a","b8","9b","5f","09","bf","f6","85","2f","a8","a4","92",
"c4","01","57","ef","8b","b9","9d","a9","56","84","c1","4d","0e","d2","15","fc","d8","b1","6f","1c","a5","81","c6","b3","4e","1f","5d","d0","4d","5c",
"db","f2","21","60","38","b4","12","47","2e","99","c0","d6","ce","d0","68","b2","28","cc","a9","01","3e","e3","2b","a8","74","25","41","7b","05","ea",
"c9","b8","12","5c","bd","52","36","21","13","6d","c3","8a","bc","1b","53","11","f2","43","b3","1c","aa","94","d6","43","92","56","b7","7c","1c","e6",
"42","99","6e","20","66","51","96","1e","57","73","e9","ae","dd","60","31","67","8c","97","56","32","90","c0","b5","18","ae","ea","04","a6","a3","8b",
"0e","df","00","e7","d8","13","e5","a8","4f","80","e8","6c","f7","26","19","9a","23","f1","6f","78","29","f3","5e","95","14","6e","67","ba","86","bf",
"c2","58","f0","ab","8e","72","93","be","1c","d8","21","20","2e","54","3a","ee","42","38","ea","bc","6b","a3","c7","e8","00","9b","e4","2e","1d","01",
"d3","a7","3b","82","7f","06","03","7d","49","45","e1","87","8d","b0","34","73","fc","70","be","ce","68","92","66","db","c1","5e","04","14","83","23",
"12","01","d6","92","8d","29","17","29","49","d8","b1","b8","c1","9f","ed","a1","17","16","25","de","bf","76","4e","bb","9e","d9","39","81","12","91",
"23","30","18","1e","00","70","91","16","4b","af","a9","36","7a","51","57","6e","b1","31","9a","b2","47","ce","38","54","d4","ec","d5","3f","a5","2d",
"48","fe","c8","b8","52","72","87","c6","72","3b","e0","7c","7a","f8","75","f0","ec","87","3a","cd","c8","1f","58","41","f0","da","56","a5","0e","4f",
"b4","f1","de","70","2b","da","d7","3f","bc","78","ad","1e","17","d6","d9","6c","ff","b2","66","fc","87","2f","64","89","4d","7a","3c","f4","2c","78",
"5c","78","2a","b5","22","f3","50","5d","a3","d2","1c","88","17","2e","2c","74","47","57","15","d7","6a","c6","39","e0","01","1e","7c","a4","45","b2",
"da","b2","b1","e6","e2","89","87","a0","73","19","33","96","b1","68","af","9a","25","2d","ae","3a","bc","33","2e","75","55","8b","59","9b","36","5e",
"5b","bd","9c","ba","34","d0","8b","7d","df","d7","27","bb","b8","fc","d6","3e","ca","84","92","0f","11","43","12","0c","95","64","d5","65","0d","8a",
"1a","89","62","be","ca","a4","34","88","71","f5","3d","d3","f4","59","47","f6","5d","3e","4f","47","e9","48","5e","11","9a","01","83","c0","ed","f1",
"89","d7","c8","23","bd","c6","59","a3","45","74","77","41","03","35","1f","8e","d2","f3","54","00","0b","ce","f4","1f","89","6e","88","d0","21","05",
"41","0a","c2","fc","3b","3a","eb","51","ca","7d","ff","e4","c2","57","16","18","69","f2","5d","ad","87","68","71","41","fb","bc","57","16","e6","d1",
"fe","68","e8","2e","26","72","a3","81","e8","2a","0c","f5","a6","e4","af","f3","03","08","b5","43","05","51","7b","42","f7","6e","b4","8d","ed","53",
"78","f8","08","68","69","b0","e3","ed","8a","a5","0a","22","66","e9","b7","9c","52","90","1b","49","3b","cf","a4","3b","fc","d2","3b","bb","16","d9",
"3b","ac","82","c0","0b","7f","f7","52","d3","67","ba","d3","72","7c","3c","6f","76","5f","c7","0f","6c","dc","9b","25","34","14","42","b8","b5","73",
"fc","5d","62","48","50","63","17","30","b4","de","c9","d1","2d","73","10","13","b2","1a","de","55","55","6c","5a","70","ce","8d","9f","4c","20","76",
"80","6c","5e","08","06","d1","84","f6","69","50","9a","58","3c","72","96","1b","be","9f","2b","90","0f","c1","e1","aa","38","7c","f6","a3","09","24",
"42","7e","61","63","81","7e","e7","05","9d","db","b0","b1","1b","b0","ea","95","d1","1b","28","ed","2c","d2","42","3f","01","78","17","49","dd","de",
"f0","97","77","7c","30","70","18","8c","dd","80","9e","3f","eb","e5","22","bd","f9","2a","ce","40","fd","91","d2","42","06","fa","a7","1e","82","47",
"ef","ef","a1","4b","d6","76","39","2f","b9","a1","e2","5b","30","7f","eb","a3","46","bc","52","ea","71","78","44","23","97","ff","df","b7","e7","5b",
"3c","ec","89","7b","4b","88","3e","01","5d","fa","1d","5c","24","8f","8d","73","32","7d","cf","39","24","d3","c4","8d","94","7a","01","04","f0","fb",
"76","54","7e","e2","d7","cf","61","6a","9e","1b","bc","cd","f5","47","e7","6b","4b","e4","d6","ed","55","6f","d3","78","ed","a5","a1","f3","70","fd",
"0b","92","52","57","b2","4e","9d","18","e2","1f","80","d1","64","89","83","fb","bf","af","1b","90","22","4c","3d","1e","31","ac","64","76","aa","f5",
"c8","75","cc","a0","1f","31","a2","d3","60","a1","9f","55","13","bf","9a","cb","cd","5a","03","ba","09","e9","7d","5f","64","fc","ee","dc","f5","29",
"b0","d5","19","b6","e5","7a","d8","56","e6","18","60","62","ca","41","74","d8","24","65","a5","9f","c5","e5","13","5a","c0","65","12","b2","40","ce",
"c2","26","4d","16","12","75","29","9a","ac","b2","1d","88","5a","39","f6","b9","7a","82","56","6a","0c","ae","2f","6a","b3","b0","6b","b9","19","d1",
"f8","8b","5b","cf","0d","9f","aa","f1","b2","e4","9b","dd","68","6c","9f","51","04","e4","75","c9","33","9a","3b","fe","05","21","a6","2d","1a","56",
"01","a2","79","ff","d6","7c","b9","f3","b8","3c","d8","85","6f","ae","c2","ca","30","a8","7c","54","c2","c0","74","3c","bf","55","80","70","b1","6b",
"a4","e6","2c","6c","a3","49","99","60","eb","a5","78","b8","94","1b","6e","60","c0","b2","b4","0a","65","df","25","49","d3","41","3f","19","1e","a4",
"40","0b","41","71","43","b3","7d","65","21","55","97","de","b8","4b","9e","d1","df","26","4f","6a","39","87","d3","fa","4b","7a","05","a7","33","b1",
"94","8b","df","e4","6b","ea","a2","74","f5","77","f3","23","79","75","4a","9c","db","0c","6e","96","94","5d","de","6e","a5","7b","e5","36","81","f4",
"98","e7","70","9c","56","14","e4","b4","26","ac","ea","f2","9e","68","71","da","25","db","bf","d3","68","ab","94","6c","78","66","fc","2a","11","6a",
"aa","b9","cd","f8","8a","d6","ed","c5","9e","6d","a8","89","75","23","f9","28","db","77","5e","be","fa","f5","02","cd","c8","d9","69","1d","37","cb",
"15","68","90","c3","2f","25","da","c8","94","25","b3","c8","82","51","76","11","55","03","c1","5b","24","b0","94","a2","c3","aa","ca","b9","73","76",
"3f","8f","d0","52","df","2f","6d","e1","77","1c","d0","32","6b","49","78","e5","83","45","58","a4","73","ec","30","40","59","1e","f7","f2","cb","d6",
"3e","d6","40","2a","e7","58","79","e8","14","95","24","49","ba","af","2e","25","58","1b","7a","e1","95","10","1e","38","fa","00","06","92","bc","21",
"22","b7","32","96","4b","d9","aa","46","02","6b","de","0b","9b","89","ba","62","94","36","70","09","c2","f0","d3","d2","a8","a1","56","58","1d","93",
"14","8c","ce","55","61","a6","b6","e1","b9","0a","17","bc","41","2c","d0","d3","89","83","d7","b3","d3","0e","30","c3","cf","19","f7","49","4b","44",
"6a","9f","d4","0c","6b","28","6d","a6","ab","77","48","ba","20","be","18","2b","e7","ce","88","2a","c9","2e","ad","ef","86","cc","f7","1f","a3","30",
"47","2f","b7","70","5d","07","36","7f","24","dd","97","c7","1c","7d","a2","17","4e","a4","97","d4","fb","ff","94","0d","18","4b","9c","8d","40","f6",
"b2","b1","2e","ae","dc","3e","e2","cf","eb","f8","13","9c","88","2b","a3","d3","e4","3e","70","24","a5","23","0c","ce","d9","5c","a3","62","74","7a",
"b4","c0","4d","5b","df","16","4b","d4","e4","8c","9b","17","a4","24","86","4b","a3","cb","e7","ee","a1","64","d2","b1","03","dd","63","18","8e","ac",
"fa","40","c8","6e","15","3a","04","0c","2f","3b","aa","73","1a","2a","55","e5","e5","57","dc","4f","cc","61","7d","90","a9","72","db","68","38","db",
"5f","98","ca","84","52","92","de","23","65","bb","b7","6b","15","69","ad","28","87","93","fa","09","30","50","f5","34","28","db","22","9a","f6","31",
"d3","8b","c4","58","de","a3","0b","1c","03","af","05","1c","97","7c","db","8b","6a","2b","bf","39","b0","79","91","cc","4a","4a","d3","90","5c","70",
"02","96","ba","c8","9b","3e","07","a1","84","f6","0e","94","ec","3c","ee","75","48","6f","a4","c1","a9","b8","03","90","d8","01","1f","d4","3d","0a",
"43","5e","0f","97","d6","95","45","10","61","fe","27","e1","74","f2","15","fb","a5","fb","4e","4d","78","57","91","34","b6","18","aa","cb","89","9c",
"cc","4c","f1","99","fc","e4","7a","48","d4","b0","90","d2","f9","ee","df","cb","19","a5","ec","56","a2","3e","85","92","0e","ea","35","cc","2e","ab",
"81","d4","33","ee","03","b6","20","d2","a7","39","40","55","e2","76","74","1e","22","9d","fb","d6","81","ce","70","94","ca","3e","a8","5e","d8","5d",
"17","b4","b2","f9","6a","22","c5","87","91","c7","0e","4a","fe","cc","a2","d7","a0","42","3f","8d","08","01","38","1b","7b","8f","0c","f0","b3","63",
"dc","43","ae","ef","21","d6","01","56","4a","3a","1b","6c","42","be","ae","dc","ef","60","d8","61","c1","2e","dd","53","57","2d","3b","87","21","78",
"86","bc","82","78","5b","55","1b","a7","84","7a","3b","da","36","8c","ff","01","27","32","14","43","b5","f1","4d","bc","52","e4","00","81","bf","ac",
"ec","fa","ef","30","9d","9b","0f","ce","63","b3","06","92","52","a5","79","e0","d4","63","ec","10","70","6d","1a","63","e6","03","6b","be","44","9a",
"d4","87","9f","79","6f","08","19","01","76","3a","eb","3f","8b","99","1e","52","8b","0e","e9","f5","ea","31","45","9e","f4","45","58","1a","ea","03",
"b8","81","5c","c1","e3","83","d6","d0","3a","fc","38","86","1b","48","f9","84","21","1a","39","37","77","81","f5","ad","02","38","92","d0","13","3d",
"0d","03","05","2b","1b","14","e1","8a","d1","66","05","cc","53","b9","48","28","10","e0","43","3d","65","70","c7","90","b3","47","69","c3","de","ac",
"20","82","25","70","a2","7d","63","0e","1b","44","8c","2f","7c","2c","ca","d2","60","54","b4","09","22","db","9b","41","66","8a","33","d5","93","a4",
"c9","18","d7","01","c7","02","6a","b0","b0","a5","b7","9d","89","c4","96","89","e7","2c","43","74","2e","c3","9a","da","b4","f9","ee","4a","43","d8",
"67","6e","b9","3e","67","ef","82","75","ca","43","8a","dd","36","45","30","11","a0","78","35","f1","a6","ed","20","38","9b","0d","95","fe","77","07",
"17","c3","7b","a9","e1","e8","d2","36","57","f6","3f","bd","12","3e","7b","12","5f","fa","5a","68","72","16","4e","93","7c","a5","39","68","57","49",
"6b","f6","54","ec","01","bd","31","36","58","5b","1c","f9","7c","36","a9","d8","e5","3a","37","7f","e5","ee","38","df","8b","52","5f","2e","df","96",
"5b","37","ea","70","46","9b","b2","20","90","38","5e","37","3b","2d","65","50","c1","05","bc","64","01","63","69","aa","0d","b4","69","93","b8","55",
"ac","77","82","11","76","4d","e4","e3","f4","48","e8","09","b1","91","37","34","d6","47","33","6b","ff","ae","51","e0","a3","87","f9","0d","0c","e2",
"77","d1","00","aa","1f","f0","a4","25","00","d9","b4","35","1b","71","0a","8f","88","5d","fb","25","f8","6b","5a","b2","f2","6f","dd","81","aa","f6",
"ee","38","09","b0","1f","e0","0a","24","93","0b","bf","11","ab","b0","19","ed","36","14","6b","6e","43","5c","1c","aa","05","42","f8","6a","db","2c",
"aa","e4","8c","65","6a","77","39","8f","3c","52","9e","80","de","76","44","b0","c5","b6","ab","3a","c7","1f","cf","73","c9","5b","19","2c","0b","6f",
"ad","59","8e","ea","6e","8f","53","68","63","fc","3d","38","b9","e7","de","5f","87","c2","e8","5b","31","61","9b","21","74","23","69","7b","34","82",
"69","c9","e2","89","14","7d","59","e7","11","bc","fc","97","ab","eb","8a","1b","7d","d5","9f","66","38","3b","a0","0f","40","17","a9","cb","5b","79",
"ef","f5","04","60","1b","2c","5f","dd","c2","56","6f","e9","c7","87","38","2e","32","a1","0d","b0","47","ea","70","ab","8d","60","cb","18","97","ec",
"32","9f","9a","51","1c","3b","ec","2b","8e","5c","4a","dd","90","40","27","bd","f3","b4","97","dc","4c","59","11","52","86","62","68","e6","26","a4",
"d2","16","62","66","b4","3a","5d","bd","86","1e","79","6d","7e","42","95","14","78","cb","67","01","7d","97","25","8a","7c","de","39","7e","c8","9b",
"40","52","60","3d","2a","27","d8","ab","04","8f","b0","6d","3a","ce","05","20","bc","50","b3","11","85","e4","2a","9e","a4","e2","33","6b","62","24",
"8a","9e","55","48","59","08","70","53","60","04","5e","06","ed","14","d8","e8","6a","c8","3e","3e","43","60","df","2a","c9","c3","9a","2a","f9","1d",
"c7","22","d1","64","61","e2","cd","95","d0","0c","2f","de","ce","bc","5c","18","6e","5c","49","7f","dd","02","30","70","9d","9f","6a","c0","98","65",
"3f","e6","9e","7a","23","d7","6b","34","e7","bb","ac","f6","b4","f0","02","29","ba","a9","35","05","10","69","3a","2f","6f","aa","8e","5a","b9","27",
"aa","83","00","2c","db","af","1b","37","b9","66","9d","35","ce","15","57","38","fc","94","07","6e","85","1e","97","77","b1","d9","c2","82","4a","01",
"64","f8","b5","95","dc","82","be","8b","69","70","96","22","7f","51","08","1d","64","b0","6b","82","09","14","62","f1","f0","7a","df","b2","cc","a8",
"00","b9","82","e0","eb","ff","c1","c8","ba","1d","db","4b","85","47","ce","80","5a","21","f1","c3","bd","62","dd","d7","09","21","f1","76","61","f6",
"71","ab","b3","e7","fe","7f","ce","d9","ef","2c","ab","26","79","ba","cc","ec","fa","53","7d","b2","1f","cf","e8","f8","71","6f","d6","4f","40","c1",
"8d","19","76","69","64","01","93","f1","ca","2e","89","23","dd","e1","11","a9","1f","78","73","ca","de","90","8d","f9","a9","70","7d","0c","ba","0b",
"ec","de","4b","67","d4","12","55","10","e1","f2","7f","9d","98","ef","95","f6","10","64","b9","ac","90","c3","b1","03","92","a4","37","e0","ea","62",
"fb","5e","db","00","bb","ba","b6","e5","c2","7d","35","1e","af","d5","82","ad","8e","52","01","e9","e6","13","1a","56","21","4b","d8","16","b9","98",
"72","e5","82","ab","c7","81","7b","d8","10","b5","e5","37","dd","80","9a","43","83","4f","40","aa","4a","46","86","a5","5c","57","96","02","7d","4a",
"90","cd","7d","43","9e","4d","5c","7a","5f","ed","77","e9","8e","b6","f5","04","dc","5d","35","94","c5","b4","7e","b2","de","2d","ef","04","a1","4f",
"4b","34","ba","3f","c7","f6","09","25","05","91","49","53","d5","85","fe","d3","67","39","bd","fe","1a","90","cc","53","f8","72","14","33","eb","80",
"de","fd","0f","71","2e","29","9f","1c","5c","8c","4a","75","3c","4d","9c","97","22","a9","31","b2","9a","1f","87","9f","3c","6c","55","d5","37","32",
"ff","75","74","81","87","96","3d","87","42","9c","3a","61","72","a3","69","8e","48","bc","c8","cf","c8","41","b6","60","ee","37","82","e5","4a","53",
"1a","e5","16","63","29","08","b4","d4","ef","9f","c2","4c","32","a2","c6","c5","34","c1","f6","bb","c8","ce","61","cf","41","73","20","cc","69","9d",
"94","20","be","5c","e9","13","70","fd","67","f8","05","10","4f","b5","05","78","08","32","fe","7b","f0","68","8a","13","fb","eb","45","a5","5b","4a",
"2e","ec","4a","fb","45","fd","ff","47","66","e2","01","e8","0f","0f","8f","1d","fa","bb","00","66","51","5b","22","27","fc","53","1d","91","93","55",
"ac","21","eb","54","87","e5","4b","5b","7e","a1","ac","c8","7d","1b","cf","59","22","9c","9c","34","11","91","78","77","50","6c","54","ca","d2","ba",
"61","ca","04","6b","e8","cd","f6","ac","ea","46","33","03","74","97","5d","4b","21","9c","bc","80","75","aa","af","a1","92","4d","94","21","74","46",
"bd","f3","a0","39","e6","aa","eb","0e","79","cc","a1","90","0c","d8","35","41","af","ca","ce","fb","63","c3","a2","73","1f","94","0b","55","7b","a8",
"f3","c7","f8","ca","63","7e","7f","69","91","13","06","66","3e","70","7c","37","43","6a","70","42","01","af","3b","23","6d","f6","05","01","63","bf",
"c8","51","42","12","8f","af","b8","0f","d4","12","90","aa","f2","37","c9","4f","9b","d7","a2","b8","83","a6","99","b0","7a","1a","41","7e","74","60",
"40","64","ed","e4","1d","15","2c","4f","55","d9","bc","e8","f1","75","6e","96","c5","9e","cd","55","c3","7c","85","11","87","d1","ce","d2","16","66",
"1a","9f","d7","9f","d4","97","8d","7e","43","4d","2e","fc","c8","fc","95","17","a1","dd","cb","93","5f","4a","fb","9e","0a","4a","bd","9b","4a","a9",
"d6","0f","43","14","10","ef","32","52","c6","7b","63","bb","01","d9","0a","31","ee","e3","51","7a","4f","86","a0","67","6b","ee","8c","6c","53","31",
"c4","43","b1","ec","30","5b","a9","f6","dd","32","3e","08","3d","ce","7b","9f","13","e1","0f","61","92","22","c0","02","67","fd","d8","8b","a8","99",
"3d","98","d0","01","2e","ce","12","ff","98","8c","07","ac","ce","a0","19","ac","cd","40","a3","e0","44","0b","66","ba","a3","58","a6","84","83","38",
"82","4e","8d","ae","64","06","df","3a","c5","66","a6","6f","1d","b8","5a","76","2f","62","6b","5c","c4","79","8e","64","42","9a","7f","72","d6","b9",
"f9","8b","a0","34","92","68","02","fc","c9","b7","d6","dd","00","62","17","51","76","5c","63","94","72","ad","b5","e0","4f","dd","5a","92","1f","a8",
"c5","69","1a","20","1d","ce","b7","43","71","5a","9b","49","98","36","29","bc","70","45","cf","0d","d7","a7","fd","21","a7","51","47","3b","c4","87",
"14","e0","2e","c6","4f","11","9d","26","0e","a4","e2","a3","46","64","22","67","25","0e","b6","75","4c","85","05","94","38","9a","94","08","72","f2",
"be","2b","50","f5","38","c6","69","d4","a0","c7","eb","6b","14","1e","0f","43","e0","9c","8d","b2","1d","37","f6","4c","48","cb","cd","3d","3a","49",
"39","55","3d","b7","d5","29","9f","9e","cf","83","8b","b2","b4","9d","94","0d","2b","d5","0d","67","9e","58","16","4e","3d","08","46","28","e9","e3",
"ae","55","f6","ac","36","e4","96","ca","c5","f6","3b","4a","b9","60","6a","aa","d7","4b","dc","cb","18","c6","aa","55","0d","7d","f4","96","9e","50",
"d5","69","cc","78","5f","23","45","e8","05","b5","a6","2c","41","18","62","cf","8a","d3","e3","23","d1","0c","c0","a5","0e","29","3a","aa","f3","ac",
"a6","c8","7d","3d","1c","2a","4b","71","a2","ea","6b","cf","af","23","ad","90","ec","fd","83","f3","32","57","36","da","97","4c","dc","b8","c8","24",
"31","81","dd","35","86","25","d7","13","3c","96","71","0f","f0","0d","33","96","9a","9b","75","57","fa","b2","5b","bf","5a","5b","81","0f","07","43",
"60","26","e2","ed","dc","6f","2b","dd","e5","5b","58","2e","2e","30","24","9a","d7","51","bd","4a","21","df","76","f5","78","82","3d","fb","1a","ce",
"0d","48","1a","b3","ea","bb","24","b6","03","73","75","ac","e1","5f","18","43","1e","f9","ed","ed","68","a0","e3","41","b3","09","bd","2f","4a","4b",
"b0","6e","5d","c3","92","87","fd","6f","ea","b0","40","d0","d5","e4","23","b3","09","4b","fb","ff","5c","11","7e","b3","2f","74","93","5e","e4","3a",
"0f","31","b9","83","fa","3e","7e","21","dc","c5","08","65","b6","ca","24","03","37","fd","e8","77","fa","0e","21","45","1b","3b","a8","56","fe","b5",
"3e","da","87","21","16","ec","86","b3","08","3d","b2","3c","bf","25","75","a5","0e","83","d2","6a","19","cf","36","b3","e6","2c","62","28","80","f9",
"b7","35","fa","ef","7a","0b","e7","15","7b","28","9b","ed","76","bf","62","a2","fd","66","34","f0","1a","cd","d1","34","b7","a3","f3","27","58","21",
"d2","3e","58","90","1f","82","2f","a6","17","4e","74","0b","3d","d3","23","9c","13","da","1c","36","ce","41","23","90","c0","aa","3a","14","9a","49",
"f6","36","ec","a5","0e","d9","7a","7f","d7","a6","c1","7f","94","35","05","02","bb","78","d7","1a","b0","e9","23","58","14","16","3d","a1","48","c1",
"e5","bf","a6","26","34","ab","b7","c2","b8","82","47","f0","b4","73","81","38","bc","05","94","ea","72","a3","c6","80","ac","fe","39","85","d3","52",
"9e","60","34","a2","d9","2c","bb","93","ed","fc","ab","55","09","5c","f1","64","af","60","51","26","81","fe","02","e5","e0","c9","35","7b","92","94",
"0c","6a","b1","92","3a","95","2e","b7","b7","76","ac","6b","4c","a4","c8","ef","e1","36","13","95","eb","55","1e","fb","74","d6","2c","c4","cb","3d",
"2f","21","df","b2","39","a2","c6","c3","3c","a2","30","37","b6","c0","4d","c3","e7","c6","a0","92","0a","1c","57","3a","45","66","14","28","ef","e9",
"9c","ee","09","58","ce","80","75","69","05","3f","b3","4b","b5","ca","5d","0e","82","d5","eb","a5","f3","d1","fc","a3","96","6c","64","2e","c0","3e",
"aa","7d","c7","6c","37","71","da","d7","ec","5b","5e","c6","ca","43","4c","0b","91","00","9f","43","a7","06","7e","b9","0a","dc","0a","e9","87","53",
"f6","26","cc","43","dd","db","01","cd","3d","a6","fe","0e","5b","e1","0f","a2","3a","9d","7c","e3","2e","33","6d","6c","ad","2a","f8","9b","0c","66",
"a9","1d","fd","dd","99","d5","86","73","d7","2c","56","21","d6","7c","3d","3c","e4","73","8e","f0","5e","a2","de","6a","57","20","6b","c5","89","c9",
"5d","e7","01","0d","92","68","71","67","61","4f","98","c2","13","13","51","4e","81","8f","22","3a","47","fd","3d","6d","b4","6d","2e","1e","3a","a7",
"22","a4","4a","cd","4f","f1","f4","f6","ae","6b","11","04","71","23","77","ee","36","29","b7","bf","31","95","9b","6e","da","d0","97","8d","a0","15",
"51","8b","6e","09","07","8f","af","aa","dd","6f","2a","3c","58","95","57","2d","04","70","63","7c","17","f7","27","36","e2","f6","ed","b1","12","5d",
"75","e5","d4","5e","a9","0e","b9","09","be","c0","4c","34","1e","91","d5","cf","85","1d","85","4d","d9","98","cd","ee","82","8a","ec","75","b7","80",
"4f","50","c1","79","15","9d","db","36","95","f0","c8","fa","87","35","e2","77","2d","c9","dd","3b","d1","79","da","73","41","d6","77","8d","6b","b7",
"d3","a0","c6","27","c4","0b","ff","cc","7b","d7","25","f7","eb","a8","55","02","50","6f","62","61","6d","ff","ab","14","a4","a6","df","54","54","07",
"70","ee","da","b2","b1","84","3b","5e","3a","e4","c1","d5","69","3d","9c","87","8b","a2","9f","67","a0","fb","c9","04","1f","f2","54","59","d5","8b",
"12","65","57","8a","90","94","af","60","cb","36","8b","aa","ba","08","ed","b4","53","e8","0f","c0","c5","9a","6d","39","f2","ad","b4","62","f1","77",
"37","bf","88","80","6d","f7","59","a7","9c","b8","e0","0c","79","c1","e2","bc","94","6e","73","4a","28","f3","5e","9f","5d","1c","28","a2","04","65",
"a7","33","fd","d0","b0","68","32","e7","4a","39","22","82","c6","7a","e7","be","32","51","ab","c7","e0","16","5f","1f","9c","13","f5","48","8e","9c",
"18","0b","8b","00","c8","c8","80","f8","06","d5","5b","93","a6","57","83","a6","e0","ef","31","98","a5","4e","b2","5a","3e","0f","2f","56","13","fb",
"58","e5","df","e4","7d","f7","a5","93","cb","a9","e6","10","d4","85","c8","cb","1c","17","e6","8c","f1","99","a8","1b","35","e8","d8","77","99","e9",
"ac","46","a4","a0","61","65","56","73","41","93","88","18","de","f5","a5","24","6c","42","ef","68","3c","92","02","e4","ca","d3","cd","06","d9","33",
"0b","13","a2","69","f5","e1","3c","1d","ad","21","ff","b8","12","1a","b1","09","3b","28","7f","cf","5f","fe","9b","55","50","d0","43","5f","a7","bb",
"8d","58","95","33","d7","e2","0a","65","a2","bf","ee","90","1a","f9","9a","07","17","f9","7d","2a","c1","e7","ec","0d","2d","b0","36","8f","02","9c",
"eb","24","c9","c5","e4","4d","f8","48","ff","2c","44","78","87","b8","be","32","66","fb","a2","39","7a","a8","51","f2","8e","1e","da","3b","02","45",
"8f","d6","31","8d","51","7a","7e","c5","c4","0f","91","b5","c1","ff","54","81","4a","d8","b4","7d","90","f1","fc","24","e3","f0","bf","38","3a","46",
"a9","61","80","87","6a","60","df","71","06","7c","ee","01","ac","d2","9e","42","4d","55","32","c1","01","98","22","a1","62","1e","0a","a5","eb","42",
"b5","9f","7d","1c","0d","f7","73","d1","d4","da","80","72","05","27","b7","53","4a","dc","e2","41","56","ec","b1","72","b2","ab","8f","7e","41","88",
"0d","de","8d","8f","54","30","e7","05","f6","53","3d","ce","42","3b","9a","c3","42","3b","4d","a2","5b","f9","fd","13","66","1b","d9","83","e8","25",
"47","9b","36","64","15","1f","65","d4","0d","65","be","d0","13","4b","e4","29","47","64","f8","2f","34","ce","bd","ec","d7","5c","39","08","0b","1d",
"51","c4","32","5a","68","bd","c6","62","0c","4a","63","08","c5","8f","cd","96","3e","df","dc","7d","9f","3c","f3","7b","8d","8d","f9","af","99","e5",
"2d","c9","d2","eb","87","ca","40","03","f0","c4","08","19","71","dc","e1","57","3a","64","38","d0","13","4c","34","df","43","ae","6e","96","bf","15",
"4c","4f","e4","5d","cd","b1","93","c6","b1","cc","44","8c","44","8d","b1","c6","98","18","6f","67","41","72","4c","d1","ca","aa","7e","cc","8a","91",
"4e","cb","0d","36","3c","f2","04","68","d1","25","94","aa","22","99","04","be","95","78","54","a8","79","60","ac","4c","c7","96","f6","92","aa","36",
"84","2f","ee","20","2b","17","21","a8","05","25","df","6d","c1","6d","d6","74","fc","f9","7f","ed","a3","d0","67","00","2e","b3","11","93","4b","f8",
"aa","d0","6c","5d","95","ac","25","37","f4","77","bd","96","84","07","c9","3b","6f","46","dc","d3","17","a1","7b","0f","1f","23","f7","46","6a","9e",
"28","46","82","00","e9","76","3f","ed","c8","95","ff","96","b4","cf","4d","57","dd","3f","8f","f8","13","07","8a","dd","af","b7","df","45","fa","94",
"4f","f3","cb","9a","38","51","ac","30","58","47","51","49","0f","3a","dc","63","a3","6a","63","ab","e9","19","b1","d2","49","b4","00","8f","64","df",
"d6","06","ae","d3","59","f8","c3","aa","b1","54","1f","66","e1","87","1d","98","d8","e2","78","12","d8","0a","1f","5e","97","b3","35","5e","23","c1",
"76","74","f9","de","c0","4b","8a","b1","eb","12","9f","97","dc","b8","da","ef","49","86","96","05","bb","03","6d","10","92","06","a6","51","92","9c",
"4a","9c","66","b2","de","69","5e","7a","96","8c","5b","60","f7","45","f3","cc","95","6a","15","8e","4e","46","54","00","80","83","9f","40","cd","12",
"5b","93","87","54","12","76","09","d4","7f","92","25","c7","3a","0d","78","4a","8e","ca","16","10","85","16","c6","9c","22","dc","7a","2e","9a","1b",
"6d","ca","79","83","07","11","7e","14","3b","7c","41","10","96","b6","fb","ba","7f","82","39","ff","b2","40","2d","0c","70","a2","c8","87","60","77",
"90","93","7b","4f","23","2f","04","06","6c","0b","a9","07","3b","8b","80","c0","2a","66","e0","cd","d2","93","0f","a1","d9","c9","4d","a3","e8","75",
"df","f2","81","fa","83","f6","76","89","f9","32","fb","7c","df","98","2f","b6","c2","5c","31","7b","96","e3","f6","b0","bc","e5","33","45","7f","34",
"15","3a","94","6c","3f","9c","7a","ac","4e","86","ba","7d","a7","6b","e2","ec","99","66","31","d2","7f","8b","10","9e","67","bd","d9","a8","5d","22",
"e6","44","b4","fd","fb","09","d7","66","fa","49","d3","85","4c","78","d3","a2","87","76","cb","08","47","51","90","89","f2","a3","75","aa","74","ff",
"10","78","33","68","e1","91","a5","3c","01","c3","40","27","41","45","69","fa","7a","ee","bc","30","aa","87","5d","fc","ef","9e","e0","5e","72","b8",
"31","a3","d4","f8","f6","6f","16","63","75","47","03","9c","b9","9f","fd","c5","1a","8f","ed","38","6f","a5","af","cb","56","d9","da","a4","af","91",
"db","20","20","5b","88","2f","c6","f4","f5","3d","da","8c","c5","a9","35","54","fd","d1","9f","2e","83","ae","01","54","1d","35","31","22","33","d9",
"09","ae","7f","79","26","c1","40","46","a8","68","00","cf","da","b4","7f","6a","c3","08","5a","19","25","ee","ef","57","4f","4d","48","90","6b","c2",
"c5","a9","91","1c","c2","29","a7","cc","c8","85","66","98","5b","b8","75","7f","de","a5","35","d7","fa","c6","81","f8","3c","fd","92","6f","50","90",
"54","1c","b2","14","fb","c7","5f","a8","be","e9","da","a8","eb","99","f5","73","7c","aa","34","b2","b1","be","d0","7f","9b","d2","9d","09","cb","f6",
"00","bd","e9","81","18","92","79","1d","6e","06","00","18","65","31","97","73","48","7e","72","6d","df","34","c4","fb","21","91","55","75","7c","86",
"87","60","79","7a","de","b4","b7","08","d3","e4","77","b0","2c","f6","db","e2","8b","99","20","e4","e9","13","c7","fc","4a","0c","04","92","4c","b5",
"87","fc","6d","82","d6","68","72","32","6c","d7","36","23","0b","e6","d4","bf","ca","3e","90","0d","7c","6b","5e","6b","d3","c3","cf","9a","5a","20",
"c7","f4","f7","c3","78","9d","c7","7f","7b","b8","55","68","50","6e","94","e0","e9","93","20","58","63","ec","3c","53","28","7f","9e","9e","c3","b6",
"7d","78","2a","08","af","64","95","82","27","eb","d4","c1","1e","95","b1","f9","47","58","60","c7","6e","10","0a","41","13","5c","ee","ee","b5","3c",
"44","ad","ce","35","ac","a1","0f","ee","10","08","c4","bf","a7","e9","99","0c","d0","9f","2d","de","ea","3f","f5","a3","59","5f","c8","18","54","8d",
"60","b5","96","e5","c1","3f","6b","30","be","35","c6","43","04","50","55","f0","c6","1f","b9","84","42","12","f1","6f","a8","3a","3e","24","0d","26",
"8a","b3","67","60","4f","dc","a9","cd","da","f2","b6","1c","1a","67","a3","7c","c0","94","ac","c1","fc","e1","69","87","89","09","6b","a1","3d","0a",
"c3","5c","ac","11","12","db","93","ab","25","b3","a9","0f","4c","41","7f","12","02","c3","61","06","83","c1","be","1d","39","e0","85","c1","c3","6c",
"ed","47","cb","40","e4","94","78","d8","03","6a","03","0b","dd","06","c6","e3","f6","b8","53","18","99","f1","a8","22","3c","6c","73","a6","e7","22",
"99","ba","36","46","1c","6c","16","1c","92","80","53","a0","e7","9b","b5","f5","27","50","b6","ae","9c","53","12","00","14","ff","35","2d","81","c9",
"f0","ba","3f","4e","23","67","e6","e5","26","17","d7","80","ef","84","17","e4","39","78","61","6a","23","a0","30","2f","23","9f","67","2a","07","2f",
"28","55","f2","3b","ca","0a","c1","8b","54","f4","dd","52","cb","ba","48","68","4f","dd","a5","e7","60","87","be","a8","2a","4a","48","08","30","47",
"bf","db","81","68","eb","fb","0d","e4","6d","e5","c3","f1","19","9f","cf","fd","2a","c7","1f","54","c4","87","69","3d","4d","11","2f","b4","18","66",
"5d","c4","01","cd","6b","cf","96","b4","b8","b4","8d","84","0d","22","5e","7b","3e","37","de","3d","a6","4d","40","6c","27","cb","e3","f7","c8","49",
"dd","02","d9","c2","32","c6","7e","d7","54","58","05","02","4e","e7","f6","b5","d1","80","42","24","7d","71","63","8b","82","1e","5d","71","3e","a4",
"6b","09","86","99","57","27","da","87","5e","ef","71","da","a1","81","7f","31","3b","47","cd","90","1d","ed","db","29","fd","5b","f9","23","9d","01",
"f5","39","bf","fb","a5","49","8b","16","b1","f8","d5","a8","ba","5e","e3","c7","88","02","4d","8a","7a","e3","9e","e9","e0","30","60","96","00","00",
"3b","10","c4","4a","3e","92","06","16","51","99","f2","81","b9","f9","4e","97","f5","bf","4d","8e","a2","e5","77","8d","82","68","9f","39","3d","41",
"6d","40","7b","50","b9","07","74","b5","2b","4c","c6","4d","55","6d","28","50","87","05","e9","4d","82","59","09","ff","fe","64","96","18","75","ef",
"23","7c","ca","75","8d","4c","19","02","42","66","ba","74","f6","e2","9a","50","a7","b3","08","b5","4b","60","51","83","d5","21","00","94","af","bd",
"bb","4b","7a","da","60","e0","b5","5d","d6","d4","f1","c7","4a","1b","1b","b7","cb","f0","bd","f2","ef","bd","73","58","e3","aa","aa","ba","79","e2",
"9d","34","aa","73","e5","0a","e4","35","cc","d8","59","20","06","75","12","f5","67","5f","7b","4b","67","db","81","61","39","87","e0","13","71","b9",
"c9","42","50","70","96","ec","81","ee","fa","12","df","5f","4a","87","ae","28","35","43","ee","f9","0c","6e","c2","6e","77","25","22","36","cd","79",
"ab","aa","cc","1d","7b","40","a3","91","a3","c3","67","72","ab","5d","17","69","e7","ca","be","79","fd","9c","cf","a1","72","5c","e6","95","23","70",
"f0","55","02","2a","22","dd","2a","90","d7","f4","f6","03","12","95","f8","e0","95","f1","b0","8d","a6","60","1b","d8","fe","1b","3a","ca","41","55",
"29","2f","c8","8d","e9","55","72","72","8a","5f","f6","1d","bc","89","84","0e","de","04","d5","bb","db","f1","d3","b3","e2","fc","3c","d1","43","aa",
"d3","ab","e4","c2","45","ab","ea","c2","70","5f","36","00","00","14","8e","a7","72","98","d3","03","c4","28","7b","20","7b","33","5c","52","f7","f7",
"00","54","6f","da","ba","db","f1","cc","aa","8d","df","4a","3d","b7","43","b9","e4","8d","f4","12","e6","f7","11","76","7d","2c","c9","d0","8a","49",
"67","fe","34","26","76","c8","95","cb","29","50","95","59","15","29","f5","68","f8","de","c4","e9","1a","79","24","d2","25","a1","b6","89","0f","4d",
"b7","0c","63","67","81","9b","42","e1","40","75","4a","46","60","a5","4e","d7","4f","77","8a","35","15","c9","a6","b5","ef","06","61","4b","a9","8a",
"2c","29","c2","73","6b","eb","82","ed","aa","2a","10","e5","ff","60","94","7d","4f","b0","d0","bb","6e","9a","1e","65","eb","d5","58","05","d4","3e",
"9b","f0","34","55","22","b5","6e","61","62","d9","da","2c","b5","db","49","4c","70","16","29","b5","68","7a","12","71","5d","d5","7d","9b","7a","f9",
"4c","de","3e","20","83","b3","eb","39","93","85","69","a4","fe","ec","a6","51","d9","77","2f","5e","19","b8","0f","9d","3f","fe","ab","24","eb","eb",
"fa","5d","48","a0","bc","79","95","72","92","6c","36","f5","e9","49","25","de","07","63","09","2c","5b","a2","c7","db","3d","73","b8","b8","53","18",
"69","5c","7e","9f","99","3d","6d","19","cb","95","b4","8b","b5","0e","e6","82","50","13","28","95","da","ed","84","82","4c","cc","d7","bc","c1","17",
"33","22","71","04","eb","3b","53","dd","bc","57","38","7e","77","12","f6","a7","9a","33","9d","14","0b","f7","50","46","20","8d","75","55","4d","ac",
"17","0e","7b","96","05","72","1a","83","6c","79","b6","66","40","4d","40","02","79","17","8f","8f","60","45","a8","7e","b7","33","6a","10","59","ca",
"b5","97","ac","1e","ab","6c","e4","76","9d","ad","46","f5","3c","39","38","15","b4","6a","c2","b9","c8","8f","40","70","73","f2","6b","f3","3f","c2",
"5e","bf","5a","42","cd","b0","51","3b","a9","ad","ce","b0","14","97","5d","46","cd","af","31","a5","af","bf","99","5e","9c","2b","8b","33","67","15",
"e8","91","a7","81","f1","63","43","85","10","1a","6b","32","08","f4","f3","d8","4e","17","e6","fc","ce","5c","bf","1d","6b","4b","73","23","19","87",
"49","54","d7","58","15","56","20","74","a5","7d","cb","2f","7e","8f","96","84","7a","47","ba","9c","1d","2f","d1","87","2b","cb","e6","4d","47","9d",
"f1","9c","69","fc","12","1e","26","2c","9d","81","90","6b","4c","a8","8d","7d","d3","5f","83","62","db","fb","1f","92","aa","7d","38","c2","13","f0",
"c6","34","a9","31","6d","c5","76","8a","00","c5","f9","9c","55","ad","ba","da","51","f1","53","40","3e","6e","72","fd","04","56","ad","4f","2c","b5",
"98","9b","49","95","c1","55","e1","84","4b","96","f2","84","2c","fd","c7","f4","33","bc","cc","9e","41","1a","1c","58","90","07","7e","3b","6b","f8",
"dc","77","60","90","07","bb","16","5c","97","87","d3","30","62","15","d3","0b","a4","c9","1a","8e","9c","bd","fe","c5","6d","e1","97","e8","cb","ff",
"ef","84","c2","3b","b7","f6","d6","22","76","b3","64","3a","ba","70","11","94","63","5d","e2","ad","ca","9e","4d","31","8e","4a","27","88","5c","ff",
"98","74","a5","68","c9","45","16","87","31","63","f7","10","2e","01","90","dd","fc","a1","5c","8a","e2","65","67","56","2a","df","dc","24","e1","00",
"a1","a6","f0","8b","f9","2d","f7","b5","c1","f5","5b","16","06","07","4c","d0","f0","80","e9","c7","9f","2a","33","e1","26","ce","a1","53","f4","76",
"3c","0c","b7","9d","c9","9c","8d","4c","95","cf","85","45","a7","e8","bf","bc","d4","7e","46","dd","92","e4","77","0e","5b","e6","24","b4","86","b8",
"21","71","ac","c3","ed","58","98","73","f5","60","fd","14","5f","ee","30","66","f1","83","c8","38","f5","11","bf","b6","e8","9e","5c","c1","b4","38",
"cd","2b","4b","2b","29","fa","fd","55","a2","6e","8a","5d","6b","f5","66","5c","f9","3d","6c","c6","85","ef","9a","0f","1b","71","5c","34","38","a7",
"f2","5d","8f","02","14","9f","5a","84","1b","7a","9d","3d","31","ab","ed","a9","ba","2f","c5","2d","4d","60","19","96","bc","9b","1f","f2","df","13",
"a1","e2","8d","ed","84","68","89","f6","04","86","8c","a0","39","a6","aa","ce","d8","38","05","20","06","63","9d","75","dc","bd","05","15","16","a6",
"e0","2f","5e","9d","e6","f9","fc","99","c3","a0","98","07","12","6a","9c","7e","00","d4","ec","ca","80","87","02","78","d5","23","67","4c","68","7e",
"80","f6","54","32","58","dc","12","e6","fe","72","0e","92","e2","97","00","5a","8d","b5","b3","8b","07","e8","8c","ca","b5","bf","df","e2","6a","6e",
"49","81","14","fa","72","bf","41","45","6c","c1","7b","05","46","53","0e","f7","ee","1b","90","ac","25","7d","3d","ec","8d","fd","b2","4a","ec","07",
"a7","7c","52","b0","e0","0e","ea","48","19","b7","17","e1","86","96","00","05","7d","97","ea","c2","bd","6d","c5","fb","70","bc","df","3e","e9","9a",
"ae","ef","29","b7","d3","e2","78","c6","5d","1c","90","9f","f2","3b","22","e9","2c","a5","e7","51","5a","54","04","f7","dd","d6","46","31","23","55",
"b2","f7","44","86","36","4e","cd","57","12","3f","99","67","2d","ea","4d","54","8c","c4","8e","fe","d1","45","79","d3","d1","24","76","7e","dc","43",
"fd","f3","fe","54","6e","4a","64","3a","29","fa","af","76","bb","2b","61","71","3c","00","a8","ac","a6","c7","2f","cf","1e","16","11","9c","0f","7c",
"b1","7f","1c","c0","14","41","c3","99","6a","b7","b2","5a","35","ff","2f","28","04","23","4b","d2","a3","4f","d2","f5","a8","c8","fa","9b","fa","ef",
"96","ff","68","37","d2","01","71","78","96","8d","ef","a2","a4","e5","7c","2b","48","c8","8e","61","e3","91","32","a2","50","7a","83","cd","33","76",
"c1","a4","97","c9","eb","00","46","57","52","21","71","16","10","7f","26","48","51","8e","8f","9e","0c","18","77","49","44","ec","1c","dd","ad","9e",
"20","07","cd","aa","cb","56","bc","fb","29","1d","63","b0","91","49","4d","8f","76","15","9b","8e","e6","48","fa","84","65","46","36","75","e9","16",
"f2","08","3d","21","73","2e","0d","a0","9d","2f","b2","f8","a9","ed","49","0f","2d","9f","9d","10","22","7a","70","c1","eb","bb","b8","a5","b4","94",
"9b","4a","17","5b","20","9c","5c","a2","e1","42","67","ad","54","fd","5f","fb","e6","52","3a","fc","9f","3b","e4","36","b4","cd","0a","96","a2","93",
"96","61","7d","a0","b1","fd","88","8d","c8","36","f7","fc","a3","07","a1","93","61","05","61","96","d3","7f","33","88","87","d9","c0","e4","27","cc",
"14","31","5b","67","61","95","e8","25","20","a9","f6","c6","b1","00","44","a3","b8","a1","dd","81","19","80","e9","71","96","ec","bf","ef","8b","ad",
"1e","00","77","58","bc","7b","8d","ed","e0","10","82","0a","6d","45","52","d3","fa","7e","4d","ea","62","6b","09","fc","7d","da","b6","ab","b9","fb",
"23","a5","13","d0","98","84","ce","81","a8","23","44","a5","6c","8d","bd","0c","07","df","25","ca","8c","d4","7b","32","65","1d","7b","d8","48","eb",
"d8","43","ec","3f","21","85","40","8f","61","55","54","33","7e","18","e4","2b","7e","4f","f1","c9","25","41","2d","f6","bd","64","f5","fe","f5","39",
"bf","3e","06","5e","69","81","a5","96","1e","71","69","cd","de","0b","f8","5b","ae","c2","c6","e6","2b","76","3b","5b","96","ee","b5","76","7c","c7",
"cc","15","d2","8c","0e","52","a5","47","39","97","64","93","bd","b3","d7","84","72","12","04","a1","45","75","f5","e9","5b","a2","af","77","e9","7c",
"fc","94","b6","94","f0","ad","42","73","9d","f2","5e","cc","47","1e","38","26","74","d2","dc","a7","c2","a2","83","af","be","50","83","c1","3e","1a",
"42","79","d4","ea","bb","1a","4f","b6","fb","c8","6e","ed","60","0c","07","13","aa","8f","16","2b","f2","77","41","fc","d3","41","db","f9","f4","3b",
"bd","8c","e0","28","f4","04","84","25","84","ac","59","b9","4f","ca","c2","1a","09","cd","d7","a2","c8","24","52","6b","ae","ca","22","da","da","3b",
"ab","02","dd","39","6f","44","9e","d9","05","da","80","a1","b9","d1","b7","e8","c1","cf","bc","34","c3","b8","8b","25","07","b5","bc","13","48","3d",
"74","30","95","3a","22","47","79","29","da","2b","9b","99","cb","f3","fc","24","92","2d","f9","6e","5e","b2","e3","35","77","cf","01","68","ca","0d",
"21","f6","92","30","05","da","53","c8","3f","c7","16","0e","4e","d5","01","89","49","7b","83","8a","ad","48","59","2b","e6","2a","1d","7b","f6","26",
"96","3d","63","f0","63","e5","ff","05","52","34","b2","dd","ec","65","7e","96","55","06","1f","37","fe","b7","3e","7f","91","9b","04","04","65","db",
"56","3f","dc","2e","ff","d0","4b","3d","6a","d5","9c","3b","c7","77","a5","cc","08","14","1e","c3","b8","75","df","fb","60","75","f8","ed","07","94",
"27","01","68","0d","11","b2","b4","71","95","21","80","de","8f","32","a0","5a","dd","a0","90","5d","c8","2b","95","96","53","07","3d","09","db","fc",
"fa","58","0e","f8","b3","1a","88","60","82","a7","e8","81","c7","ca","57","2c","4f","25","25","85","14","bb","be","5f","a3","f8","11","36","60","9f",
"4c","b3","41","42","5d","bf","3c","db","d1","f8","13","01","83","ae","76","be","83","93","18","4a","40","ee","d9","e7","22","98","cc","2c","d9","cc",
"25","ae","c4","24","38","34","ca","6f","8a","03","6a","88","3a","aa","68","95","44","8c","35","81","3c","e3","cc","b3","ae","11","a4","da","64","a7",
"cf","80","41","9b","a6","be","54","78","ba","a7","cb","b0","2a","09","44","12","df","4b","32","d9","0c","26","53","62","36","df","9e","c8","98","eb",
"83","cb","da","62","01","f2","4d","fa","f1","1b","c9","fe","dc","b4","fe","d6","f6","2b","5e","a8","93","cc","29","55","ef","40","8a","b5","af","c5",
"de","9f","08","6a","bd","e8","50","c8","a9","b1","6e","61","bb","9f","4c","5d","30","00","84","72","42","25","6b","59","da","ff","9f","43","03","fc",
"60","22","27","ed","cd","6c","1d","38","4a","d1","f8","7f","c9","f8","33","57","99","72","93","a6","6e","42","4d","0e","43","79","00","c1","bf","49",
"c1","67","d0","2c","e0","40","b6","cc","a1","d7","4d","a9","6e","27","06","e7","05","81","a2","01","8f","9d","be","0e","7d","5f","f4","9a","e6","1f",
"e0","84","82","8c","e9","41","7a","a6","3f","2c","6d","10","66","dc","2f","5c","4e","26","90","b8","19","a5","ef","b7","b3","77","bb","94","0c","a7",
"a4","0a","2f","35","e9","ee","08","e6","e1","a2","27","0d","10","23","36","c0","11","5b","eb","1c","59","f6","4b","57","21","0e","8e","e2","d7","6e",
"d0","33","1f","1b","37","45","40","e8","71","7c","8e","eb","8f","7a","dd","78","32","18","e2","ce","c7","45","0d","11","32","fc","3c","17","4f","2c",
"0e","b0","a4","d5","d2","f7","7b","80","3e","95","6a","ca","8d","e0","e7","6a","65","30","7f","55","b9","9f","a0","c4","07","2e","19","10","34","a4",
"00","3e","3a","8f","0e","73","07","8e","98","a4","7b","27","4b","22","26","25","95","1b","f6","65","c1","41","b5","86","5e","98","01","35","12","ca",
"93","44","38","54","19","94","dc","1b","62","ba","06","70","3c","47","56","7a","7f","25","af","ae","bc","91","97","28","67","78","6f","d8","75","15",
"c6","73","8d","64","33","50","d7","f5","a5","14","aa","3d","d9","b3","9c","e6","6d","b7","6b","33","7b","9f","a5","6f","06","88","41","7d","00","ad",
"36","72","38","21","66","45","33","59","aa","49","96","d5","9d","eb","2a","ca","9d","b9","41","8d","78","00","09","a4","59","af","91","2f","63","a8",
"75","80","a6","40","a8","e2","81","74","63","d8","7e","6f","47","d9","7d","e8","3b","ee","86","73","1f","69","24","77","c5","00","eb","3d","7e","1e",
"17","69","e2","a8","b3","14","e8","85","a7","e2","a3","b0","b2","d6","9f","8d","4d","59","a4","7a","10","3a","3a","ed","2b","01","af","a1","ff","02",
"bb","a9","db","cc","b7","a9","22","3c","dc","b0","6a","bc","1e","89","1f","87","99","81","39","fd","a3","93","fe","32","14","12","ee","b9","9d","63",
"b0","5a","90","bc","18","22","3d","67","1c","e3","c1","73","69","68","33","10","ad","dd","19","5c","16","c9","9f","3b","1e","cd","57","0f","d1","43",
"3b","b2","2a","63","4e","ea","a5","0b","5e","07","00","76","88","a0","28","25","57","a1","3b","e2","06","ff","50","ff","77","c4","1b","16","17","2c",
"a3","89","d4","01","15","e1","c1","46","09","04","05","58","a9","ed","23","09","f8","f0","e4","c5","b5","53","f2","2f","ba","fa","5b","1b","78","89",
"3f","35","c3","b1","b8","99","92","17","a4","fa","1a","6b","1f","fa","93","78","72","62","23","c1","77","91","a9","04","2d","32","6c","85","b4","3c",
"22","04","d8","24","3d","8a","cc","bc","33","7c","82","34","c2","0b","47","f4","6d","bd","28","e5","31","35","f3","b6","c2","9f","ca","ff","ec","fd",
"3b","70","ee","b0","e9","8f","d0","22","42","db","36","fe","97","3b","65","8f","d4","be","05","1a","eb","71","32","76","a3","a5","e3","bf","2a","26",
"72","90","d5","f0","aa","d0","bb","eb","bd","a8","55","75","52","ae","b1","53","b1","b0","48","43","8b","09","59","8e","3d","ee","70","71","9f","10",
"36","0b","b8","ca","c8","62","1f","53","6d","98","fd","c3","ef","2a","bf","57","86","24","4d","16","16","67","ca","a2","d3","d0","65","49","ec","36",
"1b","7e","5c","fd","7e","03","91","ec","3d","2e","d3","97","81","6c","e1","0d","01","6c","7a","2a","71","d3","55","a1","7f","cd","98","88","0c","e1",
"f1","c9","7c","c5","78","fb","a2","13","20","32","c0","be","82","e7","ae","a8","c5","be","23","d1","1a","79","19","bb","4b","33","f7","f5","84","17",
"5c","64","6e","8d","b3","1a","ad","43","81","c8","a3","b7","2e","20","b2","d0","a5","9f","33","30","86","26","3e","cd","17","66","86","66","05","1d",
"e0","84","43","26","2c","89","4b","bd","a8","04","84","3b","61","a9","94","33","58","5c","6c","a8","d5","74","b6","2c","d4","0e","fe","bf","92","21",
"22","09","b8","be","d3","4f","19","d1","ee","f7","fc","c8","8a","ef","2f","e9","65","f2","80","0a","99","03","0a","f2","76","88","c5","79","a9","d0",
"2a","8d","f3","39","88","db","33","a3","51","80","9b","b3","24","ce","b4","d6","25","da","35","2e","0f","f8","2d","c1","c1","33","bc","0f","82","9c",
"14","22","bf","63","d9","e2","25","38","cc","a3","32","35","88","ad","ca","51","fb","fc","6e","9f","e4","19","60","76","cd","2e","79","10","9f","62",
"1a","71","0b","46","bd","d5","9e","22","13","e0","e4","7e","db","5f","6d","cf","9e","2c","5c","72","79","21","e0","d0","a3","86","cb","16","81","93",
"bc","d2","5d","e1","04","c1","b1","76","c6","db","50","b5","59","49","f5","c8","33","c9","5f","20","ab","4c","60","58","b2","09","a3","d2","0d","ef",
"8b","99","3c","0a","a9","fb","66","c8","73","ac","b9","5c","23","03","d9","10","23","6f","99","fe","cd","be","0a","08","f7","d7","55","8d","b6","cc",
"30","96","f5","ed","67","3d","d0","a6","b8","8f","85","5c","d8","4f","7e","d0","91","5d","cb","32","40","0f","1b","29","56","e2","a4","94","2a","d1",
"56","e9","f1","56","53","8d","90","28","dc","d9","42","a2","c2","40","a3","35","90","12","20","b8","86","f5","8d","83","ca","23","8f","8d","ed","7c",
"5b","8a","1d","3e","22","69","0d","e4","be","8a","30","cb","14","27","76","4d","94","7d","f1","e5","a9","e3","6b","1e","9d","9f","04","9b","fd","18",
"5c","66","a3","aa","a4","94","58","25","7d","d5","e2","bd","d8","48","ef","77","78","eb","70","92","44","91","e2","73","41","8c","b4","3c","6b","6b",
"f7","45","f2","20","cb","af","00","3e","2a","9f","55","05","90","18","01","13","2f","24","83","4c","e9","24","62","41","a0","55","f7","ee","8e","47",
"f9","3b","7d","01","27","4e","7f","f4","74","50","e6","e6","ca","5e","61","6b","1f","30","d4","42","18","17","3d","d4","75","97","37","ca","07","b4",
"29","d2","32","b4","2b","28","f3","5d","ca","27","1f","7b","43","3c","7c","63","a8","fa","13","30","20","c3","21","1b","75","02","60","76","1a","cb",
"b0","72","ee","59","17","00","db","a8","6d","90","20","ce","7c","30","f3","9e","36","84","b9","cc","7c","f3","37","cb","af","ba","c2","30","5e","b3",
"a8","56","2c","f0","a6","09","f3","fe","53","97","9d","1f","2a","d1","c0","1a","03","1b","d8","87","c7","48","d8","2e","b2","15","13","72","7d","09",
"c3","dd","03","83","a8","76","77","19","30","a5","ff","b0","74","ac","14","b0","df","84","4e","f4","98","fd","10","54","e7","13","a8","03","c4","24",
"ec","19","81","c1","fb","47","83","b5","05","da","f9","57","1b","a0","4c","06","22","1b","7a","cd","4e","08","b3","b0","95","97","39","54","9e","df",
"32","73","ef","ef","f8","9f","11","0a","6c","24","b8","cc","2b","35","ad","eb","4d","61","a4","1c","ae","e1","36","e8","ae","a3","8d","d0","f7","60",
"cb","94","b6","f5","44","29","46","8a","bb","02","c8","be","ff","87","c8","e4","c2","4c","e2","63","c0","86","7b","e1","8e","64","48","fd","92","42",
"32","6b","a6","d0","fe","74","e5","5f","76","1b","bd","a5","49","9f","b8","2c","68","73","c0","38","55","f6","f5","f3","6e","27","51","b8","c5","6d",
"d9","85","bd","bc","0b","d0","16","56","a9","95","27","f1","48","42","b5","69","d5","9c","20","b2","43","66","b8","4a","2f","60","a6","c1","c8","81",
"9a","3a","2b","74","e6","52","4b","60","0a","a4","05","ef","24","b0","ff","7c","cd","e7","45","d1","39","b5","3f","77","b6","c5","5e","37","26","d8",
"b4","0c","a8","c6","c7","eb","8f","49","3c","5b","84","bf","9a","29","2d","ed","e3","1b","31","82","4a","24","61","a9","8d","ea","e3","e0","4a","3a",
"56","4e","87","2a","14","5f","13","a3","76","4b","63","30","ed","72","47","93","c9","38","71","6a","3a","fb","72","59","b7","f8","4b","44","13","21",
"76","ab","e2","cf","5d","00","af","c3","f6","dc","55","61","70","4b","e8","dd","f0","fd","63","6e","5a","26","e6","d7","0f","13","c1","bd","86","cf",
"b7","4f","64","0f","89","23","31","46","7c","62","7b","ed","f1","5f","ff","84","f8","76","ad","14","f4","08","de","b5","9d","cb","37","b6","e1","5c",
"42","e8","92","49","3b","8a","36","9a","f8","26","48","fd","c9","b5","e9","b5","3b","c7","69","41","dc","9e","04","a9","0a","b5","eb","6c","d6","80",
"61","24","ac","a0","5e","7a","58","13","2e","4c","76","9a","9c","88","ff","1a","94","52","b6","f5","85","80","cb","8d","3f","ab","76","82","55","11",
"c2","06","2f","62","a4","13","e4","bb","93","52","0f","f1","e1","f8","83","2a","9b","f5","f8","b1","f7","b1","62","bf","eb","9d","92","ed","f4","60",
"90","ed","6f","7e","43","96","23","3c","aa","77","40","3f","19","39","24","bd","fe","13","a2","ff","91","6b","f5","b7","46","66","78","d7","ac","dc",
"a4","40","4f","c0","c1","3d","da","cb","7c","2f","03","61","6f","93","72","58","c9","d0","f2","0e","fd","9c","96","a0","cf","cb","b3","50","0e","90",
"78","37","74","89","57","ea","89","b3","bf","12","cf","76","b7","15","c1","e5","d7","19","2e","01","c4","02","51","7f","9c","13","fb","df","38","ef",
"97","86","36","35","ef","9b","fd","1d","5f","37","6d","d5","bf","cb","55","7c","e7","af","b0","53","c9","16","ef","94","bd","19","e9","07","88","9c",
"33","7a","30","fa","43","b8","07","76","e6","0d","e7","df","bc","9b","cf","ed","d3","06","22","34","c6","40","01","1a","f5","e0","17","0b","dc","72",
"04","53","ad","c9","ab","59","db","b9","6c","d6","24","35","5f","38","cd","3a","35","27","c2","9e","a0","7d","52","aa","7a","c2","ee","be","70","a6",
"f4","f1","9e","2d","1f","3d","a9","38","2c","57","d2","9e","d3","7b","ac","8c","f5","60","a5","b8","51","46","e0","7b","88","e1","02","02","f5","ab",
"c0","2b","73","a3","70","58","54","1a","e6","c6","31","89","04","3a","71","2a","12","78","9e","40","c5","d5","66","71","ee","58","23","4e","f3","44",
"d9","c5","55","73","51","28","75","d6","50","a7","d6","91","42","da","5f","ea","08","37","3a","82","8a","f0","c2","58","ed","a8","48","30","a0","86",
"73","b8","f8","33","fa","f4","05","cb","29","2f","a7","83","e8","6c","c5","3a","ad","88","34","61","4b","23","b7","1e","ce","a8","4f","fb","d7","f4",
"e2","ab","47","99","62","14","17","1a","2b","b7","05","bb","5a","ab","ff","09","2f","44","06","f9","ef","10","48","86","fb","7b","52","6a","61","97",
"c5","d4","40","d3","b9","4f","b4","37","f1","6f","58","8d","dc","04","60","2e","5b","b1","7c","24","3d","28","21","ff","d4","9c","ad","4a","4c","53",
"39","5a","70","87","5f","f6","c4","5f","0f","8a","90","8f","d3","87","0c","b4","0b","cc","51","6e","d7","36","a0","0d","cc","b4","aa","ef","f4","d1",
"d6","60","e1","07","13","b5","97","33","7e","46","78","0a","d8","cb","f8","c6","00","3f","b5","34","b8","14","17","e7","b8","b1","d6","a4","37","cc",
"fe","44","74","38","aa","17","f4","dc","46","28","cd","23","72","a0","06","3a","fc","ea","c1","bc","45","ff","08","8d","6c","22","83","46","a5","9a",
"35","51","b2","0f","5d","ba","6d","59","fb","e2","0c","6a","db","20","fa","f3","b6","15","91","32","bc","38","d3","87","ae","37","d2","19","a8","0d",
"6f","67","07","e4","e5","81","a1","56","40","59","d5","82","87","d2","b1","77","7b","c1","86","8f","a3","04","25","f2","55","10","6c","7e","ff","68",
"aa","22","21","0f","e7","bd","fe","6c","f8","d7","e1","2e","a3","18","e3","e6","44","b0","75","ee","9f","3f","04","06","52","84","bb","97","3d","70",
"70","33","1d","45","6e","ac","61","a9","11","1a","f5","90","70","71","cd","eb","74","03","84","f0","4b","bc","2a","b3","59","0b","cf","fd","40","ed",
"82","af","b8","d2","f6","cf","d5","48","45","27","77","b4","55","df","f6","e4","3c","b9","36","63","62","35","33","50","cc","43","00","7f","79","8a",
"6a","b2","87","59","14","5e","54","2b","a0","4c","ea","ae","4e","57","d8","e7","ae","b8","46","42","3f","a7","4d","69","92","b0","4d","35","03","11",
"f2","8d","1d","d8","16","9a","5c","a0","b9","6f","ee","48","44","f1","71","b2","7b","83","c3","d0","35","a2","0b","9b","ec","25","f4","51","90","08",
"3f","d1","e3","4a","f5","cf","95","1e","1e","fe","41","66","d2","87","30","2e","db","75","6d","68","d2","3c","85","7c","c4","cb","30","2c","fd","64",
"e5","95","c7","8e","6c","12","84","3d","50","63","0d","9a","e9","1d","b8","7d","64","14","9f","07","01","a8","d6","4f","d4","2a","48","5a","31","7a",
"c1","9c","7d","47","0e","a5","1d","90","aa","2b","db","b3","b9","1b","98","b1","fe","5a","4d","bf","0c","17","5e","04","de","46","c0","d1","5f","28",
"66","91","69","3a","02","ba","62","fd","2d","f7","ad","91","44","be","da","44","17","4c","95","06","48","ef","19","77","96","45","36","4d","a2","6d",
"a7","4b","8f","71","b9","12","68","b8","6e","88","9b","9d","25","f6","45","06","8f","09","af","d4","eb","bb","2a","64","1f","5b","3f","4e","6e","fe",
"0f","b7","e8","e7","a4","96","4e","d7","cb","98","87","11","03","e6","b0","b0","c0","67","b5","8a","37","86","a0","e6","aa","5c","06","70","14","17",
"db","ae","88","3a","54","21","58","7e","36","c1","47","92","8b","2d","dc","1a","88","0e","07","f7","c6","e7","e4","81","ed","32","85","4a","4c","35",
"47","67","ec","0a","f0","eb","e6","a0","e2","47","4d","79","05","51","ff","c7","0e","eb","61","32","4f","34","5d","2f","da","16","87","8b","e7","ba",
"b0","10","26","d4","db","4c","9f","b7","e6","16","ea","ed","01","be","37","63","29","84","83","bb","b4","9d","97","21","d7","e8","65","d8","b7","9a",
"13","70","b8","93","79","76","3b","cc","23","4a","11","d5","15","b1","17","6e","89","f4","d6","1a","d1","63","57","05","cb","60","f1","25","9b","f8",
"9a","8f","21","16","41","a7","1b","18","82","84","ef","bc","a6","22","d4","b1","66","01","65","31","56","62","da","08","fb","7d","c1","ba","23","09",
"17","b7","df","b3","6b","48","09","5a","d6","d7","63","39","52","09","ad","c8","62","ef","f6","7c","67","c3","1a","26","63","3e","18","14","90","ef",
"72","84","3d","29","18","b4","64","66","00","28","48","d8","82","27","ef","a5","b5","0e","8d","48","22","59","0f","86","ba","c5","f9","3e","7d","39",
"76","31","48","27","6a","93","3a","9f","95","21","2d","03","55","33","bc","aa","11","79","38","1d","5f","07","76","f5","5a","00","f6","57","42","5a",
"66","49","d5","3e","cb","55","a0","64","73","25","70","14","26","e2","00","aa","ee","8f","23","c4","81","7d","f2","1d","a7","6d","bf","74","fb","0e",
"3b","58","d1","61","fa","b6","88","34","87","3d","74","4c","b3","2e","67","6d","d3","55","5e","0c","93","53","80","96","4d","99","f1","e3","22","ca",
"74","e1","9f","bd","13","fc","87","95","66","5d","76","d6","1b","46","39","4e","fb","9c","49","2a","b6","0b","89","c9","40","53","33","2a","24","11",
"63","a3","3d","cd","4f","de","74","a7","38","7d","a5","8c","f8","e8","7a","3c","d4","2f","98","f9","41","ea","12","94","10","ce","cd","16","f9","64",
"69","bc","15","ff","75","1a","27","c6","01","e4","f1","74","a6","8b","07","fd","93","95","49","e4","55","6d","f9","85","2c","9b","b0","bc","55","d1",
"0b","01","72","cc","c0","76","e2","86","bb","85","c6","46","d6","35","57","65","e0","b2","94","80","86","2e","be","6e","75","4d","da","3e","eb","46",
"ba","c8","30","37","99","c8","81","78","30","b5","f3","66","3e","67","47","f9","e9","61","df","e8","fa","45","07","64","6a","68","e8","cc","92","4f",
"92","f4","a0","da","62","9f","33","d5","5b","e1","d9","42","8c","b5","75","95","cc","5f","83","69","d6","db","2c","30","12","d6","d3","d9","29","51",
"36","1e","7e","a3","61","d3","55","e7","d2","7c","dc","53","02","89","4c","ba","bd","ae","4d","a6","d3","f2","bf","18","66","84","2a","ac","db","50",
"42","2e","2d","57","b9","2e","aa","5b","55","8a","62","88","01","e0","4e","94","11","44","ec","fa","fb","0c","dc","44","19","8d","be","4d","7c","b9",
"6f","57","f8","28","35","35","1a","e7","a6","cd","25","5d","cd","45","4a","0a","71","7d","a1","22","ad","62","cf","29","fc","06","65","40","27","51",
"a8","53","b1","a2","b8","c2","42","92","7a","55","8a","19","19","6f","24","c6","30","cf","69","e8","0a","de","40","f5","01","47","54","9d","a0","8d",
"c5","b7","9b","7d","e8","3f","00","fe","eb","7e","51","56","84","56","8e","cb","e7","8e","74","82","58","66","cd","81","6d","df","54","4f","86","12",
"90","04","03","d4","e3","62","19","99","6b","5a","33","2c","5a","a3","a5","31","55","70","06","27","03","2a","36","54","9b","34","53","f1","d2","5a",
"be","99","1f","9a","fa","7c","5f","87","29","da","27","6b","3e","02","bd","22","94","83","27","c0","a6","28","65","85","80","64","a0","16","7a","db",
"a1","06","f0","68","bf","99","cb","cd","a3","e2","69","c0","18","31","9d","60","a6","c0","a7","be","32","69","e6","6f","69","2c","ea","64","fa","d6",
"1b","43","61","4e","ca","95","94","b6","dc","92","d0","ff","d4","2a","f2","49","b3","4b","c4","bf","d5","f6","5b","bc","ae","54","9b","40","b6","66",
"a3","96","7c","3a","55","a9","c7","b1","ff","90","a9","1b","e8","da","7e","0a","ba","48","c7","ef","ab","8f","f3","8d","89","9b","3f","59","49","f0",
"84","44","c5","ec","94","e4","89","87","7f","b5","d6","38","ff","8f","87","41","f1","cf","c0","6c","f1","12","20","55","d5","cc","38","25","a9","c3",
"7a","33","d8","6c","48","44","1d","3b","ca","4c","65","5c","de","8d","41","44","66","c4","1a","6b","81","11","b1","a0","01","59","1f","b1","0b","0c",
"76","81","f6","79","5f","3c","33","9d","ce","d5","34","ed","a2","75","8b","4c","84","a5","7e","1a","17","40","13","08","27","6f","32","ee","e7","b7",
"3c","34","06","49","cf","2d","d5","f7","ae","96","1a","e2","64","9f","d3","a1","e6","66","64","a7","33","78","e9","9c","f1","02","b9","41","1f","88",
"92","f6","e5","15","c9","19","09","42","5b","a6","5b","c7","04","fb","af","4e","fa","f5","70","48","cc","6b","8b","ec","13","7d","ee","97","85","57",
"e9","16","fd","98","78","d7","9b","0d","44","ac","7b","2c","6f","f3","b7","99","44","c6","cd","16","b1","3d","fe","c9","f7","3e","3c","ec","5e","1c",
"3f","8c","89","06","f4","de","6d","bb","e6","a2","6b","ff","4a","c8","84","53","a5","04","4e","24","f9","47","ea","de","75","9d","66","25","a3","f0",
"6d","ae","c3","b1","e5","86","3d","66","0b","d3","28","6b","ba","30","01","f7","6d","7b","23","7f","6f","71","f7","b4","cc","be","db","ce","b4","94",
"45","5e","a2","44","28","58","95","dd","78","d6","c2","a6","09","b1","2b","85","8b","73","a4","cc","e0","68","12","03","2d","a1","40","79","37","79",
"b5","08","cd","66","e5","3c","df","56","0e","1b","4a","39","3f","d1","0a","1b","67","67","f8","12","73","36","4b","85","10","c3","81","41","0a","52",
"e8","5c","86","1f","29","82","49","45","56","6c","d0","42","52","97","db","05","36","b2","30","93","76","06","43","c6","39","95","0e","41","14","a5",
"c1","1c","45","e8","97","ad","ba","cc","12","71","b6","96","88","d2","df","e8","af","42","4c","2f","24","d4","8b","26","cf","02","b4","86","aa","45",
"a0","dd","80","6f","5b","eb","a5","67","b3","30","68","a6","e0","b4","2a","fc","a8","13","61","2b","88","cc","0a","29","ec","14","3a","3a","48","08",
"95","96","b1","f5","88","56","9f","d4","df","23","ea","79","e4","bb","31","bd","90","31","e0","74","97","38","34","ad","11","00","aa","f8","2d","dd",
"cf","dd","81","c5","93","cf","0f","6b","a0","f7","a9","e2","17","8c","f1","c3","6e","cf","d7","c1","60","a8","9a","84","fe","bc","9f","c8","1a","25",
"52","48","37","53","cd","13","86","0e","3d","7e","50","a9","bc","1c","ec","46","49","f0","ae","11","5f","75","08","b7","50","a9","1b","ea","09","b2",
"aa","4f","4c","b2","5d","e0","26","16","4b","eb","59","7d","70","80","e0","0f","94","a1","81","c9","f2","ae","b1","39","d0","3a","be","cf","bb","5d",
"94","8e","37","c5","85","e1","62","24","2f","28","eb","3b","81","16","3d","1a","34","5b","3f","1c","68","af","9b","6e","a6","6d","c0","6f","8a","ba",
"3e","db","7c","79","b1","dc","9d","a2","56","61","9a","94","ae","f6","c8","40","2d","59","2d","b8","d2","af","7e","9d","e9","bc","ab","31","3f","df",
"46","10","5a","e1","5b","63","20","a9","f8","5e","5f","2f","2e","c6","56","8b","d6","4e","14","3a","61","1e","e0","d2","d4","13","b3","47","4a","71",
"39","51","df","b2","d9","50","7f","fe","e9","9f","44","6d","e5","1a","08","ba","67","15","b1","15","da","98","5c","a6","23","e3","e3","6e","b8","e3",
"56","21","be","21","93","02","2a","05","4d","8e","80","09","13","7c","b3","6b","0e","47","20","d3","0a","04","7d","4e","3a","41","ff","e1","42","68",
"31","e3","b9","2b","4b","6c","fd","4e","bd","6d","e3","75","91","8d","79","2e","da","92","5d","ed","55","49","cf","9b","21","a2","e9","7c","c1","18",
"ac","ca","37","66","c7","3e","b6","be","6a","3e","aa","92","86","66","b6","61","6e","4a","f0","74","bf","4f","ee","e5","2b","3a","7b","e9","f6","0d",
"c5","82","6a","e2","12","24","22","70","00","a7","c0","e8","1d","ec","54","b9","eb","7a","23","ba","e2","65","be","fe","ce","76","3a","dd","4a","d1",
"1e","17","50","5f","b9","c4","54","ea","38","5c","47","0d","76","cd","2e","15","6e","6b","60","04","7a","e8","84","18","b2","23","c0","61","1a","17",
"e1","72","50","72","a0","be","80","a5","04","63","66","b3","9c","c0","a7","5f","fd","9b","4f","45","91","88","09","ce","32","10","df","02","6e","db",
"55","9f","9b","29","d1","b1","6b","a5","69","70","78","fc","7a","26","c4","88","80","2e","58","8b","6c","76","a1","60","9e","24","14","fe","e2","f4",
"7e","ff","44","bb","60","9d","3e","12","07","e1","5f","2b","c5","4c","f0","8d","ca","26","0b","53","6f","d8","35","b6","36","8f","5e","39","be","cc",
"b2","cc","ed","48","b0","de","43","43","19","34","b8","70","49","54","fc","7d","65","57","43","1e","ed","5f","0c","82","85","46","6c","80","d9","d6",
"70","27","da","68","12","c9","c2","10","75","5e","e5","96","23","b4","41","03","45","d9","90","11","07","8f","42","62","3f","4c","0a","4f","64","d2",
"45","6b","ed","ac","4e","98","ae","0e","5b","cb","88","79","21","ee","f4","e3","4f","27","56","30","46","17","9f","c6","ca","76","13","e2","dd","03",
"51","d3","22","44","bc","cf","f6","c5","78","b9","a3","5e","85","8f","06","7b","e5","6b","7d","8c","9e","73","74","f6","f5","3f","4c","2a","d7","83",
"b9","f2","00","80","ac","d7","56","af","73","4b","b7","d7","8a","09","86","7b","1b","8c","63","ba","00","b6","08","0a","c2","ec","e9","3e","66","71",
"02","42","84","2a","94","2c","e2","f6","27","f3","10","b6","76","93","f4","32","10","e3","ec","7b","7b","19","cc","a0","1b","8a","d6","50","1e","ad",
"eb","5c","16","a5","fa","f0","87","dc","26","01","d0","2e","c8","6e","f6","5c","67","35","1b","36","3a","3f","1e","3b","66","66","30","b8","d7","18",
"7b","30","c0","af","02","b7","66","f7","6c","70","fa","9c","e6","0f","99","dd","ee","e5","00","75","51","48","01","40","8f","16","15","5b","8a","52",
"38","3f","ba","4a","1f","06","f6","85","80","88","0f","56","f6","9b","1c","2b","9a","72","51","bd","bd","d6","ae","b2","8c","a4","7c","05","cc","e2",
"5e","14","aa","03","53","f4","ab","4b","cf","41","5a","cf","e4","39","79","21","75","0b","14","2b","a1","61","9e","e6","02","4f","66","b8","88","7a",
"56","fe","15","3b","54","e3","4a","2b","60","db","b9","24","61","a7","56","06","c2","ca","73","0e","b5","ac","da","5e","fe","94","51","d6","26","f2",
"c1","cf","8f","1c","9a","86","fa","cb","25","b1","ec","f9","33","dc","be","b8","cc","96","c4","1a","3f","97","e4","74","b2","b0","07","42","d6","4f",
"87","6b","b0","36","ac","d3","29","53","13","a6","92","56","86","82","e2","c5","aa","67","5f","63","c4","52","d2","9e","34","12","e6","a9","d2","31",
"7f","93","24","57","87","51","a9","a9","8b","62","20","61","63","15","db","54","ff","64","75","22","b3","23","1c","21","02","05","10","80","d4","77",
"9c","35","67","5c","c7","59","0e","b5","9e","ff","27","55","a3","ae","9c","2d","e2","e8","49","be","01","57","44","58","30","2c","a6","0c","d4","82",
"d0","31","a7","5b","1a","84","ca","91","db","24","a8","0b","60","40","1c","93","6b","65","1e","2c","99","04","6d","be","e6","5a","94","5f","bb","4d",
"72","c4","51","5c","da","e8","dc","ab","5e","76","dc","89","93","c0","73","17","f5","10","29","9a","28","89","c4","5a","4c","c1","25","e4","36","d5",
"b7","fd","90","20","07","da","f0","43","e1","ed","f2","2a","31","a7","df","f9","a9","51","b2","25","06","3e","d4","38","e9","53","ed","f0","a4","a9",
"f1","8e","da","8c","eb","c8","00","c5","f7","f8","c9","7f","91","31","52","a8","b9","87","c4","eb","56","e6","4f","0b","75","02","57","a9","b0","e2",
"15","de","22","f4","45","3b","4d","11","e6","34","5c","29","4c","4c","2b","b4","57","a4","27","8c","25","f9","1f","4f","93","08","21","00","30","a2",
"6c","74","d2","85","9d","f3","6d","4b","f7","db","5e","b0","a1","08","9c","70","11","7f","1a","e2","26","47","27","fc","a8","ed","03","16","24","37",
"f9","18","c4","8b","a4","b6","b5","5d","bd","d7","72","1c","29","b6","da","25","27","92","89","e5","a0","5e","12","82","80","cf","e4","70","73","97",
"32","78","f0","df","f1","72","d8","34","f9","1d","74","6b","9b","56","a3","24","90","d6","05","84","61","61","a5","06","53","0a","78","4e","3d","90",
"d3","d6","57","83","62","67","5c","94","c8","db","03","71","f6","15","40","59","8e","61","75","a1","a3","b1","14","db","13","d9","6d","42","35","df",
"20","f6","f4","e1","67","b1","70","b1","5d","4e","18","a1","c4","83","2d","02","bf","cf","d4","9e","75","30","17","1b","9e","05","77","fe","68","89",
"e2","2d","d5","d4","cc","51","c1","53","94","4b","09","c0","90","fc","87","59","72","df","03","f8","62","f2","d6","32","a3","16","77","76","1c","20",
"38","a2","a3","78","76","a2","aa","0c","c0","21","37","4c","0d","c0","ee","1d","52","53","e3","ab","9c","54","5b","5b","be","34","fb","ca","88","7b",
"8b","9b","72","71","75","5d","d4","9c","70","79","38","6a","8c","42","a8","b8","fe","a6","b7","06","ac","9d","af","14","46","cf","39","d0","ae","f3",
"21","c1","25","8b","80","7b","fb","26","e4","44","52","6b","42","64","72","95","30","06","65","c7","5b","68","86","16","74","8a","93","fb","0a","ba",
"91","0f","70","91","80","61","b4","39","79","ad","c9","0d","24","3a","5e","2f","7e","bc","db","ae","ce","49","ce","a5","0f","ab","72","39","b3","f9",
"03","68","95","47","19","5f","7b","d5","2e","f6","9a","60","d4","63","32","b9","47","a8","e9","e7","cc","d1","ba","67","f4","9b","4a","53","3d","ee",
"e8","86","13","4a","8b","ed","49","b5","06","22","2e","4e","4a","1b","b5","ff","4c","3d","b6","57","34","24","90","cb","a5","d3","15","d6","5b","82",
"c7","42","fb","2b","4e","34","a0","7a","02","cc","93","03","6f","7b","53","3b","3f","c9","2e","31","87","13","05","eb","31","9f","27","97","e6","bc",
"e9","07","25","33","6f","ae","ad","cc","57","62","f3","aa","aa","db","28","38","01","40","ae","c4","0b","6a","ee","11","ca","75","0b","8b","c3","be",
"7d","99","4c","81","1b","b6","65","ad","38","0b","96","2f","c7","9f","1e","1c","d4","d5","80","a8","73","6f","93","d2","20","ab","32","61","fa","2f",
"15","58","13","38","b9","56","c9","e6","11","75","55","04","51","21","5e","7b","1b","42","a1","4e","64","61","af","a7","ba","91","d0","3f","09","0a",
"06","51","99","a5","98","a3","90","c5","ff","36","1e","6d","c8","6f","75","b1","56","9d","49","cc","8d","74","9c","fa","4a","9f","e2","70","67","50",
"8e","ca","0f","16","bb","a5","e1","d3","d3","d3","61","5e","a2","e4","59","46","5f","69","b4","0d","4f","78","2c","f8","75","9a","b3","3b","ed","dd",
"5f","1a","65","55","aa","43","ee","bb","0a","68","9e","46","c6","93","f3","d1","46","a4","85","97","24","f0","67","2e","7e","22","f4","e0","d4","62",
"72","eb","ed","54","93","4f","7b","57","16","73","4d","dc","84","49","da","b4","a8","51","71","44","34","6f","ab","db","04","c9","19","a6","02","0c",
"cf","df","59","af","0c","8e","74","82","6f","9c","c6","b2","b8","5b","8c","7d","3b","f8","69","d0","7f","16","28","24","ef","66","9a","02","45","2c",
"4c","69","71","1d","b2","45","21","b8","55","88","ff","b0","a2","e1","6a","84","cd","e4","6c","3c","af","65","be","87","1d","6e","62","1b","67","47",
"d7","79","11","62","d5","94","d0","6f","8d","56","33","3b","cf","20","f7","e4","3a","84","d0","b3","6f","c2","c0","77","6c","c1","17","5a","6d","5c",
"2e","6e","b5","d5","d4","3f","23","4f","5c","6e","6b","7b","bf","91","52","76","7a","b6","48","d5","de","db","98","87","43","a3","56","c7","35","74",
"ca","c7","f8","b2","10","bc","61","19","25","67","d1","f0","fb","36","00","92","79","63","f4","11","89","e4","dd","95","f7","b4","d3","48","aa","45",
"fa","f5","7b","66","86","b4","5a","2f","ee","b5","f6","08","56","c8","c2","2f","38","37","74","f8","90","89","46","aa","c4","46","d8","e7","78","e7",
"d6","bc","eb","5b","52","00","a9","42","32","13","c1","fd","fd","a9","db","90","09","a0","f9","23","29","1b","93","13","c2","b9","3f","a6","54","a1",
"5d","a1","b5","fb","1c","44","56","44","62","13","3a","88","28","bb","de","8e","6a","e7","a4","fd","79","ea","62","3e","59","6c","5a","f4","dc","70",
"43","7b","00","56","1f","d4","98","3e","cb","52","4e","78","28","e2","a0","11","4c","15","ca","3a","93","bd","f5","b5","2b","12","28","40","de","8f",
"20","fb","8f","c6","81","c3","33","47","fd","51","00","af","1f","81","bc","44","50","9e","dc","e3","c6","24","54","42","2d","fa","86","bd","bb","57",
"f8","6c","ed","6b","24","2c","30","8d","c4","ba","c6","a0","34","32","dc","f7","c9","76","b6","96","3f","50","32","ac","82","83","35","0a","d3","22",
"82","15","28","7e","59","dd","14","6b","24","8f","0f","dd","da","2a","72","aa","0b","8e","1d","3a","fd","f0","dd","14","d4","3d","80","01","f8","35",
"5a","fa","74","b2","67","8c","67","c8","82","58","4e","8f","72","74","bb","8d","04","29","04","8d","04","4e","37","16","f9","04","c4","9d","83","19",
"8b","73","fe","b1","7f","bb","14","28","30","3d","55","75","1e","53","c2","d3","87","ba","bf","1e","f2","90","6b","5f","bb","8b","e3","03","b2","6c",
"01","6c","4f","59","48","8a","56","b9","a0","2a","cc","f6","41","ec","13","be","ee","2c","e1","d9","7d","d6","17","62","b7","b1","7c","82","39","96",
"8d","06","e6","1a","79","51","fb","1b","d2","2e","e2","a5","e8","2e","f6","aa","d1","ba","41","75","33","2b","2a","74","d4","c5","d1","a9","df","52",
"d6","d8","c4","a7","da","30","93","62","79","88","53","2e","b1","a7","d3","3b","a1","8c","5b","bd","38","94","7e","d6","1a","c5","68","e4","14","b4",
"f8","fb","72","2d","08","f5","ec","9f","a4","a2","e2","d7","1e","e7","1e","a3","a2","af","ed","7d","72","67","3b","31","e3","8b","48","74","27","3d",
"fc","3d","da","f5","c0","a2","66","b4","24","dd","6e","77","7e","08","b3","4e","b1","2b","a1","f4","da","06","f6","8f","37","da","eb","18","1d","ad",
"0f","23","96","c8","66","f4","c2","0f","95","2f","a2","01","e0","82","a9","6e","60","6f","d4","f0","92","90","39","a3","7d","2b","64","89","55","c2",
"2c","7e","ca","ff","d9","a3","15","ec","c4","43","cc","f2","c7","c1","4a","62","70","81","77","c6","ac","db","e6","cf","32","e7","66","c2","f2","80",
"3a","6c","58","02","81","62","bc","2d","6c","5f","5f","9b","ec","67","14","fa","49","c0","89","09","b1","93","09","80","0a","c5","2c","ce","71","f4",
"26","9f","e9","2c","2c","9c","d7","54","49","06","fd","36","c9","c2","20","93","1c","22","a3","72","c2","60","19","56","e8","31","ea","d2","21","42",
"b9","a1","9c","11","a2","39","f3","6f","5e","c8","b6","91","30","91","19","5f","d5","fa","4d","45","50","ad","fd","a7","ef","a7","04","53","cc","cb",
"d3","e7","44","6a","6d","5d","ac","ce","5b","5b","95","cb","b6","6f","d2","38","ef","77","7e","57","75","0a","66","38","6c","35","6e","61","ca","77",
"2a","e1","41","b5","07","0f","a0","2b","8d","03","7b","76","c0","69","51","d9","62","f5","d8","f4","e5","6e","dd","0a","53","94","a1","df","70","ce",
"d4","40","30","da","a7","da","fb","c7","26","ee","94","80","af","ac","88","a1","ff","58","19","24","65","13","a6","c2","21","98","2a","24","25","ef",
"1a","c2","5e","c8","a6","bd","53","bd","12","fe","e7","75","f8","e3","eb","32","c0","14","6f","61","ca","c6","bf","a0","f7","c3","fb","3c","d5","63",
"b9","74","fd","ce","a5","d5","35","12","a3","3d","b3","24","7d","50","46","f2","c2","1d","de","8e","ea","35","f9","ea","02","2f","c0","6f","34","6b",
"5e","3a","ec","c7","57","b6","7b","48","10","0e","29","bd","dd","42","5b","92","2e","50","07","4c","2f","f7","b7","a9","98","a3","9b","1b","ad","23",
"a1","13","ec","0a","64","fb","53","7f","c5","c4","90","4d","c8","e1","30","35","c1","91","97","c0","b0","4b","ca","d1","5a","65","60","45","c4","ee",
"b6","b9","15","90","6f","3b","f7","48","5c","74","79","dc","7b","51","56","d7","d4","0f","43","48","5a","91","fc","72","99","fe","31","a1","32","f6",
"38","08","5b","ab","e0","75","91","fd","c5","b9","11","cb","41","ad","7c","5c","cf","69","a1","31","7c","74","cb","ce","d9","75","ed","0c","cb","04",
"0f","47","eb","91","28","11","ed","90","e7","fb","f5","d7","62","27","a0","25","d1","88","1a","47","7e","9e","69","7f","42","ff","74","3a","ca","10",
"9c","15","05","60","01","8d","02","91","e8","d7","08","9e","3f","ea","d5","48","04","70","dd","63","b5","ad","3e","4a","18","66","33","00","66","ed",
"c9","17","9e","fe","f5","1f","cd","d8","7c","f5","70","44","21","03","b5","c3","21","1e","fb","4d","d8","c0","ab","e8","f3","1a","0e","9b","d1","c2",
"b0","5a","8d","d4","4d","55","92","cd","0d","2f","5d","91","d0","d8","d9","0a","1f","30","88","6c","f4","81","74","0c","80","4d","f1","35","8d","26",
"8b","86","63","0b","a0","1f","c1","8f","c3","f2","3a","b2","a7","0c","ad","69","7a","86","43","a8","c8","a5","94","72","ac","cf","66","66","f4","e0",
"c4","95","02","d7","06","e2","68","4d","19","60","71","12","47","65","98","57","74","7f","39","ca","75","62","38","2a","ac","37","e7","ae","b5","83",
"d9","f9","a2","59","0d","31","14","15","73","3b","61","06","4f","23","ab","87","07","8a","c1","71","84","ca","d2","16","87","1d","5f","03","b9","f6",
"d9","a8","64","4c","fc","ca","ba","fc","e4","d9","24","ef","f8","6b","3b","70","f1","f4","82","19","84","d7","dc","02","2c","69","2f","1d","0c","ab",
"47","b9","08","32","04","7f","05","53","8a","07","ff","49","a9","df","d9","cb","a5","6f","d0","a9","ae","29","29","a7","37","09","5f","3a","5a","6b",
"91","6f","f1","34","87","7d","a2","69","31","5d","f0","ec","09","1e","ab","3c","a1","ed","e2","ba","e6","72","b3","b2","f0","df","0c","f9","bb","cd",
"4e","ff","d3","e2","98","77","8d","5e","2d","39","01","6c","65","52","e1","a5","3c","2a","4e","ef","2e","54","2a","59","3e","96","4c","4f","b7","e6",
"ff","61","bc","f9","60","47","21","42","f6","27","9e","a8","0e","a1","f6","ed","75","bc","e3","08","84","3f","44","33","28","e9","75","eb","80","04",
"59","cb","7f","38","41","ea","99","c8","0c","c7","c2","66","47","4d","3a","05","ad","b4","fe","35","93","f5","7c","d5","7c","6e","c1","ab","fe","c0",
"a3","6c","e6","57","30","0e","bf","d2","7f","70","6b","0b","a8","bc","62","c6","b5","c6","ee","ef","69","40","18","78","98","73","f8","25","b0","0e",
"47","59","15","db","d3","51","53","72","0f","2e","de","45","7f","00","eb","c4","a7","aa","32","97","8a","6a","e6","3f","15","d3","da","3f","a7","ed",
"5b","82","08","6e","04","e9","84","fc","a5","b4","32","bc","68","0c","87","43","4c","6c","6a","42","f2","69","14","74","6a","d7","67","b1","23","e2",
"9a","a0","9e","f9","db","e0","f1","19","72","98","66","ae","0a","dd","1b","a4","1a","46","34","85","7d","75","16","94","6c","f7","14","80","01","a3",
"13","8a","60","63","df","87","96","37","cb","3b","61","6a","da","93","d7","0e","a4","51","2b","fd","88","a3","1e","11","68","81","47","36","02","f9",
"08","4d","95","a9","de","20","23","54","e8","f1","ab","ff","68","00","46","f9","52","1d","67","a8","a8","87","70","8a","16","0f","e7","d8","51","ef",
"c7","04","2d","20","d0","68","fb","19","25","58","cd","9d","21","50","44","55","a2","77","3f","fa","f0","cb","97","46","b2","7b","80","38","3e","93",
"9f","6d","25","c6","14","72","93","84","9b","d4","81","7c","c2","2e","63","83","5b","ce","09","9f","b5","35","48","fc","58","ca","5e","ca","71","7b",
"a2","00","f9","ec","af","f9","93","b4","dc","46","60","da","24","53","4b","0d","fe","67","49","13","20","b1","5b","58","1c","7c","04","b2","5d","59",
"77","ca","2a","d8","08","6f","ce","02","0b","7a","d8","8d","7e","a9","c1","68","12","a0","92","dc","40","d9","bc","c7","03","26","6a","cb","a1","fd",
"1e","ab","9e","2c","44","1f","1e","36","1d","d6","47","83","80","e0","18","3c","ff","0f","6d","a0","2a","26","97","45","8f","07","41","e6","c2","28",
"94","a5","81","64","b6","a8","43","3d","05","7e","96","76","4f","a6","cd","b3","5f","2a","b1","e8","00","5d","e3","8c","3e","c8","6d","0c","43","a9",
"c1","33","bf","76","17","ed","83","d2","3d","21","9b","ba","72","8f","7f","a4","1a","e2","57","b4","ce","35","66","04","88","80","91","27","60","9f",
"24","1e","81","0a","e2","a4","6b","32","ba","b5","27","0a","70","26","db","70","d8","d1","af","8e","0a","8a","8b","a8","52","b4","89","7c","6d","d3",
"32","a4","58","d3","58","86","dd","e7","fb","4e","cb","1e","31","f1","a5","24","63","84","44","28","f5","79","cd","99","60","ca","8c","e0","71","24",
"d4","79","22","58","d9","bf","5b","85","7c","7f","0b","1d","18","fa","d0","9e","23","d4","25","8e","7c","0f","a8","31","2b","9c","59","e2","8b","ef",
"c3","a9","a7","bb","b2","39","9d","17","16","fe","ce","0f","d7","44","fa","db","ec","b7","18","da","51","59","c5","b1","f4","b0","2c","ca","3d","39",
"d5","72","f3","79","5a","c5","ad","65","a1","47","05","39","d4","d6","6a","15","38","97","a4","be","bf","e7","09","75","4b","5c","2d","b6","fb","27",
"92","0d","3b","06","63","4b","c6","28","a3","3c","f4","3f","49","b4","e7","8d","2f","7a","69","61","e6","f7","ae","9c","cd","c8","92","65","1d","09",
"40","78","eb","2c","72","f3","ab","3f","bf","22","1b","5f","d2","8e","19","58","5d","87","6c","e5","87","f4","55","15","cb","4e","d6","c7","e4","64",
"81","f2","15","a9","00","d6","82","5d","38","9d","02","86","f3","e1","1f","4f","e7","83","1a","82","a2","55","7d","4b","38","5a","dc","73","44","12",
"48","92","38","ad","c0","b1","ce","18","f2","6c","a7","da","6c","d7","2a","0b","1b","3c","74","4b","2a","fe","b2","d7","4a","f4","28","40","54","6f",
"a1","77","3c","6e","5b","af","02","ef","3c","ee","26","13","4f","ac","4d","29","28","f7","d8","22","8d","f7","66","2f","74","d5","26","b2","0c","be",
"f9","b0","7f","57","77","d3","61","28","70","4c","e5","99","11","db","42","03","8a","5b","a2","c6","b1","39","79","44","8f","f4","fe","85","71","de",
"a4","c0","9d","57","be","aa","ea","a2","4e","60","85","ba","9a","16","c1","61","f3","9a","f9","0e","32","53","99","a4","ca","a2","ec","47","c1","75",
"de","47","1a","f7","80","ea","b1","e9","66","f0","c6","78","a3","1f","13","fa","d4","64","30","f4","fc","fe","4e","47","16","d9","3c","39","7e","30",
"89","dd","17","8c","72","21","15","30","dc","5f","de","13","45","ed","12","16","f8","9c","13","09","9e","08","1e","64","31","77","1a","b7","f2","8c",
"be","3d","47","78","40","fa","8d","46","6f","f4","c4","31","f9","a0","e1","85","1d","d7","07","dc","cb","42","f0","24","14","1c","b9","c2","79","8b",
"9c","ee","f4","6a","fb","96","dc","43","c3","f9","da","96","f1","ba","a1","e0","f5","be","99","ad","d6","95","2a","b3","cb","90","31","1c","0c","b4",
"fd","c0","5f","0a","d4","72","36","ca","6f","b8","f1","12","b0","f2","e1","e7","f7","4c","b2","bb","fd","b6","c1","93","ba","a9","87","10","e5","b0",
"c5","e6","03","49","b7","8f","cc","11","98","a5","c8","b5","d6","d8","a0","ea","3a","c5","c2","24","bb","e4","12","4e","0a","06","31","0a","b9","16",
"3d","5b","72","4d","8e","79","75","63","7c","cf","a1","c0","b6","6c","49","2c","72","66","ec","ba","38","4f","d7","a8","3b","ac","ac","f1","d0","40",
"16","b2","47","60","a4","7c","6c","bd","ad","28","ac","58","27","84","9b","fd","a6","0f","d7","c3","02","22","d7","04","fb","65","32","31","de","7b",
"45","6c","36","bd","a5","88","73","66","1c","8c","43","3a","1a","97","19","37","27","76","f0","26","b6","df","0c","2a","58","a9","17","ec","4e","e6",
"33","db","32","b2","19","d6","21","f8","5b","70","92","a4","a6","10","e1","fc","9a","84","d9","0c","64","39","42","7b","af","8e","be","7b","e2","fe",
"db","38","03","47","ea","2e","40","e0","1c","2d","05","41","a4","84","bf","91","22","41","89","9b","f3","25","30","e2","64","9b","64","e2","f8","6e",
"b0","9c","79","e7","fb","6e","cc","61","43","ca","ce","b9","15","7c","1e","11","b4","b3","a6","97","83","77","92","2c","55","6e","15","7f","a5","83",
"97","d6","23","32","34","54","21","ad","c0","c7","31","ce","af","d4","21","f0","01","d2","0b","10","a3","54","0d","bd","90","5d","0c","f0","af","cd",
"69","ed","d7","3d","84","69","d3","c0","58","b8","7f","06","be","01","b1","cc","55","9d","83","bd","02","67","75","b3","19","b1","a1","2c","8b","d3",
"f0","f2","6b","f5","09","50","8c","bd","b1","d8","e1","71","48","4f","0c","4a","a1","76","73","f6","09","55","12","61","d9","a6","8c","dc","7e","b5",
"68","a5","7d","3e","21","41","0b","d4","59","81","03","e3","33","ec","0b","a5","44","b5","77","84","68","fb","9d","26","d3","60","3c","74","85","35",
"ef","c7","71","c9","ce","76","a3","09","0d","7c","be","d4","16","93","8e","ec","b3","5d","39","23","ce","fa","64","b1","51","2d","96","2e","87","f1",
"52","bc","1d","70","f9","4e","3a","d4","b5","b7","78","b5","73","57","6a","79","0a","4b","f2","95","7e","41","00","05","40","b2","64","94","d6","18",
"d3","16","b2","e1","07","ba","e9","63","61","0c","76","43","a2","51","35","c9","2e","7e","07","eb","fe","1f","19","20","0e","21","71","04","f2","e0",
"96","fe","66","be","74","af","55","2d","70","08","db","c6","c7","47","0a","a2","52","60","e8","e0","1d","a8","67","2e","be","62","47","87","c3","0a",
"c3","48","75","02","fc","92","8e","76","0d","1c","34","09","c4","a2","b1","bd","ea","af","46","d8","33","0f","aa","99","86","76","bc","01","dc","03",
"14","db","a3","e0","22","25","44","1c","19","1b","29","02","1d","37","39","59","5a","48","5c","aa","be","0c","4b","41","50","bd","de","cb","47","dc",
"c6","a2","3f","35","49","46","63","55","74","57","24","1c","32","31","7d","e0","84","f8","b8","a1","0d","38","bb","d9","dd","a3","20","17","89","06",
"40","e4","15","0b","0d","1d","93","78","7e","7e","50","f8","70","80","3f","0a","f3","55","55","1d","c3","31","e5","13","25","4d","c3","78","0e","84",
"5a","9b","c7","e7","51","ed","27","9b","95","7e","c0","c2","be","f9","f9","a2","52","2f","10","23","04","35","e2","44","fa","c3","ba","71","ef","17",
"5c","98","88","95","27","4f","15","0c","a0","f5","32","d8","36","a1","0a","19","40","1e","ed","a9","39","a0","b7","b0","2d","79","d6","ee","d6","ae",
"4e","72","da","2b","1c","49","81","b4","8c","a7","04","48","69","ab","44","a4","20","b5","7d","66","ef","a4","fe","b2","49","b0","2d","04","52","8a",
"16","7a","7e","44","a8","a3","01","7f","1b","fa","2e","cf","71","2d","fd","cf","af","7e","15","27","ab","e0","33","3a","e3","8c","61","3a","66","69",
"32","af","d9","d3","15","b6","bb","ba","a0","c1","49","7f","0c","48","68","db","4d","4c","c2","ad","a5","e4","c9","8c","30","c0","bb","93","6e","5b",
"1f","5b","9d","91","f0","57","18","9b","eb","24","a7","72","e9","c2","1a","e5","e2","d1","ae","79","23","91","f0","d4","57","ed","d7","ef","25","51",
"95","13","01","8e","ac","22","10","de","77","19","7a","06","37","7c","0f","88","b9","fd","a0","5c","5d","55","7c","30","c2","5c","a8","a6","df","25",
"e1","d9","f6","c3","0a","b4","36","60","67","d6","29","07","44","94","26","b1","00","03","65","8b","0e","8d","a0","22","9b","c4","cf","af","2d","2c",
"a0","b5","04","9d","94","19","3e","78","00","b9","8c","13","3b","82","3c","9a","98","99","b4","2e","16","c1","a3","c7","f5","3d","a4","5f","c2","ea",
"55","54","fc","5c","9c","93","45","40","78","8d","6b","a5","85","04","0d","02","43","d8","46","62","33","14","d2","05","0b","37","99","e0","92","9f",
"c8","59","e0","7a","b3","68","08","d4","85","c8","1e","5d","bf","ef","5e","21","6a","4e","9a","da","77","8d","95","2a","53","c4","19","d2","76","f0",
"fd","ba","a8","23","15","be","c1","70","2f","74","53","fd","95","d7","c8","37","13","17","5f","fe","58","a1","24","1b","72","af","2f","fc","5d","31",
"da","4b","16","2d","75","5a","ac","21","de","0c","db","13","2d","05","18","e1","79","8a","85","14","6b","b9","bc","2a","01","cc","cb","c2","78","9a",
"2f","65","f0","c7","81","9e","6f","57","c5","87","69","48","3c","20","fa","f1","e8","41","dd","ed","e5","29","7f","a0","6a","72","0a","bf","0b","7a",
"50","0b","59","42","cd","9e","4d","7e","53","0a","21","40","8c","e2","fa","db","6a","84","60","0f","2a","68","9e","0f","c8","5b","07","a0","f5","1a",
"35","75","2d","89","0b","a2","14","05","c0","c4","c8","cd","1a","6a","af","e3","29","d4","c3","72","08","b2","7d","1c","28","b0","27","1b","12","07",
"4b","00","3a","a0","78","8f","d5","b6","f2","50","ef","11","fd","b6","05","00","23","ce","bb","2e","3a","95","63","83","06","56","6f","07","d9","25",
"4b","32","42","f1","c0","62","c5","fe","80","96","b5","72","e6","fe","12","77","01","3d","ea","38","50","8e","ad","a5","c8","2a","5e","7d","5c","6e",
"37","33","7c","34","ca","b3","50","20","66","20","7f","5f","94","b6","7c","3b","2b","11","22","4b","f3","ce","29","b4","88","98","52","24","45","c4",
"fb","e5","b9","4c","da","10","04","b7","3a","4c","6a","78","4d","53","4d","2a","b4","cb","6d","0c","17","d1","49","45","7e","49","a2","1b","1c","74",
"9e","1f","5e","42","5a","df","b7","60","89","ba","07","3a","92","15","a9","6c","17","1f","06","62","2e","fa","44","75","cc","1f","75","94","72","6a",
"cb","01","38","6d","24","0b","20","0d","34","66","c2","63","84","68","02","e8","ad","0c","f8","75","b4","00","40","a9","5b","03","fa","a9","d0","32",
"e1","c9","29","11","2b","3c","59","3a","5c","16","11","f8","3e","e1","31","8a","1c","28","d0","46","ab","fe","4c","b5","e4","0c","ec","e3","f5","fb",
"f4","7a","cc","4e","d5","d5","8d","90","66","f5","54","24","0a","99","18","1e","65","a4","4b","9c","75","c7","a3","84","3f","b7","dd","91","53","dc",
"3d","63","3b","1d","c7","56","d0","4b","12","4a","03","43","c2","b9","00","5a","4a","24","73","38","d2","d5","74","6b","f6","07","9d","87","e7","06",
"97","58","46","fa","a4","7d","50","00","df","4b","c7","48","0b","9e","9b","aa","8a","74","34","18","d1","eb","fb","b3","92","d9","73","c8","c0","6f",
"1c","71","c0","f9","d2","92","6f","34","40","90","48","99","72","ac","4c","70","f3","2b","11","2b","4a","58","4d","21","cd","19","54","2c","13","25",
"6c","c1","98","35","47","18","6a","b6","5b","b9","4f","a0","e7","95","8b","d8","d9","1b","3c","b9","c6","38","70","17","0b","de","d6","d9","32","4b",
"a2","d4","d1","09","1e","29","18","eb","00","a4","f1","71","77","a7","dd","c5","e6","b1","ae","b0","c1","88","e9","8b","21","f4","97","b6","20","6d",
"12","f7","b7","e8","e2","91","ab","fe","31","c9","10","e5","45","5c","43","f0","3d","23","7a","5a","d6","54","fa","29","22","f3","6d","08","73","28",
"fd","29","54","38","2b","55","0d","26","00","c0","00","c6","e2","16","cf","34","3a","a0","b9","b5","bb","52","2e","c0","30","29","d4","73","cf","70",
"f3","81","f0","bf","a8","22","9d","92","f8","f4","b7","0a","40","bc","ad","bf","c0","7f","dc","e7","f7","45","50","05","a8","ac","c0","82","e5","f7",
"a3","36","1a","ea","21","7e","dc","a9","bc","7a","a0","a2","43","f0","00","92","2b","60","67","ca","05","fd","6e","71","7a","26","9e","d4","b4","f0",
"6d","01","e8","ea","1a","f9","a3","a3","03","04","34","a9","aa","2e","5f","d1","e5","37","be","a6","5a","0c","80","5d","4a","03","78","21","7c","aa",
"02","58","cf","a7","29","6c","87","86","d1","00","cd","f9","36","66","d3","66","2b","98","07","9d","a2","00","d3","b7","c8","90","a4","fc","1b","d6",
"bc","e9","2c","db","84","82","9d","f2","cf","b8","41","d4","e4","0e","0d","a5","8e","94","25","a8","55","fd","77","30","db","1b","03","8e","64","b4",
"b8","63","59","c2","0f","d3","e9","17","cf","31","0f","63","fb","dc","82","4b","d2","c2","12","14","4d","49","33","69","22","1a","85","f7","9e","c4",
"e6","b2","68","79","04","4e","26","33","8d","0b","de","e1","17","75","a8","61","bd","b2","7f","f8","ce","a4","cc","32","5b","58","b2","a0","e8","e2",
"ad","8c","86","11","64","c2","49","76","05","a6","1f","79","6a","bc","53","8c","80","c7","8d","c2","47","01","b1","26","94","cd","27","54","33","96",
"36","d9","37","9f","b0","d3","d9","09","46","f1","ba","a0","00","07","61","b1","7a","7d","0c","ee","b4","99","b2","86","97","26","6d","be","ea","6d",
"88","af","08","b2","51","a3","19","03","65","d4","98","ae","7f","b5","cf","4c","c3","a1","7c","f1","40","58","ea","f1","70","27","88","19","cb","61",
"37","d5","6b","54","7b","18","70","24","fb","7d","7c","c5","cd","0e","be","0d","c9","27","5d","e7","8f","14","10","36","10","ea","a2","cd","5a","24",
"d9","da","e0","d3","50","87","f1","e1","4c","4c","79","e5","33","5a","97","32","a7","23","39","0c","36","b5","84","86","b3","41","0a","9e","1d","fd",
"76","bb","03","35","b1","e2","6e","20","e0","61","de","5e","21","7e","9b","d1","eb","3e","23","3b","5b","06","80","c6","b0","bd","b8","94","48","85",
"ae","83","1c","6a","88","31","c5","1f","0f","0d","a0","45","64","fb","de","08","78","58","0b","7d","2b","71","50","c2","e8","3d","de","35","6f","77",
"e2","13","44","2d","20","2c","05","b0","02","b5","79","29","74","b9","f0","0c","2e","0f","64","1e","39","85","74","46","e6","16","8f","78","3b","db",
"ba","97","27","bc","d5","62","0e","a0","57","63","b1","2f","60","3f","5e","a3","e2","c6","c2","e6","43","63","22","11","5f","ee","56","98","77","1e",
"79","37","5d","b5","c1","3d","7a","79","d3","4e","40","85","e2","c5","9b","bb","3b","00","a9","5f","5d","6b","f2","17","d9","b9","d4","3c","7c","6b",
"03","74","e0","3e","c5","83","9c","82","24","cd","ea","39","4e","40","c1","29","e8","4f","60","9d","62","3c","27","30","a0","20","25","7f","b4","95",
"bd","25","2d","6b","74","a3","01","d6","18","a0","d4","da","79","49","42","50","ff","3b","3c","44","f8","40","df","08","30","e2","5a","ce","b8","57",
"ef","32","60","0d","78","78","56","31","3b","52","dc","19","a9","8c","9f","fa","dc","a7","3d","07","51","ae","df","79","63","b6","ec","d3","25","4a",
"f3","ca","fb","c2","cb","9f","81","4f","46","6a","34","e4","ba","80","d9","4c","9a","51","4b","2a","ed","56","7b","78","6e","b0","35","8b","43","fc",
"b4","e6","75","f7","c4","2e","5c","b6","1a","c6","b5","33","79","57","2e","bf","d8","a3","1a","f4","64","48","1e","75","40","ec","97","1e","fa","c9",
"12","7e","d0","a6","f4","cf","9d","68","17","f8","39","af","d6","c8","8d","c2","a4","76","73","80","9f","2a","1e","74","3f","59","c0","a9","73","9e",
"df","cb","3d","00","bd","39","56","64","4b","e5","50","19","17","16","e2","85","bf","b5","22","14","e7","d9","ee","af","03","1a","0f","6e","66","92",
"4f","a5","0e","6c","32","35","68","e1","8f","7e","6f","f7","d7","68","73","4d","98","85","41","84","b5","b6","f4","ab","2e","47","4a","92","a6","97",
"37","6f","8f","8f","79","74","52","19","3b","2e","fc","d2","1a","46","29","d1","65","da","c2","e2","b6","aa","b3","61","e6","b0","d3","c2","3a","48",
"8b","56","b6","a5","f4","0f","99","03","78","b3","10","de","8c","64","fe","a2","d0","31","23","6c","a0","7f","5c","96","48","23","23","81","a6","6f",
"66","50","a4","7d","cc","56","e7","2f","38","ce","ff","78","f7","d8","9d","46","95","51","46","bb","14","7c","92","76","f3","ae","41","09","75","14",
"a3","b9","d5","4c","16","49","5b","22","f8","89","7a","5d","68","09","5d","26","62","4d","ed","35","5e","70","5e","91","f5","88","8b","01","14","2a",
"0b","2d","6d","07","d2","05","b5","c7","1f","2e","b8","bf","7e","42","55","87","7f","a0","a7","2f","ff","e1","a6","75","7d","72","03","dc","94","32",
"e8","e2","52","cb","9d","3b","21","0f","fb","13","05","f7","c9","23","47","a9","cc","02","c9","68","21","77","da","88","2e","bd","26","dd","29","d2",
"8b","6d","58","11","6e","3d","5f","09","05","c4","32","f4","81","8f","37","26","36","6a","68","0b","59","7c","87","9c","b4","91","02","cd","64","08",
"c5","2e","cf","4e","c6","7d","ba","77","c7","d3","3d","03","9e","ae","20","c8","1a","3e","0f","e6","0a","2e","0a","0f","69","03","66","7e","ee","7e",
"49","d1","fc","51","cb","6e","b4","01","47","35","55","57","1e","94","b1","f2","f2","97","85","16","89","20","f5","e6","93","e6","64","93","f1","7f",
"8d","cf","e6","b9","85","5b","48","80","3e","13","e5","92","62","fd","1c","b4","5d","66","bb","4a","85","b6","6c","9a","26","89","1b","78","ee","4c",
"0e","b4","79","1b","84","9e","ea","a4","26","ee","7a","cd","56","d4","bd","05","9f","40","0c","02","18","26","c2","51","80","9b","22","c5","da","4f",
"f2","23","d4","59","49","88","d6","03","52","ae","9d","ca","9f","3b","5e","a7","b1","4b","11","3e","fb","0b","2f","22","cb","c5","34","1a","77","be",
"97","5e","65","ca","33","83","db","57","44","69","80","72","bf","89","d0","c7","a3","55","a3","f3","2d","ab","59","6d","b2","1e","ae","eb","27","31",
"ad","3c","b6","ef","0a","f4","c5","21","f7","a4","70","a7","4e","73","a1","e8","53","a1","9f","c8","3a","12","89","8a","62","28","a5","2c","d2","35",
"c8","64","ed","61","22","be","0c","4f","87","33","e5","7e","3d","3f","0f","39","f8","0b","aa","0f","67","57","6a","0a","77","0b","4c","81","e2","60",
"94","79","1d","5e","00","60","c3","e3","75","ad","31","0a","3b","c0","e6","2e","41","c8","8e","b6","5c","37","ce","d7","56","df","1c","bc","f3","d7",
"4f","d0","7e","b6","4c","e8","67","4a","66","76","c5","c3","86","7e","ad","fa","95","e4","61","19","f2","c2","a5","5a","72","e6","61","05","18","fa",
"ff","2e","00","3f","d9","b9","c4","00","a5","d4","f2","e9","94","37","a7","38","68","c3","db","d5","20","39","af","08","95","87","fe","4e","74","34",
"f4","6a","95","fa","76","2f","68","16","53","72","00","8c","1c","b4","6b","e2","fc","7c","ab","3c","21","31","84","5a","46","2a","c8","d7","a9","89",
"84","f2","77","dc","de","f0","d8","2e","5a","ed","db","b7","e4","f9","a7","ab","a0","e2","fb","34","62","4d","aa","68","87","40","b0","30","c8","be",
"ab","f9","bf","a5","ea","1a","9d","c2","10","40","8c","81","38","e0","41","c4","1f","cc","2c","0c","9b","d7","51","cd","fe","38","00","a3","1c","d3",
"d2","de","74","e1","4a","5c","3c","7b","fb","58","14","c3","47","ca","12","1c","75","ba","69","d3","2f","8d","f8","d3","48","22","a4","ad","97","fd",
"28","08","e1","67","43","e7","41","f9","fb","05","9f","90","09","4b","3a","16","77","4f","de","21","95","be","c8","71","fc","be","4b","02","79","9a",
"ee","19","4e","03","e6","58","b4","02","df","72","ff","82","52","5c","19","51","38","b5","52","6a","e5","28","93","8b","02","7f","bc","cd","48","03",
"21","ec","c8","78","00","75","28","1a","16","01","59","d7","7b","fe","82","b6","78","80","6b","9d","b1","b4","96","0b","79","5e","fc","e5","00","a4",
"f3","3a","2f","4d","80","8f","fb","56","59","5f","99","bc","58","f1","3e","c3","e5","d1","99","64","6b","cb","47","57","50","e6","0f","19","3a","34",
"13","07","9a","71","5b","99","82","c2","18","6f","55","09","05","e6","80","11","95","7b","8c","e8","8d","3e","8a","5d","d9","ed","46","62","b4","cc",
"6f","8d","b6","6d","04","41","ed","95","c1","54","8f","c8","29","4a","1f","a0","49","ca","ea","65","4b","60","3f","08","6d","60","db","68","f7","16",
"88","6b","24","b2","1c","0c","57","95","df","a0","0f","4d","e4","1f","9e","fc","e2","30","fe","1f","3f","1e","42","21","e2","ca","d7","4d","15","06",
"66","7b","60","38","38","5d","17","0a","be","e0","6c","13","dc","aa","ed","6f","2c","2c","ef","88","34","fd","41","d4","ea","78","54","55","5a","f4",
"d2","16","ce","3a","23","f3","7f","81","ff","40","66","6e","4b","a0","f5","f0","f0","c6","3c","4d","e0","6f","eb","8e","9e","05","48","ff","f4","26",
"0e","34","95","44","95","6e","f1","6f","22","4a","64","3b","19","78","b8","5b","4d","f4","65","4c","8c","bc","22","a8","b3","5d","6e","ed","e4","da",
"0a","6c","06","ea","3a","56","ec","3c","55","a2","23","c4","2e","5c","84","6e","00","98","48","0d","8d","af","fa","91","e1","b9","6f","0a","84","13",
"1c","e1","a2","bb","21","21","01","18","05","72","b5","e9","e0","f0","fb","77","af","4a","fe","47","9e","c8","15","59","dc","40","30","31","8e","5c",
"9c","6f","b9","eb","da","fe","ab","32","71","9f","50","bd","c0","a0","b2","cb","35","67","ba","b8","18","60","0b","7f","85","93","89","d1","18","66",
"6a","0a","e0","3a","1d","21","74","13","a7","c2","8c","ec","4c","10","d6","14","0a","89","6a","86","d1","93","09","41","40","4e","e5","c4","8c","07",
"71","f4","9f","bb","d2","58","cc","fb","2c","52","c8","f1","4a","4b","39","0a","fb","9b","ea","dc","1d","18","f7","92","f2","1d","69","c2","d2","df",
"dc","5d","40","77","30","1f","cc","c2","53","98","a7","b0","79","05","49","d5","42","ec","15","6b","6a","bc","45","3e","54","22","ea","52","62","ca",
"9b","42","d0","36","89","d8","21","2c","2b","20","3b","2b","92","dd","e2","b8","7a","8f","e8","ef","fe","47","12","e2","dc","e8","2a","8c","40","ca",
"c0","d5","28","82","68","53","d8","52","f9","8c","15","e4","b2","07","1c","01","bc","90","f8","df","1c","ef","f9","b7","52","fa","ac","46","37","f9",
"75","bc","50","27","17","bc","34","a2","99","1a","32","c0","28","0e","1f","16","d7","d0","40","18","80","4e","0d","d4","d8","4d","da","eb","07","d6",
"5c","c3","aa","e1","cf","01","44","21","40","f4","58","be","70","88","ea","b4","ac","16","a0","50","46","ee","67","ea","9b","0a","2a","c2","0c","5d",
"02","3f","38","ab","26","94","67","25","5a","bc","4e","22","ed","0c","2e","52","81","a5","04","47","41","82","e8","7e","fb","ae","bf","f9","96","6f",
"ec","26","11","a0","ba","ec","68","f8","d0","9c","76","59","02","d8","83","a8","ae","5f","22","5c","41","ae","2f","15","bf","72","d4","f1","55","3a",
"2f","81","51","74","66","54","ca","e3","20","1e","77","cf","55","fb","7f","06","56","ac","d3","8a","56","dc","6f","08","cc","29","7a","a2","61","d5",
"37","f7","1a","2b","6e","06","4f","b3","83","c9","7d","1d","a1","15","12","0a","5a","c3","be","d1","3d","6a","90","18","ba","4c","8d","e0","da","be",
"11","02","5b","af","fc","c4","5d","91","77","e8","55","80","fa","59","78","0e","c4","61","a2","1a","c0","d1","c3","d5","85","9a","35","bc","68","85",
"93","76","89","64","61","c5","69","90","94","09","4d","80","0d","6f","1e","89","f3","5f","b5","94","67","1f","a7","44","c6","47","1e","84","68","31",
"ee","32","d8","f1","5a","73","d7","63","d4","bc","ba","1a","1f","06","c5","83","ea","93","9b","f8","3e","97","9b","68","5b","48","49","24","ea","ce",
"77","75","63","e0","e0","b1","11","07","d3","ce","31","5a","bb","0f","d4","d1","34","4d","59","12","2a","d4","74","a4","03","c4","e7","9b","22","fc",
"00","0b","b0","36","73","0d","f9","4e","4a","47","f0","55","d3","0a","1e","77","9d","b6","2a","f2","12","b0","9e","7e","29","af","c1","4a","fe","9c",
"83","a4","38","3b","0f","85","98","56","ad","60","ba","6b","98","d9","57","ee","24","db","f9","4e","1b","19","9d","81","19","83","e5","9b","b8","e4",
"d0","af","c3","fe","8f","51","d4","c2","9c","31","64","f8","9b","b8","9f","d4","b5","77","a2","ed","64","ca","4c","67","c5","24","26","13","b3","d7",
"4c","ac","ab","af","11","d0","4d","6e","78","fa","5c","cf","a7","1b","4b","ac","44","ef","20","51","de","6e","2b","39","7a","78","37","71","bb","13",
"e2","65","57","a9","82","bd","e1","9a","2a","3e","a4","ab","18","d8","cf","98","db","b0","70","6f","2c","83","d7","31","8a","39","03","17","cc","2e",
"68","25","7d","79","fb","9a","aa","8a","c1","40","a6","59","3e","48","bf","36","bd","76","4b","a5","ac","f0","16","2a","82","b7","7b","a2","3d","ea",
"1e","81","d2","d7","73","4f","44","97","8e","65","84","4e","bc","82","b7","69","c2","e8","af","15","bb","5b","02","90","11","36","af","53","48","7f",
"09","f8","2c","0f","52","f9","ad","03","89","ee","db","57","cc","d9","22","d7","38","d2","ae","eb","1c","e5","b4","df","c0","17","31","97","09","d7",
"43","66","d1","20","9d","ea","46","36","11","f8","a9","a2","08","1c","fb","f1","2a","43","c5","93","7e","1e","dc","1e","18","fe","57","8e","96","3e",
"0a","e0","0e","8b","29","87","e0","4e","40","9a","6c","09","96","b8","27","0b","18","52","ee","d0","55","e4","2d","92","fe","7f","cb","49","ee","20",
"b0","e1","eb","0f","f1","e9","4d","55","af","1d","06","78","4b","4b","b5","58","93","48","c3","a6","ac","f8","e8","4c","c2","6e","4b","70","e2","4b",
"7a","6f","e6","74","cb","8b","45","3f","cf","88","ad","3e","76","44","b3","8b","59","94","8e","8d","77","cb","b8","e6","3d","35","c0","11","d0","c3",
"e3","05","41","d9","de","65","64","7c","42","89","6f","fe","a0","3c","9d","6e","be","ad","2a","99","c8","98","0b","6a","09","16","d2","a8","e2","a7",
"17","2c","94","37","be","bb","2d","cf","5b","00","bb","f7","84","a9","be","60","36","8c","fa","d3","a8","ab","c0","b4","51","2b","26","02","e0","78",
"80","6e","46","56","c0","b9","b2","a5","ee","4c","62","b5","06","69","87","bb","08","51","0a","d2","eb","e7","b6","5e","08","61","50","98","be","db",
"a2","c3","1d","95","3e","66","0c","3c","4f","90","e3","81","d0","82","c0","36","59","78","f4","43","f2","80","f2","7d","84","42","15","5d","08","91",
"89","b2","8e","bb","fa","a3","98","22","66","ab","04","1b","b6","47","56","1b","95","11","80","5a","8a","26","56","a1","7b","26","c3","d3","e7","09",
"3c","aa","78","b3","9d","fd","fa","ff","68","81","ab","81","9f","22","aa","87","77","26","d4","86","47","3e","75","87","3a","44","b6","e0","e1","ff",
"86","51","63","3f","99","2b","98","9b","62","ce","12","f4","f4","a3","10","8e","0a","08","44","b0","02","fe","b3","3c","02","83","c7","43","7e","fc",
"0c","e2","ad","db","ab","0f","6b","be","35","18","b7","85","61","82","6e","16","35","e4","d3","17","ba","fd","02","ed","01","66","0e","3f","d1","c7",
"fa","99","3e","8d","9d","e1","11","2e","64","fb","33","ff","cd","93","73","43","1b","a8","01","f6","db","d7","b6","c0","57","84","49","e5","e1","cb",
"6d","1f","d7","90","71","e7","34","62","45","a8","42","51","50","82","6b","bf","59","a9","c2","d3","11","77","f8","33","c2","6e","f9","72","08","99",
"1d","e1","87","31","ef","2e","fc","61","b8","21","4f","db","5f","88","b0","74","1e","d9","23","0c","fd","cd","ce","20","02","2a","62","0b","87","cb",
"c6","29","cc","ae","fd","9e","6a","91","69","91","8f","24","4e","f9","1c","95","df","51","44","10","cc","92","0e","c1","11","4c","c2","0f","27","51",
"58","40","1e","ab","0e","99","e6","ce","56","dd","93","70","8e","b7","6a","18","5a","3f","bc","b4","95","e0","67","b7","66","7a","b4","4e","ba","d7",
"b6","dc","63","68","d1","8f","73","59","53","48","61","d0","8c","37","d6","92","c3","e9","8f","7e","25","8f","67","ce","4c","19","50","46","91","64",
"b9","c1","1a","dd","10","6f","3a","1f","c9","4e","17","7f","1c","0f","e0","4e","fe","c2","75","60","05","0a","5a","33","a5","0e","fe","5b","27","13",
"94","75","ff","d4","1e","65","e1","d0","b4","58","76","12","59","43","26","e4","03","62","4b","03","2a","48","d7","f6","47","5e","ea","7b","46","a9",
"c4","57","a9","a3","4c","d1","3d","cc","91","69","b2","df","c4","01","9c","cc","1a","b8","02","63","ed","cb","78","c9","da","6a","e9","c7","75","99",
"0c","41","3d","b0","e0","f7","bc","76","52","db","3e","9e","1e","95","57","6b","22","46","50","fa","11","dc","20","31","56","94","fb","eb","df","ff",
"1c","c2","50","b4","87","85","5b","8c","30","1d","ef","1e","41","75","69","ee","4a","4a","81","87","96","60","32","44","00","09","89","c5","bb","d1",
"01","93","81","f0","1b","da","af","f5","e4","c8","8f","63","10","88","87","65","b9","fd","04","1a","00","7f","44","2b","40","74","1c","77","fb","5b",
"d6","e7","f9","54","1c","34","93","49","af","75","44","f3","4b","95","29","10","86","b1","57","50","98","7e","9b","97","30","6c","5e","b8","6b","f8",
"16","ac","06","d2","b0","8d","00","0c","a9","52","be","48","f1","56","e3","1e","e5","7f","23","c9","ec","02","e3","42","a0","b7","d8","4a","a3","c3",
"ea","91","1b","49","68","a1","06","41","5c","1f","43","a2","7c","0e","05","1d","34","dd","62","a5","74","a4","3d","df","71","52","f0","98","ad","6c",
"3e","1b","55","c8","64","64","44","90","24","15","82","63","49","12","61","b4","a5","92","72","e2","b4","c3","ca","17","e2","81","3d","af","16","ff",
"70","f6","2b","25","df","ab","b9","55","5d","1f","5d","67","ce","b1","6a","00","1f","02","db","18","00","e1","85","19","24","b7","3b","50","c6","3c",
"dc","e8","93","92","8e","43","5a","a6","7b","f9","cd","ae","4c","df","6b","0e","d4","e2","8b","2d","e0","1e","37","6a","ef","72","0a","5c","a0","fa",
"08","ef","1d","fb","18","d7","0b","3a","f0","1b","99","5e","8c","49","b7","29","7c","48","0c","e9","11","72","95","30","34","96","c5","b8","37","c6",
"41","26","7d","47","c5","1f","4c","5e","0f","3f","97","57","ae","c6","1c","4e","e8","76","07","e7","5e","4a","62","9c","37","e3","bd","22","4b","94",
"9b","f1","c2","52","4f","45","6f","7f","e7","6c","d6","ae","04","d8","68","ce","10","04","60","45","88","a9","22","df","25","21","1d","b4","02","e1",
"c3","2a","a0","33","eb","54","11","eb","1b","15","70","ae","57","49","38","fe","dc","8c","ce","60","90","c8","87","b3","ee","c7","c9","23","d1","6c",
"dc","21","dc","9d","7a","a9","b7","e9","4f","e5","a0","9a","6b","74","fa","41","d2","5e","62","a4","97","a4","b5","dd","27","cb","d1","2a","32","4a",
"84","42","1e","38","71","e8","23","9a","9c","45","91","d0","17","eb","29","ed","bb","55","6e","de","e7","5f","ff","32","fe","61","97","40","49","b9",
"86","c6","8b","d0","3c","69","53","bd","a0","11","6f","6d","a8","1e","c3","4d","5a","e1","2e","f9","bd","05","8d","4f","69","c8","74","03","fb","af",
"21","d5","68","63","04","d8","83","fd","fc","50","e6","b7","84","77","36","64","3e","ba","b3","40","99","3b","09","aa","a6","ae","a2","fc","91","ce",
"8c","9c","55","b3","81","c5","9a","ef","82","6f","69","73","97","ac","d4","1b","d3","e2","d1","e3","4c","14","71","36","1c","97","e1","7c","00","bb",
"ff","70","21","9f","d6","9b","c0","8c","a5","46","25","d5","de","51","23","c8","d0","b2","38","7e","74","2a","35","7b","3e","8f","f7","e4","03","c4",
"8c","d2","2f","2b","9e","98","b9","de","c5","7d","2a","28","cf","43","1f","04","71","34","66","73","11","6b","3a","08","c3","8d","7d","4d","7c","61",
"6a","b2","c8","19","79","0c","7d","61","a7","9c","ac","e5","e1","aa","7b","7b","2a","dd","21","91","72","23","6f","41","9e","19","cb","31","58","59",
"92","88","c5","5a","e4","7e","5d","a4","d7","d9","f1","37","2a","7e","eb","35","85","b9","66","c2","c2","20","59","80","37","e5","a5","e9","d8","f9",
"b0","71","03","ff","fd","e8","eb","68","ca","b7","94","0b","df","95","ed","0e","a7","cb","da","66","3b","44","29","8f","7c","ba","2e","57","fb","54",
"ee","0c","e3","fc","98","fd","fc","ee","ad","2b","4b","07","49","2b","46","e1","ee","18","84","67","87","ed","5f","66","f7","fe","aa","97","83","e2",
"02","3a","99","a7","12","2a","e0","77","a5","27","2b","de","f2","24","e3","70","9f","4d","ef","b8","eb","c0","02","b1","bc","29","c7","4f","2a","b0",
"b8","e3","65","e3","f8","c0","f2","68","29","3b","0c","72","9f","09","f7","9a","b3","f7","b8","e3","0e","de","b2","0f","0c","0c","6b","d7","fd","31",
"96","fd","2b","57","05","98","57","1c","16","64","12","92","6c","51","f0","23","33","6f","ce","d1","df","d3","81","78","ef","d5","b0","3b","59","a5",
"0a","8f","8b","e6","19","bc","82","e2","b6","84","d8","85","a0","47","17","80","d3","f7","b2","7c","95","a1","75","89","0e","a1","71","2b","0f","ae",
"47","43","af","4e","b2","e5","a0","47","74","74","79","81","72","c9","f8","0c","f2","95","a8","2d","a6","4b","aa","15","71","c3","5f","44","c5","15",
"30","7d","dc","c9","66","1e","40","9c","a1","6c","8c","38","55","0d","d0","d8","80","b6","82","f0","34","be","8d","35","8b","fe","f3","a8","44","91",
"7e","a9","b8","0b","ca","24","44","c4","cd","bd","47","c6","9c","ac","b3","18","03","2b","ab","f5","3c","2e","56","81","bf","89","76","9d","5e","02",
"c1","21","50","3a","57","05","db","ef","89","a4","e8","79","dd","ae","ca","88","10","3c","27","84","7d","9b","60","87","d6","be","1a","f5","de","7e",
"18","9c","e6","dc","f7","6d","78","4a","80","7c","5a","8e","3f","d6","17","3a","9b","cb","86","7f","bc","e4","82","f1","8e","8d","a3","55","bf","73",
"b7","51","8e","12","f3","23","26","c8","b1","5e","61","7c","b5","14","bb","bb","bf","c7","05","91","15","38","ba","17","bd","28","89","ac","ca","3b",
"93","16","c1","28","cb","8a","09","08","06","9b","6c","64","57","32","2c","f0","7f","58","be","e2","f5","2f","d4","88","fb","7d","24","83","5b","fe",
"cd","2a","21","08","21","6a","3f","ff","de","fe","62","53","0a","92","bf","30","55","35","90","fa","b3","66","62","01","6d","fd","0c","fc","8b","00",
"de","29","a2","0f","b8","1e","c9","77","7d","cf","d8","58","cb","92","ce","09","33","af","b9","71","45","d4","df","e0","5d","0c","db","98","f3","b7",
"c4","4a","f0","a0","e7","5c","25","31","e0","2e","d5","fe","d8","58","7d","95","bf","bf","6a","41","65","f4","00","64","10","6e","05","7a","31","a2",
"f7","83","e5","56","80","05","6e","27","c5","f1","28","0b","63","8d","99","7c","f4","d3","ac","78","12","69","6f","ac","42","cf","07","56","79","d2",
"32","06","8a","f8","13"))

