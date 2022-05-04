
# source('~/repos/3_Scripts/4_MiniDOASAuswertung/RScripts/miniDOAS_functions_hac5.R')


#### read doas data
read_data <- function(folder, from, to = NULL, tz = 'Etc/GMT-1', doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL){
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
        readDOASdata(di, folder, rawdataOnly = TRUE)
        , class = 'rawdat'
        )
}

#### print specdat method
print.rawdat <- function(x, ...){
    nc <- ncol(x$RawData)
    wl <- get_wl(x)
    cat('~~~~\n')
    cat('\t', x$DOASinfo$DOASmodel, '/', x$DOASinfo$Spectrometer$Serial, '- raw data\n')
    cat('\t', nc, 'entries\n')
    cat('\t', min(wl), 'to', max(wl), 'nm', sprintf('(%s pixel)\n', length(wl)))
    cat('\t recorded between', format(x$Header$st[1]),'and', format(x$Header$et[nc]),'\n')
    cat('~~~~\n')
}

#### print avgdat method
print.avgdat <- function(x, ...) {
    has_dc <- attr(x, 'dark.corrected')
    has_sc <- attr(x, 'straylight.corrected')
    has_lc <- attr(x, 'linearity.corrected')
    x <- attr(x, 'RawData')
    nc <- ncol(x$RawData)
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
        # convert from/to to POSIXct
        from <- parse_date_time3(from, tz = tz)
        to <- parse_date_time3(to, tz = tz)
        # get indices
        ind <- which(folder$Header[['et']] > from & folder$Header[['st']] < to)
        if (length(ind)) {
            folder$RawData <- folder$RawData[, ind, drop = FALSE]
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
        folder$RawData <- folder$RawData[, ind, drop = FALSE]
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


#### average raw data
avg_spec <- function(folder, from = NULL, to = NULL, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.straylight = TRUE, correct.linearity = TRUE, dark = NULL) {
    if(inherits(folder, 'rawdat')){
        if (!is.null(from)) {
            # convert from/to to POSIXct
            from <- parse_date_time3(from, tz = tz)
            to <- parse_date_time3(to, tz = tz)
            # get indices
            ind <- which(folder$Header[['et']] > from & folder$Header[['st']] < to)
            if (length(ind)) {
                folder$RawData <- folder$RawData[, ind, drop = FALSE]
                folder$Header <- folder$Header[ind, , drop = FALSE]
                # update timerange
                folder[['DOASinfo']][['timerange']] <- c(folder[['Header']][1, 'st'], folder[['Header']][length(ind), 'et'])
            } else {
                stop('No data within specified timerange available')
            }
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
        rowMeans(folder$RawData)
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
        out$RawData <- x$RawData[ind, , drop = FALSE]
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
    x$RawData - y$data[, cnt]
}

correct_straylight <- function(x) {
    win <- getWindows(x$DOASinfo)
    if (is.data.frame(x$RawData)) {
        sweep(x$RawData, 2, colMeans(x$RawData[win$pixel_straylight, , drop = FALSE]))
    } else {
        stop('Fix correct_straylight for non-data.frames')
        x$RawData - mean(x$RawData[win$pixel_straylight, ])
    }
}

#### helper function to correct for non-linearity
correct_linearity <- function(x) {
    lin.coef <- x$DOASinfo$Spectrometer$"Linearity Coefficients"
    x$RawData / linearity.func(x$RawData, lin.coef)
}


#### get Cheng 2006 cross-section
get_Cheng <- function() {
    structure(
        list(
            data = fread('~/repos/3_Scripts/4_MiniDOASAuswertung/ReferenceSpectras/Literature_absorption_crossections/NH3_Cheng(2006)_298K_140-230nm(0.02nm).txt',
                col.names = c('wl', 'cnt'))
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
                DOASinfo = getDOASinfo(info[2, val], timerange = info[1, val], tzone = tz, Serial = Serial)
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
                DOASinfo = DOASinfo
            )
        }
    } else {
        # if file is result from getSpecSet
        if (is.null(spec)) {
            stop('argument spec is unset. spec defines which SpecSet list entry will be selected')
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
            DOASinfo = DOASinfo
            )
    }
    # correct dark
    out$Calinfo$dark.corrected = FALSE
    if (!is_dark && is.null(dark)) {
        if (!dont_warn_dark) warning('no dark spectrum provided')
    } else if(!is_dark) {
        cdark <- dark$data[, cnt]
        out$data[, cnt := cnt - cdark]
        out$Calinfo$dark.corrected = TRUE
    }
    # correct straylight
    out$Calinfo$straylight.corrected = FALSE
    if (correct.straylight) {
        win <- getWindows(out$DOASinfo)
        dark <- out$data[, mean(cnt[win$pixel_straylight])]
        if (!lin_before_dark) {
            out$data[, cnt := cnt - dark]
        }
        out$Calinfo$straylight.corrected = TRUE
    }
    # correct linearity
    out$Calinfo$linearity.corrected = FALSE
    if (correct.linearity) {
        lin.coef <- out$DOASinfo$Spectrometer$'Linearity Coefficients'
        out$data[, cnt := cnt / linearity.func(cnt, lin.coef)]
        out$Calinfo$linearity.corrected = TRUE
    }
    # correct.straylight afterwards
    if (correct.straylight && lin_before_dark) {
        out$data[, cnt := cnt - dark]
    }
    structure(
        out,
        class = 'caldat'
        , straylight.corrected = correct.straylight
        , linearity.corrected = correct.linearity
        )
}

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


#### doascurve
calc_dc <- function(meas, ref, ftype = 'BmHarris', fstrength = 25, fwin = NULL,
    fitwin = NULL, shift = NULL, correct.straylight = TRUE, correct.linearity = TRUE,
    lin_before_dark = FALSE, do_lowpass_filtering = FALSE, lp.type = 'Rect', lp.strength = 5) {
    if (is.character(meas)) {
        # file path or chen
        if (tolower(meas) == 'chen') {
            # call cheng2dc
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
# TODO: add argument molar mass
# sigma2dc inverse
dc2sigma <- function(dc, mgm3 = NULL, copy = FALSE) {
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
    dc$cnt <- dc$cnt / 1e-4 / 6.02214076e23 / mgm3 * 17e3 / 0.075
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
#' @param index      A vector of time indices used to subset
#' @param to         An optional time index providing the end of a time range. 
#'                   Argument \code{index} has to be of length 1 to provide
#'                   the start of the time range. Defaults to \code{NULL}.
#' @param including  logical. Only used if time range is provided. Should
#'                   data at the range edges be included?
#' @return A subset of the initial \code{rawdat} object
filter_time <- function(rawdat, index, to = NULL, including = TRUE) {
    # check to
    if (length(to) > 1) stop('argument "to" should be NULL or of length 1')
    rd_st <- rawdat[['Header']][['st']]
    rd_et <- rawdat[['Header']][['et']]
    rd_tz <- tzone(rd_st)
    # check from
    if (!is.null(to)) {
        if (length(index) != 1) stop('argument "index" should be of length 1 if providing "to"')
        # parse from/to
        if (is.character(index)) {
            index <- parse_date_time3(index, tz = rd_tz)
        }
        if (is.character(to)) {
            to <- parse_date_time3(to, tz = rd_tz)
        }
        # from / to
        if (including) {
           ind <- which(rd_et > index & rd_st < to) 
        } else {
           ind <- which(rd_st >= index & rd_et <= to) 
        }
    } else {
        # index only
        # parse index
        if (is.character(index)) {
            index <- parse_date_time3(index, tz = rd_tz)
        }
        # find intervals
        ind <- ibts::findI_st(as.numeric(index), as.numeric(rd_st), as.numeric(rd_et))
    }
    # return rawdat class
    structure(
        list(
            RawData = rawdat[['RawData']][, ind],
            Header = rawdat[['Header']][ind, ],
            DOASinfo = rawdat[['DOASinfo']]
            ),
        class = 'rawdat'
        )
}

# filter raw data by revolver position
filter_position <- function(rawdat, position) {
    # get indices
    ind <- which(rawdat[['Header']][['RevPos']] %in% position)
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
            RawData = rawdat[['RawData']][, ind],
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

