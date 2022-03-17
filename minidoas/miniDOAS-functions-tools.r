
require(data.table)
require(ibts)
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

#### extract raw spectra
get_specs <- function(folder, from, to = NULL, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.dark = FALSE, correct.linearity = FALSE) {
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
    if (correct.dark) {
        folder$RawData <- correct_dark(folder)
    }
    # correct linearity
    if (correct.linearity) {
        folder$RawData <- correct_linearity(folder)
    }
    structure(
        folder$RawData
        , RawData = folder
        , class = 'specs'
        , dark.corrected = correct.dark
        , linearity.corrected = correct.linearity
        )
}

#### extract single spectrum
single_specs <- function(folder, at, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.dark = FALSE, correct.linearity = FALSE) {
    # convert at to POSIXct
    at <- parse_date_time3(at, tz = tz)
    # get indices
    ind <- which(folder$Header[['et']] >= at & folder$Header[['st']] <= at)
    if (length(ind)) {
        folder$RawData <- folder$RawData[, ind, drop = FALSE]
        folder$Header <- folder$Header[ind, , drop = FALSE]
    } else {
        stop('No data at specified time available')
    }
    # correct for dark current
    if (correct.dark) {
        folder$RawData <- correct_dark(folder)
    }
    # correct linearity
    if (correct.linearity) {
        folder$RawData <- correct_linearity(folder)
    }
    structure(
        folder$RawData[[1]]
        , RawData = folder
        , class = 'single_spec'
        , dark.corrected = correct.dark
        , linearity.corrected = correct.linearity
        )
}

### get spec within wavelength range

### methods for single_spec class
print.single_spec <- function(x, lo = 200, hi = 230, ...) {
    rd <- attr(x, 'RawData')
    xi <- cut_wl(x, lo, hi)
    cat('***\nSingle spectrum:\n')
    cat('   recorded between', format(rd$Header[['st']]), 'and', format(rd$Header[['et']]), '\n')
    cat(sprintf('   I (%i to %i nm) min/avg/max: %1.0f/%1.0f/%1.0f\n***\n', lo, hi, min(xi, na.rm = TRUE), mean(xi, na.rm = TRUE), max(xi, na.rm = TRUE)))
}
plot.single_spec <- function(x, y, lo = 190, hi = 230, ylab = 'counts', xlab = 'nm', ...) {
    x <- cut_wl(x, lo, hi)
    y <- x
    class(y) <- 'numeric'
    x <- get_wl(x)
    plot(x, y, xlab = xlab, ylab = ylab, ...)
}
lines.single_spec <- function(x, lo = 190, hi = 230, ...) {
    x <- cut_wl(x, lo, hi)
    y <- x
    class(y) <- 'numeric'
    x <- get_wl(x)
    lines(x, y, ...)
}
points.single_spec <- function(x, lo = 190, hi = 230, ...) {
    x <- cut_wl(x, lo, hi)
    y <- x
    class(y) <- 'numeric'
    x <- get_wl(x)
    points(x, y, ...)
}


#### average raw data
avg_spec <- function(folder, from, to = NULL, tz = 'Etc/GMT-1', 
    doas = sub('.*(S[1-6]).*', '\\1', folder), Serial = NULL, 
    correct.dark = TRUE, correct.linearity = TRUE) {
    if(inherits(folder, 'rawdat')){
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
    } else {
        folder <- read_data(folder, from, to, tz, doas, Serial)
    }
    # correct for dark current
    if (correct.dark) {
        folder$RawData <- correct_dark(folder)
    }
    # correct linearity
    if (correct.linearity) {
        folder$RawData <- correct_linearity(folder)
    }
    structure(
        rowMeans(folder$RawData)
        , RawData = folder
        , class = 'avgdat'
        , dark.corrected = correct.dark
        , linearity.corrected = correct.linearity
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
plot.avgdat <- function(x, y = NULL, what = c('avg', 'resid'), xlim = c(190, 230), 
    ylim = NULL, type = 'l', main = NULL, ylab = NULL, ...) {
    # cut to xlim
    lo <- xlim[1] - diff(xlim) * 0.04
    hi <- xlim[2] + diff(xlim) * 0.04
    x <- cut_wl(x, lo = lo, hi = hi)
    # switch plot kind
    switch(pmatch(what[1], c('avg', 'resid'), nomatch = 3)
        # avg
        , yp <- x
        # resid
        , yp <- sweep(attr(x, 'RawData')$RawData, 1, x)
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
        plot(xp, x, type = 'n', xlim = xlim, ylim = ylim, xlab = '', ylab = ylab, main = main, ...)
        for (i in seq_len(ncol(yp))) {
            lines(xp, yp[, i], type = type, ...)
        }
    } else {
        plot(xp, yp, xlim = xlim, ylim = ylim, type = type, xlab = '', ylab = ylab, main = main, ...)
    }
}

#### lines method avgdat
lines.avgdat <- function(x, y = NULL, what = c('avg', 'resid'), ...) {
    # cut to xlim
    usr <- par('usr')[1:2]
    x <- cut_wl(x, lo = usr[1], hi = usr[2])
    # switch plot kind
    switch(pmatch(what[1], c('avg', 'resid'), nomatch = 3)
        # avg
        , yp <- x
        # resid
        , yp <- sweep(attr(x, 'RawData')$RawData, 1, x)
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
correct_dark <- function(x) {
    win <- getWindows(x$DOASinfo)
    if (is.data.frame(x$RawData)) {
        sweep(x$RawData, 2, colMeans(x$RawData[win$pixel_straylight, ]))
    } else {
        stop('Fix correct_dark for non-data.frames')
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
read_cal <- function(file, tz = 'Etc/GMT-1', Serial = NULL, spec = NULL,
    correct.dark = TRUE, correct.linearity = TRUE, lin_before_dark = FALSE) {
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
        # get entry
        spec_out <- file[[spec]]
        # get doas info
        DOASinfo <- getDOASinfo(file[[1]], timerange = spec_out$timerange)
        out <- list(
            data = data.table(wl = DOASinfo$Spectrometer$wavelength, cnt = as.numeric(spec_out[['dat.spec']])),
            Calinfo = list(
                info = data.table(var = 'timerange', val = paste0(spec_out$timerange, collapse = ' and ')),
                cuvette.gas = spec_out[['cuvette']][['cuvetteGas']],
                cuvette.conc = spec_out[['cuvette']][['cuvetteConc_mg']],
                cuvette.path = spec_out[['cuvette']][['cuvetteLength']]
                ),
            DOASinfo = DOASinfo
            )
    }
    # correct.dark
    if (correct.dark) {
        win <- getWindows(out$DOASinfo)
        dark <- out$data[, mean(cnt[win$pixel_straylight])]
        if (!lin_before_dark) {
            out$data[, cnt := cnt - dark]
        }
    }
    # correct linearity
    if (correct.linearity) {
        lin.coef <- out$DOASinfo$Spectrometer$'Linearity Coefficients'
        out$data[, cnt := cnt / linearity.func(cnt, lin.coef)]
    }
    # correct.dark afterwards
    if (correct.dark && lin_before_dark) {
        out$data[, cnt := cnt - dark]
    }
    structure(
        out,
        class = 'caldat'
        , dark.corrected = correct.dark
        , linearity.corrected = correct.linearity
        )
}

#### print caldat method
print.caldat <- function(x, ...){
    wl <- get_wl(x)
    cat('~~~~\n')
    cat('\t', x$DOASinfo$DOASmodel, '/', x$DOASinfo$Spectrometer$Serial, '- calibration spectrum\n')
    cat('\t cuvette gas:', x$Calinfo$cuvette.gas, '\n')
    cat('\t cuvette concentration (mg/m3):', x$Calinfo$cuvette.conc, '\n')
    cat('\t cuvette length (m):', x$Calinfo$cuvette.path, '\n')
    cat('\t', min(wl), 'to', max(wl), 'nm', sprintf('(%s pixel)\n', length(wl)))
    switch(as.character(nrow(x$Calinfo$info))
        , '1' = cat('\t recorded between', x$Calinfo$info[1, val], '\n')
        , '9' = cat('\t recorded between', x$Calinfo$info[14, val], '\n')
        , cat('\t recorded between', x$Calinfo$info[6, val], '\n')
    )
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
    fitwin = NULL, shift = NULL, correct.dark = TRUE, correct.linearity = TRUE,
    lin_before_dark = FALSE) {
    if (is.character(meas)) {
        # file path or chen
        if (tolower(meas) == 'chen') {
            # call chen2dc
            chen2dc(meas, ref, ftype, fstrength, fwin, fitwin, shift)
        } else {
            meas <- read_cal(meas, correct.dark = correct.dark, correct.linearity = correct.linearity, lin_before_dark = lin_before_dark)
        }
    } else if (inherits(meas, 'single_spec') || inherits(meas, 'avgdat')) {
        meas <- attr(meas, 'RawData')
    }
    if (is.character(ref)) ref <- read_cal(ref, correct.dark = correct.dark, correct.linearity = correct.linearity, lin_before_dark = lin_before_dark)
    # get counts
    m <- get_cnt(meas)
    r <- get_cnt(ref)
    # calc diffspec
    ds <- suppressWarnings(log(m / r))
    ds[!is.finite(ds)] <- NA_real_
    # calc dc
    win <- getWindows(meas$DOASinfo, filter.type = ftype, 
        filter.strength = fstrength, filter.window = fwin, fit.window = fitwin,
        tau.shift = shift)
    dc <- highpass.filter(ds, win)
    structure(
        list(
            wl = get_wl(meas)[win$pixel_filter],
            cnt = dc
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
cheng2dc <- function(dc, cheng = NULL, shift = FALSE) {
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
plot.caldat <- function(x, type = 'l', ...) {
    x$data[, plot(wl, cnt, type = type, ...)]
}
lines.caldat <- function(x, ...) {
    x$data[, lines(wl, cnt, ...)]
}

#### plot method for dc
plot.dc <- function(x, type = 'l', xlab = 'nm', ylab = 'doascurve', ...) {
    plot(x$wl, x$cnt, type = type, xlab = xlab, ylab = ylab, ...)
}
lines.dc <- function(x, ...) {
    lines(x$wl, x$cnt, ...)
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
    dcAula <- calc_dc(FileAulaNH3, FileAulaN2, correct.dark = c_dark, correct.linearity = c_lin, lin_before_dark = lin_first)
    dcOld <- calc_dc(FileOldNH3, FileOldN2, correct.dark = c_dark, correct.linearity = c_lin, lin_before_dark = lin_first)
    dsAula <- attr(dcAula, 'ds')[attr(dcAula, 'win')$pixel_filter]
    dsOld <- attr(dcOld, 'ds')[attr(dcOld, 'win')$pixel_filter]
    chAula <- cheng2dc(dcAula)
    chOld <- cheng2dc(dcOld)

    type <- 'l'
    plot(dcAula$wl, dcAula$cnt, type = type)
    # lines(dcAula$wl, chAula * (193e3 * 0.075), col = 'red', type = type)
    lines(dcOld$wl, dcOld$cnt, col = 'orange', type = type)

    # plot(dcAula$wl, dsAula, type = 'l', ylim = c(-0.8, 0))
    plot(calc_wl(attr(dcAula, 'meas'), attr(dcAula, 'win')$pixel_filter - 10), dsAula, type = 'l', ylim = c(-0.8, 0))
    lines(get_wl(nh3), cheng * 193e3 * 0.075, col = 'orange')

    dc1 <- calc_dc(FileAulaNH3, FileAulaN2, correct.dark = c_dark, correct.linearity = c_lin, lin_before_dark = FALSE)
    dc2 <- calc_dc(FileAulaNH3, FileAulaN2, correct.dark = c_dark, correct.linearity = c_lin, lin_before_dark = TRUE)
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

    Anh3 <- read_cal(FileAulaNH3, correct.dark = FALSE, correct.linearity = FALSE)
    Onh3 <- read_cal(FileOldNH3, correct.dark = FALSE, correct.linearity = FALSE)
    An2 <- read_cal(FileAulaN2, correct.dark = FALSE, correct.linearity = FALSE)
    On2 <- read_cal(FileOldN2, correct.dark = FALSE, correct.linearity = FALSE)
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

