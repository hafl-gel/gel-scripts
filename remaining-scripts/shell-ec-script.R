#!/usr/bin/env Rscript

# set user library paths?

# get docopt
suppressMessages(
    if (!require(docopt)) {
        install.packages('docopt', verbose = FALSE, quiet = TRUE)
        library(docopt)
    }
)
# get jsonlite
suppressMessages(
    if (!require(jsonlite)) {
        install.packages('jsonlite', verbose = FALSE, quiet = TRUE)
        library(jsonlite)
    }
)

# load gel
suppressMessages(library(gel))

# build shell interface
doc <- paste0('Usage: shell-ec-script.R',
    ' SONIC_FILE [-dfsFvh] [-H HT_FILE] [-L LICOR_FILE] [-M MIRO_FILE]',
    ' [-S START] [-E END] [-Z ZSONIC] [-C ZCANOPY] [-D DEVNORTH]',
    ' [ARGS...]

SONIC_FILE                          (required) path to sonic file
-H HT_FILE --ht=HT_FILE             path to ht8700 file [default: NULL]
-L LICOR_FILE --licor=LICOR_FILE    path to licor file [default: NULL]
-M MIRO_FILE --miro=MIRO_FILE       path to miro file [default: NULL]
-S START --start=START              start datetime [default: -30mins]
-E END --end=END                    end datetime [default: now]
-Z ZSONIC --zsonic=ZSONIC           sonic height in m a.g.l [default: 2]
-C CANOPY --canopy=CANOPY           canopy height in m [default: 0.2]
-D DEVNORTH --north=DEVNORTH        sonic north deviation [default: 0]
-d --dyn                            provide reduced output (turb. + dyn. fluxes, lag, damping) [default: FALSE]
-f --fix                            provide reduced output (turb. + fix fluxes, lag, damping) [default: FALSE]
-s --subintervals                   should subintervals (default n = 5) be calculated [default: FALSE]
-F --foken                          apply Foken-Wichura classification (implies -s) [default: FALSE]
-v --verbose                        be verbose when running the main function [default: FALSE]
ARGS                                (optional) additional arguments passed to process_ec_fluxes() [default: NULL]
-h, --help                          show this help text

NOTE:
Some defaults are changed compared to the original R function, namely subintervals=FALSE, create_graphs=FALSE, create_dailygraphs=FALSE and as_ibts=FALSE
')

# process
opt <- docopt(doc)

# capture help
if (opt$help) {
    invisible()
} else {

    # verbose mode?
    verbose <- opt$verbose
    
    # get dyn & fix
    dyn <- opt$dyn
    fix <- opt$fix
    opt[c('dyn', 'fix')] <- NULL

    opt

    # check ARGS
    if (length(opt$ARGS) > 0) {
        ARGS <- opt$ARGS
        # fix strings
        for (a in seq_along(ARGS)) {
            if (grepl('/', ARGS[[a]])) {
                ARGS[[a]] <- sub('^([^=]+=)(.*)$', '\\1\"\\2\"', ARGS[[a]])
            }
        }
        # parse additional arguments
        ARGS <- eval(parse(text = paste0(
            'list(', paste(ARGS, collapse = ','), ')'
        )))
    } else {
        ARGS <- list()
    }

    # check licor
    if (opt$licor == 'NULL') {
        opt$licor <- NULL
    }
    # check ht
    if (opt$ht == 'NULL') {
        opt$ht <- NULL
    }
    # check miro
    if (opt$miro == 'NULL') {
        opt$miro <- NULL
    }

    # check time zone and set default UTC if unset
    if (!('tz_user' %in% names(ARGS))) {
        ARGS$tz_user <- 'UTC'
    }

    # check start
    nw <- now()
    if (opt$start != 'first') {
        # try parsing
        start <- try(parse_date_time3(opt$start, tz = ARGS$tz_user))
        if (inherits(start, 'try-error') || is.na(start)) {
            start <- nw + parse_time_diff(opt$start)
        }
        # check start
        stopifnot(start < nw)
        opt$start <- start
    }

    # check end
    if (opt$end %in% c('now', 'last')) {
        opt$end <- nw
    } else {
        # try parsing
        end <- try(parse_date_time3(opt$end, tz = ARGS$tz_user))
        if (inherits(end, 'try-error') || is.na(end)) {
            end <- nw + parse_time_diff(opt$end)
        }
        opt$end <- end
    }

    # check foken
    if (foken <- opt$foken) {
        # we need subintervals
        opt$subintervals <- TRUE
    }

    # get arguments
    opt <- opt[c('SONIC_FILE', 'ht', 'licor', 'miro', 'start', 'end',
        'zsonic', 'canopy', 'north', 'subintervals')]
    # fix numeric
    for (v in c('zsonic', 'canopy', 'north')) {
        opt[[v]] <- as.numeric(opt[[v]])
    }
    # fix names
    names(opt) <- c('sonic_directory', 'ht_directory', 'licor_directory',
        'miro_directory', 'start_time', 'end_time', 'z_ec', 'z_canopy', 'dev_north', 'subintervals')

    # get averaging time
    if ('avg_period' %in% names(ARGS)) {
        opt$avg_period <- ARGS$avg_period
        ARGS$avg_period <- NULL
    } else if (is.character(opt$start)) {
        opt$avg_period <- '30mins'
    } else {
        opt$avg_period <- as.numeric(opt$end - opt$start, units = 'secs')
    }


    # fix declination default
    if ('declination' %in% names(ARGS)) {
        opt$declination <- ARGS$declination
        ARGS$declination <- NULL
    } else {
        opt$declination <- 0
    }
    # fix create_graphs defaults
    if ('create_graphs' %in% names(ARGS)) {
        opt$create_graphs <- ARGS$create_graphs
        ARGS$create_graphs <- NULL
    } else {
        opt$create_graphs <- FALSE
    }
    # fix create_dailygraphs defaults
    if ('create_dailygraphs' %in% names(ARGS)) {
        opt$create_dailygraphs <- ARGS$create_dailygraphs
        ARGS$create_dailygraphs <- NULL
    } else {
        opt$create_dailygraphs <- FALSE
    }
    # fix as_ibts defaults
    if ('as_ibts' %in% names(ARGS)) {
        opt$as_ibts <- ARGS$as_ibts
        ARGS$as_ibts <- NULL
    } else {
        opt$as_ibts <- FALSE
    }

    # gather arguments
    all_args <- c(opt, ARGS)

    # fix missing lazy objects
    for (a in names(all_args)) {
        assign(a, all_args[[a]])
    }

    # fix missing defaults (why ???)
    frmls <- formals(process_ec_fluxes)
    subf <- frmls[!(names(frmls) %in% c(names(all_args), '...'))]
    for (a in names(subf)) {
        all_args[[a]] <- eval(subf[[a]])
    }

    # call flux processing function
    if (verbose) {
        out <- try(do.call(process_ec_fluxes, all_args), silent = TRUE)
    } else {
        suppressMessages(nirvana <- capture.output(out <- try(do.call(process_ec_fluxes, all_args), silent = TRUE)))
    }

    # check output
    if (inherits(out, 'try-error')) {
        attr(out, 'condition')$message
    } else {
        # run foken
        if (foken) {
            out <- foken_wichura(out)
        }
        cols <- character(0)
        if (dyn) {
            cols <- c('st', 'et', 'WD', 'Ustar', 'L', 'Zo', 'sUu', 'sVu', 
                'sWu', 'U_sonic', 'T_sonic', 'Ra', 'Rb_nh3',
                grep('^(phi|alpha|beta|w_bias)$', names(out), value = TRUE),
                grep('^avg_', names(out), value = TRUE),
                grep('^lag_dyn_', names(out), value = TRUE),
                grep('^flux_dyn_', names(out), value = TRUE),
                grep('^damping_dyn_pbreg_', names(out), value = TRUE)
                # grep('^damping_dyn_deming_', names(out), value = TRUE)
            )
        }
        if (fix) {
            cols <- c(cols,
                c('st', 'et', 'WD', 'Ustar', 'L', 'Zo', 'sUu', 'sVu', 
                    'sWu', 'U_sonic', 'T_sonic', 'Ra', 'Rb_nh3',
                    grep('^(phi|alpha|beta|w_bias)$', names(out), value = TRUE),
                    grep('^avg_', names(out), value = TRUE),
                    grep('^lag_fix_', names(out), value = TRUE),
                    grep('^flux_fix_', names(out), value = TRUE),
                    grep('^damping_fix_pbreg_', names(out), value = TRUE)
                    # grep('^damping_fix_deming_', names(out), value = TRUE)
                )
            )
        }
        if (length(cols)) {
            out <- out[, .SD, .SDcols = unique(cols)]
        }
        # add Linv
        out[, Linv := 1 / L]
        # convert to json
        jsonlite::toJSON(as.data.frame(out), data.frame = 'columns',
            digits = NA)
    }

}
