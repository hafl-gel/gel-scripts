#!/usr/bin/env Rscript

# set user library paths?

# get docopt
suppressMessages(
    if (!require(docopt)) {
        install.packages('docopt', verbose = FALSE, quiet = TRUE)
        library(docopt)
    }
)

# load gel
suppressMessages(library(gel))

# build shell interface
doc <- paste0('Usage: shell-ec-script.R',
    ' SONIC_FILE [-H HT_FILE] [-L LICOR_FILE] [-M MIRO_FILE]',
    ' [-S START] [-E END] [-Z ZSONIC] [-C ZCANOPY] [-D DEVNORTH]',
    ' [ARGS...] [-h]

SONIC_FILE                          (required) path to sonic file
-H HT_FILE --ht=HT_FILE             path to ht8700 file [default: NULL]
-L LICOR_FILE --licor=LICOR_FILE    path to licor file [default: NULL]
-M MIRO_FILE --miro=MIRO_FILE       path to miro file [default: NULL]
-S START --start=START             start datetime [default: -30mins]
-E END --end=END                   end datetime [default: now]
-Z ZSONIC --zsonic=ZSONIC          sonic height in m a.g.l [default: 2]
-C CANOPY --canopy=CANOPY          canopy height in m [default: 0.2]
-D DEVNORTH --north=DEVNORTH       sonic north deviation [default: 0]
ARGS                                (optional) additional arguments passed to process_ec_fluxes() [default: NULL]
-h, --help                          show this help text

NOTE:
Some defaults are changed compared to the original R function, namely subintervals=FALSE, create_graphs=FALSE and as_ibts=FALSE
')

# process
opt <- docopt(doc)

# capture help
help <- opt$help

# check ARGS
if (length(opt$ARGS) > 0) {
    # parse additional arguments
    ARGS <- eval(parse(text = paste0(
        'list(', paste(opt$ARGS, collapse = ','), ')'
    )))
} else {
    ARGS <- NULL
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

# get arguments
opt <- opt[c('SONIC_FILE', 'ht', 'licor', 'miro', 'start', 'end',
    'zsonic', 'canopy', 'north')]
# fix numeric
for (v in c('zsonic', 'canopy', 'north')) {
    opt[[v]] <- as.numeric(opt[[v]])
}
# fix names
names(opt) <- c('sonic_directory', 'ht_directory', 'licor_directory',
    'miro_directory', 'start_time', 'end_time', 'z_ec', 'z_canopy', 'dev_north')

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
# fix subintervals defaults
if ('subintervals' %in% names(ARGS)) {
    opt$subintervals <- ARGS$subintervals
    ARGS$subintervals <- NULL
} else {
    opt$subintervals <- FALSE
}
# fix create_graphs defaults
if ('create_graphs' %in% names(ARGS)) {
    opt$create_graphs <- ARGS$create_graphs
    ARGS$create_graphs <- NULL
} else {
    opt$create_graphs <- FALSE
}
# fix as_ibts defaults
if ('as_ibts' %in% names(ARGS)) {
    opt$as_ibts <- ARGS$as_ibts
    ARGS$as_ibts <- NULL
} else {
    opt$as_ibts <- FALSE
}

if (help) {
    invisible()
} else {
    # call flux processing function
    out <- try(do.call(process_ec_fluxes, c(opt, ARGS)), silent = TRUE)
    if (inherits(out, 'try-error')) {
        attr(out, 'condition')$message
    } else {
        out
    }
}
