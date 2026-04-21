#!/usr/bin/env Rscript

# set user library paths?

# get docopt
suppressMessages(
    if (!require(docopt)) {
        install.packages('docopt', verbose = FALSE, quiet = TRUE)
        library(docopt)
    }
)

suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN

doc <- "Usage: shell-ec-script.R [-S SONIC] [-h]

-S --sonic SONIC    sonic file or path to sonic file (required)
-h --help           show this help text"

opt <- docopt(doc)

print(opt)

# if (opt$deps == "TRUE" || opt$deps == "FALSE") {
#     opt$deps <- as.logical(opt$deps)
# } else if (opt$deps == "NA") {
#     opt$deps <- NA
# }

# if (opt$repos == "NULL") {
#     opt$repos <- NULL
# }

# # check if PACKAGES is a path to a file/directory
# if (length(opt$PACKAGES) == 1 && file.exists(opt$PACKAGES)) {
#     old_path <- opt$PACKAGES
#     # try to find package top directory (max two levels up)
#     i <- 0
#     while (!(is_top <- 'NAMESPACE' %in% dir(opt$PACKAGES)) && i <= 2) {
#         i <- i + 1
#         opt$PACKAGES <- dirname(opt$PACKAGES)
#     }
#     opt$repos <- NULL
#     opt$deps <- NA
# } else {
#     is_top <- TRUE
# }

# if (!is_top) {
#     cat('cannot find package top directory from given path!\n')
# } else if (opt$error) {
#     withCallingHandlers(
#         install.packages(pkgs  = opt$PACKAGES,
#                          lib   = opt$libloc,
#                          repos = opt$repos,
#                          clean = TRUE,
#                          Ncpus = opt$ncpus,
#                          dependencies=opt$deps),
#         warning = stop)

# } else {
#     install.packages(pkgs  = opt$PACKAGES,
#                      lib   = opt$libloc,
#                      repos = opt$repos,
#                      clean = TRUE,
#                      ncpus = opt$ncpus,
#                      dependencies=opt$deps)
# }
