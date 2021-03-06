
# transition from require/library to box::use
if ('box' %in% .packages(TRUE) && length(box::name()) != 0) {
    # although considered bad praxis, we want all
    # data.table functions available in .GlobalEnv:
    evalq(box::use(data.table[...]), .GlobalEnv)
    # load ibts/data.table package for module
    cat(paste0('*** ibts namespace not available ***\n-> Consider adding:\n\n',
        'box::use(ibts[...])\n\n',
        'to your script...\n*** ~~~~~~~~~~~~~~~~~~ ***\n'))
    # import necessary functions
    box::use(
        data.table[...], 
        ibts[as.ibts], 
        lubridate[fast_strptime],
        utils[unzip],
        graphics[par]
        )
} else {
    library(lubridate)
    library(data.table)
    library(ibts)
}

#' Print IDA data summary 
#'
#' This function prints a short summary of a idaweb data set.
#'
#' @param x A `data.table` containing idaweb data
#' @return A `data.table` with a data summary (one row per id)
#' @export
smry_ida <- function(x) {
    nms0 <- names(x[, id:et])
    nms1 <- names(x)[!(names(x) %in% nms0)]
    x[, c(
            list(
                st.min = st[1],
                et.max = et[.N],
                N = .N
            ),
            lapply(mget(nms1), function(x) {
                vls <- sum(!is.na(x)) / .N * 100
                sprintf('%s%% (%s, %s)',
                    as.character(round(vls, 1)),
                    as.character(ifelse(vls > 0, round(min(x, na.rm = TRUE), 2), NA)),
                    as.character(ifelse(vls > 0, round(max(x, na.rm = TRUE), 2), NA))
                )
            })
        ), by = .(id, stn, name)]
}

#' Check IDA data coverage
#'
#' This function prints the percentage of available data in
#' an idaweb data set.
#'
#' @param x A `data.table` containing idaweb data
#' @return A `data.table` with data coverage summary (one row per id)
#' @export
check_ida <- function(x) {
    nms0 <- names(x[, id:et])
    nms1 <- names(x)[!(names(x) %in% nms0)]
    x[, 
        lapply(mget(nms1), function(x) {
            sum(!is.na(x)) / .N * 100
            }) , by = id]
}

#' Plot IDA data time series
#'
#' This function plots the columns of an idaweb data set as 
#' time series plots.
#'
#' @param x A `data.table` containing idaweb data
#' @export
plot_ida <- function(x) {
    ch <- check_ida(x)
    nms1 <- names(ch[, -1])
    nt <- sum(ch[, mget(nms1)] > 0)
    nr <- floor(sqrt(nt))
    nc <- ceiling(nt / nr)
    pars <- attr(x, 'parameters')
    rownames(pars) <- pars[, 1]
    par(mfrow = c(nr, nc))
    x[, {
        y <- as.ibts(mget(nms1), st = st, et = et, tz = 'UTC')
        for (nm in nms1) {
            if (ch[id %in% .BY$id, get(nm) > 0]) plot(y[, nm], 
                main = paste(.BY$stn, pars[nm, 3], sep = ' - '),
                ylab = paste0(nm, ' (', pars[nm, 2], ')'),
                xlab_fmt = '%y-%m-%d')
        }
        }, by = .(id, stn, name)]
    invisible(NULL)
}

#' Read IDA data from an idaweb zip file 
#'
#' This function reads data from an idaweb zip file.
#'
#' @param File character. The path to an idaweb zip file.
#' @return A `data.table` containing idaweb data
#' @export
read_ida <- function(File) {

    # get zipped file names
    zip_names <- unzip(File, list = TRUE)$Name

    # read zipped legend
    leg <- readLines(
        con <- unz(File, grep('legend', zip_names, value = TRUE))
        ); close(con)
    # convert encoding
    leg <- iconv(leg, 'latin1')

    # get Stationen and Parameter
    iEmpty <- grep('^$', leg)
    iSf <- grep('Stationen', leg) + 2
    iPf <- grep('^Parameter$', leg) + 2
    iSt <- iEmpty[which(iEmpty > iSf)[1]] - 1
    iPt <- iEmpty[which(iEmpty > iPf)[1]] - 1

    # read Parameter
    Pars <- fread(text = gsub('[ ]{2, }', ',', leg[iPf:iPt]), sep = ',', header = TRUE)
    setnames(Pars, 'V1', 'Variable')

    # read Stationen
    Stats <- unique(
        fread(text = gsub('([ ]{2, }|[[]km[]])', '\\1,', leg[iSf:iSt]), sep = ',',
            check.names = TRUE)[, c('Datenquelle', 'Parameter') := NULL]
        )

    # read zipped data
    dat <- readLines(
        con <- unz(File, grep('data', zip_names, value = TRUE))
        ); close(con)
    dat <- iconv(dat, 'latin1')

    # get empty lines to see different data chunks
    iEmpty <- unique(c(grep('^$', dat), length(dat)))

    # read different chunks
    Data <- rbindlist(
        lapply(seq_len(length(iEmpty) - 1), 
            function(x) {
                out <- fread(text = dat[iEmpty[x]:iEmpty[x + 1]], na.strings = '-', colClasses = c(time = 'character'))
                # fix time
                nc <- out[1, nchar(time)]
                if (nc == 10L) {
                    # hourly data
                    out[, granularity := '1hours']
                    out[, et := fast_strptime(time, '%Y%m%d%H', tz = 'UTC', lt = FALSE)][,
                        st := et - 3600]
                } else if (nc == 12L) {
                    # 10-minute data
                    out[, granularity := '10mins']
                    out[, et := fast_strptime(time, '%Y%m%d%H%M', tz = 'UTC', lt = FALSE)][,
                        st := et - 600]
                } else {
                    # other granularity
                    stop('granularity not yet implemented!')
                }

            })
        , use.names = TRUE, fill = TRUE)

    # merge Data and Stats
    Stat2 <- Stats[, .(stn, name = Name, 
        id = paste0('id', seq_len(.N)),
        ch.x = as.numeric(sub('([0-9]{6})[/].*', '\\1', Koordinaten..km.)),
        ch.y = as.numeric(sub('[0-9]{6}[/]([0-9]{6})', '\\1', Koordinaten..km.)),
        # m.asl = H??he.??..M...m.
        m.asl = get('H\u00f6he.\u00fc..M...m.')
        )]
    out <- merge(Data, Stat2, all = TRUE)[, time := NULL][]

    # set key column
    setkey(out, id)

    # prepare output columns
    nms0 <- c('id', 'stn', 'name', 'ch.x', 'ch.y', 'm.asl', 'granularity', 'st', 'et')
    nms1 <- names(out)[!(names(out) %in% nms0)]

    # change column order
    setcolorder(out, c(nms0, nms1))

    # add Parameter attribute
    setattr(out, 'parameters', as.data.frame(Pars))

    # return
    out
}

