
library(lubridate)
library(data.table)
library(ibts)

####### functions

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

check_ida <- function(x) {
    nms0 <- names(x[, id:et])
    nms1 <- names(x)[!(names(x) %in% nms0)]
    x[, 
        lapply(mget(nms1), function(x) {
            sum(!is.na(x)) / .N * 100
            }) , by = id]
}

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
            function(x) fread(text = dat[iEmpty[x]:iEmpty[x + 1]], na.strings = '-', colClasses = c(time = 'character')))
        , use.names = TRUE, fill = TRUE)

    # fix time
    nc <- Data[1, nchar(time)]
    if (nc == 10L) {
        # hourly data
        Data[, et := fast_strptime(time, '%Y%m%d%H', tz = 'UTC', lt = FALSE)][,
            st := et - 3600]
    } else if (nc == 12L) {
        # 10-minute data
        Data[, st := fast_strptime(time, '%Y%m%d%H%M', tz = 'UTC', lt = FALSE)][,
            et := st + 600]
    } else {
        # other granularity
        stop('granularity not yet implemented!')
    }

    # merge Data and Stats
    Stat2 <- Stats[, .(stn, name = Name, 
        id = paste0('id', seq_len(.N)),
        ch.x = as.numeric(sub('([0-9]{6})[/].*', '\\1', Koordinaten..km.)),
        ch.y = as.numeric(sub('[0-9]{6}[/]([0-9]{6})', '\\1', Koordinaten..km.)),
        m.asl = Höhe.ü..M...m.
        )]
    out <- merge(Data, Stat2, all = TRUE)[, time := NULL]

    # set key column
    setkey(out, id)

    # prepare output columns
    nms0 <- c('id', 'stn', 'name', 'ch.x', 'ch.y', 'm.asl', 'st', 'et')
    nms1 <- names(out)[!(names(out) %in% nms0)]

    # change column order
    setcolorder(out, c(nms0, nms1))

    # add Parameter attribute
    setattr(out, 'parameters', as.data.frame(Pars))

    # return
    out
}
