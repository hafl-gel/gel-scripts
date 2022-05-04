
# formulas & parameters:
# https://geodesy.geo.admin.ch/declination/doc/geomagnetic-survey-ch-en.pdf

require(data.table)

# calculate declination
declination <- function(X, Y, Date = Sys.time(), 
    type = c('declination', 'inclination', 'field')[1],
    folder = './declination', correct.anomaly = type == 'declination',
    correct.merid = type == 'declination') {
    # convert Date to POSIXt
    if (is.character(Date)) {
        require(ibts)
        Date <- parse_date_time3(Date)
    }
    # convert to lt
    Datelt <- as.POSIXlt(Date)
    # get year as number
    year <- Datelt$year + 1900L
    # get fraction
    Datelt$year <- Datelt$year + 1
    Tyear <- year + Datelt$yday / as.numeric(Datelt - Date, units = 'days')
    # fix X & Y
    if (X < 10000) {
        X <- X * 1000
    }
    if (Y < 10000) {
        Y <- Y * 1000
    }
    # check year
    file_paras <- file.path(folder, 'parameters')
    if (!file.exists(file_paras)) {
        stop('File "', file_paras, '" does not exist.')
    }
    # read parameters 
    paras_raw <- fread(file_paras, fill = TRUE)
    setnames(paras_raw, c('lab', names(paras_raw)[-length(paras_raw)]))
    # fix type
    i_type <- pmatch(type[1], c('declination', 'inclination', 'field'))
    if (is.na(i_type)) stop('Argument type should be one of "declination", "inclination" or "field"!')
    # get parameters for correct year
    paras <- paras_raw[, {
        ind_paras <- which(lab == as.character(min(2021, year))) + i_type
        if (!length(ind_paras)) stop('parameters for year ', floor(Tyear), ' missing!')
        setNames(as.numeric(unlist(.SD[ind_paras, -1])), names(paras_raw)[-1])
    }]
    # get zeros
    zeros <- c(X0 = 200, Y0 = 600, T0 = 2003.5)
    # calc d_x etc.
    d_loc <- c(X = X / 1000, Y = Y / 1000, T = Tyear) - zeros
    # calculate result
    out <- D_xyt(d_loc[['X']], d_loc[['Y']], d_loc[['T']], paras)
    # correct meridian convergence?
    if (correct.merid) {
        # get meridian convergence (see excel file)
        out <- out - calc_merid(d_loc[['X']], d_loc[['Y']])
    }
    # correct anomaly?
    if (correct.anomaly) {
        # read file
        anomaly_file <- file.path(folder, 'anomalies.grd')
        # anomalies header
        header <- fread(anomaly_file, nrows = 6)[, setNames(V2, V1)]
        # anomalies values
        values <- as.matrix(fread(anomaly_file, skip = 6, na.strings = as.character(header['nodata_value'])))
        # get raw row/col
        row <- header[['nrows']] - ((X - header[['yllcorner']]) / header[['cellsize']])
        col <- (Y - header[['xllcorner']]) / header[['cellsize']] + 1
        # get rows/cols
        rows <- c(floor(row), ceiling(row))
        cols <- c(floor(col), ceiling(col))
        # # print
        # print(values[rows, cols])
        # get anomaly
        out + mean(values[rows, cols]) / 60
    } else {
        out
    }
}

# helper functions
f_t <- function(d_yr, ps = paras) {
    (d_yr + ps[['beta']] * d_yr ^ 2) / (1 + ps[['gamma']] * d_yr)
}
D_xyt <- function(d_x, d_y, d_yr, ps = paras) {
    (
        ps[['A']] + 
        ps[['B']] * d_x + 
        ps[['C']] * d_y
    ) + 
    (
        ps[['alpha']] * 
            (1 + ps[['b']] * d_x + ps[['c']] * d_y) 
    ) * f_t(d_yr, ps)
}
calc_merid <- function(x, y, ps = c(0.0096, 0.0000016)) {
    ps[1] * y + ps[2] * x * y
}

# # checks with https://www.swisstopo.admin.ch/en/maps-data-online/calculation-services/deklination.html
# declination(200, 600, '21.01.2021')
# # should read 2.66425 (reads: 2.084385)
# declination(200, 600, '21.01.2021', correct.anomaly = FALSE)
# # should read 2.68508 (reads: 2.684385)
# declination(200, 600, '21.03.2022')
# # should read 2.86293 (reads: 2.283059)
# declination(200, 600, '21.03.2022', correct.anomaly = FALSE)
# # should read 2.88377 (reads: 2.883059)

# declination(250, 600, '21.03.2022')
# declination(250, 600, '21.03.2022', correct.anomaly = FALSE)


# declination(100, 620, '21.01.2021') # soll: 1.76; ist: 2.13

# declination(100, 620, '21.01.2021') # soll: 1.76; ist: 2.13
# declination(90, 670, '21.01.2021')  # soll: 2.86; ist: 3.54

# declination(100, 620, '21.01.2021', correct.anomaly = FALSE) # soll: 2.62; ist: 2.81
# declination(90, 670, '21.01.2021', correct.anomaly = FALSE)  # soll: 2.32; ist: 2.98



