
# formulas & parameters:
# https://geodesy.geo.admin.ch/declination/doc/geomagnetic-survey-ch-en.pdf

require(data.table)

# functions
f_t <- function(d_yr, ps = paras) {
    (d_yr + ps['beta'] * d_yr ^ 2) / (1 + ps['gamma'] * d_yr)
}
D_xyt <- function(d_x, d_y, d_yr, ps = paras) {
    ps['A'] + ps['B'] * d_x + ps['C'] * d_y + 
        ps['alpha'] * (1 + ps['b'] * d_x + ps['c'] * d_y) * 
        f_t(d_yr, ps)
}

declination <- function(X, Y, Date = Sys.time(), 
    type = c('declination', 'inclination', 'field')[1],
    folder = './declination') {
    # convert Date to POSIXt
    if (is.character(Date)) {
        require(ibts)
        Date <- parse_date_time3(Date)
    }
    # convert to lt
    Datelt <- as.POSIXlt(Date)
    # get year
    year <- Datelt$year + 1900
    # get fraction
    Datelt$year <- Datelt$year + 1
    Tyear <- year + Datelt$yday / as.numeric(Datelt - Date, units = 'days')
    # fix X & Y
    if (X > 10000) {
        X <- X / 1000
    }
    if (Y > 10000) {
        Y <- Y / 1000
    }
    # check year
    file_paras <- file.path(folder, 'parameters')
    if (!file.exists(file_paras)) {
        stop('File "', file_paras, '" does not exist.')
    }
    # read parameters 
    paras_raw <- fread(file_paras, fill = TRUE)
    setnames(paras_raw, c('coef', 'declination', 'inclination', 'field'))
    # get parameters
    type <- c('declination', 'inclination', 'field')[pmatch(type[1], c('declination', 'inclination', 'field'))]
    if (is.na(type)) stop('Argument type should be one of "declination", "inclination" or "field"!')
    paras <- paras_raw[grepl('[^0]$', coef), setNames(get(type), coef)]
    # get zeros
    zeros <- paras_raw[grepl('0$', coef), setNames(declination, coef)]
    # calculate result
    D_xyt((X - zeros['X0']), (Y - zeros['Y0']), (Tyear - zeros['T0']), paras)
}

# declination(200, 600, '21.01.2021')

