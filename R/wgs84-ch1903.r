## ~~~~ new functions using EPSG ~~~~ ##

## input data format:
# 1.) blsmodelr -> Sources and Sensors
# 2.) matrix/data.frame/list with x & y coords and other entries

## functions:
# coord_transf
#   -> crs: map, local, crs
#   -> local or map: check origin not missing!!!
#   -> crs_to/crs_from == Map?
# local_to_...
# ch1903_to_...
# wgs84_to_...
# map_to_...

# check EPSG codes on https://epsg.io
# wgs84: 4326
# ch1903/LV03: 21781
# ch1903/LV95: 2056

## ~~~~~~~~~~~~~~~ functions

fix_crs <- function(crs = NULL,
    new_origin_at = NULL) {
    if (is.null(crs)) {
        return(NULL)
    }
    if (is.character(crs) && tolower(crs[1]) %in% 
        c('wgs84', 'ch1903/lv03', 'ch03', 'lv03', 'ch1903/lv95', 'lv95', 'ch95')) {
        crs <- switch(tolower(crs[1])
            , 'wgs84' = 'EPSG:4326'
            , 'ch1903/lv03' = 
                , 'ch03' = 
                , 'lv03' = 'EPSG:21781'
            , 'ch1903/lv95' = 
                , 'ch95' = 
                , 'lv95' = 'EPSG:2056'
            )
    }
    crs <- sf::st_crs(crs)
    if (!is.null(new_origin_at) && length(new_origin_at) > 0) {
        base_origin <- get_origin(crs)
        crs <- modify_crs(crs, new_origin = base_origin - as.numeric(new_origin_at))
    }
    crs
}
get_crs <- function(obj, fix = FALSE) {
    out <- list(
        crs_gel = attr(obj, 'crs_gel'),
        origin = attr(obj, 'origin')
        )
    if (fix) {
        out <- fix_crs(out$crs_gel, out$origin)
    }
    out
}
set_crs <- function(obj, crs = NULL,
    new_origin_at = NULL) {
    validate_crs <- fix_crs(crs, new_origin_at)
    if (inherits(obj, 'data.table')) {
        setattr(obj, 'crs_gel', crs)
        setattr(obj, 'origin', new_origin_at)
    } else {
        attr(obj, 'crs_gel') <- crs
        attr(obj, 'origin') <- new_origin_at
    }
    obj
}



# ## ~~~ general helpers
get_names <- function(obj) UseMethod('get_names')
get_names.default <- function(obj) names(obj)
get_names.matrix <- function(obj) colnames(obj)
# get_length <- function(obj) UseMethod('get_length')
# get_length.default <- function(obj) length(obj)
# get_length.matrix <- function(obj) ncol(obj)
# names_exist <- function(obj, col) UseMethod('names_exist')
# setMethod('names_exist', 
#     signature(obj = 'ANY', col = 'numeric'),
#     function(obj, col) {
#         col <= get_length(obj)
#     })
# setMethod('names_exist', 
#     signature(obj = 'ANY', col = 'character'),
#     function(obj, col) {
#         col %in% get_names(obj)
#     })
# setMethod('names_exist', 
#     signature(obj = 'ANY', col = 'NULL'),
#     function(obj, col) {
#         TRUE
#     })

get_col_classes <- function(obj) UseMethod('get_col_classes')
setMethod('get_col_classes',
    signature(obj = 'ANY'),
    function(obj) {
        unlist(lapply(obj, \(x) class(x)[[1]]))
    })
setMethod('get_col_classes',
    signature(obj = 'matrix'),
    function(obj) {
        setNames(
            rep(class(obj[[1]])[[1]], ncol(obj)),
            colnames(obj)
            )
    })
check_numeric <- function(obj, name) UseMethod('check_numeric')
setMethod('check_numeric',
    signature(obj = 'ANY', name = 'ANY'),
    function(obj, name) {
        cc <- get_col_classes(obj)[name] 
        setNames(
            cc %in% c('integer', 'numeric'),
            names(cc)
            )
    })
setMethod('check_numeric',
    signature(obj = 'ANY', name = 'missing'),
    function(obj, name) {
        setNames(
            get_col_classes(obj) %in% c('integer', 'numeric'),
            get_names(obj)
            )
    })

guess_coord_x <- function(obj, coord_system = NULL) {
    isn <- check_numeric(obj)
    if (sum(isn) < 2) {
        stop('coordinates must be provided as numeric values!')
    }
    nms <- get_names(obj)
    # check null
    if (is.null(nms)) {
        if (is.character(obj)) {
            nms <- obj
        } else {
            nms <- NULL
        }
    }
    nnms <- nms[isn]
    wisn <- which(isn)
    # check coord_system
    if (!is.null(coord_system)) {
        if (is.list(coord_system)) {
            coord_system <- switch(sub('.*,(\\d+)[]]*$', '\\1', coord_system$wkt)
                , '4326' = 'wgs84'
                , '21781' = 'lv03'
                , '2056' = 'lv95'
            )
        } else {
            coord_system <- sub('ch1903/', '', match.arg(tolower(coord_system), 
                    c('wgs84', 'lv03', 'lv95', 'ch1903/lv03', 'ch1903/lv95')), fixed = TRUE)
        }
    }
    # check x
    x_nm <- integer(0)
    if (is.null(coord_system) || coord_system %in% c('lv03', 'lv95')) {
        x_nm <- grep('^((x|X)(.*(c|C)oord.*)?)$|^((c|C)oord.*(x|X))$', nnms)
    }
    if ((is.null(coord_system) && length(x_nm) == 0) || (!is.null(coord_system) && coord_system == 'wgs84')) {
        # check lon
        x_nm <- grep('^(l|L)o?n((g|gitude).*)?$', nnms)
    }
    # check length
    if (length(x_nm) == 0) {
        # find numeric
        if (is.matrix(obj)) {
            obj <- as.data.frame(obj)
        }
        if (is.null(coord_system)) {
            x_nm <- which(sapply(wisn, \(x) {
                # wgs84
                all(
                    (obj[[x]] > 4 & obj[[x]] < 12) |
                    (obj[[x]] > 2400000 & obj[[x]] < 2900000) |
                    (obj[[x]] > 400000 & obj[[x]] < 900000)
                )
            }))
        } else if (coord_system == 'wgs84') {
            x_nm <- which(sapply(wisn, \(x) {
                # wgs84
                all(obj[[x]] > 4 & obj[[x]] < 12)
            }))
        } else if (coord_system == 'lv95') {
            x_nm <- which(sapply(wisn, \(x) {
                # lv95
                all(obj[[x]] > 2400000 & obj[[x]] < 2900000)
            }))
        } else {
            x_nm <- which(sapply(wisn, \(x) {
                # lv03
                all(obj[[x]] > 400000 & obj[[x]] < 900000)
            }))
        }
    }
    switch(as.character(length(x_nm))
        # 0
        , '0' = {
            stop('x-coordinates could not be guessed successfully!')
        }
        # 1
        , '1' = {}
        # >1
        , {
            warning('more than one potential x coordinate found. returning first finding')
        }
    )
    wisn[x_nm[1]]
}

guess_coord_y <- function(obj, coord_system = NULL) {
    isn <- check_numeric(obj)
    if (sum(isn) < 2) {
        stop('coordinates must be provided as numeric values!')
    }
    nms <- get_names(obj)
    # check null
    if (is.null(nms)) {
        if (is.character(obj)) {
            nms <- obj
        } else {
            nms <- NULL
        }
    }
    nnms <- nms[isn]
    wisn <- which(isn)
    # check coord_system
    if (!is.null(coord_system)) {
        if (is.list(coord_system)) {
            coord_system <- switch(sub('.*,(\\d+)[]]*$', '\\1', coord_system$wkt)
                , '4326' = 'wgs84'
                , '21781' = 'lv03'
                , '2056' = 'lv95'
            )
        } else {
            coord_system <- sub('ch1903/', '', match.arg(tolower(coord_system), 
                    c('wgs84', 'lv03', 'lv95', 'ch1903/lv03', 'ch1903/lv95')), fixed = TRUE)
        }
    }
    # check y
    y_nm <- integer(0)
    if (is.null(coord_system) || coord_system %in% c('lv03', 'lv95')) {
        y_nm <- grep('^((y|Y)(.*(c|C)oord.*)?)$|^((c|C)oord.*(y|Y))$', nnms)
    }
    if ((is.null(coord_system) && length(y_nm) == 0) || (!is.null(coord_system) && coord_system == 'wgs84')) {
        # check lat
        y_nm <- grep('^(l|L)at(itude.*)?$', nnms)
    }
    # check length
    if (length(y_nm) == 0) {
        # find numeric
        if (is.matrix(obj)) {
            obj <- as.data.frame(obj)
        }
        if (is.null(coord_system)) {
            y_nm <- which(sapply(wisn, \(x) {
                # wgs84
                all(
                    (obj[[x]] > 44 & obj[[x]] < 49) |
                    (obj[[x]] > 1050000 & obj[[x]] < 1400000) |
                    (obj[[x]] > 50000 & obj[[x]] < 400000)
                )
            }))
        } else if (coord_system == 'wgs84') {
            y_nm <- which(sapply(wisn, \(x) {
                # wgs84
                all(obj[[x]] > 44 & obj[[x]] < 49)
            }))
        } else if (coord_system == 'lv95') {
            y_nm <- which(sapply(wisn, \(x) {
                # lv95
                all(obj[[x]] > 1050000 & obj[[x]] < 1400000)
            }))
        } else {
            y_nm <- which(sapply(wisn, \(x) {
                # lv03
                all(obj[[x]] > 50000 & obj[[x]] < 400000)
            }))
        }
    }
    switch(as.character(length(y_nm))
        # 0
        , '0' = {
            stop('y-coordinates could not be guessed successfully!')
        }
        # 1
        , '1' = {}
        # >1
        , {
            warning('more than one potential y coordinate found. returning first finding')
        }
    )
    wisn[y_nm[1]]
}

# TODO: add option to indicate coord system
guess_coords <- function(obj, coord_system = NULL) {
    c(
        guess_coord_x(obj, coord_system = coord_system)[1], 
        guess_coord_y(obj, coord_system = coord_system)[1]
    )
}

parse_wkt <- function(crs) {
    # crs <- sf::st_crs('EPSG:2056')
    wkt_a <- gsub('\\s*\\n\\s*', '', crs$wkt)
    # replace ,[] inside string with {1}, {2}[ {3}]
    wkt_a2 <- gsub('([,]|\\[|\\])(?=([^"]*["][^"]*["])*[^"]*$)', '{\\1}', wkt_a, perl = TRUE)
    wkt_b <- unlist(strsplit(wkt_a2, split = '(?<=[{](\\[|,)[}])', perl = TRUE))
    wkt_c <- sub('^([^0-9^"]+)\\{,}', 'as.name("\\1"),', wkt_b)
    wkt_d <- sub('{[}', '=list(', wkt_c, fixed = TRUE)
    wkt_e <- gsub('{]}', ')', wkt_d, fixed = TRUE)
    wkt_f <- gsub('{,}', ',', wkt_e, fixed = TRUE)
    eval(parse(text = paste0('list(', paste(wkt_f, collapse = ''), ')')))
}
deparse_wkt <- function(x) {
    dep_crs <- deparse(x)
    dep2 <- gsub('(\\))(?=([^"]*["][^"]*["])*[^"]*$)', '{\\1}', dep_crs, perl = TRUE)
    paste(sub('(^list\\(|\\]$)', '', gsub('{)}', ']', gsub(' = list(', '[', dep2, fixed = TRUE),
            fixed = TRUE)), collapse = '')
}

traverse_list <- function(x) .traverse_list(x)
.traverse_list <- function(x, i = numeric(0)) {
    if (is.list(x)) {
        out <- mapply(.traverse_list, x, i = lapply(seq_along(x), function(a) c(i, a)), SIMPLIFY = FALSE)
        do.call(rbind, out)
    } else {
        cbind(list(i), x)
    }
}

# find EPSG code:
.get_indices <- function(t_wkt, code) {
    i_epsg <- which(t_wkt[, 2] %in% 'EPSG')
    ind <- which(t_wkt[i_epsg + 1, 2] %in% code)
    if (length(ind) > 1) {
        stop('several codes ', paste(t_wkt[i_epsg[ind] + 1, 2], collapse = ' + '),
            ' (n != 1) were found in wkt string')
    }
    index <- t_wkt[[i_epsg[ind] + 1]]
    c(index[-length(index) + (0:1)], 2)
    # y[[index]]
}

modify_crs <- function(base_crs, new_origin = NULL, new_angle = NULL) {
    if (missing(base_crs)) {
        stop('base crs is missing')
    }
    if (is.character(base_crs)) {
        base_crs <- try(sf::st_crs(base_crs), silent = TRUE)
    }
    if (!inherits(base_crs, 'crs')) {
        stop('base crs argument is not a valid crs')
    }
    if (is.null(new_origin) && is.null(new_angle)) {
        return(base_crs)
    }
    # parse wkt
    p_wkt <- parse_wkt(base_crs)
    # traverse wkt
    i_wkt <- traverse_list(p_wkt)
    # false easting/northing -> check epsg code in Table F.3 (https://docs.ogc.org/is/18-010r7/18-010r7.html#106)
    tabf.3 <- list(
        false_easting = c(8806, 8816, 8826),
        false_northing = c(8807, 8817, 8827),
        angle_northing = 8814
        )
    if (!is.null(new_origin)) {
        if (!is.numeric(new_origin)) {
            stop('argument new_origin must be numeric')
        }
        if (length(new_origin) != 2) {
            stop('argument new_origin must be of length 2')
        }
        # get false easting (x-axis)
        i_x <- .get_indices(i_wkt, tabf.3[['false_easting']])
        # get false northing (y-axis)
        i_y <- .get_indices(i_wkt, tabf.3[['false_northing']])
        # replace values
        p_wkt[[i_x]] <- new_origin[1]
        p_wkt[[i_y]] <- new_origin[2]
    }
    if (!is.null(new_angle)) {
        if (!is.numeric(new_angle)) {
            stop('argument new_angle must be numeric')
        }
        if (length(new_angle) != 1) {
            stop('argument new_angle must be of length 1')
        }
        # get angle northing
        i_n <- .get_indices(i_wkt, tabf.3[['angle_northing']])
        # replace value
        p_wkt[[i_n]] <- new_angle
    }
    # replace name
    p_wkt[[1]][[1]] <- "User"
    sf::st_crs(deparse_wkt(p_wkt))
}
get_origin <- function(base_crs) {
    if (missing(base_crs)) {
        stop('base crs is missing')
    }
    if (is.character(base_crs)) {
        base_crs <- try(sf::st_crs(base_crs), silent = TRUE)
    }
    if (!inherits(base_crs, 'crs')) {
        stop('base crs argument is not a valid crs')
    }
    # parse wkt
    p_wkt <- parse_wkt(base_crs)
    # traverse wkt
    i_wkt <- traverse_list(p_wkt)
    # false easting/northing -> check epsg code in Table F.3 (https://docs.ogc.org/is/18-010r7/18-010r7.html#106)
    tabf.3 <- list(
        false_easting = c(8806, 8816, 8826),
        false_northing = c(8807, 8817, 8827),
        angle_northing = 8814
        )
    # get false easting (x-axis)
    i_x <- .get_indices(i_wkt, tabf.3[['false_easting']])
    # get false northing (y-axis)
    i_y <- .get_indices(i_wkt, tabf.3[['false_northing']])
    as.numeric(c(p_wkt[[i_x]], p_wkt[[i_y]]))
}

guess_ch <- function(x, y = NULL) {
    if (is.null(y)) {
        crs_gel <- attr(x, 'crs_gel')
        if (!is.null(crs_gel)) {
            return(crs_gel)
        }
        cnms <- guess_coords(x)
        if (is.matrix(x)) {
            y <- x[, cnms[2]]
            x <- x[, cnms[1]]
        } else {
            y <- x[[cnms[2]]]
            x <- x[[cnms[1]]]
        }
    }
    if (any(y > x)) {
        stop('ch coordinates should be provided as x pointing',
            ' towards east and y pointing towards north\n',
            'This is contrary to the official axis naming of "CH1903 / LV03" and "CH1903 / LV95"')
    }
    if (all(xbelow <- x < 1e6)) {
        return('CH1903/LV03')
    }
    if (any(xbelow)) stop('ch coordinate systems "CH1903 / LV03" and "CH1903 / LV95" cannot be mixed!')
    return('CH1903/LV95')
}


.coord_transf <- function(x, y, crs_from = NULL, crs_to = NULL, as_list = FALSE) {
    out <- sf::sf_project(crs_from, crs_to, cbind(x, y))
    if (as_list) {
        out <- as.list(as.data.table(out))
    }
    out
}

coord_transf <- function(x, crs_to,
    y = NULL, crs_from = NULL, 
    x_column = guess_coord_x(x, coord_system = crs_from), 
    y_column = guess_coord_y(x, coord_system = crs_from),
    as_list = FALSE, new_origin_at = NULL,
    old_origin_at = NULL, add_crs = inherits(x, 'data.frame'), append = TRUE) {
    # copy original x
    x_in <- copy(x)
    # check crs_to
    if (missing(crs_to)) {
        stop('argument crs_to is missing')
    }
    # convert to numeric
    if (!is.null(new_origin_at)) {
        new_origin_at <- as.numeric(new_origin_at)
    }
    # check crs from
    if (inherits(crs_from, 'staticMap')) {
        x <- map_to_wgs(crs_from, x, x_column = x_column,
            y_column = y_column)
        crs_from <- 'EPSG:4326'
    }
    if (is.null(crs_from)) {
        crs_from <- get_crs(x, TRUE)
    } else {
        crs_from <- fix_crs(crs_from, old_origin_at)
    }
    # check crs_from
    if (is.null(crs_from)) stop('crs of argument x is unknown')
    # fix crs to
    if (post_convert <- inherits(crs_to, 'staticMap')) {
        map <- crs_to
        crs_to <- 'EPSG:4326'
        new_origin_at <- NULL
    }
    # check crs_to
    if (!inherits(crs_to, 'crs') && length(crs_to) > 1) {
        stop('argument crs_to is not a valid crs input')
    }
    # save crs + offset
    crs_out <- list(crs = crs_to, new_origin_at = new_origin_at)
    crs_to <- fix_crs(crs_to, new_origin_at)
	if (inherits(x, "Sources") && ncol(x) == 4) {
        append <- FALSE
		out <- x
		out[, 2:3] <- .coord_transf(x[, 2], x[, 3], 
            crs_from = crs_from, crs_to = crs_to)
	} else if (inherits(x, "Sensors") && ncol(x) >= 7) {
        append <- FALSE
		out <- convert(x)
        out[, c('x-Coord (m)', 'y-Coord (m)')] <- .coord_transf(
            x[, 'x-Coord (m)'], x[, 'y-Coord (m)'], 
            crs_from = crs_from, crs_to = crs_to)
	} else {
		if (is.null(y)) {
            if (inherits(x, 'data.table')) {
                out <- data.table::copy(x)
            } else {
                out <- x
            }
			if (is.matrix(x)) {
                out[, c(x_column, y_column)] <- .coord_transf(
                    x[, x_column], x[, y_column],
                    crs_from = crs_from, crs_to = crs_to)
                colnames(out)[c(x_column, y_column)] <- c('x', 'y')
            } else if (inherits(x, 'area_sources')) {
                out[, c('x', 'y') :=  coord_transf(
                    x, y, crs_from = crs_from,
                    crs_to = crs_to, x_column = 'x',
                    y_column = 'y', as_list = TRUE)]
                attr(out, 'cadastre')[, c('x', 'y') := coord_transf(
                    attr(..x, 'cadastre')[, .(x, y)], crs_from = crs_from,
                    crs_to = crs_to, x_column = 'x',
                    y_column = 'y', as_list = TRUE)]
            } else if (inherits(x, 'data.table')) {
                if (is.numeric(x_column)) x_column <- names(out)[x_column]
                if (is.numeric(y_column)) y_column <- names(out)[y_column]
                out[, c(x_column, y_column) :=  .coord_transf(
                    get(x_column), get(y_column), crs_from = crs_from,
                    crs_to = crs_to, as_list = TRUE)]
            } else if (is.data.frame(x) || is.list(x)) {
                out <- x
				y <- x[[y_column]]
				x <- x[[x_column]]
                coords <- .coord_transf(x, y, 
                    crs_from = crs_from, crs_to = crs_to)
                colnames(coords) <- c('x', 'y')
                out[[x_column]] <- coords[, 'x']
                out[[y_column]] <- coords[, 'y']
			} else {
				y <- x[[y_column]]
				x <- x[[x_column]]
                coords <- .coord_transf(x, y, 
                    crs_from = crs_from, crs_to = crs_to)
                colnames(coords) <- c('x', 'y')
                out[[x_column]] <- coords[, 'x']
                out[[y_column]] <- coords[, 'y']
			}
		} else {
            out <- .coord_transf(x, y, 
                crs_from = crs_from, crs_to = crs_to)
        }
	}
    # add coords?
    if (append) {
        cx <- function(.x, ...) {
            if (inherits(.x, c('data.frame', 'matrix'))) {
                cbind(.x, ...)
            } else {
                c(.x, ...)
            }
        }
        if (is_wgs(crs_from)) {
            if (any(c('x', 'y') %in% names(out))) {
                out <- cx(x_in,
                    xnew = getElement(out, x_column),
                    ynew = getElement(out, y_column)
                )
            } else {
                out <- cx(x_in,
                    x = getElement(out, x_column),
                    y = getElement(out, y_column)
                )
            }
        } else {
            out <- cx(x_in,
                lat = getElement(out, y_column),
                lon = getElement(out, x_column)
            )
        }
    }
    # add crs
    if (add_crs) {
        out <- set_crs(out, crs = crs_out$crs, new_origin_at = crs_out$new_origin_at)
    }
    # check if we need to convert
    # from wgs84 to map
    if (post_convert) {
        out <- wgs_to_map(map, out, x_column = x_column,
            y_column = y_column)
    }
    if (as_list) {
        out <- as.list(as.data.frame(out))
    }
	return(out)
}

is_wgs <- function(x) {
    if (is.list(x)) {
        if ('wkt' %in% names(x)) {
            sub('.*,(\\d+)[]]*$', '\\1', x$wkt) == '4326'
        } else {
            stop('not able to guess coordinate system from the provided list')
        }
    } else if (is.numeric(x)) {
        x == 4326
    } else {
        grepl('^ *wgs[ _-]?84 *$', x, ignore.case = TRUE)
    }
}

##
# FIXME: add append argument
wgs_to_map <- function(MyMap, lon, lat = NULL, zoom,
    x_column = guess_coord_x(lon, coord_system = 'wgs84'), 
    y_column = guess_coord_y(lon, coord_system = 'wgs84'), append = TRUE) {
    require(RgoogleMaps, quietly = TRUE)
    if (missing(zoom)) zoom <- MyMap$zoom
	if (inherits(lon, "Sources") && ncol(lon) == 4) {
		out <- lon
		dummy <- LatLon2XY.centered(MyMap, lon[, 3], lon[, 2], zoom)
		out[, 2:3] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} else if (inherits(lon, "Sensors") && ncol(lon) >= 7) {
		out <- convert(lon)
		dummy <- LatLon2XY.centered(MyMap, lon[, "y-Coord (m)"], lon[, "x-Coord (m)"], zoom)
		out[, c("x-Coord (m)", "y-Coord (m)")] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} else {
		if (is.null(lat)) {
            if (inherits(lon, 'data.table')) {
                out <- copy(lon)
            } else {
                out <- lon
            }
            if (inherits(lon, 'area_sources')) {
                out[, c('x', 'y') := {
                    LatLon2XY.centered(MyMap, y, x, zoom)
                }]
                attr(out, 'cadastre')[, c('x', 'y') := {
                    LatLon2XY.centered(MyMap, y, x, zoom)
                }]
            } else if (inherits(lon, 'data.table')) {
                out[, c(x_column, y_column) := {
                    LatLon2XY.centered(MyMap, get(y_column), get(x_column), zoom)
                }]
            } else if (!is.matrix(lon) & !is.data.frame(lon)) {
                stop('fix me in wgs_to_map()')
			} else {
                coords <- LatLon2XY.centered(MyMap, 
                    lon[, y_column], lon[, x_column], zoom)
                out[, x_column] <- coords$newX
                out[, y_column] <- coords$newY
			}
		}
        return(out)
	}
}

##
map_to_wgs <- function(MyMap, x, y = NULL, zoom,
    x_column = guess_coord_x(x), 
    y_column = guess_coord_y(x)) {
    require(RgoogleMaps, quietly = TRUE)
    if (missing(zoom)) zoom <- MyMap$zoom
    # Sources?
	if (inherits(x, "Sources") && ncol(x) == 4) {
		out <- x
		dummy <- XY2LatLon(MyMap, x[, 3], x[, 2], zoom)
		out[, 2:3] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} 
    # Sensors?
    if (inherits(x, "Sensors") && ncol(x) >= 7) {
		out <- convert(x)
		dummy <- XY2LatLon(MyMap, x[, "y-Coord (m)"], x[, "x-Coord (m)"], zoom)
		out[, c("x-Coord (m)", "y-Coord (m)")] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} 
    # other
    if (is.null(y)) {
        if (inherits(x, 'data.table')) {
            out <- copy(x)
        } else {
            out <- x
        }
        if (inherits(x, 'area_sources')) {
            out[, c('x', 'y') := {
                XY2LatLon(MyMap, y, x, zoom)
            }]
            attr(out, 'cadastre')[, c('x', 'y') := {
                XY2LatLon(MyMap, y, x, zoom)
            }]
        } else if (inherits(x, 'data.table')) {
            out[, c(x_column, y_column) := {
                XY2LatLon(MyMap, get(y_column), get(x_column), zoom)
            }]
        } else if (!is.matrix(x) & !is.data.frame(x)) {
            stop('fix me in map_to_wgs()')
        } else {
            coords <- XY2LatLon(MyMap,
                x[, x_column], x[, y_column], zoom)
            out[, x_column] <- coords[, 'lon']
            out[, y_column] <- coords[, 'lat']
        }
    } else {
        stop('fix me in map_to_wgs() argument y is not null')
    }
    # return
    out
}

## 
ch_to_wgs <- function(x, y = NULL, append = TRUE) {
    crs_ch <- guess_ch(x, y)
    coord_transf(x, y, 
        crs_from = crs_ch, crs_to = 4326, append = append)
}
ch_to_map <- function(MyMap, x, y = NULL, ...) {
	WGS84 <- ch_to_wgs(x, y)
	wgs_to_map(MyMap, WGS84, ...)
} 
ch_to_user <- function(x, new_origin_at = NULL, y = NULL, append = TRUE) {
    crs_from <- guess_ch(x)
    coord_transf(x, y = y, crs_from = crs_from,
        crs_to = crs_from, new_origin_at = new_origin_at,
        append = append
    )
}

user_to_ch <- function(x, crs_from, origin_at, lv95 = TRUE,
    y = NULL, x_column = 'x', y_column = 'y', append = TRUE) {
    # check crs
    crs_x <- get_crs(x)
    if (!is.null(crs_x)) {
        if (missing(origin_at)) origin_at <- crs_x$origin_at
        if (missing(crs_from)) crs_from <- crs_x$crs_from
    } else if (missing(origin_at) || missing(crs_from)) {
        stop('crs_from and origin_at are both required!')
    }
    coord_transf(x, y, 
        crs_from = crs_from, crs_to = if (lv95) 2056 else 21781, 
        old_origin_at = origin_at, x_column = x_column, 
        y_column = y_column, append = append)
}
user_to_wgs <- function(x, crs_from, origin_at, y = NULL,
    x_column = 'x', y_column = 'y', append = TRUE) {
    # check crs
    crs_x <- get_crs(x)
    if (!is.null(crs_x)) {
        if (missing(origin_at)) origin_at <- crs_x$origin_at
        if (missing(crs_from)) crs_from <- crs_x$crs_from
    } else if (missing(origin_at) || missing(crs_from)) {
        stop('crs_from and origin_at are both required!')
    }
    coord_transf(x, y, x_column = x_column, y_column = y_column,
        crs_from = crs_from, crs_to = 4326, 
        old_origin_at = origin_at, append = append)
}
user_to_map <- function(MyMap, x, crs_from, origin_at,
    x_column = 'x', y_column = 'y', y = NULL, ...) {
	WGS84 <- user_to_wgs(x, y, crs_from = crs_from, 
        origin_at = origin_at, x_column = x_column, 
        y_column = y_column)
	wgs_to_map(MyMap, WGS84, x_column = x_column, y_column = y_column, ...)
}

wgs_to_user <- function(lon, crs_to, new_origin_at = NULL, lat = NULL, 
    append = TRUE) {
    coord_transf(lon, crs_to = crs_to, new_origin_at = new_origin_at,
        crs_from = 'wgs84', y = lat, append = append)
}
wgs_to_ch <- function(lon, lv95 = TRUE, lat = NULL, append = TRUE) {
    coord_transf(lon, lat, 
        crs_from = 4326, crs_to = if (lv95) 2056 else 21781,
        append = append
    )
}
# wgs_to_map -> see above

map_to_user <- function(MyMap, x, crs_to,
    new_origin_at = NULL, y = NULL, zoom,
    x_column = guess_coord_x(x),
    y_column = guess_coord_y(x)) {
    WGS84 <- map_to_wgs(MyMap, x, y, zoom, x_column, y_column)
    wgs_to_user(WGS84, crs_to, new_origin_at)
}
map_to_ch <- function(MyMap, x, lv95 = TRUE, y = NULL, zoom,
    x_column = guess_coord_x(x),
    y_column = guess_coord_y(x)) {
    WGS84 <- map_to_wgs(MyMap, x, y, zoom, x_column, y_column)
    wgs_to_ch(WGS84, lv95)
}
# map_to_wgs -> see above

## RgoogleMaps convenience wrappers
get_map <- function(loc, file = NULL, zoom = 16, 
    maptype = 'satellite', type = 'google',
    token = Sys.getenv('R_GOOGLE_MAPS'),
    crs_from = NULL, crs_origin_at = NULL, ...) {
    require(RgoogleMaps, quietly = TRUE)
    # get temporary file
    if (is.null(file)) {
        file <- tempfile('staticMap', fileext = '.png')
    }
    # convert loc
    loc_wgs <- coord_transf(loc, crs_to = 4326, 
        crs_from = crs_from, old_origin_at = crs_origin_at)
    colnames(loc_wgs) <- c('lon', 'lat')
    # check token
    if (type != 'google') token <- ''
    if (token == '' && type == 'google') {
        type <- 'osm'
        if (maptype == 'satellite') maptype <- 'roadmap'
    }
    # check center or bbox
    if (nrow(loc_wgs) == 1) {
        fu <- function(...) GetMap(loc_wgs[, c('lat', 'lon')], destfile = file, zoom = zoom,
            maptype = maptype, type = type, ...)
        cat('center: lon:', 
            round(loc_wgs[, 1], 1), 
            ' lat:', 
            round(loc_wgs[, 2], 1), 
            '\n')
    # get map
    } else {
        # get range
        loc_bbox <- apply(loc_wgs, 2, range, na.rm = TRUE)
        # get map
        fu <- function(...) GetMap.bbox(loc_bbox[, 1], loc_bbox[, 2], destfile = file, 
            maptype = maptype, type = type, ...)
        loc_ctr <- colMeans(loc_bbox)
        cat('bbox:\n  lon:', 
            round(loc_bbox[, 1], 1), 
            '\n  lat:', 
            round(loc_bbox[, 2], 1), 
            '\n')
    }
    if (token != '') {
        fu(API_console_key = token, ...)
    } else {
        fu(...)
    }
}
plot.staticMap <- function(x, y, ...) {
    PlotOnStaticMap(x, NEWMAP = FALSE, ...)
}

get_map_wgs <- function(loc, ...) get_map(loc, crs_from = 4326, ...)
get_map_ch <- function(loc, ...) get_map(loc, crs_from = guess_ch(loc), ...)

# https://gis.stackexchange.com/a/87159
metric_wgs84 <- function(lon, lat) {
    sf::st_crs(
        paste0('+proj=laea +lat_0=',
            lat, ' +lon_0=', lon, 
            ' +ellps=WGS84 +units=m +no_defs'
            )
        )
}

## ~~~~~~~~~~~~~ Tests

if (FALSE) {

    library(gel)
    library(RgoogleMaps)

    ## read gps data
    gps <- read.table('~/repos/4_Projects/6_EVEMBI/01_Waedenswil/MK2/Point_2.csv', fill = TRUE, sep = '\t', header = TRUE)
    rgmap <- ReadMapTile("~/repos/4_Projects/6_EVEMBI/01_Waedenswil/MK1/RSaves/Christoph/Waedi.png")
    gps_m <- as.matrix(gps[77:100, c('x', 'y')])
    gps_l <- as.list(gps)
    # TODO:
    # test Sources/Sensors/gral classes!

    ### ~~~~~~~~ coord_transf
    args(coord_transf)
    # matrix
    gm_ch03 <- coord_transf(gps_m, 'ch03', crs_from = 'wgs84')
    gm_ch95 <- coord_transf(gm_ch03, 'ch95')
    gm_ch03_user <- coord_transf(gm_ch95, 'ch03', new_origin_at = colMeans(gm_ch03))
    gm_user2 <- set_crs(gm_ch03_user, NULL)
    par(mfrow = c(2, 2))
    plot(gm_ch03, type = 'b')
    plot(gm_ch95, type = 'b')
    plot(gm_ch03_user, type = 'b')

    # data.frame
    coord_transf(gps, 'ch95', crs_from = 'wgs84')

    # list
    coord_transf(gps_l, 'ch03', crs_from = 'wgs84')

    # x + y
    coord_transf(gps$x, gps$y, crs_from = 'wgs84') # fails
    coord_transf(gps$x, y = gps$y, crs_from = 'wgs84') # fails
    coord_transf(gps$x, 'ch03', y = gps$y, crs_from = 'wgs84')

    # test convenience wrappers
    # ch_to_*
    ch_to_wgs(gm_ch03)
    ch_to_user(gm_ch95, c(2.6e6, 1.2e6))
    gm_map <- ch_to_map(rgmap, gm_ch03)
    PlotOnStaticMap(rgmap)
    points(gm_map, cex = 1.5, col = 'orange', lwd = 2)
    # wgs_to_*
    wgs_to_ch(gps_m)
    wgs_to_ch(gps_m, FALSE)
    args(wgs_to_user)
    wgs_to_user(gps_m, 'lv95', c(2.6e6, 1.2e6))
    gm_map2 <- wgs_to_map(rgmap, gps_m)
    PlotOnStaticMap(rgmap)
    points(gm_map2, cex = 1.5, col = 'orange', lwd = 2)
    # map_to_*
    map_to_wgs(rgmap, gm_map)
    map_to_ch(rgmap, gm_map)
    map_to_ch(rgmap, gm_map, FALSE)
    map_to_user(rgmap, gm_map, 'lv03', c(6e5, 2e5))
    # user_to_*
    # part crs defined by obj
    user_to_wgs(gm_ch03_user)
    user_to_ch(gm_ch03_user)
    user_to_ch(gm_ch03_user, lv95 = FALSE)
    user_to_map(rgmap, gm_ch03_user)
    # part crs undefined by obj
    user_to_wgs(gm_user2, 'ch03', attr(gm_ch03_user, 'origin'))
    user_to_ch(gm_user2, 'ch03', attr(gm_ch03_user, 'origin'))
    user_to_ch(gm_user2, lv95 = FALSE, 'ch03', attr(gm_ch03_user, 'origin'))
    user_to_map(rgmap, gm_user2, 'ch03', attr(gm_ch03_user, 'origin'))

}
