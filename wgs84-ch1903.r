## ~~~~ new functions using EPSG ~~~~ ##

## input data format:
# 1.) blsmodelr -> Sources and Sensors
# 2.) matrix/data.frame/list with x & y coords and other entries

## functions:
# coord_transf
#   -> crs: map, local, crs
#   -> local or map: check origin not missing!!!
#   -> crs_to/crs_from == Map?
# ev. swap_coords (ch1903 <-> wgs84)
# local_to_...
# ch1903_to_...
# wgs84_to_...
# map_to_...
#
# - 1 function to change in-between all of them
# - use S3 methods for different classes
# - use attributes to store offset & coordinates
# - use st_point for offset & st_multipoint for resid
# - more partial functions

# check EPSG codes on https://epsg.io
# wgs84: 4326
# ch1903/LV03: 21781
# ch1903/LV95: 2056

library(sf, quietly = TRUE)


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
    crs <- st_crs(crs)
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
    attr(obj, 'crs_gel') <- crs
    attr(obj, 'origin') <- new_origin_at
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
        unlist(lapply(obj, class))
    })
setMethod('get_col_classes',
    signature(obj = 'matrix'),
    function(obj) {
        setNames(
            rep(class(obj[[1]]), ncol(obj)),
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
guess_coord_x <- function(obj, value = TRUE) {
    isn <- check_numeric(obj)
    nms <- get_names(obj)[isn]
    # check null
    if (is.null(nms)) {
        if (is.character(obj)) {
            nms <- obj
        } else {
            return(which(isn)[1])
        }
    }
    # check x
    x_nm <- grep('^((x|X)(.*(c|C)oord.*)?)$|^((c|C)oord.*(x|X))$', nms, value = value)
    if (length(x_nm) == 0) {
        # check lon
        x_nm <- grep('^(l|L)on((g|gitude).*)?$', nms, value = value)
        if (length(x_nm) == 0) return(which(isn)[1])
    }
    x_nm
}
guess_coord_y <- function(obj, value = TRUE) {
    isn <- check_numeric(obj)
    nms <- get_names(obj)[isn]
    # check null
    if (is.null(nms)) {
        if (is.character(obj)) {
            nms <- obj
        } else {
            return(which(isn)[2])
        }
    }
    # check y
    y_nm <- grep('^((y|Y)(.*(c|C)oord.*)?)$|^((c|C)oord.*(y|Y))$', nms, value = value)
    if (length(y_nm) == 0) {
        # check lon
        y_nm <- grep('^(l|L)at(itude.*)?$', nms, value = value)
        if (length(y_nm) == 0) return(which(isn)[2])
    }
    y_nm
}
guess_coords <- function(obj, value = TRUE) {
    c(
        guess_coord_x(obj)[1], 
        guess_coord_y(obj)[1]
    )
}


# ## ~~~ mess with coordinate attributes
# coord_attr_name <- function(obj, what, which_coord = NULL) {
#     attr_name <- switch(what
#         # coordinate column names
#         , 'coord_name' = {
#             stopifnot(length(which_coord) == 1 && which_coord %in% c('x', 'y'))
#             paste0(which_coord, '_coord')
#         }
#         # else
#         , stop('attribute ', what, ' does not exist')
#         )
# }
# coord_attr <- function(obj, what, which_coord = NULL) {
#     attr(obj, coord_attr_name(obj, what, which_coord = which_coord))
# }
# 'coord_attr<-' <- function(obj, what, which_coord = NULL, value) {
#     attr(obj, coord_attr_name(obj, what, which_coord = which_coord)) <- value
#     obj
# }
# 'coord_name_x<-' <- function(obj, value) {
#     # check if name exists
#     if (length(value) > 1 || !names_exist(obj, value)) {
#         stop('value with length > 1 cannot be assigned to coord_name x')
#     }
#     # assign
#     coord_attr(obj, 'coord_name', 'x') <- value
#     obj
# }
# 'coord_name_y<-' <- function(obj, value) {
#     # check if name exists
#     if (length(value) > 1 || !names_exist(obj, value)) {
#         stop('value with length > 1 cannot be assigned to coord_name y')
#     }
#     # assign
#     coord_attr(obj, 'coord_name', 'y') <- value
#     obj
# }
# 'coord_names<-' <- function(obj, value) {
#     nms <- get_names(value)
#     if (!is.null(nms)) {
#         stopifnot(all(nms %in% c('x', 'y')))
#     } else {
#         names(value) <- c('x', 'y')
#     }
#     coord_name_x(obj) <- value[['x']]
#     coord_name_y(obj) <- value[['y']]
#     obj
# }
# coord_name_x <- function(obj) {
#     coord_attr(obj, 'coord_name', 'x')
# }
# coord_name_y <- function(obj) {
#     coord_attr(obj, 'coord_name', 'y')
# }
# coord_names <- function(obj) {
#     c(
#         x = coord_name_x(obj),
#         y = coord_name_y(obj)
#         )
# }

# ## guess coordinate columns

# ## automatically assign/fix missing coordinate names
# fix_coord_name_x <- function(obj) {
#     # check if coord name x/y is missing
#     cn <- coord_name_x(obj)
#     # if existing
#     if (!is.null(cn)) {
#         # check if column exists and if numeric
#         names_exist(cn)
#             # if inexistent -> missing
#             # else return cn
#     }
#     # coord name not ok
#     # search for x/y or lat/lon
#     cn <- guess_coord_x(obj)
#     # check numeric columns (1st -> x, 2nd -> y)
# }
# fix_coord_name_y
# fix_coord_names

# ## additional functions:
# # get coordinate values (format matrix? list?)
# coords <- function(obj)
# 'coords<-' <- function(obj)
# coord_x
# coord_y
# 'coord_x<-'
# 'coord_y<-'

## --- local -> map ---
## local -> crs (proj) -> wgs84 -> map

# ## ~~~ local x/y -> ch1903/crs_proj
# local_to_ch1903/crs_proj:
#     - obj_local, origin (in ch1903/crs_proj), LV03 = FALSE

# local_to_ch1903 <- function(obj_local, origin = NULL, LV03 = FALSE) {
#     if (LV03) {
#         crs_to <- 'EPSG:21781'
#     } else {
#         crs_to <- 'EPSG:2056'
#     }
#     transform(obj_local, crs_from = NULL, crs_to = crs_to)
# }

# ## ~~~ ch1903/crs_proj -> wgs84
# ch1903_to_wgs84:
#     - obj_ch1903

# ## ~~~ wgs84 -> map
# wgs84_to_map:
#     - obj_wgs84, map

# ## ~~~ map x/y -> local
# map_to_local:
#     - obj_map, map, origin (in crs_local), crs_local != NULL

.change_coords <- function(x, y, crs_from = NULL, crs_to = NULL) {
    sf_project(crs_from, crs_to, cbind(x, y))
}

change_coords <- function(x, crs_to,
    y = NULL, crs_from = NULL, 
    x_column = guess_coord_x(x), 
    y_column = guess_coord_y(x),
    as_list = FALSE, new_origin_at = NULL,
    old_origin_at = NULL, add_crs = TRUE) {
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
		out <- x
		out[, 2:3] <- .change_coords(x[, 2], x[, 3], 
            crs_from = crs_from, crs_to = crs_to)
	} else if (inherits(x, "Sensors") && ncol(x) >= 7) {
		out <- convert(x)
        out[, c('x-Coord (m)', 'y-Coord (m)')] <- .change_coords(
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
                out[, c(x_column, y_column)] <- .change_coords(
                    x[, x_column], x[, y_column],
                    crs_from = crs_from, crs_to = crs_to)
            } else if (inherits(x, 'area_sources')) {
                out[, c('x', 'y') :=  change_coords(
                    x, y, crs_from = crs_from,
                    crs_to = crs_to, x_column = 'x',
                    y_column = 'y', as_list = TRUE)]
                attr(out, 'cadastre')[, c('x', 'y') := change_coords(
                    attr(..x, 'cadastre')[, .(x, y)], crs_from = crs_from,
                    crs_to = crs_to, x_column = 'x',
                    y_column = 'y', as_list = TRUE)]
            } else if (inherits(x, 'data.table')) {
                if (is.numeric(x_column)) x_column <- names(out)[x_column]
                if (is.numeric(y_column)) y_column <- names(out)[y_column]
                out[, c(x_column, y_column) :=  change_coords(
                    x, y, crs_from = crs_from,
                    crs_to = crs_to, as_list = TRUE)]
            } else if (is.data.frame(x) || is.list(x)) {
                out <- x
				y <- x[[y_column]]
				x <- x[[x_column]]
                coords <- .change_coords(x, y, 
                    crs_from = crs_from, crs_to = crs_to)
                colnames(coords) <- c('x', 'y')
                out[[x_column]] <- coords[, 'x']
                out[[y_column]] <- coords[, 'y']
			} else {
				y <- x[[y_column]]
				x <- x[[x_column]]
                out <- .change_coords(x, y, 
                    crs_from = crs_from, crs_to = crs_to)
                colnames(out) <- c('x', 'y')
                if (as_list) {
                    out <- as.list(as.data.frame(out))
                }
			}
		} else {
            out <- .change_coords(x, y, 
                crs_from = crs_from, crs_to = crs_to)
            colnames(out) <- c('x', 'y')
            if (as_list) {
                out <- as.list(as.data.frame(out))
            }
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
	return(out)
}

##
wgs_to_map <- function(MyMap, lon, lat = NULL, zoom,
    x_column = guess_coord_x(lon), 
    y_column = guess_coord_y(lon)) {
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
    if (y > x) {
        stop('ch coordinates should be provided as x pointing',
            ' towards east and y pointing towards north\n',
            'This is contrary to the official axis naming of "CH1903 / LV03" and "CH1903 / LV95"')
    }
    if (x < 1e6) {
        return('CH1903/LV03')
    }
    return('CH1903/LV95')
}

## 
ch_to_wgs <- function(x, y = NULL) {
    crs_ch <- guess_ch(x, y)
    change_coords(x, y, 
        crs_from = crs_ch, crs_to = 4326)
}
ch_to_map <- function(MyMap, x, y = NULL, ...) {
	WGS84 <- ch_to_wgs(x, y)
	wgs_to_map(MyMap, WGS84, ...)
} 
ch_to_user <- function(x, new_origin_at = NULL, y = NULL) {
    crs_from <- guess_ch(x)
    change_coords(x, y = y, crs_from = crs_from,
        crs_to = crs_from, new_origin_at = new_origin_at)
}

user_to_ch <- function(x, crs_from, origin_at, lv95 = TRUE,
    y = NULL, x_column = 'x', y_column = 'y') {
    # check crs
    crs_x <- get_crs(x)
    if (!is.null(crs_x)) {
        if (missing(origin_at)) origin_at <- crs_x$origin_at
        if (missing(crs_from)) crs_from <- crs_x$crs_from
    } else if (missing(origin_at) || missing(crs_from)) {
        stop('crs_from and origin_at are both required!')
    }
    change_coords(x, y, 
        crs_from = crs_from, crs_to = if (lv95) 2056 else 21781, 
        old_origin_at = origin_at, x_column = x_column, 
        y_column = y_column)
}
user_to_wgs <- function(x, crs_from, origin_at, y = NULL,
    x_column = 'x', y_column = 'y') {
    # check crs
    crs_x <- get_crs(x)
    if (!is.null(crs_x)) {
        if (missing(origin_at)) origin_at <- crs_x$origin_at
        if (missing(crs_from)) crs_from <- crs_x$crs_from
    } else if (missing(origin_at) || missing(crs_from)) {
        stop('crs_from and origin_at are both required!')
    }
    change_coords(x, y, x_column = x_column, y_column = y_column,
        crs_from = crs_from, crs_to = 4326, 
        old_origin_at = origin_at)
}
user_to_map <- function(MyMap, x, crs_from, origin_at,
    x_column = 'x', y_column = 'y', y = NULL, ...) {
	WGS84 <- user_to_wgs(x, y, crs_from = crs_from, 
        origin_at = origin_at, x_column = x_column, 
        y_column = y_column)
	wgs_to_map(MyMap, WGS84, x_column = x_column, y_column = y_column, ...)
}

wgs_to_user <- function(lon, crs_to, new_origin_at = NULL, lat = NULL) {
    change_coords(lon, crs_to = crs_to, new_origin_at = new_origin_at,
        crs_from = 'wgs84', y = lat)
}
wgs_to_ch <- function(lon, lv95 = TRUE, lat = NULL) {
    change_coords(lon, lat, 
        crs_from = 4326, crs_to = if (lv95) 2056 else 21781)
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
    loc_wgs <- change_coords(loc, crs_to = 4326, 
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

parse_wkt <- function(crs) {
    # crs <- st_crs('EPSG:2056')
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

traverse_list <- function(x, i = numeric(0)) {
    if (is.list(x)) {
        out <- mapply(traverse_list, x, i = lapply(seq_along(x), function(a) c(i, a)), SIMPLIFY = FALSE)
        # out <- mapply(traverse_list, x, i = paste(i, seq_along(x), sep = ','), SIMPLIFY = FALSE)
        do.call(rbind, out)
        # do.call(rbind, lapply(x, traverse_list, i = paste(i, seq_along(x), sep = ':')))
    } else {
        cbind(list(i), x)
    }
}

# find EPSG code:
get_index <- function(t_wkt, code) {
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
        base_crs <- try(st_crs(base_crs), silent = TRUE)
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
        i_x <- get_index(i_wkt, tabf.3[['false_easting']])
        # get false northing (y-axis)
        i_y <- get_index(i_wkt, tabf.3[['false_northing']])
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
        i_n <- get_index(i_wkt, tabf.3[['angle_northing']])
        # replace value
        p_wkt[[i_n]] <- new_angle
    }
    # replace name
    p_wkt[[1]][[1]] <- "User"
    st_crs(deparse_wkt(p_wkt))
}
get_origin <- function(base_crs) {
    if (missing(base_crs)) {
        stop('base crs is missing')
    }
    if (is.character(base_crs)) {
        base_crs <- try(st_crs(base_crs), silent = TRUE)
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
    i_x <- get_index(i_wkt, tabf.3[['false_easting']])
    # get false northing (y-axis)
    i_y <- get_index(i_wkt, tabf.3[['false_northing']])
    as.numeric(c(p_wkt[[i_x]], p_wkt[[i_y]]))
}

## ~~~~~~~~~~~~~ Tests

if (FALSE) {

    library(RgoogleMaps)

    ## read gps data
    gps <- read.table('~/repos/4_Projects/6_EVEMBI/01_Waedenswil/MK2/Point_2.csv', fill = TRUE, sep = '\t', header = TRUE)
    rgmap <- ReadMapTile("~/repos/4_Projects/6_EVEMBI/01_Waedenswil/MK1/RSaves/Christoph/Waedi.png")
    gps_m <- as.matrix(gps[77:100, c('x', 'y')])
    gps_l <- as.list(gps)
    # TODO:
    # test Sources/Sensors/gral classes!

    ### ~~~~~~~~ change_coords
    args(change_coords)
    # matrix
    gm_ch03 <- change_coords(gps_m, 'ch03', crs_from = 'wgs84')
    gm_ch95 <- change_coords(gm_ch03, 'ch95')
    gm_ch03_user <- change_coords(gm_ch95, 'ch03', new_origin_at = colMeans(gm_ch03))
    gm_user2 <- set_crs(gm_ch03_user, NULL)
    par(mfrow = c(2, 2))
    plot(gm_ch03, type = 'b')
    plot(gm_ch95, type = 'b')
    plot(gm_ch03_user, type = 'b')

    # data.frame
    change_coords(gps, 'ch95', crs_from = 'wgs84')

    # list
    change_coords(gps_l, 'ch03', crs_from = 'wgs84')

    # x + y
    change_coords(gps$x, gps$y, crs_from = 'wgs84') # fails
    change_coords(gps$x, y = gps$y, crs_from = 'wgs84') # fails
    change_coords(gps$x, 'ch03', y = gps$y, crs_from = 'wgs84')

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
