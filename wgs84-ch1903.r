## ~~~~ new functions using EPSG ~~~~ ##

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

## ~~~ general helpers
get_names <- function(obj) UseMethod('get_names')
get_names.default <- function(obj) names(obj)
get_names.matrix <- function(obj) colnames(obj)
get_length <- function(obj) UseMethod('get_length')
get_length.default <- function(obj) length(obj)
get_length.matrix <- function(obj) ncol(obj)
names_exist <- function(obj, col) UseMethod('names_exist')
setMethod('names_exist', 
    signature(obj = 'ANY', col = 'numeric'),
    function(obj, col) {
        col <= get_length(obj)
    })
setMethod('names_exist', 
    signature(obj = 'ANY', col = 'character'),
    function(obj, col) {
        col %in% get_names(obj)
    })
setMethod('names_exist', 
    signature(obj = 'ANY', col = 'NULL'),
    function(obj, col) {
        TRUE
    })

## ~~~ mess with coordinate attributes
coord_attr_name <- function(obj, what, which_coord = NULL) {
    attr_name <- switch(what
        # coordinate column names
        , 'coord_name' = {
            stopifnot(length(which_coord) == 1 && which_coord %in% c('x', 'y'))
            paste0(which_coord, '_coord')
        }
        # else
        , stop('attribute ', what, ' does not exist')
        )
}
coord_attr <- function(obj, what, which_coord = NULL) {
    attr(obj, coord_attr_name(obj, what, which_coord = which_coord))
}
'coord_attr<-' <- function(obj, what, which_coord = NULL, value) {
    attr(obj, coord_attr_name(obj, what, which_coord = which_coord)) <- value
    obj
}
'coord_name_x<-' <- function(obj, value) {
    # check if name exists
    if (length(value) > 1 || !names_exist(obj, value)) {
        stop('value with length > 1 cannot be assigned to coord_name x')
    }
    # assign
    coord_attr(obj, 'coord_name', 'x') <- value
    obj
}
'coord_name_y<-' <- function(obj, value) {
    # check if name exists
    if (length(value) > 1 || !names_exist(obj, value)) {
        stop('value with length > 1 cannot be assigned to coord_name y')
    }
    # assign
    coord_attr(obj, 'coord_name', 'y') <- value
    obj
}
'coord_names<-' <- function(obj, value) {
    nms <- get_names(value)
    if (!is.null(nms)) {
        stopifnot(all(nms %in% c('x', 'y')))
    } else {
        names(value) <- c('x', 'y')
    }
    coord_name_x(obj) <- value[['x']]
    coord_name_y(obj) <- value[['y']]
    obj
}
coord_name_x <- function(obj) {
    guess coord name here if missing?
    coord_attr(obj, 'coord_name', 'x')
}
coord_name_y <- function(obj) {
    coord_attr(obj, 'coord_name', 'y')
}
coord_names <- function(obj) {
    c(
        x = coord_name_x(obj),
        y = coord_name_y(obj)
        )
}

guess_coord_x <- function(obj, value = TRUE) {
    nms <- get_names(obj)
    # check null
    if (is.null(nms)) {
        if (is.character(obj)) {
            nms <- obj
        } else {
            return(NULL)
        }
    }
    # check x
    x_nm <- grep('^((x|X)(.*(c|C)oord.*)?)$|^((c|C)oord.*(x|X))$', nms, value = value)
    if (length(x_nm) == 0) {
        # check lon
        x_nm <- grep('^(l|L)on((g|gitude).*)?$', nms, value = value)
    }
    x_nm
}
guess_coord_y <- function(obj, value = TRUE) {
    nms <- get_names(obj)
    # check null
    if (is.null(nms)) {
        if (is.character(obj)) {
            nms <- obj
        } else {
            return(NULL)
        }
    }
    # check y
    y_nm <- grep('^((y|Y)(.*(c|C)oord.*)?)$|^((c|C)oord.*(y|Y))$', nms, value = value)
    if (length(y_nm) == 0) {
        # check lon
        y_nm <- grep('^(l|L)at(itude.*)?$', nms, value = value)
    }
    y_nm
}
guess_coords <- function(obj, value = TRUE) {
    c(
        guess_coord_x(obj)[1], 
        guess_coord_y(obj)[1]
    )
}

## automatically assign/fix missing coordinate names
fix_coord_col <- function(obj) {
    # check if coord name x/y is missing
    # if not missing
        # check if column exists and if numeric
            # if inexistent -> missing
            # else ok
    # if missing
        # search for x/y or lat/lon
        # check numeric columns (1st -> x, 2nd -> y)
}

## additional functions:
# get coordinate values (format matrix? list?)
coords <- function(obj)
'coords<-' <- function(obj)
coord_x
coord_y
'coord_x<-'
'coord_y<-'

Also: make only certain functions visible for final version

## --- local -> map ---
## local -> crs (proj) -> wgs84 -> map

## ~~~ local x/y -> ch1903/crs_proj
local_to_ch1903/crs_proj:
    - obj_local, origin (in ch1903/crs_proj), LV03 = FALSE

local_to_ch1903 <- function(obj_local, origin = NULL, LV03 = FALSE) {
    if (LV03) {
        crs_to <- 'EPSG:21781'
    } else {
        crs_to <- 'EPSG:2056'
    }
    transform(obj_local, crs_from = NULL, crs_to = crs_to)
}

## ~~~ ch1903/crs_proj -> wgs84
ch1903_to_wgs84:
    - obj_ch1903

## ~~~ wgs84 -> map
wgs84_to_map:
    - obj_wgs84, map

## ~~~ map x/y -> local
map_to_local:
    - obj_map, map, origin (in crs_local), crs_local != NULL

.change_coords <- function(x, y, crs_from = NULL, crs_to = NULL, swap_xy = FALSE,
    offset = c(0, 0)) {
    # insist on numeric
    x <- as.numeric(x) + offset[1]
    y <- as.numeric(y) + offset[2]
    # swap y > x?
    if (swap_xy && (
        (y[1] < 100 && y[1] < x[1]) || 
        (y[1] > 100 && y[1] > x[1])
        ))
    {
        xcp <- x
        x <- y
        y <- xcp
    }
    # check xy crs
    if (is.null(crs_from)) {
        if (y[1] < 90) {
            # wgs84
            crs_from <- 'EPSG:4326'
        } else if (y[1] > 5e5) {
            # LV95
            crs_from <- 'EPSG:2056'
        } else {
            # LV03
            crs_from <- 'EPSG:21781'
        }
    }
    # check output crs
    if (is.null(crs_to)) {
        if (crs_from == 'EPSG:4326' || crs_from == 4326) {
            # LV95
            crs_to <- 'EPSG:2056'
        } else {
            # WGS84
            crs_to <- 'EPSG:4326'
        }
    }
    obj <- st_sfc(st_multipoint(cbind(x, y)))
    st_crs(obj) <- crs_from
    out <- st_transform(obj, crs_to)
    as.matrix(out[[1]])
}

change_coords <- function(x, y = NULL, crs_from = NULL, 
    crs_to = NULL, swap_xy = FALSE, x_col = 1, y_col = 2,
    as_list = FALSE, offset = c(0, 0)) {
    offset <- as.numeric(offset)
	if (inherits(x, "Sources") && ncol(x) == 4) {
		out <- x
		out[, 2:3] <- .change_coords(x[, 2], x[, 3], 
            crs_from = crs_from, crs_to = crs_to, 
            swap_xy = swap_xy, offset = offset)
	} else if (inherits(x, "Sensors") && ncol(x) >= 7) {
		out <- convert(x)
        out[, c('x-Coord (m)', 'y-Coord (m)')] <- .change_coords(
            x[, 'x-Coord (m)'], x[, 'y-Coord (m)'], 
            crs_from = crs_from, crs_to = crs_to, 
            swap_xy = swap_xy, offset = offset)
	} else {
		if (is.null(y)) {
            out <- copy(x)
			if (is.matrix(x)) {
                out[, c(x_col, y_col)] <- .change_coords(
                    x[, x_col], x[, y_col],
                    crs_from = crs_from, crs_to = crs_to, 
                    swap_xy = swap_xy, offset = offset)
            } else if (inherits(x, 'area_sources')) {
                out[, c('x', 'y') :=  change_coords(
                    x, y, crs_from = crs_from,
                    crs_to = crs_to, swap_xy = swap_xy, x_col = 'x',
                    y_col = 'y', as_list = TRUE, offset = offset)]
                attr(out, 'cadastre')[, c('x', 'y') := change_coords(
                    attr(..x, 'cadastre')[, .(x, y)], crs_from = crs_from,
                    crs_to = crs_to, swap_xy = swap_xy, x_col = 'x',
                    y_col = 'y', as_list = TRUE, offset = offset)]
            } else if (is.data.table(x)) {
                if (is.numeric(x_col)) x_col <- names(out)[x_col]
                if (is.numeric(y_col)) y_col <- names(out)[y_col]
                out[, c(x_col, y_col) :=  change_coords(
                    x, y, crs_from = crs_from,
                    crs_to = crs_to, swap_xy = swap_xy, 
                    as_list = TRUE, offset = offset)]
            } else if (is.data.frame(x) || is.list(x)) {
                stop('Fix me in change_coords()!')
			} else {
				y <- x[[y_col]]
				x <- x[[x_col]]
                out <- .change_coords(x, y, 
                    crs_from = crs_from, crs_to = crs_to, 
                    swap_xy = swap_xy, offset = offset)
                colnames(out) <- c('x', 'y')
                if (as_list) {
                    out <- as.list(as.data.frame(out))
                }
			}
		} else {
            out <- .change_coords(x, y, 
                crs_from = crs_from, crs_to = crs_to, 
                swap_xy = swap_xy, offset = offset)
            colnames(out) <- c('x', 'y')
            if (as_list) {
                out <- as.list(as.data.frame(out))
            }
        }
	}
	return(out)
}

## 
wgs_to_ch <- function(lon, lat = NULL) {
    change_coords(lon, lat, 
        crs_from = 4326, crs_to = 2056,
        swap_xy = TRUE)
}
ch_to_wgs <- function(x, y = NULL) {
    change_coords(x, y, 
        crs_from = NULL, crs_to = 4326,
        swap_xy = TRUE)
}

##
wgs_to_map <- function(MyMap, lat, lon = NULL, zoom,
    x_col = 1, y_col = 2) {
    require(RgoogleMaps, quietly = TRUE)
    if (missing(zoom)) zoom <- MyMap$zoom
	if (inherits(lat, "Sources") && ncol(lat) == 4) {
		out <- lat
		dummy <- LatLon2XY.centered(MyMap, lat[, 3], lat[, 2], zoom)
		out[, 2:3] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} else if (inherits(lat, "Sensors") && ncol(lat) >= 7) {
		out <- convert(lat)
		dummy <- LatLon2XY.centered(MyMap, lat[, "y-Coord (m)"], lat[, "x-Coord (m)"], zoom)
		out[, c("x-Coord (m)", "y-Coord (m)")] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} else {
		if (is.null(lon)) {
            out <- copy(lat)
            if (inherits(lat, 'area_sources')) {
                out[, c('x', 'y') := {
                    LatLon2XY.centered(MyMap, y, x, zoom)
                }]
                attr(out, 'cadastre')[, c('x', 'y') := {
                    LatLon2XY.centered(MyMap, y, x, zoom)
                }]
            } else if (is.data.table(lat)) {
                out[, c(x_col, y_col) := {
                    LatLon2XY.centered(MyMap, get(y_col), get(x_col), zoom)
                }]
            } else if (!is.matrix(lat) & !is.data.frame(lat)) {
                stop('fix me in wgs_to_map()')
			} else {
                stop('fix me in wgs_to_map() data.frame/matrix')
			}
		}
        return(out)
	}
}

ch_to_map <- function(MyMap, x, y = NULL, ...) {
	WGS84 <- ch_to_wgs(x, y)
	wgs_to_map(MyMap, WGS84, ...)
} 
ch_to_wgs <- function(x, y = NULL) {
    change_coords(x, y, 
        crs_from = NULL, crs_to = 4326,
        swap_xy = TRUE)
}

xy_to_ch <- function(x, y = NULL, crs_from, offset,
    x_col = 'x', y_col = 'y') {
    if (missing(offset) || missing(crs_from)) {
        stop('crs_from and offset between xy and crs_from are both required!')
    }
    change_coords(x, y, 
        crs_from = crs_from, crs_to = 2056, 
        swap_xy = TRUE, offset = offset, x_col = x_col, y_col = y_col)
}
xy_to_wgs <- function(x, y = NULL, crs_from, offset, x_col = 'x', y_col = 'y') {
    if (missing(offset) || missing(crs_from)) {
        stop('crs_from and offset between xy and crs_from are both required!')
    }
    change_coords(x, y, x_col = x_col, y_col = y_col,
        crs_from = crs_from, crs_to = 4326, 
        swap_xy = TRUE, offset = offset)
}
xy_to_map <- function(MyMap, x, y = NULL, crs_from, offset,
    x_col = 'x', y_col = 'y', ...) {
	WGS84 <- xy_to_wgs(x, y, crs_from = crs_from, offset = offset, x_col = x_col, y_col = y_col)
	wgs_to_map(MyMap, WGS84, x_col = x_col, y_col = y_col, ...)
}

## RgoogleMaps convenience wrappers
get_map <- function(loc, file = NULL, zoom = 16, 
    maptype = 'satellite', type = 'google',
    token = Sys.getenv('R_GOOGLE_MAPS'), crs_from = NULL,
    fix_xy = TRUE, ...) {
    require(RgoogleMaps, quietly = TRUE)
    # get temporary file
    if (is.null(file)) {
        file <- tempfile('staticMap', fileext = '.png')
    }
    # convert loc
    loc_wgs <- change_coords(loc, crs_to = 4326, swap_xy = fix_xy, 
        crs_from = crs_from)
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

# TODO:
# use S3 methods
# user coord (m) -> crs (m)
# user coord (m) <- crs (m)
# crs (m) -> crs/wgs84
# crs (m) <- crs/wgs84
# wgs84 -> RgoogleMaps
# wgs84 <- RgoogleMaps

