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

## guess coordinate columns
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

obj <- data.frame(x = 1:3, y = 1:3)

## automatically assign/fix missing coordinate names
fix_coord_name_x <- function(obj) {
    # check if coord name x/y is missing
    cn <- coord_name_x(obj)
    # if existing
    if (!is.null(cn)) {
        # check if column exists and if numeric
        names_exist(cn)
            # if inexistent -> missing
            # else return cn
    }
    # coord name not ok
    # search for x/y or lat/lon
    cn <- guess_coord_x(obj)
    # check numeric columns (1st -> x, 2nd -> y)
}
fix_coord_name_y
fix_coord_names

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

# sf package
library(sf)

obj <- data.frame(lon = 7.4, lat = 47)
o <- st_sfc(st_point(as.matrix(obj)))
st_crs(o) <- 'EPSG:4326'
cat(st_as_text(st_crs(o), projjson = TRUE))

o_ch <- st_transform(o, crs = 'EPSG:21781')
cat(st_as_text(st_crs(o_ch), projjson = TRUE))

crs_ch <- st_as_text(st_crs(o_ch))
## user defined
crs_user <- sub('CH1903 / LV03', 'User', crs_ch)
crs_user <- sub('600000', '0', crs_user)

st_transform(o, crs = crs_ch)
st_transform(o, crs = crs_user)

st_as_text(st_crs('EPSG:21781'))
st_as_text(st_crs('EPSG:2056'))
st_as_text(st_crs('EPSG:2057'))
# TODO: create local/user-defined crs from projection with 
#   - UNIT[\"metre\",1] -> check if true
#   - PARAMETER[\"false_northing\",200000] -> check for false_northing
#   - PARAMETER[\"false_easting\",600000] -> check for false_easting
#   - what about scale_factor != 1??
crs_test <- st_as_text(st_crs('EPSG:2057'))
#   - UNIT[\"metre\",1] -> check if true
grepl('UNIT["metre",1]', crs_test, fixed = TRUE)
#   - PARAMETER[\"false_northing\",3044969.194] -> check for false_northing
sub('.*PARAMETER\\["false_northing",(\\d(\\d|[.])+)\\].*', '\\1', crs_test)
#   - PARAMETER[\"false_easting\",658377.437] -> check for false_easting
sub('.*PARAMETER\\["false_easting",(\\d(\\d|[.])+)\\].*', '\\1', crs_test)
#   - what about scale_factor != 1?? -> don't mess with it :)
test_ch <- st_sfc(st_multipoint(cbind(x = c(600000, 610000), y = c(200000, 210000))))
st_crs(test_ch) <- 'EPSG:21781'
# st_transform(test_ch, crs = 'EPSG:4326')
crs_ch <- st_as_text(st_crs(o_ch))
## user defined
crs_user <- sub('CH1903 / LV03', 'User', crs_ch)
crs_user <- sub('600000', '0', crs_user)
crs_user <- sub('200000', '0', crs_user)
test_2 <- st_transform(test_ch, crs = crs_user)
crs_user2 <- sub('scale_factor",1', 'scale_factor",0.5', crs_user, fixed = TRUE)
test_2_scaled <- st_transform(test_ch, crs = crs_user2)
# false easting/northing -> check epsg code in Table F.3 (https://docs.ogc.org/is/18-010r7/18-010r7.html#106)
tabf.3 <- list(
    false_easting = c(8806, 8816, 8826),
    false_northing = c(8807, 8817, 8827),
    angle_northing = 8814
    )

st_crs('+proj=laea +lat_0=<LAT> +lon_0=<LON> +ellps=WGS84 +units=m +no_defs')
crs_wgs84_m <- st_crs('+proj=laea +lat_0=47 +lon_0=7.4 +ellps=WGS84 +units=m +no_defs')
crs_wgs84 <- st_crs('EPSG:4326')
cat(st_as_text(crs_wgs84_m, projjson = TRUE))
require(jsonlite)
crs_json1 <- parse_json(st_as_text(crs_wgs84_m, projjson = TRUE))
crs_json2 <- parse_json(st_as_text(crs_wgs84, projjson = TRUE))

test_wgs84_m <- st_transform(test_ch, crs = crs_wgs84_m)

apply(st_coordinates(test_ch), 2, diff)
apply(st_coordinates(test_wgs84_m), 2, diff)

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

x <- parse_wkt(crs_wgs84_m)
y <- parse_wkt(st_crs('EPSG:2056'))

lapply(x[['PROJCRS']][['CONVERSION']], names)
lapply(y[['PROJCRS']][['CONVERSION']], names)

y[['PROJCRS']][['CONVERSION']][['PARAMETER']]

traverse <- function(x, i = numeric(0)) {
    if (is.list(x)) {
        out <- mapply(traverse, x, i = lapply(seq_along(x), function(a) c(i, a)), SIMPLIFY = FALSE)
        # out <- mapply(traverse, x, i = paste(i, seq_along(x), sep = ','), SIMPLIFY = FALSE)
        do.call(rbind, out)
        # do.call(rbind, lapply(x, traverse, i = paste(i, seq_along(x), sep = ':')))
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

modify_epsg <- function(base_crs, new_center = NULL, new_angle = NULL) {
    if (missing(base_crs)) {
        stop('base crs is missing')
    }
    if (is.character(base_crs)) {
        base_crs <- try(st_crs(base_crs), silent = TRUE)
    }
    if (!inherits(base_crs, 'crs')) {
        stop('base crs argument is not a valid crs')
    }
    if (is.null(new_center) && is.null(new_angle)) {
        return(base_crs)
    }
    # parse wkt
    p_wkt <- parse_wkt(base_crs)
    # traverse wkt
    i_wkt <- traverse(p_wkt)
    # false easting/northing -> check epsg code in Table F.3 (https://docs.ogc.org/is/18-010r7/18-010r7.html#106)
    tabf.3 <- list(
        false_easting = c(8806, 8816, 8826),
        false_northing = c(8807, 8817, 8827),
        angle_northing = 8814
        )
    if (!is.null(new_center)) {
        if (!is.numeric(new_center)) {
            stop('argument new_center must be numeric')
        }
        if (length(new_center) != 2) {
            stop('argument new_center must be of length 2')
        }
        # get false easting (x-axis)
        i_x <- get_index(i_wkt, tabf.3[['false_easting']])
        # get false northing (y-axis)
        i_y <- get_index(i_wkt, tabf.3[['false_northing']])
        # replace values
        p_wkt[[i_x]] <- new_center[1]
        p_wkt[[i_y]] <- new_center[2]
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
    deparse_wkt(p_wkt)
}

crs_user3 <- modify_epsg('EPSG:2056', c(0, 0), 60)

## check angle!
rect_ch <- st_sfc(st_polygon(list(cbind(
                c(6e5, 6.1e5, 6.1e5, 6.05e5, 6e5), 
                c(200000, 2e5, 2.1e5, 2.1e5, 2e5)
                ))))
st_crs(rect_ch) <- crs_ch
# transform
rect_a25 <- st_transform(rect_ch, crs = crs_user3)

par(mfrow = c(2, 1))
plot(rect_ch, asp = 1, axes = TRUE)
plot(rect_a25, asp = 1, axes = TRUE)
sqrt(rowSums(apply(st_coordinates(rect_ch)[, 1:2], 2, function(x) diff(x)) ^ 2))
sqrt(rowSums(apply(st_coordinates(rect_a25)[, 1:2], 2, function(x) diff(x)) ^ 2))


### sf blsmodelr -> TODO: change to sf in blsmodelr!!!
poly <- st_polygon(list(cbind(
            x = c(0, 0, 10, 10, 0),
            y = c(0, 10, 10, 0, 0)
            )))
pts <- st_multipoint(cbind(-5:5, -5:5))
st_intersection(pts, poly)
