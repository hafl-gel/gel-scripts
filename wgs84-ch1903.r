## ~~~~ new functions using EPSG ~~~~ ##

# WGS_to_CH
# CH_to_WGS
# WGS_to_map
# CH_to_map
# # map_to_WGS
# # map_to_CH

# check EPSG codes on https://epsg.io
# wgs84: 4326
# ch1903/LV03: 21781
# ch1903/LV95: 2056

library(sf, quietly = TRUE)
st_sfc
st_transform
st_crs

x <- data.frame(lat = 47, lon = 7.5)

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
		dummy <- LatLon2XY.centered(MyMap, lat[, "y-Coord (m) "], lat[, "x-Coord (m) "], zoom)
		out[, c("x-Coord (m) ", "y-Coord (m) ")] <- cbind(dummy$newX, dummy$newY)
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

