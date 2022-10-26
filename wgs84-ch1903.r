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

.change_coords <- function(x, y, crs_from = NULL, crs_to = NULL, swap_xy = FALSE) {
    # insist on numeric
    x <- as.numeric(x)
    y <- as.numeric(y)
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
            crs_from <- 'EPSG:21781'
        } else {
            # LV03
            crs_from <- 'EPSG:2056'
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
    crs_to = NULL, swap = FALSE) {
	if (inherits(x, "Sources") && ncol(x) == 4) {
		out <- x
		out[, 2:3] <- .change_coords(x[, 2], x[, 3], 
            crs_from = crs_from, crs_to = crs_to, 
            swap_xy = swap_xy)
	} else if (inherits(x, "Sensors") && ncol(x) >= 7) {
		out <- convert(x)
        out[, c('x-Coord (m)', 'y-Coord (m)')] <- .change_coords(
            y[, 'x-Coord (m)'], y[, 'y-Coord (m)'], 
            crs_from = crs_from, crs_to = crs_to, 
            swap_xy = swap_xy)
	} else {
		if (is.null(y)) {
			if (is.matrix(x)) {
				y <- x[, 2]
				x <- x[, 1]
			} else {
				y <- x[[2]]
				x <- x[[1]]
			}
		}
        out <- .change_coords(x, y, 
            crs_from = crs_from, crs_to = crs_to, 
            swap_xy = swap_xy)
        colnames(out) <- c('x', 'y')
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
wgs_to_map <- function(MyMap, lat, lon = NULL, zoom) {
    require(RgoogleMaps, quietly = TRUE)
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
			if (!is.matrix(lat) & !is.data.frame(lat)) {
				lon <- lat[2]
				lat <- lat[1]
			} else {
				lon <- lat[, 2]
				lat <- lat[, 1]
			}
		}
		out <- LatLon2XY.centered(MyMap, lat, lon, zoom)
		return(cbind(x = out$newX, y = out$newY))
	}
}

ch_to_map <- function(MyMap, x, y = NULL, ...) {
	WGS84 <- ch_to_wgs(x, y)
	wgs_to_map(MyMap, WGS84, ...)
} 

