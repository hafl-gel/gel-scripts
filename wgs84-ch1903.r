# options(digits=20)
## The MIT License (MIT)
## 
## Copyright (c) 2014 Federal Office of Topography swisstopo, Wabern, CH
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
##  of this software and associated documentation files (the "Software"), to deal
##  in the Software without restriction, including without limitation the rights
##  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##  copies of the Software, and to permit persons to whom the Software is
##  furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in
##  all copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
##  THE SOFTWARE.
## 
## 
## Source: http://www.swisstopo.admin.ch/internet/swisstopo/en/home/topics/survey/sys/refsys/projections.html (see PDFs under "Documentation")
## Updated 9 dec 2014
## Please validate your results with NAVREF on-line service: http://www.swisstopo.admin.ch/internet/swisstopo/en/home/apps/calc/navref.html (difference ~ 1-2m)

## Convert WGS lat/long (° dec) to CH y
WGS.to.CH.y <- function(lat, lng){
	## Converts decimal degrees to sexagesimal seconds
	lat <- DEC.to.SEX(lat)
	lng <- DEC.to.SEX(lng)

	## Auxiliary values (% Bern)
	lat_aux <- (lat - 169028.66)/10000
	lng_aux <- (lng - 26782.5)/10000
  
	## Process Y
 	y <- {600072.37 +
	211455.93 * lng_aux -
	10938.51 * lng_aux * lat_aux -
	0.36 * lng_aux * (lat_aux^2) -
	44.54 * (lng_aux^3)}
     
	return(y)
}

## Convert WGS lat/long (° dec) to CH x
WGS.to.CH.x <- function(lat, lng){

	## Converts decimal degrees to sexagesimal seconds
	lat <- DEC.to.SEX(lat)
	lng <- DEC.to.SEX(lng)
  
	## Auxiliary values (% Bern)
	lat_aux <- (lat - 169028.66)/10000
	lng_aux <- (lng - 26782.5)/10000

	## Process X
  	x <- {200147.07 +
	308807.95 * lat_aux + 
	3745.25 * (lng_aux^2) +
	76.63 * (lat_aux^2) -
	194.56 * (lng_aux^2) * lat_aux +
	119.79 * (lat_aux^3)}

	return(x)
}


## Convert CH y/x to WGS lat
CH.to.WGS.lat <- function (y, x){

	## Converts military to civil and  to unit = 1000km
	## Auxiliary values (% Bern)
	y_aux <- (y - 600000)/1000000
	x_aux <- (x - 200000)/1000000
  
	## Process lat
	lat <- {16.9023892 +
	3.238272 * x_aux -
	0.270978 * (y_aux^2) -
	0.002528 * (x_aux^2) -
	0.0447   * (y_aux^2) * x_aux -
	0.0140   * (x_aux^3)}
    
	## Unit 10000" to 1 " and converts seconds to degrees (dec)
	lat <- lat * 100/36
  
  	return(lat)  
}

## Convert CH y/x to WGS long
CH.to.WGS.lng <- function (y, x){

	## Converts military to civil and  to unit = 1000km
	## Auxiliary values (% Bern)
	y_aux <- (y - 600000)/1000000
	x_aux <- (x - 200000)/1000000
  
	## Process long
	lng <- {2.6779094 +
	4.728982 * y_aux +
	0.791484 * y_aux * x_aux +
	0.1306   * y_aux * (x_aux^2) -
	0.0436   * (y_aux^3)}
     
	## Unit 10000" to 1 " and converts seconds to degrees (dec)
  	lng <- lng * 100/36

	return(lng)
}


## Convert decimal degrees to sexagesimal seconds
DEC.to.SEX <- function(angle){

	## Extract DMS
	angle_chr <- as.character(angle)
	deg <- as.numeric(strsplit(angle_chr, "\\.")[[1]][1])
	min <- as.numeric(strsplit(as.character((angle-deg)*60), "\\.")[[1]][1])
	sec <- (((angle-deg)*60) - min) * 60

	## Result in seconds
  	return(sec + min*60 + deg*3600)

}


## 
WGS.to.CH <- function(lat,lng=NULL){
	if(inherits(lat,"Sources")&&ncol(lat)==4){
		out <- lat
		out[,2:3] <- cbind(WGS.to.CH.x(lat[,3],lat[,2]),WGS.to.CH.y(lat[,3],lat[,2]))
	} else if(inherits(lat, "Sensors") && ncol(lat) >= 7){
		out <- convert(lat)
		out[, c("x-Coord (m)", "y-Coord (m)")] <- cbind(
      WGS.to.CH.x(lat[, "y-Coord (m)"], lat[, "x-Coord (m)"]),
      WGS.to.CH.y(lat[, "y-Coord (m)"], lat[, "x-Coord (m)"])
      )
	} else {
		if(is.null(lng)){
			if(!is.matrix(lat)&!is.data.frame(lat)){
				lng <- lat[2]
				lat <- lat[1]
			} else {
				lng <- lat[,2]
				lat <- lat[,1]
			}
		}
		out <- cbind(y=WGS.to.CH.y(lat,lng),x=WGS.to.CH.x(lat,lng))
	}
	return(out)
}

## 
CH.to.WGS <- function(x,y=NULL){
	if(inherits(x,"Sources")&&ncol(x)==4){
		out <- x
		out[,2:3] <- cbind(CH.to.WGS.lng(x[,2],x[,3]),CH.to.WGS.lat(x[,2],x[,3]))
	} else if(inherits(x,"Sensors") && ncol(x) >= 7){
		out <- convert(x)
		out[, c("x-Coord (m)", "y-Coord (m)")] <- cbind(
      CH.to.WGS.lng(x[, "x-Coord (m)"],x[, "y-Coord (m)"]),
      CH.to.WGS.lat(x[, "x-Coord (m)"],x[, "y-Coord (m)"])
      )
	} else {
		if(is.null(y)){
			if(!is.matrix(x)&!is.data.frame(x)){
				y <- x[2]
				x <- x[1]
			} else {
				y <- x[,2]
				x <- x[,1]
			}
		}
		out <- cbind(lat=CH.to.WGS.lat(x,y),lng=CH.to.WGS.lng(x,y))
	}
	return(out)
}


##
WGS.to.map <- function(MyMap, lat, lon=NULL, zoom){
    require(RgoogleMaps, quietly = TRUE)
	if(inherits(lat,"Sources")&&ncol(lat)==4){
		out <- lat
		dummy <- LatLon2XY.centered(MyMap, lat[,3], lat[,2], zoom)
		out[,2:3] <- cbind(dummy$newX,dummy$newY)
		return(out)
	} else if(inherits(lat,"Sensors") && ncol(lat) >= 7){
		out <- convert(lat)
		dummy <- LatLon2XY.centered(MyMap, lat[, "y-Coord (m)"], lat[, "x-Coord (m)"], zoom)
		out[, c("x-Coord (m)", "y-Coord (m)")] <- cbind(dummy$newX, dummy$newY)
		return(out)
	} else {
		if(is.null(lon)){
			if(!is.matrix(lat)&!is.data.frame(lat)){
				lon <- lat[2]
				lat <- lat[1]
			} else {
				lon <- lat[,2]
				lat <- lat[,1]
			}
		}
		out <- LatLon2XY.centered(MyMap, lat, lon, zoom)
		return(cbind(x=out$newX,y=out$newY))
	}
}

##
CH.to.map <- function(MyMap,x,y=NULL,...){
	WGS84 <- CH.to.WGS(x,y,...)
	WGS.to.map(MyMap,WGS84)
} 

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
wgs_to_map <- WGS.to.map
ch_to_map <- function(MyMap, x, y = NULL, ...) {
	WGS84 <- ch_to_wgs(x, y)
	wgs_to_map(MyMap, WGS84, ...)
} 

