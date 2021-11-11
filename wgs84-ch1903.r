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
WGS.to.xy <- function(MyMap, lat, lon=NULL, zoom){
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
CH.to.xy <- function(MyMap,x,y=NULL,...){
	WGS84 <- CH.to.WGS(x,y,...)
	WGS.to.xy(MyMap,WGS84)
} 


  # rosavent2 <- function (frec, fnum = 4, fint = 5, flab = 2, ang = 3 * pi/16, 
  #     col = rainbow(10, 0.5, 0.92, start = 0.33, end = 0.2), margen = c(0, 
  #         0, 0, 0), key = TRUE, uni = "m/s",num.rast=4,dRfrac=0.8,polyOnly=FALSE,farb = "black",farbT = farb,cex.uleg=1, ...){
  #     # old.par <- par(no.readonly = TRUE)
  #     # on.exit(par(old.par))
  #     dR <- max(diff(par("usr")[1:2]),diff(par("usr")[3:4]))/2*dRfrac
  #     if(is.table(frec)){
  #       if(length(dim(frec)) > 1){
  #         class(frec) <- "matrix"
  #       } else {
  #         frec <- as.numeric(frec)
  #       }
  #     } 
  #     # if (is.matrix(frec)) 
  #     #     frec <- as.data.frame(frec)
  #     if (is.vector(frec)) {
  #         ndir <- length(frec)
  #         nr <- 1
  #     }
  #     else {
  #         ndir <- length(frec[1, ])
  #         nr <- nrow(frec)
  #     }
  #     fmax <- fnum * fint
  #     tot <- sum(frec)
  #     fr <- 100 * frec/tot
  #     key <- (nr > 1) && key
  #     if (key) 
  #         mlf <- 3
  #     else mlf <- 1
  #     # par(mar = margen, new = FALSE, bg = dev.bg, cex=cex)
  #     fx <- cos(pi/2 - (2 * pi/ndir * 0:(ndir - 1)))
  #     fy <- sin(pi/2 - (2 * pi/ndir * 0:(ndir - 1)))
  #     # plot(fx, fy, xlim = c(-fmax - mlf * fint, fmax + fint), ylim = c(-fmax - 
  #     #     fint, fmax + fint), xaxt = "n", yaxt = "n", xlab = "", 
  #     #     ylab = "", bty = "n", asp = 1, type = "n", ...)
  #     if (nr == 1) {
  #         cx <- fx * fr
  #         cy <- fy * fr
  #     }
  #     else {
  #       f <- colSums(fr)
  #         cx <- fx * f
  #         cy <- fy * f
  #         for (i in nr:2) {
  #           # browser()
  #             f <- f - fr[i, ]
  #             cx <- c(cx, NA, fx * f)
  #             cy <- c(cy, NA, fy * f)
  #         }
  #     }
  #     # browser()
  #     polygon(cx/fmax*dR, cy/fmax*dR, col = col[nr:1])
  #     if(!polyOnly){
  #       symbols(c(0 * 1:fnum), c(0 * 1:fnum), circles = c(fint * 
  #           1:fnum)/fmax*dR, inches = FALSE, add = TRUE,fg=farb)
  #       segments(0 * 1:num.rast, 0 * 1:num.rast, dR * cos(pi/2 - (2 * pi/num.rast * 0:(num.rast - 1))), dR * sin(pi/2 - (2 * pi/num.rast * 0:(num.rast - 1))),col=farb)
  #       fmaxi <- (fmax + fint/4)/fmax*dR
  #       text(0, fmaxi, "N",col=farbT)
  #       text(0, -fmaxi, "S",col=farbT)
  #       text(fmaxi, 0, "E",col=farbT)
  #       text(-fmaxi, 0, "W",col=farbT)
  #       if (flab == 2) 
  #           for (i in 1:fnum) text(i * fint * cos(ang)/fmax*dR, i * fint * 
  #               sin(ang)/fmax*dR, paste(i * fint, "%"),col=farbT)
  #       else if (flab == 1) 
  #           text(dR * cos(ang), dR * sin(ang), paste(fmax, "%"),col=farbT)
  #       if (key) {
  #         # browser()
  #           legend("topright",fill = col, legend = rownames(frec),title=uni,cex=cex.uleg)
  #           # text(-fmaxi - 0.4 * fint/fmax*dR, fmaxi + 0.9 * fint/fmax*dR, uni)
  #       }
  #   }
  #     invisible()
  #   }

