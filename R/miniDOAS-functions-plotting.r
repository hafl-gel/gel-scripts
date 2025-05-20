
### function to draw (and annotate) (deca-)logarithmic ticks on R standard plots (plot() has to be plotted with x/yaxt="n" beforehand) (input: axis.loc = axis location: 1,2,3 or 4 for c(bottom, left, top, right), axis.values = vector of the data plotted on respective axis, tick.dim = c(major tick length, minor tick length) - eg. c(-0.5,-0.25), cex.labels = cex for axis labels)
### ****************************************************************************** 
log10.ticks <- function(axis.loc, axis.values, tick.dim, cex.labels) {
	axis.values <- ifelse(axis.values <= 0, NA, axis.values)
	val.range <- floor(log10(range(axis.values, na.rm=TRUE, finite=TRUE)))
	pow <- seq(val.range[1],val.range[2]+1)
	ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
	axis(axis.loc,pow,labels=10^pow,tcl=tick.dim[1],cex.axis=cex.labels)
	axis(axis.loc,log10(ticksat),labels=NA,tcl=tick.dim[2],lwd=0,lwd.ticks=1)
}



### plot raw reference (with and without cuvette), dark, and calibration spectra
### ****************************************************************************** 
plot.calibration.spectra <- function(calref.cols, wavelengths, lamp.reference, dark, lamp.N2, dark.N2, NH3, SO2, NO) {
	cex.annotations <- 1.25
	y.limes <- log10(range(c(lamp.reference,dark,dark.N2,lamp.N2,NO,SO2,NH3), na.rm=TRUE))
	plot(wavelengths, log10(lamp.reference), type="l", lty=2, lwd=1, yaxt="n", xlab="wavelength [nm]", ylab="averaged signal [counts]", main="miniDOAS reference, calibration, & dark spectra", cex.axis=cex.annotations, cex.lab=cex.annotations, ylim=y.limes)
	lines(wavelengths, log10(dark), lty=1, lwd=1, col="gray30")
	lines(wavelengths, log10(dark.N2), lty=2, lwd=1, col="gray50")
	lines(wavelengths, log10(lamp.N2), lty=1, lwd=1, col=calref.cols[4])
	lines(wavelengths, log10(NO), lty=1, lwd=1, col=calref.cols[3])
	lines(wavelengths, log10(SO2), lty=1, lwd=1, col=calref.cols[2])
	lines(wavelengths, log10(NH3), lty=1, lwd=1.5, col=calref.cols[1])
	log10.ticks(2, c(lamp.reference,dark,dark.N2,lamp.N2,SO2,NO,NH3), c(-0.5,-0.25), cex.annotations)
	legend("topleft", legend=c(expression(italic("reference")),expression(italic("dark")),expression(italic(paste(N[2]," dark",sep=""))),expression(italic(N[2])),expression(italic("NO")),expression(italic(paste("SO"[2],sep=""))),expression(italic(paste("NH"[3],sep="")))), bty="n", col=c("black","gray50","gray80",rev(calref.cols)), lty=c(2,1,2,1,1,1,1), cex=0.75, ncol=2)
}



### plot calibration DOAS curves
### ****************************************************************************** 
plot.calibration.DOAScurves <- function(calref.cols, wavelengths1, wavelengths2, NH3.dc, SO2.dc, NO.dc) {
	cex.annotations <- 1.25
	par(mar=c(5,4,2,2)+0.6)
	plot(wavelengths1, wavelengths1, type="n", ylim=range(c(NH3.dc,SO2.dc,NO.dc)), xlab="wavelength [nm]", ylab=expression(paste("differential absorption cross section [",m^2,"/",mu,"g]",sep="")), main="calibration DOAS curves", cex.axis=cex.annotations, cex.lab=cex.annotations)
	abline(h=0,lty=2,col="gray60")
	lines(wavelengths2, NO.dc, lty=1, lwd=1.5, col=calref.cols[3])
	lines(wavelengths2, SO2.dc, lty=1, lwd=1.5, col=calref.cols[2])
	lines(wavelengths2, NH3.dc, lty=1, lwd=1.5, col=calref.cols[1])
	legend("bottomright", legend=c(expression(italic(paste("NH"[3],sep=""))),expression(italic(paste("SO"[2],sep=""))),expression(italic("NO"))), lty=c(1,1,1), cex=0.75, col=calref.cols[1:3], bty="n")
}	


### plot calibration diffspec
### ****************************************************************************** 
plot.calibration.diffspecs <- function(calref.cols, wavelengths1, NH3.diffspec, NO.diffspec, SO2.diffspec, NH3.dc, NO.dc, SO2.dc) {
	cex.annotations <- 1.25
	par(mar=c(5,4,2,2)+0.6)
	# plot(wavelengths1, wavelengths1, type="n", ylim=c(0,max(c(NH3.diffspec,NO.diffspec,SO2.diffspec))), xlab="wavelength [nm]", ylab="I/I0", main="calibration diffspec", cex.axis=cex.annotations, cex.lab=cex.annotations)
  plot(wavelengths1, wavelengths1, type="n", ylim=range(c(NH3.diffspec,NO.diffspec,SO2.diffspec)), xlab="wavelength [nm]", ylab="I/I0", main="calibration diffspec", cex.axis=cex.annotations, cex.lab=cex.annotations)
	abline(h=1,lty=2,col="gray60")
  lines(wavelengths1, SO2.diffspec/exp(SO2.dc), lty=1, lwd=1.5, col="black")
  lines(wavelengths1, SO2.diffspec/exp(SO2.dc), lty=2, col=calref.cols[2])
	lines(wavelengths1, SO2.diffspec, lty=1, lwd=1.5, col=calref.cols[2])
  lines(wavelengths1, NO.diffspec/exp(NO.dc), lty=1, lwd=1.5, col="black")
  lines(wavelengths1, NO.diffspec/exp(NO.dc), lty=2, col=calref.cols[3])
	lines(wavelengths1, NO.diffspec, lty=1, lwd=1.5, col=calref.cols[3])
  lines(wavelengths1, NH3.diffspec/exp(NH3.dc), lty=1, lwd=1.5, col="black")
  lines(wavelengths1, NH3.diffspec/exp(NH3.dc), lty=2, col=calref.cols[1])
	lines(wavelengths1, NH3.diffspec, lty=1, lwd=1.5, col=calref.cols[1])
	legend("bottomright", legend=c(expression(italic(paste("NH"[3],sep=""))),expression(italic(paste("SO"[2],sep=""))),expression(italic("NO"))), lty=c(1,1,1), cex=0.75, col=calref.cols[1:3], bty="n")
}	

