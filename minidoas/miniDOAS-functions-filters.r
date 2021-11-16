library(IDPmisc)

### neue Filter Funktionen:

# Blackman-Harris Windowing:
winBmHarris <- function(n, ...){
  N <- n - 1
  x <- 0.35875 - 
    0.48829 * cos(2 * pi * seq.int(0, N) / N) + 
    0.14128 * cos(4 * pi * seq.int(0, N) / N) - 
    0.01168 * cos(6 * pi * seq.int(0, N) / N)
  x / sum(x)
}

# Blackman-Nuttall Windowing:
winBmNuttall <- function(n, ...){
  N <- n - 1
  x <- 0.3635819 - 
    0.4891775 * cos(2 * pi * seq.int(0, N) / N) + 
    0.1365995 * cos(4 * pi * seq.int(0, N) / N) - 
    0.0106411 * cos(6 * pi * seq.int(0, N) / N)
  x / sum(x)
}

winTukey <- function(n, p1 = 0.5, ...){
  stopifnot(0 < p1 & 1 >= p1)
  N <- n - 1
  n1 <- seq.int(0, p1 * N / 2)
  n2 <- rev(seq.int(N, N * (1 - p1 / 2)))
  x <- c(
    0.5 * (1 + cos(pi * (2 * n1 / (p1 * N) - 1))),
    rep(1, length(seq.int(0, N)[-c(n1, n2)])),
    0.5 * (1 + cos(pi * (2 * n2 / (p1 * N) - 2 / p1 + 1)))
    )
  x / sum(x)
}

winPoisson <- function(n, p1 = n / 2 / 6.9, ...){
  N <- n - 1
  x <- exp(-abs(seq.int(0, N) - N / 2) / p1)
  x / sum(x)
}

winExp <- function(n, p1 = n / 2, ...){
  N <- (n - 1) / 2
  x <- exp(p1 * sqrt(1 - (seq.int(-N, N) / N) ^ 2)) / exp(p1)
  x / sum(x)
}

winExpHamming <- function(n, p1 = n / 2, ...){
  N <- (n - 1) / 2
  x <- 0.5 * exp(p1 * sqrt(1 - (seq.int(-N, N) / N) ^ 2)) / exp(p1) + 
    0.27 - 0.23 * cos(2 * pi * (seq.int(-N, N) / N / 2 + 0.5))
  x / sum(x)
}

chebwin <- function(n, at){
    if (!(length(n) == 1 && (n == round(n)) && (n > 0))) 
        stop("n has to be a positive integer")
    if (!(length(at) == 1 && (at == Re(at)))) 
        stop("at has to be a real scalar")
    if (n == 1) 
        w <- 1
    else {
        gamma <- 10 ^ (-at / 20)
        beta <- cosh(1 / (n - 1) * acosh(1 / gamma))
        k <- 0:(n - 1)
        x <- beta * cos(pi * k / n)
        p <- cheb(n - 1, x)
        if (n %% 2) {
            w <- Re(fft(p))
            M <- (n + 1) / 2
            w <- w[1:M] / w[1]
            w <- c(w[M:2], w)
        }
        else {
            p <- p * exp((0+1i) * pi / n * (0:(n - 1)))
            w <- Re(fft(p))
            M <- n / 2 + 1
            w <- w / w[2]
            w <- c(w[M:2], w[2:M])
        }
    }
    w
}

cheb <- function(n, x){
    if (!(is.numeric(n) && (n == round(n)) && (n >= 0))) 
        stop("n has to be a positive integer")
    Tt <- numeric(length(x))
    ind <- x <= 1
    if (any(ind)) 
        Tt[ind] <- cos(n * acos(as.complex(x[ind])))
    ind <- x > 1
    myacosh <- function(x) log(x + sqrt(x ^ 2 - 1))
    if (any(ind)) 
        T[ind] <- cosh(n * myacosh(as.complex(x[ind])))
    Re(Tt)
}

# Rectangular Windowing (Moving Average):
winRect <- function(n, ...){
  rep(1 / n, n)
}
# Cosine/Sine Windowing:
winSin <- function(n, ...){
  N <- n - 1 
  x <- sin(pi * seq.int(0, N) / N)
  x / sum(x)
}
# Gaussian Windowing:
winGauss <- function(n, p1 = 0.4, ...){
  N <- n - 1 
  x <- exp(-0.5 * ((seq.int(0, N) - N / 2) / (N / 2 * p1)) ^ 2)
  x / sum(x)
}
# Hann Windowing:
winHann <- function(n, ...){
  N <- n - 1 
  sin(pi * seq.int(0, N) / N) ^ 2 / (0.5 * N)
}

# Hamming Windowing (p1=25/46,p2=21/46)/(p1=0.54,p2=0.46)/(p1=0.53836,p2=0.46164):
winHamming <- function(n, p1 = 25 / 46, p2 = 21 / 46, ...){
  N <- n - 1
  x <- p1 - p2 * cos(2 * pi * seq.int(0, N) / N)
  x / sum(x)
}

# Blackman Windowing (p1=0.42,p2=0.5,p3=0.08)/(p1=7938/18608,p2=9240/18608,p3=1430/18608):
winBlackman <- function(n, p1 = 0.42, p2 = 0.5, p3 = 0.08, ...){
  N <- n - 1
  x <- p1 - 
    p2 * cos(2 * pi * seq.int(0, N) / N) + 
    p3 * cos(4 * pi * seq.int(0, N) / N)
  x/sum(x)
}

# Flat-Top Windowing:
winFlatTop <- function(n, p0 = 0.21557895, p1 = 0.41663158, 
  p2 = 0.277263158, p3 = 0.083578947, p4 = 0.006947368){
  N <- n - 1
  x <- p0 - 
    p1 * cos(2 * pi * seq.int(0, N) / N) + 
    p2 * cos(4 * pi * seq.int(0, N) / N) - 
    p3 * cos(6 * pi * seq.int(0, N) / N) + 
    p4 * cos(8 * pi * seq.int(0, N) / N)
  x/sum(x)
}

# Dolph-Chebyshev Windowing:
winDolphChebyshev <- function(n, p1 = 100, ...){
  x <- chebwin(n, p1)
  x / sum(x)
}

# Kaiser Windowing alpha = 2:
winKaiser <- function(n, p1 = 2 * pi, ...){
  x <- kaiser(n, p1)
  x / sum(x)
}

kaiser <- function(n, beta){
    if (!(length(n) == 1 && (n == round(n)) && (n > 0))) 
        stop("kaiser:  n has to be a positive integer")
    if (!(length(beta) == 1 && (beta == as.double(beta)))) 
        stop("kaiser:  beta has to be a real scalar")
    if (n == 1) 
        w = 1
    else {
        m = n - 1
        k = 0:m
        k = 2 * beta / m * sqrt(k * (m - k))
        w = besselI(k, 0) / besselI(beta, 0)
    }
    w
}

double.filter <- function(x.dat, filter.strength, winFUN, double = FALSE, ...){
  filt <- winFUN(filter.strength, ...)
  if(double){
    x.dat - rev(filter(
      rev(filter(x.dat, filt, "convolution", 2, circular = FALSE)), 
      filt, "convolution", 2, circular = FALSE))
  } else {
    x.dat - filter(x.dat, filt, "convolution", 2, circular = FALSE)
  }
}

### low pass filter based on a combination of loess (gaussian) and rfbaseline function to fit a baseline through the high-pass filtered spectrum
### ****************************************************************************** 
special.filter <- function(x.dat, filter.strength_rebs = 0.175, 
  filter.strength_loess = 0.2, b, Scale, delta, fam, maxit = c(10, 0)) {  
  loess(
    -IDPmisc::rfbaseline(x = seq_along(x.dat), y = -x.dat, 
      span = filter.strength_rebs, maxit = maxit, b = b, Scale = Scale,
      delta = delta)$fit ~ seq_along(x.dat), 
    span = filter.strength_loess, family = fam)$fitted
}


### low pass filter based on R's 'loess' function (local polynomial fitting)
### ****************************************************************************** 
loess.filter <- function(x.dat, filter.strength, fam){
  loess(x.dat ~ seq_along(x.dat), span = filter.strength, family = fam)$fitted
}


### low pass filter based on R-package IDPmisc's 'rfbaseline' function (asymmetric local polynomial fitting)
### ****************************************************************************** 
rfbaseline.filter <- function(x.dat, filter.strength, b, Scale, delta,
  maxit = c(10, 0)) {
  -IDPmisc::rfbaseline(x = seq_along(x.dat), y = -x.dat, 
    span = filter.strength, maxit = maxit, b = b, Scale = Scale, 
    delta = delta)$fit
}



### high pass filter (for broadband absorption exclusion) (simply makes use of "lowpass.filter"; y.dat = data to be filtered, x.data = data to be filtered against - in case of moving averages x.dat should be set = y.dat -, ...)
### ****************************************************************************** 
highpass.filter <- function(dat, DOAS.win, ...){
  dat <- dat[DOAS.win$pixel1]
  switch(DOAS.win$filter.type,
    "Rect" = double.filter(dat,DOAS.win$filter.strength, winRect, DOAS.win$double, ...),
    "Hann" = double.filter(dat,DOAS.win$filter.strength, winHann, DOAS.win$double, ...),
    "Hamming" = double.filter(dat,DOAS.win$filter.strength, winHamming, DOAS.win$double, ...),
    "Blackman" = double.filter(dat,DOAS.win$filter.strength, winBlackman, DOAS.win$double, ...),
    "BmNuttall" = double.filter(dat,DOAS.win$filter.strength, winBmNuttall, DOAS.win$double, ...),
    "FlatTop" = double.filter(dat,DOAS.win$filter.strength, winFlatTop, DOAS.win$double, ...),
    "Sin" = double.filter(dat,DOAS.win$filter.strength, winSin, DOAS.win$double, ...),
    "Gauss" = double.filter(dat,DOAS.win$filter.strength, winGauss, DOAS.win$double, ...),
    "Kaiser" = double.filter(dat,DOAS.win$filter.strength, winKaiser, DOAS.win$double, ...),
    "DolphChebyshev" = double.filter(dat,DOAS.win$filter.strength, winDolphChebyshev, DOAS.win$double, ...),
    "BmHarris" = double.filter(dat,DOAS.win$filter.strength, winBmHarris, DOAS.win$double, ...),
    "Tukey" = double.filter(dat,DOAS.win$filter.strength, winTukey, DOAS.win$double, ...),
    "Poisson" = double.filter(dat,DOAS.win$filter.strength, winPoisson, DOAS.win$double, ...),
    "Exp" = double.filter(dat,DOAS.win$filter.strength, winExp, DOAS.win$double, ...),
    "ExpHamming" = double.filter(dat,DOAS.win$filter.strength, winExpHamming, DOAS.win$double, ...),
    "loess" = dat - loess(dat ~ seq_along(dat), span = DOAS.win$filter.strength, family = DOAS.win$special.control$fam)$fitted,
    "special" = {
      x <- dat - special.filter(dat, filter.strength_rebs = DOAS.win$filter.strength, 
        filter.strength_loess = DOAS.win$special.control$filter.strength_loess, 
        b = DOAS.win$special.control$b, Scale = DOAS.win$special.control$Scale,
        delta = DOAS.win$special.control$delta, fam = DOAS.win$special.control$fam,
        maxit = DOAS.win$special.control$maxit)
      x - special.filter(x, filter.strength_rebs = DOAS.win$filter.strength * 
        DOAS.win$special.control$filter.strength_multiplication, 
        filter.strength_loess = DOAS.win$special.control$filter.strength_loess,
        b = DOAS.win$special.control$b, Scale = DOAS.win$special.control$Scale,
        delta = DOAS.win$special.control$delta, fam = DOAS.win$special.control$fam,
        maxit = DOAS.win$special.control$maxit) # Joerg
    }
  )
}


# ######## testing:
# fsAVG <- 25
# fs <- 0.25
# DOAS.win <- getWindows(DOAS.info, "special", timerange, 
# straylight.window, c(202,229.1), fit.window, 
# filter.strength, tau.shift)
# # DOAS.win3 <- getWindows(DOAS.info, "specialAVG", timerange, 
# # DOAS.win3 <- getWindows(DOAS.info, "specialMED", timerange, 
# # DOAS.win3 <- getWindows(DOAS.info, "specialHamming", timerange, 
# straylight.window, c(202,229.1), fit.window, 
# fs, tau.shift)