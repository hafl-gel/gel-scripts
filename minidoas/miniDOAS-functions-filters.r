
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

filter.functions <- list(
    "Rect" = winRect,
    "Hann" = winHann,
    "Hamming" = winHamming,
    "Blackman" = winBlackman,
    "BmNuttall" = winBmNuttall,
    "Sin" = winSin,
    "Gauss" = winGauss,
    "Kaiser" = winKaiser,
    "DolphChebyshev" = winDolphChebyshev,
    "BmHarris" = winBmHarris,
    "Tukey" = winTukey,
    "Poisson" = winPoisson,
    "Exp" = winExp,
    "ExpHamming" = winExpHamming
)


highpass.filter <- function(dat, DOAS.win, ...){
    filt <- filter.functions[[DOAS.win$filter.type]](DOAS.win$filter.strength, ...)
    dat - filter(dat[DOAS.win$pixel_filter], filt, 'convolution', 2, circular = FALSE)
}
