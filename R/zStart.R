.onLoad <- function(libname, pkgname) {
    ### neue Filter Funktionen:
    options(
        md.filter.function.list = list(
            # Blackman-Harris:
            "BmHarris" = function(n, ...){
              N <- n - 1
              x <- 0.35875 - 
                0.48829 * cos(2 * pi * seq.int(0, N) / N) + 
                0.14128 * cos(4 * pi * seq.int(0, N) / N) - 
                0.01168 * cos(6 * pi * seq.int(0, N) / N)
              x / sum(x)
            },
            # Hann:
            "Hann" = function(n, ...){
              N <- n - 1 
              sin(pi * seq.int(0, N) / N) ^ 2 / (0.5 * N)
            },
            # Hamming:
            "Hamming" = function(n, p1 = 25 / 46, p2 = 21 / 46, ...){
              N <- n - 1
              x <- p1 - p2 * cos(2 * pi * seq.int(0, N) / N)
              x / sum(x)
            },
            # Blackman:
            "Blackman" = function(n, p1 = 0.42, p2 = 0.5, p3 = 0.08, ...){
              N <- n - 1
              x <- p1 - 
                p2 * cos(2 * pi * seq.int(0, N) / N) + 
                p3 * cos(4 * pi * seq.int(0, N) / N)
              x/sum(x)
            },
            # Blackman-Nuttall:
            "BmNuttall" = function(n, ...){
              N <- n - 1
              x <- 0.3635819 - 
                0.4891775 * cos(2 * pi * seq.int(0, N) / N) + 
                0.1365995 * cos(4 * pi * seq.int(0, N) / N) - 
                0.0106411 * cos(6 * pi * seq.int(0, N) / N)
              x / sum(x)
            },
            # Rectangular (Moving Average):
            "Rect" = function(n, ...){
              rep(1 / n, n)
            },
            # Exponential:
            "Exp" = function(n, p1 = n / 2, ...){
              N <- (n - 1) / 2
              x <- exp(p1 * sqrt(1 - (seq.int(-N, N) / N) ^ 2)) / exp(p1)
              x / sum(x)
            },
            # Poisson:
            "Poisson" = function(n, p1 = n / 2 / 6.9, ...){
              N <- n - 1
              x <- exp(-abs(seq.int(0, N) - N / 2) / p1)
              x / sum(x)
            }
        )
    )
    # set minidoas options (except filter functions)
    options(
        # # old evaluation (using BmHarris)
        # md.filter.type = 'BmHarris',
        # md.filter.strength = 25,
        # md.filter.window = c(202.2, 230.6),
        # md.fit.window = c(204.5, 228.4)
        # new evaluation (using Exp filter)
        md.filter.type = 'Exp',
        md.filter.strength = 51,
        md.filter.window = c(200, 232),
        md.fit.window = c(204.4, 227.7)
    )
}
