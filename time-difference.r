#### estimate time lag between x1 and x2 (add time lag to x2)
# library(ibts)
# x <- seq(0, 2*pi, length.out = 200)
# y <- sin(x) + rnorm(200, sd = 0.2)
# tim <- floor_time(Sys.time()) + seq_along(x)
# ind1 <- 20:80
# x1 <- as.ibts(y[ind1], st = tim[ind1], et = tim[ind1] + 1, tz = "UTC")
# ind2 <- ind1 +5
# x2 <- as.ibts(y[ind2], st = tim[ind1], et = tim[ind1] + 1, tz = "UTC")
# xy <- na.omit(x1[x2])
# dT <- checkTimeDiff(xy, dyn = c(-25, 25), Plot = TRUE)
# plot(x1)
# lines(x2, col = "blue", lty = 2)
# x3 <- x2
# st(x3) <- st(x2) + dT
# et(x3) <- et(x2) + dT
# lines(x3, col = "red")
# dT

#### check time shifts:
cov_fun <- function(x1,x2){
  n <- length(x1)
  x1 <- fft(x1)/n
  x2 <- fft(x2)/n
  if(n%%2){
    # n ungerade:
    Re(fft(Conj(x2) * x1, inverse=TRUE))[c(((n+1)/2+1):n,1:((n+1)/2))]*n/(n-1)
  } else {
    # n gerade:
    Re(fft(Conj(x2) * x1, inverse=TRUE))[c((n/2+1):n,1:(n/2))]*n/(n-1)
  }
}

find_dynlag <- function(x,dyn){
  n <- length(x)
  if(n%%2){
    # ungerade
    m <- (n+1)/2
  } else {
    # gerade
    m <- n/2 + 1
  }
  ind <- seq(dyn[1],dyn[2]) + m
  # find max:
  maxis <- ind[which.max(abs(x[ind]))]
  maxis - m
}

checkTimeDiff <- function(xy,dyn = c(-1, 1)*1E4, Plot = FALSE, xlim = dl + c(-1, 1)*50, ...){
  dt_xy <- as.numeric(median(diff(st(xy))), units = "secs") 
  cv <- cov_fun(xy[[1]], xy[[2]])
  n <- length(cv)
  if((n %% 2) == 1) { 
    xx <- seq(-(n - 1)/2,(n - 1)/2,by=dt_xy) * dt_xy
  } else {
    xx <- seq(-n/2 + 1,n/2,by=dt_xy) * dt_xy
  }
  dl <- find_dynlag(cv, dyn)
  if(Plot){
    plot(xx, cv, type = "l", xlim = xlim, ...)
  }
  dl * dt_xy
}
