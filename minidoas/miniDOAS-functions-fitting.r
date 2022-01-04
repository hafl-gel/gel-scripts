

### multiple linear least square fit (see Stutz & Platt, 1996, Applied Optics 30):
### fit calibration curves to measured DOAS curve, for parallel computing...
### determine best fit after shifting over given tau range and optionally consider fixed (residual) pattern beforehand 
### ******************************************************************************
fit.curves.ARIMA2 <- function(meas.doascurve, x4, dyn.fixed.pattern, Xreg, fit.weights, tau.shift, order2, path.length) {
  aicc <- numeric(length(tau.shift))
  for(i in seq_along(tau.shift)){
    aicc[i] <- AICc(lm(as.numeric(meas.doascurve[x4 + tau.shift[i]] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights))
  }
  tau.best <- tau.shift[which.min(aicc)]
  delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[x4] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights)) - min(aicc)
  fitcurves <- apply(order2,1, function(x){
    out <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),xreg=Xreg, order=x),silent=TRUE)
    if(inherits(out,"try-error"))out <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),xreg=Xreg, order=x, method="ML"),silent=TRUE)
      return(out)
    })
  index <- which(!sapply(fitcurves,function(x)inherits(x,"try-error")))
  nams <- apply(order2,1,paste,collapse="/")[index]
  if(length(index)){
    fitcurves <- fitcurves[index]
    aicc <- sapply(fitcurves,"[[","aicc")
    ord <- order(aicc)
    best.aicc <- aicc[ord[1]]
    deltaAICc <- aicc - best.aicc
    cofs <- t(sapply(fitcurves,function(x)coef(x)[c("NH3","SO2","NO")]))
    mll <- exp(-0.5*deltaAICc)
    wts <- mll/sum(mll)
    aicctab <- data.frame(Order=nams,cofs/path.length,AICc=aicc,delta=deltaAICc,modLik=mll,wts=wts)[ord,]

    best.order <- nams[ord[1]]
    Var <- sapply(fitcurves,function(x){
      diag(x$var.coef)[c("NH3","SO2","NO")]
    })
    coeffs <- c(NH3=sum(cofs[,1]*wts),SO2=sum(cofs[,2]*wts),NO=sum(cofs[,3]*wts))
    se <- sqrt(sapply(1:3,function(x){sum(wts*(Var[x,]+(coeffs[x]-cofs[,x])^2))}))

    ### for record fitted.doascurve.best over all averaging intervals
    fitted.doascurve <- colSums(coeffs*t(Xreg))
    ### corresponding residual spectrum
    residual.best <- meas.doascurve[x4 + tau.best] - dyn.fixed.pattern - fitted.doascurve
    out <- list(tau.best, delta.AICc.zero, best.order, aicctab, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
  } else {
    out <- NULL
  }
  return(out)
} 

fit.curves.ARIMA2.rob <- function(meas.doascurve, x4, dyn.fixed.pattern, Xreg, fit.weights, tau.shift, order2, path.length) {
  require(robustbase)
  aicc <- numeric(length(tau.shift))
  for(i in seq_along(tau.shift)){
    aicc[i] <- AICc(lm(as.numeric(meas.doascurve[x4 + tau.shift[i]] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights))
  }
  tau.best <- tau.shift[which.min(aicc)]
  delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[x4] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights)) - min(aicc)
  # fitcurves <- lmrob(as.numeric(meas.doascurve[x4 + tau.best] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], method="MM")
  fitcurves <- lmrob(as.numeric(meas.doascurve[x4 + tau.best] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights, setting="KS2014")
  res <- meas.doascurve[x4 + tau.best] - colSums(coef(fitcurves)[2:4]*t(Xreg))
  fitres <- apply(order2,1, function(x){
    out <- try(forecastArima(res, order=x),silent=TRUE)
    if(inherits(out,"try-error"))out <- try(forecastArima(res, order=x, method="ML"),silent=TRUE)
      return(out)
    })
  fitcurves <- lapply(seq(nrow(order2)),function(x){
    out <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),transform.pars = FALSE,xreg=Xreg, order=as.numeric(order2[x,]),fixed=c(coef(fitres[[x]]),rep(NA,3))),silent=TRUE)
    if(inherits(out,"try-error"))out <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),transform.pars = FALSE,xreg=Xreg, order=as.numeric(order2[x,]),fixed=c(coef(fitres[[x]]),rep(NA,3)), method="ML"),silent=TRUE)
      return(out)
    })
  index <- which(!sapply(fitcurves,function(x)inherits(x,"try-error")))
  nams <- apply(order2,1,paste,collapse="/")[index]
  if(length(index)){
    fitcurves <- fitcurves[index]
    aicc <- sapply(fitcurves,"[[","aicc")
    ord <- order(aicc)
    best.aicc <- aicc[ord[1]]
    deltaAICc <- aicc - best.aicc
    cofs <- t(sapply(fitcurves,function(x)coef(x)[c("NH3","SO2","NO")]))
    mll <- exp(-0.5*deltaAICc)
    wts <- mll/sum(mll)
    aicctab <- data.frame(Order=nams,cofs/path.length,AICc=aicc,delta=deltaAICc,modLik=mll,wts=wts)[ord,]

    best.order <- nams[ord[1]]
    Var <- sapply(fitcurves,function(x){
      diag(x$var.coef)[c("NH3","SO2","NO")]
    })
    coeffs <- c(NH3=sum(cofs[,1]*wts),SO2=sum(cofs[,2]*wts),NO=sum(cofs[,3]*wts))
    se <- sqrt(sapply(1:3,function(x){sum(wts*(Var[x,]+(coeffs[x]-cofs[,x])^2))}))

    ### for record fitted.doascurve.best over all averaging intervals
    fitted.doascurve <- colSums(coeffs*t(Xreg))
    ### corresponding residual spectrum
    residual.best <- meas.doascurve[x4 + tau.best] - dyn.fixed.pattern - fitted.doascurve
    out <- list(tau.best, delta.AICc.zero, best.order, aicctab, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
  } else {
    out <- NULL
  }
  return(out)
} 

fit.curves.ARIMA1 <- function(meas.doascurve, x4, dyn.fixed.pattern, Xreg, fit.weights, tau.shift, order2, path.length) {
  aicc <- numeric(length(tau.shift))
  for(i in seq_along(tau.shift)){
    aicc[i] <- AICc(lm(as.numeric(meas.doascurve[x4 + tau.shift[i]] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights))
  }
  tau.best <- tau.shift[which.min(aicc)]
  delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[x4] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights)) - min(aicc)

  fitcurves <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),xreg=Xreg, order=as.numeric(order2)),silent=TRUE)
  if(inherits(fitcurves,"try-error"))fitcurves <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),xreg=Xreg, order=as.numeric(order2), method="ML"),silent=TRUE)

  if(!inherits(fitcurves,"try-error")){
    best.order <- paste(order2,collapse="/")
    aicctab <- NULL
    coeffs <- coefficients(fitcurves)[c("NH3","SO2","NO")]
    se <- sqrt(diag(fitcurves$var.coef)[c("NH3","SO2","NO")])
    ### for record fitted.doascurve.best over all averaging intervals
    fitted.doascurve <- colSums(coeffs*t(Xreg))
    ### corresponding residual spectrum
    residual.best <- meas.doascurve[x4 + tau.best] - dyn.fixed.pattern - fitted.doascurve
    out <- list(tau.best, delta.AICc.zero, best.order, aicctab, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
  } else {
    out <- NULL
  }
  return(out)
}

fit.curves.ARIMA1.rob <- function(meas.doascurve, x4, dyn.fixed.pattern, Xreg, fit.weights, tau.shift, order2, path.length) {
  aicc <- numeric(length(tau.shift))
  for(i in seq_along(tau.shift)){
    aicc[i] <- AICc(lm(as.numeric(meas.doascurve[x4 + tau.shift[i]] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights))
  }
  tau.best <- tau.shift[which.min(aicc)]
  delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[x4] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights)) - min(aicc)
  fitcurves <- lmrob(as.numeric(meas.doascurve[x4 + tau.best] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights, setting="KS2014")
  res <- meas.doascurve[x4 + tau.best] - colSums(coef(fitcurves)[2:4]*t(Xreg))
  fitres <- try(forecastArima(res, order=order2),silent=TRUE)
  if(inherits(fitres,"try-error"))fitres <- try(forecastArima(res, order=order2, method="ML"),silent=TRUE)

  fitcurves <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),transform.pars = FALSE,xreg=Xreg, order=as.numeric(order2),fixed=c(coef(fitres),rep(NA,3))),silent=TRUE)
  if(inherits(fitcurves,"try-error"))fitcurves <- try(forecastArima((meas.doascurve[x4 + tau.best] - dyn.fixed.pattern),transform.pars = FALSE,xreg=Xreg, order=as.numeric(order2),fixed=c(coef(fitres),rep(NA,3)), method="ML"),silent=TRUE)

  if(!inherits(fitcurves,"try-error")){
    best.order <- paste(order2,collapse="/")
    aicctab <- NULL
    coeffs <- coefficients(fitcurves)[c("NH3","SO2","NO")]
    se <- sqrt(diag(fitcurves$var.coef)[c("NH3","SO2","NO")])
    ### for record fitted.doascurve.best over all averaging intervals
    fitted.doascurve <- colSums(coeffs*t(Xreg))
    ### corresponding residual spectrum
    residual.best <- meas.doascurve[x4 + tau.best] - dyn.fixed.pattern - fitted.doascurve
    out <- list(tau.best, delta.AICc.zero, best.order, aicctab, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
  } else {
    out <- NULL
  }
  return(out)
}

fit.curves.OLS <- function(meas.doascurve, x4, dyn.fixed.pattern, Xreg, fit.weights, tau.shift, order2, path.length){
  fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[x4 + x] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], weights=fit.weights, model=FALSE))
  aicc <- sapply(fitcurves,AICc)
  index <- which.min(aicc)
  tau.best <- tau.shift[index]
  delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[x4] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights)) - aicc[index]

  fitcurves <- fitcurves[[index]]
  best.order <- aicctab <- NA
  coeffs <- coefficients(fitcurves)[2:4]
  se <- sqrt(diag(vcov(fitcurves))[2:4])
  ### for record fitted.doascurve.best over all averaging intervals
  fitted.doascurve <- colSums(coeffs*t(Xreg))
  ### corresponding residual spectrum
  residual.best <- meas.doascurve[x4 + tau.best] - dyn.fixed.pattern - fitted.doascurve
  return(
    list(tau.best, delta.AICc.zero, best.order, aicctab, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
 )
}   

fit.curves.OLS.rob <- function(meas.doascurve, x4, dyn.fixed.pattern, Xreg, fit.weights, tau.shift, order2, path.length){
  require(MASS)
  fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[x4 + x] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], weights=fit.weights, model=FALSE))
  aicc <- sapply(fitcurves,AICc)
  index <- which.min(aicc)
  tau.best <- tau.shift[index]
  delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[x4] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights)) - aicc[index]

    # fitcurves <- try(rlm(as.numeric(meas.doascurve[x4 + tau.best] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights, method="MM"),silent=TRUE)
    fitcurves <- try(lm(as.numeric(meas.doascurve[x4 + tau.best] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights),silent=TRUE)
  # if(inherits(fitcurves,"try-error") || !fitcurves$converged){
  #   cat("Trying rlm 'M'\n")
  #   fitcurves <- try(rlm(as.numeric(meas.doascurve[x4 + tau.best] - dyn.fixed.pattern) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],weights=fit.weights),silent=TRUE)
  # }
  # if(inherits(fitcurves,"try-error") || !fitcurves$converged){
  #   out <- NULL
  # } else {  
    best.order <- aicctab <- NA
    coeffs <- coefficients(fitcurves)[2:4]
    se <- sqrt(diag(vcov(fitcurves))[2:4])
    ### for record fitted.doascurve.best over all averaging intervals
    fitted.doascurve <- colSums(coeffs*t(Xreg))
    ### corresponding residual spectrum
    residual.best <- meas.doascurve[x4 + tau.best] - dyn.fixed.pattern - fitted.doascurve
    out <- list(tau.best, delta.AICc.zero, best.order, aicctab, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
  # }
  return(
    out
 )
}     



forecastArima <- function (x, order = c(0, 0, 0), seasonal = c(0, 0, 0), xreg = NULL, 
    include.mean = TRUE, include.drift = FALSE, include.constant, 
    lambda = model$lambda, transform.pars = TRUE, fixed = NULL, 
    init = NULL, method = c("CSS-ML", "ML", "CSS"), n.cond, optim.control = list(maxit=1000), 
    kappa = 1e+06, model = NULL, SSinit = "Rossignol2011"){
    series <- deparse(substitute(x))
    origx <- x
    if (!is.null(lambda)) 
        x <- BoxCox(x, lambda)
    if (!is.null(xreg)) {
        nmxreg <- deparse(substitute(xreg))
        xreg <- as.matrix(xreg)
        if (ncol(xreg) == 1 & length(nmxreg) > 1) 
            nmxreg <- "xreg"
        if (is.null(colnames(xreg))) 
            colnames(xreg) <- if (ncol(xreg) == 1) 
                nmxreg
            else paste(nmxreg, 1:ncol(xreg), sep = "")
    }
    if (!is.list(seasonal)) {
        if (frequency(x) <= 1) 
            seasonal <- list(order = c(0, 0, 0), period = NA)
        else seasonal <- list(order = seasonal, period = frequency(x))
    }
    if (!missing(include.constant)) {
        if (include.constant) {
            include.mean <- TRUE
            if ((order[2] + seasonal$order[2]) == 1) 
                include.drift <- TRUE
        }
        else {
            include.mean <- include.drift <- FALSE
        }
    }
    if ((order[2] + seasonal$order[2]) > 1 & include.drift) {
        warning("No drift term fitted as the order of difference is 2 or more.")
        include.drift <- FALSE
    }
    if (!is.null(model)) {
      # cat("hac5: include model part of Arima needs to be written...\n")
        tmp <- forecast:::arima2(x, model, xreg = xreg, method = method)
        xreg <- tmp$xreg
    }
    else {
        if (include.drift) {
            drift <- 1:length(x)
            xreg <- cbind(drift = drift, xreg)
        }
        if (is.null(xreg)) 
        #     suppressWarnings(tmp <- stats::arima(x = x, order = order, 
        #         seasonal = seasonal, include.mean = include.mean, 
        #         transform.pars = transform.pars, fixed = fixed, 
        #         init = init, method = method, n.cond = n.cond, 
        #         optim.control = optim.control, kappa = kappa, SSinit = SSinit))
        # else suppressWarnings(tmp <- stats::arima(x = x, order = order, 
        #     seasonal = seasonal, xreg = xreg, include.mean = include.mean, 
        #     transform.pars = transform.pars, fixed = fixed, init = init, 
        #     method = method, n.cond = n.cond, optim.control = optim.control, 
        #     kappa = kappa, SSinit = SSinit))
            tmp <- stats::arima(x = x, order = order, 
                seasonal = seasonal, include.mean = include.mean, 
                transform.pars = transform.pars, fixed = fixed, 
                init = init, method = method, n.cond = n.cond, 
                optim.control = optim.control, kappa = kappa, SSinit = SSinit)
        else tmp <- stats::arima(x = x, order = order, 
            seasonal = seasonal, xreg = xreg, include.mean = include.mean, 
            transform.pars = transform.pars, fixed = fixed, init = init, 
            method = method, n.cond = n.cond, optim.control = optim.control, 
            kappa = kappa, SSinit = SSinit)
    }
    npar <- length(tmp$coef) + 1
    nstar <- length(tmp$residuals) - tmp$arma[6] - tmp$arma[7] * 
        tmp$arma[5]
    tmp$aicc <- tmp$aic + 2 * npar * (nstar/(nstar - npar - 1) - 
        1)
    tmp$bic <- tmp$aic + npar * (log(nstar) - 2)
    tmp$series <- series
    tmp$xreg <- xreg
    tmp$call <- match.call()
    tmp$lambda <- lambda
    tmp$x <- origx
    return(structure(tmp, class = c("ARIMA", "Arima")))
}

AICc <- function(x){
  ll <- logLik(x)
  no <- attr(ll, "nobs")
  df <- attr(ll, "df")
  (-2 * c(ll)) + (2 * df) * (1 + ((df + 1)/(no - 
  df - 1)))
}

