

fit.curves <- function(meas.doascurve, ind_fit, Xreg, tau.shift, path.length){
    # fit lm
    if (length(tau.shift) > 1) {
        fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[ind_fit + x]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model=FALSE))
        aicc <- sapply(fitcurves,AICc)
        index <- which.min(aicc)
        tau.best <- tau.shift[index]
        delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[ind_fit]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3])) - aicc[index]
        fitcurves <- fitcurves[[index]]
    } else {
        fitcurves <- lm(as.numeric(meas.doascurve[ind_fit]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model=FALSE)
        tau.best <- tau.shift
        delta.AICc.zero <- 0
    }
    # get coefs
    coeffs <- coefficients(fitcurves)[2:4]
    se <- sqrt(diag(vcov(fitcurves))[2:4])
    ### for record fitted.doascurve.best over all averaging intervals
    fitted.doascurve <- colSums(coeffs*t(Xreg))
    ### corresponding residual spectrum
    residual.best <- meas.doascurve[ind_fit + tau.best] - fitted.doascurve
    return(
        list(tau.best, delta.AICc.zero, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
    )
}   

fit.curves.rob <- function(meas.doascurve, ind_fit, Xreg, tau.shift, path.length){
    # tau shift?
    if (length(tau.shift) > 1) {
        fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[ind_fit + x]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],  model=FALSE))
        aicc <- sapply(fitcurves,AICc)
        index <- which.min(aicc)
        tau.best <- tau.shift[index]
        delta.AICc.zero <- AICc(lm(as.numeric(meas.doascurve[ind_fit]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3])) - aicc[index]
    } else {
        tau.best <- tau.shift
        delta.AICc.zero <- 0
    }
    # fit robust
    fitcurves <- try(lmrob(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], setting="KS2014"),silent=TRUE)
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(lmrob(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], setting="KS2011"),silent=TRUE)
    }
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(rlm(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], method="MM"),silent=TRUE)
    }
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(rlm(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3]),silent=TRUE)
    }
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        out <- NULL
    } else {  
        # get coefs
        coeffs <- coefficients(fitcurves)[2:4]
        se <- sqrt(diag(vcov(fitcurves))[2:4])
        ### for record fitted.doascurve.best over all averaging intervals
        fitted.doascurve <- colSums(coeffs*t(Xreg))
        ### corresponding residual spectrum
        residual.best <- meas.doascurve[ind_fit + tau.best] - fitted.doascurve
        out <- list(tau.best, delta.AICc.zero, coeffs/path.length, se/path.length, fitted.doascurve, residual.best)
    }
    return(
        out
    )
}     



AICc <- function(x){
    ll <- logLik(x)
    no <- attr(ll, "nobs")
    df <- attr(ll, "df")
    (-2 * c(ll)) + (2 * df) * (1 + ((df + 1)/(no - 
                df - 1)))
}

