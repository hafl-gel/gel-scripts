

fit.curves <- function(meas.doascurve, ind_fit, Xreg, tau.shift, path.length, all_coefs = FALSE, return_resid = FALSE){
    # fit lm
    if (length(tau.shift) > 1) {
        fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[ind_fit + x]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model=FALSE))
        aicc <- sapply(fitcurves,AICc)
        index <- which.min(aicc)
        tau.best <- tau.shift[index]
        fitcurves <- fitcurves[[index]]
    } else {
        fitcurves <- lm(as.numeric(meas.doascurve[ind_fit]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model=FALSE)
        tau.best <- tau.shift
    }
    # return result
    out <- as.list(c(
        # coefficients NH3/SO2/NO
        coefficients(fitcurves)[2:4] / path.length,
        # standard errors NH3/SO2/NO
        sqrt(diag(vcov(fitcurves))[2:4]) / path.length,
        tau.best,
        # intercept
        if (all_coefs) coefficients(fitcurves)[1],
        use.names = FALSE
    ))
    if (return_resid) {
        out <- c(out, list(resid(fitcurves)))
    }
    out
}   

fit.curves.rob <- function(meas.doascurve, ind_fit, Xreg, tau.shift, path.length, all_coefs = FALSE, return_resid = FALSE){
    # tau shift?
    if (length(tau.shift) > 1) {
        fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[ind_fit + x]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3],  model=FALSE))
        aicc <- sapply(fitcurves,AICc)
        index <- which.min(aicc)
        tau.best <- tau.shift[index]
    } else {
        tau.best <- tau.shift
    }
    # fit robust
    fitcurves <- try(lmrob(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], setting = "KS2014", model = FALSE), silent = TRUE)
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(lmrob(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], setting = "KS2011", model = FALSE), silent = TRUE)
    }
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(rlm(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], method = "MM", model = FALSE), silent = TRUE)
    }
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(rlm(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model = FALSE), silent = TRUE)
    }
    # return result
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        # was not able to fit - return NAs
        out <- as.list(c(rep(NA_real_, 6), NA_integer_))
        if (return_resid) {
            out <- c(out, list(rep(NA_real_, length(ind_fit))))
        }
    } else {  
        out <- as.list(c(
            # coefficients NH3/SO2/NO
            coefficients(fitcurves)[2:4] / path.length,
            # standard errors NH3/SO2/NO
            sqrt(diag(vcov(fitcurves))[2:4]) / path.length,
            tau.best,
            # intercept
            if (all_coefs) coefficients(fitcurves)[1],
            use.names = FALSE
        ))
        if (return_resid) {
            out <- c(out, list(resid(fitcurves)))
        }
    }
    out
}     

AICc <- function(x){
    ll <- logLik(x)
    no <- attr(ll, "nobs")
    df <- attr(ll, "df")
    (-2 * c(ll)) + (2 * df) * (1 + ((df + 1)/(no - 
                df - 1)))
}

