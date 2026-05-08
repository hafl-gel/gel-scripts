

fit.curves <- function(meas.doascurve, ind_fit, Xreg, tau.shift, path.length, all_coefs = FALSE, return_resid = FALSE){
    # fit lm
    if (length(tau.shift) > 1) {
        fitcurves <- lapply(tau.shift, function(x) lm(as.numeric(meas.doascurve[ind_fit + x]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model=FALSE))
        aicc <- sapply(fitcurves,AICc)
        index <- which.min(aicc)
        tau.best <- tau.shift[index]
        fitcurves <- fitcurves[[index]]
    } else {
        fitcurves <- lm(as.numeric(meas.doascurve[ind_fit + tau.shift]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], model=FALSE)
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


prep_robust_fit <- function(Xreg) {
    # Xreg: n x 3 matrix of predictors (without intercept)
    X <- cbind(1.0, Xreg)  # Add intercept: n x 4
    n <- nrow(X)
    stopifnot(ncol(X) == 4)
    stopifnot(n < 2000)  # Hard limit for stack buffers in C
    # 1. QR decomposition for fast OLS screening of tau.shift values
    qrX <- qr(X)
    # 2. Cholesky of X'X for C code (packed upper triangular)
    XtX <- crossprod(X)  # 4x4
    R_full <- chol(XtX)  # Upper triangular 4x4
    # Pack into 10 elements: (0,0), (0,1), (0,2), (0,3), (1,1), (1,2), (1,3), (2,2), (2,3), (3,3)
    R_packed <- c(R_full[1,1], R_full[1,2], R_full[1,3], R_full[1,4],
                  R_full[2,2], R_full[2,3], R_full[2,4],
                  R_full[3,3], R_full[3,4],
                  R_full[4,4])
    # 3. Workspace for C: n (weights) + 4 (xty) + n (resid) + buffer
    work_size <- n + 4 + n + 100
    list(
        X = X,                    # Constant n x 4 design matrix
        qrX = qrX,                # For tau optimization
        R_chol = R_packed,        # 10 elements for C
        n = as.integer(n),
        work = double(work_size),
        # Settings for KS2014 (primary)
        KS2014 = list(
            settings = c(1.4734061, 0.9822707, 1.5, 0.5),  # c(b, c, s, bb)
            nResample = 1000L,
            k_fast = 2L,
            best_r = 20L,
            tol = 1e-7,
            scale_tol = 1e-10,
            maxit_scale = 200L
        ),
        # Settings for KS2011 (fallback)
        KS2011 = list(
            settings = c(0.4015457, 0.2676971, 1.5, 0.5),
            nResample = 500L,      # Can be lower for fallback
            k_fast = 1L,
            best_r = 10L,
            tol = 1e-7,
            scale_tol = 1e-10,
            maxit_scale = 200L
        )
    )
}

fit.curves.rob.fast <- function(meas.doascurve, ind_fit, X_prep, tau.shift, 
                                path.length, return_resid = FALSE) {
    # X_prep: output from prep_robust_fit()
    n <- X_prep$n
    X <- X_prep$X
    qrX <- X_prep$qrX
    # --- Step 1: Optimize tau using pre-computed QR (fast OLS) ---
    if (length(tau.shift) > 1) {
        # Vectorized RSS calculation for all tau candidates
        rss <- vapply(tau.shift, function(tau) {
            y <- as.numeric(meas.doascurve[ind_fit + tau])
            sum(qr.resid(qrX, y)^2)
        }, numeric(1))
        n_obs <- length(ind_fit)
        k <- 4  # parameters (intercept + 3 slopes)
        aicc <- n_obs * log(rss / n_obs) + 2 * k + 2 * k * (k + 1) / (n_obs - k - 1)
        tau.best <- tau.shift[which.min(aicc)]
    } else {
        tau.best <- tau.shift[1]
    }
    # Extract response
    y <- as.numeric(meas.doascurve[ind_fit + tau.best])
    # --- Step 2: Robust fit with KS2014 ---
    res <- try(.C("R_lmrob_fast_p4",
        X = as.double(X),
        y = as.double(y),
        n = X_prep$n,
        settings = as.double(X_prep$KS2014$settings),
        nResample = as.integer(X_prep$KS2014$nResample),
        k_fast = as.integer(X_prep$KS2014$k_fast),
        best_r = as.integer(X_prep$KS2014$best_r),
        tol = as.double(X_prep$KS2014$tol),
        scale_tol = as.double(X_prep$KS2014$scale_tol),
        maxit_scale = as.integer(X_prep$KS2014$maxit_scale),
        R_chol = as.double(X_prep$R_chol),
        coeffs = double(4),
        scale = double(1),
        ses = double(3),
        converged = integer(1),
        work = as.double(X_prep$work),
        PACKAGE = "gel"
    ), silent = TRUE)
    # --- Step 3: Fallback to KS2011 if needed ---
    if (inherits(res, "try-error") || !res$converged) {
        res <- .C("R_lmrob_fast_p4",
            X = as.double(X),
            y = as.double(y),
            n = X_prep$n,
            settings = as.double(X_prep$KS2011$settings),
            nResample = as.integer(X_prep$KS2011$nResample),
            k_fast = as.integer(X_prep$KS2011$k_fast),
            best_r = as.integer(X_prep$KS2011$best_r),
            tol = as.double(X_prep$KS2011$tol),
            scale_tol = as.double(X_prep$KS2011$scale_tol),
            maxit_scale = as.integer(X_prep$KS2011$maxit_scale),
            R_chol = as.double(X_prep$R_chol),
            coeffs = double(4),
            scale = double(1),
            ses = double(3),
            converged = integer(1),
            work = as.double(X_prep$work),
            PACKAGE = "gel"
        )
    }
    # --- Step 4: Format output ---
    if (inherits(res, "try-error") || !res$converged) {
        res <- try(rlm(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[, 1] + Xreg[, 2] + Xreg[, 3], method = "MM", model = FALSE), silent = TRUE)
    }
    if (inherits(res, "try-error") || !res$converged) {
        res <- try(rlm(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[, 1] + Xreg[, 2] + Xreg[, 3], model = FALSE), silent = TRUE)
    }
    # return result
    if (inherits(res, "try-error") || !res$converged) {
        out <- list(
            coef = rep(NA_real_, 3),
            se = rep(NA_real_, 3),
            tau = tau.best,
            intercept = NA_real_,
            scale = NA_real_
        )
        if (return_resid) out$resid <- rep(NA_real_, n)
        return(out)
    }
    out <- list(
        coef = res$coeffs[2:4] / path.length,    # 3 slopes only
        se = res$ses / path.length,              # Corresponding SEs
        tau = tau.best,
        intercept = res$coeffs[1],
        scale = res$scale
    )
    if (return_resid) {
        # Compute residuals if requested
        y_hat <- X %*% res$coeffs
        out$resid <- as.numeric(y - y_hat)
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

# old fit.curves.rob function
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
    fitcurves <- try(robustbase::lmrob(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], setting = "KS2014", model = FALSE), silent = TRUE)
    if(inherits(fitcurves,"try-error") || !fitcurves$converged){
        fitcurves <- try(robustbase::lmrob(as.numeric(meas.doascurve[ind_fit + tau.best]) ~ Xreg[,1] + Xreg[,2] + Xreg[,3], setting = "KS2011", model = FALSE), silent = TRUE)
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

