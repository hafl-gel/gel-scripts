#include <Rcpp.h>
using namespace Rcpp;

// Declare the C function (defined in robust_fit.c)
extern "C" void R_lmrob_fast_p4(double *X, double *y, int *n, double *settings,
                     int *nResample, int *k_fast, int *best_r,
                     double *tol, double *scale_tol, int *maxit_scale,
                     double *R_chol, double *coeffs, double *scale, 
                     double *ses, int *converged, double *work);

// [[Rcpp::export]]
List lmrob_fast_p4(NumericVector X, NumericVector y, int n, 
                   NumericVector settings,
                   int nResample, int k_fast, int best_r,
                   double tol, double scale_tol, int maxit_scale,
                   NumericVector R_chol, 
                   NumericVector work) {

    // Output containers
    NumericVector coeffs(4);
    NumericVector ses(3);
    double scale = 0.0;
    int converged = 0;

    // Call the C function
    R_lmrob_fast_p4(&X[0], &y[0], &n, &settings[0],
                    &nResample, &k_fast, &best_r,
                    &tol, &scale_tol, &maxit_scale,
                    &R_chol[0], &coeffs[0], &scale, &ses[0], 
                    &converged, &work[0]);

    return List::create(
        Named("coeffs") = coeffs,
        Named("ses") = ses,
        Named("scale") = scale,
        Named("converged") = converged == 1
    );
}

