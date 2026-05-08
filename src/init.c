#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

// Declare the C function signature
void R_lmrob_fast_p4(double *X, double *y, int *n, double *settings,
                     int *nResample, int *k_fast, int *best_r,
                     double *tol, double *scale_tol, int *maxit_scale,
                     double *R_chol, double *coeffs, double *scale, 
                     double *ses, int *converged, double *work);

// Register the callable C routines
static const R_CMethodDef CEntries[] = {
    {"R_lmrob_fast_p4", (DL_FUNC) &R_lmrob_fast_p4, 16},
    {NULL, NULL, 0}
};

// Package initialization function
// attribute_visible ensures the symbol is exported on all platforms
void attribute_visible R_init_gel(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);  // Optional but safer: forces use of registered symbols
}

