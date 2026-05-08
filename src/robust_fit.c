#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

// Hardcoded LQQ parameters structure
typedef struct {
    double b, c, s;
    double rho_norm;
} LQQ_PARS;

// Forward declaration
static double fmedian_abs(double *x, int n);

// ============================================================================
// LQQ FUNCTIONS (Hardcoded for KS2014/KS2011)
// ============================================================================

static inline double lqq_rho(double x, const LQQ_PARS *k) {
    const double ax = fabs(x);
    const double b = k->b, c = k->c, s = k->s;
    const double bc = b + c;
    const double a = k->rho_norm;

    if (ax <= c) {
        return (3.0 * s - 3.0) * x * x / a;
    }
    if (ax <= bc) {
        const double d = ax - c;
        return ((3.0 * s - 3.0) * c * c + 
                (3.0 * ax + 6.0 * c - 3.0 * (ax + 2.0 * c) * d / b) * d * d) / a;
    }
    const double s5 = s - 1.0;
    const double s6 = -2.0 * bc * s5;
    const double limit = bc - s6 / s5;
    if (ax < limit) {
        const double d = ax - bc;
        const double t = s6 / s5;
        return 1.0 + s5 * s5 * d * d * (3.0 * t + 2.0 * d) / (3.0 * s6 * s6 * ax);
    }
    return 1.0;
}

static inline double lqq_psi(double x, const LQQ_PARS *k) {
    const double ax = fabs(x);
    const double b = k->b, c = k->c, s = k->s;
    const double bc = b + c;

    if (ax <= c) return x;
    if (ax <= bc) {
        const double d = ax - c;
        return copysign(ax - s * d * d / (2.0 * b), x);
    }
    const double s5 = s - 1.0;
    const double s6 = -2.0 * bc + b * s;
    const double limit = bc - s6 / s5;
    if (ax < limit) {
        const double d = ax - bc;
        return copysign(-s6 / 2.0 - s5 * s5 / s6 * (d * d / 2.0 + s6 / s5 * d), x);
    }
    return 0.0;
}

static inline double lqq_wgt(double x, const LQQ_PARS *k) {
    const double ax = fabs(x);
    const double b = k->b, c = k->c, s = k->s;

    if (ax <= c) return 1.0;
    const double bc = b + c;
    if (ax <= bc) {
        const double d = ax - c;
        return 1.0 - s * d * d / (2.0 * ax * b);
    }
    const double s5 = s - 1.0;
    const double s6 = -2.0 * bc + b * s;
    const double limit = bc - s6 / s5;
    if (ax < limit) {
        const double d = ax - bc;
        return -(s6 / 2.0 + s5 * s5 / s6 * d * (d / 2.0 + s6 / s5)) / ax;
    }
    return 0.0;
}

// ============================================================================
// OPTIMIZED LINEAR ALGEBRA (P=4)
// ============================================================================

// Solve (R'R) beta = xty where R is upper triangular (packed storage)
// R: 10 elements [R00, R01, R02, R03, R11, R12, R13, R22, R23, R33]
static void solve_chol_4(const double *R, const double *xty, double *beta) {
    // Forward substitution: solve R' z = xty
    double z[4];
    z[0] = xty[0] / R[0];
    z[1] = (xty[1] - R[1]*z[0]) / R[4];
    z[2] = (xty[2] - R[2]*z[0] - R[5]*z[1]) / R[7];
    z[3] = (xty[3] - R[3]*z[0] - R[6]*z[1] - R[8]*z[2]) / R[9];

    // Backward substitution: solve R beta = z
    beta[3] = z[3] / R[9];
    beta[2] = (z[2] - R[8]*beta[3]) / R[7];
    beta[1] = (z[1] - R[5]*beta[2] - R[6]*beta[3]) / R[4];
    beta[0] = (z[0] - R[1]*beta[1] - R[2]*beta[2] - R[3]*beta[3]) / R[0];
}

// Weighted crossproduct X'Wy for p=4 (manual unrolling)
static void crossprod_XWy(const double *X, const double *w, const double *y, 
                          int n, double *xty) {
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    for (int i = 0; i < n; i++) {
        const double wy = w[i] * y[i];
        s0 += X[i] * wy;
        s1 += X[i + n] * wy;
        s2 += X[i + 2*n] * wy;
        s3 += X[i + 3*n] * wy;
    }
    xty[0] = s0; xty[1] = s1; xty[2] = s2; xty[3] = s3;
}

// ============================================================================
// SCALE FINDING
// ============================================================================

static int find_scale_bisect(double *scale, const double *resid, int n,
                             const LQQ_PARS *k, double bb, 
                             double tol, int maxit) {
    double lo = 0.0, hi = 0.0, mid;

    for (int i = 0; i < n; i++) {
        if (fabs(resid[i]) > hi) hi = fabs(resid[i]);
    }
    if (hi == 0.0) { *scale = 1e-10; return 0; }

    // Check if solution exists in [0, hi]
    double rho_hi = 0.0;
    for (int i = 0; i < n; i++) rho_hi += lqq_rho(resid[i] / hi, k);
    rho_hi /= n;

    if (rho_hi >= bb) {
        *scale = hi;
        return 0;
    }

    for (int iter = 0; iter < maxit; iter++) {
        mid = (lo + hi) / 2.0;
        if (mid < 1e-20) { *scale = hi; return -1; }

        double rho_mid = 0.0;
        for (int i = 0; i < n; i++) rho_mid += lqq_rho(resid[i] / mid, k);
        rho_mid /= n;

        if (fabs(rho_mid - bb) < tol * bb) {
            *scale = mid;
            return 0;
        }
        if (rho_mid > bb) lo = mid; else hi = mid;
    }
    *scale = mid;
    return -1;
}

// Median absolute deviation (insertion sort for n ~ 300-400)
static double fmedian_abs(double *x, int n) {
    for (int i = 0; i < n; i++) x[i] = fabs(x[i]);

    // Insertion sort
    for (int i = 1; i < n; i++) {
        double key = x[i];
        int j = i - 1;
        while (j >= 0 && x[j] > key) {
            x[j+1] = x[j];
            j--;
        }
        x[j+1] = key;
    }
    if (n % 2) return x[n/2];
    return (x[n/2 - 1] + x[n/2]) / 2.0;
}

// ============================================================================
// IRWLS REFINEMENT
// ============================================================================

static double refine_p4(const double *X, const double *y, int n,
                        const double *R_chol,
                        double *beta,
                        const LQQ_PARS *k,
                        double *work) {
    double *w = work;
    double *xty = work + n;
    double *resid = work + n + 4;  // Temporary storage

    for (int iter = 0; iter < 50; iter++) {
        // Compute weights
        for (int i = 0; i < n; i++) {
            const double r = y[i] - (X[i]*beta[0] + X[i+n]*beta[1] + 
                                   X[i+2*n]*beta[2] + X[i+3*n]*beta[3]);
            w[i] = lqq_wgt(r, k);
        }

        crossprod_XWy(X, w, y, n, xty);

        double beta_new[4];
        solve_chol_4(R_chol, xty, beta_new);

        double max_diff = 0.0;
        for (int j = 0; j < 4; j++) {
            if (beta[j] != 0.0) {
                const double diff = fabs(beta_new[j] - beta[j]) / fabs(beta[j]);
                if (diff > max_diff) max_diff = diff;
            }
            beta[j] = beta_new[j];
        }
        if (max_diff < 1e-7) break;
    }

    // Return MADN of residuals
    for (int i = 0; i < n; i++) {
        resid[i] = y[i] - (X[i]*beta[0] + X[i+n]*beta[1] + 
                          X[i+2*n]*beta[2] + X[i+3*n]*beta[3]);
    }
    return 1.4826 * fmedian_abs(resid, n);
}

// ============================================================================
// MAIN FAST-S ALGORITHM
// ============================================================================

void fast_s_p4(const double *X, const double *y, int n,
               const LQQ_PARS *k, double bb,
               int nResample, int k_fast, int best_r,
               double tol, double scale_tol, int maxit_scale,
               const double *R_chol,
               double *coeffs, double *scale_out, double *ses,
               int *converged, double *work) {

    const int p = 4;
    const double INF = 1e20;
    double best_scales[20];
    double best_betas[80];  // 20 * 4

    for (int i = 0; i < best_r; i++) best_scales[i] = INF;

    // Stack buffer for residuals (n <= 400)
    double resid_buf[400];

    GetRNGstate();

    for (int iter = 0; iter < nResample; iter++) {
        // Sample 4 distinct indices
        int idx[4];
        idx[0] = (int)(unif_rand() * n);
        do { idx[1] = (int)(unif_rand() * n); } while (idx[1] == idx[0]);
        do { idx[2] = (int)(unif_rand() * n); } while (idx[2] == idx[0] || idx[2] == idx[1]);
        do { idx[3] = (int)(unif_rand() * n); } while (idx[3] == idx[0] || idx[3] == idx[1] || idx[3] == idx[2]);

        // Build 4x4 system (column-major for LAPACK)
        double A[16], b_vec[4];
        for (int j = 0; j < 4; j++) {
            const int row = idx[j];
            A[0*4 + j] = X[row];
            A[1*4 + j] = X[row + n];
            A[2*4 + j] = X[row + 2*n];
            A[3*4 + j] = X[row + 3*n];
            b_vec[j] = y[row];
        }

        int ipiv[4], info;
        F77_CALL(dgesv)(&p, &(int){1}, A, &p, ipiv, b_vec, &p, &info);
        if (info) continue;

        // Quick refinements
        double beta_cand[4];
        memcpy(beta_cand, b_vec, 4 * sizeof(double));

        for (int ref = 0; ref < k_fast; ref++) {
            refine_p4(X, y, n, R_chol, beta_cand, k, work);
        }

        // Compute scale
        for (int i = 0; i < n; i++) {
            resid_buf[i] = y[i] - (X[i]*beta_cand[0] + X[i+n]*beta_cand[1] + 
                                  X[i+2*n]*beta_cand[2] + X[i+3*n]*beta_cand[3]);
        }

        double sc;
        if (find_scale_bisect(&sc, resid_buf, n, k, bb, scale_tol, maxit_scale) != 0) 
            continue;

        // Insert into best list (maintain sorted)
        if (sc < best_scales[best_r - 1]) {
            int pos = best_r - 1;
            while (pos > 0 && sc < best_scales[pos - 1]) pos--;
            for (int i = best_r - 1; i > pos; i--) {
                best_scales[i] = best_scales[i-1];
                memcpy(best_betas + i*4, best_betas + (i-1)*4, 4*sizeof(double));
            }
            best_scales[pos] = sc;
            memcpy(best_betas + pos*4, beta_cand, 4*sizeof(double));
        }
    }

    PutRNGstate();

    if (best_scales[0] >= INF/2) {
        *converged = 0;
        return;
    }

    // Full refinement of best candidates
    double best_scale = INF;
    double best_beta[4];

    for (int i = 0; i < best_r; i++) {
        if (best_scales[i] >= INF/2) break;

        double beta[4];
        memcpy(beta, best_betas + i*4, 4*sizeof(double));

        // Full convergence
        for (int ref = 0; ref < 200; ref++) {
            double sc = refine_p4(X, y, n, R_chol, beta, k, work);
            (void)sc; // Could use for convergence check
            // Check convergence implicitly in refine_p4
        }

        // Final scale
        for (int j = 0; j < n; j++) {
            resid_buf[j] = y[j] - (X[j]*beta[0] + X[j+n]*beta[1] + 
                                  X[j+2*n]*beta[2] + X[j+3*n]*beta[3]);
        }
        double sc;
        find_scale_bisect(&sc, resid_buf, n, k, bb, scale_tol, maxit_scale);

        if (sc < best_scale) {
            best_scale = sc;
            memcpy(best_beta, beta, 4*sizeof(double));
        }
    }

    memcpy(coeffs, best_beta, 4*sizeof(double));
    *scale_out = best_scale;

    // Compute diagonal covariance for slopes only
    double xwx[16] = {0};
    for (int i = 0; i < n; i++) {
        const double r = y[i] - (X[i]*coeffs[0] + X[i+n]*coeffs[1] + 
                               X[i+2*n]*coeffs[2] + X[i+3*n]*coeffs[3]);
        const double wi = lqq_wgt(r / best_scale, k);
        const double xi0 = X[i], xi1 = X[i+n], xi2 = X[i+2*n], xi3 = X[i+3*n];

        xwx[0] += xi0*xi0*wi; xwx[1] += xi0*xi1*wi; xwx[2] += xi0*xi2*wi; xwx[3] += xi0*xi3*wi;
        xwx[5] += xi1*xi1*wi; xwx[6] += xi1*xi2*wi; xwx[7] += xi1*xi3*wi;
        xwx[10] += xi2*xi2*wi; xwx[11] += xi2*xi3*wi;
        xwx[15] += xi3*xi3*wi;
    }
    // Symmetrize
    xwx[4] = xwx[1]; xwx[8] = xwx[2]; xwx[12] = xwx[3];
    xwx[9] = xwx[6]; xwx[13] = xwx[7];
    xwx[14] = xwx[11];

    int info;
    char uplo = 'U';
    F77_CALL(dpotrf)(&uplo, &(int){4}, xwx, &(int){4}, &info FCONE);
    if (info == 0) {
        F77_CALL(dpotri)(&uplo, &(int){4}, xwx, &(int){4}, &info FCONE);
        // Return SEs for slopes only (indices 1,2,3)
        const double scale2 = best_scale * best_scale;
        ses[0] = sqrt(xwx[5] * scale2);   // Var(beta_1)
        ses[1] = sqrt(xwx[10] * scale2);  // Var(beta_2)
        ses[2] = sqrt(xwx[15] * scale2);  // Var(beta_3)
    } else {
        ses[0] = ses[1] = ses[2] = NA_REAL;
    }

    *converged = 1;
}

// ============================================================================
// R INTERFACE
// ============================================================================

void R_lmrob_fast_p4(double *X, double *y, int *n, double *settings,
                     int *nResample, int *k_fast, int *best_r,
                     double *tol, double *scale_tol, int *maxit_scale,
                     double *R_chol,  // 10 elements
                     double *coeffs, double *scale, double *ses,
                     int *converged, double *work) {

    LQQ_PARS k;
    k.b = settings[0];
    k.c = settings[1];
    k.s = settings[2];
    const double bc = k.b + k.c;
    k.rho_norm = (3.0 * k.s - 3.0) * k.c * k.c + 3.0 * bc * bc - 
                 2.0 * (k.s - 1.0) * bc * bc * bc / k.b;
    double bb = settings[3];

    fast_s_p4(X, y, *n, &k, bb,
              *nResample, *k_fast, *best_r,
              *tol, *scale_tol, *maxit_scale,
              R_chol, coeffs, scale, ses, converged, work);
}

