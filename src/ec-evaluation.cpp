
#include <Rcpp.h>
using namespace Rcpp;

// C++ helper function used in hard limits function
// [[Rcpp::export]]
Rcpp::List find_window(NumericVector x, NumericVector y1, NumericVector y2)
	{
	Rcpp::List Out(y1);
	int lenx = x.size();
	int leny = y1.size();
	int run = 0;
	int j = 0;
	for (int i = 0; i < leny; i++) {
        // pass j to runner & set j to 0
        run = j;
        j = 0;
        LogicalVector LogVec = rep(false, lenx);
		if ((x[lenx - 1] < y1[i]) || (x[run] > y2[i])) {
			Out[i] = LogVec;
		} else {
            // goto y1[i]
			while ((x[run] < y1[i]) && (run < lenx)) {
				run += 1;
			}
            // goto y2[i]
			while ((run < lenx) && (x[run] < y2[i])) {
                if (i < leny && j == 0 && x[run] >= y1[i + 1]) {
                    j = run;
                }
                LogVec[run] = true;
				run += 1; 
			}
            // check j
            if (j == 0) {
                j = run;
            } else if (j == lenx) {
                j -= 1;
            }
            Out[i] = LogVec;
		}
	}
	return(Out);
}

// function to fit theoretical ogive shape to measurement
// [[Rcpp::export]]
double fit_ogive(const NumericVector paras, const NumericVector ogive, const NumericVector f,
    const int ilo, const int ihi) {
    const double len = ogive.size();
    const double m = 3.0 / 4.0;
    const double fx = paras[0];
    const double mu = paras[1];
    const double A0 = paras[2];
    const double A0_fx = A0 / fx;
    const double m_mu_pow = (m + 1.0) / (2.0 * mu * m);
    double last_value = 0.0;
    double ss = 0.0;
    // loop in reverse over ogive and get cumsum
    for (int i = len - 1; i >= (ilo - 1); i--) {
        // get cospec value devided by f
        last_value = last_value + 
            A0_fx / (
                std::pow(1.0 + m * std::pow(f[i] / fx, (2.0 * mu)), m_mu_pow)
            );
        // get difference to ogive
        if (i < ihi) {
            ss += std::fabs(ogive[i] - last_value) * std::sqrt(1 / f[i]);
        }
    }
    return ss;
}

