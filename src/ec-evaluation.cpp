
#include <Rcpp.h>
using namespace Rcpp;

// C++ helper function used in hard limits function
// [[Rcpp::export]]
Rcpp::List find_window(NumericVector x, NumericVector y1, NumericVector y2,
        const int maximum_size)
	{
	Rcpp::List Out(y1);
	int lenx = x.size();
	int leny = y1.size();
	int run = 0;
    int run_next = 0;
    int mody = leny / 10;
    int j;
    Rcout << "[";
	for (int i = 0; i < leny; i++) {
        if ((i % mody) < 1) Rcout << "=";
        // pass previous position to runner
        run = run_next;
		if (x[lenx - 1] < y1[i] || x[run] > y2[i]) {
			Out[i] = IntegerVector::create();
		} else {
            // goto y1[i]
			while (run < lenx && x[run] < y1[i]) {
				run += 1;
			}
            // assign run to run_next (here should next loop start)
            // Rcout << "\r00000000000000000\r" << i << "/" << run << "\n";
            run_next = run;
            IntegerVector IntVec(maximum_size);
            j = 0;
            // goto y2[i]
			while (run < lenx && x[run] < y2[i] && j < maximum_size) {
                IntVec[j] = run;
				run += 1; 
                j += 1;
			}
            Out[i] = IntVec;
		}
	}
    Rcout << "] 100%";
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

