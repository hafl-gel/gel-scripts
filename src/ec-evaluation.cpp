
#include <Rcpp.h>
using namespace Rcpp;

// C++ helper function used in hard limits function
// [[Rcpp::export]]
Rcpp::List find_window(NumericVector x, NumericVector y1, NumericVector y2)
	{
	Rcpp::List Out(y1);
	int lenx = x.size();
    // LogicalVector LogVec = rep(false, lenx);
	int leny = y1.size();
	int run = 0;
    int run_next = 0;
	// int j = 0;
	for (int i = 0; i < leny; i++) {
        // pass j to runner & set j to 0
        run = run_next;
        // j = 0;
        LogicalVector LogVec = rep(false, lenx);
		if (x[lenx - 1] < y1[i] || x[run] > y2[i]) {
			Out[i] = LogVec;
		} else {
            // goto y1[i]
			while (run < lenx && x[run] < y1[i]) {
				run += 1;
			}
            // assign run to run_next (here should next loop start)
            // Rcout << run << "\n";
            run_next = run;
            // goto y2[i]
			while (run < lenx && x[run] < y2[i]) {
                // if (i < (leny - 1) && j == 0 && x[run] >= y1[i + 1]) {
                //     j = run;
                // }
                LogVec[run] = true;
				run += 1; 
			}
            // // check overflowing run
            // if (run == lenx) {
            //     run -= 1;
            // }
            // // check j
            // if (j == 0) {
            //     j = run;
            // }
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
            // ss += std::fabs(ogive[i] - last_value) * std::sqrt(1 / f[i]);
            ss += std::log(std::fabs(ogive[i] - last_value) + 1.0) * std::sqrt(1 / f[i]);
        }
    }
    return ss;
}

// function to fit theoretical ogive shape to measurement
// [[Rcpp::export]]
double fit_ogive_binned(const NumericVector paras, const NumericVector ogive, const NumericVector f,
    const int jlo, const int jhi, const IntegerVector jbin) {
    const double len = f.size();
    const double m = 3.0 / 4.0;
    const double fx = paras[0];
    const double mu = paras[1];
    const double A0 = paras[2];
    const double hf_offset = paras[3];
    const double A0_fx = A0 / fx;
    const double m_mu_pow = (m + 1.0) / (2.0 * mu * m);
    double last_value = 0.0;
    double ss = 0.0;
    double bin_sum = 0.0;
    double bin_avg = 0.0;
    int i = len - 1;
    int j = jbin[i];
    int sum_len = 0;
    // loop in reverse over ogive and get cumsum
    while(i >= 0 && j <= jhi) {
        // get cospec value devided by f
        last_value = last_value + 
            A0_fx / (
                std::pow(1.0 + m * std::pow(f[i] / fx, (2.0 * mu)), m_mu_pow)
            );
        // above jlo?
        if (jbin[i] >= jlo) {
            // add to bin sum
            bin_sum += last_value;
            sum_len++;
        }
        // get j
        j = jbin[i];
        // decrement
        i--;
        // check if next i new bin
        if (i < 0 || (sum_len > 0 && i >= 0 && jbin[i] != j)) {
            // get bin avg
            bin_avg = bin_sum / sum_len;
            // add to sum
            // ss += std::fabs(hf_offset + ogive[j - 1] - bin_avg);
            ss += std::log(std::fabs(hf_offset + ogive[j - 1] - bin_avg) + 1.0);
            // reset
            bin_sum = 0.0;
            sum_len = 0;
        }
    }
    return ss;
}
