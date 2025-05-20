
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <map>
#include <string>
#include <functional>

// even
double med_even(std::vector<double> vIn, int nIn) {
    return (vIn[nIn / 2 - 1] + vIn[nIn / 2]) / 2;
}

// uneven
double med_uneven(std::vector<double> vIn, int nIn) {
    return vIn[(nIn - 1) / 2];
}

// TODO:
// compare both?: 
//   1) sort first and get ranks (efficient for large values of slices/total)
//   2) get slice and sort afterwards (efficient for few slices in many data points)

// [[Rcpp::export]]
std::vector<double> calc_penalty(
        Rcpp::DataFrame mddataIn, Rcpp::List weightsIn,
        std::vector<int> startsIn, NumericVector minimaIn,
        int reftimeIn) {

    // initialize
    int nrow = mddataIn.nrows();
    int nstarts = startsIn.size();
    NumericVector vcol(nrow);
    NumericVector wabs = weightsIn["abs"];
    NumericVector wdiff = weightsIn["diff"];
    NumericVector wmad = weightsIn["mad"];
    double i_max = weightsIn["i_max"];
    CharacterVector colnames = mddataIn.names();
    std::vector<double> out(nstarts, 0.0);
    std::vector<double> sub(reftimeIn);
    std::vector<double> absdiff(reftimeIn);
    std::string mcall, cname;
    std::map<std::string, std::function<double(std::vector<double>, int)>> median =
        {
            {"even", med_even},
            {"uneven", med_uneven}
        };
    double med, subsum;

    // which median function to call
    if (reftimeIn % 2 == 0) {
        mcall =  "even";
    } else {
        mcall =  "uneven";
    }

    // loop over columns
    for (int col = 0; col < mddataIn.size(); col++) {

        // get column
        cname = colnames[col];
        vcol = mddataIn[cname];

        // loop over startsIn
        for (int i = 0; i < nstarts; i++) {

            // get slice (from/to inclusive)
            sub = std::vector<double>(vcol.begin() + startsIn[i] - 1, vcol.begin() + startsIn[i] + reftimeIn - 1);

            // sort vector subset
            std::sort(sub.begin(), sub.end());

            // 2) add penalty for max diff (second list entry)
            out[i] = out[i] + (sub[reftimeIn - 1] - sub[0]) * wdiff[cname];

            // get median
            med = median[mcall](sub, reftimeIn);

            // set sum-up value to zero
            subsum = 0.0;

            // loop over sorted subset
            for (int s = 0; s < reftimeIn; s++) {

                // sum up
                subsum = subsum + sub[s];

                // get absolut difference
                absdiff[s] = std::fabs(sub[s] - med);

            }

            // sort absdiff
            std::sort(absdiff.begin(), absdiff.end());

            // 3) add penalty for mad
            out[i] = out[i] + median[mcall](absdiff, reftimeIn) * wmad[cname] * 1.4826;

            // 1) add penalty for absolut level
            if (cname == "i_max") {
                // i_max
                out[i] = out[i] + (i_max - subsum / reftimeIn) * wabs[cname];
            } else {
                out[i] = out[i] + (subsum / reftimeIn - minimaIn[cname]) * wabs[cname];
            }

        }

    }

    // return penalties
    return out;
}

// [[Rcpp::export]]
// find continuous series within given range
std::vector<int> cont_within_range (
        std::vector<double> xIn, std::vector<bool> bIn, 
        const double dxIn, const int refLengthIn
        ) {

    // vars
    const int n = xIn.size();
    std::vector<int> out(n);
    int count, j;

    // loop over vector
    for (int i = 0; i < n - 1; i++) {

        // check NA
        if (bIn[i]) {

            // set count to zero
            count = 0;

            // scan forward
            j = i + 1;
            while (j < std::min(n, i + refLengthIn) && bIn[j] && std::fabs(xIn[j] - xIn[i]) <= dxIn) {
                j++;
                count++;
            }

            // assign counts
            out[i] = count;

        } else {
            out[i] = 0;
        }
    }

    // add last element
    out[n - 1] = 0;

    // return result
    return out;
}
