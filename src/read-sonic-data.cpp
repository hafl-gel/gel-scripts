
#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include <zlib.h>
using namespace Rcpp;

// C++ helper function for HS data
// [[Rcpp::export]]
Rcpp::List hs_read_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    // create output
    CharacterVector col1_time(9e5);
    NumericVector col4_u(9e5);
    NumericVector col5_v(9e5);
    NumericVector col6_w(9e5);
    NumericVector col7_T(9e5);
    int cline = 0;
    int n_fields = 10 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    while (input.get(c)) {
        // check for comma
        if (c == \',\') {
            // add s to current line vector
            line[field] = s;
            // increase field counter
            field += 1;
            // reset s
            s.clear();
        } else if (c == \'\\n\') {
        // check for newline -> newline
            // add s to current line vector
            line[field] = s;
            // check field counter
            if (field == n_fields) {
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                col4_u[cline] = std::stod(line[3]);
                col5_v[cline] = std::stod(line[4]);
                col6_w[cline] = std::stod(line[5]);
                col7_T[cline] = std::stod(line[6]);
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // increase line counter
            cline += 1;
        } else if (field < 7 && field != 1 && field != 2) {
            // append to string or new line
            s += c;
        } else if (field > n_fields) {
            // scan to newline without consuming newline
            // this might fail
            char sp[256];
            input.get(sp, 256, \'\\n\');
        }
    }
    return Rcpp::List::create(
		_["time_string"] = col1_time,
		_["u_string"] = col4_u,
		_["v_string"] = col5_v,
		_["w_string"] = col6_w,
		_["t_string"] = col7_T
    );
}

// C++ helper function for gzipped HS data
// [[Rcpp::export]]
Rcpp::List hs_read_cpp_gzip(Rcpp::String filename) {
    // open file
    gzFile input = gzopen(filename.get_cstring(), "rb");
    if (input == NULL) {
        Rcpp::Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    // create output
    Rcpp::CharacterVector col1_time(9e5);
    Rcpp::NumericVector col4_u(9e5);
    Rcpp::NumericVector col5_v(9e5);
    Rcpp::NumericVector col6_w(9e5);
    Rcpp::NumericVector col7_T(9e5);
    int cline = 0;
    int n_fields = 10 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    while (gzread(input, &c, 1) > 0) {
        if (c == \'\\n\') {
        // check for newline -> newline
            // check field counter
            if (field == n_fields) {
                // add s to current line vector
                line[field] = s;
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                col4_u[cline] = std::stod(line[3]);
                col5_v[cline] = std::stod(line[4]);
                col6_w[cline] = std::stod(line[5]);
                col7_T[cline] = std::stod(line[6]);
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // increase line counter
            cline += 1;
        } else if (field <= n_fields) {
            // check for comma
            if (c == \',\') {
                // add s to current line vector
                line[field] = s;
                // increase field counter
                field += 1;
                // reset s
                s.clear();
            } else if (field < 7 && field != 1 && field != 2) {
                // append to string or new line
                s += c;
            }
        }
        // else ignore all characters up to newline
    }
    // close properly
    if (gzclose(input) != Z_OK) {
        Rcpp::Rcout << "Failed to close file\n";
        return R_NilValue;
    }
    return Rcpp::List::create(
		Rcpp::_["time_string"] = col1_time,
		Rcpp::_["u_string"] = col4_u,
		Rcpp::_["v_string"] = col5_v,
		Rcpp::_["w_string"] = col6_w,
		Rcpp::_["t_string"] = col7_T
    );
}
