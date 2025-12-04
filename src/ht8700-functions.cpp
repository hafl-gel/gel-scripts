
#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <Rcpp.h>
using namespace Rcpp;

// C++ helper function
// [[Rcpp::export]]
Rcpp::List ht8700_read_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    int max_lines = 870000;
    // create output
    CharacterVector col1_time(max_lines, NA_STRING);
    CharacterVector col17(max_lines, NA_STRING);
    CharacterVector col18(max_lines, NA_STRING);
    CharacterVector col20(max_lines, NA_STRING);
    IntegerVector col2(max_lines, NA_INTEGER);
    IntegerVector col19(max_lines, NA_INTEGER);
    NumericVector col3(max_lines, NA_REAL);
    NumericVector col4(max_lines, NA_REAL);
    NumericVector col5(max_lines, NA_REAL);
    NumericVector col6(max_lines, NA_REAL);
    NumericVector col7(max_lines, NA_REAL);
    NumericVector col8(max_lines, NA_REAL);
    NumericVector col9(max_lines, NA_REAL);
    NumericVector col10(max_lines, NA_REAL);
    NumericVector col11(max_lines, NA_REAL);
    NumericVector col12(max_lines, NA_REAL);
    NumericVector col13(max_lines, NA_REAL);
    NumericVector col14(max_lines, NA_REAL);
    NumericVector col15(max_lines, NA_REAL);
    NumericVector col16(max_lines, NA_REAL);
    int cline = 0;
    int n_fields = 20 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    while (input.get(c)) {
        // check for comma
        if (c == ',') {
            // add s to current line vector
            line[field] = s;
            // increase field counter
            field += 1;
            // reset s
            s.clear();
        } else if (c == '\n') {
        // check for newline -> newline
            // add s to current line vector
            line[field] = s;
            // check field counter
            if (field == n_fields) {
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                try {
                    col2[cline] = std::stoi(line[1]);
                } catch (...) {
                    col2[cline] = NA_INTEGER;
                }
                col3[cline] = std::stod(line[2]);
                col4[cline] = std::stod(line[3]);
                col5[cline] = std::stod(line[4]);
                col6[cline] = std::stod(line[5]);
                col7[cline] = std::stod(line[6]);
                col8[cline] = std::stod(line[7]);
                col9[cline] = std::stod(line[8]);
                col10[cline] = std::stod(line[9]);
                col11[cline] = std::stod(line[10]);
                col12[cline] = std::stod(line[11]);
                col13[cline] = std::stod(line[12]);
                col14[cline] = std::stod(line[13]);
                col15[cline] = std::stod(line[14]);
                col16[cline] = std::stod(line[15]);
                col17[cline] = line[16];
                col18[cline] = line[17];
                col19[cline] = std::stoi(line[18]);
                col20[cline] = line[19];
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // increase line counter
            cline += 1;
        } else if (field <= n_fields) {
            // append to string or new line
            s += c;
        } else if (field > n_fields) {
            // scan to newline without consuming newline
            // this might fail
            char sp[256];
            input.get(sp, 256, '\n');
        }
    }
    return Rcpp::List::create(
		_["time_string"] = col1_time,
        _["sn"] = col2,
        _["nh3_ppb"] = col3,
        _["nh3_ugm3"] = col4,
        _["rh_int"] = col5,
        _["temp_int"] = col6,
        _["temp_amb"] = col7,
        _["press_amb"] = col8,
        _["oss"] = col9,
        _["peak_pos"] = col10,
        _["temp_leaser_chip"] = col11,
        _["temp_leaser_housing"] = col12,
        _["temp_mct"] = col13,
        _["temp_mct_housing"] = col14,
        _["laser_current"] = col15,
        _["ref_road_2f"] = col16,
        _["alarm_lower_bit"] = col17,
        _["alarm_upper_bit"] = col18,
        _["cleaning_flag"] = col19,
        _["notused"] = col20
    );
}

// gzip version
// [[Rcpp::export]]
Rcpp::List ht8700_read_cpp_gzip(Rcpp::String filename) {
    // open file
    gzFile input = gzopen(filename.get_cstring(), "rb");
    if (input == NULL) {
        Rcpp::Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    int max_lines = 870000;
    // create output
    CharacterVector col1_time(max_lines, NA_STRING);
    CharacterVector col17(max_lines, NA_STRING);
    CharacterVector col18(max_lines, NA_STRING);
    CharacterVector col20(max_lines, NA_STRING);
    IntegerVector col2(max_lines, NA_INTEGER);
    IntegerVector col19(max_lines, NA_INTEGER);
    NumericVector col3(max_lines, NA_REAL);
    NumericVector col4(max_lines, NA_REAL);
    NumericVector col5(max_lines, NA_REAL);
    NumericVector col6(max_lines, NA_REAL);
    NumericVector col7(max_lines, NA_REAL);
    NumericVector col8(max_lines, NA_REAL);
    NumericVector col9(max_lines, NA_REAL);
    NumericVector col10(max_lines, NA_REAL);
    NumericVector col11(max_lines, NA_REAL);
    NumericVector col12(max_lines, NA_REAL);
    NumericVector col13(max_lines, NA_REAL);
    NumericVector col14(max_lines, NA_REAL);
    NumericVector col15(max_lines, NA_REAL);
    NumericVector col16(max_lines, NA_REAL);
    int cline = 0;
    int n_fields = 20 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    while (gzread(input, &c, 1) > 0) {
        if (c == '\n') {
        // check for newline -> newline
            // check field counter
            if (field == n_fields) {
                // add s to current line vector
                line[field] = s;
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                try {
                    col2[cline] = std::stoi(line[1]);
                } catch (...) {
                    col2[cline] = NA_INTEGER;
                }
                col3[cline] = std::stod(line[2]);
                col4[cline] = std::stod(line[3]);
                col5[cline] = std::stod(line[4]);
                col6[cline] = std::stod(line[5]);
                col7[cline] = std::stod(line[6]);
                col8[cline] = std::stod(line[7]);
                col9[cline] = std::stod(line[8]);
                col10[cline] = std::stod(line[9]);
                col11[cline] = std::stod(line[10]);
                col12[cline] = std::stod(line[11]);
                col13[cline] = std::stod(line[12]);
                col14[cline] = std::stod(line[13]);
                col15[cline] = std::stod(line[14]);
                col16[cline] = std::stod(line[15]);
                col17[cline] = line[16];
                col18[cline] = line[17];
                col19[cline] = std::stoi(line[18]);
                col20[cline] = line[19];
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
            if (c == ',') {
                // add s to current line vector
                line[field] = s;
                // increase field counter
                field += 1;
                // reset s
                s.clear();
            } else {
                // append to string or new line
                s += c;
            }
        }
        // else advance until newline or eof
    }
    // close properly
    if (gzclose(input) != Z_OK) {
        Rcpp::Rcout << "Failed to close file\n";
        return R_NilValue;
    }
    return Rcpp::List::create(
		Rcpp::_["time_string"] = col1_time,
        Rcpp::_["sn"] = col2,
        Rcpp::_["nh3_ppb"] = col3,
        Rcpp::_["nh3_ugm3"] = col4,
        Rcpp::_["rh_int"] = col5,
        Rcpp::_["temp_int"] = col6,
        Rcpp::_["temp_amb"] = col7,
        Rcpp::_["press_amb"] = col8,
        Rcpp::_["oss"] = col9,
        Rcpp::_["peak_pos"] = col10,
        Rcpp::_["temp_leaser_chip"] = col11,
        Rcpp::_["temp_leaser_housing"] = col12,
        Rcpp::_["temp_mct"] = col13,
        Rcpp::_["temp_mct_housing"] = col14,
        Rcpp::_["laser_current"] = col15,
        Rcpp::_["ref_road_2f"] = col16,
        Rcpp::_["alarm_lower_bit"] = col17,
        Rcpp::_["alarm_upper_bit"] = col18,
        Rcpp::_["cleaning_flag"] = col19,
        Rcpp::_["notused"] = col20
    );
}

// merge sonic & ht8700
// [[Rcpp::export]]
List match_times(NumericVector time1, NumericVector time2, double deltat)
{
    int len1 = time1.size();
    int len2 = time2.size();
    // iterator for time2
    int j = 0;
    // time diffs
    double dj;
    double dj_1;
    // loop over time1 only
	IntegerVector index1 = seq_len(len1);
    IntegerVector index2(len1);
    // loop over time1
    for (int i = 0; i < len1; i++) {
        // advance time2 index until time2[j] > time1[i] or j == len2 - 1
        while (time2[j] < time1[i] && j < len2 - 1) {
            j++;
        }
        // check j
        dj = time2[j] - time1[i];
        if (j == 0) {
            if (dj <= deltat) {
                index2[i] = j + 1;
            }
        } else {
            // check j and j-1
            dj_1 = time1[i] - time2[j - 1];
            if (dj_1 <= dj && std::abs(dj_1) <= deltat) {
                index2[i] = j;
            } else if (std::abs(dj) <= deltat) {
                index2[i] = j + 1;
            }
        }
        // check if index2[i] is 0
        if (index2[i] == 0) {
            index1[i] = 0;
        }
    }
    return List::create(index1, index2);
}

// alarm decoding helper function
// [[Rcpp::export]]
List decal(IntegerVector x, IntegerVector y)
{
    int len1 = x.size();
    IntegerVector index2(len1);
    List out = List(len1);
    LogicalMatrix mat(len1, 32);
    // loop over x
    for (int i = 0; i < len1; i++) {
        // create empty vector
        IntegerVector l;
        if (y[i] != 0) {
            // loop
            int a = y[i];
            for (int j = 15; j >= 0; j--) {
                int p = std::pow(2, j);
                if (a >= p) {
                    // add value of j + 1 + 16
                    l.push_front(j + 17);
                    // update mat
                    mat(i, j + 16) = TRUE;
                    // update a
                    a = a % p;
                }
            }
        }
        if (x[i] == 0) {
            out[i] = 0;
        } else {
            // loop
            int a = x[i];
            for (int j = 15; j >= 0; j--) {
                int p = std::pow(2, j);
                if (a >= p) {
                    // add value of j + 1
                    l.push_front(j + 1);
                    // update mat
                    mat(i, j) = TRUE;
                    // update a
                    a = a % p;
                }
            }
        }
        // assign vector to list entry
        if (l.size() == 0) {
            l = 0;
        }
        out[i] = l;
    }
    // add matrix attribute
    out.attr("mat") = mat;
    return out;
}
