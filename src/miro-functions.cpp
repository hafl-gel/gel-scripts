
/*
   file header local: t-stamp;H2O;CH4 wet;CH4 dry;N2O wet;N2O dry;FitWin3 bad fit;FitWin8 bad fit
   loggerbox adds loggerbox timestamp as first column and removes last two columns
*/

#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <Rcpp.h>
using namespace Rcpp;

// C++ helper function
// [[Rcpp::export]]
Rcpp::List miro_read_loggerbox_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    int max_lines = 870000;
    // create output
    CharacterVector col1_time(max_lines);
    CharacterVector col2_time(max_lines);
    NumericVector col3(max_lines);
    NumericVector col4(max_lines);
    NumericVector col5(max_lines);
    NumericVector col6(max_lines);
    NumericVector col7(max_lines);
    NumericVector col8(max_lines);
    NumericVector col9(max_lines);
    NumericVector col10(max_lines);
    NumericVector col11(max_lines);
    NumericVector col12(max_lines);
    int cline = 0;
    int n_fields_old = 7 - 1;
    int n_fields = 12 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    while (input.get(c)) {
        // check for comma (first field loggerbox) or semi-colon (miro)
        if (c == ',' || c == ';') {
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
                col2_time[cline] = line[1];
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
            } else if (field == n_fields_old) {
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                col2_time[cline] = line[1];
                col3[cline] = std::stod(line[2]);
                col4[cline] = std::stod(line[3]);
                col5[cline] = std::stod(line[4]);
                col6[cline] = std::stod(line[5]);
                col7[cline] = std::stod(line[6]);
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
		_["time_miro"] = col2_time,
        _["h2o"] = col3,
        _["ch4_wet"] = col4,
        _["ch4_dry"] = col5,
        _["n2o_wet"] = col6,
        _["n2o_dry"] = col7,
        _["ld_cold_plate"] = col8,
        _["outside_t"] = col9,
        _["optics_t"] = col10,
        _["p_cell"] = col11,
        _["p_valve_pwm"] = col12
    );
}

// gzip version
// [[Rcpp::export]]
Rcpp::List miro_read_loggerbox_cpp_gzip(Rcpp::String filename) {
    // open file
    gzFile input = gzopen(filename.get_cstring(), "rb");
    if (input == NULL) {
        Rcpp::Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    int max_lines = 870000;
    // create output
    CharacterVector col1_time(max_lines);
    CharacterVector col2_time(max_lines);
    NumericVector col3(max_lines);
    NumericVector col4(max_lines);
    NumericVector col5(max_lines);
    NumericVector col6(max_lines);
    NumericVector col7(max_lines);
    NumericVector col8(max_lines);
    NumericVector col9(max_lines);
    NumericVector col10(max_lines);
    NumericVector col11(max_lines);
    NumericVector col12(max_lines);
    int cline = 0;
    int n_fields_old = 7 - 1;
    int n_fields = 12 - 1;
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
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                col2_time[cline] = line[1];
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
            } else if (field == n_fields_old) {
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                col2_time[cline] = line[1];
                col3[cline] = std::stod(line[2]);
                col4[cline] = std::stod(line[3]);
                col5[cline] = std::stod(line[4]);
                col6[cline] = std::stod(line[5]);
                col7[cline] = std::stod(line[6]);
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // increase line counter
            cline += 1;
        } else if (field <= n_fields) {
            // check for comma or semi-colon
            if (c == ',' || c == ';') {
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
		_["time_string"] = col1_time,
		_["time_miro"] = col2_time,
        _["h2o"] = col3,
        _["ch4_wet"] = col4,
        _["ch4_dry"] = col5,
        _["n2o_wet"] = col6,
        _["n2o_dry"] = col7,
        _["ld_cold_plate"] = col8,
        _["outside_t"] = col9,
        _["optics_t"] = col10,
        _["p_cell"] = col11,
        _["p_valve_pwm"] = col12
    );
}

/*
    adapt these functions below to read local txt files
    filenames: "YYYY-mm-dd MGA SN??.txt"
// [[Rcpp::export]]
Rcpp::List miro_read_loggerbox_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    int max_lines = 870000;
    // create output
    CharacterVector col1_time(max_lines);
    CharacterVector col2_time(max_lines);
    IntegerVector col8(max_lines);
    IntegerVector col9(max_lines);
    NumericVector col3(max_lines);
    NumericVector col4(max_lines);
    NumericVector col5(max_lines);
    NumericVector col6(max_lines);
    NumericVector col7(max_lines);
    int cline = 0;
    int n_fields = 9 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    while (input.get(c)) {
        // check for comma (first field loggerbox) or semi-colon (miro)
        if (c == ',' || c == ';') {
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
                col2_time[cline] = line[1];
                col3[cline] = std::stod(line[2]);
                col4[cline] = std::stod(line[3]);
                col5[cline] = std::stod(line[4]);
                col6[cline] = std::stod(line[5]);
                col7[cline] = std::stod(line[6]);
                col8[cline] = std::stoi(line[7]);
                col9[cline] = std::stoi(line[8]);
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
		_["time_miro"] = col2_time,
        _["h2o"] = col3,
        _["ch4_wet"] = col4,
        _["ch4_dry"] = col5,
        _["n2o_wet"] = col6,
        _["n2o_dry"] = col7,
        _["fitwin3_bad_fit"] = col8,
        _["fitwin8_bad_fit"] = col9
    );
}

// gzip version
// [[Rcpp::export]]
Rcpp::List miro_read_loggerbox_cpp_gzip(Rcpp::String filename) {
    // open file
    gzFile input = gzopen(filename.get_cstring(), "rb");
    if (input == NULL) {
        Rcpp::Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    int max_lines = 870000;
    // create output
    CharacterVector col1_time(max_lines);
    CharacterVector col2_time(max_lines);
    IntegerVector col8(max_lines);
    IntegerVector col9(max_lines);
    NumericVector col3(max_lines);
    NumericVector col4(max_lines);
    NumericVector col5(max_lines);
    NumericVector col6(max_lines);
    NumericVector col7(max_lines);
    int cline = 0;
    int n_fields = 9 - 1;
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
                col2_time[cline] = line[1];
                col3[cline] = std::stod(line[2]);
                col4[cline] = std::stod(line[3]);
                col5[cline] = std::stod(line[4]);
                col6[cline] = std::stod(line[5]);
                col7[cline] = std::stod(line[6]);
                col8[cline] = std::stoi(line[7]);
                col9[cline] = std::stoi(line[8]);
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // increase line counter
            cline += 1;
        } else if (field <= n_fields) {
            // check for comma or semi-colon
            if (c == ',' || c == ';') {
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
		_["time_string"] = col1_time,
		_["time_miro"] = col2_time,
        _["h2o"] = col3,
        _["ch4_wet"] = col4,
        _["ch4_dry"] = col5,
        _["n2o_wet"] = col6,
        _["n2o_dry"] = col7,
        _["fitwin3_bad_fit"] = col8,
        _["fitwin8_bad_fit"] = col9
    );
}
*/
