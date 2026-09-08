
#include <iostream>
#include <fstream>
#include <string>
#include <zlib.h>
#include <Rcpp.h>
using namespace Rcpp;

// C++ helper function
// [[Rcpp::export]]
Rcpp::List licor_read_cpp(String filename, const int hertz) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    // create output (10 Hz = 870000)
    int max_lines = 87000 * hertz;
    CharacterVector col1_time(max_lines, NA_STRING);
    IntegerVector col4_DiagVal(max_lines, NA_INTEGER);
    NumericVector col6_CO2D(max_lines, NA_REAL);
    NumericVector col7_H2OD(max_lines, NA_REAL);
    NumericVector col8_Temp(max_lines, NA_REAL);
    NumericVector col9_Pres(max_lines, NA_REAL);
    NumericVector col10_Cooler(max_lines, NA_REAL);
    NumericVector col11_SFVin(max_lines, NA_REAL);
    NumericVector col12_H2OMF(max_lines, NA_REAL);
    NumericVector col13_DewPt(max_lines, NA_REAL);
    NumericVector col14_CO2SS(max_lines, NA_REAL);
    IntegerVector col15_CO2AWO(max_lines, NA_INTEGER);
    int cline = 0;
    const int max_field = 15; // we're just interested in up to this field
    const int n_fields = 17; // we need to check for complete lines
    int field = 0;
    std::vector<std::string> line(n_fields);
    // loop over lines
    char c;
    std::string s;
    bool append = true;
    while (input.get(c)) {
        if (field < n_fields) {
            // check for comma (first col)
            if (c == ',') {
                // add s to current line vector
                line[field] = s;
                // increase field counter
                field += 1;
                // stop appending
                append = false;
            } else if (c == ')') {
            // match end of field
                // add s to current line vector
                line[field] = s;
                // stop appending
                append = false;
            } else if (c == ' ') {
            // match start of field
                // increase field counter
                field += 1;
                // start appending
                append = true;
                // reset s
                s.clear();
            } else if (c == '\n') {
                // check n_fields for complete line
                if (field == (n_fields - 1)) {
                    // line, number of fields is ok
                    // assign to vectors
                    col1_time[cline] = line[0];
                    col4_DiagVal[cline] = std::stoi(line[3]);
                    col6_CO2D[cline] = std::stod(line[5]);
                    col7_H2OD[cline] = std::stod(line[6]);
                    col8_Temp[cline] = std::stod(line[7]);
                    col9_Pres[cline] = std::stod(line[8]);
                    col10_Cooler[cline] = std::stod(line[9]);
                    col11_SFVin[cline] = std::stod(line[10]);
                    col12_H2OMF[cline] = std::stod(line[11]);
                    col13_DewPt[cline] = std::stod(line[12]);
                    col14_CO2SS[cline] = std::stod(line[13]);
                    col15_CO2AWO[cline] = std::stoi(line[14]);
                }
                // else premature newline
                // reset field counter
                field = 0;
                // reset s
                s.clear();
                // start appending for first column
                append = true;
                // increase line counter
                cline += 1;
            } else if (append && field < max_field) {
                // append to string
                s += c;
            }
            // else ignore all characters until start of new field
        } else if (c == '\n') {
            // newline => line not ok (too many fields)
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // start appending for first column
            append = true;
            // increase line counter
            cline += 1;
        }
        // else ignore all characters up to newline
    }
    return Rcpp::List::create(
		_["time_string"] = col1_time,
        _["DiagVal"] = col4_DiagVal,
        _["CO2D"] = col6_CO2D,
        _["H2OD"] = col7_H2OD,
        _["Temp"] = col8_Temp,
        _["Pres"] = col9_Pres,
        _["Cooler"] = col10_Cooler,
        _["SFVin"] = col11_SFVin,
        _["H2OMF"] = col12_H2OMF,
        _["DewPt"] = col13_DewPt,
        _["CO2SS"] = col14_CO2SS,
        _["CO2AWO"] = col15_CO2AWO
    );
}

// gzip version to read raw data
// [[Rcpp::export]]
Rcpp::List licor_read_cpp_gzip(Rcpp::String filename, const int hertz) {
    // open file
    gzFile input = gzopen(filename.get_cstring(), "rb");
    if (input == NULL) {
        Rcpp::Rcout << "Could not read file: " << filename.get_cstring() << "\n";
        return R_NilValue;
    }
    // create output (10 Hz = 870000)
    int max_lines = 87000 * hertz;
    CharacterVector col1_time(max_lines, NA_STRING);
    IntegerVector col4_DiagVal(max_lines, NA_INTEGER);
    NumericVector col6_CO2D(max_lines, NA_REAL);
    NumericVector col7_H2OD(max_lines, NA_REAL);
    NumericVector col8_Temp(max_lines, NA_REAL);
    NumericVector col9_Pres(max_lines, NA_REAL);
    NumericVector col10_Cooler(max_lines, NA_REAL);
    NumericVector col11_SFVin(max_lines, NA_REAL);
    NumericVector col12_H2OMF(max_lines, NA_REAL);
    NumericVector col13_DewPt(max_lines, NA_REAL);
    NumericVector col14_CO2SS(max_lines, NA_REAL);
    IntegerVector col15_CO2AWO(max_lines, NA_INTEGER);
    int cline = 0;
    const int max_field = 15; // we're just interested in up to this field
    const int n_fields = 17; // we need to check for complete lines
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    bool append = true;
    while (gzread(input, &c, 1) > 0) {
        if (field < n_fields) {
            // check for comma (first col)
            if (c == ',') {
                // add s to current line vector
                line[field] = s;
                // increase field counter
                field += 1;
                // stop appending
                append = false;
            } else if (c == ')') {
            // match end of field
                // add s to current line vector
                line[field] = s;
                // stop appending
                append = false;
            } else if (c == ' ') {
            // match start of field
                // increase field counter
                field += 1;
                // start appending
                append = true;
                // reset s
                s.clear();
            } else if (c == '\n') {
                // check n_fields for complete line
                if (field == (n_fields - 1)) {
                    // line, number of fields is ok
                    // assign to vectors
                    col1_time[cline] = line[0];
                    col4_DiagVal[cline] = std::stoi(line[3]);
                    col6_CO2D[cline] = std::stod(line[5]);
                    col7_H2OD[cline] = std::stod(line[6]);
                    col8_Temp[cline] = std::stod(line[7]);
                    col9_Pres[cline] = std::stod(line[8]);
                    col10_Cooler[cline] = std::stod(line[9]);
                    col11_SFVin[cline] = std::stod(line[10]);
                    col12_H2OMF[cline] = std::stod(line[11]);
                    col13_DewPt[cline] = std::stod(line[12]);
                    col14_CO2SS[cline] = std::stod(line[13]);
                    col15_CO2AWO[cline] = std::stoi(line[14]);
                }
                // else premature newline
                // reset field counter
                field = 0;
                // reset s
                s.clear();
                // start appending for first column
                append = true;
                // increase line counter
                cline += 1;
            } else if (append && field < max_field) {
                // append to string
                s += c;
            }
            // else ignore all characters until start of new field
        } else if (c == '\n') {
            // newline => line not ok (too many fields)
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // start appending for first column
            append = true;
            // increase line counter
            cline += 1;
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
        Rcpp::_["DiagVal"] = col4_DiagVal,
        Rcpp::_["CO2D"] = col6_CO2D,
        Rcpp::_["H2OD"] = col7_H2OD,
        Rcpp::_["Temp"] = col8_Temp,
        Rcpp::_["Pres"] = col9_Pres,
        Rcpp::_["Cooler"] = col10_Cooler,
        Rcpp::_["SFVin"] = col11_SFVin,
        Rcpp::_["H2OMF"] = col12_H2OMF,
        Rcpp::_["DewPt"] = col13_DewPt,
        Rcpp::_["CO2SS"] = col14_CO2SS,
        Rcpp::_["CO2AWO"] = col15_CO2AWO
    );
}
