
library(data.table)
library(lubridate)
library(Rcpp)



col_names <- c('Time', 'DiagVal', 'CO2D', 'H2OD', 'Temp', 'Pres', 'Cooler', 'SFVin', 'H2OMF',
    'DewPt', 'CO2SS', 'CO2AWO')
columns <- c(1, 3, 5:14)

col_classes <- c('POSIXct', 'character', 'numeric', 'character', 
    rep('numeric', 10), 'character', 'character')

f <- '~/LFE/08_gelhub/wauwilermoos/licor/py_fnf_01_licor_2025_02_19.csv'
# f <- '~/LFE/08_gelhub/wauwilermoos/licor/py_fnf_01_licor_2025_02_20.csv'
for (i in 1:10) gc()
system.time({
    xx <- readLines(f)
    yy <- gsub('[,)]*[(][^ ]+ ', ',', xx)
    bzz <- na.omit(fread(text = yy, select = columns, col.names = col_names, 
        colClasses = col_classes, fill = TRUE))
})
for (i in 1:10) gc()
system.time({
    czz <- read_licor_data(f)
})



f <- '~/LFE/08_gelhub/wauwilermoos/licor/py_fnf_01_licor_2025_02_19.csv'
xx <- readLines(f, n = 1)


cppFunction('
#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
Rcpp::List licor_read_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\\n";
        return R_NilValue;
    }
    // create output
    int max_lines = 9e5;
    CharacterVector col1_time(max_lines);
    CharacterVector col2_DiagVal(max_lines);
    CharacterVector col3_CO2D(max_lines);
    CharacterVector col4_H2OD(max_lines);
    CharacterVector col5_Temp(max_lines);
    CharacterVector col6_Pres(max_lines);
    CharacterVector col7_Cooler(max_lines);
    CharacterVector col8_SFVin(max_lines);
    CharacterVector col9_H2OMF(max_lines);
    CharacterVector col10_DewPt(max_lines);
    CharacterVector col11_CO2SS(max_lines);
    CharacterVector col12_CO2AWO(max_lines);
    int cline = 0;
    int n_fields = 17 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    bool append = true;
    std::vector<int> skip_field_chars = {0, 14, hier weiter !!! };
    while (input.get(c)) {
        if (field > 0) {
        // all fields with opening parenthesis as separator
            // state: c == "("
            // skip predefined characters until space
            for (int i = 0; i < skip_field_chars[field]; i++) {
                // read next char (skip initial opening parenthesis)
                input.get(c);
                // check for newline
                if (c == "\\n") {
                    // reset field
                    field = 0;
                    // reset s
                    s.clear();
                    break
                }
            }
            // process next line if newline has been found
            if (field == 0) continue
            // state: c == " "
            // skip space
            input.get(c);
            // state: c == first char of field value
            // run until we find the closing parenthesis
            while (c != ")") {
                // check for newline
                if (c == "\\n") {
                    // reset field
                    field = 0;
                    // reset s
                    s.clear();
                    break
                }
                // append to s
                s += c;
                // read next char
                input.get(c);
            }
            // process next line if newline has been found
            if (field == 0) continue
            // state: c == ")"
            // ~~~ > field ok here
            // assign to line
            line[field] = s;
            // reset s
            s.clear();
            // skip to next char (either "(" or "\n")
            input.get(c)
            // state: c == "(" or c == "\n"
            // check for newline
            if (c == "\\n") {
                if (field == n_fields) {
                    // assign to vectors
                    col1_time[cline] = line[0];
                    col2_DiagVal[cline] = line[1];
                    col3_CO2D[cline] = line[2];
                    col4_H2OD[cline] = line[3];
                    col5_Temp[cline] = line[4];
                    col6_Pres[cline] = line[5];
                    col7_Cooler[cline] = line[6];
                    col8_SFVin[cline] = line[7];
                    col9_H2OMF[cline] = line[8];
                    col10_DewPt[cline] = line[9];
                    col11_CO2SS[cline] = line[10];
                    col12_CO2AWO[cline] = line[11];
                }
                // reset field
                field = 0;
            } else {
                // increase field
                field++;
            }
        } else {
        // first field (time) (least frequent -> place last)
            // increase field (to use in tests for newline)
            field++;
            // run for 24 characters (time string)
            for (int i = 0; i < 24; i++) {
                // check for newline
                if (c == "\\n") {
                    // reset field
                    field = 0;
                    // reset s
                    s.clear();
                    break
                }
                // otherwise append to s
                s += c;
                // read next character
                input.get(c);
            }
            // process next line if newline has been found
            if (field == 0) continue
            // skip comma + next 6 characters until opening parenthesis
            for (int i = 0; i < 7; i++) {
                // check for newline
                if (c == "\\n") {
                    // reset field
                    field = 0;
                    // reset s
                    s.clear();
                    break
                }
                // read next char
                input.get(c);
            }
            // process next line if newline has been found
            if (field == 0) continue
            // first field ok here
            // assign to line
            line[0] = s;
        }
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
')

cppFunction('
#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
Rcpp::List licor_read_cpp2(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\\n";
        return R_NilValue;
    }
    // create output
    CharacterVector col1_time(9e5);
    CharacterVector col4_DiagVal(9e5);
    CharacterVector col6_CO2D(9e5);
    CharacterVector col7_H2OD(9e5);
    CharacterVector col8_Temp(9e5);
    CharacterVector col9_Pres(9e5);
    CharacterVector col10_Cooler(9e5);
    CharacterVector col11_SFVin(9e5);
    CharacterVector col12_H2OMF(9e5);
    CharacterVector col13_DewPt(9e5);
    CharacterVector col14_CO2SS(9e5);
    CharacterVector col15_CO2AWO(9e5);
    int cline = 0;
    int n_fields = 17 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    std::string s;
    bool append = true;
    while (input.get(c)) {
        // check for comma (first col)
        if (c == \',\') {
            // add s to current line vector
            line[field] = s;
            // increase field counter
            field += 1;
            // reset s
            s.clear();
            // stop appending
            append = false;
        } else if (c == \')\') {
        // stop appending
            append = false;
        } else if (c == \' \') {
        // check for space (all other columns)
            // add s to current line vector
            line[field] = s;
            // increase field counter
            field += 1;
            // reset s
            s.clear();
            // start appending
            append = true;
        } else if (c == \'\\n\') {
        // check for newline -> newline
            // add s to current line vector
            line[field] = s;
            // check field counter
            if (field == n_fields) {
                // line ok
                // assign to vectors
                col1_time[cline] = line[0];
                col4_DiagVal[cline] = line[3];
                col6_CO2D[cline] = line[5];
                col7_H2OD[cline] = line[6];
                col8_Temp[cline] = line[7];
                col9_Pres[cline] = line[8];
                col10_Cooler[cline] = line[9];
                col11_SFVin[cline] = line[10];
                col12_H2OMF[cline] = line[11];
                col13_DewPt[cline] = line[12];
                col14_CO2SS[cline] = line[13];
                col15_CO2AWO[cline] = line[14];
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // start appending for first column
            append = true;
            // reset c for appending
            c = \'\\0\';
            // increase line counter
            cline += 1;
        } else if (field > n_fields) {
            // check field counter > number of columns -> newline
            // reset field counter
            field = 0;
            // drop readings
            // reset s
            s.clear();
            // start appending for first column
            append = true;
            // reset c for appending
            c = \'\\0\';
            // increase line counter
            cline += 1;
        }
        if (append && c) {
            // append to string or new line
            s += c;
        }
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
')

read_licor_data <- function(FilePath) {
    raw_list <- licor_read_cpp(normalizePath(FilePath))
    out <- as.data.table(raw_list)
    # convert to integer
    int_cols <- c('DiagVal', 'CO2AWO')
    suppressWarnings(out[, (int_cols) := lapply(.SD, as.integer), .SDcols = int_cols])
    # convert to numeric
    num_cols <- names(raw_list)[-c(1, 2, 12)]
    suppressWarnings(out[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols])
    # convert time
    out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    # remove time_string
    out[, time_string := NULL]
    # place Time first
    setcolorder(out, 'Time')
    # remove NA entries
    na.omit(out)
}

