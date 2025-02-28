
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
# f <- '~/LFE/08_gelhub/wauwilermoos/licor/py_fnf_01_licor_2025_02_24.csv'
rm(xx, yy)
for (i in 1:10) gc()
system.time({
    xx <- readLines(f)
    yy <- gsub('[,)]*[(][^ ]+ ', ',', xx)
    bzz <- na.omit(fread(text = yy, select = columns, col.names = col_names, 
        colClasses = col_classes, fill = TRUE))
})
rm(xx, yy)
for (i in 1:10) gc()
system.time({
    xx <- fread(f, sep = ',')
    dzz <- na.omit(xx[, c(list(V1), 
        fread(text = gsub('[)]{,2}[(][^ ]+ ', ',', V2),
            select = columns[-1], col.names = col_names[-1], 
            colClasses = c('character', col_classes[-1]), fill = TRUE)
        )])
})
rm(xx, yy)
for (i in 1:10) gc()
system.time({
    czz <- read_licor_data(f)
})

as.data.table(licor_read_cpp(normalizePath(f)))

f_test <- 'test'
uu <- licor_read_cpp(normalizePath(f_test))
as.data.table(uu)

cppFunction('
#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
Rcpp::List licor_read_cpp(String filename) {
    // open file
    std::ifstream finput{filename};
    if (!finput.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\\n";
        return R_NilValue;
    }
    // create output
    int max_lines = 870000;
    std::vector<std::string> col1_time(max_lines);
    std::vector<int> col4_DiagVal(max_lines);
    std::vector<double> col6_CO2D(max_lines);
    std::vector<double> col7_H2OD(max_lines);
    std::vector<double> col8_Temp(max_lines);
    std::vector<double> col9_Pres(max_lines);
    std::vector<double> col10_Cooler(max_lines);
    std::vector<double> col11_SFVin(max_lines);
    std::vector<double> col12_H2OMF(max_lines);
    std::vector<double> col13_DewPt(max_lines);
    std::vector<double> col14_CO2SS(max_lines);
    std::vector<int> col15_CO2AWO(max_lines);
    int cline = 0;
    int n_fields = 17 - 1;
    int field = 0;
    std::vector<std::string> line(n_fields + 1);
    // loop over lines
    char c;
    bool append = true;
    // create field chars + dummy
    char s[256];
    while (finput.read(s, 14) && cline < 2) {
        // fixed read time (14 chars) -> check comma
        // check comma -> otherwise read to eol
        finput.get(c);
        if (c != \',\') {
        // erroneous line
            // check s for eol
            // if eol -> rewind => loop finput.putpack()
            // else getline until eol
            continue;
        }
        // terminate string
        s[14] = \'\\0\';
        col1_time[cline] = s;
        // fixed skip up to DiagVal -> check space
        // fixed read DiagVal -> check closing parenthesis
        // etc.
        // fixed read eol (newline or carriage return)
        // increment line
        cline++;
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
Rcpp::List licor_read_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\\n";
        return R_NilValue;
    }
    // create output
    int max_lines = 870000;
    std::vector<std::string> col1_time(max_lines);
    std::vector<int> col4_DiagVal(max_lines);
    std::vector<double> col6_CO2D(max_lines);
    std::vector<double> col7_H2OD(max_lines);
    std::vector<double> col8_Temp(max_lines);
    std::vector<double> col9_Pres(max_lines);
    std::vector<double> col10_Cooler(max_lines);
    std::vector<double> col11_SFVin(max_lines);
    std::vector<double> col12_H2OMF(max_lines);
    std::vector<double> col13_DewPt(max_lines);
    std::vector<double> col14_CO2SS(max_lines);
    std::vector<int> col15_CO2AWO(max_lines);
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
            // stop appending
            append = false;
        } else if (c == \')\') {
        // match end of field
            // add s to current line vector
            line[field] = s;
            // stop appending
            append = false;
        } else if (c == \' \') {
        // match start of field
            // increase field counter
            field += 1;
            // start appending
            append = true;
            // reset s
            s.clear();
        } else if (c == \'\\n\') {
        // check for newline -> newline
            // check field counter
            if (field == n_fields) {
                // line ok
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
                // else drop readings
            }
            // reset field counter
            field = 0;
            // reset s
            s.clear();
            // start appending for first column
            append = true;
            // increase line counter
            cline += 1;
        } else if (field > n_fields) {
            // scan to newline without consuming newline
            // this might fail
            char sp[256];
            input.get(sp, 256, \'\\n\');
        } else if (append) {
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
    # convert time
    out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    # remove time_string
    out[, time_string := NULL]
    # place Time first
    setcolorder(out, 'Time')
    # remove NA entries
    na.omit(out)
}
