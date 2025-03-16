
library(Rcpp)
library(data.table)
library(ibts)

## main function to read sonic data
read_sonic <- function(FilePath) {
	bn <- basename(FilePath)
    # check if provided as qs or rds
    if (grepl('[.]qs$', bn)) {
        if (!require(qs)) {
            stop('data is provided as *.qs file -> install qs library',
                ' running "install.packages("qs")"')
        }
        qs::qread(FilePath)
    } else if (grepl('[.]rds$', bn)) {
        readRDS(FilePath)
    } else if (grepl("^(py_)?fnf_", bn)) {
        # new data format -> TODO: check HS vs Windmaster
        read_hs_ascii(FilePath)
        # read_windmaster_ascii(FilePath)
    } else {
        read_windmaster_ascii(FilePath)
    }
}

#### To Do:
# - neue readWindMaster Routine
# - Plotting & anderes Gheu aufräumen
read_windmaster_ascii <- function(FilePath){
	### get Date
	bn <- basename(FilePath)
	if(!grepl("^data_", bn)){
		# run old script
		return(read_windmaster_old_ascii(FilePath))
	}
    if (grepl('[.]gz$', bn)) {
        if (!any(grepl('R.utils', installed.packages()[, 'Package']))) {
            stop('package "R.utils" must be installed to process gz files!')
        }
    }
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	Date <- gsub("^data_.*_([0-9]{8})_.*", "\\1", bn)
	### read File
    raw <- readLines(FilePath, warn = FALSE)
    # filter out erroneous multibyte strings
    raw <- raw[grepl('^\\d{2}[0-9.:]+,[^,]+,([0-9.+-]+,){3}M,([0-9.+-]+,){2}[^,]+$', raw, 
        useBytes = TRUE)]
    out <- fread(text = raw, header = FALSE, na.strings = '999.99', showProgress = FALSE, 
        blank.lines.skip = TRUE)
    # check if file is empty
	if(nrow(out) == 0){
        cat('no valid data\n')
		return(NULL)
	}
    # check which columns to convert columns if necessary
    vnums <- paste0('V', c(3, 4, 5, 7))
    is.char <- out[, sapply(.SD, is.character), .SDcols = vnums]
    # convert to numeric
    out[, vnums[is.char] := {
        lapply(.SD, as.numeric)
    }, .SDcols = vnums[is.char]]
    # remove NA lines that come from conversion
    out <- na.omit(out)
    if (out[, .N == 0]) {
        cat('no valid data\n')
        return(NULL)
    }
    # be verbose and print sonic names:
    sonic_label <- out[, sub('^[^A-Z]*', '', unique(V2))]
    sonic_label <- unique(sonic_label[!(sonic_label == '')])
    if (length(sonic_label) == 0) stop('sonic label not available!')
    if (length(sonic_label) > 1) stop('more than one unique sonic label!')
    cat(paste0("data recorded by sonic-", tolower(sonic_label), "\n"))
    sonic_file <- sub("data_(.*)_[0-9]{8}_[0-9]{6}([.]gz)?$", "\\1", bn)
    if(sonic_label != toupper(sub("sonic-", "", sonic_file))){
        warning(paste0("Sonic label '", sonic_label, "', and hostname '", sonic_file, "' don't match!"), call. = FALSE)
    }
    # check units
    if(out[, V6[1]] != "M"){
        stop("Units of recorded data not compatible with evaluation script! Column 6 should contain 'M' for m/s!")
    }
    # fix time etc.
    out[, c(
        # remove columns
        'V1', 'V2', 'V6', 'V8', 'V9',
        # add columns
        'sonic', 'Time', 'Hz',
        # replace °C by K
        'V7'
        ) := {
        # set times correctly
        st.dec <- fast_strptime(paste(Date, V1), lt = FALSE, format = "%Y%m%d %H:%M:%OS", tz = "Etc/GMT-1")
        # get Hz (faster than as.factor(sub('[.].*', '', V1)) !
        Hz <- round(median(tabulate(trunc(as.numeric(st.dec) - as.numeric(st.dec[1])))), -1)
        if (
            hour(st.dec[.N]) == 0 && 
            (hour(st.dec[1]) != 0 || .N > (3 * Hz))
        ) {
            ## fix last 3 * Hz entries, where hour == 0
            sub.st <- st.dec[.N - seq_len(3 * Hz) + 1]
            hr <- hour(sub.st)
            st.dec[.N - seq_len(3 * Hz) + 1] <- fifelse(hr == 0, sub.st + 24 * 3600, sub.st)
        }
        # return list
        list(
            # remove columns
            NULL, NULL, NULL, NULL, NULL,
            # sonic
            sonic_label,
            # Time
            st.dec,
            # Hz
            Hz,
            # T
            V7 + 273.15
            )
    }]
    ### set Output names and order
    setnames(out, c("u", "v", "w", "T", "sonic", "Time","Hz"))
    setcolorder(out,c("Time","Hz","u","v","w","T", "sonic"))
    # check if GillBug affected sonic
    setattr(out, 'GillBug', sonic_label %in% c('C', 'D'))
    # return
    out
}

read_windmaster_old_ascii <- function(FilePath){
	### get Date
	bn <- basename(FilePath)
	Date <- gsub("^..._([0-9]{6})_.*","\\1",bn)
    ### check if rg is available
    use_rg <- try(system('rg -V', intern = TRUE), silent = TRUE)
    # use_rg <- length(system('command -v rg', intern = TRUE)) > 0
	### read File
	# browser()
    if (inherits(use_rg, 'try-error')) {
        suppressWarnings(out <- fread(cmd=paste0("grep -v -e ',,' -e '[A-Za-z]' '",path.expand(FilePath),"'"),fill = TRUE,blank.lines.skip = TRUE))
    } else {
        suppressWarnings(out <- fread(cmd=paste0("rg -v -e ',,' -e '[A-Za-z]' '",path.expand(FilePath),"'"),fill = TRUE,blank.lines.skip = TRUE))
    }
	if(nrow(out) == 0){
		cat("File empty:",path.expand(FilePath),"\n")
		return(NULL)
	}
	# remove first (empty) column
    if ('V1' %in% names(out)) out[, V1 := NULL]
    # remove NAs
 	out <- na.omit(out)
    ## get Hz
    Hz <- out[, .N, by = V2][, round(median(N), -1)]
    out[, Hz := Hz]
	## set times
	out[, st.dec := fast_strptime(paste0(Date,V2),lt = FALSE,format = "%y%m%d%H:%M:%S",tz = "Etc/GMT-1")+V3][,c("V2","V3"):=NULL]
    ## fix end of day
    if (out[, 
        hour(st.dec[.N]) == 0 && 
        (hour(st.dec[1]) != 0 || .N > (3 * Hz[1]))
        ]) {
        ## fix last 3 * Hz entries, where hour == 0
        out[.N - seq_len(3 * Hz[1]) + 1, st.dec := {
            hr <- hour(st.dec)
            fifelse(hr == 0, st.dec + 24 * 3600, st.dec)
        }]
    }
	### set Output names and order
	setnames(out,c("u", "v", "w", "T", "Hz", "Time"))
	setcolorder(out,c("Time","Hz","u","v","w","T"))
	### remove 999.99 entries
	out <- out[!(u%in%999.99|v%in%999.99|w%in%999.99|T%in%999.99),]
	### change units from °C to K
	out[, T := T + 273.15]
    setattr(out, "GillBug", TRUE)
	out
}

readSonicEVS_csv <- function(FilePath, tz = "Etc/GMT-1"){
	### read File
	suppressWarnings(out <- fread(cmd=paste0("grep -v -e ',,' -e [A-Za-z] '", 
		path.expand(FilePath), "'"), fill = TRUE, blank.lines.skip = TRUE))
	### remove rows with NA
 	out <- na.omit(out)
	### set times
	out[,st.dec := fast_strptime(paste(V1, V2), lt = FALSE, format = "%d.%m.%Y %H.%M.%OS",
		tz = tz)][,c("V1","V2"):=NULL]
	# get start time
	start_time <- out[, as.POSIXct(trunc(st.dec[1]))]
	out[, dt := trunc(as.numeric(st.dec - start_time, units = "secs"))]
	# correct new day
	out[(seq_len(.N) > 30) & dt == 0, dt := dt + 24 * 3600]
	out[dt < 0, dt := dt + 24 * 3600 - 1]
	# add st column
	out[, st := start_time + dt]
	# add Hz column
	Hz <- out[, .N, by = trunc(dt)][, round(median(N), -1)]
	out[, Hz := Hz]
	# remove columns
	out[, c("st.dec", "dt") := NULL]
	### set Output names and order
	setnames(out,c("u", "v", "w", "T", "Time","Hz"))
	setcolorder(out,c("Time","Hz","u","v","w","T"))
	### remove 999.99 entries
	out <- out[!(u%in%999.99|v%in%999.99|w%in%999.99|T%in%999.99),]
	### change units from °C to K
	out[,T := T + 273.15]
	out
}

# HS Sonic Agroscope Wauwilermoos
read_hs_ascii_old <- function(FilePath){
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
	### read File
    raw <- readLines(FilePath, warn = FALSE)
    # filter out erroneous multibyte strings
    raw <- raw[grepl('^\\d{4}-\\d{2}-\\d{2}T\\d{2}[0-9.:]+Z,([^,]+,){2}([0-9.+-]+,){4}([0-9.+-]+,){2}[^,]+$', raw, 
        useBytes = TRUE)]
    out <- fread(text = raw, header = FALSE, na.strings = '999.99', showProgress = FALSE, 
        blank.lines.skip = TRUE, tz = 'UTC')
    # check if file is empty
	if(nrow(out) == 0){
        cat('no valid data\n')
		return(NULL)
	}
    # check if first column has class POSIXct
    if (out[, inherits(class(V1), 'POSIXct')]) {
        out[, V1 := fast_strptime(V1, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE, tz = 'UTC')]
    }
    # check which columns to convert columns if necessary
    vnums <- paste0('V', 4:7)
    is.char <- out[, sapply(.SD, is.character), .SDcols = vnums]
    # convert to numeric
    out[, vnums[is.char] := {
        lapply(.SD, as.numeric)
    }, .SDcols = vnums[is.char]]
    # remove NA lines that come from conversion
    out <- na.omit(out)
    if (out[, .N == 0]) {
        cat('no valid data\n')
        return(NULL)
    }
    # be verbose and print sonic names:
    sonic_label <- 'HS'
    cat(paste0("data recorded by sonic-", tolower(sonic_label), "\n"))
    # fix time etc.
    out[, c(
        # remove columns
        'V1', 'V2', 'V3', 'V8', 'V9', 'V10',
        # add columns
        'sonic', 'Time', 'Hz'
        ) := {
        # get Hz (faster than as.factor(sub('[.].*', '', V1)) !
        Hz <- round(median(tabulate(trunc(as.numeric(V1) - as.numeric(V1[1])))), -1)
        # return list
        list(
            # remove columns
            NULL, NULL, NULL, NULL, NULL, NULL,
            # sonic
            sonic_label,
            # Time
            V1,
            # Hz
            Hz
            )
    }]
    ### set Output names and order
    setnames(out, c("u", "v", "w", "T", "sonic", "Time","Hz"))
    setcolorder(out,c("Time","Hz","u","v","w","T", "sonic"))
    # check if GillBug affected sonic
    setattr(out, 'GillBug', sonic_label %in% c('C', 'D'))
    # return
    out
}

# C++ helper function for HS data
cppFunction('
#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
Rcpp::List hs_read_cpp(String filename) {
    // open file
    std::ifstream input{filename};
    if (!input.is_open()) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\\n";
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
')

# C++ helper function for gzipped HS data
cppFunction('
#include <zlib.h>
#include <string>
#include <Rcpp.h>
Rcpp::List hs_read_cpp_gzip(String filename) {
    // open file
    gzFile input = gzopen(filename.get_cstring(), "rb");
    if (input == NULL) {
        Rcout << "Could not read file: " << filename.get_cstring() << "\\n";
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
    while (gzread(input, &c, 1) > 0) {
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
    // close properly
    if (gzclose(input) != Z_OK) {
        Rcpp::Rcout << "Failed to close file\\n";
        return R_NilValue;
    }
    return Rcpp::List::create(
		_["time_string"] = col1_time,
		_["u_string"] = col4_u,
		_["v_string"] = col5_v,
		_["w_string"] = col6_w,
		_["t_string"] = col7_T
    );
}
')

# HS Sonic Agroscope Wauwilermoos (new function using C++)
read_hs_ascii <- function(FilePath) {
    # be verbose
    cat("File:", path.expand(FilePath), "- ")
    if (grepl('[.]csv$', basename(FilePath)) {
        # uncompressed
        raw <- hs_read_cpp(normalizePath(FilePath, mustWork = FALSE))
    } else {
        # gzipped data
        raw <- hs_read_cpp_gzip(normalizePath(FilePath, mustWork = FALSE))
    }
    if (length(raw) == 0) {
        stop('File path: "', FilePath, '" is not accessible!')
    }
    out <- as.data.table(raw)
    out[, Time := fast_strptime(time_string, '%Y-%m-%dT%H:%M:%OSZ', lt = FALSE)]
    out <- na.omit(out)
    cat(paste0("data recorded by HS sonic\n"))
    out[, {
        .(Time, Hz = round(1 / median(as.numeric(diff(Time)))), u = u_string, v = v_string, 
            w = w_string, T = t_string, sonic = 'HS')
    }]
}
