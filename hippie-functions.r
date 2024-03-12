
library(data.table)
library(ibts)



read_hippie <- function(file, as_ibts = TRUE) {
    dat_raw <- readLines(file)
    hdr_lines <- grep('^Date', dat_raw)

    if (length(hdr_lines) > 1 || hdr_lines != 1) {
        cat('check header lines\n')
        browser()
        # -> instrument restart
    }

    hdr <- unlist(read.table(text = dat_raw[1], nrows = 1, sep = ';'))
    dat <- fread(text = dat_raw[-hdr_lines], header = FALSE)

    setnames(dat, sub('\\s+.*$', '', hdr))

    # TODO: 
    #   - switch HIPPIE to UTC+1
    #   - capture: 
    #       a) no data (tair, pair) 
    #       b) different values tair (-> always take minimum?)
    dat[, c('st', 'et', 'Date_Time') := {
        Time <- fast_strptime(Date_Time, format = '%Y%m%dT%H%M%SZ', lt = FALSE, tz = 'Etc/GMT-1')
        dT <- median(diff(Time), na.rm = TRUE)
        list(
            Time - dT,
            Time,
            NULL
            )
    }]
    dat[, c('nFlow_1', 'nFlow_2', 'Flow_1', 'Flow_2') := {
        Tair <- rowMeans(
            cbind(Temperature_outside_left, Temperature_outside_right), 
            na.rm = TRUE)
        list(
            Flow_1,
            Flow_2,
            Flow_1 * 101325 / 273.15 * (Tair + 273.15) / Air_pressure,
            Flow_2 * 101325 / 273.15 * (Tair + 273.15) / Air_pressure
        )
    }]

    # get averages
    dat_avg <- dat[, .(
        st = st[1],
        et = et[.N],
        # mt = st[1] + (et[.N] - st[1]) / 2,
        duration_secs = as.numeric(et[.N] - st[1], units = 'secs'),
        N_int = .N,
        flow1_lmin = mean(Flow_1, na.rm = TRUE) / 1e3,
        flow2_lmin = mean(Flow_2, na.rm = TRUE) / 1e3
        ), by = .(Valve_Is, Cycle_Set)
    ]

    # as ibts
    if (as_ibts) {
        hi_res <- as.ibts(dat)
        lo_res <- as.ibts(dat_avg)
    }

    structure(lo_res, hi_res = hi_res)
}

# test_file <- '~/LFE/02_Daten/6-hippie/20240215.TXT'
# read_hippie(test_file)


# x11(width = 10)
# # png('flow.png', width = 480 / 7 * 10)
# plot(hi_res[, 'nFlow_1'], col = 'lightgrey', ylim = c(0, 1.2e3))
# lines(hi_res[, 'Flow_1'])
# dat_avg [, 
#     text(mt, flow1_lmin, labels = paste(Cycle_Set, Valve_Is, sep = '/'), pos = 3)
# ]
# # lines(lo_res[, 'flow1_lmin'], lwd = 2, col = 'red')
# # dev.off()

# x11(width = 10)
# # png('temp.png', width = 480 / 7 * 10)
# plot(hi_res[, 'Temperature_inside'], col = 'orange', ylim = c(20, 55), 
#     ylab = 'temperature (°C)')
# lines(hi_res[, 'Temperature_outside_left'], col = '#BBA529')
# lines(hi_res[, 'Temperature_outside_right'], col = '#BB5B29')
# legend('topleft', legend = c('inside', 'outside left', 'outside right'),
#     col = c('orange', '#BBA529', '#BB5B29'), lty = 1, bty = 'n')
# # dev.off()

# x11(width = 10)
# plot(hi_res[, 'Air_pressure'], col = 'orange')

