
read_hippie <- function(file, as_ibts = TRUE, tz_data = 'UTC', tz_out = 'UTC', 
    flatten = FALSE) {
    if (length(file) > 1) {
        if (all(!grepl(',|;', file))) {
            return(lapply(file, read_hippie, as_ibts = as_ibts, tz_data = tz_data, 
                    tz_out = tz_out, flatten = flatten))
        } else {
            dat_raw <- file
        }
    } else {
        dat_raw <- readLines(file)
    }
    hdr_lines <- grep('^(Date|SDcard)', dat_raw)

    # fix any consecutive lines
    if (any(dhl <- diff(hdr_lines) == 1)) {
        dhli <-c(dhl, FALSE)
        # exclude these lines
        dat_raw <- dat_raw[-hdr_lines[dhli]]
        hdr_lines <- grep('^(Date|SDcard)', dat_raw)
    }

    # detect restarts in node-red file
    if (grepl('^SDcard_present', dat_raw[hdr_lines[1]])) {
        cycles <- suppressWarnings(
            as.integer(sub('^([^;]+;){4}(\\d).*', '\\2', dat_raw))
        )
        # get header text
        hdr_txt <- dat_raw[hdr_lines[1]]
        # get restarts
        i_restart <- which(diff(cycles) < 0) + 1
        # fix missing data -> restarts
        remain <- suppressWarnings(
            as.integer(sub('.+;', '', dat_raw))
        )
        r_check <- which(c(diff(remain) > 0, FALSE))
        r_fix <- r_check[which(remain[r_check] > 30)]
        if (length(r_fix)) {
            i_restart <- sort(unique(c(i_restart, r_fix)))
        }
        # update header lines
        hdr_lines <- c(
            hdr_lines,
            i_restart + seq_along(i_restart) - 1
        )
        # fix consecutive headers
        i_remove <- which(c(FALSE, diff(hdr_lines) == 1))
        if (length(i_remove)) {
            i_restart <- i_restart[!(i_restart %in% hdr_lines[i_remove])]
            hdr_lines <- hdr_lines[-i_remove]
        }
        # fix restarts
        for (i in rev(i_restart)) {
            # update raw data
            dat_raw <- c(
                dat_raw[1:(i - 1)],
                hdr_txt,
                dat_raw[i:length(dat_raw)]
            )
        }
    }

    if (length(hdr_lines) > 1 || hdr_lines != 1) {
        # -> instrument restarts
        # break into several results
        out <- mapply(\(from, to) {
            read_hippie(dat_raw[from:to], as_ibts = FALSE)
        }, from = hdr_lines, to = c(hdr_lines[-1] - 1, length(dat_raw)),
        SIMPLIFY = FALSE)
        # check times
        if (any(is_na <- sapply(out, \(x) x[, is.na(st[1])]))) {
            # get time interval
            if (any(!is_na)) {
                dT <- out[[which(!is_na)]][, 
                    as.numeric(et[1] - st[1], units = 'secs')]
            } else {
                dT <- 30
            }
            # fix times
            out <- lapply(seq_along(out), \(i) {
                if (is_na[i]) {
                    out[[i]][, st := et - dT]
                }
                out[[i]]
            })
        }
        # fix ibts output
        if (as_ibts) {
            out <- lapply(out, as.ibts)
        }
        if (flatten) {
            out <- do.call(rbind, out)
        }
        return(out)
    }

    # get header
    hdr <- unlist(read.table(text = dat_raw[hdr_lines[1]], nrows = 1, sep = ';'))
    # get data without header
    dat <- na.omit(fread(text = dat_raw[-hdr_lines], header = FALSE, fill = TRUE))

    setnames(dat, sub('\\s+.*$', '', hdr))

    # remove erroneous valve positions other than 1 to 8
    dat <- dat[Valve_Is %chin% as.character(1:8)]

    # TODO: 
    #   - switch HIPPIE to UTC+1
    #   - capture: 
    #       a) no data (tair, pair) 
    #       b) different values tair (-> always take minimum?)
    dat[, c('st', 'et', 'Date_Time') := {
        Time <- with_tz(fast_strptime(Date_Time, format = '%Y%m%dT%H%M%SZ', lt = FALSE, 
                tz = tz_data), tz_out)
        dT <- diff(Time)
        mdT <- median(dT, na.rm = TRUE)
        dT[is.na(dT) | dT > mdT * 1.5] <- mdT
        list(
            Time - c(mdT, dT),
            Time,
            NULL
            )
    }]
    dat[, c('Tair', 'nFlow_1', 'nFlow_2', 'Flow_1', 'Flow_2') := {
        Tair <- rowMeans(
            cbind(Temperature_outside_left, Temperature_outside_right), 
            na.rm = TRUE)
        list(
            Tair,
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
        duration_mins = as.numeric(et[.N] - st[1], units = 'mins'),
        n_data = .N,
        flow1_lmin = mean(Flow_1, na.rm = TRUE) / 1e3,
        flow2_lmin = mean(Flow_2, na.rm = TRUE) / 1e3,
        flow1_Ln_min = mean(nFlow_1, na.rm = TRUE) / 1e3,
        flow2_Ln_min = mean(nFlow_2, na.rm = TRUE) / 1e3,
        flow1_mol_min = mean(nFlow_1, na.rm = TRUE) * 101325 / 8.3144598 / 273.14 / 1e6,
        flow2_mol_min = mean(nFlow_2, na.rm = TRUE) * 101325 / 8.3144598 / 273.14 / 1e6,
        t_outside = mean(Tair) + 273.14,
        p_outside = mean(Air_pressure),
        sn = as.character(Serial[1])
        ), by = .(position = Valve_Is, cycle = Cycle_Set)
    ]

    # reorder columns
    setcolorder(dat_avg, c('st', 'et', 'sn', 'cycle', 'position'))

    # as ibts
    if (as_ibts) {
        dat <- as.ibts(dat)
        dat_avg <- as.ibts(dat_avg)
    }

    structure(dat_avg, data = dat)
}


read_vials <- function(file, sheet = 1, filter_na = TRUE, as_data_table = TRUE, return_all = FALSE) {
    if (requireNamespace('openxlsx')) {
        dat <- openxlsx::read.xlsx(file, sheet = sheet, sep.names = ' ')
    }
    # get column names
    nms <- names(dat)
    # find concentration columns
    conc_cols <- grep('NH4.*meas', nms, value = TRUE)
    # find related times
    time_cols <- grep('time.*meas', nms, value = TRUE)
    # find start
    start_col <- grep('analysis.*(start|time)', nms, value = TRUE)
    # fix start time
    dat[, start_col] <- parse_date_time3(dat[[start_col]])
    # find Cycle/Position/etc.
    cycle_cols <- grep('(C|c)ycle|(P|p)osition|(C|c)hannel|(R|r)ank', nms, value = TRUE)
    # find full/empty weights
    full_col <- grep('full', nms, value = TRUE)
    empty_col <- grep('empty', nms, value = TRUE)
    # calculate net ml
    net_col <- 'net_ml'
    dat[, net_col] <- dat[, full_col] - dat[, empty_col]
    # fix volume if weighted after analysis
    add_vol_col <- grep('add.*vol.*ml', nms, value = TRUE)
    suppressWarnings({
        add <- as.numeric(dat[, add_vol_col])
        add[is.na(add)] <- 0
        dat[, net_col] <- dat[, net_col] + add
    })
    # find vial nr.
    vial_col <- grep('vial.*nr', nms, value = TRUE)
    # find measurement date
    date_col <- grep('meas.*(D|d)ate', nms, value = TRUE)
    # find dilution and cuvette type
    dil_col <- grep('dilut(e|ion)', nms, value = TRUE)
    cuvette_col <- grep('cuvette.+type', nms, value = TRUE)
    # fix numeric values
    suppressWarnings(dat[, conc_cols] <- lapply(dat[, conc_cols], as.numeric))
    # add mg in Solution and time difference
    for (i in seq_along(conc_cols)) {
        # total NH4-N
        dat[, paste0('nh4.n.tot.mg.', i)] <- dat[, conc_cols[i]] * dat[, net_col] /
            1000 / dat[, dil_col]
        # reaction time
        dat[, paste0('reaction.time.', i)] <- difftime(parse_date_time3(dat[[time_cols[i]]]), dat[, start_col], units = 'mins')
    }
    # get max NH4-N
    dat[, 'nh4_n_tot_mg'] <- apply(dat[, paste0('nh4.n.tot.mg.', 
            seq_along(conc_cols))], 1, max, na.rm = TRUE)
    # get limits
    dat[, c('nh4_n_mg_min', 'nh4_n_mg_max')] <- list(
        c('304' = 0.015, '305' = 1.0, '303' = 2.0, '302' = 47, '502' = 100)[
            sub('^[^0-9]*([0-9]+)$', '\\1', dat[[cuvette_col]])],
        c('304' = 2.0, '305' = 12.0, '303' = 47.0, '302' = 130, '502' = 1800)[
            sub('^[^0-9]*([0-9]+)$', '\\1', dat[[cuvette_col]])]
    )
    # build output
    cols <- c(date_col, cycle_cols, net_col, 'nh4_n_tot_mg')
    if (return_all) {
        resid_cols <- nms[!(nms %in% cols)]
        out <- dat[, c(cols, resid_cols)]
    } else {
        out <- dat[, cols]
    }
    # filter NA values in cycle, position, ...
    setDT(out)
    if (return_all) {
        # fix analysis start & weights
        setnames(out, c(start_col, empty_col, full_col), 
            c('analysis.start', 'empty.weight.g', 'full.weight.g'))
    }
    # fix meas.date
    setnames(out, date_col, 'meas.date')
    # fix format
    out[, meas.date := as.Date(parse_date_time3(meas.date))]
    # fix names
    setnames(out, make.names(tolower(names(out))))
    if (filter_na) {
        out <- na.omit(out, cols = c('cycle', 'position', 'channel'))
    }
    if (!as_data_table) {
        return(as.data.frame(out))
    }
    out
}
# file <- '~/LFE/02_Daten/6-hippie/lab/Vials_blau_41-80_20240221.xlsx'
# xx <- read_vials(file)
# read_vials(file, 2)

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

