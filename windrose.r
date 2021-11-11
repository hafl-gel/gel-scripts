
hexadec <- function(x){
    hd_string <- c(0:9, letters[1:6])
    r16 <- function(r) {
        f <- (r %% 16) + 1L
        new <- floor(r / 16)
        if (new == 0) {
            hd_string[f]
        } else {
            paste0(r16(new), hd_string[f])
        }
    }
    r16(x)
}


.polar2xy <- function(wd, ws, scale = 1, center = c(0, 0), asp = 1){
    z <- scale * ws * exp((90 - wd) / 180 * pi * 1i)
    data.frame(x = Re(z) + center[1], y = Im(z) * asp + center[2])
}

windrose <- function(
    data, wd = "wd", ws = "ws", type = c("wedge", "spider"), ws_breaks = 1,
    delta_wd = 22.5, delta_freq = NULL, max_freq = NULL, circ_freq = delta_freq/4,
    scale = 1, center = c(0,0), add = FALSE, border = "black", width = 0.9,
    draw_wd = c(0, 360), exclude = FALSE, start = 0, mirror = FALSE, alpha = 0.5,
    colors = rainbow(max(10, length(breaks) - 1), end = 0.7),
    grid.angle = 45, grid.lty = 3, grid.col = "black", draw.grid = TRUE, 
    lab.col = grid.col, lab.angle = 30, key = TRUE, unit = "m/s",
    legend.cex = 0.7, legend.x = "topright", legend.y = NULL, legend.bty = "o",
    legend.text.col = par("col"), legend.bg = par("bg"), ...
    ){

    if(!requireNamespace("reshape2")) stop("install.packages('reshape2')")

    if(missing(data)){
        data <- data.frame(wd = wd, ws = ws)
        wd <- "wd"
        ws <- "ws"
    } else if(!all(c(wd, ws) %in% names(data))) stop("please provide valid column names")

    if (is.data.table(data)) data <- as.data.frame(data)

    if((360 / grid.angle) %% 1 != 0) stop("360 degrees needs to be a multiplier of 'grid.angle'")
    if(length(delta_wd) != 1) stop("argument 'delta_wd' must be of length 1")
    if((360 / delta_wd) %% 1 != 0) stop("360 degrees needs to be a multiplier of 'delta_wd'")
    type <- type[1]
    if(!(type %in% c("spider","wedge")))stop("argument 'type' must be either 'spider' or 'wedge'")
    if(width > 1 || width < 0)stop("argument width must be between 0 and 1")
    if(length(draw_wd) != 2) stop("argument 'draw_wd' must be of length 2 (exclude wd *from* *to*, clockwise)")

    df <- as.data.frame(data[,c(wd, ws)])
    df_nona <- na.exclude(df)

    ind0 <- df_nona[,ws] == 0
    if(any(ind0)) {
        cat(sprintf("%i cases (%1.1f%%) with wind speed == 0\n",
            sum(ind0),sum(ind0)/length(ind0)*100))
    }

    ## wind speed
    if (length(ws_breaks) == 1 && ws_breaks <= max(df_nona[!ind0, ws])) {
        breaks <- seq(0, max(df_nona[!ind0, ws]) + ws_breaks, by = ws_breaks)
    } else {
        # -> warning if data outside?
        breaks <- ws_breaks
    }
    wind_speed_cut <- cut(df_nona[!ind0, ws], breaks = breaks)
    # -> warning if ws == 0?
    lvls_ws <- levels(wind_speed_cut)

    ## wind direction
    # browser()
    dst <- start %% delta_wd
    wd_num <- (seq(0, 360 - delta_wd, by = delta_wd) + dst) %% 360
    wd_labels <- paste0("(",
            (wd_num - delta_wd/2) %% 360,
            ",",
            (wd_num + delta_wd/2) %% 360,
            "]")
    names(wd_labels) <- wd_num
    wind_dir_cut <- cut((df_nona[!ind0, wd] + delta_wd/2 + dst + 
        if(mirror) 180 else 0) %% 360, 
        breaks = seq(0, 360, by = delta_wd),
        labels = wd_labels
        )
    tab <- table(wind_dir=wind_dir_cut,wind_speed=wind_speed_cut)
    freq_tab <- tab/sum(tab)*100
    # margin.table(freq_tab,1)
    # margin.table(freq_tab,2)

    #### freqeuncies
    if(is.null(max_freq)) max_freq <- max(margin.table(freq_tab,1))
    if(is.null(delta_freq)){
        freqs <- pretty(c(0,max_freq))
        freqs <- freqs[freqs <= max_freq]
        delta_freq <- diff(freqs[1:2])
    } else {
        freqs <- seq(0, max_freq, delta_freq)
    }
    max_f <- max_freq + circ_freq
    
    if(!add){
        plot(c(-1.04,1.04)*max_f,c(-1.04,1.04)*max_f,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n", asp = 1, ...)
        usr <- par("usr")
        pin <- par("pin")
        asp <- pin[1]/pin[2]/diff(usr[1:2])*diff(usr[3:4])
    } else {
        usr <- par("usr")
        dx <- diff(usr[1:2])
        pin <- par("pin")
        asp <- pin[1]/pin[2]/diff(usr[1:2])*diff(usr[3:4])
        scale <- scale / (max_f * 1.04) * dx / 2.08 / asp
    }


    #### angles
    angles <- seq(0, 360 - grid.angle/2, by = grid.angle)
    plot_freqs <- freqs + circ_freq
    draw_freqs <- plot_freqs[plot_freqs!=0]

    if(draw_wd[1] != 0 || draw_wd[2] != 360){
        if(exclude){
            dwd <- -delta_wd/2
        } else {
            dwd <- delta_wd/2
        }
        if(draw_wd[1] > draw_wd[2]){
            excl <- which(wd_num + dwd < draw_wd[1] & 
                wd_num - dwd > draw_wd[2])
            fromto <- as.character(c(
                wd_num[excl[length(excl)] %% length(wd_num) + 1] - dwd, 
                wd_num[(excl[1] - 2) %% length(wd_num) + 1] + dwd) %% 360)
        } else {
            excl <- which(wd_num + dwd < draw_wd[1] | 
                wd_num - dwd > draw_wd[2])
            fromto <- as.character(c(
                min(wd_num[-excl]) - dwd, 
                max(wd_num[-excl]) + dwd) %% 360)
        }
        freq_tab[excl,] <- 0
        cat(sprintf("wind direction drawn only from %s to %s (%1.1f%% of data)\n",
            fromto[1],
            fromto[2],
            sum(freq_tab)))
    }

    ## do calc
    mfreq <- reshape2::melt(freq_tab)

    # cumulative frequencies
    cumfreq <- do.call(rbind,lapply(names(wd_labels),function(x){
        out <- mfreq[mfreq[,1] == wd_labels[x],]
        out$wd <- as.numeric(x)
        out$c_value <- cumsum(out$value) + circ_freq
        out
    }))
    # add polar coordinates
    cumfreq[,c("x","y")] <- .polar2xy(cumfreq$wd,cumfreq$c_value, scale = scale, center = center, asp = asp)

    # add alpha to colors
    if (!is.na(alpha[1])) {
        add <- unlist(lapply(alpha * 255, function(x) {
            out <- hexadec(x)
            if (nchar(out) == 1) out <- paste0('0', out)
            out
            }))
        colors <- paste0(ibts::gplots_col2hex(colors), add)
    }

    switch(type
        , "spider" = {
            # by wind_speed
            cum_ws <- by(cumfreq,cumfreq$wind_speed,identity)
            #
            xy_0 <- .polar2xy(seq(0, 360, length.out = 1E3), circ_freq, scale = scale, center = center, asp = asp)
            xy_0$y <- xy_0$y
            ord_0 <- c(seq.int(nrow(xy_0)), 1)
            ord <- c(seq_along(wd_labels), 1)
            ord_rev <- c(length(wd_labels):1,length(wd_labels))
            ## first polygon fill
            polygon(
                c(xy_0[ord_0, 1], cum_ws[[1]][ord_rev, "x"]),
                c(xy_0[ord_0, 2], cum_ws[[1]][ord_rev, "y"]),
                col = colors[1],
                border = NA
                )
            for(i in 2:length(lvls_ws)){
                ## polygon fill
                polygon(
                    c(cum_ws[[i - 1]][ord, "x"], cum_ws[[i]][ord_rev, "x"]),
                    c(cum_ws[[i - 1]][ord, "y"], cum_ws[[i]][ord_rev, "y"]),
                    col = colors[i],
                    border = NA
                    )
            }
            ## inner circle
            polygon(
                xy_0[ord_0, 1],
                xy_0[ord_0, 2],
                border = border
                )
            for(i in seq_along(lvls_ws)){
                ## polygon border
                polygon(
                    cum_ws[[i]][ord, "x"],
                    cum_ws[[i]][ord, "y"],
                    border = border
                    )
            }

        }
        , "wedge" = {
            cumfreq$wedge_l <- cumfreq$wd - delta_wd / 2 * width
            cumfreq$wedge_r <- cumfreq$wd + delta_wd / 2 * width
            # by wind_dir
            cum_wd <- by(cumfreq,cumfreq$wind_dir,function(x){
                # browser()
                xy0 <- list(.polar2xy(seq(x[1, "wedge_l"], x[1, "wedge_r"], length.out = 1E3),
                    circ_freq, scale = scale, center = center, asp = asp))
                xy <- lapply(seq_along(x$wind_speed),function(y){
                    .polar2xy(seq(x[y, "wedge_l"], x[y, "wedge_r"], length.out = 1E3),
                        x[y, "c_value"], scale = scale, center = center, asp = asp)
                })
                c(xy0,xy)
            })
            # outer loop
            ind_rev <- 1E3:1
            for(i in seq_along(cum_wd)){
                # inner loop
                for(j in seq_along(lvls_ws)){
                    # draw polygon
                    polygon(
                        c(cum_wd[[i]][[j]][,"x"], cum_wd[[i]][[j + 1]][ind_rev,"x"]),
                        c(cum_wd[[i]][[j]][,"y"], cum_wd[[i]][[j + 1]][ind_rev,"y"]),
                        col = colors[j],
                        border = border
                        )            
                }

            }
        }
        )

    if(draw.grid){
        symbols(rep(center[1],length(draw_freqs)),rep(center[2],length(draw_freqs)),circles = draw_freqs * scale, 
            lty = grid.lty, col = grid.col, add=TRUE, inches = FALSE)
        # for(ang in angles) lines(.polar2xy(ang, plot_freqs[c(1,length(plot_freqs))]), lty = grid.lty)
        line_min <- .polar2xy(angles, plot_freqs[1], scale = scale, center = center, asp = asp)
        line_max <- .polar2xy(angles, plot_freqs[length(plot_freqs)], scale = scale, center = center, asp = asp)
        segments(x0=line_min[,1],y0=line_min[,2],x1=line_max[,1],y1=line_max[,2],
            lty = grid.lty, col = grid.col)
        xy_lab <- .polar2xy(lab.angle, draw_freqs[draw_freqs != circ_freq], scale = scale, center = center, asp = asp)
        text(xy_lab, sprintf("%1.0f%%", draw_freqs[draw_freqs != circ_freq] - circ_freq), col = lab.col)
        xy_NESW <- .polar2xy(seq(0,270,90), draw_freqs[length(draw_freqs)] * 1.05, scale = scale, center = center, asp = asp)
        xy_NESESWNW <- .polar2xy(seq(45,315,90), draw_freqs[length(draw_freqs)] * 1.08, scale = scale, center = center, asp = asp)
        text(xy_NESW[1, 1], xy_NESW[1, 2], "N",col = lab.col)
        text(xy_NESW[2, 1], xy_NESW[2, 2], "E",col = lab.col)
        text(xy_NESW[3, 1], xy_NESW[3, 2], "S",col = lab.col)
        text(xy_NESW[4, 1], xy_NESW[4, 2], "W",col = lab.col)
        text(xy_NESESWNW[1, 1], xy_NESESWNW[1, 2], "NE",col = lab.col)
        text(xy_NESESWNW[2, 1], xy_NESESWNW[2, 2], "SE",col = lab.col)
        text(xy_NESESWNW[3, 1], xy_NESESWNW[3, 2], "SW",col = lab.col)
        text(xy_NESESWNW[4, 1], xy_NESESWNW[4, 2], "NW",col = lab.col)
    }

    if (key) {
        legend(x = legend.x, y= legend.y,fill = colors, legend = lvls_ws, title = unit, cex = legend.cex, bty = legend.bty,
            text.col = legend.text.col, bg = legend.bg)
    }

    invisible(freq_tab)
}
