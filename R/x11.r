
# scale x11
x11 <- function(ratio = NULL, frac_winbar = 1 / 12, dpi = 96, 
    monitor = 1, ...) {
    if (is.null(ratio) || !is.numeric(ratio)) {
        grDevices::x11(...)
    } else {
        if (length(ratio) == 1L) {
            width_ratio <- ratio
            height_ratio <- ratio
        } else if (length(ratio) == 2L) {
            width_ratio <- ratio[1]
            height_ratio <- ratio[2]
        } else {
            return(grDevices::x11(...))
        }
        # call xrandr
        res <- try(system("xrandr | grep ' connected'", intern = TRUE), silent = TRUE)
        if (inherits(res, 'try-error')) {
            grDevices::x11(...)
        } else {
            # Parse dimension & position, e.g. "DP-1 connected primary 3840x2160+0+0 ..."
            dims <- strsplit(gsub(".* connected( primary)? ([0-9]+x[0-9]+\\+[0-9]+\\+[0-9]+) .*", "\\2", res), "x|\\+")
            if (length(dims) > 1) {
                # get x-positions order
                xord <- order(as.integer(sapply(dims, '[', 3)))
                # pick monitor from left
                ind <- xord[monitor]
                if (is.na(ind)) {
                    # default to 0/0 position monitor (most left)
                    ind <- xord[1]
                }
                # select
                dims <- dims[ind]
            }
            dims <- unlist(dims)
            sw_px <- as.integer(dims[1])
            sh_px <- as.integer(dims[2])
            xpos <- as.integer(dims[3])
            ypos <- as.integer(dims[4])
            grDevices::x11(
                width  = sw_px / dpi * width_ratio,
                height = sh_px / dpi * height_ratio * (1 - frac_winbar),
                xpos   = xpos,
                ypos   = ypos,
                ...
            )
        }
    }
}
max11 <- function(ratio = 1, ...) {
    x11(ratio = ratio, ...)
}
# x11(ratio = 1)
# x11(ratio = 1, monitor = 2)
# x11(ratio = c(0.5, 0.7))

