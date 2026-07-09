
# # highpass filter function
# highpass.filter <- function(dat, DOAS.win, ...){
#     filt <- getOption('md.filter.function.list')[[DOAS.win$filter.type]](DOAS.win$filter.strength, ...)
#     dat[DOAS.win$pixel_filter] - filter(dat[DOAS.win$pixel_filter], filt, 'convolution', 2, circular = FALSE)
# }

# # complementary lowpass filter function
# lowpass.filter <- function(dat, filter.type = 'Rect', filter.strength = 5, ...){
#     filt <- getOption('md.filter.function.list')[[filter.type]](filter.strength, ...)
#     filter(dat, filt, 'convolution', 2, circular = FALSE)
# }

# new double filter
highpass.filter2 <- function(dat, filt) {
    C_cfilter <- getFromNamespace('C_cfilter', 'stats')
    dat - (
        .Call(C_cfilter, dat, filt, 2L, FALSE) +
            rev(.Call(C_cfilter, rev(dat), filt, 2L, FALSE))
        ) / 2
}

# highpass filter function
highpass.filter <- function(dat, DOAS.win, ...){
    C_cfilter <- getFromNamespace('C_cfilter', 'stats')
    filt <- getOption('md.filter.function.list')[[DOAS.win$filter.type]](DOAS.win$filter.strength, ...)
    dat[DOAS.win$pixel_filter] - (
        .Call(C_cfilter, dat[DOAS.win$pixel_filter], filt, 2L, FALSE) +
            rev(.Call(C_cfilter, rev(dat[DOAS.win$pixel_filter]), filt, 2L, FALSE))
        ) / 2
}

# complementary lowpass filter function
lowpass.filter <- function(dat, filter.type = 'Rect', filter.strength = 5, ...){
    C_cfilter <- getFromNamespace('C_cfilter', 'stats')
    filt <- getOption('md.filter.function.list')[[filter.type]](filter.strength, ...)
    (
        .Call(C_cfilter, dat, filt, 2L, FALSE) +
            rev(.Call(C_cfilter, rev(dat), filt, 2L, FALSE))
        ) / 2
}
