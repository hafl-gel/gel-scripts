
# highpass filter function
highpass.filter <- function(dat, DOAS.win, ...){
    filt <- getOption('md.filter.function.list')[[DOAS.win$filter.type]](DOAS.win$filter.strength, ...)
    dat[DOAS.win$pixel_filter] - filter(dat[DOAS.win$pixel_filter], filt, 'convolution', 2, circular = FALSE)
}

# complementary lowpass filter function
lowpass.filter <- function(dat, filter.type = 'Rect', filter.strength = 5, ...){
    filt <- getOption('md.filter.function.list')[[filter.type]](filter.strength, ...)
    filter(dat, filt, 'convolution', 2, circular = FALSE)
}
