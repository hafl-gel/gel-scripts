# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

find_window <- function(x, y1, y2) {
    .Call(`_gel_find_window`, x, y1, y2)
}

fit_ogive <- function(paras, ogive, f, ilo, ihi) {
    .Call(`_gel_fit_ogive`, paras, ogive, f, ilo, ihi)
}

ht8700_read_cpp <- function(filename) {
    .Call(`_gel_ht8700_read_cpp`, filename)
}

ht8700_read_cpp_gzip <- function(filename) {
    .Call(`_gel_ht8700_read_cpp_gzip`, filename)
}

match_times <- function(time1, time2, deltat) {
    .Call(`_gel_match_times`, time1, time2, deltat)
}

decal <- function(x, y) {
    .Call(`_gel_decal`, x, y)
}

licor_read_cpp <- function(filename) {
    .Call(`_gel_licor_read_cpp`, filename)
}

licor_read_cpp_gzip <- function(filename) {
    .Call(`_gel_licor_read_cpp_gzip`, filename)
}

calc_penalty <- function(mddataIn, weightsIn, startsIn, minimaIn, reftimeIn) {
    .Call(`_gel_calc_penalty`, mddataIn, weightsIn, startsIn, minimaIn, reftimeIn)
}

cont_within_range <- function(xIn, bIn, dxIn, refLengthIn) {
    .Call(`_gel_cont_within_range`, xIn, bIn, dxIn, refLengthIn)
}

nakai_correction_2012 <- function(uIn, vIn, wIn) {
    .Call(`_gel_nakai_correction_2012`, uIn, vIn, wIn)
}

hs_read_cpp <- function(filename) {
    .Call(`_gel_hs_read_cpp`, filename)
}

hs_read_cpp_gzip <- function(filename) {
    .Call(`_gel_hs_read_cpp_gzip`, filename)
}

