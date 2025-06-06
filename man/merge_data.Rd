% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ht8700-functions.r
\name{merge_data}
\alias{merge_data}
\title{Merge Sonic and HT8700 Data Based on Time}
\usage{
merge_data(basis_sonic, draw_ht = NULL, draw_licor = NULL)
}
\arguments{
\item{basis_sonic}{A data.table containing the sonic data to be used as the basis for merging.}

\item{draw_ht}{A data.table containing the HT8700 data to be merged. Default is NULL.}

\item{draw_licor}{A data.table containing the Licor data to be merged. Default is NULL.}
}
\value{
A data.table containing the merged data with the same times as `basis`.
}
\description{
This function merges sonic and HT8700 data based on time. The output contains the same times as the `basis` input. Values from `draw` will be repeated or dropped to match `basis` times. Licor data is optional and must be provided by `draw_licor`.
}
