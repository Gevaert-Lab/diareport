% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_function.R
\name{parse_comparison}
\alias{parse_comparison}
\title{parse_comparison}
\usage{
parse_comparison(input_comparison, variable_names, pe)
}
\arguments{
\item{input_comparison}{input comparison label with possible}

\item{variable_names}{list of the variable name in the formula}

\item{pe}{Q-feature object related to the analysis.}
}
\value{
list with left and right parsed value of the comparison
}
\description{
This function split the input  comparison label and extract variables and
their validity with respect to the ColData() of the current Q-feature object
}
\author{
Andrea Argentini
}
