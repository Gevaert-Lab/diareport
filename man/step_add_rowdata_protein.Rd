% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step_pipeline.R
\name{step_add_rowdata_protein}
\alias{step_add_rowdata_protein}
\title{step_add_rowdata_protein}
\usage{
step_add_rowdata_protein(input)
}
\arguments{
\item{input}{A list containing pe and design.}
}
\value{
The updated input list with pe modified at proteinRS level.
}
\description{
Step function that adds row data detection at the proteinRS level using add_rowdata_detection.
}
\author{
Andrea Argentini
}
