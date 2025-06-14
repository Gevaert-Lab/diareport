% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_function.R
\name{add_ensembl_annotation}
\alias{add_ensembl_annotation}
\title{add_ensembl_annotation}
\usage{
add_ensembl_annotation(pe_, params)
}
\arguments{
\item{pe_}{q_feature object where to add rowdata}

\item{params}{EDF data frame}
}
\value{
status int 0 non error  / 1 error found

q_feat q-feature created / modified

error error message
}
\description{
This function adds Ensembl annotation to the assay proteinRS in the q-feature object.
Ensembl annotation file is csv file, that must include hgnc_symbol columns. Annotation columns indicated in
ensembl_col are included in the rowData(pe_[['proteinRS']]), joining the gene symbol present in DIA-NN report.
hgnc_symbol information should be the same of the Genes column present in the DIA-NN, otherwise an error is throw.
}
\author{
Andrea Argentini
}
