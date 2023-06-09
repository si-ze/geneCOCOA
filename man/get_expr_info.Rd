% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneCOCOA.R
\name{get_expr_info}
\alias{get_expr_info}
\title{Create Gene-COCOA input from expression data.}
\usage{
get_expr_info(expr, GOI, clean = TRUE)
}
\arguments{
\item{expr}{This can either be a string holding a path to a table with expression data (tsv / csv) or a data frame already loaded into R.Columns should correspond to sample names, rows to gene names.}

\item{GOI}{String specifying the gene symbol of the GOI}

\item{clean}{Bool specifying whether the expression data should be cleaned (see \code{\link{clean_df}}) or not.}
}
\value{
Returns a list of two data frames with $GOIexpr = expression of gene of interest, $expr_df = expression info on all other genes.
}
\description{
Create Gene-COCOA input from expression data.
}
\examples{
\dontrun{
# let CAD_disease: data frame holding expression data in coronary artery disease with rows=genes, columns=samples

expr_info <- get_expr_info(expr=CAD_disease, GOI="PDGFD")
}
}
