% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneCOCOA.R
\name{compute_x_RMSEs}
\alias{compute_x_RMSEs}
\title{Gets n random gene samples from a given gene set and computes regression models.}
\usage{
compute_x_RMSEs(
  GOI,
  GOI_expr,
  my_t.tpms,
  samplesize = 10,
  nsims = 1000,
  geom_mean_cors = FALSE
)
}
\arguments{
\item{geom_mean_cors}{}
}
\value{
Returns a list with $RMSEs=root mean square errors of all regression models, $geom_mean_cors=geometric mean correlation of each sample with GOI
}
\description{
Gets n random gene samples from a given gene set and computes regression models.
}
\keyword{Internal}
