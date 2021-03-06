% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{result_names}
\alias{result_names}
\title{Get the result_names}
\usage{
result_names(fit, ...)
}
\arguments{
\item{fit}{the fit to get the result_names from}

\item{...}{additional arguments used by the concrete implementation}
}
\value{
a character vector
}
\description{
Get the result_names
}
\examples{
  syn_data <- generate_synthetic_data(n_proteins = 10)
  fit <- proDA(syn_data$Y, design = syn_data$groups)
  result_names(fit)

}
