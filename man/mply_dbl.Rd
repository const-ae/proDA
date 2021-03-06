% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{mply_dbl}
\alias{mply_dbl}
\alias{msply_dbl}
\title{apply function that always returns a numeric matrix}
\usage{
mply_dbl(x, FUN, ncol = 1, ...)

msply_dbl(x, FUN, ...)
}
\arguments{
\item{x}{a vector that will be passed to `vapply` or a matrix that will be
passed to apply with \code{MARGIN=1}.}

\item{FUN}{the function that returns a vector of length ncol}

\item{ncol}{the length of the vector returned by `FUN`.}

\item{...}{additional arguments to FUN}
}
\value{
a matrix of size \code{length(x) x ncol}
}
\description{
The function is modeled after `vapply`, but always returns a matrix
with one row for each iteration. You need to provide the number
of elements each function call produces beforehand (i.e. the number of
resulting columns). For a more flexible version where you don't need to
provide the number of columns see \code{\link{msply_dbl}}
}
\section{Functions}{
\itemize{
\item \code{mply_dbl}: apply function that always returns a numeric matrix

\item \code{msply_dbl}: flexible version that automatically infers the number
of columns
}}

\examples{
  # Behaves similar to sapply(), but it always returns a matrix
  t(sapply(1:5, function(i) c(i - i/3, i, i + i/3)))
  proDA:::mply_dbl(1:5, function(i) c(i - i/3, i, i + i/3), ncol=3)

  # Which can avoid some bad surprises
  t(sapply(1:5, identity))
  proDA:::mply_dbl(1:5, identity)


  # Works also with matrix input
  mat <- matrix(1:20, ncol=4)
  mat
  proDA:::msply_dbl(mat, function(i) rep(i, each=2))


}
\keyword{internal}
