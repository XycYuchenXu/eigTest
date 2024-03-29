% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expmPartSchur.R
\name{expmPartSchur}
\alias{expmPartSchur}
\title{Optimize Common Invariant Orthogoonal Sub-spaces}
\usage{
expmPartSchur(A, k, warmup = FALSE, iter = 5000, tol = 10^(-12))
}
\arguments{
\item{A}{The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.}

\item{k}{The number of common Schur components. Must be an integer within (0, \code{d}).}

\item{warmup}{Logical, whether use \code{partSchur} for a warm-up initial value, default to \code{warmup = FALSE} \insertCite{Flury86}{eigTest}.}

\item{iter}{The maximum iteration number.}

\item{tol}{The tolerance error for iteration termination.}
}
\value{
The orthogonal matrix \code{Q} of dimension \code{d}-\code{d} with the first \code{k} columns to be the estimated common components.
}
\description{
Tangent space optimization for partially joint Schur decomposition using matrix exponential parameterizations.
See \insertCite{tensor;textual}{eigTest}.
}
\examples{
expmPartSchur(countryCoeff, 2)
}
\references{
\insertAllCited{}
}
