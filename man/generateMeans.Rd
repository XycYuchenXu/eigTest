% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateMeans.R
\name{generateMeans}
\alias{generateMeans}
\title{Generate mean matrices for simulation samples}
\usage{
generateMeans(
  d,
  p,
  k = d,
  V = matrix(0, ncol = d, nrow = d),
  v = rep(0, d),
  snr = 10,
  control.g = FALSE,
  nonneg = FALSE
)
}
\arguments{
\item{d}{Size of matrices}

\item{p}{Number of matrices}

\item{k}{Number of shared components}

\item{V}{Input of eigenvector matrix}

\item{v}{Input of stationary distribution}

\item{snr}{Signal to noice variance ratio (SNR), can be a vector}

\item{control.g}{Whether the control group of samples with perturbed common
eigenvectors should be output or not}

\item{nonneg}{Whether the generated matrix should be nonnegative as a transition
probability matrix.}
}
\value{
Array of matrices: p-by-q-by-d-by-d, where q is the number of SNRs.
        If \code{control.g = TRUE}, \code{q = length(snr) + 1} otherwise \code{q = 1}.
}
\description{
Generate mean matrices for simulation samples
}
\examples{
generateMeans(5,8,3)
}
