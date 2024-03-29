% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vec.test.R
\name{vec.test}
\alias{vec.test}
\title{Calculate the p-value for normal random vectors with specified covariance matrices}
\usage{
vec.test(V.arr, cn, testType, cov.arr = NULL, eps = NULL, param.out = FALSE)
}
\arguments{
\item{V.arr}{The array of vectors, with dimension \code{p}-\code{L}, where \code{p} is the number of vectors, \code{L} is the vector length.}

\item{cn}{The convergence rate(s) to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators. If \code{length(cn) < p}, all vectors share the same rate \code{cn[1]}, otherwise \code{cn = cn[1:p]}.}

\item{testType}{The test methods, can be exact chi-squared test \code{testType = 'chi'}, and/or approximated gamma test \code{testType = 'gam'}.}

\item{cov.arr}{The array of covariance matrices corresponding to the vectors with dimension \code{p}-\code{L^2}-\code{L^2}, default will use identity matrices when \code{is.null(cov.arr) = TRUE}.}

\item{eps}{The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.}

\item{param.out}{Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output P-value.}
}
\value{
A P-value when \code{param.out=FALSE} or a list of test information when \code{param.out = TRUE}.
\itemize{
\item \code{testType}: The test methods. Either exact chi-squared test \code{'chi'}, or approximated gamma test \code{'gam'}.
\item \code{statistic}: The test statistic.
\item \code{df}: The degrees of freedom for chi-squared distribution when \code{testType = 'chi'}.
\item \code{shape}: The shape parameter in gamma distribution when \code{testType = 'gam'}.
\item \code{rate}: The rate parameter in gamma distribution when \code{testType = 'gam'}.
\item \code{pvalue}: The P-value.
}
}
\description{
Calculate the p-value for normal random vectors with specified covariance matrices
}
\keyword{internal}
