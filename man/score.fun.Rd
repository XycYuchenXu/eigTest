% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.fun.R
\name{score.fun}
\alias{score.fun}
\title{Score function (off-diagonal squared norm or first semi-column squared norm)}
\usage{
score.fun(A, Q = NULL, k = NULL)
}
\arguments{
\item{A}{The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.}

\item{Q}{Optional. The orthogonal matrix of dimension \code{d}-\code{d} with the first \code{k} columns evaluated for common. Default (\code{is.null(Q) = TRUE}) will ignore \code{Q} for off-diagonal norm.}

\item{k}{Optional. The number of common orthogonal components to be evaluated. Only applicable when \code{is.null(Q) = FALSE}. Default (\code{is.null(k) = TRUE}) will use \code{k = 1}.}
}
\value{
Either the sum of squared off-diagonal norms of \code{Ai} or the sum of squared norms of \code{t(Q[,(k+1):d]) Ai Q[,1:k]}.
}
\description{
Calculate either the summed squared off-diagonal 'F' norms of \code{Ai} or the summed squared 'F' norms of \code{Bi[(j+1):n,j]} where \code{Bi = t(Q) Ai  Q} when \code{Q} is supplied.
}
\keyword{internal}
