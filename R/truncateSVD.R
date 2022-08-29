#' Compute truncated SVD of a matrix.
#'
#' @param A Matrix
#' @param eps The truncation threshold of singular values for SVD
#' @param tr.approx Logical whether approximate the truncated version of \code{A}
#'
#' @return A list of information about truncated SVD.
#' \itemize{
#' \item Id The truncated approximation of \code{A}.
#' \item r The rank of \code{A} after truncation.
#' \item ginv The general inverse of the truncated SVD approximation.
#' }
#'
#' @keywords internal
#'
truncateSVD = function(A, eps, tr.approx = FALSE){

  s = svd(A)
  A.eps = s$u %*% diag((s$d > eps)*s$d) %*% t(s$v)
  r = sum(s$d > eps)
  d.inv = (s$d > eps)/s$d
  d.inv[is.na(d.inv)] = 0
  ginv.A = s$v %*% diag(d.inv) %*% t(s$u)

  if (tr.approx) {
    output = list(Id = s$u %*% diag((s$d > eps)*s$d) %*% t(s$v), r = r, ginv = ginv.A)
  } else {
    output = list(r = r, ginv = ginv.A)
  }
  return(output)
}
