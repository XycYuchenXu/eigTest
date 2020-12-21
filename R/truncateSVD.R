#' Compute truncated SVD of a matrix.
#'
#' @param A Matrix
#' @param eps The truncation threshold of singular values for SVD
#'
#' @return A list of information about truncated SVD.
#' \itemize{
#' \item Id The truncated approximation of \code{A}.
#' \item r The rank of \code{A} after truncation.
#' \item ginv The general inverse of the truncated SVD approximation.
#' }
#' @export
#'
#' @keywords internal
#'
truncateSVD = function(A, eps){

  s = svd(A)
  A.eps = s$u %*% diag((s$d > eps)*s$d) %*% t(s$v)
  r = sum(s$d > eps)
  ginv.A = s$v %*% diag((s$d > eps)/s$d) %*% t(s$u)

  output = list(Id = A.eps, r = r, ginv = ginv.A)
  return(output)
}
