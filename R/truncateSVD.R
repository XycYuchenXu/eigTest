#' Compute truncated SVD of a matrix.
#'
#' @param A The square matrix of dimension \code{d}-\code{d} to be decomposed.
#' @param eps The truncation threshold of singular values for SVD.
#' @param tr.approx Logical, whether approximate the truncated version of \code{A}
#'
#' @return A list of information about truncated SVD.
#' \itemize{
#' \item \code{Id}: The truncated approximation of \code{A}. Included only when \code{tr.approx = TRUE}.
#' \item \code{r}: The rank of \code{A} after truncation.
#' \item \code{ginv}: The general inverse of the truncated SVD approximation.
#' }
#'
#' @keywords internal
#'
#' @import svd
#'
truncateSVD = function(A, eps, tr.approx = FALSE){
  gotit = F
  try( {s = svd(A); gotit = T}, silent = T)
  if (!gotit) {s = propack.svd(A)}
  r = sum(s$d > eps)
  d.inv = (s$d > eps)/s$d
  d.inv[is.na(d.inv)] = 0
  ginv.A = tcrossprod(tcrossprod(s$v, diag(d.inv)), s$u)

  if (tr.approx) {
    output = list(Id = tcrossprod(tcrossprod(s$u, diag((s$d > eps)*s$d)), s$v), r = r, ginv = ginv.A)
  } else {
    output = list(r = r, ginv = ginv.A)
  }
  return(output)
}
