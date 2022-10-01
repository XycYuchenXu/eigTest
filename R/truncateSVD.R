#' Compute truncated SVD of a (positive semi-definite) matrix.
#'
#' @param A The square matrix of dimension \code{d}-\code{d} to be decomposed.
#' @param eps The truncation threshold of singular values for SVD.
#' @param tr.approx Logical, whether approximate the truncated version of \code{A}
#'
#' @return A list of information about truncated SVD.
#' \itemize{
#' \item \code{Id}: The truncated approximation of \code{A}. Included only when \code{tr.approx = TRUE}.
#' \item \code{r}: The rank of \code{A} after truncation.
#' \item \code{rootdinv.u} The left square root matrix of \code{A}, i.e., the truncated generalized inverse \eqn{A^+ \approx}\code{tcrossprod(rootdinv.u)}.
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
#  ginv.A = tcrossprod(tcrossprod(s$v, diag(d.inv)), s$u)

  if (tr.approx) {
    output = list(Id = tcrossprod(tcrossprod(s$u, diag((s$d > eps)*s$d)), s$v), r = r, rootdinv.u = t(sqrt(d.inv) * t(s$u)))
  } else {
    output = list(r = r, rootdinv.u = t(sqrt(d.inv) * t(s$u)))
  }
  return(output)
}
