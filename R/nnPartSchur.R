#' Joint Schur decomposition for the non-negative leading component
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#'
#' @return The orthogonal matrix \code{Q} of dimension \code{d}-\code{d} with the first column to be the estimated nonnegative common component.
#' @export
#'
#' @importFrom 'pracma' quadprog
#'
#' @examples nnPartSchur(countryCoeff)
nnPartSchur = function(A){
  d = dim(A)[2]; p = dim(A)[1]
  mid = matrix(0, d, d)
  for (i in 1:p) {
    mid = mid + tcrossprod(A[i,,] - diag(d))
  }
  Q = quadprog(mid, d = rep(0, d), Aeq = rep(1, d), beq = 1, lb = rep(0, d))$xmin
  Q = cbind(Q, matrix(rnorm(d*(d-1)), nrow = d))
  Q = qr.Q(qr(Q))
  return(Q)
}
