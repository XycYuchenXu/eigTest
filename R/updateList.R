#' Update the array of matrices
#'
#' Calculate the array of matrices after orthogonal transformation by matrix \code{Q}, i.e., \code{t(Q) Ai Q}.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param Q The orthogonal matrix of dimension \code{d}-\code{d}.
#'
#' @return The array of transformed matrices with dimension \code{p}-\code{d}-\code{d}.
#'
#' @keywords internal
updateList = function(A, Q){
  p = dim(A)[1]
  d = dim(Q)[2]
  for (i in 1:p) {
    A[i,1:d,1:d] = t(Q) %*% A[i,,] %*% Q
  }
  A = A[,1:d,1:d]
  return(A)
}
