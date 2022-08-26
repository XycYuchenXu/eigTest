#' Update the list of matrices
#'
#' Calculate the list of matrices after orthogonal transformation by matrix \code{Q}
#'
#' @param A Original array of matrices
#' @param Q Orthogonal matrix
#'
#' @return Array of transformed matrices
#' @export
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
