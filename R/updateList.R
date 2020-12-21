#' Update the list of matrices
#'
#' Calculate the list of matrices after orthogonal transformation by matrix \code{Q}
#'
#' @param A Original list of matrices
#' @param Q Orthogonal matrix
#'
#' @return LIst of transformed matrices
#' @export
#'
#' @keywords internal
updateList = function(A, Q){
  p = length(A)
  for (i in 1:p) {
    A[[i]] = t(Q) %*% A[[i]] %*% Q
  }
  return(A)
}
