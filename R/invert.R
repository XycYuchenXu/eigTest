#' Matrix inverse wrapper
#'
#' Prefer \code{solve}. When singular error is encountered, use \code{ginv} instead.
#'
#' @param X The square matrix to be inverted.
#'
#' @return The matrix inverse of \code{X}.
#'
#' @importFrom 'MASS' ginv
#' @keywords internal
#'
invert = function(X){
  try( {X_inv = solve(X); return(X_inv)}, silent = T )
  return(ginv(X))
}
