#' Score function (first semi-column squared norm)
#'
#' Calculate the summation of squared 'F' norms of \code{Bi[(j+1):n,j]} where \code{Bi = t(Q) Ai Q}.
#'
#' @param A List of matrices
#' @param Q The orthogonal matrix
#' @param j The column to be summed
#'
#' @return The sum of squared norms
#'
#' @keywords internal
#'
score.fun = function(A, Q){
  score.fun = 0
  n = ncol(A[[1]])
  p = length(A)

  for (i in 1:p) {
    Ai = crossprod(Q, A[[i]] %*% Q)
    vi = Ai[-1,1]
    score.fun = score.fun + as.double(crossprod(vi))
  }
  return(score.fun)
}
