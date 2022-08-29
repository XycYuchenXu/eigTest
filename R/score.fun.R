#' Score function (first semi-column squared norm)
#'
#' Calculate the summation of squared 'F' norms of \code{Bi[(j+1):n,j]} where \code{Bi = t(Q) Ai  Q}.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param Q The orthogonal matrix of dimension \code{d}-\code{d} with the first column evaluated for common.
#'
#' @return The sum of squared norms of \code{t(Q[,-1]) Ai Q[,1]}.
#'
#' @keywords internal
#'
score.fun = function(A, Q){

  score.fun = 0
  d = dim(A)[2]
  p = dim(A)[1]

  for (i in 1:p) {
    score.fun = score.fun +
      norm(crossprod(Q[,-1], A[i,,] %*% Q[,1]), 'F')^2
  }
  return(score.fun)
}
