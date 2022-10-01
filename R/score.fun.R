#' Score function (off-diagonal squared norm or first semi-column squared norm)
#'
#' Calculate either the summed squared off-diagonal 'F' norms of \code{Ai} or the summed squared 'F' norms of \code{Bi[(j+1):n,j]} where \code{Bi = t(Q) Ai  Q} when \code{Q} is supplied.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param Q Optional. The orthogonal matrix of dimension \code{d}-\code{d} with the first column evaluated for common. Default (\code{is.null(Q) = TRUE}) will be ignored for off-diagonal norm.
#'
#' @return Either the sum of squared off-diagonal norms of \code{Ai} or the sum of squared norms of \code{t(Q[,-1]) Ai Q[,1]}.
#'
#' @keywords internal
#'
score.fun = function(A, Q = NULL){

  score.fun = 0
  p = dim(A)[1]

  if (is.null(Q)) {
    for (i in 1:p) {
      Ai = A[i,,]; diag(Ai) = 0
      score.fun = score.fun + sum(Ai^2)
    }
  } else {
    for (i in 1:p) {
      score.fun = score.fun +
        sum(crossprod(Q[,-1], crossprod(t(A[i,,]), Q[,1]))^2)
    }
  }

  return(score.fun)
}
