#' Score function (off-diagonal squared norm or first semi-column squared norm)
#'
#' Calculate either the summed squared off-diagonal 'F' norms of \code{Ai} or the summed squared 'F' norms of \code{Bi[(j+1):n,j]} where \code{Bi = t(Q) Ai  Q} when \code{Q} is supplied.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param Q Optional. The orthogonal matrix of dimension \code{d}-\code{d} with the first \code{k} columns evaluated for common. Default (\code{is.null(Q) == TRUE}) will ignore \code{Q} for off-diagonal norm.
#' @param k Optional. The number of common orthogonal components to be evaluated. Only applicable when \code{is.null(Q) == FALSE}. Default (\code{is.null(k) == TRUE}) will use \code{k = 1}.
#'
#' @return Either the sum of squared off-diagonal norms of \code{Ai} or the sum of squared norms of \code{t(Q[,(k+1):d]) Ai Q[,1:k]}.
#'
#' @keywords internal
#'
score.fun = function(A, Q = NULL, k = NULL){

  score.fun = 0
  p = dim(A)[1]

  if (is.null(Q)) {
    for (i in 1:p) {
      Ai = A[i,,]; diag(Ai) = 0
      score.fun = score.fun + sum(Ai^2)
    }
  } else {
    d = dim(A)[2]
    if (is.null(k)) {k = 1}
    for (i in 1:p) {
      score.fun = score.fun +
        sum(crossprod(Q[,(k+1):d], crossprod(t(A[i,,]), Q[,1:k]))^2)
    }
  }

  return(score.fun)
}
