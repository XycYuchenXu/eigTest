#' Initiate orthogonal matrix \code{Q}.
#'
#' @param A List of matrices
#'
#' @return Best initial condition for function \code{partSchur}.
#' @export
#'
#' @keywords internal
#'
init.Q = function(A){

  n = ncol(A[[1]])
  p = length(A)

  if (n > 3) {
    scores = c()
    Q = array(0, dim = c(p,n,n))
    Q.temp = diag(n)

    for (i in 1:p) {
      Q[i,,] = qr.Q(qr(Re(eigen(A[[i]])$vectors)))
      Qi = Q[i,,]
      scores = c(scores, score.fun(A, Qi))
    }
    k = which.min(scores)[1]
    Q = Q[k,,]
  } else {
    Q = diag(n)
  }

  return(Q)
}
