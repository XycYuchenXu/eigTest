#' Initiate orthogonal matrix \code{Q}.
#'
#' @param A Array of matrices
#'
#' @return Best initial condition for function \code{partSchur}.
#' @export
#'
#' @keywords internal
#'
init.Q = function(A){

  d = dim(A)[2]
  p = dim(A)[1]

  if (d > 3) {
    scores = c()
    Q = array(0, dim = c(p,d,d))
    Q.temp = diag(d)

    for (i in 1:p) {
      Q[i,,] = qr.Q(qr(Re(eigen(A[i,,])$vectors)))
      Qi = Q[i,,]
      scores = c(scores, score.fun(A, Qi))
    }
    k = which.min(scores)[1]
    Q = Q[k,,]
  } else {
    Q = diag(d)
  }

  return(Q)
}
