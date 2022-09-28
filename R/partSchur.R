#' Partial Joint Schur Decomposition
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param k The number of common components to be tested. \code{k} must be an integer within (0, \code{d}).
#' @param iter The maximum iteration number.
#' @param tol The tolerance error for iteration termination.
#'
#' @return The orthogonal matrix \code{Q} of dimension \code{d}-\code{d} with the first \code{k} columns to be the estimated common components.
#'
#' @keywords internal
#'
partSchur = function(A, k, iter = 5000, tol = 10^(-16)){

  d = dim(A)[2]; p = dim(A)[1]
  if (k <= 0 || k >= d || k != round(k)) {
    print('Unavailable setup. The number of components must be an integer within (0, d).')
    return()
  }

  for (i in 1:p) {
    A[i,,] = t(A[i,,])
  }
  if (k >= d) {k = d}

  oneSchur = function(A, d = dim(A)[2], Q = diag(d)){
    Qr = cbind(Q[,1], Q[,2])

    score.old = score.fun(A, Q)

    for (i in 1:iter) {
      for (r in d:2) {
        Qr[,2] = Q[,r]

        mid = matrix(0, nrow = 2, ncol = 2)
        for (j in 1:p) {
          Mj = t(crossprod(Qr, A[j,,]))
          mid = mid + tcrossprod(crossprod(Mj, Q[,1]))
          for (k in 2:d) {
            mid = mid - tcrossprod(crossprod(Mj, Q[,k]))
          }
        }

        Qr = crossprod(t(Qr), eigen(mid, symmetric = TRUE)$vectors)
        Q[,c(1,r)] = Qr
      }

      score.new = score.fun(A, Q)
      if (abs(score.new - score.old) < tol) {break()}
      score.old = score.new
    }
    return(Q)
  }

  orderUpdate = function(A, Q){
    d = dim(A)[2]
    Qi = Q

    for (j in 1:k) {
      if (j == d) {break()}
      Ai = updateList(A, Qi)
      D = diag(d - j + 1)
      minS = score.fun(Ai, D); minD = D
      for (i in (j+1):d) {
        Di = D
        Di[c(1,i-j+1), c(1,i-j+1)] = matrix(c(0,1,1,0), ncol = 2)
        tempS = score.fun(Ai, Di)
        if (tempS < minS) {
          minS = tempS; minD = Di
        }
      }
      Q[,j:d] = crossprod(t(Q[,j:d]), oneSchur(Ai, Q = minD))
      Qi = Q[,-(1:j)]
    }
    return(Q)
  }

  Q = orderUpdate(A, diag(d))
  Q = orderUpdate(A, Q)
  return(Q)

}
