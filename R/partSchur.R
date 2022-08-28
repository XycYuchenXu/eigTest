#' Partial Joint Schur Decomposition
#'
#' @param A Array of matrices
#' @param k Number of Schur components
#' @param d Size of matrices
#' @param p Number of matrices
#' @param iter Maximum iteration number
#' @param tol Tolerance error
#' @param nonneg Logical whether the eigenvector elements are nonnegative
#'
#' @return Orthogonal matrix \code{Q}
#' @export
#'
#' @importFrom Matrix Schur
#' @importFrom pracma quadprog
#'
#' @examples partSchur(countryCoeff, k = 2)
partSchur = function(A, k, d = dim(A)[2], p = dim(A)[1], iter = 5000, tol = 10^(-16), nonneg = FALSE){

  if (nonneg == TRUE) {
    mid = matrix(0, d, d)
    for (i in 1:p) {
      mid = mid + tcrossprod(A[i,,] - diag(d))
    }
    Q = quadprog(mid, d = rep(0, d), Aeq = rep(1, d), beq = 1, lb = rep(0, d))$xmin
    Q = cbind(Q, matrix(rnorm(d*(d-1)), nrow = d))
    Q = qr.Q(qr(Q))
    return(Q)
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
          Mj = crossprod(crossprod(Qr, A[j,,]))
          mid = mid + crossprod(Q[,1], Mj %*% Q[,1])
          for (k in 2:d) {
            mid = mid - crossprod(Q[,k], Mj %*% Q[,k])
          }
        }

        Qr = Qr %*% eigen(mid, symmetric = TRUE)$vectors
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
      minD = D
      minS = score.fun(Ai, D)
      for (i in (j+1):d) {
        Di = D
        Di[c(1,i-j+1), c(1,i-j+1)] = matrix(c(0,1,1,0), ncol = 2)
        tempS = score.fun(Ai, Di)
        if (tempS < minS) {
          minD = Di; minS = tempS
        }
      }
      Q[,j:d] = Q[,j:d] %*% oneSchur(Ai, Q = minD)
      Qi = Q[,-(1:j)]
    }
    return(Q)
  }

  Q = orderUpdate(A, diag(d))
  Q = orderUpdate(A, Q)
  return(Q)

}
