#' Partial Joint Schur Decomposition
#'
#' @param A List of matrices
#' @param k Number of Schur components
#' @param n Size of matrices
#' @param p Number of matrices
#' @param iter Maximum iteration number
#' @param tol Tolerance error
#' @param nonneg Logical whether the eigenvector elements are nonnegative
#'
#' @return Orthogonal matrix \code{Q}
#' @export
#'
#' @importFrom Matrix Schur
#' @importFrom MASS ginv
#' @importFrom pracma quadprog
#'
#' @examples partSchur(countryCoeff, k = 2)
partSchur = function(A, k, n = ncol(A[[1]]), p = length(A), iter = 5000, tol = 10^(-16), nonneg = FALSE){

  if (nonneg == TRUE) {
    mid = matrix(0, ncol = n, nrow = n)
    for (i in 1:p) {
      AiI = A[[i]] - diag(n)
      mid = mid + tcrossprod(AiI)
    }
    Q = quadprog(mid, d = rep(0,n), Aeq = rep(1, n), beq = 1, lb = rep(0,n))$xmin
    Q = cbind(Q, matrix(rnorm(n*(n-1)), nrow = n))
    Q = qr.Q(qr(Q))
    return(Q)
  }

  for (i in 1:p) {
    A[[i]] = t(A[[i]])
  }
  if (k >= n) {k = n}

  oneSchur = function(A, n = ncol(A[[1]]), Q = diag(n)){
    Qr = cbind(Q[,1], Q[,2])

    score.old = score.fun(A, Q)

    for (i in 1:iter) {
      seq.order = n:2
      for (r in seq.order) {
        Qr[,2] = Q[,r]

        mid = matrix(0, nrow = 2, ncol = 2)
        for (j in 1:p) {
          Aj = A[[j]]
          mid = mid + tcrossprod(crossprod(Qr, Aj %*% Q[,1]))
          for (k in 2:n) {
            mid = mid - tcrossprod(crossprod(Aj %*% Qr, Q[,k]))
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
    n = ncol(A[[1]])
    Qi = Q

    for (j in 1:k) {
      if (j == n) {break()}
      Ai = updateList(A, Qi)
      D = diag(n - j + 1)
      scores = score.fun(Ai, D)
      for (i in (j+1):n) {
        Di = D
        Di[c(1,i-j+1), c(1,i-j+1)] = matrix(c(0,1,1,0), ncol = 2)
        scores = c(scores, score.fun(Ai, Di))
      }
      kk = which.min(scores)[1]
      Di = D
      if (kk != 1) {
        Di[c(1,kk), c(1,kk)] = matrix(c(0,1,1,0), ncol = 2)
      }
      Q[,j:n] = Q[,j:n] %*% oneSchur(Ai, Q = Di)
      Qi = Q[,-(1:j)]
    }
    return(Q)
  }

  Ai = A
  Q = orderUpdate(A, diag(n))
  Q = orderUpdate(A, Q)
  return(Q)

}
