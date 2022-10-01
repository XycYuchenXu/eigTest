#' JDTE algorithm.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param iter The maximum iteration times.
#' @param tol The tolerance error for iteration termination.
#'
#' @return The eigenvector matrix \code{U} with dimension \code{d}-\code{d}.
#' @export
#'
#' @importFrom 'MASS' ginv
#'
#' @examples JDTE(countryCoeff)
JDTE = function(A, iter = 5000, tol = 10^(-16)){

  d = dim(A)[2]; p = dim(A)[1]

  sumA = A[1,,]
  for (i in 2:p) {
    sumA = sumA + A[i,,]
  }
  U = eigen(sumA/p, symmetric = TRUE)$vectors
  V = solve(U)
  for (i in 1:p) {
    A[i,,] = tcrossprod(crossprod(t(V), A[i,,]), t(U))
  }

  tempA = A
  score.old = score.fun(A)
  for (i in 1:iter) {

    Z = matrix(0, d, d)

    for (r in 1:d) {
      for (s in 1:d) {
        if (r == s) {next()}
        denum = 0
        num = 0
        for (j in 1:p) {
          Aj = tempA[j,,]
          denum = denum + (Aj[r,r] - Aj[s,s])^2
          num = num + Aj[r,s]*(Aj[r,r] - Aj[s,s])
        }
        Z[r,s] = - num/denum
      }
    }

    wnum = 0; wdenum = 0
    for (j in 1:p) {
      Oj = t(tempA[j,,])
      Cj = tcrossprod(Z, Oj) - crossprod(Oj, Z); diag(Cj) = 0
      diag(Oj) = 0
      wdenum = wdenum + sum(Cj^2)
      wnum = wnum + sum(t(Oj) * Cj)
    }
    Z = diag(d) + Z * wnum / wdenum

    for (j in 1:p) {
      tempA[j,,] = tcrossprod(tcrossprod(ginv(Z), t(tempA[j,,])), t(Z))
    }

    U = tcrossprod(U, t(Z))
    score.new = score.fun(tempA)
    if (abs(score.new - score.old) < tol) {
      break()
    }
    score.old = score.new
  }
  U
}
