#' JDTE algorithm.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param iter The maximum iteration times.
#' @param tol The tolerance error for iteration termination.
#'
#' @return The eigenvector matrix \code{U} with dimension \code{d}-\code{d}.
#' @export
#'
#' @examples JDTE(countryCoeff)
JDTE = function(A, iter = 500, tol = 10^(-8)){

  d = dim(A)[2]; p = dim(A)[1]
  U = eigen(colSums(A), symmetric = T)$vectors
  tempA = A
  for (i in 1:p) {tempA[i,,] = tcrossprod(crossprod(U, tempA[i,,]), t(U))}

  score.old = score.fun(A)
  for (i in 1:iter) {

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
        if (denum > 0) {
          coeff = - num/denum
          U[,s] = U[,s] - num/denum * U[,r]
          tempA[,,s] = tempA[,,s] + coeff * tempA[,,r]
          tempA[,r,] = tempA[,r,] - coeff * tempA[,s,]
          tempA[,r,s] = tempA[,r,s] + coeff^2 * tempA[,s,r]
        }
      }
    }

    score.new = score.fun(tempA)
    if (abs(score.new - score.old) < tol) {break()}
    score.old = score.new
  }
  U
}
