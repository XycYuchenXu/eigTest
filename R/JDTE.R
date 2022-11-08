#' JDTE algorithm.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param iter The maximum iteration times.
#' @param tol The tolerance error for iteration termination.
#'
#' @return The eigenvector matrix \code{U} with dimension \code{d}-\code{d}.
#'
#' @importFrom 'Rdpack' reprompt
#' @export
#'
#' @references
#'     \insertRef{xu2021testing}{eigTest}
#'     \insertRef{andre}{eigTest}
#'
#' @examples JDTE(countryCoeff)
JDTE = function(A, iter = 500, tol = 10^(-8)){

  d = dim(A)[2]; p = dim(A)[1]
  U = eigen(colSums(A), symmetric = T)$vectors
  tempA = A
  for (i in 1:p) {tempA[i,,] = tcrossprod(crossprod(U, tempA[i,,]), t(U))}

  score.old = score.fun(tempA)
  for (i in 1:iter) {
    Z = matrix(0, d, d)
    for (r in 1:d) {
      for (s in 1:d) {
        if (r == s) {next()}
        denum = sum((tempA[,r,r] - tempA[,s,s])^2)
        if (denum > 0) {
          Z[r,s] = - sum(tempA[,r,s] * (tempA[,r,r] - tempA[,s,s]))/denum
        }
      }
    }

    wdenum = 0; wnum = 0
    for (j in 1:p) {
      Oj = t(tempA[j,,])
      Cj = crossprod(Oj, Z) - tcrossprod(Z, Oj); diag(Cj) = 0
      wdenum = wdenum + sum(Cj^2)
      wnum = wnum + sum(t(Oj) * Cj)
    }
    if (wdenum > 0) {
      coeff = wnum / wdenum; M = 1 / max(colSums(abs(Z)), rowSums(abs(Z)))
      if (abs(coeff) > M) {coeff = M * sign(coeff)}
      Z = diag(d) - coeff * Z; Zinv = solve(Z)
      U = crossprod(t(U), Z)
      for (j in 1:p) {tempA[j,,] = crossprod(t(Zinv), tcrossprod(tempA[j,,], t(Z)))}
    }

    score.new = score.fun(tempA)
    if (abs(score.new - score.old) < tol) {break()}
    score.old = score.new
  }
  U
}
