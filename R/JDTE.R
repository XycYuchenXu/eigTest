#' JDTE algorithm.
#'
#' @param A Array of matrices.
#' @param d Size of matrices.
#' @param p Number of matrices.
#' @param iter Maximum iteration times.
#' @param tol Tolerance error.
#'
#' @return Matrix \code{U}.
#' @export
#'
#' @importFrom 'MASS' ginv
#'
#' @examples JDTE(countryCoeff)
JDTE = function(A, d = dim(A)[2], p = dim(A)[1], iter = 5000, tol = 10^(-16)){

  score.fun = function(A, p){
    sc = 0
    for (i in 1:p) {
      Ai = A[i,,]
      sc = sc + norm(Ai - diag(diag(Ai)), type = 'f')^2
    }
    sc
  }

  sumA = A[1,,]
  for (i in 2:p) {
    sumA = sumA + A[i,,]
  }
  U = eigen(sumA/p, symmetric = TRUE)$vectors
  V = solve(U)
  for (i in 1:p) {
    A[i,,] = V %*% A[i,,] %*% U
  }

  tempA = A
  score.old = score.fun(A, p)
  for (i in 1:iter) {

    Z = diag(d)

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

    for (j in 1:p) {
      tempA[j,,] = ginv(Z) %*% tempA[j,,] %*% Z
    }

    U = U %*% Z
    #U = U/norm(U, type = '2')
    score.new = score.fun(tempA, p)
    if (abs(score.new - score.old) < tol) {
      break()
    }
    score.old = score.new
  }
  U
}
