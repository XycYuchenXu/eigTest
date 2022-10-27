#' Tangent space optimization for partially joint Schur decomposition using matrix exponential
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param k The number of common Schur components. Must be an integer within (0, \code{d}).
#' @param warmup Logical, whether use \code{partSchur} for a warm-up initial value, default to \code{FALSE}.
#' @param iter The maximum iteration number.
#' @param tol The tolerance error for iteration termination.
#'
#' @return The orthogonal matrix \code{Q} of dimension \code{d}-\code{d} with the first \code{k} columns to be the estimated common components.
#' @export
#'
#' @importFrom 'Matrix' expm
#'
#' @examples expmPartSchur(countryCoeff, 2)
expmPartSchur = function(A, k, warmup = FALSE, iter = 5000, tol = 10^(-12)){

  d = dim(A)[2]; p = dim(A)[1]
  if (k <= 0 || k >= d || k != round(k)) {
    print('Unavailable setup. The number of components must be an integer within (0, d).')
    return()
  }

  if (warmup) {
    Ui = partSchur(A, k)
  } else {
    Ui = diag(d)
  }

  ZeroM = matrix(0, d, d)

  gridNodes = seq(0,1,0.02)

  CY = c()
  for (i in 1:d) {
    for (j in 1:d) {
      if (i != j) {
        ZeroM[i,j] = 1
        ZeroM[j,i] = -1
      }
      CY = cbind(CY, as.vector(ZeroM))
      ZeroM[i,j] = 0; ZeroM[j,i] = 0
    }
  }

  UL = ZeroM
  UL[1:k, (k+1):d] = 1
  UL = diag(as.vector(UL))

  listAB = function(U){
    matA = 0
    vecB = 0
    for (i in 1:p) {
      Mi = crossprod(U, tcrossprod(A[i,,], t(U)))
      Ti = kronecker(diag(d), Mi) - kronecker(t(Mi), diag(d))
      matA = matA + crossprod(Ti, crossprod(UL, Ti))
      vecB = vecB + crossprod(Ti, crossprod(UL, as.vector(Mi)))
    }
    matA = as.matrix(matA)
    return(matrix(crossprod(invert(crossprod(CY, crossprod(matA, CY))),
                            crossprod(CY, vecB)), ncol = d))
  }

  scoresUL = function(U){
    S = 0
    for (i in 1:p) {
      Mi = crossprod(U, tcrossprod(A[i,,], t(U)))
      S = S + sum(Mi[1:k, (k+1):d]^2)
    }
    return(S)
  }

  resultAB = function(U, X){
    minS = Inf; Z = X
    for (i in 1:length(gridNodes)) {
      Zi = as.matrix(crossprod(t(U), expm(gridNodes[i] * X)))
      tempS = scoresUL(Zi)
      if (tempS < minS) {
        minS = tempS; Z = Zi
      }
    }
    return(list(Z, minS))
  }

  score.old = scoresUL(Ui)
  for (i in 1:iter) {
    Xi = listAB(Ui)
    solAB = resultAB(Ui, Xi)
    Ui = solAB[[1]]
    score.new = solAB[[2]]
    if (abs(score.new - score.old) < tol) {break()}
    score.old = score.new
  }
  return(Ui)
}
