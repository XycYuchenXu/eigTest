#' Tangent space optimization for partial Schur decomposition
#'
#' @param A Array of matrices
#' @param k Number of Schur components
#' @param d Size of matrices
#' @param p Number of matrices
#' @param iter Maximum iteration number
#' @param tol Tolerance error
#' @param nn Logical whether the eigenvector elements are nonnegative
#'
#' @return Orthogonal matrix \code{Q}
#' @export
#'
#' @importFrom Matrix expm
#' @importFrom 'MASS' ginv
#'
#' @examples expmPartSchur(countryCoeff, 2)
expmPartSchur = function(A, k, d = dim(A)[2], p = dim(A)[1], iter = 5000, tol = 10^(-12), nn = FALSE){

  if (k <= 0) {k = d}
  if (k >= d) {return(JDTE(A))}

  Ui = partSchur(A, k, nonneg = nn)
  if (nn) {return(Ui)}
  ZeroM = matrix(0, d, d)

  gridNodes = seq(0,1,0.02)

  CY = c()
  for (i in 1:d) {
    for (j in 1:d) {
      Ci = ZeroM
      if (i != j) {
        Ci[i,j] = 1
        Ci[j,i] = -1
      }
      CY = cbind(CY, as.vector(Ci))
    }
  }

  ul = ZeroM
  ul[1:k, (k+1):d] = 1
  UL = diag(as.vector(ul))

  listAB = function(U){
    matA = 0
    vecB = 0
    for (i in 1:p) {
      Mi = crossprod(U, A[i,,] %*% U)
      Ti = kronecker(diag(d), Mi) - kronecker(t(Mi), diag(d))
      matA = matA + crossprod(Ti, UL %*% Ti)
      vecB = vecB + crossprod(Ti, UL %*% as.vector(Mi))
    }
    return(list(as.matrix(matA), vecB))
  }

  scoresUL = function(U){
    S = 0
    for (i in 1:p) {
      Mi = crossprod(U, A[i,,] %*% U)
      S = S + norm(UL %*% as.vector(Mi), 'F')^2
    }
    return(S)
  }

  resultAB = function(U, X){
    scores = c()
    for (i in 1:length(gridNodes)) {
      newU = as.matrix(U %*% expm(gridNodes[i] * X))
      scores = c(scores, scoresUL(newU))
    }
    minT = which.min(scores)
    finalU = as.matrix(U %*% expm(gridNodes[minT] * X))
    return(list(finalU, scores[minT]))
  }

  score.old = scoresUL(Ui)
  for (i in 1:iter) {
    AB = listAB(Ui)
    Xi = ginv(crossprod(CY, AB[[1]] %*% CY)) %*% crossprod(CY, AB[[2]])
    Xi = matrix(Xi, ncol = d)
    solAB = resultAB(Ui, Xi)
    Ui = solAB[[1]]
    score.new = solAB[[2]]
    if (abs(score.new - score.old) < tol) {break()}
    score.old = score.new
  }
  return(Ui)
}
