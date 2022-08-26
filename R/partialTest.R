#' Test Subset of Common Eigenvectors
#'
#' Test whether a subset of common eigenvectors are shared by the list of random matrices.
#'
#' @param A Array of original matrices
#' @param cov.arr Array of original covariance matrices, default will use identity matrices
#' @param k Size of the block matrices to be tested
#' @param cn Constant for convergence
#' @param Q Schur components for transformation
#' @param d Size of original matrices
#' @param p Number of matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
#' @param param.out The parameters of limiting distribution should be output or not
#' @param nn Logical whether the eigenvector elements are nonnegative
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Must be supplied when \code{testType = 'chi'}
#'
#' @return P-value or a list of test information.
#' \itemize{
#' \item testType Test method, \code{'chi'} or \code{'gam'}.
#' \item statistic Test statistic.
#' \item df Degrees of freedom for chi-squared distribution.
#' \item shape Shape parameter in gamma distribution.
#' \item rate Rate parameter in gamma distribution.
#' \item pvalue P-value.
#' }
#' @export
#'
#' @importFrom 'MASS' ginv
#'
#' @examples partialTest(countryCoeff, countryCovar, k = 2, cn = sqrt(112), testType = 'gam')
partialTest = function(A, cov.arr = NULL, k, cn, eps, nn = FALSE,
                       Q = expmPartSchur(A, k, nn = nn), d = dim(A)[2],
                       p = dim(A)[1], testType = c('chi', 'gam'),
                       param.out = FALSE){

  if (is.null(cov.arr)) {
    cov.arr = array(0, c(p, d^2, d^2))
    for (i in 1:p) {
      cov.arr[i,,] = diag(d^2)
    }
  }

  if (k >= d) {k = d}
  if (k <= 0) {k = d}
  if (k == d) {
    return(eigTest(A, cov.arr, cn, eps, testType = testType, param.out = param.out))
  }
  if (k == 1) {
    return(schurTest(A, cov.arr, k, cn, eps, testType = testType, param.out = param.out, nn = nn))
  }

  matB = matrix(0, d, d)
  matC = matrix(0, d, d)
  matB[1:k, 1:k] = 1
  SB = diag(as.double(matB))
  SB = SB[which(as.double(matB) == 1),]
  matC[1:k, (k+1):d] = 1
  SC = diag(as.double(matC))
  SC = SC[which(as.double(matC) == 1),]
  SBC = rbind(SB, SC)

  SQ = SBC %*% kronecker(t(Q), t(Q))

  B = array(0, c(p,k,k))
  for (i in 1:p) {
    Ai = t(Q) %*% A[i,,] %*% Q
    B[i,,] = Ai[1:k, 1:k]
  }
  V = JDTE(B)
  SV = (diag(k^2) - diag(as.vector(diag(k)))) %*% kronecker(t(V), ginv(V))

  SW = matrix(0, nrow = k*d, ncol = d*k)
  SW[1:(k^2), 1:(k^2)] = SV
  SW[(k^2+1):(k*d), (k^2+1):(k*d)] = diag(k*(d-k))

  Varr = array(0, c(p, k*d))
  for (i in 1:p) {
    Varr[i,] = SW %*% SQ %*% as.double(A[i,,])
    cov.arr[i, 1:(k*d), 1:(k*d)] = tcrossprod(SW %*% SQ, SW %*% SQ %*% cov.arr[i,,])
  }
  cov.arr = cov.arr[,1:(k*d),1:(k*d)]

  if (length(testType) == 2) {
    output = c(vec.test(Varr, cov.arr, cn, eps, 'chi')$pvalue, vec.test(Varr, cov.arr, cn, eps, 'gam')$pvalue)
    return(output)
  } else {
    testResult = vec.test(Varr, cov.arr, cn, eps, testType)
    if (param.out) {return(testResult)}
    return(testResult$pvalue)
  }

}
