#' Test Whole Set of Eigenvectors
#'
#' @param A Array of matrices
#' @param V Eigenvector to be tested
#' @param cov.arr List of covariance matrices, default with use identity matrices
#' @param d Size of matrices
#' @param p Number of matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test.
#' @param cn Constant for convergence
#' @param param.out The parameters of limiting distribution should be output or not
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
#' @examples eigTest(countryCoeff, countryCovar, cn = sqrt(112), testType = 'gam')
eigTest = function(A, cov.arr = NULL, cn, eps=NULL, V = JDTE(A), d = dim(A)[2],
                   p = dim(A)[1], testType = c('chi', 'gam'), param.out = FALSE){

  if (is.null(cov.arr)) {
    cov.arr = array(0, c(p, d^2, d^2))
    for (i in 1:p) {
      cov.arr[i,,] = diag(d^2)
    }
  }
  S = diag(d^2) - diag(as.vector(diag(d)))
  S = S[-which(as.double(diag(d)) == 1),]
  SV = S %*% kronecker(t(V), ginv(V))

  Varr = array(0, c(p, d^2 - d))
  for (i in 1:p) {
    Varr[i,] = SV %*% as.double(A[i,,])
    cov.arr[i, 1:(d^2 - d), 1:(d^2 - d)] = tcrossprod(SV, SV %*% cov.arr[i,,])
  }
  cov.arr = cov.arr[, 1:(d^2 - d), 1:(d^2 - d)]

  if (length(testType) == 2) {
    output = c(vec.test(Varr, cov.arr, cn, eps, 'chi')$pvalue, vec.test(Varr, cov.arr, cn, eps, 'gam')$pvalue)
    return(output)
  } else {
    testResult = vec.test(Varr, cov.arr, cn, eps, testType)

    if (param.out) {return(testResult)}
    return(testResult$pvalue)
  }

}
