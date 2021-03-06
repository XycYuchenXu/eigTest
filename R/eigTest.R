#' Test Whole Set of Eigenvectors
#'
#' @param A List of matrices
#' @param V Eigenvector to be tested
#' @param covList List of covariance matrices, default with use identity matrices
#' @param n Size of matrices
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
eigTest = function(A, covList = list(), cn, eps=NULL, V = JDTE(A), n = ncol(A[[1]]),
                   p = length(A), testType = c('chi', 'gam'), param.out = FALSE){

  if (length(covList) == 0) {
    covList = vector('list', p)
    for (i in 1:p) {
      covList[[i]] = diag(n^2)
    }
  }
  S = diag(n^2) - diag(as.vector(diag(n)))
  S = S[-which(as.double(diag(n)) == 1),]
  SV = S %*% kronecker(t(V), ginv(V))

  Vlist = vector('list', p)
  for (i in 1:p) {
    Vlist[[i]] = SV %*% as.double(A[[i]])
    covList[[i]] = tcrossprod(SV, SV %*% covList[[i]])
  }

  if (length(testType) == 2) {
    output = c(vec.test(Vlist, covList, cn, eps, 'chi')$pvalue, vec.test(Vlist, covList, cn, eps, 'gam')$pvalue)
    return(output)
  } else {
    testResult = vec.test(Vlist, covList, cn, eps, testType)

    if (param.out) {return(testResult)}
    return(testResult$pvalue)
  }

}
