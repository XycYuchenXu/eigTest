#' Calculate the p-value for normal random vectors with specified covariance matrices
#'
#' @param Vlist List of vectors
#' @param covList List of covariance matrices corresponding to the vectors, default will use identity matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
#' @param cn The constant for convergence, or a vector of constants.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Must be supplied when \code{testType = 'chi'}
#'
#' @return A list of test information.
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
#' @import 'stats'
#'
#' @examples p = length(countryCoeff)
#' vlist = vector('list', p)
#' for (i in 1:p) {
#'   vlist[[i]] = as.double(countryCoeff[[i]])
#' }
#' vec.test(vlist, countryCovar, cn = sqrt(112), eps = 112^(-1/3), testType = 'chi')
vec.test = function(Vlist, covList = list(), cn, eps=NULL, testType = c('chi', 'gam')){

  p = length(Vlist)
  if (length(covList) == 0) {
    n = length(Vlist[[1]])
    covList = vector('list', p)
    for (i in 1:p) {
      covList[[i]] = diag(n)
    }
  }

  if (length(cn) != p) {
    cn = rep(cn[1], p)
  }

  testVal = 0

  if (testType == 'chi') {
    r = 0
    for (i in 1:p) {
      vi = Vlist[[i]]
      s = truncateSVD(covList[[i]], eps)
      testVal = testVal + crossprod(vi, s$ginv %*% vi) * cn[i]^2
      r = r + s$r
    }
    testResult = list(testType, testVal, r, 1 - pchisq(testVal, r))
    names(testResult) = c('testType', 'statistic', 'df', 'pvalue')
  } else {
    me = 0
    va = 0
    for (i in 1:p) {
      vi = Vlist[[i]]
      testVal = testVal + norm(vi, 'F')^2 * cn[i]^2
      sig.i = covList[[i]]
      me = me + sum(diag(sig.i))
      va = va + 2*sum(diag(crossprod(sig.i)))
    }
    alphaP = me^2/va
    betaP = me/va
    testResult = list(testType, testVal, alphaP, betaP, 1 - pgamma(testVal, shape = alphaP, rate = betaP))
    names(testResult) = c('testType', 'statistic', 'shape', 'rate', 'pvalue')
  }

  return(testResult)
}
