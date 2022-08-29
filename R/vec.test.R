#' Calculate the p-value for normal random vectors with specified covariance matrices
#'
#' @param V.arr Array of vectors
#' @param cn The constant for convergence, or a vector of constants.
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
#' @param cov.arr Array of covariance matrices corresponding to the vectors, default will use identity matrices
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied
#' @param param.out The parameters of limiting distribution should be output or not
#'
#' @return A list of test information.
#' \itemize{
#' \item statistic Test statistic.
#' \item df Degrees of freedom for chi-squared distribution.
#' \item shape Shape parameter in gamma distribution.
#' \item rate Rate parameter in gamma distribution.
#' \item pvalue P-value.
#' }
#'
#' @keywords internal
#'
#' @import 'stats'
#'
vec.test = function(V.arr, cn, testType, cov.arr = NULL,
                    eps=NULL, param.out = FALSE){

  if (testType != 'chi' && testType != 'gam') {
    print('Unavailable setup. Only support "chi" for chi-test and "gam" for gamma-test.')
    return()
  }

  p = dim(V.arr)[1]
  if (is.null(cov.arr)) {
    n = dim(V.arr)[2]
    cov.arr = array(0, c(p, n, n))
    for (i in 1:p) {
      cov.arr[i,,] = diag(n)
    }
  }

  if (length(cn) < p) {
    cn = rep(cn[1], p)
  } else {
    cn = cn[1:p]
  }
  if (is.null(eps)) {eps = cn^(-2/3)}
  else if (length(eps) < p) {eps = rep(eps[1], p)}
  else {eps = eps[1:p]}

  testVal = 0
  if (testType == 'chi') {
    r = 0
    for (i in 1:p) {
      vi = V.arr[i,]
      s = truncateSVD(cov.arr[i,,], eps[i])
      testVal = testVal + crossprod(vi, s$ginv %*% vi) * cn[i]^2
      r = r + s$r
    }
    testResult = list(testVal, r, 1 - pchisq(testVal, r))
    names(testResult) = c('statistic', 'df', 'pvalue')
  } else {
    me = 0
    va = 0
    for (i in 1:p) {
      vi = V.arr[i,]
      testVal = testVal + norm(vi, '2')^2 * cn[i]^2
      sig.i = cov.arr[i,,]
      me = me + sum(diag(sig.i))
      va = va + 2*sum(diag(crossprod(sig.i)))
    }
    alphaP = me^2/va
    betaP = me/va
    testResult = list(testVal, alphaP, betaP, 1 - pgamma(testVal, shape = alphaP, rate = betaP))
    names(testResult) = c('statistic', 'shape', 'rate', 'pvalue')
  }

  if (param.out) {return(testResult)}
  return(testResult$pvalue)
}
