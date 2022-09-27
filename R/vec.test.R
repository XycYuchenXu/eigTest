#' Calculate the p-value for normal random vectors with specified covariance matrices
#'
#' @param V.arr The array of vectors, with dimension \code{p}-\code{L}, where \code{p} is the number of vectors, \code{L} is the vector length.
#' @param cn The convergence rate(s) to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators. If \code{length(cn) < p}, all vectors share the same rate \code{cn[1]}, otherwise \code{cn = cn[1:p]}.
#' @param testType The test methods, can be exact chi-squared test \code{testType = 'chi'}, and/or approximated gamma test \code{testType = 'gam'}.
#' @param cov.arr The array of covariance matrices corresponding to the vectors with dimension \code{p}-\code{L^2}-\code{L^2}, default will use identity matrices when \code{is.null(cov.arr)}.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.
#' @param param.out Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output P-value.
#'
#' @return A P-value when \code{param.out=FALSE} or a list of test information when \code{param.out = TRUE}.
#' \itemize{
#' \item \code{testType}: The test methods. Either exact chi-squared test \code{'chi'}, or approximated gamma test \code{'gam'}.
#' \item \code{statistic}: The test statistic.
#' \item \code{df}: The degrees of freedom for chi-squared distribution when \code{testType = 'chi'}.
#' \item \code{shape}: The shape parameter in gamma distribution when \code{testType = 'gam'}.
#' \item \code{rate}: The rate parameter in gamma distribution when \code{testType = 'gam'}.
#' \item \code{pvalue}: The P-value.
#' }
#'
#' @keywords internal
#'
#' @import stats
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
    testResult = list(testType, testVal, r, as.numeric(1 - pchisq(testVal, r)))
    names(testResult) = c('testType', 'statistic', 'df', 'pvalue')
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
    testResult = list(testType, testVal, alphaP, betaP, as.numeric(1 - pgamma(testVal, shape = alphaP, rate = betaP)))
    names(testResult) = c('testType', 'statistic', 'shape', 'rate', 'pvalue')
  }

  if (param.out) {return(testResult)}
  return(testResult$pvalue)
}
