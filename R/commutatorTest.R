#' Test on the commutator of two matrices
#'
#' @param mat.arr The array of two matrices to be tested
#' @param cn Constant for convergence
#' @param cov.arr The array of covariance matrices of the two matrices, default will use identity matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied
#' @param param.out The parameters of limiting distribution should be output or not
#'
#' @return A named vector of P-value(s) when \code{param.out = FALSE}, or a named list of test information when \code{param.out = TRUE}, with name(s) to be \code{'chi'} and/or \code{'gam'}. Test information is a sub-list with elements including:
#' \itemize{
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
#' @examples means = generateMeans(5,2)
#' samples = simuSamples(means, sqrt(100), 1)
#' commutatorTest(samples[[1]][1,1,1,,,], sqrt(400))
commutatorTest = function(mat.arr, cn, cov.arr = NULL, testType = c('chi', 'gam'),
                          eps = NULL, param.out = FALSE) {

  d = dim(mat.arr)[2]
  if (length(testType) == 2){
    testType = 'chi'
  }

  if (is.null(cov.arr)) {
    covA = diag(d^2)
    covB = diag(d^2)
  } else {
    covA = cov.arr[1,,]
    covB = cov.arr[2,,]
  }

  A = mat.arr[1,,]
  B = mat.arr[2,,]

  y = array(A %*% B - B %*% A, dim = c(1, d^2))

  QA = kronecker(diag(d), A) - kronecker(t(A), diag(d))
  QB = kronecker(diag(d), B) - kronecker(t(B), diag(d))

  sigma.y = array(tcrossprod(QA, QA %*% covB) + tcrossprod(QB, QB %*% covA), dim = c(1, d^2, d^2))

  output = c();
  for (tp in testType) {
    temp = vec.test(y, cn, tp, sigma.y, eps, param.out)
    if (!is.null(temp)) {
      output = c(output, temp); names(output)[length(output)] = tp
    }
  }
  if (is.null(output)) {
    print('Unavailable test type. Please try using "chi" and/or "gam".')
    return()
  }
  return(output)
}
