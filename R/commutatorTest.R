#' Test on the Commutator for Two Matrices
#'
#' Two matrices with common eigenvectors commute.
#' See \insertCite{xu2021testing;textual}{eigTest}.
#'
#' @param mat.arr The array of two matrices to be tested with dimension \code{2}-\code{d}-\code{d}. If \code{dim(mat.arr)[1] > 2}, only the first two are used.
#' @param cn The convergence rate(s) to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators. If \code{length(cn) < p}, all matrices share the same rate \code{cn[1]}, otherwise \code{cn = cn[1:p]}.
#' @param cov.arr The array of covariance matrices of the two matrices with dimension \code{2}-\code{d^2}-\code{d^2}, default will use identity matrices when \code{is.null(cov.arr)}. If \code{dim(cov.arr)[1] > 2}, only the first two are used.
#' @param testType The test methods, can be exact chi-squared test \code{testType = 'chi'}, and/or approximated gamma test \code{testType = 'gam'}.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.
#' @param param.out Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output p-value.
#'
#' @return A named vector of p-value(s) when \code{param.out = FALSE}, or a named list of test information when \code{param.out = TRUE}, with name(s) to be \code{'chi'} and/or \code{'gam'}. Test information is a sub-list with elements including:
#' \itemize{
#' \item \code{testType}: The test methods. Either exact chi-squared test \code{'chi'}, or approximated gamma test \code{'gam'}.
#' \item \code{statistic}: The test statistic.
#' \item \code{df}: The degrees of freedom for chi-squared distribution when \code{testType = 'chi'}.
#' \item \code{shape}: The shape parameter in gamma distribution when \code{testType = 'gam'}.
#' \item \code{rate}: The rate parameter in gamma distribution when \code{testType = 'gam'}.
#' \item \code{pvalue}: The p-value.
#' }
#' @export
#'
#' @import stats
#' @importFrom 'Rdpack' reprompt
#'
#' @references
#' \insertAllCited{}
#'
#' @examples means = generateMeans(5,2)
#' samples = simuSamples(means, sqrt(100), 1)
#' commutatorTest(samples[[1]]$mu.bar, cn = samples[[1]]$CovRate)
commutatorTest = function(mat.arr, cn, cov.arr = NULL, testType = c('chi', 'gam'),
                          eps = NULL, param.out = FALSE) {

  d = dim(mat.arr)[2]
  if (length(testType) == 2){
    testType = 'chi'
  }

  A = mat.arr[1,,]
  tB = t(mat.arr[2,,])

  y = matrix(tcrossprod(A, tB) - crossprod(tB, A), 1, d^2)

  if (is.null(cov.arr)) {
    sigma.y = array(
      kronecker(diag(d), tcrossprod(A) + crossprod(tB)) +
        kronecker(crossprod(A) + tcrossprod(tB), diag(d)) -
        kronecker(A, A) - kronecker(t(A), t(A)) -
        kronecker(tB, tB) - kronecker(t(tB), t(tB)),
      dim = c(1, d^2, d^2))
  } else {
    QA = kronecker(diag(d), A) - kronecker(t(A), diag(d))
    QB = kronecker(diag(d), t(tB)) - kronecker(tB, diag(d))
    sigma.y = array(tcrossprod(QA, tcrossprod(QA, cov.arr[2,,])) +
                      tcrossprod(QB, tcrossprod(QB, cov.arr[1,,])),
                    dim = c(1, d^2, d^2))
  }

  if (param.out) {output = list()}
  else {output = c()}
  for (tp in testType) {
    temp = vec.test(y, cn, tp, sigma.y, eps, param.out)
    if (!is.null(temp)) {
      if (param.out) {output[[length(output) + 1]] = temp}
      else {output = c(output, temp)}
      names(output)[length(output)] = tp
    }
  }
  if (is.null(output)) {
    print('Unavailable test type. Please try using "chi" and/or "gam".')
    return()
  }
  return(output)
}
