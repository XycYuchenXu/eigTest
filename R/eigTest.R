#' Test Whole Set of Eigenvectors
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param cn The convergence rate(s) to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators. If \code{length(cn) < p}, all matrices share the same rate \code{cn[1]}, otherwise \code{cn = cn[1:p]}.
#' @param cov.arr The array of covariance matrices corresponding to the vectors with dimension \code{p}-\code{d^2}-\code{d^2}, default will use identity matrices when \code{is.null(cov.arr)}.
#' @param V The eigenvectors to be tested with dimension \code{d}-\code{d}. Default will use \code{V = JDTE(A)} when \code{is.null(V)}.
#' @param testType The test methods, can be exact chi-squared test \code{testType = 'chi'}, and/or approximated gamma test \code{testType = 'gam'}.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.
#' @param param.out Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output P-value.
#'
#' @return A named vector of P-value(s) when \code{param.out = FALSE}, or a named list of test information when \code{param.out = TRUE}, with name(s) to be \code{'chi'} and/or \code{'gam'}. Test information is a sub-list with elements including:
#' \itemize{
#' \item \code{testType}: The test methods. Either exact chi-squared test \code{'chi'}, or approximated gamma test \code{'gam'}.
#' \item \code{statistic}: The test statistic.
#' \item \code{df}: The degrees of freedom for chi-squared distribution when \code{testType = 'chi'}.
#' \item \code{shape}: The shape parameter in gamma distribution when \code{testType = 'gam'}.
#' \item \code{rate}: The rate parameter in gamma distribution when \code{testType = 'gam'}.
#' \item \code{pvalue}: The P-value.
#' }
#' @export
#'
#' @importFrom 'MASS' ginv
#'
#' @examples eigTest(countryCoeff, cn = sqrt(112), cov.arr = countryCovar, testType = 'gam')
eigTest = function(A, cn, cov.arr = NULL, V = NULL, testType = c('chi', 'gam'),
                   eps = NULL, param.out = FALSE){

  d = dim(A)[2]; p = dim(A)[1]
  if (is.null(cov.arr)) {
    cov.arr = array(0, c(p, d^2, d^2))
    for (i in 1:p) {
      cov.arr[i,,] = diag(d^2)
    }
  }
  if (is.null(V)) {V = JDTE(A)}

  S = diag(d^2) - diag(as.vector(diag(d)))
  S = S[-which(as.double(diag(d)) == 1),]
  SV = S %*% kronecker(t(V), ginv(V))

  Varr = array(0, c(p, d^2 - d))
  for (i in 1:p) {
    Varr[i,] = SV %*% as.double(A[i,,])
    cov.arr[i, 1:(d^2 - d), 1:(d^2 - d)] = tcrossprod(SV, SV %*% cov.arr[i,,])
  }
  cov.arr = cov.arr[, 1:(d^2 - d), 1:(d^2 - d)]

  if (param.out) {output = list()}
  else {output = c()}
  for (tp in testType) {
    temp = vec.test(Varr, cn, tp, cov.arr, eps, param.out)
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
