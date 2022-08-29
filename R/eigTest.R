#' Test Whole Set of Eigenvectors
#'
#' @param A Array of matrices
#' @param cn Constant for convergence
#' @param cov.arr List of covariance matrices, default with use identity matrices
#' @param V Eigenvector to be tested
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test.
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

  output = c();
  for (tp in testType) {
    temp = vec.test(Varr, cn, tp, cov.arr, eps, param.out)
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
