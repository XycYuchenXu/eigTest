#' Test of common Schur components
#'
#' @param A Array of matrices
#' @param cn Constant for convergence
#' @param cov.arr Array of covariance matrices, default will use identity matrices
#' @param nn Logical whether the eigenvector elements are nonnegative
#' @param k Number of components to be tested, can be \code{NULL} when \code{nn=TRUE}
#' @param warmup Logical whether use \code{partSchur} for a warm-up \code{Q}
#' @param Q Orthogonal components to be tested
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
#' @examples schurTest(countryCoeff, cn = sqrt(112), cov.arr = countryCovar, k = 2, testType = 'chi')
schurTest = function(A, cn, cov.arr = NULL, nn = FALSE, k = NULL, warmup = FALSE, Q = NULL,
                     testType = c('gam', 'chi'), eps = NULL, param.out = FALSE){

  d = dim(A)[2]; p = dim(A)[1]
  if (nn) {
    if (!is.null(k) && k != 1) {
      print('Unavailable setup. The number of components must be 1 for nonnegative test.')
      return()
    }
    k = 1
  }
  if (is.null(k) || k >= d || k <= 0 || k != round(k)) {
    print('Unavailable setup. The number of components must be an integer within (0, d) for Schur test.')
    return()
  }

  if (is.null(Q)) {
    if (nn) {Q = nnPartSchur(A)}
    else {Q = expmPartSchur(A, k, warmup = warmup)}
  }

  if (is.null(cov.arr)) {
    cov.arr = array(0, c(p, d^2, d^2))
    for (i in 1:p) {
      cov.arr[i,,] = diag(d^2)
    }
  }

  Mat = matrix(0, d, d)
  Mat[1:k,(k+1):d] = 1
  Mat = as.double(Mat)
  select = diag(Mat)
  S = select[which(Mat == 1),]
  SV = S %*% kronecker(t(Q), t(Q))

  Varr = array(0, c(p, k*(d-k)))
  for (i in 1:p) {
    Varr[i,] = S %*% as.double(crossprod(Q, A[i,,] %*% Q))
    cov.arr[i,1:(k*(d-k)),1:(k*(d-k))] = tcrossprod(SV, SV %*% cov.arr[i,,])
  }
  cov.arr = cov.arr[,1:(k*(d-k)),1:(k*(d-k))]

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
