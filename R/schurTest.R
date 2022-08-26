#' Test Schur Components
#'
#' @param A Array of matrices
#' @param cov.arr Array of covariance matrices, default will use identity matrices
#' @param k Number of components to be tested
#' @param cn Constant for convergence
#' @param Q Orthogonal components to be tested
#' @param d Size of matrices
#' @param p Number of matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
#' @param param.out The parameters of limiting distribution should be output or not
#' @param nn Logical whether the eigenvector elements are nonnegative
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
#' @examples schurTest(countryCoeff, countryCovar, k = 2, cn = sqrt(112), eps = 112^(-1/3), testType = 'chi')
schurTest = function(A, cov.arr = NULL, k, cn, eps=NULL, nn = FALSE, Q = NULL, d = dim(A)[2], p = dim(A)[1], testType = c('gam', 'chi'), param.out = FALSE){

  if (is.null(Q)) {Q = partSchur(A, k, nonneg = nn)}
  if (is.null(cov.arr)) {
    cov.arr = array(0, c(p, d^2, d^2))
    for (i in 1:p) {
      cov.arr[i,,] = diag(d^2)
    }
  }
  if (k >= d) {return(eigTest(A, cov.arr, cn, eps, testType, param.out))}
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

  if (length(testType) == 2) {
    output = c(vec.test(Varr, cov.arr, cn, eps, 'chi')$pvalue, vec.test(Varr, cov.arr, cn, eps, 'gam')$pvalue)
    return(output)
  } else {
    testResult = vec.test(Varr, cov.arr, cn, eps, testType)

    if (param.out) {return(testResult)}
    return(testResult$pvalue)
  }

}
