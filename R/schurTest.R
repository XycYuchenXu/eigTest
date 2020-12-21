#' Test Schur Components
#'
#' @param A List of matrices
#' @param covList List of covariance matrices, default will use identity matrices
#' @param k Number of components to be tested
#' @param cn Constant for convergence
#' @param Q Orthogonal components to be tested
#' @param n Size of matrices
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
schurTest = function(A, covList = list(), k, cn, eps=NULL, nn = FALSE, Q = NULL, n = ncol(A[[1]]), p = length(A), testType = c('gam', 'chi'), param.out = FALSE){

  if (is.null(Q)) {Q = partSchur(A, k, nonneg = nn)}
  if (length(covList) == 0) {
    covList = vector('list', p)
    for (i in 1:p) {
      covList[[i]] = diag(n^2)
    }
  }
  if (k >= n) {return(eigTest(A, covList, cn, eps, testType, param.out))}
  Mat = matrix(0, nrow = n, ncol = n)
  Mat[1:k,(k+1):n] = 1
  Mat = as.double(Mat)
  select = diag(Mat)
  S = select[which(Mat == 1),]
  SV = S %*% kronecker(t(Q), t(Q))

  Vlist = vector('list', p)
  for (i in 1:p) {
    Vlist[[i]] = S %*% as.double(crossprod(Q, A[[i]] %*% Q))
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
