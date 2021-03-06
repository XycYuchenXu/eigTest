#' Test on the commutator of two matrices
#'
#' @param matList The list of two matrices to be tested
#' @param cn Constant for convergence
#' @param d Size of matrices
#' @param covList The list of covariance matrices of the two matrices, default will use identity matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
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
#' @examples means = generateMeans(5,2)
#' samples = simuSamples(means, sqrt(100), 1)
#' commutatorTest(samples[[1]][[1]][[1]][[1]], sqrt(400), 400^(-1/3))
commutatorTest = function(matList, cn, eps=NULL, d = nrow(A),
                          covList = list(diag(d^2), diag(d^2)),
                          testType = c('chi', 'gam')) {

  if (length(testType) == 2){
    testType = 'chi'
  }

  A = matList[[1]]
  B = matList[[2]]
  covA = covList[[1]]
  covB = covList[[2]]

  y = matrix(data = A %*% B - B %*% A, ncol = 1)

  QA = kronecker(diag(d), A) - kronecker(t(A), diag(d))
  QB = kronecker(diag(d), B) - kronecker(t(B), diag(d))

  sigma.y = tcrossprod(QA, QA %*% covB) + tcrossprod(QB, QB %*% covA)

  return(vec.test(list(y), list(sigma.y), cn, eps, testType))
}
