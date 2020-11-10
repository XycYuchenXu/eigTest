#' Test Whole set of Eigenvectors of Upper Left Block after Schur Decomposition
#'
#' @param A List of original matrices
#' @param covList List of original covariance matrices
#' @param k Size of the block matrices to be tested
#' @param cm Constant for convergence
#' @param Q Schur components for transformation
#' @param n Size of original matrices
#' @param p Number of matrices
#' @param testType The test methods. Either for exact chi-squared test, or approximated gamma test
#' @param param.out The parameters of limiting distribution should be output or not
#' @param nn Logical whether the eigenvector elements are nonnegative
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
#' @import 'MASS'
#'
#' @examples schurBlockTest(countryCoeff, countryCovar, k = 2, cm = 102, testType = 'gam')
schurBlockTest = function(A, covList = list(), k, cm, nn = FALSE,
                          Q = expmPartSchur(A,k, nn = nn), n = nrow(A[[1]]),
                          p = length(A), testType = c('chi', 'gam'),
                          param.out = FALSE){

  if (length(covList) == 0) {
    covList = vector('list', p)
    for (i in 1:p) {
      covList[[i]] = diag(n^2)
    }
  }
  if (k >= n) {k = n}
  if (k <= 0) {k = n}
  if (k == n) {
    return(eigTest(A, covList, cm, testType = testType, param.out = param.out))
  }
  if (k == 1) {
    return(schurTest(A, covList, k, cm, testType = testType, param.out = param.out, nn = nn))
  }

  matB = matrix(0, ncol = n, nrow = n)
  matC = matrix(0, ncol = n, nrow = n)
  matB[1:k, 1:k] = 1
  SB = diag(as.double(matB))
  SB = SB[which(as.double(matB) == 1),]
  matC[1:k, (k+1):n] = 1
  SC = diag(as.double(matC))
  SC = SC[which(as.double(matC) == 1),]
  SBC = rbind(SB, SC)

  SQ = SBC %*% kronecker(t(Q), t(Q))

  B = vector('list', p)
  for (i in 1:p) {
    Ai = t(Q) %*% A[[i]] %*% Q
    B[[i]] = Ai[1:k, 1:k]
  }
  V = JDTE(B)
  SV = (diag(k^2) - diag(as.vector(diag(k)))) %*% kronecker(t(V), ginv(V))

  SW = matrix(0, nrow = k*n, ncol = n*k)
  SW[1:(k^2), 1:(k^2)] = SV
  SW[(k^2+1):(k*n), (k^2+1):(n*k)] = diag(k*(n-k))

  Vlist = vector('list', p)
  for (i in 1:p) {
    Vlist[[i]] = SW %*% SQ %*% as.double(A[[i]])
    covList[[i]] = tcrossprod(SW %*% SQ, SW %*% SQ %*% covList[[i]])
  }

  if (length(testType) == 2) {
    output = c(vec.test(Vlist, covList, cm, 'chi')$pvalue, vec.test(Vlist, covList, cm, 'gam')$pvalue)
    return(output)
  } else {
    testResult = vec.test(Vlist, covList, cm, testType)
    if (param.out) {return(testResult)}
    return(testResult$pvalue)
  }

}
