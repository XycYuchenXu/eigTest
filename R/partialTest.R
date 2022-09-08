#' Test Subset of Common Eigenvectors
#'
#' Test whether a subset of common eigenvectors are shared by the list of random matrices.
#'
#' @param A The array of matrices with dimension \code{p}-\code{d}-\code{d}, where \code{p} is the number of matrices, \code{d} is the dimension of the matrices.
#' @param cn The convergence rate(s) to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators. If \code{length(cn) < p}, all matrices share the same rate \code{cn[1]}, otherwise \code{cn = cn[1:p]}.
#' @param cov.arr The array of covariance matrices corresponding to the vectors with dimension \code{p}-\code{d^2}-\code{d^2}, default will use identity matrices when \code{is.null(cov.arr)}.
#' @param nn Logical, whether the matrices \code{A} are regarded as nonnegative transition probability matrices and the the eigenvector elements should be nonnegative.
#' @param k The number of common components to be tested. When \code{nn=TRUE}, \code{k} can be \code{NULL} or \code{k = 1}, otherwise \code{k} must be an integer within (0, \code{d}]. Will call \code{eigTest} instead if \code{k = d}.
#' @param warmup Logical, whether use \code{partSchur} for a warm-up initial value, default to \code{FALSE}.
#' @param Q The orthogonal components to be tested with dimension \code{d}-\code{d}. Default (when \code{is.null(Q)}) will use \code{Q = nnPartSchur(A)} if \code{nn} else \code{Q = expmPartSchur(A, k, warmup = warmup)}.
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
#' @examples partialTest(countryCoeff, cn = sqrt(112), cov.arr = countryCovar, k = 2, testType = 'gam')
partialTest = function(A, cn, cov.arr = NULL, nn = FALSE, k=NULL, warmup = FALSE,
                       Q = NULL, testType = c('chi', 'gam'),
                       eps = NULL, param.out = FALSE){

  d = dim(A)[2]; p = dim(A)[1]
  if (is.null(cov.arr)) {
    cov.arr = array(0, c(p, d^2, d^2))
    for (i in 1:p) {
      cov.arr[i,,] = diag(d^2)
    }
  }

  if (nn) {
    if (is.null(k) || k == 1) {
      return(schurTest(A, cn, cov.arr = cov.arr, nn = nn, Q = Q,
                       testType = testType, eps = eps, param.out = param.out))
    } else {
      print('Unavailable setup. The number of components must be 1 for nonnegative test.')
      return()
    }
  }

  if (is.null(k) || k > d || k <= 0 || k != round(k)) {
    print('Unavailable setup. The number of common components must be an integer within (0, d] for partial test.')
    return()
  }

  if (k == 1) {
    return(schurTest(A, cn, cov.arr = cov.arr, k = k, warmup = warmup, Q = Q,
                     testType = testType, eps = eps, param.out = param.out))
  }
  if (k == d) {
    return(eigTest(A, cn, cov.arr = cov.arr, testType = testType,
                   eps = eps, param.out = param.out))
  }

  if (is.null(Q)) {Q = expmPartSchur(A, k, warmup = warmup)}

  matB = matrix(0, d, d)
  matC = matrix(0, d, d)
  matB[1:k, 1:k] = 1
  SB = diag(as.double(matB))
  SB = SB[which(as.double(matB) == 1),]
  matC[1:k, (k+1):d] = 1
  SC = diag(as.double(matC))
  SC = SC[which(as.double(matC) == 1),]
  SBC = rbind(SB, SC)

  SQ = SBC %*% kronecker(t(Q), t(Q))

  B = array(0, c(p,k,k))
  for (i in 1:p) {
    Ai = t(Q) %*% A[i,,] %*% Q
    B[i,,] = Ai[1:k, 1:k]
  }
  V = JDTE(B)
  SV = (diag(k^2) - diag(as.vector(diag(k)))) %*% kronecker(t(V), ginv(V))

  SW = matrix(0, nrow = k*d, ncol = d*k)
  SW[1:(k^2), 1:(k^2)] = SV
  SW[(k^2+1):(k*d), (k^2+1):(k*d)] = diag(k*(d-k))

  Varr = array(0, c(p, k*d))
  for (i in 1:p) {
    Varr[i,] = SW %*% SQ %*% as.double(A[i,,])
    cov.arr[i, 1:(k*d), 1:(k*d)] = tcrossprod(SW %*% SQ, SW %*% SQ %*% cov.arr[i,,])
  }
  cov.arr = cov.arr[,1:(k*d),1:(k*d)]

  output = c();
  for (tp in testType) {
    temp = vec.test(Varr, cn, tp, cov.arr, eps, param.out)
    if (!is.null(temp)) {
      if (param.out) {output = list(output, temp)}
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
