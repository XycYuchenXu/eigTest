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
#' @param Q The orthogonal components to be tested with dimension \code{d}-\code{d}. Default (when \code{is.null(Q)}) will use \code{Q = nnPartSchur(A)} if \code{nn = TRUE} else \code{Q = expmPartSchur(A, k, warmup = warmup)}.
#' @param V The \code{k}-\code{k} eigenvector matrix to be tested for the upper diagonal block after orthogonally transformed by \code{Q}. Only applicable when \code{nn = FALSE}. Default (when \code{is.null(V)}) will use \code{V = JDTE(A)} if \code{k = d} else \code{V = JDTE(B)} where \code{B[i,,] = (t(Q) A[i,,] Q)[1:k,1:k]}.
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
#'
#' @importFrom 'Rdpack' reprompt
#' @export
#'
#' @references
#' \insertRef{xu2021testing}{eigTest}
#' \insertRef{Flury86}{eigTest}
#' \insertRef{tensor}{eigTest}
#' \insertRef{andre}{eigTest}
#'
#' @examples partialTest(countryCoeff, cn = sqrt(112), cov.arr = countryCovar, k = 2, testType = 'gam')
partialTest = function(A, cn, cov.arr = NULL, nn = FALSE, k=NULL, warmup = FALSE,
                       Q = NULL, V = NULL, testType = c('chi', 'gam'),
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

  if (k == d){
    if (is.null(V)) {V = JDTE(A)}
  } else {
    if (is.null(Q)) {Q = expmPartSchur(A, k, warmup = warmup)}
    if (is.null(V)) {
      B = array(0, c(p,k,k))
      for (i in 1:p) {
        B[i,,] = tcrossprod(crossprod(Q[,1:k], A[i,,]), t(Q[,1:k]))
      }
      V = JDTE(B)
    }
  }

  if (k == d) {
    return(eigTest(A, cn, cov.arr = cov.arr, V = V, testType = testType,
                   eps = eps, param.out = param.out))
  }

  matB = matrix(0, d, d)
  matC = matrix(0, d, d)
  matB[1:k, 1:k] = 1
  SB = diag(as.double(matB))
  SB = SB[which(as.double(matB) == 1),]
  SVB = tcrossprod(crossprod((diag(k^2) - diag(as.vector(diag(k)))),
                             kronecker(t(V), invert(V))),
                   t(SB))

  matC[(k+1):d, 1:k] = 1
  SC = diag(as.double(matC))
  SC = SC[which(as.double(matC) == 1),]

  SWQ = tcrossprod(rbind(SVB, SC), kronecker(Q, Q))
  Varr = array(0, c(p, k*d))
  for (i in 1:p) {
    Varr[i,] = crossprod(t(SWQ), as.double(A[i,,]))
    cov.arr[i, 1:(k*d), 1:(k*d)] = tcrossprod(SWQ, tcrossprod(SWQ, cov.arr[i,,]))
  }
  cov.arr = cov.arr[,1:(k*d),1:(k*d)]

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
