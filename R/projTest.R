#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A The array of two matrices to be tested with dimension \code{2}-\code{d}-\code{d}. If \code{dim(A)[1] > 2}, only the first two are used.
#' @param cn The convergence rate to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators.
#' @param cov.arr The array of covariance matrices of the two matrices with dimension \code{2}-\code{d^2}-\code{d^2}, default will use identity matrices when \code{is.null(cov.arr)}. If \code{dim(cov.arr)[1] > 2}, only the first two are used.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.
#' @param refMat The array of reference matrices with dimension \code{2}-\code{d}-\code{d} that have the same eigenvectors. Default (when \code{is.null(refMat)}) is to use \code{A[2,,]} as reference for \code{A[1,,]} and vice versa. If \code{dim(refMat)[1] = 1}, the only reference matrix is shared.
#' @param param.out Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output P-value.
#'
#' @return A P-value when \code{param.out = FALSE} or a list of test information when \code{param.out = TRUE}.
#' \itemize{
#' \item \code{statistic}: The test statistic.
#' \item \code{df}: The degrees of freedom for chi-squared distribution.
#' \item \code{pvalue} The P-value.
#' }
#' @importFrom 'MASS' ginv
#' @export
#'
#' @description Only valid for nonsingular covariance matrices \code{cov.arr}.
#'
#' @examples projTest(countryCoeff, cn = sqrt(112), countryCovar)
projTest = function(A, cn, cov.arr = NULL, eps = NULL,
                    refMat = NULL, param.out = FALSE){

  p = 2; d = dim(A)[2]

  if (is.null(cov.arr)) {
    cov1 = diag(d^2)
    cov2 = diag(d^2)
  } else {
    cov1 = cov.arr[1,,]
    cov2 = cov.arr[2,,]
  }

  if (is.null(refMat)) {
    C = A[1,,]; D = A[2,,]
  } else if (dim(refMat)[1] == 1) {
    C = refMat[1,,]; D = refMat[1,,]
  } else {
    C = refMat[1,,]; D = refMat[2,,]
  }

  if (is.null(eps)) {eps = cn^(-2/3)}

  X = A[1,,]; Y = A[2,,]

  C.temp = diag(d)
  D.temp = diag(d)
  SP.C = as.vector(C.temp)/norm(C.temp, 'F')
  SP.D = as.vector(D.temp)/norm(D.temp, 'F')
  for (i in 1:(d-1)) {
    C.temp = C.temp %*% C
    D.temp = D.temp %*% D
    C.temp = C.temp/norm(C.temp, 'F')
    D.temp = D.temp/norm(D.temp, 'F')
    SP.C = cbind(SP.C, as.vector(C.temp))
    SP.D = cbind(SP.D, as.vector(D.temp))
  }
  cov1.svd = truncateSVD(cov1, eps); cov2.svd = truncateSVD(cov2, eps)
  inner1 = ginv(crossprod(SP.D, cov1.svd$ginv %*% SP.D)); inner2 = ginv(crossprod(SP.C, cov2.svd$ginv %*% SP.C))
  Q1 = cov1.svd$ginv - cov1.svd$ginv %*% SP.D %*% inner1 %*% t(SP.D) %*% cov1.svd$ginv
  Q2 = cov2.svd$ginv - cov2.svd$ginv %*% SP.C %*% inner2 %*% t(SP.C) %*% cov2.svd$ginv

  r = cov1.svd$r + cov2.svd$r
  testVal = crossprod(as.double(X), Q1 %*% as.double(X)) + crossprod(as.double(Y), Q2 %*% as.double(Y))
  testVal = testVal * cn^2

  if (param.out) {
    testResult = list(testVal, r, 1 - pchisq(testVal, r))
    names(testResult) = c('statistic', 'df', 'pvalue')
  } else {
    testResult = 1 - pchisq(testVal, r)
  }

  return(testResult)
}
