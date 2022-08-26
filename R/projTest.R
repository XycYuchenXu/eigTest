#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A Array of matrices to be tested. Require the length to be 2.
#' @param cov.arr List of covariance matrices corresponding to the random matrices, default will use identity matrices
#' @param cn Constant for convergence
#' @param d Matrix dimension
#' @param param.out Logical. Whether the parameters need to be included in output
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices
#' @param refMat The list of reference matrices with the same eigenvectors, default is to use the estimated \code{A}
#'
#' @return The p-value or a list of test information.
#' \itemize{
#' \item statistic Test statistic.
#' \item df Degrees of freedom for chi-squared distribution.
#' \item pvalue P-value.
#' }
#' @importFrom MASS ginv
#' @export
#'
#' @description Only valid for nonsingular covariance matrices \code{cov.arr}.
#'
#' @examples projTest(countryCoeff, countryCovar, cn = sqrt(112), eps = 112^(-1/3))
projTest = function(A, cov.arr = NULL, refMat = A, cn, eps, d = dim(A)[2], param.out = FALSE){

  p = 2

  if (is.null(cov.arr)) {
    cov1 = diag(d^2)
    cov2 = diag(d^2)
  } else {
    cov1 = cov.arr[1,,]
    cov2 = cov.arr[2,,]
  }

  if (dim(refMat)[1] == 1) {
    C = refMat[1,,]; D = refMat[1,,]
  } else {
    C = refMat[1,,]; D = refMat[2,,]
  }

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
