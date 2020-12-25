#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A List of matrices to be tested. Require the length to be 2.
#' @param covList List of covariance matrices corresponding to the random matrices, default will use identity matrices
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
#' @description Only valid for nonsingular covariance matrices \code{covList}.
#'
#' @examples projTest(countryCoeff, countryCovar, cn = sqrt(112), eps = 112^(-1/3))
projTest = function(A, covList = list(), refMat = A, cn, eps, d = ncol(A[[1]]), param.out = FALSE){

  p = 2

  if (length(covList) == 0) {
    cov1 = diag(d^2)
    cov2 = diag(d^2)
  } else {
    cov1 = covList[[1]]
    cov2 = covList[[2]]
  }

  if (length(refMat) == 1) {
    refMat = list(refMat[[1]], refMat[[1]])
  }

  X = A[[1]]; Y = A[[2]]

  C = refMat[[1]]; D = refMat[[2]]

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
  cov1.r = ginv(cov1); cov2.r = ginv(cov2)
  inner1 = ginv(crossprod(SP.D, cov1.r %*% SP.D)); inner2 = ginv(crossprod(SP.C, cov2.r %*% SP.C))
  Q1 = cov1.r - cov1.r %*% SP.D %*% inner1 %*% t(SP.D) %*% cov1.r
  Q2 = cov2.r - cov2.r %*% SP.C %*% inner2 %*% t(SP.C) %*% cov2.r

  r = 2*d*(d-1)
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
