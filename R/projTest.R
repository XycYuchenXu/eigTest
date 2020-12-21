#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A List of matrices to be tested. Require the length to be 2.
#' @param covList List of covariance matrices corresponding to the random matrices, default will use identity matrices
#' @param cn Constant for convergence
#' @param d Matrix dimension
#' @param param.out Logical. Whether the parameters need to be included in output
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices
#'
#' @return The p-value or a list of test information.
#' \itemize{
#' \item statistic Test statistic.
#' \item df Degrees of freedom for chi-squared distribution.
#' \item pvalue P-value.
#' }
#' @export
#'
#' @examples projTest(countryCoeff, countryCovar, cn = sqrt(112), eps = 112^(-1/3))
projTest = function(A, covList = list(), cn, eps, d = ncol(A[[1]]), param.out = FALSE){

  p = 2

  if (length(covList) == 0) {
    cov1 = diag(d^2)
    cov2 = diag(d^2)
  } else {
    cov1 = covList[[1]]
    cov2 = covList[[2]]
  }

  C = A[[1]]
  D = A[[2]]

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
  tSVD1 = truncateSVD(cov1, eps); tSVD2 = truncateSVD(cov2, eps)
  inner1 = truncateSVD(crossprod(SP.D, tSVD1$ginv %*% SP.D), eps)
  inner2 = truncateSVD(crossprod(SP.C, tSVD2$ginv %*% SP.C), eps)
  Q1 = tSVD1$ginv - tSVD1$ginv %*% SP.D %*% inner1$ginv %*% t(SP.D) %*% tSVD1$ginv
  Q2 = tSVD2$ginv - tSVD2$ginv %*% SP.C %*% inner2$ginv %*% t(SP.C) %*% tSVD2$ginv

  r = tSVD1$r + tSVD2$r - inner1$r - inner2$r
  testVal = crossprod(as.double(C), Q1 %*% as.double(C)) + crossprod(as.double(D), Q2 %*% as.double(D))
  testVal = testVal * cn^2

  if (param.out) {
    testResult = list(testVal, r, 1 - pchisq(testVal, r))
    names(testResult) = c('statistic', 'df', 'pvalue')
  } else {
    testResult = 1 - pchisq(testVal, r)
  }

  return(testResult)
}
