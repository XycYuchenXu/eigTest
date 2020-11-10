#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A List of matrices to be tested. Require the length to be 2.
#' @param covList List of covariance matrices corresponding to the random matrices
#' @param cm Constant for convergence
#' @param n Matrix dimension
#' @param param.out Logical. Whether the parameters need to be included in output
#'
#' @return The p-value or a list of test information.
#' \itemize{
#' \item statistic Test statistic.
#' \item df Degrees of freedom for chi-squared distribution.
#' \item pvalue P-value.
#' }
#' @export
#'
#' @importFrom MASS ginv
#'
#' @examples projTest(countryCoeff, countryCovar, cm = 102)
projTest = function(A, covList = list(), cm, n = ncol(A[[1]]), param.out = FALSE){

  p = 2

  if (length(covList) == 0) {
    cov1 = diag(n^2)
    cov2 = diag(n^2)
  } else {
    cov1 = covList[[1]]
    cov2 = covList[[2]]
  }

  C = A[[1]]
  D = A[[2]]

  C.temp = diag(n)
  D.temp = diag(n)
  SP.C = as.vector(C.temp)/norm(C.temp, 'F')
  SP.D = as.vector(D.temp)/norm(D.temp, 'F')
  for (i in 1:(n-1)) {
    C.temp = C.temp %*% C
    D.temp = D.temp %*% D
    C.temp = C.temp/norm(C.temp, 'F')
    D.temp = D.temp/norm(D.temp, 'F')
    SP.C = cbind(SP.C, as.vector(C.temp))
    SP.D = cbind(SP.D, as.vector(D.temp))
  }
  Q1 = ginv(cov1) - ginv(cov1) %*% SP.D %*% ginv(crossprod(SP.D, ginv(cov1) %*% SP.D)) %*% t(SP.D) %*% ginv(cov1)
  Q2 = ginv(cov2) - ginv(cov2) %*% SP.C %*% ginv(crossprod(SP.C, ginv(cov2) %*% SP.C)) %*% t(SP.C) %*% ginv(cov2)

  r = 2*n^2-2*n
  testVal = crossprod(as.double(C), Q1 %*% as.double(C)) + crossprod(as.double(D), Q2 %*% as.double(D))
  testVal = testVal * cm^2

  if (param.out) {
    testResult = list(testVal, r[1], 1 - pchisq(testVal, r))
    names(testResult) = c('statistic', 'df', 'pvalue')
  } else {
    testResult = 1 - pchisq(testVal, r)
  }

  return(testResult)
}
