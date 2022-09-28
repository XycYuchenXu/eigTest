#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A The array of two matrices to be tested with dimension \code{2}-\code{d}-\code{d}. If \code{dim(A)[1] > 2}, only the first two are used.
#' @param cn The convergence rate to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators.
#' @param cov.arr The array of covariance matrices of the two matrices with dimension \code{2}-\code{d^2}-\code{d^2}, default will use identity matrices when \code{is.null(cov.arr)}. If \code{dim(cov.arr)[1] > 2}, only the first two are used.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.
#' @param refMat The array of reference matrices with dimension \code{2}-\code{d}-\code{d} that have the same eigenvectors. Default (when \code{is.null(refMat)}) is to use \code{A[2,,]} as reference for \code{A[1,,]} and vice versa. If \code{dim(refMat)[1] = 1}, the only reference matrix is shared.
#' @param poly.sp Logical, whether the space matrix is generated from polynomial basis. Default \code{poly.sp = TRUE} to use Legendre polynomials. Otherwise \code{poly.sp = FALSE}, the space matrix is generated with estimated common eigenvectors from \code{JDTE(refMat)}.
#' @param param.out Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output P-value.
#'
#' @return A P-value when \code{param.out = FALSE} or a list of test information when \code{param.out = TRUE}.
#' \itemize{
#' \item \code{statistic}: The test statistic.
#' \item \code{df}: The degrees of freedom for chi-squared distribution.
#' \item \code{pvalue} The P-value.
#' }
#' @importFrom 'MASS' ginv
#' @importFrom 'abind' abind
#' @export
#'
#' @description Only valid for nonsingular covariance matrices \code{cov.arr}.
#'
#' @examples projTest(countryCoeff, cn = sqrt(112), countryCovar)
projTest = function(A, cn, cov.arr = NULL, eps = NULL, refMat = NULL,
                    poly.sp = TRUE, param.out = FALSE){

  p = 2; d = dim(A)[2]

  if (is.null(cov.arr)) {
    cov1 = diag(d^2)
    cov2 = diag(d^2)
  } else {
    cov1 = cov.arr[1,,]
    cov2 = cov.arr[2,,]
  }

  X = A[1,,]; Y = A[2,,]
  if (is.null(eps)) {eps = cn^(-2/3)}

  if (is.null(refMat)) {
    refMat = A
  } else if (dim(refMat)[1] == 1) {
    refMat = abind(refMat[1,,], refMat[1,,], along = 0)
  }

  if (poly.sp) {
    C = t(refMat[1,,]); D = t(refMat[2,,])
    C0 = diag(d); D0 = diag(d); C1 = t(C); D1 = t(D)
    C1 = C1 / norm(C1); D1 = D1 / norm(D1)
    SP.C = cbind(as.vector(C0), as.vector(C1))
    SP.D = cbind(as.vector(D0), as.vector(D1))
    if (d > 2){
      for (i in 1:(d-2)) {
        C.temp = ((2*i + 1) * crossprod(C, C1) - i * C0) / (i + 1)
        D.temp = ((2*i + 1) * crossprod(D, D1) - i * D0) / (i + 1)
        #C.temp = C.temp / norm(C.temp); D.temp = D.temp / norm(D.temp)
        SP.C = cbind(SP.C, as.vector(C.temp / norm(C.temp)))
        SP.D = cbind(SP.D, as.vector(D.temp / norm(D.temp)))
        C0 = C1; D0 = D1; C1 = C.temp; D1 = D.temp
      }
    }
  } else {
    CV = JDTE(refMat)
    CV_inv = ginv(CV)
    SP.C = matrix(0, d^2, d)
    for (i in 1:d) {
      SP.C[,i] = tcrossprod(CV[,i], CV_inv[i,])
    }
    SP.D = SP.C
  }

  cov1.svd = truncateSVD(cov1, eps); cov2.svd = truncateSVD(cov2, eps)
  inner1 = crossprod(SP.D, crossprod(cov1.svd$ginv, SP.D))
  inner2 = crossprod(SP.C, crossprod(cov2.svd$ginv, SP.C))
  ginvSPD = crossprod(cov1.svd$ginv, SP.D); ginvSPC = crossprod(cov2.svd$ginv, SP.C)
  Q1 = cov1.svd$ginv - tcrossprod(ginvSPD, tcrossprod(ginvSPD, ginv(inner1)))
  Q2 = cov2.svd$ginv - tcrossprod(ginvSPC, tcrossprod(ginvSPC, ginv(inner2)))

  r = cov1.svd$r + cov2.svd$r - sum(svd(inner1)$d > 1e-10) - sum(svd(inner2)$d > 1e-10)
  testVal = crossprod(as.double(X), crossprod(Q1, as.double(X))) + crossprod(as.double(Y), crossprod(Q2, as.double(Y)))
  testVal = testVal * cn^2

  if (param.out) {
    testResult = list(testVal, r, as.numeric(1 - pchisq(testVal, r)))
    names(testResult) = c('statistic', 'df', 'pvalue')
  } else {
    testResult = as.numeric(1 - pchisq(testVal, r))
  }

  return(testResult)
}
