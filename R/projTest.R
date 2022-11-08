#' Calculate the p-value for projected MLE with specified covariance matrices
#'
#' @param A The array of two matrices to be tested with dimension \code{2}-\code{d}-\code{d}. If \code{dim(A)[1] > 2}, only the first two are used.
#' @param cn The convergence rate to normality. Assume \code{n} is the sample size, usually CLT indicates \code{cn = sqrt(n)} for consistent estimators.
#' @param cov.arr The array of covariance matrices of the two matrices with dimension \code{2}-\code{d^2}-\code{d^2}, default will use identity matrices when \code{is.null(cov.arr)}. If \code{dim(cov.arr)[1] > 2}, only the first two are used.
#' @param eps The threshold of eigenvalues when compute general inverse of covariance matrices. Required when \code{testType = 'chi'} but with default \code{cn^(-2/3)} when unsupplied.
#' @param refMat The array of reference matrices with dimension \code{2}-\code{d}-\code{d} that have the same eigenvectors. Default (when \code{is.null(refMat)}) is to use \code{A[2,,]} as reference for \code{A[1,,]} and vice versa. If \code{dim(refMat)[1] = 1}, the only reference matrix is shared.
#' @param poly.sp Logical, whether the space matrix is generated from polynomial basis. Default \code{poly.sp = TRUE} to use Legendre polynomials. Otherwise \code{poly.sp = FALSE}, the space matrix is generated with common eigenvectors.
#' @param CV Matrix of dimension \code{d}-by-\code{d} only functionable when \code{poly.sp = FALSE}, as the supplied reference common eigenvector matrix. Default (when \code{is.null(V)}) is to call \code{V = JDTE(refMat)}.
#' @param param.out Logical, whether the parameters of limiting distribution should be output or not. Default \code{param.out = FALSE} to only output P-value.
#'
#' @return A P-value when \code{param.out = FALSE} or a list of test information when \code{param.out = TRUE}.
#' \itemize{
#' \item \code{statistic}: The test statistic.
#' \item \code{df}: The degrees of freedom for chi-squared distribution.
#' \item \code{pvalue} The P-value.
#' }
#' @importFrom 'abind' abind
#' @importFrom 'Rdpack' reprompt
#' @export
#'
#' @description Only valid for nonsingular covariance matrices \code{cov.arr}.
#'
#' @references
#' \insertRef{xu2021testing}{eigTest}
#' \insertRef{andre}{eigTest}
#'
#' @examples projTest(countryCoeff, cn = sqrt(112), countryCovar)
projTest = function(A, cn, cov.arr = NULL, eps = NULL, refMat = NULL,
                    poly.sp = TRUE, CV = NULL, param.out = FALSE){

  p = 2; d = dim(A)[2]

  X = A[1,,]; Y = A[2,,]
  if (is.null(eps)) {eps = cn^(-2/3)}

  if (is.null(refMat)) {
    if (poly.sp || is.null(CV)) {refMat = A}
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
    if (is.null(CV)) {CV = JDTE(refMat)}
    CV_inv = invert(CV)
    SP.C = matrix(0, d^2, d)
    for (i in 1:d) {
      SP.C[,i] = tcrossprod(CV[,i], CV_inv[i,])
    }
    SP.D = SP.C
  }

  if (is.null(cov.arr)) {
    inner1 = crossprod(SP.D); inner2 = crossprod(SP.C)
    UDX = crossprod(SP.D, as.double(X)); UCY = crossprod(SP.C, as.double(Y))

    r = 2 * d^2 - sum(svd(inner1)$d > 1e-10) - sum(svd(inner2)$d > 1e-10)
    testVal = sum(X^2 + Y^2) -
      crossprod(UDX, crossprod(invert(inner1), UDX)) -
      crossprod(UCY, crossprod(invert(inner2), UCY))
    testVal = as.numeric(testVal * cn^2)
  } else {
    cov1.svd = truncateSVD(cov.arr[1,,], eps); cov2.svd = truncateSVD(cov.arr[2,,], eps)
    UD = crossprod(cov1.svd$rootdinv.u, SP.D); inner1 = crossprod(UD)
    UC = crossprod(cov2.svd$rootdinv.u, SP.C); inner2 = crossprod(UC)
    UX = crossprod(cov1.svd$rootdinv.u, as.double(X))
    UY = crossprod(cov2.svd$rootdinv.u, as.double(Y))
    UDX = crossprod(UD, UX); UCY = crossprod(UC, UY)

    r = cov1.svd$r + cov2.svd$r - sum(svd(inner1)$d > 1e-10) - sum(svd(inner2)$d > 1e-10)
    testVal = sum(UX^2) + sum(UY^2) -
      crossprod(UDX, crossprod(invert(inner1), UDX)) -
      crossprod(UCY, crossprod(invert(inner2), UCY))
    testVal = as.numeric(testVal * cn^2)
  }

  if (param.out) {
    testResult = list(testVal, r, 1 - pchisq(testVal, r))
    names(testResult) = c('statistic', 'df', 'pvalue')
  } else {
    testResult = 1 - pchisq(testVal, r)
  }

  return(testResult)
}
