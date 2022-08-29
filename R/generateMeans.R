#' Generate mean matrices for simulation samples
#'
#' @param d Size of matrices
#' @param p Number of matrices
#' @param k Number of shared components
#' @param snr Signal to noice variance ratio (SNR), can be a vector
#' @param control.g Whether the control group of samples with perturbed common
#'                  eigenvectors should be output or not
#' @param V Input of eigenvector matrix
#' @param v Input of stationary distribution
#' @param nonneg Whether the generated matrix should be nonnegative as a transition
#'               probability matrix.
#'
#' @return Array of matrices: p-by-q-by-d-by-d, where q is the number of SNRs.
#'         If \code{control.g = TRUE}, \code{q = length(snr) + 1} otherwise \code{q = 1}.
#' @export
#'
#' @importFrom 'MASS' ginv
#' @import gtools
#'
#' @examples generateMeans(5,8,3)
generateMeans = function(d, p, k = d, snr = 10, control.g = FALSE,
                         V = NULL, v = NULL, nonneg = FALSE) {

  if (k <= 0 || k > d) {k = d}
  means.groups = 1
  SNR = 0
  if (control.g) {
    means.groups = length(snr) + means.groups
    for (l in 1:(means.groups - 1)) {
      if (snr[l] <= 0) {
        snr[l] = 0
      } else {
        snr[l] = 1/sqrt(snr[l])
      }
    }
    SNR = c(SNR, snr)
  }

  mu = array(0, dim = c(p, means.groups, d, d))
  dimnames(mu) = list(paste0('Mat ID: i=', 1:p),
                      paste0('1/SNR=', SNR^2),
                      NULL, NULL)

  if (!nonneg) {
    if (is.null(V)) {
      orth = qr.Q(qr(matrix(rnorm(d^2), ncol = d)))
    } else {
      orth = qr.Q(qr(orth))
    }
    groups = sample(d, k)

    if (is.null(V) || k < d){
      V = matrix(0, ncol = d, nrow = d)
      V[1:k,] = matrix(runif(k^2, -2, 2), nrow = k) %*% orth[groups,]
    }
  } else {
    if (is.null(v)) {v = rdirichlet(1, rep(1,d))}
  }

  for (i in 1:p) {
    if (!nonneg) {
      Vi = V
      di = diag(runif(d, -1, 1))
      if (k < d) {Vi[(k+1):d,] = matrix(runif((d-k)^2, -2, 2), ncol = d-k) %*% orth[-groups,]}

    }
    for (l in 1:means.groups) {
      if (nonneg) {
        vv = v + SNR[l]*rdirichlet(1, rep(1,d))
        MarkovMat2 = matrix(0, nrow = d, ncol = d)
        for (j in 1:d) {
          for (l in 1:d) {
            if (j != l) {
              randNum = abs(rexp(1))
              MarkovMat2[j,l] = vv[l]*randNum
              MarkovMat2[l,j] = vv[j]*randNum
            }
          }
        }
        sumsMat2 = rowSums(MarkovMat2)
        total = max(sumsMat2) + abs(rnorm(1))
        for (j in 1:d) {
          MarkovMat2[j,j] =  total - sumsMat2[j]
        }
        mu[i,l,,] = MarkovMat2/total
      } else {
        Vi.perturb = Vi + matrix(rnorm(d^2), nrow = d, ncol = d) * SNR[l]
        mu[i,l,,] = ginv(Vi.perturb) %*% di %*% Vi.perturb
      }
    }
  }

  return(mu)

}
