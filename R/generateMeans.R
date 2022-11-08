#' Generate mean matrices for simulation samples
#'
#' @param d The dimension of matrices.
#' @param p The number of matrices.
#' @param k The number of common Schur components. Must be an integer within (0, \code{d}), otherwise set \code{k = d}.
#' @param snr The positive signal to noise variance ratio (SNR), can be a vector. Default \code{is.null(snr) = TRUE} will only generate means with non-perturbed eigenvectors.
#' @param V The input of eigenvector matrix with dimension \code{d}-\code{d}, needed when \code{nn = FALSE}. Default will use random sampling when \code{is.null(V) = TRUE}.
#' @param v The input of stationary distribution with length \code{d}, needed when \code{nn = TRUE}. Default will use Dirichlet random sampling when \code{is.null(v) = TRUE}.
#' @param nn Logical, whether the generated matrices should be nonnegative as transition probability matrices.
#'
#' @return An array of mean matrices: \code{p}-by-\code{q}-by-\code{d}-by-\code{d}, where \code{q} is the number of SNRs: \code{q = length(snr) + 1}.
#' @export
#'
#' @import gtools
#'
#' @examples generateMeans(5,8,3)
generateMeans = function(d, p, k = d, snr = NULL, V = NULL, v = NULL, nn = FALSE) {

  if (k <= 0 || k > d || k != round(k)) {k = d}
  means.groups = 1
  SNR = 0
  if (is.null(snr) == FALSE) {
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

  if (!nn) {
    if (is.null(V)) {
      orth = qr.Q(qr(matrix(rnorm(d^2), ncol = d)))
    } else {
      orth = qr.Q(qr(V))
    }
    groups = sample(d, k)

    if (is.null(V) || k < d){
      V = matrix(0, ncol = d, nrow = d)
      coefM = matrix(runif(k^2, -1, 1), nrow = k); coefM[!upper.tri(coefM)] = 0
      V[,1:k] = tcrossprod(orth[,groups], diag(k) + coefM)
    }
  } else {
    if (is.null(v)) {v = rdirichlet(1, rep(1,d))}
  }

  for (i in 1:p) {
    if (!nn) {
      Vi = V
      di = diag(runif(d, 0.5, d) * sample(c(-1,1), d, replace = T))
      if (k < d) {
        coefM = matrix(runif((d-k)^2, -2, 2), ncol = d-k)
        Vi[,(k+1):d] = tcrossprod(orth[,-groups], coefM)
      }

    }
    for (l in 1:means.groups) {
      if (nn) {
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
        mu[i,l,,] = t(MarkovMat2/total)
      } else {
        Vi.perturb = Vi + matrix(rnorm(d^2), nrow = d, ncol = d) * SNR[l]
        mu[i,l,,] = tcrossprod(Vi.perturb, crossprod(invert(Vi.perturb), di))
      }
    }
  }

  return(mu)

}
