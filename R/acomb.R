#' Combine the simulated samples, used in the foreach loop
#'
#' @param ... Any number of lists, each list is either length one that only includes the array of mean matrices, or length two with both array of mean matrices and array of covariance matrices.
#' @param along The dimension along which to bind the arrays in the list.
#'
#' @return A list of either length one that only includes the array of mean matrices, or length two with both array of mean matrices and array of covariance matrices.
#'
#' @importFrom 'abind' abind
#' @keywords internal
acomb <- function(..., along) {
  argl = list(...)
  m = argl[[1]][[1]]
  if (length(argl[[1]]) > 1) {v = argl[[1]][[2]]}
  if (length(argl) >= 2) {
    for (i in 2:length(argl)) {
      temp = argl[[i]]
      m = abind(m, temp[[1]], along = along)
      if (length(temp) > 1) {v = abind(v, temp[[2]], along = along)}
    }
  }
  if (exists('v')) {return(list(mu.bar = m, cov.bar = v))}
  return(list(mu.bar = m))
}

