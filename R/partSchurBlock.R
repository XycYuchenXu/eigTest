partSchurBlock = function(A, k, n = ncol(A[[1]]), p = length(A), iter = 5000, tol = 10^(-8)) {

  if (k <= 0) {k = n}
  if (k >= n) {return(JDTE(A, iter = iter, tol = tol))}

  score.fun = function(A, Q){
    score.fun = 0
    n = ncol(A[[1]])

    for (i in 1:p) {
      Ai = crossprod(Q, A[[i]] %*% Q)
      vi = Ai[1:k, (k+1):n]
      score.fun = score.fun + norm(vi, 'F')^2
    }
    return(score.fun)
  }

  Q = Schur(A[[1]])$Q
  Qrs = cbind(Q[,1], Q[,2])

  score.old = score.fun(A, Q)

  for (i in 1:iter) {
    for (r in 1:k) {
      for (s in (k+1):n) {
        param.2 = 0
        param.1 = 0
        for (j in 1:p) {
          Aj = A[[j]]
          param.2 = param.2 - 2 * crossprod(Q[,r], Aj %*% Q[,s]) * crossprod(Q[,s], Aj %*% Q[,r])
          for (w in 1:n) {
            if (w <= k) {
              param.2 = param.2 + (crossprod(Q[,w], Aj %*% Q[,r]))^2
              param.1 = param.1 + crossprod(Q[,w], Aj %*% Q[,r]) * crossprod(Q[,w], Aj %*% Q[,s])
            } else {
              param.2 = param.2 + (crossprod(Q[,s], Aj %*% Q[,w]))^2
              param.1 = param.1 - crossprod(Q[,r], Aj %*% Q[,w]) * crossprod(Q[,s], Aj %*% Q[,w])
            }
          }
        }
        theta = -param.1/param.2
        print(theta)
        Q[,c(s,r)] = Q[,c(s,r)] %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2)
      }
    }
    score.new = score.fun(A, Q)
    print(score.new)
    if (abs(score.new - score.old) < tol) {break()}
    score.old = score.new
  }

  return(Q)

}
