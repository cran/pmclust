### This file contains functions for log density of MVN.
### These will majorly update W.spmd.

logdmvnorm <- function(PARAM, i.k){
#  for(i.k in 1:PARAM$K){
#    U <- chol(PARAM$SIGMA[[i.k]])
    U <- PARAM$U[[i.k]]
    logdet <- sum(log(abs(diag(U)))) * 2
#    B <- t.X.spmd - PARAM$MU[, i.k]
#    A <- backsolve(U, B, upper.tri = TRUE, transpose = TRUE)
#    distval <- colSums(A * A)
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow(X.spmd), ncol(X.spmd))
    B <- B %*% backsolve(U, diag(PARAM$p))
    distval <- rowSums(B * B)
    W.spmd[, i.k] <<- -(p.times.logtwopi + logdet + distval) * 0.5
#  }
} # End of logdmvnorm().

