### This is an independent function to provide logL.
### Note that Inf, -Inf, NA, NaN is drop from the summation.

indep.logL <- function(PARAM){
  nrow <- nrow(X.worker)
  ncol <- ncol(X.worker)

  log.a <- log(2 * pi) * (-0.5 * ncol)

  ret <- matrix(0, nrow = nrow, ncol = PARAM$K)
  for(i.k in 1:PARAM$K){
    tmp.X <- W.plus.y(X.worker, -PARAM$MU[, i.k], nrow, ncol)

    tmp.S <- PARAM$SIGMA[[i.k]]
    log.b <- -0.5 * log(abs(det(tmp.S)))

    tmp.S <- solve(tmp.S)
    log.c <- -0.5 * rowSums((tmp.X %*% tmp.S) * tmp.X)

    ret[, i.k] <- log.c + log.b + log.a + PARAM$log.ETA[i.k]
  }

  ret <- rowSums(exp(ret))

  if(CONTROL$debug > 10){
    catmpi("  >>Not finite: ", sep = "")
    for(i.rank in 0:(COMM.SIZE - 1)){
      if(i.rank == COMM.RANK){
        cat(COMM.RANK, ":", sum(!is.finite(ret)), " ", sep = "")
      }
      invisible(mpi.barrier())
    }
    catmpi("\n", sep = "")
  }

  ret <- sum(log(ret[is.finite(ret)]))
  ret <- mpi.allreduce(ret, type = 2, op = "sum")
  ret
} # End of indep.logL().

