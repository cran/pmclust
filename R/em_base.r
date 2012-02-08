### This file contains major functions for EM iterations.

### E-step.
e.step.worker <- function(PARAM, update.logL = TRUE){
  for(i.k in 1:PARAM$K){
    logdmvnorm(PARAM, i.k)
  }

  update.expectation(PARAM, update.logL = update.logL)
} # End of e.step.worker().

### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
update.expectation <- function(PARAM, update.logL = TRUE){
  N <- nrow(X.worker)
  K <- PARAM$K

  U.worker <<- W.plus.y(W.worker, PARAM$log.ETA, N, K)
  Z.worker <<- exp(U.worker)

  tmp.id <- rowSums(U.worker < CONTROL$exp.min) == K |
            rowSums(U.worker > CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.worker <- U.worker[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.worker) - CONTROL$exp.max / K
    } else{
      tmp.scale <- apply(tmp.worker, 1, max) - CONTROL$exp.max / K
    }
    Z.worker[tmp.id,] <<- exp(tmp.worker - tmp.scale)
  }

  W.worker.rowSums <<- rowSums(Z.worker)
  Z.worker <<- Z.worker / W.worker.rowSums

  ### For semi-supervised clustering.
  # if(SS.clustering){
  #   Z.worker[SS.id.worker,] <<- SS.Z.worker
  # }

  Z.colSums <<- colSums(Z.worker)
  Z.colSums <<- mpi.allreduce(Z.colSums, type = 2, op = "sum")

  if(update.logL){
    W.worker.rowSums <<- log(W.worker.rowSums)
    if(tmp.flag){
      W.worker.rowSums[tmp.id] <<- W.worker.rowSums[tmp.id] + tmp.scale
    }
  }
} # End of update.expectation().


### M-step.
m.step.worker <- function(PARAM){
  ### MLE For ETA
  PARAM$ETA <- Z.colSums / sum(Z.colSums)
  PARAM$log.ETA <- log(PARAM$ETA)

  p <- PARAM$p
  for(i.k in 1:PARAM$K){
    ### MLE for MU
    tmp.MU <- colSums(X.worker * Z.worker[, i.k]) / Z.colSums[i.k]
    PARAM$MU[, i.k] <- mpi.allreduce(tmp.MU, type = 2, op = "sum")

    ### MLE for SIGMA
    if(PARAM$U.check[[i.k]]){
      B <- W.plus.y(X.worker, -PARAM$MU[, i.k],
                    nrow(X.worker), ncol(X.worker)) *
           sqrt(Z.worker[, i.k] / Z.colSums[i.k])
      tmp.SIGMA <- crossprod(B)
      tmp.SIGMA <- mpi.allreduce(tmp.SIGMA, type = 2, op = "sum") 
      dim(tmp.SIGMA) <- c(p, p)

      tmp.U <- decompsigma(tmp.SIGMA)
      PARAM$U.check[[i.k]] <- tmp.U$check
      if(tmp.U$check){
        PARAM$U[[i.k]] <- tmp.U$value
        PARAM$SIGMA[[i.k]] <- tmp.SIGMA
      }
    } else{
      if(CONTROL$debug > 2){
        catmpi("  SIGMA[[", i.k, "]] is fixed.\n", sep = "")
      }
    }
  }

  PARAM
} # End of m.step.worker().


### log likelihood.
logL.step <- function(){
  tmp.logL <- sum(W.worker.rowSums)
  mpi.allreduce(tmp.logL, type = 2, op = "sum")
} # End of logL.step().


### Check log likelihood convergence.
check.em.convergence <- function(PARAM.org, PARAM.new, i.iter){
  abs.err <- PARAM.new$logL - PARAM.org$logL
  rel.err <- abs.err / abs(PARAM.org$logL)
  convergence <- 0

  if(abs.err < 0){
    convergence <- 4
  } else if(any(PARAM.new$ETA < PARAM.new$min.N.CLASS / PARAM.new$N)){
    convergence <- 3
  } else if(i.iter > CONTROL$max.iter){
    convergence <- 2
  } else if(rel.err < CONTROL$rel.err){
    convergence <- 1
  }

  if(CONTROL$debug > 1){
    catmpi("  check.em.convergence:",
           " abs: ", abs.err,
           ", rel: ", rel.err,
           ", conv: ", convergence, "\n", sep = "")
  }

  list(method = CHECK$method,
       iter = i.iter, abs.err = abs.err, rel.err = rel.err,
       convergence = convergence)
} # End of check.em.convergence().


### EM-step.
em.step.worker <- function(PARAM.org){
  CHECK <<- list(method = "em", i.iter = 0, abs.err = Inf, rel.err = Inf,
                 convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .GlobalEnv)){
      SAVE.param <<- NULL
      SAVE.iter <<- NULL
      CLASS.iter.org <<- unlist(apply(Z.worker, 1, which.max))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- try(em.onestep.worker(PARAM.org))
    if(class(PARAM.new) == "try-error"){
      catmpi("Results of previous iterations are returned.\n")
      CHECK$convergence <<- 99
      PARAM.new <- PARAM.org
      break
    }

    CHECK <<- check.em.convergence(PARAM.org, PARAM.new, i.iter)
    if(CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      SAVE.param <<- c(SAVE.param, PARAM.new)
      CLASS.iter.new <- unlist(apply(Z.worker, 1, which.max))
      tmp <- as.double(sum(CLASS.iter.new != CLASS.iter.org))
      tmp <- mpi.allreduce(tmp, type = 2, op = "sum")
      tmp.all <- c(tmp / PARAM$N, PARAM.new$logL,
                   PARAM.new$logL - PARAM.org$logL,
                   (PARAM.new$logL - PARAM.org$logL) / PARAM.org$logL)
      SAVE.iter <<- rbind(SAVE.iter, c(tmp, tmp.all, tmp.time))
      CLASS.iter.org <<- CLASS.iter.new
    }

    PARAM.org <- PARAM.new
    i.iter <- i.iter + 1
  }

  PARAM.new
} # End of em.step.worker().

em.onestep.worker <- function(PARAM){
#  if(COMM.RANK == 0){
#    Rprof(filename = "em.Rprof", append = TRUE)
#  }

  PARAM <- m.step.worker(PARAM)
  e.step.worker(PARAM)

#  if(COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step()

  if(CONTROL$debug > 0){
    catmpi(">>em.onestep: ", format(Sys.time(), "%H:%M:%S"),
           ", iter: ", CHECK$iter, ", logL: ",
                       sprintf("%-30.15f", PARAM$logL), "\n", sep = "")
    if(CONTROL$debug > 4){
      logL <- indep.logL(PARAM)
      catmpi("  >>indep.logL: ", sprintf("%-30.15f", logL), "\n", sep = "")
    }
    if(CONTROL$debug > 20){
      mb.print(PARAM, CHECK)
    }
  }

  PARAM
} # End of em.onestep.worker().


### Obtain classifications.
em.update.class.worker <- function(){
  CLASS.worker <<- unlist(apply(Z.worker, 1, which.max))
} # End of em.update.class.worker().

