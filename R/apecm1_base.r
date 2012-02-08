### This file contains major functions for EM iterations.

### E-step.
ape1.step.worker <- function(PARAM){
  for(i.k in 1:PARAM$K){
    logdmvnorm(PARAM, i.k)
  }

  ape1.update.expectation(PARAM)
} # End of ape1.step.worker().

ape1.step.worker.k <- function(PARAM, i.k, update.logL = TRUE){
  logdmvnorm(PARAM, i.k)
  ape1.update.expectation.k(PARAM, i.k, update.logL)
} # End of ape1.step.worker.k().


### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
ape1.update.expectation <- function(PARAM, update.logL = TRUE){
  N <- nrow(X.worker)
  K <- PARAM$K

  W.worker <<- W.plus.y(W.worker, PARAM$log.ETA, N, K)
  U.worker <<- exp(W.worker)
  Z.worker <<- U.worker

  tmp.id <- rowSums(W.worker < CONTROL$exp.min) == K |
            rowSums(W.worker > CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.worker <- W.worker[tmp.id,]

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
} # End of update.expectation().

ape1.update.expectation.k <- function(PARAM, i.k, update.logL = TRUE){
  N <- nrow(X.worker)
  K <- PARAM$K

  W.worker[, i.k] <<- W.plus.y.k(W.worker, PARAM$log.ETA, N, K, i.k)
  U.worker[, i.k] <<- exp(W.worker[, i.k])
  Z.worker <<- U.worker

  tmp.id <- rowSums(W.worker < CONTROL$exp.min) == K |
            rowSums(W.worker > CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.worker <- W.worker[tmp.id,]

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
} # End of ap1.update.expectation().


### APECM1-step.
apecm1.step.worker <- function(PARAM.org){
  CHECK <<- list(method = "apecm1", i.iter = 0, abs.err = Inf, rel.err = Inf,
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

    PARAM.new <- try(apecm1.onestep.worker(PARAM.org))
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
} # End of apecm1.step.worker().

apecm1.onestep.worker <- function(PARAM){
#  if(COMM.RANK == 0){
#    Rprof(filename = "apecm1.Rprof", append = TRUE)
#  }

  ### Update ETA
  PARAM <- cm.step.worker.ETA(PARAM)
  ape1.step.worker(PARAM)

  ### Update MU and SIGMA
  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.worker.MU.SIGMA.k(PARAM, i.k)
    ape1.step.worker.k(PARAM, i.k,
                       update.logL = ifelse(i.k == PARAM$K, TRUE, FALSE))
  }

#  if(COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step()

  if(CONTROL$debug > 0){
    catmpi(">>apecm1.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of apecm1.onestep.worker().

