### This file contains major functions for EM iterations.

### E-step.
ape1.step.spmd <- function(PARAM){
  for(i.k in 1:PARAM$K){
    logdmvnorm(PARAM, i.k)
  }

  ape1.update.expectation(PARAM)
} # End of ape1.step.spmd().

ape1.step.spmd.k <- function(PARAM, i.k, update.logL = TRUE){
  logdmvnorm(PARAM, i.k)
  ape1.update.expectation.k(PARAM, i.k, update.logL)
} # End of ape1.step.spmd.k().


### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
ape1.update.expectation <- function(PARAM, update.logL = TRUE){
  N <- nrow(X.spmd)
  K <- PARAM$K

  W.spmd <<- W.plus.y(W.spmd, PARAM$log.ETA, N, K)
  U.spmd <<- exp(W.spmd)
  Z.spmd <<- U.spmd

  tmp.id <- rowSums(W.spmd < CONTROL$exp.min) == K |
            rowSums(W.spmd > CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.spmd <- W.spmd[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.spmd) - CONTROL$exp.max / K
    } else{
      tmp.scale <- apply(tmp.spmd, 1, max) - CONTROL$exp.max / K
    }
    Z.spmd[tmp.id,] <<- exp(tmp.spmd - tmp.scale)
  }

  W.spmd.rowSums <<- rowSums(Z.spmd)
  Z.spmd <<- Z.spmd / W.spmd.rowSums

  ### For semi-supervised clustering.
  # if(SS.clustering){
  #   Z.spmd[SS.id.spmd,] <<- SS.Z.spmd
  # }

  Z.colSums <<- colSums(Z.spmd)
  Z.colSums <<- mpi.allreduce(Z.colSums, type = 2, op = "sum")
} # End of update.expectation().

ape1.update.expectation.k <- function(PARAM, i.k, update.logL = TRUE){
  N <- nrow(X.spmd)
  K <- PARAM$K

  W.spmd[, i.k] <<- W.plus.y.k(W.spmd, PARAM$log.ETA, N, K, i.k)
  U.spmd[, i.k] <<- exp(W.spmd[, i.k])
  Z.spmd <<- U.spmd

  tmp.id <- rowSums(W.spmd < CONTROL$exp.min) == K |
            rowSums(W.spmd > CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.spmd <- W.spmd[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.spmd) - CONTROL$exp.max / K
    } else{
      tmp.scale <- apply(tmp.spmd, 1, max) - CONTROL$exp.max / K
    }
    Z.spmd[tmp.id,] <<- exp(tmp.spmd - tmp.scale)
  }

  W.spmd.rowSums <<- rowSums(Z.spmd)
  Z.spmd <<- Z.spmd / W.spmd.rowSums

  ### For semi-supervised clustering.
  # if(SS.clustering){
  #   Z.spmd[SS.id.spmd,] <<- SS.Z.spmd
  # }

  Z.colSums <<- colSums(Z.spmd)
  Z.colSums <<- mpi.allreduce(Z.colSums, type = 2, op = "sum")

  if(update.logL){
    W.spmd.rowSums <<- log(W.spmd.rowSums)
    if(tmp.flag){
      W.spmd.rowSums[tmp.id] <<- W.spmd.rowSums[tmp.id] + tmp.scale
    }
  }
} # End of ap1.update.expectation().


### APECM1-step.
apecm1.step.spmd <- function(PARAM.org){
  CHECK <<- list(method = "apecm1", i.iter = 0, abs.err = Inf, rel.err = Inf,
                 convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .GlobalEnv)){
      SAVE.param <<- NULL
      SAVE.iter <<- NULL
      CLASS.iter.org <<- unlist(apply(Z.spmd, 1, which.max))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- try(apecm1.onestep.spmd(PARAM.org))
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
      CLASS.iter.new <- unlist(apply(Z.spmd, 1, which.max))
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
} # End of apecm1.step.spmd().

apecm1.onestep.spmd <- function(PARAM){
#  if(COMM.RANK == 0){
#    Rprof(filename = "apecm1.Rprof", append = TRUE)
#  }

  ### Update ETA
  PARAM <- cm.step.spmd.ETA(PARAM)
  ape1.step.spmd(PARAM)

  ### Update MU and SIGMA
  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.spmd.MU.SIGMA.k(PARAM, i.k)
    ape1.step.spmd.k(PARAM, i.k,
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
} # End of apecm1.onestep.spmd().

