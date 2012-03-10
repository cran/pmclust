### This file contains major functions for EM iterations.

### E-step.
ape.step.spmd.k <- function(PARAM, i.k, update.logL = TRUE){
  logdmvnorm(PARAM, i.k)
  update.expectation(PARAM, update.logL = update.logL)
} # End of ape.step.spmd.k().

### CM-step
cm.step.spmd.ETA.MU.SIGMA.k <- function(PARAM, i.k){
  ### MLE For ETA
  PARAM$ETA <- Z.colSums / sum(Z.colSums)
  PARAM$log.ETA <- log(PARAM$ETA)

  ### MLE for MU and SIGMA
  PARAM <- cm.step.spmd.MU.SIGMA.k(PARAM, i.k)

  PARAM
} # End of cm.step.spmd.ETA.MU.SIGMA.k().

cm.step.spmd.MU.SIGMA.k <- function(PARAM, i.k){
  ### MLE for MU
  tmp.MU <- colSums(X.spmd * Z.spmd[, i.k]) / Z.colSums[i.k]
  PARAM$MU[, i.k] <- mpi.allreduce(tmp.MU, type = 2, op = "sum")

  ### MLE for SIGMA
  if(PARAM$U.check[[i.k]]){
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow(X.spmd), ncol(X.spmd)) *
         sqrt(Z.spmd[, i.k] / Z.colSums[i.k])
    tmp.SIGMA <- crossprod(B)
    tmp.SIGMA <- mpi.allreduce(tmp.SIGMA, type = 2, op = "sum") 
    dim(tmp.SIGMA) <- c(PARAM$p, PARAM$p)

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

  PARAM
} # End of cm.step.spmd.MU.SIGMA.k().


### APECM-step.
apecm2.step.spmd <- function(PARAM.org){
  CHECK <<- list(method = "apecm2", i.iter = 0, abs.err = Inf, rel.err = Inf,
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

    PARAM.new <- try(apecm2.onestep.spmd(PARAM.org))
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
} # End of apecm2.step.spmd().

apecm2.onestep.spmd <- function(PARAM){
#  if(COMM.RANK == 0){
#    Rprof(filename = "apecm2.Rprof", append = TRUE)
#  }

  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.spmd.ETA.MU.SIGMA.k(PARAM, i.k)
    ape.step.spmd.k(PARAM, i.k,
                      update.logL = ifelse(i.k == PARAM$K, TRUE, FALSE))
  }

#  if(COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step()

  if(CONTROL$debug > 0){
    catmpi(">>apecm2.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of apecm2.onestep.spmd().

