
### Assign sample size to all N.worker and N.allworkers.
# assign.N.sample <- function(total.sample = 5000, N.org.worker,
#     var.table.worker = NULL){
assign.N.sample <- function(total.sample = 5000, N.org.worker){
  N.new.worker <- N.org.worker

  ### Check for SS.
  # if(!is.null(var.table.worker)){
  #   N.new.worker <- N.new.worker - nrow(var.table.worker)
  #   var.table.worker <- var.table.worker[order(var.table.worker$id),]
  # }
  my.size <- mpi.comm.size()
  N.org.allworkers <- mpi.allgather(as.integer(N.org.worker), type = 1,
                                    integer(my.size))
  N.new.allworkers <- mpi.allgather(as.integer(N.new.worker), type = 1,
                                    integer(my.size))
  if(any(N.new.allworkers <= 0)){
    stop("N.org.worker is too small.")
  }

  ### Sampling start.
  if(any(N.new.allworkers < total.sample / my.size)){
    N.sample.worker <- N.org.worker
    N.sample.allworkers <- N.org.allworkers
    ID.sample.worker <- 1:N.org.worker

    ### Check for SS.
    # if(!is.null(var.table.worker)){
    #   var.table.worker$id.sample <-
    #     which(ID.sample.worker %in% var.table.worker$id)
    # }
  } else{
    N.sample.allworkers <- floor(N.org.allworkers / sum(N.org.allworkers) *
                                 total.sample)
    remainder <- total.sample - sum(N.sample.allworkers)
    if(remainder > 0){
      N.sample.allworkers[(my.size - remainder + 1):my.size] <-
        N.sample.allworkers[(my.size - remainder + 1):my.size] + 1 
    }

    N.sample.worker <- N.sample.allworkers[mpi.comm.rank() + 1]

    ### Check for SS.
    # if(!is.null(var.table.worker)){
    #   ID.sample.worker <- sample((1:N.org.worker)[-var.table.worker$id],
    #                              N.sample.worker)
    #   ID.sample.worker <- sort(c(var.table.worker$id, ID.sample.worker))

    #   var.table.worker$id.sample <-
    #     which(ID.sample.worker %in% var.table.worker$id)
    # } else{
      ID.sample.worker <- sort(sample(1:N.org.worker, N.sample.worker))
    # }
  }

  N.sample <- sum(N.sample.allworkers)

  list(N = N.sample, N.worker = N.sample.worker,
       N.allworkers = N.sample.allworkers,
       ID.worker = ID.sample.worker)
#       SS.table.worker = var.table.worker)
} # End of assign.N.sample().

