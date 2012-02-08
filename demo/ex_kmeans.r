### Setup mpi environment.
library(Rmpi)
library(pmclust)
invisible(mpi.comm.dup(0, 1))

### Generate an example data.
N.allworkers <- rep(5000, mpi.comm.size())
N.worker <- 5000
N.K.worker <- c(2000, 3000)
N <- 5000 * mpi.comm.size()
p <- 2
K <- 2
seed <- 123 + mpi.comm.rank()
data.worker <- generate.basic.worker(N.allworkers, N.worker, N.K.worker,
                                     N, p, K, seed)
X.worker <- data.worker$X.worker

### Run clustering.
PARAM.org <- set.global(K = K)
PARAM.org <- initial.center.worker(PARAM.org)
PARAM.new <- kmeans.step.worker(PARAM.org)
kmeans.update.class.worker()
mb.print(PARAM.new, CHECK)

### Get results.
N.CLASS <- get.N.CLASS(K)
catmpi("# of class:", N.CLASS, "\n")

### Print run time and quit Rmpi.
# printmpi(proc.time())
mpi.quit()
