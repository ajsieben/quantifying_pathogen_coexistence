# Code for transmission of 2 pathogen system that tracks within-host dynamics.
# Individual based between-host transmission model with modified Lotka-Volterra competition model for w/n host dynamics.
library("Rmpi")

# Weird workaround for using Rmpi on Teton. Load Rmpi library from the "system" and then specify the path to the directory
# with the other package libraries.

.libPaths("/project/coexistence/pathogen_coexistence/library")
library("parallel")
library("snow")
library("deSolve")

source("decomposition.R")
source("model_functions.R")
source("within_host_model.R")
source("between_host_model.R")

nProc <- 512
numSims <- 511


simFunc <- function(v) {
  runModel(v = v)
}


cl <- makeCluster(nProc - 1, type = "MPI", outfile = "")

clusterEvalQ(cl, library("deSolve"))

clusterEvalQ(cl, source("decomposition.R"))
clusterEvalQ(cl, source("model_functions.R"))
clusterEvalQ(cl, source("within_host_model.R"))
clusterEvalQ(cl, source("between_host_model.R"))

clusterApply(cl, x = 1:numSims, fun = simFunc)

# stopCluster(cl)
