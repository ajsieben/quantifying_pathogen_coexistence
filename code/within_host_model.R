#### WITHIN-HOST MODEL ####  -----------------------------------------------------------------------------
# library("deSolve")
# library("ggplot2")

# source("./code/model_functions.R")
source("model_functions.R")

# Include event definition that represents pathogens being cleared from the host. Otherwise, an infected host will always remain infected,
# and host immunity will reach equilibrium and never decline. This condition is assessed at every time step.
eventfun <- function(t, y, parms) {
  with (as.list(y), {
    P1 <- ifelse(P1 < 0.01, 0, P1)
    P2 <- ifelse(P2 < 0.01, 0, P2)

    return(c(P1, P2, I1, I2, M1, M2, R))
  })
}

# Series of equations that describe within-host immune dynamics that is used by the ODE solver.
deriv <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dR <- theta * (1 - (R/Rmax)) - (delta1 * P1 * R) - (delta2 * P2 * R)

    # dP1 <- (alpha1 * delta1 * P1 * R) - (beta1 * I1 * P1)
    # dP2 <- (alpha2 * delta2 * P2 * R) - (beta2 * I2 * P2)
    
    dP1 <- (alpha1 * delta1 * P1 * R) - (beta1 * I1 * P1) - (beta2 * I2 * P1)  # General immunity scenario.
    dP2 <- (alpha2 * delta2 * P2 * R) - (beta2 * I2 * P2) - (beta1 * I1 * P2)  # General immunity scenario.

    # dI1 <- (epsilon * beta1 * I1 * P1) + (sigma1 * M1 * P1) + (tau * M2 * P1) - (gamma * I1)
    # dI2 <- (epsilon * beta2 * I2 * P2) + (sigma2 * M2 * P2) + (tau * M1 * P2) - (gamma * I2)
    
    dI1 <- (epsilon * beta1 * I1 * P1) + (epsilon * beta2 * I2 * P1) + (sigma1 * M1 * P1) + (sigma2 * M2 * P1) + (tau * M2 * P1) + (tau * M1 * P1) - (gamma * I1) # General immunity scenario.
    dI2 <- (epsilon * beta2 * I2 * P2) + (epsilon * beta1 * I1 * P2) + (sigma2 * M2 * P2) + (sigma1 * M1 * P2) + (tau * M1 * P2) + (tau * M2 * P2) - (gamma * I2) # General immunity scenario.
    
    # dM1 <- (rho * I1) - (mu * M1)
    # dM2 <- (rho * I2) - (mu * M2)
    
    dM1 <- (rho * I1) + (rho * I2) - (mu * M1)   # General immunity scenario.
    dM2 <- (rho * I2) + (rho * I1) - (mu * M2)   # General immunity scenario.
    
    return(list(c(P1 = dP1, P2 = dP2, I1 = dI1, I2 = dI2, M1 = dM1, M2 = dM2, R = dR)))
  })
}

# Function that resolves within-host dynamics.
withinHost <- function(temp, 
                       N, 
                       times,
                       theta,
                       Rmax,
                       delta1,
                       delta2,
                       alpha1,
                       alpha2,
                       beta1,
                       beta2,
                       epsilon,
                       sigma1,
                       sigma2,
                       tau,
                       gamma,
                       rho,
                       mu) {
  for (i in 1:N) {
    parms <- c("theta" = theta,
               "Rmax" = Rmax,
               "delta1" = delta1,
               "delta2" = delta2,
               "alpha1" = alpha1, 
               "alpha2" = alpha2, 
               "beta1" = beta1, 
               "beta2" = beta2, 
               "epsilon" = epsilon, 
               "sigma1" = sigma1,
               "sigma2" = sigma2,
               "tau" = tau, 
               "gamma" = gamma, 
               "rho" = rho, 
               "mu" = mu)
    init <- c(P1 = temp[i, "P1"], 
              P2 = temp[i, "P2"], 
              I1 = temp[i, "I1"], 
              I2 = temp[i, "I2"], 
              M1 = temp[i, "M1"], 
              M2 = temp[i, "M2"],
              R = temp[i, "R"])
    
    wnhost_results <- lsoda(y = init, times = times, func = deriv, parms = parms, events = list(func = eventfun,
                                                                                                time = times))
    temp[i,] <- wnhost_results[2, c("P1", "P2", "I1", "I2", "M1", "M2", "R")]
    # The output of lsoda is a matrix with the first column = time variable and rows are times. Only take results from the simulated step forward (i.e. second row).
  }
  
  return(temp)
  
}

##########################################################################

##### TESTING CODE ########

# N <- 2
# t <- 1000
# initP1.load <- 100    # Initial pathogen load for those infected with pathogen 1.
# initP2.load <- 100    # Initial pathogen load for those infected with pathogen 2.
# initI1.load <- 1
# initI2.load <- 1
# initM1.load <- 0
# initM2.load <- 0
# initR.load <- 100
# 
# 
# status.categ <- c("P1", "P2", "I1", "I2", "M1", "M2", "R")
# host.status <- matrix(0, nrow = N, ncol = length(status.categ),
#                       dimnames = list(seq(N), status.categ))
# host.status[, "I1"] <- initI1.load
# host.status[, "I2"] <- initI2.load
# 
# within_test <- array(host.status, dim = c(N, length(status.categ), t), dimnames = list(c(seq(N)), status.categ, seq(t)))
# within_test[,"P1",1] <- initP1.load
# within_test[,"P2",1] <- initP2.load
# within_test[,"M1",1] <- initM1.load
# within_test[,"M2",1] <- initM2.load
# within_test[,"R",1] <- initR.load
# 
# times <- seq(0, 1, by = 1)  # The number of time steps conducted in the within-host model for every time step in the between-host model.
# 
# theta <- 10
# Rmax <- 100
# delta1 <- 0.001
# delta2 <- 0.001
# alpha1 <- 100
# alpha2 <- 100
# beta1 <- 0.01
# beta2 <- 0.01
# epsilon <- 0.00001
# sigma1 <- 0.005
# sigma2 <- 0.005
# tau <- 0.0025
# gamma <- 0.05
# rho <- 0.0005
# mu <- 0.0005
# 
# 
# for (i in 1:(t-1)) {
#   temp <- within_test[,,i]
# 
#   temp <- withinHost(temp = temp,
#                      N = N,
#                      times = times,
#                      theta = theta,
#                      Rmax = Rmax,
#                      delta1 = delta1,
#                      delta2 = delta2,
#                      alpha1 = alpha1,
#                      alpha2 = alpha2,
#                      beta1 = beta1,
#                      beta2 = beta2,
#                      epsilon = epsilon,
#                      sigma1 = sigma1,
#                      sigma2 = sigma2,
#                      gamma = gamma,
#                      rho = rho,
#                      tau = tau,
#                      mu = mu)
#   within_test[,,i+1] <- temp
# }
# 
# graphWithin(t = t, results = within_test); within_test[,,t]
