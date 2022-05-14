#### INITIALIZING FUNCTIONS ####

# Function to draw parameters for between-host and within-host models.
generateParams <- function() {  
  # Input the parameters to generate depending on the scenario you wish to model (i.e. neutral, comp/col, sensitivity analysis).
  
  #----------------------------
  
  ##### NEUTRAL PATHOGENS #####
  N <- 1000   # See decomposition.R for descriptions of each parameter.          
  t <- 2000           
  
  sd.var <- 0
  
  # Between-host model parameters ##
  kappa <- 2
  bottle <- 0.05
  eta <- 0.001
  
  lo.n <- 0
  med.n <- 5000
  hi.n <- 10000
  
  lo.k.p1 <- 0
  lo.k.p2 <- 0
  med.k.p1 <- 0.025
  med.k.p2 <- 0.025
  hi.k.p1 <- 0.05
  hi.k.p2 <- 0.05
  
  ## Within-host model parameters ##
  times <- seq(0, 1, by = 1)
  
  theta <- 10
  Rmax <- 100
  delta1 <- 0.001
  delta2 <- 0.001
  alpha1 <- 100
  alpha2 <- 100
  beta1 <- 0.01
  beta2 <- 0.01
  epsilon <- 0.00001
  sigma1 <- 0.005
  sigma2 <- 0.005
  tau <- 0.0025
  gamma <- 0.05
  rho <- 0.0005
  mu <- 0.0005
  
  # Coexistence decomposition parameters.
  equil.times <- c(1900:2000)
  spat.equil.times <- 10
  
  #----------------------------
  
  parameter.list <- list("N" = N,
                         "t" = t,
                         "sd.var" = sd.var,
                         "kappa" = kappa,
                         "bottle" = bottle,
                         "eta" = eta,
                         "lo.n" = lo.n,
                         "med.n" = med.n,
                         "hi.n" = hi.n,
                         "lo.k.p1" = lo.k.p1,
                         "lo.k.p2" = lo.k.p2,
                         "med.k.p1" = med.k.p1,
                         "med.k.p2" = med.k.p2,
                         "hi.k.p1" = hi.k.p1,
                         "hi.k.p2" = hi.k.p2,
                         "times" = times,
                         "theta" = theta,
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
                         "gamma" = gamma,
                         "rho" = rho,
                         "tau" = tau,
                         "mu" = mu,
                         "equil.times" = equil.times,
                         "spat.equil.times" = spat.equil.times)
  
  return(parameter.list)
                         
}


# Function that takes initial conditions and enters them into object holding host states.
initialCondit <- function(results, 
                          initP1, 
                          initP1.load, 
                          initP2, 
                          initP2.load) { 
  results[1:initP1,"P1"] <- ifelse(initP1 != 0, initP1.load, 0)                        # Set hosts with initial P1 pathogen loads.
  results[(initP1 + 1):(initP1 + initP2),"P2"] <- ifelse(initP2 != 0, initP2.load, 0)  # Set following hosts with initial P2 pathogen loads.
  
  return(results)
  
}


# Function to create dummy dataset and generate logit function to determine probability of between-host transmission.
calcLogit <- function(lo.k, lo.n, med.k, med.n, hi.k, hi.n) {
  
  # Construct dummy dataset with two variables: pathogen load (continuous; XX.n) and whether they transmit (binary; XX.k).
  # Dummy dataset has three levels of pathogen load: low, medium, and high.
  a.k <- rbinom(n = 10000, size = 1, prob = lo.k)   # Number of successful transmitions for low pathogen load.
  a.n <- rep(x = lo.n, times = 10000)               # Low pathogen loads for corresponding transmition events 'a.k'.
  
  b.k <- rbinom(n = 10000, size = 1, prob = med.k)  # Number of successful transmitions for medium pathogen load.
  b.n <- rep(x = med.n, times = 10000)              # Medium pathogen loads for corresponding transmition events 'b.k'.
  
  c.k <- rbinom(10000, size = 1, prob = hi.k)       # Number of successful transmitions for high pathogen load.
  c.n <- rep(x = hi.n, times = 10000)               # high pathogen loads for corresponding transmition events 'c.k'.
  
  temp1 <- append(a.k, values = c(b.k, c.k))
  temp2 <- append(a.n, values = c(b.n, c.n))
  
  temp3 <- data.frame("trans" = temp1, "load" = temp2)
  
  # Conduct logistic regression on daummy dataset and pull intercept and coefficient from regression output.
  logit <- glm(trans ~ load, family = binomial(link = "logit"), data = temp3)
  coef.logit <- coef(logit)
  
  return(coef.logit)
  
}


# Function to convert log odds to probability for between-host transmission.
logitToProb <- function(coef.logit, path.load) {
  odds <- exp(coef.logit[1] + coef.logit[2] * path.load)
  prob <- odds / (1 + odds)
  return(prob)
}
    

#############################################################################################

#### DECOMPOSITION FUNCTIONS ####

# Function to determine the spatial equilibrium of invading pathogen and calculate host states at t0 and t1 with invading pathogen abundances at spatial equilibrium.
calculateLDGR <- function(resident.equilibrium,
                         host.status,
                         status.categ,
                         invader, 
                         equil.times, 
                         spat.equil.times, 
                         iabun,
                         invade.abundance,
                         t,
                         host.contacts.eq,
                         coef.logit,
                         host.trans.probs.eq,
                         N,
                         times,
                         bottle, 
                         eta,
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
                         gamma, 
                         rho, 
                         tau, 
                         mu) {
  
  # Create objects that will hold spatial invader equilibrium calculation results and host states at t0 and t1.
  space.equil <- array(NA, dim = c(N, 1, spat.equil.times), dimnames = list(c(seq(N)), invader, seq(spat.equil.times)))  # Empty array to hold invader abundances as we determine spatial equilibrium results.
  space.equil[,,1] <- invade.abundance                                                                                   # Initialize with initial invader abundance calculated from P1/P2 run.
  
  spat.equil.check <- array(NA, dim = c(N, spat.equil.times, length(equil.times)))                                       # Object to hold % differences as invader abundances are updated. Checks spatial equilibrium has been reached.
  
  resident.invader.t0 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)),                       # Empty object to hold host states at t_k0.
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  resident.invader.t1 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)),                       # Empty object to hold host states at t_k1.
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  
  # Calculate invader abundances at spatial equilibrium for each time point the resident is at spatial equilibrium, t_k, and then calculate host states at t_k1.
  for (y in 1:length(equil.times)) {
    invade.abundance.t0 <- invade.abundance   # Use initial invader abundance to start calculating invader spatial equilibrium for each t_k.
    
    # Calculate spatial equilibrium for current time by updating invader abundance 'spat.equil.times' times.
    for (r in 1:(spat.equil.times-1)) {                 
      temp.space.equib <- resident.equilibrium[,,y]
      temp.space.equib[,invader] <- invade.abundance.t0
      
      model.outcome <- betweenHost(host.statuses.input = temp.space.equib, 
                                   host.contacts = host.contacts.eq[[y]], 
                                   coef.logit = coef.logit,
                                   host.trans.probs = host.trans.probs.eq[[y]],
                                   N = N,
                                   times = times,
                                   bottle = bottle, 
                                   eta = eta,
                                   theta = theta,
                                   Rmax = Rmax,
                                   delta1 = delta1,
                                   delta2 = delta2,
                                   alpha1 = alpha1,
                                   alpha2 = alpha2,
                                   beta1 = beta1,
                                   beta2 = beta2,
                                   epsilon = epsilon,
                                   sigma1 = sigma1,
                                   sigma2 = sigma2,
                                   gamma = gamma,
                                   rho = rho,
                                   tau = tau,
                                   mu = mu)
      
      
      invade.abundance.t1 <- model.outcome[,invader] / sum(model.outcome[,invader]) * iabun   # Calculate updated invader abundance.
      
      space.equil[,,r+1] <- invade.abundance.t1   # Input updated invader abundance to 'space.equil' object.
      spat.equil.check[,r,y] <- round(((invade.abundance.t1 - invade.abundance.t0) * 100 / invade.abundance.t0), digits = 3)   # Calculate % difference to check spatial equilibrium is eventually reached.
      
      invade.abundance.t0 <- invade.abundance.t1  # Update invader abundance and repeat for 'spat.equil.times' times.
      
    }
    
    # Using the invading abundance at spatial equilibrium, calculate one more time step.
    resident.invader.t0[,,y] <- resident.equilibrium[,,y]                # Setup t_k0 host states without invader.
    resident.invader.t0[,invader,y] <- space.equil[,,spat.equil.times]   # Invade using the invader abundances (presumably now at spatial equilibrium) from the last run in the 'space.equil' object.
    resident.invader.t1[,,y] <- betweenHost(host.statuses.input = resident.invader.t0[,,y],   # Simulate one time step forward to calculate host states at t_k1.
                                            host.contacts = host.contacts.eq[[y]], 
                                            coef.logit = coef.logit,
                                            host.trans.probs = host.trans.probs.eq[[y]],
                                            N = N,
                                            times = times,
                                            bottle = bottle, 
                                            eta = eta,
                                            theta = theta,
                                            Rmax = Rmax,
                                            delta1 = delta1,
                                            delta2 = delta2,
                                            alpha1 = alpha1,
                                            alpha2 = alpha2,
                                            beta1 = beta1,
                                            beta2 = beta2,
                                            epsilon = epsilon,
                                            sigma1 = sigma1,
                                            sigma2 = sigma2,
                                            gamma = gamma,
                                            rho = rho,
                                            tau = tau,
                                            mu = mu)
    
    
  }
  
  return.list <- list("t0" = resident.invader.t0, "t1" = resident.invader.t1, "check" = spat.equil.check)
  return(return.list)
  
}

#########################

# Function to calculate host states at t0 and t1 with different forms of variation in host states removed, depending on decomposition of interest.
calculateCoexist <- function(resident.equilibrium,
                             host.status,
                             status.categ,
                             equil.times, 
                             spat.equil.times, 
                             t,
                             host.contacts.eq,
                             coef.logit,
                             host.trans.probs.eq,
                             N,
                             times,
                             bottle, 
                             eta,
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
                             gamma, 
                             rho, 
                             tau, 
                             mu) {
  
  # Create objects that will hold host states at t0 and t1.
  resident.invader.t0 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)),    # Empty object to hold host states at t_k0. 
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  resident.invader.t1 <- array(host.status, dim = c(N, length(status.categ), length(equil.times)),    # Empty object to hold host states at t_k1. 
                               dimnames = list(c(seq(N)), status.categ, c(1:(length(equil.times)))))
  
  
  # Calculate host states at t_k1. 
  # NOTE: invader abundances are still at spatial equilibrium.
  for (y in 1:length(equil.times)) {
    # Calculate one more time step.
    resident.invader.t0[,,y] <- resident.equilibrium[,,y]
    resident.invader.t1[,,y] <- betweenHost(host.statuses.input = resident.invader.t0[,,y], 
                                            host.contacts = host.contacts.eq[[y]], 
                                            coef.logit = coef.logit,
                                            host.trans.probs = host.trans.probs.eq[[y]],
                                            N = N,
                                            times = times,
                                            bottle = bottle, 
                                            eta = eta,
                                            theta = theta,
                                            Rmax = Rmax,
                                            delta1 = delta1,
                                            delta2 = delta2,
                                            alpha1 = alpha1,
                                            alpha2 = alpha2,
                                            beta1 = beta1,
                                            beta2 = beta2,
                                            epsilon = epsilon,
                                            sigma1 = sigma1,
                                            sigma2 = sigma2,
                                            gamma = gamma,
                                            rho = rho,
                                            tau = tau,
                                            mu = mu)
    
  }
  
  return.list <- list("t0" = resident.invader.t0, "t1" = resident.invader.t1)
  return(return.list)
  
}


# Function to calculate growth rate.
calculateR <- function(pop.t0, pop.t1) {
  tempt0 <- apply(pop.t0, c(2,3), sum)
  tempt1 <- apply(pop.t1, c(2,3), sum)
  
  tempR <- log(tempt1/tempt0)
  R.bar.result <- apply(tempR, 1, mean)

  return(R.bar.result)
}



#############################################################################################

#### DIAGNOSTIC FUNCTIONS #####


# Function that graphs  model results.
graphBetween <- function(t, results, i) {
  numSus <- numP1 <- numP2 <- numCoI <- c()
  times <- c(1:t)
  
  for (a in 1:t) {
    numP1[a] <- sum(results[,"P1",a] > 100)
    numP2[a] <- sum(results[,"P2",a] > 100)
    numSus[a] <- sum(results[,"P1",a] == 0 & results[,"P2",a] == 0)
    numCoI[a] <- sum(results[,"P1",a] > 0 & results[,"P2",a] > 100)
}
  
  numPop <- as.data.frame(cbind(numSus, numP1, numP2, times))
  
  ggplot(numPop, aes(x = times)) +
    geom_line(aes(y = numSus), color = "blue", size = 1) +
    geom_line(aes(y = numP1), color = "red", size = 1) +
    geom_line(aes(y = numP2), color = "green", size = 1) +
    geom_line(aes(y = numCoI), color = "black", size = 1) +
    # ggtitle(paste("SIMULATION #", i, sep = "")) +
    theme_classic()
  
}

graphWithin <- function(t, results) {
  totP1 <- totP2 <- totI1 <- totI2 <- totM1 <- totM2 <- totR <- c()
  times <- c(1:t)
  
  for (a in 1:t) {
    totP1[a] <- sum(results[,"P1",a])
    totP2[a] <- sum(results[,"P2",a])
    totI1[a] <- sum(results[,"I1",a])
    totI2[a] <- sum(results[,"I2",a])
    totM1[a] <- sum(results[,"M1",a])
    totM2[a] <- sum(results[,"M2",a])
    totR[a] <- sum(results[,"R",a])
  }
  
  totWithin <- as.data.frame(cbind(totP1, totP2, totI1, totI2, totM1, totM2, totR, times))
  
  ggplot(totWithin, aes(x = times)) +
    geom_line(aes(y = totP1), color = "red", size = 1) +
    geom_line(aes(y = totI1), color = "green", size = 1) +
    geom_line(aes(y = totM1), color = "blue", size = 1) +
    geom_line(aes(y = totP2), color = "red", linetype = 2, size = 1) +
    geom_line(aes(y = totI2), color = "green", linetype = 2, size = 1) +
    geom_line(aes(y = totM2), color = "blue", linetype = 2, size = 1) +
    # geom_line(aes(y = totR), color = "black", linetype = 1, size = 1) +
    # ggtitle(paste("SIMULATION #", i, sep = "")) +
    theme_classic()
  
}

wnHostDiag <- function(N, host, times, results) {

  diagHost <- t(results[host,,times])
  diagHost <- as.data.frame(diagHost)
  
  ggplot(diagHost, aes(x = times)) +
    geom_line(aes(y = P1), color = "orange", size = 1) +
    # geom_line(aes(y = I1), color = "green", size = 1) +
    # geom_line(aes(y = M1), color = "blue", size = 1) +
    geom_line(aes(y = P2), color = "blue", linetype = 1, size = 1) +
    # geom_line(aes(y = I2), color = "green", linetype = 2, size = 1) +
    # geom_line(aes(y = M2), color = "blue", linetype = 2, size = 1) +
    # geom_line(aes(y = R), color = "black", linetype = 1, size = 1) +
    ggtitle(paste("HOST #", host, sep = "")) +
    theme_classic()
}
  
    
