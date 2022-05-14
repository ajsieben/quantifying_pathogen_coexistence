#### COMPLETE MODEL WITH DECOMPOSITIONS ####  -----------------------------------------------------------------------------
# library("deSolve")
# library("ggplot2")
# 
# 
# source("./code/model_functions.R")
# source("./code/within_host_model.R")
# source("./code/between_host_model.R")

source("model_functions.R")
source("within_host_model.R")
source("between_host_model.R")


runModel <- function(v) {

  # residents.established <- FALSE
  # while (residents.established == FALSE) {
    
    #### DEFINE PARAMETERS ####
    parameter.list <- generateParams()
    
    ## Between-host model parameters ##
    N <- parameter.list$N                          # The total population size.
    t <- parameter.list$t                          # The amount of time the model is simulated.
    sd.var <- parameter.list$sd.var                # The amount of variation introduced for the sensitivity analysis.

    kappa <- parameter.list$kappa                  # Mean number of contacts.
    bottle <- parameter.list$bottle                # Probability a pathogen will colonize the new host upon successful transmission..
    eta <- parameter.list$eta                      # Probability a host will die during each time step.

    
    lo.n <- parameter.list$lo.n                    # Pathogen abundance defined as "low".
    med.n <- parameter.list$med.n                  # Pathogen abundance defined as "medium".
    hi.n <- parameter.list$hi.n                    # Pathogen abundance defined as "high".
    
    lo.k.p1 <- parameter.list$lo.k.p1              # Probability of transmission at low pathogen abundance for P1.
    lo.k.p2 <- parameter.list$lo.k.p2              # Probability of transmission at low pathogen abundance for P2.
    med.k.p1 <- parameter.list$med.k.p1            # Probability of transmission at medium pathogen abundance for P1.
    med.k.p2 <- parameter.list$med.k.p2            # Probability of transmission at medium pathogen abundance for P2.
    hi.k.p1 <- parameter.list$hi.k.p1              # Probability of transmission at high pathogen abundance for P1.
    hi.k.p2 <- parameter.list$hi.k.p2              # Probability of transmission at high pathogen abundance for P2.
    
    ## Within-host model parameters ##
    times <- parameter.list$times                  # The number of within-host time steps for each between-host time step.
    
    theta <- parameter.list$theta                  # Base rate of host resource recovery.
    Rmax <- parameter.list$Rmax                    # Maximum host resource.
    delta1 <- parameter.list$delta1                # Resource conversion rate for P1.
    delta2 <- parameter.list$delta2                # Resource conversion rate for P2.
    alpha1 <- parameter.list$alpha1                # P1 rate of increase.
    alpha2 <- parameter.list$alpha2                # P2 rate of increase.
    beta1 <- parameter.list$beta1                  # Rate of consumption by I1.
    beta2 <- parameter.list$beta2                  # Rate of consumption by I2.
    epsilon <- parameter.list$epsilon              # Conversion factor for immune cell growth.
    sigma1 <- parameter.list$sigma1                # Rate of immune system activation by memory cells from same serotype.
    sigma2 <- parameter.list$sigma2                # Rate of immune system activation by memory cells from same serotype.
    gamma <- parameter.list$gamma                  # Rate of immune system decay.
    rho <- parameter.list$rho                      # Rate of memory system stimulation.
    tau <- parameter.list$tau                      # Rate of immune system activation by memory cells from different serotype (i.e. cross-immunity).
    mu <- parameter.list$mu                        # Rate of memory system decay.
    
    ## Decomposition parameters ##
    equil.times <- parameter.list$equil.times             # Object specifying the time range that residents are at equilibrium.
    spat.equil.times <- parameter.list$spat.equil.times   # The number of iterations for the invader to reach spatial equilibrium. 
    
    ##########
    
    ## Setup ##
    # Generate a matrix to hold all host statuses for one time step.
    initP1.load <- 100    # Initial pathogen load for those infected with P1.
    initP2.load <- 100    # Initial pathogen load for those infected with P2.
    initI1.load <- 1      # Initial I1 load.
    initI2.load <- 1      # Initial I2 load.
    initR.load <- 100     # Initial R load.
    
    status.categ <- c("P1", "P2", "I1", "I2", "M1", "M2", "R")

    host.status <- matrix(0, nrow = N, ncol = length(status.categ), 
                          dimnames = list(seq(N), status.categ))  
    host.status[, "I1"] <- initI1.load
    host.status[, "I2"] <- initI2.load
    host.status[, "R"] <- initR.load

    
    # -----------------------------------------------------------------
    
    # Determine host contacts for all time steps.
    host.contacts <- list(NULL)
    for (i in 1:t) {
      num.contacts <- rpois(N, kappa)                 # Generate vector specifying how many other hosts each host contacts.
      tot.contacts <- sum(num.contacts)
      
      # Create matrix with first column listing contactor and second column listing contacted.
      temp.host.contacts <- matrix(data = NA, nrow = tot.contacts, ncol = 2, dimnames = list(seq(tot.contacts), 
                                                                                             c("contacter", "contacted")))
      
      vector.contacters <- rep(seq(N), num.contacts)  # Fill first column with contactors.
      list.contacted <- list(NULL)
        
      for (j in 1:N) {
        contactables <- seq(N)
        contacted <- sample(contactables[-j], num.contacts[j], replace = FALSE)   # Draw contacted hosts without replacement (also remove contactor).
        
        if (length(contacted) != 0) {
          list.contacted[[j]] <- contacted                                        # If respective host successfully contacts others, add to list.
        }
      }
      
      vector.contacted <- unlist(list.contacted)
      
      temp.host.contacts[,"contacter"] <- vector.contacters   # Combine contactor and contacted vectors.
      temp.host.contacts[,"contacted"] <- vector.contacted
      
      host.contacts[[i]] <- temp.host.contacts                # Add to 'host.contact' object, which will have a different contact list for each t.
      
    }
    
    
    # Calculate logit functions for P1 and P2 used to determine pathogen transmission.
    coef.logit.p1 <- calcLogit(lo.k = lo.k.p1, 
                              lo.n = lo.n, 
                              med.k = med.k.p1, 
                              med.n = med.n, 
                              hi.k = hi.k.p1, 
                              hi.n = hi.n)
    
    coef.logit.p1["(Intercept)"] <- rnorm(1, mean = coef.logit.p1["(Intercept)"], sd = abs(coef.logit.p1["(Intercept)"] * sd.var))   # For sensitivity analysis, vary logit outcomes.
    coef.logit.p1["load"] <- rnorm(1, mean = coef.logit.p1["load"], sd = coef.logit.p1["load"] * sd.var)
    
    coef.logit.p2 <- calcLogit(lo.k = lo.k.p2, 
                               lo.n = lo.n, 
                               med.k = med.k.p2, 
                               med.n = med.n, 
                               hi.k = hi.k.p2, 
                               hi.n = hi.n)
    
    coef.logit.p2["(Intercept)"] <- rnorm(1, mean = coef.logit.p2["(Intercept)"], sd = abs(coef.logit.p2["(Intercept)"] * sd.var))
    coef.logit.p2["load"] <- rnorm(1, mean = coef.logit.p2["load"], sd = coef.logit.p2["load"] * sd.var)
    
    coef.logit <- list("P1" = coef.logit.p1, "P2" = coef.logit.p2)
    
    # Calculate random uniform numbers (for b/w transmission probability) and store in object with length of contact list and time.
    host.trans.probs <- list(NULL)
    for (i in 1:t) {
      host.trans.probs[[i]] <- matrix(data = NA, nrow = nrow(host.contacts[[i]]), ncol = 2, dimnames = list(seq(nrow(host.contacts[[i]])), 
                                                                                                  c("P1", "P2")))
      host.trans.probs[[i]][,"P1"] <- runif(n = nrow(host.contacts[[i]]), min = 0, max = 1)
      host.trans.probs[[i]][,"P2"] <- runif(n = nrow(host.contacts[[i]]), min = 0, max = 1)
    }
   

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------
    # -----------------------------------------------------------------
    
    
    
    #### 1. #### DETERMINE PATHOGEN EQUILIBRIUMS WHEN BOTH ARE PRESENT IN THE SYSTEM. #### -------------------------------
    # Specify the initial conditions.
    initP1 <- 10          # Initial number of hosts infected with pathogen 1.
    initP2 <- 10          # Initial number of hosts infected with pathogen 2.
    
    # Array to hold the within-host variables and set the starting conditions.
    results.p1.p2.equilibrium <- array(host.status, dim = c(N, length(status.categ), t), dimnames = list(c(seq(N)), status.categ, seq(t)))
    results.p1.p2.equilibrium[,,1] <- initialCondit(results = results.p1.p2.equilibrium[,,1],
                                                    initP1 = initP1, 
                                                    initP1.load = initP1.load, 
                                                    initP2 = initP2, 
                                                    initP2.load = initP2.load)
    
    ## Model ##
    for (a in 1:(t-1)) {
      results.p1.p2.equilibrium[,,a+1] <- betweenHost(host.statuses.input = results.p1.p2.equilibrium[,,a],
                                                      host.contacts = host.contacts[[a]],
                                                      coef.logit = coef.logit,
                                                      host.trans.probs = host.trans.probs[[a]],
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

    # graphBetween(t = t, results = results.p1.p2.equilibrium)
    # graphWithin(t = t, results = results.p1.p2.equilibrium)
    
    # wnHostDiag(N = N, host = 1, times = c(1:t), results = results.p1.p2.equilibrium)
    
  #### 2. #### RESIDENT EQUILIBRIUMS ####  -----------------------------------------------------------------------------
  
    ### P1 as resident ###
    # Specify the initial conditions.
    initP1 <- 10          # Initial number of hosts infected with pathogen 1.
    initP2 <- 0           # Initial number of hosts infected with pathogen 2.
    
    # Array to hold the within-host variables and set the starting conditions.
    results.p1.resident <- array(host.status, dim = c(N, length(status.categ), t), dimnames = list(c(seq(N)), status.categ, seq(t)))
    results.p1.resident[,,1] <- initialCondit(results = results.p1.resident[,,1],
                                              initP1 = initP1, 
                                              initP1.load = initP1.load, 
                                              initP2 = initP2, 
                                              initP2.load = initP2.load)
    
    ## Model ##
    
    for (a in 1:(t-1)) {
      results.p1.resident[,,a+1] <- betweenHost(host.statuses.input = results.p1.resident[,,a],
                                                host.contacts = host.contacts[[a]],
                                                coef.logit = coef.logit,
                                                host.trans.probs = host.trans.probs[[a]],
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
    
    # graphBetween(t = t,  results = results.p1.resident)
    # graphWithin(t = t, results = results.p1.resident)
    

    ##########################################################
    
    ### P2 as resident ###
    # Specify the initial conditions.
    initP1 <- 0          # Initial number of hosts infected with pathogen 1.
    initP2 <- 10          # Initial number of hosts infected with pathogen 2.
    
    # Array to hold the within-host variables and set the starting conditions.
    results.p2.resident <- array(host.status, dim = c(N, length(status.categ), t), dimnames = list(c(seq(N)), status.categ, seq(t)))
    results.p2.resident[,,1] <- initialCondit(results = results.p2.resident[,,1],
                                              initP1 = initP1, 
                                              initP1.load = initP1.load, 
                                              initP2 = initP2, 
                                              initP2.load = initP2.load)
    
    ## Model ##
    
    for (a in 1:(t-1)) {
      results.p2.resident[,,a+1] <- betweenHost(host.statuses.input = results.p2.resident[,,a],
                                                host.contacts = host.contacts[[a]],
                                                coef.logit = coef.logit,
                                                host.trans.probs = host.trans.probs[[a]],
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

    # graphBetween(t = t, results = results.p2.resident)
    # graphWithin(t = t, results = results.p2.resident)
    
    
  #   residents.established <- ifelse(sum(results.p1.resident[,"P1", max(equil.times)]) > 100 | sum(results.p2.resident[,"P2", max(equil.times)]) > 100,
  #                                 TRUE, FALSE)
  # 
  # 
  # }
  
  param.results <- list("parameter.list" = parameter.list,
                        "coef.logit" = coef.logit)
  
  raw.results <- list("results.p1.p2.equilibrium" = results.p1.p2.equilibrium,
                      "results.p1.resident" = results.p1.resident,
                      "results.p2.resident" = results.p2.resident)


  save(param.results, file = paste("param_results_", v, ".RData", sep = ""))   # Save parameter data for sensitivity analysis and for checks.
  save(raw.results, file = paste("raw_results_", v, ".RData", sep = ""))
  

  #### 3. #### CALCULATE LOW DENSITY GROWTH RATES WITHOUT PARTITIONING #### ----------------------------------
  
  ### Setup ###
  # Pull results from time frame in which pathogens are at equilibrium.
  p1.p2.equilibrium <- results.p1.p2.equilibrium[,,equil.times, drop = FALSE]
  p1.resident.equilibrium <- results.p1.resident[,,equil.times, drop = FALSE]
  p2.resident.equilibrium <- results.p2.resident[,,equil.times, drop = FALSE]
  
  host.contacts.eq <- host.contacts[equil.times]
  host.trans.probs.eq <- host.trans.probs[equil.times]
  
  
  # Calculate initial invasion abundances for P1 and P2.
  iabun <- 0.005 * (sum(p1.p2.equilibrium[,c("P1","P2"),])/(2*length(equil.times)))
  init.abun <- rep(iabun/N, N)
  
  # If initial invader abundance is too low, manually specify.
  init.abun <- ifelse(init.abun <= 0.02, 0.02, init.abun)          # I set the manual initial invader abundance at twice the threshold for pathogen removal. 
  iabun <- ifelse(init.abun[1] == 0.02, N * 0.02, iabun)           # Recalculate total invader abundance after manual reset.
  
  
  
  ### Calculate GRWR with P1 as resident, P2 as invader ###
  model.final.outcome <- calculateLDGR(resident.equilibrium = p1.resident.equilibrium,
                                       host.status = host.status,
                                       status.categ = status.categ,
                                       invader = "P2", 
                                       equil.times = equil.times, 
                                       spat.equil.times = spat.equil.times, 
                                       iabun = iabun,
                                       invade.abundance = init.abun,
                                       t = t,
                                       host.contacts.eq = host.contacts.eq,
                                       coef.logit = coef.logit,
                                       host.trans.probs.eq = host.trans.probs.eq,
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
  
  
  p1.resident.p2.invader.t0 <- model.final.outcome$t0
  p1.resident.p2.invader.t1 <- model.final.outcome$t1
  spat.equil.check.r.bar.p2 <- model.final.outcome$check
  
  r.bar.p2 <- calculateR(pop.t0 = p1.resident.p2.invader.t0,
                         pop.t1 = p1.resident.p2.invader.t1)
  
  r.bar.p2.invader <- r.bar.p2["P2"]
  r.bar.p1.resident <- r.bar.p2["P1"]
  
  
  ### P2 as resident, P1 as invader ###
  model.final.outcome <- calculateLDGR(resident.equilibrium = p2.resident.equilibrium,
                                       host.status = host.status,
                                       status.categ = status.categ,
                                       invader = "P1", 
                                       equil.times = equil.times, 
                                       spat.equil.times = spat.equil.times, 
                                       iabun = iabun,
                                       invade.abundance = init.abun,
                                       t = t,
                                       host.contacts.eq = host.contacts.eq,
                                       coef.logit = coef.logit,
                                       host.trans.probs.eq = host.trans.probs.eq,
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
  
  p2.resident.p1.invader.t0 <- model.final.outcome$t0
  p2.resident.p1.invader.t1 <- model.final.outcome$t1
  spat.equil.check.r.bar.p1 <- model.final.outcome$check
  
  
  r.bar.p1 <- calculateR(pop.t0 = p2.resident.p1.invader.t0,
                         pop.t1 = p2.resident.p1.invader.t1)
  
  r.bar.p1.invader <- r.bar.p1["P1"]
  r.bar.p2.resident <- r.bar.p1["P2"]
  
  
  
  print("COMPLETED STEP 3")

  
  #### 4. #### PARTITION GROWTH RATE WITH EACH PATHOGEN AS INVADER #### ----------------------------------
  #### 4.1 ### No variation in fitness or density. ###
  
  # Average host states across entire population. Specific averaged conditions will be used from this object, depending on the decomposition.
  p1.resident.p2.invader.all.averaged <- colMeans(p1.resident.p2.invader.t0)
  p2.resident.p1.invader.all.averaged <- colMeans(p2.resident.p1.invader.t0)

  
  # Create object to hold results for this decomposition and pull results from time frame in which pathogens are at equilibrium.
  p1.resident.equilibrium.no.var <- p1.resident.p2.invader.t0
  p2.resident.equilibrium.no.var <- p2.resident.p1.invader.t0
  
  
  ### P1 as resident, P2 as invader ###
  # Overwrite each component abundance (other than P2 invader) to reflect averaged state.
  # NOTE: Invader abundances are not overwritten, thereby remaining at spatial equilibrium for each t_k.
  for (i in 1:length(equil.times)) {
    for (k in c("P1", "I1", "I2", "M1", "M2", "R")) {
      p1.resident.equilibrium.no.var[,k,i] <- p1.resident.p2.invader.all.averaged[k,i]
    }
  }
  
  # Calculate GRWR for each t_k under these decomposition conditions.
  model.final.outcome <- calculateCoexist(resident.equilibrium = p1.resident.equilibrium.no.var,
                                          host.status = host.status,
                                          status.categ = status.categ,
                                          equil.times = equil.times, 
                                          spat.equil.times = spat.equil.times, 
                                          t = t,
                                          host.contacts.eq = host.contacts.eq,
                                          coef.logit = coef.logit,
                                          host.trans.probs.eq = host.trans.probs.eq,
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
  
  p1.resident.p2.invader.no.var.t0 <- model.final.outcome$t0
  p1.resident.p2.invader.no.var.t1 <- model.final.outcome$t1

  
  r.bar.no.var.p2 <- calculateR(pop.t0 = p1.resident.p2.invader.no.var.t0,
                                pop.t1 = p1.resident.p2.invader.no.var.t1)
  
  epsilon.0.p2.invader <- r.bar.no.var.p2["P2"]
  epsilon.0.p1.resident <- r.bar.no.var.p2["P1"]
  
  
  
  ### P2 as resident, P1 as invader ###
  # Overwrite each component abundance (other than P1 invader) to reflect averaged state.
  for (i in 1:length(equil.times)) {
    for (k in c("P2", "I1", "I2", "M1", "M2", "R")) {
      p2.resident.equilibrium.no.var[,k,i] <- p2.resident.p1.invader.all.averaged[k,i]
    }
  }
  
  # Calculate GRWR for each t_k under these decomposition conditions.
  model.final.outcome <- calculateCoexist(resident.equilibrium = p2.resident.equilibrium.no.var,
                                          host.status = host.status,
                                          status.categ = status.categ,
                                          equil.times = equil.times, 
                                          spat.equil.times = spat.equil.times, 
                                          t = t,
                                          host.contacts.eq = host.contacts.eq,
                                          coef.logit = coef.logit,
                                          host.trans.probs.eq = host.trans.probs.eq,
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
  
  p2.resident.p1.invader.no.var.t0 <- model.final.outcome$t0
  p2.resident.p1.invader.no.var.t1 <- model.final.outcome$t1
  
  r.bar.no.var.p1 <- calculateR(pop.t0 = p2.resident.p1.invader.no.var.t0,
                                pop.t1 = p2.resident.p1.invader.no.var.t1)
  
  epsilon.0.p1.invader <- r.bar.no.var.p1["P1"]
  epsilon.0.p2.resident <- r.bar.no.var.p1["P2"]
  
  
  ####################################################################
  
  
  
  #### 4.2 ### Variation in fitness, constant density. ###
  # Create object to hold results for this decomposition and pull results from time frame in which pathogens are at equilibrium.
  p1.resident.equilibrium.fitness.var <- p1.resident.p2.invader.t0
  p2.resident.equilibrium.fitness.var <- p2.resident.p1.invader.t0
  

  ### P1 as resident, P2 as invader ###
  # Overwrite density components (other than P2 invader) to reflect variation in fitness only.
  for (i in 1:length(equil.times)) {
    for (k in c("P1")) {
      p1.resident.equilibrium.fitness.var[,k,i] <- p1.resident.p2.invader.all.averaged[k,i]
    }
  }
  
  # Calculate GRWR for each t_k under these decomposition conditions.
  model.final.outcome <- calculateCoexist(resident.equilibrium = p1.resident.equilibrium.fitness.var,
                                          host.status = host.status,
                                          status.categ = status.categ,
                                          equil.times = equil.times, 
                                          spat.equil.times = spat.equil.times, 
                                          t = t,
                                          host.contacts.eq = host.contacts.eq,
                                          coef.logit = coef.logit,
                                          host.trans.probs.eq = host.trans.probs.eq,
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
  
  p1.resident.p2.invader.fitness.var.t0 <- model.final.outcome$t0
  p1.resident.p2.invader.fitness.var.t1 <- model.final.outcome$t1

  
  r.bar.fitness.p2 <- calculateR(pop.t0 = p1.resident.p2.invader.fitness.var.t0,
                               pop.t1 = p1.resident.p2.invader.fitness.var.t1)
  
  epsilon.fitness.p2.invader <- r.bar.fitness.p2["P2"] - epsilon.0.p2.invader
  epsilon.fitness.p1.resident <- r.bar.fitness.p2["P1"] - epsilon.0.p1.resident
  
  
  
  ### P2 as resident, P1 as invader ###
  # Overwrite density components (other than P1 invader) to reflect variation in fitness only.
  for (i in 1:length(equil.times)) {
    for (k in c("P2")) {
      p2.resident.equilibrium.fitness.var[,k,i] <- p2.resident.p1.invader.all.averaged[k,i]
    }
  }
  
  # Calculate GRWR for each t_k under these decomposition conditions.
  model.final.outcome <- calculateCoexist(resident.equilibrium = p2.resident.equilibrium.fitness.var,
                                          host.status = host.status,
                                          status.categ = status.categ,
                                          equil.times = equil.times, 
                                          spat.equil.times = spat.equil.times, 
                                          t = t,
                                          host.contacts.eq = host.contacts.eq,
                                          coef.logit = coef.logit,
                                          host.trans.probs.eq = host.trans.probs.eq,
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
  
  p2.resident.p1.invader.fitness.var.t0 <- model.final.outcome$t0
  p2.resident.p1.invader.fitness.var.t1 <- model.final.outcome$t1

  
  r.bar.fitness.p1 <- calculateR(pop.t0 = p2.resident.p1.invader.fitness.var.t0,
                               pop.t1 = p2.resident.p1.invader.fitness.var.t1)
  
  epsilon.fitness.p1.invader <- r.bar.fitness.p1["P1"] - epsilon.0.p1.invader
  epsilon.fitness.p2.resident <- r.bar.fitness.p1["P2"] - epsilon.0.p2.resident
  
  
  ####################################################################
  
  
  #### 4.3 ### Variation in density, constant fitness. ###
  # Create object to hold results for this decomposition and pull results from time frame in which pathogens are at equilibrium.
  p1.resident.equilibrium.density.var <- p1.resident.p2.invader.t0
  p2.resident.equilibrium.density.var <- p2.resident.p1.invader.t0
  
  ### P1 as resident, P2 as invader ###
  # Overwrite fitness components to reflect variation in density only.
  for (i in 1:length(equil.times)) {
    for (k in c("I1","I2","M1","M2","R")) {
      p1.resident.equilibrium.density.var[,k,i] <- p1.resident.p2.invader.all.averaged[k,i]
    }
  }
  
  # Calculate GRWR for each t_k under these decomposition conditions.
  model.final.outcome <- calculateCoexist(resident.equilibrium = p1.resident.equilibrium.density.var,
                                          host.status = host.status,
                                          status.categ = status.categ,
                                          equil.times = equil.times, 
                                          spat.equil.times = spat.equil.times, 
                                          t = t,
                                          host.contacts.eq = host.contacts.eq,
                                          coef.logit = coef.logit,
                                          host.trans.probs.eq = host.trans.probs.eq,
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
  
  p1.resident.p2.invader.density.var.t0 <- model.final.outcome$t0
  p1.resident.p2.invader.density.var.t1 <- model.final.outcome$t1

  
  r.bar.density.p2 <- calculateR(pop.t0 = p1.resident.p2.invader.density.var.t0,
                               pop.t1 = p1.resident.p2.invader.density.var.t1)
  
  epsilon.density.p2.invader <- r.bar.density.p2["P2"] - epsilon.0.p2.invader
  epsilon.density.p1.resident <- r.bar.density.p2["P1"] - epsilon.0.p1.resident
  
  
  
  
  ### P2 as resident, P1 as invader ###
  # Overwrite fitness components to reflect variation in density only.
  for (i in 1:length(equil.times)) {
    for (k in c("I1","I2","M1","M2","R")) {
      p2.resident.equilibrium.density.var[,k,i] <- p2.resident.p1.invader.all.averaged[k,i]
    }
  }
  
  # Calculate GRWR for each t_k under these decomposition conditions.
  model.final.outcome <- calculateCoexist(resident.equilibrium = p2.resident.equilibrium.density.var,
                                          host.status = host.status,
                                          status.categ = status.categ,
                                          equil.times = equil.times, 
                                          spat.equil.times = spat.equil.times, 
                                          t = t,
                                          host.contacts.eq = host.contacts.eq,
                                          coef.logit = coef.logit,
                                          host.trans.probs.eq = host.trans.probs.eq,
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
  
  p2.resident.p1.invader.density.var.t0 <- model.final.outcome$t0
  p2.resident.p1.invader.density.var.t1 <- model.final.outcome$t1

  
  r.bar.density.p1 <- calculateR(pop.t0 = p2.resident.p1.invader.density.var.t0,
                               pop.t1 = p2.resident.p1.invader.density.var.t1)
  
  epsilon.density.p1.invader <- r.bar.density.p1["P1"] - epsilon.0.p1.invader
  epsilon.density.p2.resident <- r.bar.density.p1["P2"] - epsilon.0.p2.resident
  
  
  
  ##############################################################################
  
  # Calculate interaction term. #
  epsilon.interaction.p2.invader <- r.bar.p2.invader - (epsilon.0.p2.invader + epsilon.fitness.p2.invader + epsilon.density.p2.invader)
  epsilon.interaction.p1.invader <- r.bar.p1.invader - (epsilon.0.p1.invader + epsilon.fitness.p1.invader + epsilon.density.p1.invader)
  
  epsilon.interaction.p2.resident <- r.bar.p2.resident - (epsilon.0.p2.resident + epsilon.fitness.p2.resident + epsilon.density.p2.resident)
  epsilon.interaction.p1.resident <- r.bar.p1.resident - (epsilon.0.p1.resident + epsilon.fitness.p1.resident + epsilon.density.p1.resident)
  
  
  ##############################################################################
  
  
  #### 4.4 ### Calculate invader/resident comparisons. ####
  
  # Invader only results #
  p1.invader.decomp.indiv <- c(r.bar.p1.invader, 
                               epsilon.0.p1.invader, 
                               epsilon.fitness.p1.invader, 
                               epsilon.density.p1.invader,
                               epsilon.interaction.p1.invader)
  
  p2.invader.decomp.indiv <- c(r.bar.p2.invader, 
                               epsilon.0.p2.invader, 
                               epsilon.fitness.p2.invader, 
                               epsilon.density.p2.invader,
                               epsilon.interaction.p2.invader)
  
  
  names(p1.invader.decomp.indiv) <- c("LDGR", "e0", "eL", "eD", "eL*D")
  names(p2.invader.decomp.indiv) <- c("LDGR", "e0", "eL", "eD", "eL*D")
  
  # Invader-resident comparisons #
  p1.invader.decomp.compare <- c(r.bar.p1.invader - r.bar.p2.resident, 
                                 (epsilon.0.p1.invader - epsilon.0.p2.resident),
                                 (epsilon.fitness.p1.invader - epsilon.fitness.p2.resident),
                                 (epsilon.density.p1.invader - epsilon.density.p2.resident),
                                 (epsilon.interaction.p1.invader - epsilon.interaction.p2.resident))
  
  p2.invader.decomp.compare <- c(r.bar.p2.invader - r.bar.p1.resident,
                                 (epsilon.0.p2.invader - epsilon.0.p1.resident),
                                 (epsilon.fitness.p2.invader - epsilon.fitness.p1.resident),
                                 (epsilon.density.p2.invader - epsilon.density.p1.resident),
                                 (epsilon.interaction.p2.invader - epsilon.interaction.p1.resident))
  
  names(p1.invader.decomp.compare) <- c("LDGR", "e0", "eL", "eD", "eL*D")
  names(p2.invader.decomp.compare) <- c("LDGR", "e0", "eL", "eD", "eL*D")
  
  
  #### 5. #### GENERATING OUTPUT #### ----------------------------------
  
  decomp.results <- list("p1.invader.decomp.indiv" = p1.invader.decomp.indiv,
                         "p1.invader.decomp.compare" = p1.invader.decomp.compare,
                         "p2.invader.decomp.indiv" = p2.invader.decomp.indiv,
                         "p2.invader.decomp.compare" = p2.invader.decomp.compare)


  save(decomp.results, file = paste("decomp_results_", v, ".RData", sep = ""))

}

##########################################


##### TESTING CODE ########

# results.p1.p2.equilibrium <- raw.results$results.p1.p2.equilibrium
# results.p1.resident <- raw.results$results.p1.resident
# results.p2.resident <- raw.results$results.p2.resident
# 
# p1.invader.decomp.compare <- decomp.results$p1.invader.decomp.compare
# p2.invader.decomp.compare <- decomp.results$p2.invader.decomp.compare
# 
# parameter.list <- param.results$parameter.list
# coef.logit <- param.results$coef.logit
# 
# ###
# 
# ldgr.iden <- as.factor(c(1,0,0,0,0))
# decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")
# temp1.p1 <- as.data.frame(p1.invader.decomp.compare)
# temp1.p2 <- as.data.frame(p2.invader.decomp.compare)
# 
# temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
# temp2.p2 <- data.frame(t(temp1.p2), row.names = NULL)
# 
# temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
# temp.mean.p2 <- unlist(lapply(temp2.p2, mean))
# 
# neutral.specific.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, ldgr.iden)
# neutral.specific.decomp.p2 <- data.frame("decomp" = decomp, "mean" = temp.mean.p2, ldgr.iden)
# 
# neutral.specific.decomp.p1$decomp <- factor(neutral.specific.decomp.p1$decomp, levels = neutral.specific.decomp.p1$decomp)
# neutral.specific.decomp.p2$decomp <- factor(neutral.specific.decomp.p2$decomp, levels = neutral.specific.decomp.p2$decomp)
# 
# 
# ggplot(data = neutral.specific.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
#   geom_bar(stat = "identity", position = "dodge", colour = "black") +
#   scale_fill_manual(values = c("gray80", "gray40")) +
#   # scale_y_continuous(limits = c(-3,9), breaks = c(-3, 0, 3, 6, 9)) +
#   theme_classic()
# 
# 
# ggplot(data = neutral.specific.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
#   geom_bar(stat = "identity", position = "dodge", colour = "black") +
#   scale_fill_manual(values = c("gray80", "gray40")) +
#   # scale_y_continuous(limits = c(-3,9), breaks = c(-3, 0, 3, 6, 9)) +
#   theme_classic()