library("ggplot2")
library("gridExtra")
library("tidyverse")
library("reshape2")
library("here")

source(here("code/model_functions.R"))


#----------------------------------------------------------------------------------------------------------------



analyzeDecomps <- function(numSims) {
  
  tot.p1.invader.decomp.indiv <- 
    tot.p1.invader.decomp.compare <- 
    tot.p2.invader.decomp.indiv <-
    tot.p2.invader.decomp.compare <- 
    list()
  
  for (b in 1:numSims) {
    if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
      load(paste("decomp_results_", b, ".RData", sep = ""))
      load(paste("param_results_", b, ".RData", sep = ""))
      
      tot.parameter.list[[b]] <- param.results$parameter.list
      tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
      tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
      tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
      tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
      
    }
    
  }
  
  temp1.p1 <- as.data.frame(plyr::compact(tot.p1.invader.decomp.compare))
  temp1.p2 <- as.data.frame(plyr::compact(tot.p2.invader.decomp.compare))
  
  temp2.p1 <- data.frame(t(temp1.p1), row.names = NULL)
  temp2.p2 <- data.frame(t(temp1.p2), row.names = NULL)
  
  temp.mean.p1 <- unlist(lapply(temp2.p1, mean))
  temp.mean.p2 <- unlist(lapply(temp2.p2, mean))
  
  
  std.err <- function(x) sd(x)/sqrt(length(x))
  temp.se.p1 <- unlist(lapply(temp2.p1, sd))
  temp.se.p2 <- unlist(lapply(temp2.p2, sd))
  
  ldgr.iden <- as.factor(c(1,0,0,0,0))
  decomp <- c("LDGR", "e0", "eL", "eD", "eL*D")
  
  temp.decomp.p1 <- data.frame("decomp" = decomp, "mean" = temp.mean.p1, "se" = temp.se.p1, ldgr.iden)
  temp.decomp.p2 <- data.frame("decomp" = decomp, "mean" = temp.mean.p2, "se" = temp.se.p2, ldgr.iden)
  temp.decomp.p1$decomp <- factor(temp.decomp.p1$decomp, levels = temp.decomp.p1$decomp)
  temp.decomp.p2$decomp <- factor(temp.decomp.p2$decomp, levels = temp.decomp.p2$decomp)
  
  return(list("temp.decomp.p1" = temp.decomp.p1, "temp.decomp.p2" = temp.decomp.p2))
  
}
  
  
###### FIGURE 2: NEUTRAL SCENARIO #######
#### NEUTRALITY WITH GENERAL IMMUNITY #####

numSims <- 500
setwd("./results/neutral_general")

output <- analyzeDecomps(numSims = numSims)
neutral.general.decomp.p1 <- output$temp.decomp.p1
neutral.general.decomp.p2 <- output$temp.decomp.p2


setwd("./../..")



#### NEUTRALITY WITH SPECIFIC IMMUNITY #####

numSims <- 500
setwd("./results/neutral_specific")

output <- analyzeDecomps(numSims = numSims)
neutral.specific.decomp.p1 <- output$temp.decomp.p1
neutral.specific.decomp.p2 <- output$temp.decomp.p2

setwd("./../..")

#################

neutral.general.p1 <- ggplot(data = neutral.general.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
                          geom_bar(stat = "identity", position = "dodge", colour = "black") +
                          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                          scale_fill_manual(values = c("gray80", "gray40")) +
                          scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
                          theme_classic()


neutral.general.p2 <- ggplot(data = neutral.general.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
                          geom_bar(stat = "identity", position = "dodge", colour = "black") +
                          geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                          scale_fill_manual(values = c("gray80", "gray40")) +
                          scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
                          theme_classic()



neutral.specific.p1 <- ggplot(data = neutral.specific.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
                           geom_bar(stat = "identity", position = "dodge", colour = "black") +
                           geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                           scale_fill_manual(values = c("gray80", "gray40")) +
                           scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
                           theme_classic()
  
neutral.specific.p2 <- ggplot(data = neutral.specific.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
                            geom_bar(stat = "identity", position = "dodge", colour = "black") +
                            geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                            scale_fill_manual(values = c("gray80", "gray40")) +
                            scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
                            theme_classic()


grid.arrange(neutral.general.p1, neutral.general.p2,
             neutral.specific.p1, neutral.specific.p2, ncol = 2)


#-----------------------------------------------------------------------------------------

###### FIGURE 3: COMPETITION/COLONIZATION SCENARIO #######
#### COMPETITION/COLONIZATION WITH GENERAL IMMUNITY #####

numSims <- 500
setwd("./results/comp_col_general")

output <- analyzeDecomps(numSims = numSims)
comp.col.general.decomp.p1 <- output$temp.decomp.p1
comp.col.general.decomp.p2 <- output$temp.decomp.p2

setwd("./../..")


#### COMPETITION/COLONIZATION WITH SPECIFIC IMMUNITY #####

numSims <- 500
setwd("./results/comp_col_specific")

output <- analyzeDecomps(numSims = numSims)
comp.col.specific.decomp.p1 <- output$temp.decomp.p1
comp.col.specific.decomp.p2 <- output$temp.decomp.p2

setwd("./../..")

#################

comp.col.general.p1 <- ggplot(data = comp.col.general.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
  theme_classic()

comp.col.general.p2 <- ggplot(data = comp.col.general.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
  theme_classic()



comp.col.specific.p1 <- ggplot(data = comp.col.specific.decomp.p1, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
  theme_classic()

comp.col.specific.p2 <- ggplot(data = comp.col.specific.decomp.p2, aes(x = decomp, y = mean, fill = ldgr.iden)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_fill_manual(values = c("gray80", "gray40")) +
  scale_y_continuous(limits = c(-10,11), breaks = c(-10, -5, 0, 5, 10)) +
  theme_classic()


grid.arrange(comp.col.general.p1, comp.col.general.p2,
             comp.col.specific.p1, comp.col.specific.p2, ncol = 2)




#---------------------------------------------------------------------------------------------------


##### SUPPLEMENTAL FIGURE 3: SENSITIVITY ANALYSIS ######

tot.p1.invader.decomp.indiv <- 
  tot.p1.invader.decomp.compare <- 
  tot.p2.invader.decomp.indiv <-
  tot.p2.invader.decomp.compare <- 
  list()

tot.parameter.list <- 
  tot.coef.logit <-
  list()

numSims <- 511
setwd("./results/sens_analysis")


for (b in 1:numSims) {
  if (file.exists(paste("decomp_results_", b, ".RData", sep = "")) == TRUE) {
    load(paste("decomp_results_", b, ".RData", sep = ""))
    load(paste("param_results_", b, ".RData", sep = ""))
    
    tot.p1.invader.decomp.indiv[[b]] <- decomp.results$p1.invader.decomp.indiv
    tot.p1.invader.decomp.compare[[b]] <- decomp.results$p1.invader.decomp.compare
    tot.p2.invader.decomp.indiv[[b]] <- decomp.results$p2.invader.decomp.indiv
    tot.p2.invader.decomp.compare[[b]] <- decomp.results$p2.invader.decomp.compare
    
    param.results$parameter.list$equil.times <- length(param.results$parameter.list$equil.times)
    tot.parameter.list[[b]] <- param.results$parameter.list
    tot.coef.logit[[b]] <- param.results$coef.logit
    
  }
  
}

tot.parameter.list <- lapply(tot.parameter.list, unlist)
tot.coef.logit <- lapply(tot.coef.logit, unlist)

tot.parameter.list <- plyr::compact(tot.parameter.list)
tot.coef.logit <- plyr::compact(tot.coef.logit)
tot.p1.invader.decomp.indiv <- plyr::compact(tot.p1.invader.decomp.indiv)
tot.p1.invader.decomp.compare <- plyr::compact(tot.p1.invader.decomp.compare)
tot.p2.invader.decomp.indiv <- plyr::compact(tot.p2.invader.decomp.indiv)
tot.p2.invader.decomp.compare <- plyr::compact(tot.p2.invader.decomp.compare)

tot.parameter.list <- Map(c, tot.parameter.list, tot.coef.logit)

param.names <- append(names(tot.parameter.list[[1]]), names(tot.p1.invader.decomp.compare[[1]]))


sens.results <- matrix(nrow = length(tot.parameter.list), ncol = length(param.names), dimnames = list(seq(length(tot.parameter.list)), 
                                                                                                      param.names))

for (c in 1:length(tot.parameter.list)) {
  sens.results[c,] <- append(tot.parameter.list[[c]], tot.p1.invader.decomp.compare[[c]])
}

sens.results <- as.data.frame(sens.results)
colnames(sens.results)[colnames(sens.results) %in% c("P1.(Intercept)", "P2.(Intercept)")] <- c("P1.inter", "P2.inter")

sens.results <- sens.results[sens.results$LDGR != Inf,]
sens.results <- sens.results[!is.na(sens.results$LDGR),]

#####

fit <- glm(LDGR ~ theta +
             Rmax +
             delta1 +
             delta2 +
             alpha1 +
             alpha2 +
             beta1 +
             beta2 +
             epsilon +
             sigma1 +
             sigma2 +
             gamma +
             rho +
             tau +
             mu +
             P1.inter +
             P1.load +
             P2.inter +
             P2.load +
             kappa +
             bottle +
             eta,
           
           data = sens.results)



fit.coef <- summary(fit)$coef[-1,1]
sd.err <- summary(fit)$coef[-1,2]

sx <- as.numeric(lapply(fit$model[-1], sd))
sy <- as.numeric(lapply(fit$model[1], sd))
beta <- fit.coef * sx/sy
stand.sd.err <- sd.err * sx/sy
stand.sd.dev <- stand.sd.err * sqrt(nrow(sens.results))

params <- names(beta)
coefficients <- as.numeric(unname(beta))
fit.stand <- data.frame("params" = params, "coefficients" = coefficients, "se" = stand.sd.err)
fit.stand$params <- factor(fit.stand$params, levels = fit.stand$params)


ggplot(data = fit.stand, aes(x = params, y = coefficients)) + 
  geom_bar(stat = "identity", fill = "gray80", colour = "black") +
  geom_errorbar(aes(ymin = coefficients - se, ymax = coefficients + se), width = 0.2) +
  labs(title = "Sensitivity Analysis w/ SE") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


setwd("./../..")
# write.csv(fit.stand, file = "sens_analysis.csv")
