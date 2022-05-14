##### SENSITIVITY ANALYSIS #####
N <- 1000   # See decomposition.R for descriptions of each parameter.
t <- 2000           

sd.var <- 0.25

# Between-host model parameters ##
kappa <- rnorm(1, mean = 2, sd = 2 * sd.var)
bottle <- rnorm(1, mean = 0.05, sd = 0.05 * sd.var)
eta <- rnorm(1, mean = 0.001, sd = 0.001 * sd.var)

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

theta <- rnorm(1, mean = 10, sd = 10 * sd.var)
Rmax <- rnorm(1, mean = 100, sd = 100 * sd.var)
delta1 <- rnorm(1, mean = 0.001, sd = 0.001 * sd.var)
delta2 <- rnorm(1, mean = 0.001, sd = 0.001 * sd.var)
alpha1 <- rnorm(1, mean = 100, sd = 100 * sd.var)
alpha2 <- rnorm(1, mean = 100, sd = 100 * sd.var)
beta1 <- rnorm(1, mean = 0.01, sd = 0.01 * sd.var)
beta2 <- rnorm(1, mean = 0.01, sd = 0.01 * sd.var)
epsilon <- rnorm(1, mean = 0.00001, sd = 0.00001 * sd.var)
sigma1 <- rnorm(1, mean = 0.005, sd = 0.005 * sd.var)
sigma2 <- rnorm(1, mean = 0.005, sd = 0.005 * sd.var)
tau <- rnorm(1, mean = 0.0025, sd = 0.0025 * sd.var)
gamma <- rnorm(1, mean = 0.05, sd = 0.05 * sd.var)
rho <- rnorm(1, mean = 0.0005, sd = 0.0005 * sd.var)
mu <- rnorm(1, mean = 0.0005, sd = 0.0005 * sd.var)

# Coexistence decomposition parameters.
equil.times <- c(1900:2000)
spat.equil.times <- 10


########################


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


##############################


##### COMPETITION/COLONIZATION TRADEOFF #####
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
med.k.p2 <- 0.05
hi.k.p1 <- 0.05
hi.k.p2 <- 0.1


## Within-host model parameters ##
times <- seq(0, 1, by = 1)

theta <- 10
Rmax <- 100
delta1 <- 0.002
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