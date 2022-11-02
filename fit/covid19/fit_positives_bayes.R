# Packages ----------------------------------------------------------------

packs <- c("tidyverse", "magrittr", "rstan")
lapply(packs, require, character.only = T)
rm(list = ls())

Sys.setlocale(locale = "C") 

# Stan options
mc.cores = parallel::detectCores()
rstan_options(auto_write = TRUE)

# Load richards functions
stan_model <- stan_model("richards_bayes/RichSTCar_NoCovsNew.stan")

# Data preparation --------------------------------------------------------
load("data/dati_omicronfive.RData")
load("data/adjacency_italy.RData")
load("data/offset_italianRegions.RData")


# Input model preparation --------------------------------------------------------

Y <- as.integer(dati$NP)
N <- as.integer(length(Y))
Nreg <- as.integer(length(unique(dati$denominazione_regione)))
timeIdx <- as.integer(dati$WW) + 1
regIdx <- as.integer(droplevels(dati$denominazione_regione))
t <- unique(dati$WW) + 1
Ntimes <- as.integer(length(t))

lOff <- rep(logE, each = Ntimes)

# Prepare data
dat1 <- list(
  "N" = N, 
  "W" = adj,
  "W_n" = sum(adj[upper.tri(adj)]>0),
  "Y" = Y,
  "Ntimes" = Ntimes,
  "t" = t,
  "tId" = timeIdx,
  "Nreg" = Nreg, 
  "rId" = regIdx,
  "lOff" = lOff
)

# Model fitting -----------------------------------------------------------
# Chains
n_chains <- 2
M <- 25000
n_cores <- mc.cores - 2

# Fit
fit_Stan <- sampling(stan_model, data = dat1, chains = n_chains, iter = M, cores = n_cores) # 2 hours

summary(fit_Stan, pars = c("logbase", "logr", "logh", "c", "logs", "alpha", "rho"))