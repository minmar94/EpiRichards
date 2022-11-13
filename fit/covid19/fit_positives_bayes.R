# Packages ----------------------------------------------------------------

packs <- c("tidyverse", "magrittr", "rstan")
lapply(packs, require, character.only = T)
rm(list = ls())

Sys.setlocale(locale = "C") 

# Stan options
mc.cores = parallel::detectCores()
rstan_options(auto_write = TRUE)

# Load richards functions
stan_model <- stan_model("richards_bayes/RichSTCar_NoCovsNew2.stan")

# Data preparation --------------------------------------------------------
load("data/dati_omicronfive.RData")
load("data/adjacency_italy.RData")
load("data/offset_italianRegions.RData")

dati <- dati %>% arrange(WW, denominazione_regione)
# Input model preparation --------------------------------------------------------

Y <- as.integer(dati$NP)
N <- as.integer(length(Y))
Nreg <- as.integer(length(unique(dati$denominazione_regione)))
timeIdx <- as.integer(dati$WW) + 1
regIdx <- as.integer(droplevels(dati$denominazione_regione))
t <- unique(dati$WW) + 1
Ntimes <- as.integer(length(t))

loffScale <- log(1000)
lOff <- rep(logE+loffScale, times = Ntimes)
rGuesses <- dati %>% group_by(denominazione_regione) %>% summarise(sNP=sum(NP)) %>% mutate(rGuess=sNP/exp(logE+loffScale)) %$% rGuess

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
  "lOff" = lOff,
  "ttilde" = Ntimes/2
)

# Model fitting -----------------------------------------------------------
# Chains
n_chains <- 2
M <- 10000
n_cores <- mc.cores - 2

init <- function(chain_id = 1)
{
  list(# Coefficients
    logr = rnorm(1, log(mean(rGuesses)), 1),
    logh = rnorm(1, 0, 1),
    c = rnorm(1, max(dati$WW+1)/2, 1),
    logs = rnorm(1, 0, 1),
    # AR
    rho = runif(1, .2, .8),
    # CAR
    alpha = runif(1, .3, .7),
    # Random effects
    sigma = abs(rnorm(1, 0, .01)))
}

# Fit
fit_Stan <- sampling(stan_model, data = dat1, chains = n_chains, iter = M, cores = n_cores, init = init,
                     control = list(max_treedepth=8, adapt_delta = 0.9)) 

#summary(fit_Stan, pars = c("logbase", "logr", "logh", "c", "logs", "alpha", "rho"))
#save(fit_Stan, file = "WS/COVID_omicronfive.RData")
