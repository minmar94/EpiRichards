# Packages ----------------------------------------------------------------
packs <- c("tidyverse", "magrittr", "rstan", "surveillance")
lapply(packs, require, character.only = T)
rm(list = ls())

Sys.setlocale(locale = "C") 

# Load data
load("data/measles.rda")

# Stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load richards functions
stan_model <- stan_model("richards_bayes/RichSTCar_NoCovsNew2.stan")

# Data preparation --------------------------------------------------------
dat <- measles@observed %>% as_tibble() %>% 
  mutate(year = rep(2005:2018, each = 26), biweek = rep(paste("W", 1:26, sep = "_"), 14))

# Compute adjacency matrix
adj <- poly2adjmat(measles@map)

dat <- dat %>% gather(County, Count, -year, -biweek) %>% 
  filter(year == 2015) %>% 
  mutate(off = rep(measles@map$POP2015, each = 26),
         ti = unclass(factor(biweek, levels = unique(.$biweek), ordered = T))) %>%
  arrange(County, ti)

# Offset
lOff <- log(dat$off/1e04)

rGuesses <- dat %>% group_by(County) %>% summarise(sNP=sum(Count), lOff = first(log(off/1e04))) %>% 
  mutate(rGuess=sNP/exp(lOff)) %$% rGuess

# Input model preparation --------------------------------------------------------

Y <- as.integer(dat$Count)
N <- as.integer(length(Y))
Nreg <- as.integer(length(unique(dat$County)))
timeIdx <- as.integer(dat$ti)
regIdx <- as.integer(droplevels(factor(dat$County)))
t <- unique(dat$ti)
Ntimes <- as.integer(length(t))

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
#n_cores <- mc.cores - 2

init <- function(chain_id = 1)
{
  list(# Coefficients
    logr = rnorm(1, log(mean(rGuesses)), 1),
    logh = rnorm(1, 0, 1),
    c = rnorm(1, max(dat$ti+1)/2, 1),
    logs = rnorm(1, 0, 1),
    # AR
    rho = runif(1, .2, .8),
    # CAR
    alpha = runif(1, .3, .7),
    # Random effects
    sigma = abs(rnorm(1, 0, .01)))
}


# Fit
fit_Stan <- sampling(stan_model, data = dat1, chains = n_chains, iter = M, 
                     control = list(max_treedepth = 10)) # 10 minutes

#summary(fit_Stan, pars = c("logbase", "logr", "logh", "c", "logs", "alpha", "rho"))
#save(fit_Stan, file = "WS/Measles_germany.RData")
