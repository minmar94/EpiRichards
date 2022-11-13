
# Packages ----------------------------------------------------------------

packs <- c("tidyverse", "magrittr")
lapply(packs, require, character.only = T)
rm(list = ls())

# Load richards functions
source("richards_mle/Fun_DRichFit_Off.R")

# Load data
load("data/monkey_data.RData")

# Data preparation --------------------------------------------------------

monkey_data_formodel <- monkey_data %>% 
  group_by(location, week = cut.Date(date, breaks = "1 week", labels = FALSE)) %>% 
  summarise(NP = sum(new_cases), date = last(date), ndays = n()) %>% 
  filter(ndays == 7)  %>%
  mutate(ti = as.numeric(unclass(factor(date)))) %>%
  rename(ti_orig = date, region = location) %>%
  ungroup()

# Fit Richards MLE ----------------------------------------------------------
fam <- "Negative Binomial" # or Poisson
horizon <- 3 # forecast horizon
ti_orig_out <- monkey_data_formodel$ti_orig %>% unique() # original dates

# Covariates
X <- cbind(rep(0, nrow(monkey_data_formodel)))

# Name of the variable to model
varname <- "NP"
whidx <- which(colnames(monkey_data_formodel)==varname)

off <- c(60000000, 332403650)/1e06 # OFFSET

richout <- parsout <- list()

for(ctry in 1:length(unique(monkey_data_formodel$region))){
  
  dmodel <- monkey_data_formodel %>% filter(region == unique(monkey_data_formodel$region)[ctry])
  
  pc_allobs <- dmodel[, whidx, drop = T] # observed data
  pc_fit <- pc_allobs
  # observed times
  ti_fit <- dmodel %$% ti 
  ti_orig_fit <- dmodel %$% ti_orig
  
  # forecast horizons
  ti_orig_hor <- seq(min(ti_orig_fit), max(ti_orig_fit) + 7*horizon, 7) # weekly
  mti <- max(ti_fit)
  timax <- mti+horizon
  
  # Covariates
  X_fit <- as.matrix(X[ti_fit,])
  
  # Fit
  fit <- tryCatch(
    growthGLM(di = pc_fit, ti = ti_fit, alpha = .05, family = fam, 
              tPred = timax, X = X_fit, maxiter = 1000, runs = 2000, nBoot = 5000,
              off = off[ctry]),
    error = function(e){
      return("Error")
    }
  )
  
  # Output
  fitTib <- tibble(ti_orig=ti_orig_fit, pc=pc_fit)
  horTib <- tibble(ti_orig=ti_orig_hor, ly=fit$hlowdiff, yhat=c(fit$linPredDiff), uy=fit$hupdiff)
  
  
  richout[[ctry]] <- fitTib %>% 
    full_join(dmodel %>% dplyr::select(ti_orig,  varname), by="ti_orig") %>% 
    full_join(horTib, by="ti_orig") %>% 
    set_colnames(value = c("x1", "pc", "pc_all", "ly", "y", "uy"))
  
  
  # Richards parameters
  parsout[[ctry]] <- list()
  parsout[[ctry]][[1]] <- fit$pars
  parsout[[ctry]][[2]] <- fit$hubasyVar
  
}
names(richout) <- names(parsout) <- unique(monkey_data_formodel$region)

#save(richout, parsout, monkey_data_formodel, file = "WS/Monkeypox_fit.RData")
