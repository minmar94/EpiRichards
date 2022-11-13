# Packages ----------------------------------------------------------------

packs <- c("tidyverse", "magrittr", "surveillance")
lapply(packs, require, character.only = T)
rm(list = ls())

Sys.setlocale(locale = "C") 

# Load richards functions
source("richards_mle/Fun_DRichcFit_Off.R")

# Data preparation --------------------------------------------------------
data("campyDE")
dat_campy <- campyDE %>% dplyr::select(date, case) %>% 
  filter(lubridate::year(date) == 2011) %>% 
  mutate(ti = as.numeric(unclass(factor(date)))) %>%
  rename(ti_orig = date)

# Fit Richards MLE ----------------------------------------------------------
fam <- "Negative Binomial" # or Poisson
horizon <- 0 # forecast horizon
ti_orig_out <- dat_campy$ti_orig %>% unique() # original dates

# Covariates
X <- cbind(rep(0, nrow(dat_campy)))

# Name of the variable to model
varname <- "case"
whidx <- which(colnames(dat_campy)==varname)

# observed data
pc_allobs <- dat_campy[, whidx, drop = T] 
pc_fit <- pc_allobs
# observed times
ti_fit <- dat_campy %$% ti 
ti_orig_fit <- dat_campy %$% ti_orig 

# forecast horizons
ti_orig_hor <- seq(min(ti_orig_fit), max(ti_orig_fit) + horizon*7, 7) # weekly
mti <- max(ti_fit)
timax <- mti+horizon

# Covariates
X_fit <- as.matrix(X[ti_fit,])

# Fit
fit <- tryCatch(
  growthGLM(di = pc_fit, ti = ti_fit, alpha = .05, family = fam, tPred = timax, X = X_fit,
            maxiter = 1000, runs = 2000, nBoot = 5000, off = 80274984/1e07),
  error = function(e){
    return("Error")
  }
)

# Sistemo l'output
fitTib <- tibble(ti_orig=ti_orig_fit, pc=pc_fit)
horTib <- tibble(ti_orig=ti_orig_hor, ly=fit$hlowdiff, yhat=c(fit$linPredDiff), uy=fit$hupdiff)

# mi faccio ritornare una tibble con tutte le date, gli osservati (pc, pc_all), lower bound (ly), stima (y), upper bound (uy)
richout <- fitTib %>% 
  full_join(dat_campy %>% dplyr::select(ti_orig,  varname), by="ti_orig") %>% 
  full_join(horTib, by="ti_orig") %>% 
  set_colnames(value = c("x1", "pc", "pc_all", "ly", "y", "uy")) 

#save(richout, fit, dat_campy, file = "WS/Campylobacter_germany.RData")
