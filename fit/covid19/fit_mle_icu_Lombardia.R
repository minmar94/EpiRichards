
# Packages ----------------------------------------------------------------

packs <- c("tidyverse", "magrittr")
lapply(packs, require, character.only = T)
rm(list = ls())

Sys.setlocale(locale = "C") 

# Load richards functions
source("richards_mle/Fun_DRichFit_Off.R")

# Data preparation --------------------------------------------------------

dati_reg <- data.table::fread("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv") %>% as_tibble() 

dati_reg <- dati_reg %>% 
  group_by(denominazione_regione) %>% 
  mutate(data = as.Date(data),
         isMonday = ifelse(lubridate::wday(data, label = T) == "Mon", 1, 0)) %>% 
  filter(data >= as.Date("2021-02-07"), data <= as.Date("2021-07-04"), denominazione_regione == "Lombardia") %>% 
  mutate(ti = as.numeric(unclass(factor(data)))) %>%
  rename(ti_orig = data, region = stato) %>%
  ungroup()

# Fit Richards MLE ----------------------------------------------------------
fam <- "Negative Binomial" # or Poisson
horizon <- 0 # forecast horizon
ti_orig_out <- dati_reg$ti_orig %>% unique() # original dates

# Covariates
X <- cbind(rep(1, nrow(dati_reg)), dati_reg$isMonday)

# Name of the variable to model
varname <- "ingressi_terapia_intensiva"
whidx <- which(colnames(dati_reg)==varname)

# observed data
pc_allobs <- dati_reg[, whidx, drop = T] 
pc_fit <- pc_allobs
# observed times
ti_fit <- dati_reg %$% ti 
ti_orig_fit <- dati_reg %$% ti_orig 

# forecast horizons
ti_orig_hor <- seq(min(ti_orig_fit), max(ti_orig_fit) + horizon, 1) # daily
mti <- max(ti_fit)
timax <- mti+horizon

# Covariates
X_fit <- as.matrix(X[ti_fit,])

# Fit
fit <- tryCatch(
  growthGLM(di = pc_fit, ti = ti_fit, alpha = .05, family = fam, tPred = timax, X = X_fit,
            maxiter = 1000, runs = 2000, nBoot = 5000, off = 18.10),
  error = function(e){
    return("Error")
  }
)

# Sistemo l'output
fitTib <- tibble(ti_orig=ti_orig_fit, pc=pc_fit)
horTib <- tibble(ti_orig=ti_orig_hor, ly=fit$hlowdiff, yhat=c(fit$linPredDiff), uy=fit$hupdiff)

# mi faccio ritornare una tibble con tutte le date, gli osservati (pc, pc_all), lower bound (ly), stima (y), upper bound (uy)
reg_out <- fitTib %>% 
  full_join(dati_reg %>% dplyr::select(ti_orig,  varname), by="ti_orig") %>% 
  full_join(horTib, by="ti_orig") %>% 
  set_colnames(value = c("x1", "pc", "pc_all", "ly", "y", "uy")) 
