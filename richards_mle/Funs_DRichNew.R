# Packages ----------------------------------------------------------------

library(MASS)
require(tidyverse)
require(magrittr)


# Functions - Richards -----------------------------------------------------

# Linear Predictor
lRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # Denominator
  l1pehptInf <- log(1+exp(h*(p-ti)))
  l1pehpt <- ifelse(is.infinite(l1pehptInf), h*(p-ti), l1pehptInf)
  logden <- -s*l1pehpt
  
  # Out
  lout <- logden
  
  return(lout)
}

# Linear Predictor monotone
ldiffRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # All times
  tAll <- c(ti[1]-1, ti)
  
  # Auxiliary
  l1pehptInf <- log(1+exp(h*(p-tAll)))
  l1pehpt <- ifelse(is.infinite(l1pehptInf), h*(p-tAll), l1pehptInf)
  
  # Terms
  # 1
  laux1 <- log(diff(exp( -s * l1pehpt)))
  
  # Out  
  lout <- laux1
  
  return(lout)
}

# Derivate prime
d1diffRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # All times
  tAll <- c(ti[1]-1, ti)
  
  # Auxiliary
  l1pehptInf <- log(1+exp(h*(p-tAll)))
  l1pehpt <- ifelse(is.infinite(l1pehptInf), h*(p-tAll), l1pehptInf)
  lehpt <- h*(p-tAll)
  
  # Terms
  #1
  laux1 <- log(diff(exp( -s * l1pehpt)))
  
  # First derivatives
  
  # h
  dlh <- -s * h *               # Jacobian
    diff((p-tAll) * exp(lehpt - (s+1) * l1pehpt))
  
  # p
  dlp <- -s * h *
    diff(exp(lehpt - (s+1)*l1pehpt))
  
  # s  
  dls <- -diff(exp(log(l1pehpt) - s*l1pehpt))
  
  out <- matrix(c(dlh, dlp, dls), nrow=length(ti), ncol=3)
  
  return(out)
}

# Derivate seconde
d2diffRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # All times
  tAll <- c(ti[1]-1, ti)
  
  # First derivatives
  d1lt <- d1diffRich(ti = ti, pars = pars)
  dlh <- d1lt[, 1]
  dlp <- d1lt[, 2]
  dls <- d1lt[, 3]
  
  # Auxiliary
  l1pehptInf <- log(1+exp(h*(p-tAll)))
  l1pehpt <- ifelse(is.infinite(l1pehptInf), h*(p-tAll), l1pehptInf)
  lehpt <- h*(p-tAll)
  
  # Terms
  #1
  laux1 <- log(diff(exp( -s * l1pehpt)))
  
  # Second derivatives
  # h
  dlhh <- dlh - s* h^2 *                      # Jacobian
    diff( exp(log((p-tAll)^2) + lehpt + (-s-1)*l1pehpt) * 
            (1 - exp(log(s+1) - l1pehpt + lehpt)) )
  dlhp <- -s * h *                      # Jacobian
    diff( exp(lehpt + (-s-1)*l1pehpt) * 
            (1 + (p-tAll)*h*
               (1 - exp(log(s+1) - l1pehpt + lehpt))) )
  dlhs <- - h *                       # Jacobian
    diff( (p-tAll)*exp(lehpt + (-s-1) * l1pehpt) * (1-exp(log(s)+log(l1pehpt))))
  
  dlpp <- - s * h^2 *
    diff( exp((-s-1)*l1pehpt + lehpt) * (1-exp(log(s+1)+lehpt-l1pehpt)))
  dlps <- - h *
    diff( exp(lehpt + (-s-1) * l1pehpt) * (1-exp(log(s)+log(l1pehpt))))
  
  dlss <- diff(exp(2*log(l1pehpt) - s * l1pehpt))
  
  # Output
  out <- list(matrix(c(dlhh, dlhp, dlhs), nrow=length(ti), ncol=3),
              matrix(c(dlpp, dlps), nrow=length(ti), ncol=2),
              matrix(c(dlss), nrow=length(ti), ncol=1))
  
  return(out)
}