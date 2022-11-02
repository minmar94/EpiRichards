# Packages ----------------------------------------------------------------

library(MASS)
library(rmutil)
library(numDeriv)
require(tidyverse)
require(magrittr)
source("richards_mle/Funs_DRichNew.R")

# Functions - Poisson -----------------------------------------------------

# Likelihood
likPois <- function(pars, ti, di, X, off) 
{
  npars <- length(pars)
  ncovs <- ncol(X)
  richpars <- pars[1:(npars-ncovs-1)]
  betapars <- pars[(npars-ncovs):(npars-1)]
  logbase <- pars[npars]
  
  linPred <- off * (1e-19+exp(logbase)+exp(X%*%betapars + ldiffRich(richpars, ti)))
  
  out <- sum(dpois(di, linPred, log=T))
  
  return(out)
}

# Gradient Richard-Poisson
PoisRichGradient <- function(pars, ti, di, X, off)
{
  # Parameters
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- (npars-ncovs-1)
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-1)]
  logbase <- pars[npars]
  loff <- log(off)
  
  # Linear predictor
  lambdat <- exp(ldiffRich(richpars, ti))
  llambdat <- log(lambdat)
  Xb <- X%*%betapars
  eXb <- exp(Xb)
  mut <- off*(1e-19+exp(logbase)+exp(Xb+llambdat))
  lmut <- log(mut)
  
  # First derivatives
  d1lt <- d1diffRich(pars, ti)
  d1eXb <- c(eXb)*X
  d1logb <- rep(exp(logbase), length(ti))
  d1mut <- off*cbind(c(eXb)*d1lt, c(lambdat)*d1eXb, d1logb)
  
  # Auxiliary
  aux1 <- exp(log(di)-lmut)
  
  # Gradient computation
  out <- colSums(c(aux1-1)*d1mut)
  
  return(out)
}

# Hessian Richard-Poisson
PoisRichHessian <- function(pars, ti, di, X, off)
{
  # Parameters
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- (npars-ncovs-1)
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-1)]
  logbase <- pars[npars]
  loff <- log(off)
  
  # Linear predictor
  lambdat <- exp(ldiffRich(richpars, ti))
  llambdat <- log(lambdat)
  Xb <- X%*%betapars
  eXb <- exp(Xb)
  mut <- off*(1e-19+exp(logbase)+exp(Xb+llambdat))
  lmut <- log(mut)
  
  # First derivatives
  d1lt <- d1diffRich(pars, ti)
  d1eXb <- c(eXb)*X
  d1logb <- rep(exp(logbase), length(ti))
  
  # Second derivatives
  d2lt <- d2diffRich(richpars, ti)
  d2eXb <- list()
  for(i in 1:ncovs)
  {
    d2eXb[[i]] <- d1eXb[,i]*as.matrix(X[,i:ncovs])
  }
  d2logb <- rep(exp(logbase), length(ti))
  
  # Auxiliary
  aux1 <- (di-mut)/mut
  aux2 <- exp(log(di)-2*lmut)
  
  jumped <- c(0,1,3,6,10,15,21,28)
  out <- matrix(NA, npars, npars)
  
  # Hessian computation
  for (i in 1:npars)
  {
    if (i<=nrich)
    {
      out[i, i:npars] <- c(colSums(c(aux1*eXb*off)*d2lt[[i]]-c(aux2*eXb^2*off^2)*d1lt[,i]*d1lt[,i:nrich]), 
                           colSums(c(aux1*off)*d1lt[,i]*d1eXb-c(aux2*exp(llambdat+Xb)*off^2)*d1lt[,i]*d1eXb),
                           sum(-c(aux2)*(c(eXb)*off*d1lt[,i])*off*d1logb))
      out[i:npars, i] <- out[i, i:npars]
    }
    
    if (i>nrich & i<npars)
    {
      k <- i-nrich
      
      out[i, i:npars] <- c(colSums(c(aux1*lambdat*off)*d2eXb[[k]]-c(aux2*lambdat^2*off^2)*d1eXb[,k]*d1eXb[,k:ncovs]),
                           sum(-c(aux2)*(lambdat*d1eXb[,k]*off)*off*d1logb))
      out[i:npars, i] <- out[i, i:npars]
    }
    if (i==npars)
    {
      out[npars, npars] <- sum(c(aux1)*off*d2logb-c(aux2)*off^2*d1logb^2)
    }
  }
  
  return(out)
}

