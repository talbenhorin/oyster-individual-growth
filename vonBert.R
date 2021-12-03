rm(list=ls(all=TRUE))

# Download JAGS at https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
# Load these packages the first time your run R
#install.packages("RTools")
#install.packages("rjags")
#install.packages("coda")
#install.packages("R2jags")
#install.packages("hdi")
#install.packages("MCMCvis")
#install.packages("patchwork")
#install.packages("ggdistribute")

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)
library(ggplot2)
library(patchwork)
library(scales)
library(ggdistribute)

dat <- read.csv("dlvb.csv", fill = FALSE, header = TRUE)
m1 <-list(dl=dat$dl,line=dat$line,dt=dat$days,ind=dat$ind,L=dat$length) #data string

# individual and group von Bertalanffy growth model 
cat(
  "model{
    for (i in 1:XX) {
    #
    dl[i] ~ dnorm(mu[i],tau)
    mu[i] <- (Linf - L[i])*(1-exp(-K[ind[i]]*dt))
    }
    for (l in 1:XX) {
      K[l] ~ dnorm(mu_line[l],tau_line[l])
    }
    mu_line ~ dnorm(0,0.1)
    tau ~ dgamma(0.1,0.1)
    tau_line ~ dgamma(0.1,0.1)
  ",
  file="m1.jag"
)

cat(
  "model{
    for (i in 1:5320) {
      # Observation model across serial dilutions
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-MPN[i]*v[i])
      
      # Biological model for microbial abundance
      MPN[i] ~ dpois(lambda[i])
      log(lambda[i]) <- b0 + b1*gear[i] + b2*tide[i] + b3*gear[i]*tide[i] + b4*mod[i] + b5*hi[i] + U[samp[i]] +V[time[i]]
    }


mfit <- jags(data = m1,
               inits = m1.inits,
               parameters.to.save = parameters,
               model.file = "m1.jag",
               n.chains = 3,
               n.iter = 4000,
               n.burnin = 1000,
               n.thin = 3) 