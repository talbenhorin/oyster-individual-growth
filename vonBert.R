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
#library(ggdistribute)

dat <- read.csv("dlvb.csv", fill = FALSE, header = TRUE)
m1 <-list(dl=dat$dl,line=dat$lineNum,L=dat$d2) #data string

# individual and group von Bertalanffy growth model 
cat(
  "model{
    # von Bertalanffy growth model (Linf from Kraeuter et al 2007 JSR)
    for (i in 1:981) {
    
    dl[i] ~ dnorm(mu[i],tau)
    mu[i] <- (130 - L[i])*(1-exp(-K[i]*30))
    K[i] ~ dnorm(K_line[line[i]],tau_line[line[i]])
    }
    for (l in 1:6) {
      K_line[l] ~ dnorm(0,0.1)
      tau_line[l] ~ dgamma(0.1,0.1)
    }
    tau ~ dgamma(0.1,0.1)
  }",
  file="m1.jag"
)

m1.inits <- list(list("K"=numeric(981),"K_line"=numeric(6),"tau"=0.2,"tau_line"=c(0.1,0.1,0.1,0.1,0.1,0.1)),
                 list("K"=numeric(981),"K_line"=numeric(6),"tau"=0.01,"tau_line"=c(0.1,0.1,0.1,0.1,0.1,0.1)),
                 list("K"=numeric(981),"K_line"=numeric(6),"tau"=0.1,"tau_line"=c(0.1,0.1,0.1,0.1,0.1,0.1)))

parameters <- c("K_line","tau_line","tau")

mfit <- jags(data = m1,
               inits = m1.inits,
               parameters.to.save = parameters,
               model.file = "m1.jag",
               n.chains = 3,
               n.iter = 4000,
               n.burnin = 1000,
               n.thin = 3) 