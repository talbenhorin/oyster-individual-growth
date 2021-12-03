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
m1 <-list(dl=dat$dl,line=dat$line,t=dat$days) #data string

# individual and group von Bertalanffy growth model 
cat(
  "model{
    for (i in 1:XX) {

  }",
  file="m1.jag"
)

mfit <- jags(data = m1,
               inits = m1.inits,
               parameters.to.save = parameters,
               model.file = "m1.jag",
               n.chains = 3,
               n.iter = 4000,
               n.burnin = 1000,
               n.thin = 3) 