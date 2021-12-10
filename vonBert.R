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
m1 <-list(dl=dat$dl,line=dat$stocknum,L=dat$lt,t=dat$t,id=dat$fid) #data string

# individual and group von Bertalanffy growth model 
cat(
  "model{
    # von Bertalanffy growth model (Linf from Kraeuter et al 2007 JSR)
    for (i in 1:1937) {
    
    dl[i] ~ dnorm(mu[i],tau)
    mu[i] <- (130 - L[i])*(1-exp(-K[id[i]]*t[i]))
    
    }
    for (ind in 1:990) {
      K[ind] ~ dnorm(K_line[line[ind]],tau_line[line[ind]])
    }
    
    for (l in 1:6) {
      K_line[l] ~ dnorm(0,0.1)
      tau_line[l] ~ dgamma(0.1,0.1)
      sigma_line[l] <- sqrt(1/tau_line[l])
    }
    tau ~ dgamma(0.1,0.1)
  }",
  file="m1.jag"
)

m1.inits <- list(list("K"=numeric(990),"K_line"=numeric(6),"tau"=0.2,"tau_line"=c(0.1,0.1,0.1,0.1,0.1,0.1)),
                 list("K"=numeric(990),"K_line"=numeric(6),"tau"=0.01,"tau_line"=c(0.1,0.1,0.1,0.1,0.1,0.1)),
                 list("K"=numeric(990),"K_line"=numeric(6),"tau"=0.1,"tau_line"=c(0.1,0.1,0.1,0.1,0.1,0.1)))

parameters <- c("sigma_line","K_line")

mfit <- jags(data = m1,
               inits = m1.inits,
               parameters.to.save = parameters,
               model.file = "m1.jag",
               n.chains = 3,
               n.iter = 10000,
               n.burnin = 1000,
               n.thin = 3) 

# Sort output for plot
out <- data.frame(o1=mfit$BUGSoutput$sims.list$sigma_line)
mod1 <- stack(out,select=c('o1.1'))
names(mod1) <- c("K","line")
mod2 <- stack(out,select=c('o1.2'))
names(mod2) <- c("K","line")
mod3 <- stack(out,select=c('o1.3'))
names(mod3) <- c("K","line")
mod4 <- stack(out,select=c('o1.4'))
names(mod4) <- c("K","line")
mod5 <- stack(out,select=c('o1.5'))
names(mod5) <- c("K","line")
mod6 <- stack(out,select=c('o1.6'))
names(mod6) <- c("K","line")

# K plot
p1 <- ggplot(mod1,aes(x=K))+
  geom_density(color="darkblue",fill="lightblue",alpha=0.4)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(0.025,0.045)+
  labs(title=expression("Narrow River"),x = "", y = "")+
  geom_vline(xintercept = 0.03519011, linetype="dotted", 
             color = "black", size=1)
p2 <- ggplot(mod2,aes(x=K))+
  geom_density(color="darkblue",fill="lightblue",alpha=0.4)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(0.025,0.045)+
  labs(title=expression("Fisher's Island"),x = "", y = "")+
  geom_vline(xintercept = 0.03535833, linetype="dotted", 
             color = "black", size=1)
p3 <- ggplot(mod3,aes(x=K))+
  geom_density(color="darkblue",fill="lightblue",alpha=0.4)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(0.025,0.045)+
  labs(title=expression("Green Hill"),x = "", y = "")+
  geom_vline(xintercept = 0.03515213, linetype="dotted", 
             color = "black", size=1)
p4 <- ggplot(mod4,aes(x=K))+
  geom_density(color="darkblue",fill="lightblue",alpha=0.4)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(0.025,0.045)+
  labs(title=expression("NEH"),x = "", y = "")+
  geom_vline(xintercept = 0.03547868, linetype="dotted", 
             color = "black", size=1)
p5 <- ggplot(mod5,aes(x=K))+
  geom_density(color="darkblue",fill="lightblue",alpha=0.4)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(0.025,0.045)+
  labs(title=expression("Martha's Vinyard"),x = "", y = "")+
  geom_vline(xintercept = 0.03537808, linetype="dotted", 
             color = "black", size=1)
p6 <- ggplot(mod6,aes(x=K))+
  geom_density(color="darkblue",fill="lightblue",alpha=0.4)+
  theme_classic()+
  theme(legend.position="none")+
  xlim(0.025,0.045)+
  labs(title=expression("Connecticut"),x = "Variance in estimated Brody growth coefficients", y = "")+
  geom_vline(xintercept = 0.03589763, linetype="dotted", 
             color = "black", size=1)
p1/p2/p3/p4/p5/p6