# Programa Monitora
# Data analysis
# Fit multi-species N-mixture model
# tss (time for space substitution)
# multi-region, multi-species model
# Elildo Carvalho Jr @ICMBio/CENAP

rm(list = ls())

#----- load packages
library(here)
library(tidyverse)
library(jagsUI)


# load data
mydata <- readRDS(here("data", "data_long.rds"))
data_wide <- readRDS(here("data", "data_wide.rds"))
y_wide <- data_wide$y
rm(data_wide)

# create indicator w[pa,site] that tells if site was ever sampled in pa
#w <- apply(y_wide, c(5,1), max, na.rm=TRUE)
#w[w == Inf | w == -Inf ] <- 0
#w[w != 0] <- 1
#head(w)


# bundle data for jagsUI
jags_data <- list(y = mydata$y,
                  pa = dense_rank(mydata$pa),
                  #sp = dense_rank(mydata$sp),
                  pop = mydata$pop,
                  site = mydata$site,
                  #site_pa = mydata$indexes$nsite,
                  year = mydata$year-2013,
                  #rep = mydata$rep,
                  #w = w,
                  
                  nobs = mydata$nobs,
                  npa = mydata$npa,
                  #nsp = mydata$nsp,
                  npop = mydata$npop,
                  nsite = mydata$nsite,
                  nyear = mydata$nyear
                  #nrep = mydata$nrep
                  )
rm(mydata) # save space
str(jags_data)

# write model in jags
sink(here("tss", "tss.jags"))
cat("
model {

##----- global-level priors

# abundance
mu.a0.global ~ dnorm(0,0.01)
mu.a1.global ~ dnorm(0,0.01)
sigma.a0.global ~ dunif(0,10)
sigma.a1.global ~ dunif(0,10)
tau.a0.global <- pow(sigma.a0.global,-2)
tau.a1.global <- pow(sigma.a1.global,-2)


# detection
mu.b0.global ~ dnorm(0, 0.01)
sigma.b0.global ~ dunif(0,10)
tau.b0.global <- pow(sigma.b0.global,-2)

#----- PA-level priors
    for (pa in 1:npa) {
    
      # abundance
      mu.a0.pa[pa] ~ dnorm(mu.a0.global,tau.a0.global)
      mu.a1.pa[pa] ~ dnorm(mu.a1.global,tau.a1.global)
      sigma.a0.pa[pa] ~ dunif(0,10)
      sigma.a1.pa[pa] ~ dunif(0,10)
      tau.a0.pa[pa] <- pow(sigma.a0.pa[pa],-2)
      tau.a1.pa[pa] <- pow(sigma.a1.pa[pa],-2)

      # detection
      mu.b0.pa[pa] ~ dnorm(mu.b0.global,tau.b0.global)
      sigma.b0.pa[pa] ~ dunif(0,10)
      tau.b0.pa[pa] <- pow(sigma.b0.pa[pa],-2)

      #----- population-level priors
      for (pop in 1:npop) {
      
      # inclusion prob
      #omega[pa,pop] ~ dunif(0,1)
      
      # abundance
      a0[pa,pop] ~ dnorm(mu.a0.pa[pa],tau.a0.pa[pa])
      a1[pa,pop] ~ dnorm(mu.a1.pa[pa],tau.a1.pa[pa])
      
      # zero-inflation
      #for (site in 1:nsite) {
      #  w[pa,pop,site] ~ dbern(omega[pa,pop])
      #}#site

      # detection
      b0[pa,pop] ~ dnorm(mu.b0.pa[pa],tau.b0.pa[pa])
      
      # random yesr effects on detection
      tau.year[pa,pop] <- 1 / (sd.year[pa,pop]*sd.year[pa,pop])
      sd.year[pa,pop] ~ dunif(0, 1)

      }#pop
    }#pa


##----- likelihood

# ecological submodel
for (pa in 1:npa) {
  for(pop in 1:npop){
      #for(site in 1:nsite){
        for(year in 1:nyear){
      
        log(lambda[pa,pop,year]) <- a0[pa,pop] + a1[pa,pop]*(year[year]-1) # minus 1 to set 1st to zero
        
        for(site in 1:nsite){
          n[pa,pop,site,year] ~ dpois(lambda[pa,pop,year])
        }#site

   }#year
  }#pop
}#pa


# observation submodel
for (pa in 1:npa) {
  for(pop in 1:npop){
    for(year in 1:nyear){

    eps_year[pa,pop,year] ~ dnorm(0, tau.year[pa,pop])
    logit(p[pa,pop,year]) <- b0[pa,pop] + eps_year[pa,pop,year]
    
    }#year
  }#pop
}#pa 

# observation process
for(i in 1:nobs){
  
  # observation model
  y[i] ~ dbin( p[pa[i],pop[i],year[i]], n[pa[i],pop[i],site[i],year[i]] )

}#i


}#model",fill = TRUE)
sink()



# parameters to monitor
params <- c('lambda', 'n',  'p', 'a0', 'a1')


# initialize n at maximum count ever observed at a site
#n[pa,pop,site,year]
n_empirical <- apply(y_wide, c(5,4,1,3), max)
n_empirical[2,3,,] # check
str(n_empirical)
max_values <- apply(n_empirical, c(1,2,3), max, na.rm=TRUE)
max_values[max_values == Inf | max_values == -Inf ] <- 0
str(max_values)
max_values[2,3,]
for (i in 1:dim(n_empirical)[1]) {
  for (j in 1:dim(n_empirical)[2]) {
    for (k in 1:dim(n_empirical)[3]) {
      n_empirical[i,j,k,] <- (max_values[i,j,k]+1)*2
    }
  }
}
n_empirical[2,3,,] # check
#n_empirical <- (n_empirical*2)
n_empirical[2,3,,] # check again
n_empirical[1,1,,] # check again
n_empirical[2,4,,] # check again

#z_emprirical <- apply(y_wide, c(5,4), max, na.rm=TRUE)
#z_emprirical[z_emprirical == Inf | z_emprirical == -Inf ] <- 0
#z_emprirical[z_emprirical > 0] <- 1
#z_emprirical

inits <- function() list(n = n_empirical#,
                         #z = z_emprirical
)
inits()


# mcmc parameters
ni <- 750000
nt <- 2500
nb <- 250000
nc <- 3


# call jags

mod <- jagsUI::jags(data=jags_data,
                    inits=inits,
                    #inits=NULL,
                    parameters.to.save=params,
                    model.file=here("tss", "tss.jags"),
                    n.chains=nc,
                    n.thin=nt,
                    n.iter=ni,
                    n.burnin=nb,
                    parallel = TRUE
)

saveRDS(mod, here("results", "res_tss_multiregion.rds"))


