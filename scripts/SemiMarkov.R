### RUN THE HIDDEN SEMI MARKOV MODEL.

#install.packages("plyr")
#install.packages("rjags")
library(plyr)
library(rjags)

source("dataPreparationSemiMarkov.R")


# Simulation option - assume heterogeneous exposure to bites?
rel_eir_on<-FALSE

#------------------------------------------------------------------------------
# THE MODEL.

## TWO MODEL VERSIONS - WITH AND WITHOUT heterogeneous biting.
# Specify the model in JAGS language, but save it as a string in R:
if(rel_eir_on) {
  modelString = "
# JAGS model specification begins ...
  model {
  
  #likelihood
  for (i in 1:N) { #with right and left censoring
  times[i] ~ dsum(timeToInf[i], transTime[i])
  inInterval[i] ~ dinterval(times[i], intervalLimits[i,])
  }
  
  #Some parameters
  rho <- 0.85
  age0 <- 2920 #days
  #sigma2 <- 1.67
  bh	=	0.590076		#{parameters for the immunity functions}
  bmin	=	0.5
  db	=	3650
  kb	=	2.15506
  ub	=	7.19919
  x_I = 0.4999674  ### x_I only varies by age if widths of age groups vary. Have kept constant
  IB0	= 43.8787
  
  
  #unobserved time intervals et alia
  
  for (i in 1:N) {
  median[i] <- durP[drug[i], site[i]] + 3.5 #days until patent
  lambdaTranstime[i] <- r[drug[i]]/median[i] # ok here median should be the mean.. (gamma distribution)
  transTime[i] ~ dgamma(r[drug[i]], lambdaTranstime[i]) 
  timeToInf[i] ~ dexp(indFoi[i])
  rel_eir[i] ~ dlnorm(-sigma2/2,(1/sqrt(sigma2))^2)  ## mean of log, PRECISION in BUGS, NOT SD. precision=(1/sigma)^2
  #EIR experienced by person i now.
  indEirNow[i] = rel_eir[i]*siteEir[site[i]]*(1-rho*exp(-age[i]*365/age0)) #actual eir acting on an individual of certain age
  # need to iteratively calculate current IB, summing over past exposure.
  # EIR experienced by person i in first age group
  indEirAges[i,1] = rel_eir[i]*siteEir[site[i]]*(1-rho*exp(-(1/2)*365/age0))  ## what EIR did they experience when in age group 1?
  #indEirAges[i,1] = siteEir[site[i]]*(1-rho*exp(-(1/2)*365/age0))  ## what EIR did they experience when in age group 1?
  # IB immunity of the person in first age group.
  indIBEqij[i,1]= (indEirAges[i,1]/(indEirAges[i,1]*ub+1) *x_I) / (1+x_I/db)
  
  for(j in 2:agecat[i]) {   ## compute just up to the age group of current person i.
  # EIR experienced by person i in earlier age groups.
  indEirAges[i,j] = rel_eir[i]*siteEir[site[i]]*(1-rho*exp(-(j/2)*365/age0)) 
  indIBEqij[i,j] = (indIBEqij[i, j-1] + indEirAges[i,j]/(indEirAges[i,j]*ub + 1)*x_I)/(1+x_I/db)
  }
  b[i] = bh*((1-bmin)/(1+(indIBEqij[i,agecat[i]]/IB0)^kb) + bmin)
  indFoi[i] <- (1/365)*indEirNow[i]*b[i] #Eir2Foi
  }
  
  
  #Priors
  for(i in 1:nSites) { #drugs or sites, respectively
  shapeSiteEir[i] <-  1/(0.8)^2   #1/(coeff.var)^2  ## used in previous iteration
  rateSiteEir[i] <- shapeSiteEir[i]/priorSiteEir[i]
  siteEir[i] ~ dgamma(shapeSiteEir[i], rateSiteEir[i]) #dgamma(1E-3,1E-3)##sites
  durP[1,i] ~ dgamma(siteDurP1HyperR, siteDurP1HyperLambda) #durP1 # dgamma(1E-2,1E-2) #
  durP[2,i] ~ dgamma(siteDurP2HyperR, siteDurP2HyperLambda)
  }
  
  siteDurP1HyperR ~ dgamma(1E-3,1E-3)  #1/(0.420927)^2 #1/(coeff.var)^2
  siteDurP1HyperMean ~ dgamma(1E-3,1E-3)
  siteDurP1HyperLambda <-  siteDurP1HyperR/siteDurP1HyperMean#14.81382
  siteDurP1HyperVar <- siteDurP1HyperR/(siteDurP1HyperLambda)^2
  siteDurP2HyperR ~ dgamma(1E-3,1E-3)  #1/(0.420927)^2 #1/(coeff.var)^2
  siteDurP2HyperMean ~ dgamma(1E-3,1E-3)
  siteDurP2HyperLambda <-  siteDurP2HyperR/siteDurP2HyperMean #14.81382
  siteDurP2HyperVar <- siteDurP2HyperR/(siteDurP2HyperLambda)^2
  sigma2 ~ dgamma(10000,5682)   ### extremely restrictive prior otherwise it wanders off. shape=10,000, rate=1/scale, scale=1.76/10000
  
  #free estimates: r[1] 157.17 r[2] 16.042 -> mean=86.606 , sd=99.79257 , coeffVar=
  r4r <- 1/(1.152259)^2 #1/(coeff.var)^2
  lambda4r <- r4r/86.606
  r[1] ~ dgamma(r4r, lambda4r) #dgamma(.001,.001)#
  r[2] ~ dgamma(r4r, lambda4r) #dgamma(.001,.001)#
  }
  # ... JAGS model specification ends.
  " # close quote to end modelString
  
}

if(!rel_eir_on) {
  modelString = "
# JAGS model specification begins ...
  model {
  
  #likelihood
  for (i in 1:N) { #with right and left censoring
  times[i] ~ dsum(timeToInf[i], transTime[i])
  inInterval[i] ~ dinterval(times[i], intervalLimits[i,])
  }
  
  #Some parameters
  rho <- 0.85
  age0 <- 2920 #days
  bh	=	0.590076		#{parameters for the immunity functions}
  bmin	=	0.5
  db	=	3650
  kb	=	2.15506
  ub	=	7.19919
  x_I = 0.4999674  ### 
  IB0	= 43.8787
  
  
  #unobserved time intervals et alia
  
  for (i in 1:N) {
  median[i] <- durP[drug[i], site[i]] + 3.5 #days until patent
  lambdaTranstime[i] <- r[drug[i]]/median[i] # ok here median should be the mean.. (gamma distribution)
  transTime[i] ~ dgamma(r[drug[i]], lambdaTranstime[i]) 
  timeToInf[i] ~ dexp(indFoi[i])
  #EIR experienced by person i now.
  indEirNow[i] = siteEir[site[i]]*(1-rho*exp(-age[i]*365/age0)) #actual eir acting on an individual of certain age
  # need to iteratively calculate current IB, summing over past exposure.
  # EIR experienced by person i in first age group
  indEirAges[i,1] = siteEir[site[i]]*(1-rho*exp(-(1/2)*365/age0))  ## what EIR did they experience when in age group 1?
  # IB immunity of the person in first age group.
  indIBEqij[i,1]= (indEirAges[i,1]/(indEirAges[i,1]*ub+1) *x_I) / (1+x_I/db)
  
  for(j in 2:agecat[i]) {   ## compute just up to the age group of current person i.
  # EIR experienced by person i in earlier age groups.
  indEirAges[i,j] = siteEir[site[i]]*(1-rho*exp(-(j/2)*365/age0)) # j/2=current age, at end of age gp.
  indIBEqij[i,j] = (indIBEqij[i, j-1] + indEirAges[i,j]/(indEirAges[i,j]*ub + 1)*x_I)/(1+x_I/db)
  }
  b[i] = bh*((1-bmin)/(1+(indIBEqij[i,agecat[i]]/IB0)^kb) + bmin)
  indFoi[i] <- (1/365)*indEirNow[i]*b[i] #Eir2Foi
  }
  
  
  #Priors
 
  for(i in 1:nSites) { #drugs or sites, respectively
  shapeSiteEir[i] <-  1/(0.8)^2   #1/(coeff.var)^2  ## used in previous iteration
  rateSiteEir[i] <- shapeSiteEir[i]/priorSiteEir[i]
  siteEir[i] ~ dgamma(shapeSiteEir[i], rateSiteEir[i]) #dgamma(1E-3,1E-3)##sites
  durP[1,i] ~ dgamma(siteDurP1HyperR, siteDurP1HyperLambda) #durP1 # dgamma(1E-2,1E-2) #
  durP[2,i] ~ dgamma(siteDurP2HyperR, siteDurP2HyperLambda)
  }
  
  siteDurP1HyperR ~ dgamma(1E-3,1E-3)  #1/(0.420927)^2 #1/(coeff.var)^2
  siteDurP1HyperMean ~ dgamma(1E-3,1E-3)
  siteDurP1HyperLambda <-  siteDurP1HyperR/siteDurP1HyperMean#14.81382
  siteDurP1HyperVar <- siteDurP1HyperR/(siteDurP1HyperLambda)^2
  siteDurP2HyperR ~ dgamma(1E-3,1E-3)  #1/(0.420927)^2 #1/(coeff.var)^2
  siteDurP2HyperMean ~ dgamma(1E-3,1E-3)
  siteDurP2HyperLambda <-  siteDurP2HyperR/siteDurP2HyperMean #14.81382
  siteDurP2HyperVar <- siteDurP2HyperR/(siteDurP2HyperLambda)^2
  
  #free estimates: r[1] 157.17 r[2] 16.042 -> mean=86.606 , sd=99.79257 , coeffVar=
  r4r <- 1/(1.152259)^2 #1/(coeff.var)^2
  lambda4r <- r4r/86.606
  r[1] ~ dgamma(r4r, lambda4r) #dgamma(.001,.001)#
  r[2] ~ dgamma(r4r, lambda4r) #dgamma(.001,.001)#
  }
  # ... JAGS model specification ends.
  " # close quote to end modelString
}


## precompute some parameters:
# infection blocking immunity
death_rate<- 1/(21*365)
den_age<-vector(length=69*2)  # max age in dataset is 69
age_rate<-1/(0.5*365)  ## 6 monthly age groups
den_age[1]<- 1/(1+age_rate / death_rate)
x_I<- den_age[1]/death_rate  # is constant if age_width constant so ok for now.
for(i in 2:length(den_age)) {
  den_age[i] <- age_rate* den_age[i-1]/(age_rate + death_rate)
}
d$agecat<-as.numeric(cut(d$ageyears,breaks=seq(from=0,to=69,by=0.5),include.lowest=T) )

#------------------------------------------------------------------------------
# THE DATA.
# Specify the data in a form that is compatible with model, as a list:
dataList = list(times = d$times, N = nrow(d), nSites=length(unique(d$siteNum)), 
  site=d$siteNum, drug=d$treatNum, 
   priorSiteEir=siteEIRs, inInterval=d$inInterval, age=d$ageyears, agecat=d$agecat,
  intervalLimits=cbind(d$lastNegTime, d$maxObsTime))



#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

rInit = c(15,15) #15 #
eirInit   = rep(15,length(unique(d$site)))
transTimeInit <- d$lastNegTime  #for the left censored ones (now just for all)
timeToInfInit  <- (d$maxObsTime-transTimeInit) - 0.5 #for the left censored ones (now just for all)

transTimeInit[d$inInterval==2] <- d$maxObsTime[d$inInterval==2]
timeToInfInit[d$inInterval==2] <- 10


initsList = list(siteEir=eirInit, transTime=transTimeInit, r = rInit, # , siteFoi=foiInit, 
	timeToInf=timeToInfInit)

#------------------------------------------------------------------------------
# RUN THE CHAINS

if(rel_eir_on) parameters = c("siteEir", "r", "siteDurP1HyperR", "siteDurP1HyperMean", "siteDurP1HyperLambda", "siteDurP1HyperVar",
 "siteDurP2HyperR", "siteDurP2HyperMean", "siteDurP2HyperLambda", "siteDurP2HyperVar", "durP", 
  "sigma2", "rel_eir")
if(!rel_eir_on) parameters = c("siteEir", "r", "siteDurP1HyperR", "siteDurP1HyperMean", "siteDurP1HyperLambda", "siteDurP1HyperVar",
                              "siteDurP2HyperR", "siteDurP2HyperMean", "siteDurP2HyperLambda", "siteDurP2HyperVar", "durP")
                              
quickRun<-F   ## set to TRUE to do a quick run to see if the model is working.
if(!quickRun) {
  adaptSteps = 4000 
  burnInSteps = 4000
  nChains = 2
  numSavedSteps=100000
  thinSteps=30
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
} else {
  adaptSteps = 2000 
  burnInSteps = 2000
  nChains = 2
  numSavedSteps=20 
  thinSteps=25 
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
}

t1<-Sys.time()
# Create, initialize, and adapt the model:
jagsModel = jags.model(textConnection(modelString), data=dataList, inits=initsList, 
                        n.chains=nChains, n.adapt=adaptSteps)


# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update(jagsModel, n.iter=burnInSteps )



# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel, variable.names=parameters, 
                            n.iter=nIter, thin=thinSteps)

# to save the posteriors: 
#save(codaSamples, file="codaSemiMarkov")

t2<-Sys.time()
t2-t1