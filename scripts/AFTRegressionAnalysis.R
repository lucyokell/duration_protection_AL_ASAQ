# RUN ACCELERATED FAILURE TIME REGRESSION MODELS.
# Files required:
# Clinical trial data from WWARN.
# source("scripts/regressionDataPreparation.R")


#install.packages("lmtest")
library(lmtest)
library(survival)

# function for extracting regression results.
tab<-function(log_est,se_log,dp=2) {
  paste0(round(exp(log_est),dp)," (",round(exp(log_est-1.96*se_log),dp),", ",round(exp(log_est+1.96*se_log),dp),")")
}

myAic <- function(fit){
  print(paste("loglik:", fit$loglik))
  print(paste("df:", sum(fit$df)))
  -2*fit$loglik[2] + 2*sum(fit$df)
}

myLrt <- function (obj1, obj2) {  ## simpler model is obj1,more complex is obj2.
  L0 <- obj1$loglik[2]
  L1 <- obj2$loglik[2]
  L01 <- as.vector(- 2 * (L0 - L1))
  df <- sum(obj2$df) - sum(obj1$df)
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))
}

source("regressionDataPreparation.R") # run code which reads in data, cleans and recodes for regression.

######### REGRESSION ###########

# BIVARIATE REGRESSSION ADJUSTING ALL VARIABLES FOR EIR
s <- Surv(time=d$time1Surv, time2=d$time2Surv,  type = "interval2" )
fit<-survreg(s ~ log_estEIR + frailty.gaussian(site, sparse=F), data=d,dist="lognormal")
est<-fit$coeff[2]
se<-sqrt(fit$var[2,2])
tab(est,se,2)

# drug, overall.
fit<-survreg(s ~ treatNum + log_estEIR + frailty.gaussian(site, sparse=F), data=d,dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
# drug, 20% 86Y
# drug, when mut prevalence = 20% (rescale so that 20=0 in new variable)
d$prev86_20_zero<- (d$prev86-0.2)
temp<-subset(d,!is.na(prev86))
s3 <- Surv(time=temp$time1Surv, time2=temp$time2Surv,  type = "interval2" )
fit <- survreg(s3 ~  treatNum + log_estEIR +
                 treatNum*prev86_20_zero + frailty.gaussian(site, sparse=F), data=temp)
summary(fit)  ### 
est<-fit$coeff[2]
se<-sqrt(fit$var[2,2])
tab(est, se)  # ratio of reinfection time ASAQ : AL, when prev86=0.2
# drug, 80% 86Y. Reverse so that 80%=0, to get drug effect at 80%
temp$prev86_rev<- 0.8-temp$prev86
s3 <- Surv(time=temp$time1Surv, time2=temp$time2Surv,  type = "interval2" )
fit <- survreg(s3 ~  treatNum+ log_estEIR + 
                 treatNum*prev86_rev + frailty.gaussian(site, sparse=F), data=temp)
summary(fit)  ### 
est<-fit$coeff[2]
se<-sqrt(fit$var[2,2])
tab(est, se)  # ratio of reinfection time ASAQ : AL, when prev86=0.2
# drug effect when there is 20% pfcrt
temp<-subset(d,!is.na(prev76))
d$prev76_20_zero<- (d$prev76-0.2)
s3 <- Surv(time=temp$time1Surv, time2=temp$time2Surv,  type = "interval2" )
fit <- survreg(s3 ~  treatNum + log_estEIR  + 
                 treatNum*prev76_20_zero + frailty.gaussian(site, sparse=F), data=temp)
summary(fit)  ### 
est<-fit$coeff[2]
se<-sqrt(fit$var[2,2])
tab(est, se)
# drug, 80% 76T. Reverse so that 80%=0, to get drug effect at 80%
temp<-subset(d,!is.na(prev76))
temp$prev76_rev<- 0.8-temp$prev76
s3 <- Surv(time=temp$time1Surv, time2=temp$time2Surv,  type = "interval2" )
fit <- survreg(s3 ~  treatNum + log_estEIR+ 
                 treatNum*prev76_rev + frailty.gaussian(site, sparse=F), data=temp)
summary(fit)  ### 
est<-fit$coeff[2]
se<-sqrt(fit$var[2,2])
tab(est, se)

fit<-survreg(s ~ treatNum*d$prev86 + frailty.gaussian(site, sparse=F), data=d,dist="lognormal")
summary(fit)
tab(fit$coeff[2],sqrt(fit$var[2,2]))

# age
temp<-d[which(!is.na(d$ageyears2)),]
s2 <- Surv(time=temp$time1Surv, time2=temp$time2Surv,  type = "interval2" )
fit<-survreg(s2 ~ ageyears2 + ageyears2_sq + ageyears2_cub + log_estEIR + frailty.gaussian(site, sparse=F), 
             data=temp,dist="lognormal")
summary(fit)
tab(fit$coeff[2] ,  sqrt(fit$var[2,2]))
tab(fit$coeff[3]  , sqrt(fit$var[3,3]))
tab(fit$coeff[4], sqrt(fit$var[4,4]),4)
fit2<-survreg(s2 ~ log_estEIR + frailty.gaussian(site, sparse=F), 
              data=temp,dist="lognormal")
-2*fit$loglik[2] - (-2*fit2$loglik[2]) + 2

# gender
fit<-survreg(s ~ gender + log_estEIR + frailty.gaussian(site, sparse=F), data=d,dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
# anaemic
fit<-survreg(s ~ anaemic + log_estEIR + frailty.gaussian(site, sparse=F), data=d,dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
#spleen
fit<-survreg(s ~ spleen + log_estEIR, data=d, dist="lognormal")	# reported in paper
tab(fit$coeff[2],sqrt(fit$var[2,2]))
#fever
fit<-survreg(s ~ fev + log_estEIR + frailty.gaussian(site, sparse=F), data=d, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
#underweight
fit<-survreg(s ~ underweight + log_estEIR + frailty.gaussian(site, sparse=F), data=d, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
# mg per kg AQ
dASAQ<-subset(d,treatNum==1)
sASAQ <- Surv(time=dASAQ$time1Surv, time2=dASAQ$time2Surv,  type = "interval2" )
fit <- survreg(sASAQ ~ mgPerKgAQ10 + log_estEIR + frailty.gamma(site, sparse=F), data=dASAQ, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
# mg per kg AL
dAL<-subset(d,treatNum==0)
sAL <- Surv(time=dAL$time1Surv, time2=dAL$time2Surv,  type = "interval2" )
fit<-survreg(sAL ~ mgPerKgLum10 + log_estEIR + frailty.gaussian(site, sparse=F), data=dAL,dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))

# ASAQ formulation. 
d$form2<-NA
d$form2[d$form=='fdc']<-1
d$form2[d$form=='coblistered nfdc']<-2
d$form2[d$form=='loose nfdc']<-3
d$form2<-as.factor(d$form2)
dASAQ<-subset(d,treatNum==1 & !is.na(d$form))
sASAQ <- Surv(time=dASAQ$time1Surv, time2=dASAQ$time2Surv,  type = "interval2" )
fit <- survreg(sASAQ ~ form2 + log_estEIR + frailty.gaussian(site, sparse=F), data=dASAQ, dist="lognormal",
               control = list(maxiter=100,rel.tolerance=5e-4))
summary(fit)
tab(fit$coeff[2],sqrt(fit$var[2,2]))
tab(fit$coeff[3],sqrt(fit$var[3,3]))

## mutants
# 76T AL
dAL<-subset(d,treatNum==0)
sAL <- Surv(time=dAL$time1Surv, time2=dAL$time2Surv,  type = "interval2" )
fit<-survreg(sAL ~ prev76_10 + log_estEIR + frailty.gaussian(site, sparse=F), data=dAL, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
# 76T AL
dASAQ<-subset(d,treatNum==1 & !is.na(d$form))
sASAQ <- Surv(time=dASAQ$time1Surv, time2=dASAQ$time2Surv,  type = "interval2" )
fit<-survreg(sASAQ ~ prev76_10 + log_estEIR + frailty.gaussian(site, sparse=F), data=dASAQ, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))

# 86Y
fit<-survreg(sAL ~ prev86_10 + log_estEIR + frailty.gaussian(site, sparse=F), data=dAL, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))
fit<-survreg(sASAQ ~ prev86_10 + log_estEIR + frailty.gaussian(site, sparse=F), data=dASAQ, dist="lognormal")
tab(fit$coeff[2],sqrt(fit$var[2,2]))




############ MULTIVARIATE REGRESSION #########################.
## SEPARATELY BY DRUG AND BY MUTANT
########################################################################
#### PREV 76 AND AL  (cannot look at 86 and 76 mutations together, too correlated)
########################################################################
multiAL<-d[which(!is.na(d$prev76_10) & !is.na(d$mgPerKgLum)
                  & !is.na(d$ageyears2) & d$treatNum==0),]
sAL <- Surv(time=multiAL$time1Surv, time2=multiAL$time2Surv,  type = "interval2" )
fit<-survreg(sAL ~ log_estEIR + ageyears2 + ageyears2_sq + ageyears2_cub + mgPerKgLum10 + prev76_10 +
                frailty.gaussian(site, sparse=F), 
             data=multiAL, dist="lognormal")
summary(fit)
# coefficients
tab(fit$coeff[2],sqrt(fit$var[2,2])) #logEIR
tab(fit$coeff[3],sqrt(fit$var[3,3]))  #age
tab(fit$coeff[4],sqrt(fit$var[4,4])) #age sq
tab(fit$coeff[5],sqrt(fit$var[5,5]),4)  #age cub
tab(fit$coeff[6],sqrt(fit$var[6,6]))  #mgperkglum10
tab(fit$coeff[7],sqrt(fit$var[7,7])) #prev76_10

################################################################################
## check age effect - still significant.
fit2<-survreg(sAL ~ prev76_10 + mgPerKgLum10 +
                log_estEIR + frailty.gaussian(site, sparse=F), 
             data=multiAL, dist="lognormal")
-2*fit$loglik[2]- (-2*fit2$loglik[2]) + 3

##############################################################################################
### PREV 86 AND AL 
##############################################################################################
multiAL<-d[which(!is.na(d$prev86_10) & !is.na(d$mgPerKgLum)
                 & !is.na(d$ageyears2) & d$treatNum==0),]
sAL <- Surv(time=multiAL$time1Surv, time2=multiAL$time2Surv,  type = "interval2" )
fit<-survreg(sAL ~ log_estEIR + ageyears2 + ageyears2_sq + ageyears2_cub + mgPerKgLum10 + prev86_10 +
                 frailty.gaussian(site, sparse=F), 
             data=multiAL, dist="lognormal")
summary(fit)
# coefficients
tab(fit$coeff[2],sqrt(fit$var[2,2])) #logEIR
tab(fit$coeff[3],sqrt(fit$var[3,3]))  #age
tab(fit$coeff[4],sqrt(fit$var[4,4])) #age sq
tab(fit$coeff[5],sqrt(fit$var[5,5]),4)  #age cub
tab(fit$coeff[6],sqrt(fit$var[6,6]))  #mgperkglum10
tab(fit$coeff[7],sqrt(fit$var[7,7])) #prev86_10

##############################################################################
# check age - yes important
fit2<-survreg(sAL ~ prev86_10 + mgPerKgLum10 +
                + log_estEIR + frailty.gaussian(site, sparse=F), 
              data=multiAL, dist="lognormal")
-2*fit$loglik[2]- (-2*fit2$loglik[2]) + 3
myLrt(fit2,fit)


##############################################################################################
### PREV 76 AND ASAQ 
##############################################################################################
multiAQ<-d[which(!is.na(d$ageyears2) & d$treatNum==1 & !is.na(d$prev76)),]
sAQ <- Surv(time=multiAQ$time1Surv, time2=multiAQ$time2Surv,  type = "interval2" )
fit<-survreg(sAQ ~ log_estEIR + 
               ageyears2 + ageyears2_sq + ageyears2_cub + prev76_10 + frailty.gaussian(site, sparse=F), 
             data=multiAQ, dist="lognormal")
summary(fit)
## coefficients
tab(fit$coeff[2],sqrt(fit$var[2,2])) #logEIR
tab(fit$coeff[3],sqrt(fit$var[3,3]))  #age
tab(fit$coeff[4],sqrt(fit$var[4,4])) #age sq
tab(fit$coeff[5],sqrt(fit$var[5,5]),4)  #age cub
tab(fit$coeff[6],sqrt(fit$var[6,6])) #prev76_10

# check age.
fit2<-survreg(sAQ ~ + log_estEIR + prev76_10 + frailty.gaussian(site, sparse=F), 
             data=multiAQ, dist="lognormal")
summary(fit2)
-2*fit$loglik[2]- (-2*fit2$loglik[2]) + 3

##############################################################################################
### PREV 86 AND ASAQ
##############################################################################################
multiAQ<-d[which(!is.na(d$prev86_10)  & !is.na(d$ageyears2) & d$treatNum==1),]
sAQ <- Surv(time=multiAQ$time1Surv, time2=multiAQ$time2Surv,  type = "interval2" )
fit<-survreg(sAQ ~ log_estEIR +  
               ageyears2 + ageyears2_sq + ageyears2_cub + prev86_10  +frailty.gaussian(site, sparse=F), 
             data=multiAQ, dist="lognormal")
summary(fit)

# coefficients
tab(fit$coeff[2],sqrt(fit$var[2,2])) #logEIR
tab(fit$coeff[3],sqrt(fit$var[3,3]))  #age
tab(fit$coeff[4],sqrt(fit$var[4,4])) #age sq
tab(fit$coeff[5],sqrt(fit$var[5,5]),4)  #age cub
tab(fit$coeff[6],sqrt(fit$var[6,6])) #prev86_10

## age p value
fit2<-survreg(sAQ ~ prev86_10  +  log_estEIR + frailty.gaussian(site, sparse=F), 
             data=multiAQ, dist="lognormal")
myLrt(fit2,fit)

