## Prepare data for accelerated failure time models.
# Data required
# WHO anthro files (WORLD HEALTH ORGANIZATION 2011), igrowup R package  https://www.who.int/childgrowth/software/en/  Accessed 14/12/18.
# Individual patient clinical trial data from WWARN. AL and AS-AQ trial arms from studies with WWARN IDs: SBCEE, EGYMA, GPXJK, ACHGC, JXZNZ,QRBRC,UBTXH,YGTAH


library(foreign)
library(plyr)
# safer version of 'sample' function that can deal with vectors of length 1.
resample <- function(x,size0,replace0=F) x[sample.int(length(x),size=size0,replace=replace0)]


d<-read.delim("data/regression_data.txt")

#Split pid of studies with multiple obs per pers, only use 1st
d$obsNo <- 1
#For SBCEE
SBCEE <- d$sid == "SBCEE"
pidAndObs <- ldply(strsplit(as.character(d[d$sid == "SBCEE", "pid"]),split= "\\+"))
d$pid2<-d$pid
d$pid<-as.character(d$pid2)
d[d$sid == "SBCEE", "pid"] <- pidAndObs[,1]
d[d$sid == "SBCEE", "obsNo"] <- pidAndObs[,2]
#For EGYMA
EGYMA <- d$sid == "EGYMA"
pidAndObs <- ldply(strsplit(d[d$sid == "EGYMA", "pid"],split= "\\-"))
d[d$sid == "EGYMA", "pid"] <- pidAndObs[,1]
d[d$sid == "EGYMA", "obsNo"] <- pidAndObs[,2]
#cleanup
d <- d[d$obsNo==1,]
d$obsNo <- NULL


d <- d[d$site!="multi-site",] #5 sites in Senegal, Mali, Cameroon & Madascar, but cannot match..
#remove sites with only AL
d <- d[d$site!="rwrukara",] #only AL
d <- d[d$site!="ugjinja",] #only AL
d <- d[d$site!="rwmashes",] #only AL
d <- d[d$site!="ugtororo",] #only AL
#remove sites with only ASAQ
d <- d[d$site!="mzmanhic",] #only ASAQ
d <- d[d$site!="ugmbarar",] #only ASAQ

#correct data types..
d$sid <- factor(d$sid)
d$pid <- factor(d$pid)
d$site <- factor(d$site)
d$ageyears <- as.numeric(d$ageyears)
d$weight <- as.numeric(d$weight)
d$treat <- as.factor(d$treat)
d$outcome <- as.factor(d$outcome)
d$lastday <- as.numeric(d$lastday)
#same with additional covariates
d$fever0 <- as.integer(d$fever0)
d$gam0 <- as.numeric(d$gam0)
d$gender <- as.factor(d$gender)
d$hb0 <- as.numeric(d$hb0)
d$pardens <- as.numeric(d$pardens)
d$pfmicl0 <- as.numeric(d$pfmicl0)
d$pfmicl1 <- as.numeric(d$pfmicl1)
d$pfmicl2 <- as.numeric(d$pfmicl2)
d$pfmicl3 <- as.numeric(d$pfmicl3)
d$sevmalaria <- as.numeric(d$sevmalaria)
d$ht0 <- as.numeric(d$ht0)
d$fever <- as.numeric(d$fever)
d$liver <- as.numeric(d$liver)
d$diarrhea <- as.numeric(d$diarrhea)
d$spleen <- as.numeric(d$spleen)
d$temp <- as.numeric(d$temp)
d$vomit <- as.numeric(d$vomit)

#code factors as numeric variables 
d$sidNum <- as.numeric(d$sid)
d$pidNum <- as.numeric(d$pid)
d$siteNum <- as.numeric(d$site)
d$treatNum <- as.numeric(d$treat)
d$outcomeNum <- as.numeric(d$outcome)

d$inInterval[d$outcome=="ACPR"] <- 2
d$inInterval[d$outcome=="LTFRI"] <- 1
d$times <- NA #all times estimated, only whether right or left censored is known
d$maxObsTime <- d$lastday #for the right censored ones, lastday is censoring day
precedingScreen <- Vectorize(function(lastDay){
  screenDays <- c(3, 7, 14, 21, 28, 35, 42) #for interval censoring
  int <- findInterval(lastDay-0.5, screenDays)
  screenDays[int]
})
d$lastNegTime <- precedingScreen(d$lastday) #left limit for interval censoring

#set al non-al as asaq
d[d$treatNum!=1, "treatNum"] <- 2
d <- d[order(d$siteNum),]


#add variables for the Surv function in R package "survival"
d$time1Surv <- NA 
d$time2Surv <- NA
d$statusSurv <- NA
#right censored obs
d[d$outcome =="ACPR","time1Surv"] <- d[d$outcome =="ACPR","lastday"] 
#already so: d[d$outcome =="ACPR","time2Surv"] <- NA
d[d$outcome =="ACPR","statusSurv"] <- 0  #not sure if needed
#there are no left censored ones.. 
#interval censored obs
d[d$outcome =="LTFRI","time1Surv"] <- d[d$outcome =="LTFRI","lastNegTime"] 
d[d$outcome =="LTFRI","time2Surv"] <- d[d$outcome =="LTFRI","lastday"] 
d[d$outcome =="LTFRI","statusSurv"] <- 3 #not sure if needed


#### CORRECT KISUMU AGES - THEY ARE IN MONTHS
inds<-which(d$siteNum==6)
d$ageyears[inds]<-d$ageyears[inds]/12


#remove covariates that have too many missing entries
d$fever0 <- NULL
d$sevmalaria <- NULL
d$ht0 <- NULL
d$fever <- NULL
d$diarrhea <- NULL
d$vomit <- NULL

#calculate parasite reduction ratio from days 0 and 1 of parasitaemia, factor two because asexual cycle would be 48h
d$log10PRR <- -2*(log10(d$pfmicl1+1) - log10(d$pfmicl0+1))
#these are not needed anymore
d$pfmicl0 <- NULL 
d$pfmicl1 <- NULL
d$pfmicl2 <- NULL
d$pfmicl3 <- NULL

#reduce the number of levels in the site factor to the actual levels present in this subset
d$site <- as.factor(as.character(d$site))

#write.table(d, file="data/regression_data3.txt", sep="\t", row.names=F, quote=F)

source(paste0(path,"my_R_lib/recode.R"))
source(paste0(path,"my_R_lib/egen.R"))

d<-read.delim("data/regression_data3.txt")

sites_mut<-read.csv("data/sites_all_mut2.csv")

d <- merge(d, sites_mut , by.x="siteNum", by.y="siteNo")
names(d)[names(d)=="sid.x"]<-"sid"
d$treatNum[d$treatNum==1]<-0   ### must do this otherwise regression not right!
d$treatNum[d$treatNum==2]<-1

### Delete 6 individuals with implausible WEIGHT FOR AGE
inds<-which(d$ageyears<8 & d$weight>30)
d$ageyears[inds]<-NA
d$weight[inds]<-NA
inds<-which(d$ageyears>10 & d$weight<10)
d$ageyears[inds]<-NA
d$weight[inds]<-NA

########## IMPUTE AGE FOR SITE 1: TORORO YEKA.
d$gender2<-as.numeric(d$gender)
d$gender2[d$gender=='No Gender Variable']<-NA
inds<-which(d$sid=='SBCEE')
d$ageyears[inds]<-NA
weights_ug<-unique(d$weight[inds])
for(i in weights_ug) {
  for(j in 1:2) {
    inds_i<-which(d$weight>=(i-0.5) & d$weight<(i+0.5) & !is.na(d$ageyears) & d$sid!='EGYMA' & 
                  d$sid!='SBCEE' & d$gender2==j)
    inds_i2<-which(d$sid=='SBCEE' & d$weight==i & d$gender2==j)
    if(length(inds_i)>1) d$ageyears[inds_i2]<-resample(d$ageyears[inds_i],size0=length(inds_i2))
    if(length(inds_i)==1) d$ageyears[inds_i2]<-d$ageyears[inds_i]
  }
}

### # IMPUTE WEIGHT FOR SITE 10 Sikasso Mali
inds<-which(d$sid=='EGYMA')
d$weight[inds]<-NA
ages_mali<-unique(d$ageyears[inds])
for(i in ages_mali) {
  for(j in 1:2) {
    if(i<25) inds_i<-which(d$ageyears>=(i-0.5) & d$ageyears<(i+0.5) & !is.na(d$weight) 
                           & d$sid!='EGYMA' & d$sid!='SBCEE' & d$gender2==j)
    if(i>=25) inds_i<-which(d$ageyears>=(i-5) & d$ageyears<(i+5) & !is.na(d$weight))
    inds_i2<-which(d$sid=='EGYMA' & d$ageyears==i & d$gender2==j)
    if(length(inds_i)>1) d$weight[inds_i2]<-resample(d$weight[inds_i],size0=length(inds_i2))
    if(length(inds_i)==1) d$weight[inds_i2]<-d$weight[inds_i]
  }  ## gender is coded 1-2
}


##################################################
# Add dosage info
##################################################
d$mgPerKgAQ <- NA
d$mgPerKgLum <- NA

# 4ABC study. 
# ASAQ
dosAQperKgGPXJK <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 67.5 #in mg
  AqTab2 <- 135 #in mg
  AqTab3 <- 270 #in mg
  if (weight < 9) {dos <- 3*AqTab1}
  if (weight >= 9 & weight < 18) {dos <- 3*AqTab2}
  if (weight >= 18 & weight <= 36) {dos <- 3*AqTab3}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="GPXJK" & d$treatNum==1)] <- 
  dosAQperKgGPXJK(d$weight[which(d$sid=="GPXJK" & d$treatNum==1)])

# AL
# weight = 5-14 kg: one # tablet per dose; weight= 15-24 kg: two tablets per dose;
# weight =25-34 kg: three tablets per dose. 
#### no dosing info in the dataset, but clear from the paper.
dosLumPerKgGPXJK <- Vectorize(function(weight){#mg
  dos <- NA
  LumTab1 <- 120 #in mg
  if (weight >=5 & weight < 15) {dos <- 1*LumTab1*6} ## 6 dose regimen in total.
  if (weight >= 15 & weight < 25) {dos <- 2*LumTab1*6}
  if (weight >= 25 & weight <= 35) {dos <- 3*LumTab1*6}
  if (weight >= 35 ) {dos <- 4*LumTab1*6}
  return(dos/weight)
})
d$mgPerKgLum[which(d$sid=="GPXJK" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="GPXJK" & d$treatNum==0)])


# ACHGC is Juma 2005 unpublished.
dosAQperKgACHGC  <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 200 #in mg
  if (weight < 8) {dos <- 3*0.25*AqTab1}
  if (weight >= 8 & weight < 13) {dos <- 2*0.5*AqTab1 + 1*0.25*AqTab1}
  if (weight >= 13 & weight <= 14) {dos <- 2*0.75*AqTab1 + 1*0.25*AqTab1}
  if (weight >= 15 & weight <= 17) {dos <- 2*0.75*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 18 & weight <= 22) {dos <- 2*1*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 23 & weight <= 24) {dos <- 2*1.25*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 25 & weight <= 32) {dos <- 2*1.5*AqTab1 + 1*0.75*AqTab1}
  if (weight >= 33 & weight <= 37) {dos <- 2*1.75*AqTab1 + 1*1*AqTab1}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="ACHGC" & d$treatNum==1)]<- 
  dosAQperKgACHGC(d$weight[which(d$sid=="ACHGC" & d$treatNum==1)])

d$mgPerKgLum[which(d$sid=="ACHGC" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="ACHGC" & d$treatNum==0)])


## JXZNZ. Nikiema et al
# remove the one obs in this study with no weight data
d$miss<-NA
d$miss[d$sid=="JXZNZ" & is.na(d$weight)]<-1
d<-subset(d, is.na(miss))

dosAQperKgJXZNZ <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 67.5 #in mg
  AqTab2 <- 135 #in mg
  AqTab3 <- 270 #in mg
  if (weight < 9) {dos <- 3*AqTab1}
  if (weight >= 9 & weight < 18) {dos <- 3*AqTab2}
  if (weight >= 18 & weight <= 36) {dos <- 3*AqTab3}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="JXZNZ" & d$treatNum==1)]<- 
  dosAQperKgJXZNZ(d$weight[which(d$sid=="JXZNZ" & d$treatNum==1)])

### Lum.
d$mgPerKgLum[which(d$sid=="JXZNZ" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="JXZNZ" & d$treatNum==0)])

# QRBRC - Espie 2012. 
dosAQperKgQRBRC <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 67.5 #in mg
  AqTab2 <- 135 #in mg
  if (weight < 9) {dos <- 3*AqTab1}
  if (weight >= 9 & weight < 18) {dos <- 3*AqTab2}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="QRBRC" & d$treatNum==1)]<- 
  dosAQperKgQRBRC(d$weight[which(d$sid=="QRBRC" & d$treatNum==1)])

## AL standard co-artem (same as 4ABC) 
d$mgPerKgLum[which(d$sid=="QRBRC" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="QRBRC" & d$treatNum==0)])


# SBCEE - Yeka 2014.
dosAQperKgSBCEE <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 67.5 #in mg
  AqTab2 <- 135 #in mg
  AqTab3 <- 270 #in mg
  if (weight < 9) {dos <- 3*AqTab1}
  if (weight >= 9 & weight < 18) {dos <- 3*AqTab2}
  if (weight >= 18 & weight <= 36) {dos <- 3*AqTab3}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="SBCEE" & d$treatNum==1)]<- 
  dosAQperKgSBCEE(d$weight[which(d$sid=="SBCEE" & d$treatNum==1)])

## Lum: same dosing as 4ABC, standard.
d$mgPerKgLum[which(d$sid=="SBCEE" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="SBCEE" & d$treatNum==0)])


## UBTXH - Schramm 2013 (called Schramm 2013b in wwarn AL dosing paper)
dosAQperKgUBTXH <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 67.5 #in mg
  AqTab2 <- 135 #in mg
  AqTab3 <- 270 #in mg
  if (weight < 9) {dos <- 3*AqTab1}
  if (weight >= 9 & weight < 18) {dos <- 3*AqTab2}
  if (weight >= 18 & weight < 36) {dos <- 3*AqTab3}
  if (weight >= 36 ) {dos <- 3*2*AqTab3}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="UBTXH" & d$treatNum==1)]<- 
  dosAQperKgUBTXH(d$weight[which(d$sid=="UBTXH" & d$treatNum==1)])

### Lum. same standard dosing as 4ABC.
d$mgPerKgLum[which(d$sid=="UBTXH" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="UBTXH" & d$treatNum==0)])


## YGTAH - Burkirwa 2006.
dosAQperKgYGTAH <- Vectorize(function(weight){#mg
  dos <- NA
  AqTab1 <- 200 #in mg
  if (weight < 12) {dos <- 2*0.5*AqTab1 + 1*0.25*AqTab1}
  if (weight >= 12 & weight < 14) {dos <- 3*0.5*AqTab1}
  if (weight >= 14 & weight < 16) {dos <- 1*0.75*AqTab1 + 2*0.5*AqTab1}
  if (weight >= 16 & weight < 18) {dos <- 2*0.75*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 18 & weight < 20) {dos <- 1*1*AqTab1 + 1*0.75*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 20 & weight < 22) {dos <- 2*1*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 22 & weight < 24) {dos <- 1*1.25*AqTab1 + 1*1*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 24 & weight < 26) {dos <- 2*1.25*AqTab1 + 1*0.5*AqTab1}
  if (weight >= 26 & weight < 28) {dos <- 3*1*AqTab1}
  if (weight >= 28) {dos <- 2*1.25*AqTab1 + 1*1*AqTab1}
  return(dos/weight)
})
d$mgPerKgAQ[which(d$sid=="YGTAH" & d$treatNum==1)]<- 
  dosAQperKgYGTAH(d$weight[which(d$sid=="YGTAH" & d$treatNum==1)])
### AL - same dosing as 4 ABC (GPXJK)
d$mgPerKgLum[which(d$sid=="YGTAH" & d$treatNum==0)] <- 
  dosLumPerKgGPXJK(d$weight[which(d$sid=="YGTAH" & d$treatNum==0)])


## Recode some variables for regression:
d$ageyears2<-d$ageyears
d$ageyears2[which(d$ageyears2>20)]<-20
d$ageyears2_sq<-d$ageyears2^2
d$ageyears2_cub<-d$ageyears2^3
d$mgPerKgAQ10<-d$mgPerKgAQ/10
d$mgPerKgLum10<-d$mgPerKgLum/10

formulation<-c("fdc", "fdc","fdc","fdc","fdc","loose nfdc","fdc","fdc","fdc","coblistered nfdc","loose nfdc","fdc")
d$form<- c("fdc", "fdc","fdc","fdc","fdc","loose nfdc","fdc","fdc","fdc","nfdc coblistered","loose nfdc","fdc")[as.numeric(d$site)]
d$form<-(as.factor(d$form))
dASAQ<-subset(d,treatNum==1)


# ASAQ  formulation
d$form<-"fdc"
d$form[d$desc=="Tororo UG [Bukirwa]" | d$desc=="Kisumu KE"]<-"loose nfdc"
d$form[d$desc=="Sikasso ML"]<-"coblistered nfdc"
d$form<-(as.factor(d$form))

# temp
d$fev[d$temp<=37.5]<-0
d$fev[d$temp>37.5]<-1
# hb0 ## try as anaemia.
d$anaemic[d$hb0<10]<-1
d$anaemic[d$hb0>=10]<-0
# gender
d$gender<-as.character(d$gender)
d$gender[d$gender=="No Gender Variable"]<-NA
d$gender<-as.factor(d$gender)
#### mutants on better scale for regression results.
d$prev76_10<-d$prev76*10
d$prev86_10<-d$prev86*10


## make underweight variable - analysis from WHO weight for age calculations.igrowup R package.
weianthro<-read.table("weianthro.txt",header=T,sep="",skip=0) # downloaded from WHO.

mat<-d
mat$age.days=round(mat$ageyears*365)
mat$gender<-as.character(mat$gender)
mat$sex[mat$gender=="M"]<-1
mat$sex[mat$gender=="F"]<-2
mat$zwei<-NA

for(i in 1:length(mat$age.days)) {
  
  if(!is.na(mat$age.days[i]) & !is.na(mat$weight[i]) & !is.na(mat$sex[i]) & mat$age.days[i]<max(weianthro$age)) {
    
    l.val<-weianthro$l[weianthro$age==mat$age.days[i] & weianthro$sex==mat$sex[i]]
    m.val<-weianthro$m[weianthro$age==mat$age.days[i] & weianthro$sex==mat$sex[i]]
    s.val<-weianthro$s[weianthro$age==mat$age.days[i] & weianthro$sex==mat$sex[i]]
    
    mat$zwei[i]<-(((mat$weight[i]/m.val)^l.val)-1)/(s.val*l.val)
    if(!is.na(mat$zwei[i]) & mat$zwei[i]>3) {
      sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
      sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
      mat$zwei[i]<- 3+((mat$weight[i]-sd3pos)/sd23pos)
    }
    if(!is.na(mat$zwei[i]) & mat$zwei[i]< (-3)) {
      sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
      sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
      mat$zwei[i]<- (-3)+((mat$weight[i]-sd3neg)/sd23neg)
    }
    
  } else mat$zwei[i]<-NA
}
d$zwei<-mat$zwei
d$zwei[d$zwei<(-6)]<-NA
d$zwei[d$zwei>(6)]<-NA
d$underweight[d$zwei<(-2)]<-1
d$underweight[d$zwei>=(-2)]<-0
d$underweight[is.na(d$zwei)]<-NA

d$log_estEIR<-log(d$estEIR)

#write.csv(d,file="data/regression_data6all.csv")

