############ PREPARE CLINICAL TRIAL DATA FOR ANALYSIS BY HIDDEN SEMI MARKOV MODEL
# Data required:
# Individual patient clinical trial data from WWARN. AL and AS-AQ trial arms from studies with WWARN IDs: SBCEE, EGYMA, GPXJK, ACHGC, JXZNZ,QRBRC,UBTXH,YGTAH


### INSTALL AND LOAD LIBRARIES
# install.packages("foreign")
# install.packages("plyr")
# install.packages("raster")
# install.packages("rgdal")


library(foreign)
library(plyr)
library(raster)  # for raster file extraction
library(rgdal)   # for raster file extraction

##########################################################################################
#	MAIN MATTER
##########################################################################################
d <- read.csv(file="data/WwarnDataSpreadsheet.csv",stringsAsFactors = F)
# 
# #correct data types..
d$sid <- factor(d$sid)
d$pid <- factor(d$pid)
d$site <- factor(d$site)
d$treat <- as.factor(d$treat)
d$outcome <- as.factor(d$outcome)
 
# #code factors as numeric variables since jags can't read factors
d$sidNum <- as.numeric(d$sid)
d$pidNum <- as.numeric(d$pid)
d$treatNum <- as.numeric(d$treat)
# d$outcomeNum <- as.numeric(d$outcome)

# create site number for each location
d$siteNum[which(d$sid=='SBCEE')]<-1
d$siteNum[which(d$site=='bfnanoro' & d$sid=='GPXJK')]<-2
d$siteNum[which(d$site=='Bobo Dioulasso' & d$sid=='JXZNZ')]<-3
d$siteNum[which(d$site=='gafougam' & d$sid=='GPXJK')]<-4
d$siteNum[which(d$site=='Gourcy' & d$sid=='JXZNZ')]<-5
d$siteNum[which(d$sid=='ACHGC')]<-6
d$siteNum[which(d$sid=='UBTXH')]<-7
d$siteNum[which(d$site=='ngikot' & d$sid=='GPXJK')]<-8
d$siteNum[which(d$sid=='QRBRC')]<-9
d$siteNum[which(d$sid=='EGYMA')]<-10
d$siteNum[which(d$sid=='YGTAH')]<-11
d$siteNum[which(d$site=='zmndola' & d$sid=='GPXJK')]<-12

# #columns for jags survival analysis
# #this version uses right and left censoring
d$inInterval[d$outcome=="ACPR"] <- 2
d$inInterval[d$outcome=="LTFRI"] <- 1
d$times <- NA #all times estimated, only whether right or left censored is known
d$maxObsTime <- d$lastday #for the right censored, lastday is censoring day
precedingScreen <- Vectorize(function(lastDay){
  screenDays <- c(3, 7, 14, 21, 28, 35, 42) #for interval censoring
  int <- findInterval(lastDay-0.5, screenDays)
  screenDays[int]
})
d$lastNegTime <- precedingScreen(d$lastday) #left limit for interval censoring

#set Non-AL as ASAQ
d[d$treatNum!=1, "treatNum"] <- 2
 
#sort our data by site number
d <- d[order(d$siteNum),]

uniqueSites <- unique(d$site)
sites<-read.csv("data/sites_all.csv")
names(sites)[names(sites)=="pr210"]<-"pr210_old"
### First get Malaria Atlas project prevalence for each Long Lat. 5x5 km resolution, all years 2000-2015
# download raster files of prevalence in 2-10 year olds across Africa from 2000-2014 
  #from Malaria Atlas Project website. https://map.ox.ac.uk/
LL<-data.frame(x=sites$Longitude, y=sites$Latitude)   ## lon lats of  data. Call these x and y so the raster function can understand them.
for(i in 2000:2015) {
  raster_filename<-paste(path,"MAP/rasters/MODEL43.",i,".PR.rmean.stable.tif",sep="")
  temp<-raster(raster_filename)
  prev_var_name<-paste("pfpr",i,sep="")
  sites[[prev_var_name]]<-NA
  sites[[prev_var_name]]<-extract(temp,LL,buffer=20000,fun=mean) ## buffer in METRES from the original point, ie 10k=10000
  print(i)
}

######## select slide prevalence at the time of each study ######################
### first for studies which just collected data in one year
inds<-which(sites$year_end==sites$year_start) 
sites$pr210<-NA
for(i in inds) sites$pr210[i]<-sites[[paste0("pfpr",sites$year_start[i])]][i] ### none are before 2000

### now for studies collecting data in multiple years
inds<-which(sites$year_end!=sites$year_start) 
for(i in inds) {
  years<-sites$year_start[i]:sites$year_end[i]
  prevs<-vector(length=length(years))
  for(j in 1:length(years)) prevs[j]<-sites[[paste0("pfpr",years[j])]][i]
  sites$pr210[i]<-mean(prevs)
}

eirVsPr210 <- read.delim("data/EIRvsPR2_10_BMModel.csv",stringsAsFactors = F)

siteEIRs <- approx(x=eirVsPr210$PR2_10, y = eirVsPr210$EIR, xout=sites$pr210)$y


siteDesc <- NA
siteInfo <- sites
siteInfo$eir<-siteEIRs

siteInfo$desc[siteInfo$SiteCode=="1"] <- "Tororo UG [Yeka]"  # changed to Tororo from Kampala, LO.
siteInfo[siteInfo$SiteCode=="bfnanoro", "desc"] <- "Nanoro BF"    
siteInfo[siteInfo$SiteCode=="Bobo Dioulasso", "desc"] <- "Bobo Dioulasso BF"
siteInfo[siteInfo$SiteCode=="gafougam", "desc"] <- "Fougamou GA"       
siteInfo[siteInfo$SiteCode=="Gourcy", "desc"] <-   "Gourcy BF"
siteInfo[siteInfo$SiteCode=="Kisumu", "desc"] <-   "Kisumu KE"        
siteInfo[siteInfo$SiteCode=="LIBERIA", "desc"] <- "Nimba LR"
siteInfo[siteInfo$SiteCode=="ngikot", "desc"] <- "Pamol NG"          
siteInfo[siteInfo$SiteCode=="Pweto", "desc"] <- "Pweto CD"        
siteInfo[siteInfo$SiteCode=="Sikasso MALI", "desc"] <- "Sikasso ML"  
siteInfo[siteInfo$SiteCode=="Tororo", "desc"] <- "Tororo UG [Bukirwa]"        
siteInfo[siteInfo$SiteCode=="zmndola", "desc"] <- "Ndola ZM"

siteInfo$country_full<-NA
siteInfo$country_full[siteInfo$country=='UG']<-"Uganda"
siteInfo$country_full[siteInfo$country=='BF']<-"Burkina Faso"
siteInfo$country_full[siteInfo$country=='GA']<-"Gabon"
siteInfo$country_full[siteInfo$country=='BF']<-"Burkina Faso"
siteInfo$country_full[siteInfo$country=='KE']<-"Kenya"
siteInfo$country_full[siteInfo$country=='LR']<-"Liberia"
siteInfo$country_full[siteInfo$country=='NG']<-"Nigeria"
siteInfo$country_full[siteInfo$country=='CD']<-"Democratic Republic of Congo"
siteInfo$country_full[siteInfo$country=='ML']<-"Mali"
siteInfo$country_full[siteInfo$country=='ZM']<-"Zambia"



siteInfoOrdered <- siteInfo[order(siteInfo$pr210, decreasing = F),]
siteInfoOrdered$siteNo <- as.numeric(rownames(siteInfoOrdered)) #identical to siteNum in d!!

#### AGE DATA IN KISUMU - convert from months to years
inds<-which(d$siteNum==6)
d$ageyears[inds]<-d$ageyears[inds]/12

d$ageyears[inds]<-NA
d$weight[inds]<-NA
inds<-which(d$ageyears>10 & d$weight<10)
d$ageyears[inds]<-NA
d$weight[inds]<-NA

########## IMPUTE AGE FOR SITE 1: TORORO YEKA.
d$gender2<-as.numeric(d$gender)
d$gender2[d$gender=='No Gender Variable']<-NA
inds<-which(d$sid=='SBCEE')
weights_ug<-unique(d$weight[inds])
for(i in weights_ug) {
  for(j in 1:2) {
    inds_i<-which(d$weight>=(i-0.5) & d$weight<(i+0.5) & !is.na(d$ageyears) & d$sid!='EGYMA' & 
                    d$sid!='SBCEE' & d$gender2==j)
    inds_i2<-which(d$sid=='SBCEE' & d$weight==i & d$gender2==j)
    if(length(inds_i)>1) d$ageyears[inds_i2]<-sample(d$ageyears[inds_i],length(inds_i2))
    if(length(inds_i)==1) d$ageyears[inds_i2]<-d$ageyears[inds_i]
  }
}


##### now create summary dataset for plotting
days<-c(0,7,14,21,28,42)
siteNos<-siteInfoOrdered$siteNo
for(i in 1:length(siteNos)) {
  temp<-data.frame(siteNo=rep(siteNos[i], length(days)), days=days)
  temp$fail_in_interval<-temp$N_in_interval<-temp$prop_no_reinfection<-NA
  temp$drug<-NA
  for(drug in 1:2) {
    siteDataTemp<-subset(d, siteNum==siteNos[i] & treatNum==drug)
    if(drug==1) temp$drug<-"AL"
    if(drug==2) temp$drug<-"ASAQ"
    temp$fail_in_interval[1]<-0
    for(j in 2:length(days)) temp$fail_in_interval[j]<-nrow(subset(siteDataTemp, outcome=="LTFRI" & lastday>days[j-1] & lastday<=days[j]))
    for(j in 1:length(days)) temp$N_in_interval[j]<-nrow(subset(siteDataTemp, lastday>=days[j]))
    temp$prop_fail_in_interval<-temp$fail_in_interval / temp$N_in_interval
    temp$prop_no_reinfection[1]<-1-temp$fail_in_interval[1]
    for(j in 2:nrow(temp)) {
      temp$prop_no_reinfection[j]<-temp$prop_no_reinfection[j-1]*(1-temp$prop_fail_in_interval[j])
    }
    if(i==1 & drug==1) { summ_data<-temp 
    } else {
      summ_data<-rbind(summ_data,temp)
    }
  }
}
summ_data<-subset(summ_data, N_in_interval>0)
summ_data$cumu_prop_reinfected<-1-summ_data$prop_no_reinfection
summ_data$cumu_prop_reinfected_lci<-
  binom.confint(round(summ_data$cumu_prop_reinfected*summ_data$N_in_interval), 
                summ_data$N_in_interval, methods = "exact")$lower
summ_data$cumu_prop_reinfected_uci<-
  binom.confint(round(summ_data$cumu_prop_reinfected*summ_data$N_in_interval), 
                summ_data$N_in_interval, methods = "exact")$upper
#write.csv(summ_data,"results/summ_reinfection_by_site.csv")
