# ANALYSE AND PLOT RESULTS OF THE HIDDEN SEMI MARKOV MODEL
# Files required (freely available online):
# output of hidden semi markov model 'codaSemiMarkov'
# site information: 'sites_all_mut2.csv' 
# proportions reinfected over time in each site: data. 'summ_reinfection_by_site.csv'
# proportions reinfected over time in each site: model predictions. 'predict_reinfect.csv'

#install.packages("coda")
#install.packages("vioplot")
#install.packages("binom")
#install.packages("gplots")
library("coda")
library(vioplot)
library(binom)
library(gplots)

## if rerunning the analysis...:
#source("dataPreparationSemiMarkov.R")
#source("SemiMarkov.R")

plot2file<-F   ## change to true in order to save plots.

# ...Or load results and data.
load("codaSemiMarkov")
mut<-read.csv("sites_all_mut2.csv",stringsAsFactors = F)
summ_data<-read.csv("summ_reinfection_by_site.csv",stringsAsFactors = F)
sims<-read.csv("predict_reinfect.csv",stringsAsFactors = F)

codaOneDur <- codaSamples
#combine chains
bothChains <- as.matrix(codaOneDur)
# colnames(bothChains)  ##show names of chains.

apply(bothChains, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm=T) ## updated LO
### in output, row 1 is AL and row 2 is ASAQ

## Output duration of prophylaxis in the whole dataset.
durP <- vector()
r <- vector()
durP[1] <- median(bothChains[,"siteDurP1HyperMean"])     
durP[2] <- median(bothChains[,"siteDurP2HyperMean"])     
r[1] <- median(bothChains[,"r[1]"])        
r[2] <- median(bothChains[,"r[2]"])
## scale params:
durP[1]/r[1]
durP[2]/r[2]

durP
r 

quantile(bothChains[,"siteDurP1HyperMean"],probs=c(0.025,0.5,0.975))
quantile(bothChains[,"siteDurP2HyperMean"],probs=c(0.025,0.5,0.975))

samplesAL <- bothChains[,2*(1:12)-1]
samplesASAQ <- bothChains[,2*(1:12)]

#### EXTRACT RESULTS FOR EACH SITE from MCMC chains.
siteInfoOrdered<-mut[order(mut$estEIR),]
inds<-grep("siteEir", colnames(bothChains))
estSiteEirs <- bothChains[,inds]
estEirsDat<-data.frame(siteNo=siteInfoOrdered$siteNo,siteCode=siteInfoOrdered$SiteCode, estEir=NA)
##### order the EIRs
estEirsDat<-estEirsDat[order(estEirsDat$estEir),]
##write to file to use in regression analysis.
#if(plot2file) write.csv(estEirsDat,file="estEirsDatSep2018.csv")

mut$dur_p_asaq<-mut$dur_p_al<-mut$lci_dur_p_asaq<-mut$uci_dur_p_asaq<-mut$lci_dur_p_al<-mut$uci_dur_p_al<-
  mut$estEIR<-mut$lci_estEIR<-mut$uci_estEIR<- NA
siteNos<-1:12
for(i in siteNos) {
  ind<-which(mut$siteNo==i)
  mut$dur_p_asaq[ind]<-median(samplesASAQ[,i])
  mut$lci_dur_p_asaq[ind]<-quantile(samplesASAQ[,i],probs=0.025)
  mut$uci_dur_p_asaq[ind]<-quantile(samplesASAQ[,i],probs=0.975)
  mut$lci_dur_p_al[ind]<-quantile(samplesAL[,i],probs=0.025)
  mut$uci_dur_p_al[ind]<-quantile(samplesAL[,i],probs=0.975)
  mut$dur_p_al[ind]<-median(samplesAL[,i])
  mut$estEIR[ind]<-median(estSiteEirs[,i])
  mut$lci_estEIR[ind]<-quantile(estSiteEirs[,i],probs=0.025)
  mut$uci_estEIR[ind]<-quantile(estSiteEirs[,i],probs=0.975)
} 

mut$prev76<-mut$x76t_pos/mut$x76t_N
mut$prev86<-mut$x86y_pos/mut$x86y_N
mut$prev76_lci<-mut$prev76_uci<-mut$prev86_lci<-mut$prev86_uci<-NA
inds<-which(!is.na(mut$x76t_pos))
mut$prev76_lci[inds]<-binom.confint(mut$x76t_pos[inds],mut$x76t_N[inds],methods="exact")$lower
mut$prev76_uci[inds]<-binom.confint(mut$x76t_pos[inds],mut$x76t_N[inds],methods="exact")$upper
inds<-which(!is.na(mut$x86y_pos))
mut$prev86_lci[inds]<-binom.confint(mut$x86y_pos[inds],mut$x86y_N[inds],methods="exact")$lower
mut$prev86_uci[inds]<-binom.confint(mut$x86y_pos[inds],mut$x86y_N[inds],methods="exact")$upper

### write prophylactic profiles to file.
#site 5 AL short, AQ long. # site 7 AL long, AQ short.
pAl <- Vectorize(function(x){1-pgamma(x, shape=r[1], rate=r[1]/median(samplesAL[,5]))})
x<-seq(0,30,0.2)
pp_l_short<-data.frame(x, age0=0,age1=200, pAl(x))
pp_l_short<-subset(pp_l_short,pp_l_short$pAl.x.>0.0001)
#write.table(pp_l_short,file="results/p_protect_lum_short.txt", sep="\t",row.names = F,col.names = F)

pAl <- Vectorize(function(x){1-pgamma(x, shape=r[1], rate=r[1]/median(samplesAL[,7]))})
x<-seq(0,30,0.2)
pp_l_long<-data.frame(x, age0=0,age1=200, pAl(x))
pp_l_long<-subset(pp_l_long,pp_l_long$pAl.x.>0.0001)
#write.table(pp_l_long,file="results/p_protect_lum_long.txt", sep="\t",row.names = F,col.names = F)

pAsaq <- Vectorize(function(x){1-pgamma(x, shape=r[2], rate=r[2]/median(samplesASAQ[,5]))})
x<-seq(0,30,0.2)
pp_a_long<-data.frame(x, age0=0,age1=200, pAsaq(x))
pp_a_long<-subset(pp_a_long,pp_a_long$pAsaq.x.>0.0001)
#write.table(pp_a_long,file="results/p_protect_aq_long.txt", sep="\t",row.names = F,col.names = F)

pAsaq <- Vectorize(function(x){1-pgamma(x, shape=r[2], rate=r[2]/median(samplesASAQ[,7]))})
x<-seq(0,30,0.2)
pp_a_short<-data.frame(x, age0=0,age1=200, pAsaq(x))
pp_a_short<-subset(pp_a_short,pp_a_short$pAsaq.x.>0.0001)
#write.table(pp_a_short,file="results/p_protect_aq_short.txt", sep="\t",row.names = F,col.names = F)

#visual check of pp files.
plot(pp_a_long$x,pp_a_long$pAsaq.x.,type="l")
lines(pp_a_short$x,pp_a_short$pAsaq.x.,type="l")
lines(pp_l_short$x,pp_l_short$pAl.x.,type="l",col="blue")
lines(pp_l_long$x,pp_l_long$pAl.x.,type="l",col="blue")


## Extract numbers for results table 1
test<-mut
test<-test[order(test$estEIR),]
#posteriors
res_table<-cbind(as.character(test$desc),paste0(round(test$dur_p_al,1)," (",round(test$lci_dur_p_al,1),"-",
                                                round(test$uci_dur_p_al,1),")"),
                 paste0(round(test$dur_p_asaq,1)," (",round(test$lci_dur_p_asaq,1),"-",
                        round(test$uci_dur_p_asaq,1),")"),
                 paste0(round(test$estEIR,1)," (",round(test$lci_estEIR,1),"-",
                        round(test$uci_estEIR,1),")"))
#write.csv(res_table,file="results/res_table2.csv")
#EIR priors
cbind(as.character(test$desc),test$priorEIR)


if(plot2file) write.csv(mut,file="results/sites_all_mut_posteriors.csv")


### Figure 1 Duration of post-treatment prophylaxis 
if(plot2file) {
  tiff(file="fig1_protection.tiff", compression = "lzw",res=300,width=4500,height=1500)
  layout(matrix(c(1,2,3), 1,3, byrow = TRUE))
  siteNos<-1:12
  # panel 1:Violin Plots
  x1 <- bothChains[,"siteDurP1HyperMean"]
  x2 <- bothChains[,"siteDurP2HyperMean"]
  plot(NA, NA, xlim=c(0,2), ylim=c(0,30), type="n", axes=F, xlab=NA, 
       ylab="days of protection after treatment", cex.axis=1.5, cex.lab=1.5)
  vioplot(x1, col="blue", horizontal=F, at=0.5, add=TRUE, names="AL", h=1.5)
  vioplot(x2, col="forestgreen", horizontal=F, at=1.5, add=TRUE, names="ASAQ", h=1.5)
  axis(1, at=(0:1)+0.5, labels=c("AL", "AS-AQ"), cex.axis=1.5, cex.lab=1.5)
  axis(2, cex.axis=1.5, cex.lab=1.5)
  mtext("  A", side=3, line=-3, adj=0.95 , cex=1.5)
  #panel 2: profile of prophylaxis
  pAl <- Vectorize(function(x){1-pgamma(x, shape=r[1], rate=r[1]/durP[1])})
  pAsaq <- Vectorize(function(x){1-pgamma(x, shape=r[2], rate=r[2]/durP[2])})
  x <- seq(0,30,by=0.1)
  plot(x, pAsaq(x), ylim=c(0,1), type="l", col="forestgreen", ylab="proportion protected",
       xlab="days after AS-AQ treatment", lwd=3, cex.axis=1.5, cex.lab=1.5)
  for(i in 1:length(siteNos)) {
    pAsaq <- Vectorize(function(x){1-pgamma(x, shape=r[2], rate=r[2]/median(samplesASAQ[,siteNos[i]]))})
    lines(x, pAsaq(x), col="forestgreen", lwd=1,lty=2)
  }
  mtext("  B", side=3, line=-3, adj=0.95 , cex=1.5)
  plot(x, pAl(x), ylim=c(0,1), type="l", col="blue", ylab="proportion protected",
       xlab="days after AL treatment", lwd=3, cex.axis=1.5, cex.lab=1.5)
  for(i in 1:length(siteNos)) {
    pAl <- Vectorize(function(x){1-pgamma(x, shape=r[1], rate=r[1]/median(samplesAL[,siteNos[i]]))})
    lines(x, pAl(x), col="blue", lwd=1,lty=2)
  }
  mtext("  C", side=3, line=-3, adj=0.95 , cex=1.5)
  
  dev.off()
} 


#### Figure 2 Time to reinfection after treatment and model fits	

### get predicted values of prob reinfection at each site by simulation.....
#...or read in the saved predictions (see below)
## function to simulate reinfection for one site
simulate_Lucy <- function(shapei=r[2], durPi, siteEiri, age_y_vec, nReps) {
  # parameters
  rho <- 0.85
  age0 <- 2920 #days
  bh	=	0.590076		#{parameters for the immunity functions}
  bmin	=	0.5
  db	=	3650
  kb	=	2.15506
  ub	=	7.19919
  x_I = 0.4999674  ### x_I only varies by age if widths of age groups vary. Have kept constant
  IB0	= 43.8787
  
  days_to_patent <- durPi + 3.5 # minimum days until patent. Duration of protection + time from emergence from liver to patent infection
  lambdaTranstime <- shapei/days_to_patent # calculate lambda. ok here median should be the mean.. (gamma distribution)
  transTime <- rgamma(nReps,shape=shapei, rate=lambdaTranstime) 
  indEirNow <- siteEiri*(1-rho*exp(-age_y_vec*365/age0)) #actual eir acting on an individual of certain age
  agecat<-as.numeric(cut(age_y_vec,breaks=seq(from=0,to=69,by=0.5),include.lowest=T) )
  # need to iteratively calculate current IB, summing over past exposure (assuming no change in EIR).
  # EIR experienced by person i in first age group
  # IB immunity of the person in first age group.
  indEirAges<-indIBEq<-vector(length=(69*2))
  b<-indFoid<-vector(length=length(age_y_vec))
  indEirAges[1] = siteEiri*(1-rho*exp(-(1/2-0.25)*365/age0))  ## what EIR did they experience when in age group 1?
  indIBEq[1]= (indEirAges[1]/(indEirAges[1]*ub+1) *x_I) / (1+x_I/db) # IB immunity when person i was in age group 1.
  for(j in 2:length(indEirAges)) {   ##compute for all age groups within this site.
    # EIR experienced in earlier age groups.
    indEirAges[j] = siteEiri*(1-rho*exp(-(j/2-0.25)*365/age0)) # j/2-0.25=current age, at midpoint of age gp.
    indIBEq[j] = (indIBEq[j-1] + indEirAges[j]/(indEirAges[j]*ub + 1)*x_I)/(1+x_I/db)
  }
  for(i in 1:length(age_y_vec)) {  ## compute reinfection for each person i
    b[i] = bh*((1-bmin)/(1+(indIBEq[agecat[i]]/IB0)^kb) + bmin)
    indFoid[i] <- (1/365)*indEirNow[i]*b[i] #Eir2Foi
    timeToInf <- rexp(nReps,indFoid[i])  ## possible time to reinfection for person i, after drug protection ends.
    if(i==1) totalTime<-transTime+timeToInf   ## NB transtime is the same for each person, time to reinfection is not.
    if(i>1) totalTime<-c(totalTime, transTime+timeToInf)  ## put together infection times for all people.
  }
  times<-seq(0,42,by=0.1)  
  sim_prop_reinfect<-times*NA
  for(j in 1:length(times)) {
    sim_prop_reinfect[j]<-length(totalTime[totalTime<times[j]])/length(totalTime)
  }
  return(data.frame(days=times, y=sim_prop_reinfect))
}


###### Figure 2 (continued) Time to reinfection after treatment and model fits for each site.
siteNos<-siteInfoOrdered$siteNo
# either simulate (as below) or plot predictions from sims file (see above.)
if(plot2file) {
  tiff(filename = "results/site_fits.tiff", res=300, width=2500, height=1400,compression="lzw")
  layout(matrix(c(1,2:5,1,6:9,1,10:13,0,14,14,14,14),nrow=4,ncol=5,byrow = TRUE),widths=c(1,5,5,5,5),heights=c(5,5,5,2))
  par(mar=c(0,0,0,0))
  nReps<-100  
  plot(0,0, col="white",axes=F)
  text(-0.5,0,"% reinfected",srt=90,cex=2)
  par(mar=c(1,2,1,1))
  y<-summ_data
  y$cumu_prop_reinfected<-y$cumu_prop_reinfected*100
  y$cumu_prop_reinfected_lci <-y$cumu_prop_reinfected_lci*100
  y$cumu_prop_reinfected_uci<-y$cumu_prop_reinfected_uci*100
  for(i in 1:length(siteNos)) {
    x<-subset(y, siteNo==siteNos[i] & drug=="AL")
    par(mgp=c(3,0.5,0))
    plotCI(x$days, x$cumu_prop_reinfected, ui=x$cumu_prop_reinfected_uci,
           li=x$cumu_prop_reinfected_lci, ylim=c(0,100), col="blue", sfrac=0,gap=0, pch=19,las=1,xlim=c(0,42))
    x<-subset(y, siteNo==siteNos[i] & drug=="ASAQ")
    plotCI(x$days, x$cumu_prop_reinfected, ui=x$cumu_prop_reinfected_uci, add=T,
           li=x$cumu_prop_reinfected_lci, ylim=c(0,100), col="forestgreen", sfrac=0,gap=0, pch=19,xlim=c(0,42))
    siteDrugDataTemp<-subset(d, siteNum==siteNos[i] & treatNum==1)
    sim<-simulate_Lucy(r[1], median(samplesAL[,siteNos[i]]), median(estSiteEirs[,siteNos[i]]), 
                       siteDrugDataTemp$ageyears, nReps)
    max_days<-max(x$days)
    lines(sim$days[sim$days<=max_days], 100*sim$y[sim$days<=max_days], ylab="Proportion reinfected", xlab="days", main="", col="blue")
    siteDrugDataTemp<-subset(d, siteNum==siteNos[i] & treatNum==2)
    sim<-simulate_Lucy(r[2], median(samplesASAQ[,siteNos[i]]), median(estSiteEirs[,siteNos[i]]),
                       siteDrugDataTemp$ageyears, nReps)
    lines(sim$days[sim$days<=max_days], 100*sim$y[sim$days<=max_days], col="forestgreen")
    text(1,93, siteInfoOrdered$desc[i], adj=0,cex=1.2)
  }
  par(mar=c(0,0,0,0))
  plot(0,0, col="white",axes=F)
  text(0,-0.3,"days",cex=2)
  dev.off()
}


## Figure 3 Trial-specific EIR estimates
siteInfoOrdered<-siteInfoOrdered[order(siteInfoOrdered$priorEIR),]
if(plot2file) {
  tiff(file="eirEstimates.tiff", compression = "lzw",res=300,width=1800,height=1800)
  plot(NA, NA, xlim=c(1,ncol(samplesASAQ)+1), ylim=c(min(0.5,min(estSiteEirs)),max(estSiteEirs)), type="n", 
     axes=F, xlab=NA, ylab="EIR (per person per year)", cex.axis=1.5, cex.lab=1.5 #)
     , log="y")
  for (i in 1:nrow(siteInfoOrdered)){
    vioplot(estSiteEirs[,siteInfoOrdered$siteNo[i]],  horizontal=F, col="firebrick",#"gold",
  	  at=i+0.5, add=TRUE)
    points(i+0.5, siteInfoOrdered$priorEIR[i], pch=5, col="darkorange", lwd=3)#col="red")
  }
  text((1:nrow(siteInfoOrdered))+0.5, 0.35,srt = 45, adj = 1,
      labels = siteInfoOrdered$desc, xpd = TRUE, cex=0.8)
  labs<-c(0.5,5,50,500)
  axis(2, cex.axis=1.5, cex.lab=1.5, at=labs, labels=as.character(labs))
  legend("topleft", c("estimated EIR", "MAP prediction"), lty=c(1, NA), lwd=3, col=c("firebrick", "darkorange"), #col=c("red", "yellow"), 
  	pch=c(NA, 5), bg="white")
  dev.off()
}


 
###### Figure 4. Duration of protection after treatment with (A & C) AS-AQ and (B & D) AL, according to local pfmdr1 N86Y (A & B) and pfcrt K76T mutation prevalence (C & D). 
if(plot2file) {
  par(mar=c(1,1,1,1))
  tiff(filename = "dur_by_mut.tiff", res=300, width=1800, height=2000,compression="lzw")
  layout(matrix(c(1:4),nrow=2,ncol=2,byrow = TRUE))
  plotCI(mut$prev86*100, mut$dur_p_asaq, xlab="mdr1 86Y prevalence (%)", ylab="days protection after AS-AQ",
         li=mut$lci_dur_p_asaq,ui=mut$uci_dur_p_asaq,gap=0,pch=19,xlim=c(0,100),ylim=c(5,23),sfrac=0, 
         col="forestgreen")
  plotCI(mut$prev86*100,mut$dur_p_asaq,li=mut$prev86_lci*100,
         ui=mut$prev86_uci*100,add=T,err='x',gap=0,sfrac=0,col="forestgreen")
  text(2,22.5 ,"A",cex=1.2)
  
  plotCI(mut$prev86*100, mut$dur_p_al, xlab="mdr1 86Y prevalence (%)", ylab="days protection after AL",
         li=mut$lci_dur_p_al,ui=mut$uci_dur_p_al,gap=0,pch=19,xlim=c(0,100),ylim=c(5,23),sfrac=0,col="blue")
  plotCI(mut$prev86*100,mut$dur_p_al,li=mut$prev86_lci*100,
         ui=mut$prev86_uci*100,add=T,err='x',gap=0,sfrac=0,col="blue")
  text(2,22.5 ,"B",cex=1.2)
  
  plotCI(mut$prev76*100, mut$dur_p_asaq, xlab="crt 76T prevalence (%)", ylab="days protection after AS-AQ",
         li=mut$lci_dur_p_asaq,ui=mut$uci_dur_p_asaq,gap=0,pch=19,xlim=c(0,100),ylim=c(5,23),sfrac=0, 
         col="forestgreen")
  plotCI(mut$prev76*100,mut$dur_p_asaq,li=mut$prev76_lci*100,
         ui=mut$prev76_uci*100,add=T,err='x',gap=0,sfrac=0,col="forestgreen")
  text(2,22.5 ,"C",cex=1.2)
  
  plotCI(mut$prev76*100, mut$dur_p_al, xlab="crt 76T prevalence (%)", ylab="days protection after AL",
         li=mut$lci_dur_p_al,ui=mut$uci_dur_p_al,gap=0,pch=19,xlim=c(0,100),ylim=c(5,23),sfrac=0,col="blue")
  plotCI(mut$prev76*100,mut$dur_p_al,li=mut$prev76_lci*100,
         ui=mut$prev76_uci*100,add=T,err='x',gap=0,sfrac=0,col="blue")
  text(2,22.5 ,"D",cex=1.2)
  
  dev.off()
}


#Figure S1. The duration of post-treatment prophylaxis at different trial locations in order of increasing estimated EIR.
siteInfoOrdered<-siteInfoOrdered[order(siteInfoOrdered$estEIR),]
if(plot2file) {
  tiff(file="durPBySite_IB2018.tiff", compression = "lzw",res=300,width=1500,height=1500)
  plot(NA, NA, xlim=c(1,3*ncol(samplesASAQ)), ylim=c(0,40), type="n", axes=F, xlab=NA, 
       ylab="days of protection after treatment", cex.axis=1.5, cex.lab=1.5)
  abline(v=0, lty=3, lwd=0.5)
  for (i in 1:nrow(siteInfoOrdered)){
    vioplot(samplesAL[, siteInfoOrdered$siteNo[i]] , col="blue", horizontal=F, at=3*i-2, 
            add=TRUE)
    vioplot(samplesASAQ[, siteInfoOrdered$siteNo[i]], col="forestgreen", horizontal=F, 
            at=3*i-1, add=TRUE, names="ASAQ")
    abline(v=3*i, lty=3, lwd=0.5)
  }
  text(3*(1:nrow(siteInfoOrdered))-1.5, par("usr")[3] - 0.5, srt = 45, adj = 1,
       labels = c(as.character(siteInfoOrdered$desc)), xpd = TRUE, cex=0.8)
  axis(2, cex.axis=1.5, cex.lab=1.5)
  legend("top", c("AL", "AS-AQ"), col=c("blue", "forestgreen"), lwd=3, 
         bg="white",cex=0.8)
  dev.off()
} 

### Figure S2. Correlation between pfcrt 76T prevalence and pfmdr 86Y prevalence
if(plot2file) {
  tiff(filename = "corr_mut.tiff", res=300, width=1400, height=1400,compression="lzw")
  plotCI(mut$prev86*100, mut$prev76*100, xlab="86Y prevalence (%)", ylab="76T prevalence (%)",
         li=mut$prev76_lci*100,ui=mut$prev76_uci*100,gap=0,pch=19,xlim=c(0,100),ylim=c(0,100),sfrac=0, 
         col="orange")
  plotCI(mut$prev86*100,mut$prev76*100,li=mut$prev86_lci*100,
         ui=mut$prev86_uci*100,add=T,err='x',gap=0,sfrac=0,col="orange")
  dev.off()
}

