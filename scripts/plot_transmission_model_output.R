# Files required:
# Profile of protection estimated by hidden markov analysis. created here: outputPlotsNProcessingLO_het_IB2018

## write function to calculate cumulative clinical incidence
cumu <- function(x, var, start, end) {
  clin_b<-var[x$year>start & x$year<=end]
  year<-x$year[x$year>start  & x$year<=end]
  cum_clin_b<-0
  for(i in 1:(length(clin_b)-1)) {
    cum_clin_b<- cum_clin_b + clin_b[i]*(year[i+1]-year[i])
  }
  return(cum_clin_b)
}


### Gourcy, Burkina Faso: low mutations, AQ LONG, LUM SHORT
names<-c("aq_long","lum_short")
trans<-c("5","15","50")
s<-c("s","ns")
##matrix for storing results.
summary<- data.frame(seasonal=rep(s,each=6),prev_210=rep(rep(trans,each=2),2),drug=rep(names,6))
summary$cumu_clin_all<-NA
summary$cumu_clin_0_5<-NA
for(i in names) {
  for(j in trans) {
    for(k in s) {
      ind<-which(summary$prev_210==j & summary$drug==i & summary$seasonal==k)
      x<-read.delim(paste0("run_",i,"_",k,"_",j,".txt"))
      summary$cumu_clin_all[ind]<-cumu(x,x$clin_inc_all,-0.001, 5)
      summary$cumu_clin_0_5[ind]<-cumu(x,x$clin_inc_0_5,-0.001, 5)
    }
  }
}

diff_ns_site5<- summary$cumu_clin_all[which(summary$drug=="aq_long" & summary$seasonal=="ns")] -
  summary$cumu_clin_all[which(summary$drug=="lum_short" & summary$seasonal=="ns")]
diff_s_site5<- summary$cumu_clin_all[which(summary$drug=="aq_long" & summary$seasonal=="s")] -
  summary$cumu_clin_all[which(summary$drug=="lum_short" & summary$seasonal=="s")]
diff0_5_ns_site5<- summary$cumu_clin_0_5[which(summary$drug=="aq_long" & summary$seasonal=="ns")] -
  summary$cumu_clin_0_5[which(summary$drug=="lum_short" & summary$seasonal=="ns")]
diff0_5_s_site5<- summary$cumu_clin_0_5[which(summary$drug=="aq_long" & summary$seasonal=="s")] -
  summary$cumu_clin_0_5[which(summary$drug=="lum_short" & summary$seasonal=="s")]
diff_bar_site5<-rbind(diff_ns_site5,(diff_s_site5-diff_ns_site5))
diff_bar_0_5_site5<-rbind(diff0_5_ns_site5,(diff0_5_s_site5-diff0_5_ns_site5))

ratio_ns_site5<- summary$cumu_clin_all[which(summary$drug=="aq_long" & summary$seasonal=="ns")] /
  summary$cumu_clin_all[which(summary$drug=="lum_short" & summary$seasonal=="ns")]
ratio_s_site5<- summary$cumu_clin_all[which(summary$drug=="aq_long" & summary$seasonal=="s")] /
  summary$cumu_clin_all[which(summary$drug=="lum_short" & summary$seasonal=="s")]
ratio0_5_ns_site5<- summary$cumu_clin_0_5[which(summary$drug=="aq_long" & summary$seasonal=="ns")] /
  summary$cumu_clin_0_5[which(summary$drug=="lum_short" & summary$seasonal=="ns")]
ratio0_5_s_site5<- summary$cumu_clin_0_5[which(summary$drug=="aq_long" & summary$seasonal=="s")] /
  summary$cumu_clin_0_5[which(summary$drug=="lum_short" & summary$seasonal=="s")]
ratio_bar_site5<-rbind((ratio_ns_site5-1),(ratio_s_site5-ratio_ns_site5))
ratio_bar_0_5_site5<-rbind((ratio0_5_ns_site5-1),(ratio0_5_s_site5-ratio0_5_ns_site5))

### Nimba, Liberia , aq short, lum long.
names<-c("aq_short","lum_long")
trans<-c("5","15","50")
s<-c("s","ns")
##matrix for storing results.
summary<- data.frame(seasonal=rep(s,each=6),prev_210=rep(rep(trans,each=2),2),drug=rep(names,6))
summary$cumu_clin_all<-NA
summary$cumu_clin_0_5<-NA
for(i in names) {
  for(j in trans) {
    for(k in s) {
      ind<-which(summary$prev_210==j & summary$drug==i & summary$seasonal==k)
      x<-read.delim(paste0("run_",i,"_",k,"_",j,".txt"))
      summary$cumu_clin_all[ind]<-cumu(x,x$clin_inc_all,-0.001, 5)
      summary$cumu_clin_0_5[ind]<-cumu(x,x$clin_inc_0_5,-0.001, 5)
    }
  }
}
diff_ns_site7<- summary$cumu_clin_all[which(summary$drug=="aq_short" & summary$seasonal=="ns")] -
  summary$cumu_clin_all[which(summary$drug=="lum_long" & summary$seasonal=="ns")]
diff_s_site7<- summary$cumu_clin_all[which(summary$drug=="aq_short" & summary$seasonal=="s")] -
  summary$cumu_clin_all[which(summary$drug=="lum_long" & summary$seasonal=="s")]
diff0_5_ns_site7<- summary$cumu_clin_0_5[which(summary$drug=="aq_short" & summary$seasonal=="ns")] -
  summary$cumu_clin_0_5[which(summary$drug=="lum_long" & summary$seasonal=="ns")]
diff0_5_s_site7<- summary$cumu_clin_0_5[which(summary$drug=="aq_short" & summary$seasonal=="s")] -
  summary$cumu_clin_0_5[which(summary$drug=="lum_long" & summary$seasonal=="s")]
diff_bar_site7<-rbind(diff_ns_site7,(diff_s_site7-diff_ns_site7))
diff_bar_0_5_site7<-rbind(diff0_5_ns_site7,(diff0_5_s_site7-diff0_5_ns_site7))

ratio_ns_site7<- summary$cumu_clin_all[which(summary$drug=="aq_short" & summary$seasonal=="ns")] /
  summary$cumu_clin_all[which(summary$drug=="lum_long" & summary$seasonal=="ns")]
ratio_s_site7<- summary$cumu_clin_all[which(summary$drug=="aq_short" & summary$seasonal=="s")] /
  summary$cumu_clin_all[which(summary$drug=="lum_long" & summary$seasonal=="s")]
ratio0_5_ns_site7<- summary$cumu_clin_0_5[which(summary$drug=="aq_short" & summary$seasonal=="ns")] /
  summary$cumu_clin_0_5[which(summary$drug=="lum_long" & summary$seasonal=="ns")]
ratio0_5_s_site7<- summary$cumu_clin_0_5[which(summary$drug=="aq_short" & summary$seasonal=="s")] /
  summary$cumu_clin_0_5[which(summary$drug=="lum_long" & summary$seasonal=="s")]
ratio_bar_site7<-rbind((ratio_ns_site7-1),(ratio_s_site7-ratio_ns_site7))
ratio_bar_0_5_site7<-rbind((ratio0_5_ns_site7-1),(ratio0_5_s_site7-ratio0_5_ns_site7))

##### see code in outputPlotsSemiMarkov for generation of protection profiles.
pp_aq_short<-read.delim("p_protect_aq_short.txt",header=F)
pp_aq_long<-read.delim("p_protect_aq_long.txt",header=F)
pp_lum_short<-read.delim("p_protect_lum_short.txt",header=F)
pp_lum_long<-read.delim("p_protect_lum_long.txt",header=F)
x <- seq(0,30,by=0.1)

##### plot baseline seasonal variation (select the period prior to change in treatment).
x<-read.delim("run_aq_short_s_15.txt")
x<-subset(x, year>=-1 & year<=0)
x$year2<-x$year+1
x$month<-x$year2*12
if(plot2file) {
  tiff(filename = "seasonal_variation.tiff", res=300, width=1200, height=1200,compression="lzw")
  plot(x$month,x$EIR,type="l",ylab="EIR",xlab="month",col="purple")
  dev.off()
}

# FIGURE 5: Duration of prophylaxis and impact on clinical incidence in under 5 year old children of using AS-AQ rather than AL as first-line treatment
if(plot2file) {
  tiff(filename = "transmission1_update2.tiff", res=300, width=3000, height=1800,compression="lzw")
  par(mar=c(5,6.2,1,2),las=1,mgp=c(3.5,1,0))
  layout(matrix(c(1:6),nrow=2,ncol=3,byrow = T))
  # site 5
  axis_cex=1.5
  letter_cex=2
  plot(pp_aq_long$V1, pp_aq_long$V4, ylim=c(0,1), type="l", col="forestgreen", ylab="proportion protected",
       xlab="days after treatment", lwd=3, cex.axis=1.5, cex.lab=1.5)
  lines(c(pp_lum_short$V1,seq(max(pp_lum_short$V1),30,0.2)),
        c(pp_lum_short$V4,rep(0,length(seq(max(pp_lum_short$V1),30,0.2)))),col="blue",lwd=3)
  text(28,0.97,"A",cex=letter_cex)
  legend(19,0.85,c("AL","AS-AQ"),col=c("blue","forestgreen"),lty=c(1,1),lwd=c(2,2))
  axis_label<- c("low", "medium", "high")
  barplot(diff_bar_0_5_site5, col=c("orange","red1"), names.arg=rep("",3), 
          ylab="Difference in clinical \n episodes / child after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(-1.6,0))
  lines(0:3, rep(0,4))
  text(0.4,-1.5,"B",cex=letter_cex)
  legend(1.7,-1.6,legend=c("non-seasonal","seasonal"), fill=c("orange","red"),xpd=T,cex=1.5)
  barplot(ratio_bar_0_5_site5*100, col=c("orange","red1"), names.arg=rep("",3), 
          ylab="% clinical episodes\n prevented / child after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(-15,0))
  lines(0:3, rep(0,4))
  text(0.4,-14,"C",cex=letter_cex)
   # site 7
  plot(pp_aq_short$V1, pp_aq_short$V4, ylim=c(0,1), type="l", col="forestgreen", ylab="proportion protected",
       xlab="days after treatment", lwd=3, cex.axis=1.5, cex.lab=1.5,xlim=c(0,30))
  lines(c(pp_lum_long$V1,seq(max(pp_lum_long$V1),30,0.2)),
        c(pp_lum_long$V4,rep(0,length(seq(max(pp_lum_long$V1),30,0.2)))),col="blue",lwd=3)
  text(28,0.97,"D",cex=letter_cex)
  barplot(diff_bar_0_5_site7, col=c("orange","red1"), names.arg=axis_label, 
          ylab="Difference in clinical \n episodes / child after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(0,1.6))
  lines(0:3, rep(0,4))  
  text(0.4,1.5,"E",cex=letter_cex)
  barplot(ratio_bar_0_5_site7*100, col=c("orange","red1"), names.arg=axis_label, 
          ylab="% clinical episodes\n prevented / child after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(0,15))
  lines(0:3, rep(0,4))
  text(0.4,14,"F",cex=letter_cex)
  dev.off()
}

##### Figure S4. CASE REDUCTIONS, ALL AGES.
if(plot2file) {
  tiff(filename = "results/transmission_supp_updated.tiff", res=300, width=3000, height=1800,compression="lzw")
  par(mar=c(5,6.2,1,2),las=1,mgp=c(3.5,1,0))
  layout(matrix(c(1:6),nrow=2,ncol=3,byrow = T))
  # site 5
  axis_cex=1.5
  letter_cex=2
  plot(pp_aq_long$V1, pp_aq_long$V4, ylim=c(0,1), type="l", col="forestgreen", ylab="proportion protected",
       xlab="days after treatment", lwd=3, cex.axis=1.5, cex.lab=1.5)
  lines(c(pp_lum_short$V1,seq(max(pp_lum_short$V1),30,0.2)),
        c(pp_lum_short$V4,rep(0,length(seq(max(pp_lum_short$V1),30,0.2)))),col="blue",lwd=3)
  text(28,0.97,"A",cex=letter_cex)
  legend(19,0.85,c("AL","AS-AQ"),col=c("blue","forestgreen"),lty=c(1,1),lwd=c(2,2))
  axis_label<- c("low", "medium", "high")
  barplot(diff_bar_site5, col=c("orange","red1"), names.arg=rep("",3), 
          ylab="Difference in clinical \n episodes / person after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(-1.6,0))
  lines(0:3, rep(0,4))
  text(0.4,-1.5,"B",cex=letter_cex)
  legend(1.7,-1.6,legend=c("non-seasonal","seasonal"), fill=c("orange","red"),xpd=T,cex=1.5)
  barplot(ratio_bar_site5*100, col=c("orange","red1"), names.arg=rep("",3), 
          ylab="% clinical episodes\n prevented / person after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(-15,0))
  lines(0:3, rep(0,4))
  text(0.4,-14,"C",cex=letter_cex)
  # site 7
  plot(pp_aq_short$V1, pp_aq_short$V4, ylim=c(0,1), type="l", col="forestgreen", ylab="proportion protected",
       xlab="days after treatment", lwd=3, cex.axis=1.5, cex.lab=1.5,xlim=c(0,30))
  lines(c(pp_lum_long$V1,seq(max(pp_lum_long$V1),30,0.2)),
        c(pp_lum_long$V4,rep(0,length(seq(max(pp_lum_long$V1),30,0.2)))),col="blue",lwd=3)
  text(28,0.97,"D",cex=letter_cex)
  barplot(diff_bar_site7, col=c("orange","red1"), names.arg=axis_label, 
          ylab="Difference in clinical \n episodes / person after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(0,1.6))
  lines(0:3, rep(0,4))  
  text(0.4,1.5,"E",cex=letter_cex)
  barplot(ratio_bar_site7*100, col=c("orange","red1"), names.arg=axis_label, 
          ylab="% clinical episodes\n prevented / person after 5yrs",
          cex.names=1.8, cex.lab=1.5, cex.axis=axis_cex, ylim=c(0,15))
  lines(0:3, rep(0,4))
  text(0.4,14,"F",cex=letter_cex)
  dev.off()
}


