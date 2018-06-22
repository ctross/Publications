fitDem <- result1[[3]]
 library(rethinking)
 library(Cairo)
 
  Samples <- ((Iter-Warmup)*Chains)
  
###################################################################### Check Estimates
 print(fitDem,digits_summary=3,pars=c("Beta_males","Beta_females"))
 print(fitDem,digits_summary=3,pars=c("Delta_males","Delta_females"))
 print(fitDem,digits_summary=3,pars=c("Phi_males","Phi_females"))
 print(fitDem,digits_summary=3,pars=c("Theta_males","Theta_females"))
 print(fitDem,digits_summary=3,pars=c("Zeta_males","Zeta_females"))
 print(fitDem,digits_summary=3,pars=c("Gamma_males","Gamma_females"))

##################################################################### Plot Densities
 BM<-extract(fitDem,pars="Beta_males")$Beta_males
 BF<-extract(fitDem,pars="Beta_females")$Beta_females
 x <- data.frame(Male=exp(BM[,4]),Female=exp(BF[,4]))
 library(ggplot2);library(reshape2);library(Cairo)
 data<- melt(x)
 p1<-ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.4, colour=NA)

 q5a <- quantile(x$Male,.05)
 q95a <- quantile(x$Male,.95)
 medxa <- median(x$Male)
 x.densa <- density(x$Male)
 df.densa <- data.frame(x = x.densa$x, y = x.densa$y)
 p1<-p1 + geom_area(data = subset(df.densa, x >= q5a & x <= q95a),
 aes(x=x,y=y), fill = '#253494',alpha=0.4)
 q5b <- quantile(x$Female,.05)
 q95b <- quantile(x$Female,.95)
 medxb <- median(x$Female)
 x.densb <- density(x$Female)
 df.densb <- data.frame(x = x.densb$x, y = x.densb$y)
 p1<-p1 + geom_area(data = subset(df.densb, x >= q5b & x <= q95b),
 aes(x=x,y=y), fill = '#e31a1c',alpha=0.7)   + ylab("Density")+xlab("Elasticity: Timing- and Quality-\nWeighted Years Married")+
 scale_fill_manual("Sex", values = c("#253494","#e31a1c"))+geom_vline(xintercept = 1) +
 theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),
 legend.text=element_text(size=12), legend.title=element_text(size=14,face="bold"))
# CairoPDF("ElastSY",height=8,width=8)
  p1
# dev.off()

############################################################ Self Quality Curves
SAE<-extract(fitDem,pars="MarriageAgeEffect_males")$MarriageAgeEffect_males
Delta<-extract(fitDem,pars="Delta_males")$Delta_males

MaxExposure<-60
Res1<-SAE
for( j in 1:MaxExposure){
for(i in 1:Samples){
Res1[i,j]<-logistic(Delta[i,1] + SAE[i,j])
}
}

Res1M2m<-Res1M1m<-Res1Mm<-c()

for(j in 1:MaxExposure){
Res1M1m[j]<-HPDI(Res1[,j],0.65)[1]
Res1M2m[j]<-HPDI(Res1[,j],0.65)[2]
Res1Mm[j]<-median(Res1[,j])
}


SAE<-extract(fitDem,pars="MarriageAgeEffect_females")$MarriageAgeEffect_females
Delta<-extract(fitDem,pars="Delta_females")$Delta_females


for( j in 1:MaxExposure){
for(i in 1:Samples){
Res1[i,j]<-logistic(Delta[i,1] +  SAE[i,j])
}
}

Res1M2f<-Res1M1f<-Res1Mf<-c()

for(j in 1:MaxExposure){
Res1M1f[j]<-HPDI(Res1[,j],0.65)[1]
Res1M2f[j]<-HPDI(Res1[,j],0.65)[2]
Res1Mf[j]<-median(Res1[,j])
}

spr<-0.65
plot1<-data.frame(Res1M1f=smooth.spline(Res1M1f,spar=spr)$y,Res1M2f=smooth.spline(Res1M2f,spar=spr)$y,
Res1Mf=smooth.spline(Res1Mf,spar=spr)$y,
Res1M1m=smooth.spline(Res1M1m,spar=spr)$y,Res1M2m=smooth.spline(Res1M2m,spar=spr)$y,
Res1Mm=smooth.spline(Res1Mm,spar=spr)$y,
Age=c(1:MaxExposure)+11)

plot1$Res1M1m <-ifelse(plot1$Res1M1m<0,0, plot1$Res1M1m)
plot1$Res1M2m <-ifelse(plot1$Res1M2m>1,1, plot1$Res1M2m)

plot1$Res1M1f <-ifelse(plot1$Res1M1f<0,0, plot1$Res1M1f)
plot1$Res1M2f <-ifelse(plot1$Res1M2f>1,1, plot1$Res1M2f)

p2 <- ggplot(plot1,aes(x=Age,y=Res1Mm)) +  ylab("Relative Value of Marriage to RS of Ego")+ xlab("Ego Age")+
      geom_ribbon(aes(ymin=Res1M1m,ymax=Res1M2m),alpha=0.5,fill="#253494")+
      geom_line(aes(x=Age,y=Res1Mm),color="#253494",size=2)+
      theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"))

p3 <- ggplot(plot1,aes(x=Age,y=Res1M1f)) +  ylab("Relative Value of Marriage to RS of Ego")+
      geom_ribbon(aes(ymin=Res1M1f,ymax=Res1M2f),alpha=0.5,fill="#e31a1c")+  xlab("Ego Age")+
      geom_line(aes(x=Age,y=Res1Mf),color="#e31a1c",size=2)+
      theme(axis.text=element_text(size=20), axis.title=element_text(size=22,face="bold"))

ggsave("SelfM.pdf",p2,height=8,width=8)

ggsave("SelfF.pdf",p3,height=8,width=8)



############################################################ Spouse Quality Curves
SAE_m<-extract(fitDem,pars="SpouseAgeEffect_males")$SpouseAgeEffect_males
Phi<-extract(fitDem,pars="Phi_males")$Phi_males

Res1<-SAE_m
for( j in 1:MaxExposure){
for(i in 1:Samples){
Res1[i,j]<-logistic(Phi[i,2] + SAE_m[i,j])
}
}

Res1M2m<-Res1M1m<-Res1Mm<-c()

for(j in 1:MaxExposure){
Res1M1m[j]<-HPDI(Res1[,j],0.65)[1]
Res1M2m[j]<-HPDI(Res1[,j],0.65)[2]
Res1Mm[j]<-median(Res1[,j])
}

SAE_f<-extract(fitDem,pars="SpouseAgeEffect_females")$SpouseAgeEffect_females
Phi<-extract(fitDem,pars="Phi_females")$Phi_females

Res1<-SAE_f
for( j in 1:MaxExposure){
for(i in 1:Samples){
Res1[i,j]<-logistic(Phi[i,2] + SAE_f[i,j])
}
}

Res1M2f<-Res1M1f<-Res1Mf<-c()

for(j in 1:MaxExposure){
Res1M1f[j]<-HPDI(Res1[,j],0.65)[1]
Res1M2f[j]<-HPDI(Res1[,j],0.65)[2]
Res1Mf[j]<-median(Res1[,j])
}

spr<-0.65

plot1<-data.frame(Res1M1f=smooth.spline(Res1M1f,spar=spr)$y,Res1M2f=smooth.spline(Res1M2f,spar=spr)$y,
Res1Mf=smooth.spline(Res1Mf,spar=spr)$y,
Res1M1m=smooth.spline(Res1M1m,spar=spr)$y,Res1M2m=smooth.spline(Res1M2m,spar=spr)$y,
Res1Mm=smooth.spline(Res1Mm,spar=spr)$y,
Age=c(1:MaxExposure)+11)

plot1$Res1M1m <-ifelse(plot1$Res1M1m<0,0, plot1$Res1M1m)
plot1$Res1M2m <-ifelse(plot1$Res1M2m>1,1, plot1$Res1M2m)

plot1$Res1M1f <-ifelse(plot1$Res1M1f<0,0, plot1$Res1M1f)
plot1$Res1M2f <-ifelse(plot1$Res1M2f>1,1, plot1$Res1M2f)

p4 <- ggplot(plot1,aes(x=Age,y=Res1M1m)) +  ylab("Relative Value of Spouse")+ xlab("Age")+
      geom_ribbon(aes(ymin=Res1M1m,ymax=Res1M2m),alpha=0.5,fill="#253494")+
      geom_line(aes(x=Age,y=Res1Mm),color="#253494",size=2)+
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

p5 <- ggplot(plot1,aes(x=Age,y=Res1M1f)) +  ylab("Relative Value of Spouse")+
      geom_ribbon(aes(ymin=Res1M1f,ymax=Res1M2f),alpha=0.5,fill="#e31a1c")+  xlab("Age")+
      geom_line(aes(x=Age,y=Res1Mf),color="#e31a1c",size=2)+
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

#CairoPDF("SpouseM",height=8,width=8)
  p4
# dev.off()

#  CairoPDF("SpouseF",height=8,width=8)
  p5
# dev.off()

#################################################################################################
########################################################### Zero Curves
Ceta_m<-extract(fitDem,pars="Ceta_males")$Ceta_males
Fage <- 1:45

Res0m<-matrix(NA,nrow=Samples,ncol=45)
for( j in 1:45){
for(i in 1:Samples){
Res0m[i,j]<-1-logistic(Ceta_m[i,1] + Ceta_m[i,2]*log(Fage[j]))
}
}

Res0mM<-Res0mH<-Res0mL<-c()

for(j in 1:45){
Res0mM[j]<-median(Res0m[,j])
Res0mH[j]<-HPDI(Res0m[,j],0.9)[2]
Res0mL[j]<-HPDI(Res0m[,j],0.9)[1]
}

Ceta_f<-extract(fitDem,pars="Ceta_females")$Ceta_females
Fage <- 1:45

Res0f<-matrix(NA,nrow=Samples,ncol=45)
for( j in 1:45){
for(i in 1:Samples){
Res0f[i,j]<-1-logistic(Ceta_f[i,1] + Ceta_f[i,2]*log(Fage[j]))
}
}

Res0fM<-Res0fH<-Res0fL<-c()

for(j in 1:45){
Res0fM[j]<-median(Res0f[,j])
Res0fH[j]<-HPDI(Res0f[,j],0.9)[2]
Res0fL[j]<-HPDI(Res0f[,j],0.9)[1]
}

plot1<-data.frame(Median=c(Res0mM,Res0fM),Low=c(Res0mL,Res0fL),High=c(Res0mH,Res0fH),
Age=rep(Fage+11,2),Sex=factor(c(rep("Male",45),rep("Female",45))))


p6 <- ggplot(plot1,aes(x=Age,y=Median)) +  ylab("Probability Unmarried")+ xlab("Age")+ geom_line(aes(x=Age,y=Median, colour=Sex),size=1.25)+
    geom_ribbon(aes(ymin=Low,ymax=High,fill=Sex),alpha=0.5) + scale_fill_manual(values=c("#e31a1c","#253494")) +  scale_colour_manual(values=c("#e31a1c","#253494")) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

#CairoPDF("Unmarried",height=8,width=8)
  p6
# dev.off()


########################################################### RS Unmarried
Deta_m<-extract(fitDem,pars="Deta_males")$Deta_males
Fage <- 1:45

Res0m<-matrix(NA,nrow=Samples,ncol=45)
for( j in 1:45){
for(i in 1:Samples){
Res0m[i,j]<-exp(Deta_m[i,1] + Deta_m[i,2]*log(Fage[j]))
}
}

Res0mM<-Res0mH<-Res0mL<-c()

for(j in 1:45){
Res0mM[j]<-median(Res0m[,j])
Res0mH[j]<-HPDI(Res0m[,j],0.9)[2]
Res0mL[j]<-HPDI(Res0m[,j],0.9)[1]
}

Deta_f<-extract(fitDem,pars="Deta_females")$Deta_females
Fage <- 1:45

Res0f<-matrix(NA,nrow=Samples,ncol=45)
for( j in 1:45){
for(i in 1:Samples){
Res0f[i,j]<-exp(Deta_f[i,1] + Deta_f[i,2]*log(Fage[j]))
}
}

Res0fM<-Res0fH<-Res0fL<-c()

for(j in 1:45){
Res0fM[j]<-median(Res0f[,j])
Res0fH[j]<-HPDI(Res0f[,j],0.9)[2]
Res0fL[j]<-HPDI(Res0f[,j],0.9)[1]
}

plot1<-data.frame(Median=c(Res0mM,Res0fM),Low=c(Res0mL,Res0fL),High=c(Res0mH,Res0fH),
Age=rep(Fage+11,2),Sex=factor(c(rep("Male",45),rep("Female",45))))


p7 <- ggplot(plot1,aes(x=Age,y=Median)) +  ylab("Reproductive Success")+ xlab("Age")+ geom_line(aes(x=Age,y=Median, colour=Sex))+
    geom_ribbon(aes(ymin=Low,ymax=High,fill=Sex),alpha=0.5) + scale_fill_manual(values=c("#e31a1c","#253494")) +  scale_colour_manual(values=c("#e31a1c","#253494")) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

#CairoPDF("RSUnmarried",height=8,width=8)
  p7
# dev.off()


################################################################################
library(rethinking)
MaxExposure<-86

################################################################################################################################
################################################################################## Effective Partner Years - Males - Mate Timing
################################################################################################################################
MarriageAgeEffect_males<-extract(fitDem,pars="MarriageAgeEffect_males")$MarriageAgeEffect_males
Delta_males<-extract(fitDem,pars="Delta_males")$Delta_males

SQ_SpouseYears_males<-array(NA,c(N_males,Samples))

for(k in 1:Samples){
 for (i in 1:N_males){
   Ticker3 <- 0;

 for (a in 1:Age_males[i]){
  Ticker3 <- Ticker3 + (ifelse(PartnersByAge1_males[i,a+11]==0,0,1)+ifelse(PartnersByAge2_males[i,a+11]==0,0,1)+ifelse(PartnersByAge3_males[i,a+11]==0,0,1)+ifelse(PartnersByAge4_males[i,a+11]==0,0,1)+ifelse(PartnersByAge5_males[i,a+11]==0,0,1))*inv_logit(Delta_males[k,2] +  MarriageAgeEffect_males[k,a])
  }

 SQ_SpouseYears_males[i,k] <- Ticker3;
   }
   }

SQ_SpouseYears_males_median  <- c()
for(i in 1:N_males){
SQ_SpouseYears_males_median[i]<-median(SQ_SpouseYears_males[i,])
   }

################################################################################################################################
################################################################################ Effective Partner Years - Females - Mate Timing
################################################################################################################################
MarriageAgeEffect_females<-extract(fitDem,pars="MarriageAgeEffect_females")$MarriageAgeEffect_females
Delta_females<-extract(fitDem,pars="Delta_females")$Delta_females

SQ_SpouseYears_females<-array(NA,c(N_females,Samples))

for(k in 1:Samples){
 for (i in 1:N_females){
   Ticker3 <- 0;

 for (a in 1:Age_females[i]){
   Ticker3 <- Ticker3 + ifelse(PartnersByAge1_females[i,a+11]==0,0,1)*inv_logit(Delta_females[k,2]  + MarriageAgeEffect_females[k,a]);
   }
   SQ_SpouseYears_females[i,k] <- Ticker3;
   }
   }

SQ_SpouseYears_females_median  <- c()
for(i in 1:N_females){
SQ_SpouseYears_females_median[i]<-median(SQ_SpouseYears_females[i,])
   }


###################################################################### Print Estimates
library(xtable)
library(SkewCalc)
library(rethinking)
cv<-function(x){var(x)/(mean(x)^2)}

 BM<-extract(fitDem,pars="Beta_males")$Beta_males
 BF<-extract(fitDem,pars="Beta_females")$Beta_females

 Irs.m<-c()
 Irs.f<-c()
   for(i in 1:Samples){
  Irs.m[i]<- cv(sample(RS_males[which(Age_males >45)], replace = TRUE))
  Irs.f[i]<- cv(sample(RS_females[which(Age_females >45)], replace = TRUE))
                              }

 Mrs.m<-c()
 Mrs.f<-c()
   for(i in 1:Samples){
     xm1 <- sample(c(1:length(RS_males)), replace = TRUE)
     xf1 <- sample(c(1:length(RS_females)), replace = TRUE)
   Mrs.m[i]<- (M_index(RS_males[xm1],Age_males[xm1])/sqrt(mean(RS_males[xm1])))^2
   Mrs.f[i]<- (M_index(RS_females[xf1],Age_females[xf1])/sqrt(mean(RS_females[xf1])))^2
                              }

 Ims.m<-c()
 Ims.f<-c()
      for(i in 1:Samples){
 Ims.m[i]<- cv(sample(SQ_SpouseYears_males[which(Age_males >45),i], replace = TRUE))
 Ims.f[i]<- cv(sample(SQ_SpouseYears_females[which(Age_females >45),i], replace = TRUE))
                              }

 Mms.m<-c()
 Mms.f<-c()
   for(i in 1:Samples){
     xm1 <- sample(c(1:length(RS_males)), replace = TRUE)
     xf1 <- sample(c(1:length(RS_females)), replace = TRUE)
   Mms.m[i]<- (M_index(SQ_SpouseYears_males[xm1,i],Age_males[xm1])/sqrt(mean(SQ_SpouseYears_males[xm1,i])))^2
   Mms.f[i]<- (M_index(SQ_SpouseYears_females[xf1,i],Age_females[xf1])/sqrt(mean(SQ_SpouseYears_females[xf1,i])))^2
                              }

 parM <-   round(c(
  c(median(Irs.m),HPDI(Irs.m,0.95)), # I RS Male
  c(median(Mrs.m),HPDI(Mrs.m,0.95)), # M RS Male

  c(median(Ims.m),HPDI(Ims.m,0.95)), # I MS Male
  c(median(Mms.m),HPDI(Mms.m,0.95)), # M MS Male

  c(median(BM[,1]),HPDI(BM[,1],0.95)),     # Intercept Male
  c(median(BM[,2]),HPDI(BM[,2],0.95)),     # Age Male
  c(median(BM[,3]),HPDI(BM[,3],0.95)),     # MS Male
  c(median(exp(BM[,4])),HPDI(exp(BM[,4]),0.95))      # MS Male
  ),2)

 parM <- as.vector(parM)


 parF <-   round(c(
  c(median(Irs.f),HPDI(Irs.f,0.95)), # I RS Female
  c(median(Mrs.f),HPDI(Mrs.f,0.95)), # M RS Female

  c(median(Ims.f),HPDI(Ims.f,0.95)), # I MS Female
  c(median(Mms.f),HPDI(Mms.f,0.95)), # M MS Female

  c(median(BF[,1]),HPDI(BF[,1],0.95)),
  c(median(BF[,2]),HPDI(BF[,2],0.95)),
  c(median(BF[,3]),HPDI(BF[,3],0.95)),
  c(median(exp(BF[,4])),HPDI(exp(BF[,4]),0.95))
     ),2)
  parF <-  as.vector(parF)


   parD <- round(c(
  c(median(log(Irs.m/Irs.f)),HPDI(log(Irs.m/Irs.f),0.95)), # dI RS
  c(median(log(Mrs.m/Mrs.f)),HPDI(log(Mrs.m/Mrs.f),0.95)), # dM RS

  c(median(log(Ims.m/Ims.f)),HPDI(log(Ims.m/Ims.f),0.95)), # dI RS
  c(median(log(Mms.m/Mms.f)),HPDI(log(Mms.m/Mms.f),0.95)), # dM RS

  c(median(BM[,1]-BF[,1]),HPDI(BM[,1]-BF[,1],0.95)),
  c(median(BM[,2]-BF[,2]),HPDI(BM[,2]-BF[,2],0.95)),
  c(median(BM[,3]-BF[,3]),HPDI(BM[,3]-BF[,3],0.95)),
  c(median(exp(BM[,4])-exp(BF[,4])),HPDI(exp(BM[,4])-exp(BF[,4]),0.95))
     ),2)
  parD <-  as.vector(parD)



########################################################################## WAIC
Beta_males<-extract(fitDem,pars="Beta_males")$Beta_males
Sigma_males<-extract(fitDem,pars="Sigma_males")$Sigma_males

Beta_females<-extract(fitDem,pars="Beta_females")$Beta_females
Sigma_females<-extract(fitDem,pars="Sigma_females")$Sigma_females

A_males<-array(NA,c(N_males,Samples))
B_males<-array(NA,c(N_males,Samples))

A_females<-array(NA,c(N_females,Samples))
B_females<-array(NA,c(N_females,Samples))

###### Male Link Functions
for(k in 1:Samples){
 for (i in 1:N_males){
  if(Zeros_males[i]==0){
   B_males[i,k]<-1/Sigma_males[k,1];
   A_males[i,k] <- exp(Beta_males[k,1] + Beta_males[k,2]*log(Age_males[i]) + exp(Beta_males[k,4])*log(SQ_SpouseYears_males[i,k]));
   } else{
   B_males[i,k]<- -1
   A_males[i,k] <- -1
   }
   }  }

###### Female Link Functions
for(k in 1:Samples){
 for (i in 1:N_females){
 if(Zeros_females[i]==0){
   B_females[i,k] <- 1/Sigma_females[k,1];
   A_females[i,k] <- exp(Beta_females[k,1] + Beta_females[k,2]*log(Age_females[i]) + exp(Beta_females[k,4])*log(SQ_SpouseYears_females[i,k]));
    } else{
   B_females[i,k] <- -1
   A_females[i,k] <- -1
   }
   } }


# Set Range of Samples to Use for WAIC Computation
Q1 <-1
Q <- Samples

# Then for each data point, for each sample, compute log-likelihood
# pD is sum over data of variance of log-lik
pD <- 0

for ( i in 1:N_males ) {
 if(  A_males[i,1]>0){
    ll<-c()
# compute variance of log-lik for point i over posterior
 for(q in Q1:Q){
    ll[q] <- dgampois(RS_males[i], A_males[i,q],(1/B_males[i,q]), log=TRUE )
               }
    pD <- pD + var(ll,na.rm=T)
}}

# Then for each data point, for each sample, compute log of average likelihood
# lppd is log posterior likelihood
lppd <- 0
for ( i in 1:N_males ){
 if(  A_males[i,1]>0){
  ll<-c()
    # compute log of average likelihood
    for(q in Q1:Q){
    ll[q] <- dgampois(RS_males[i], A_males[i,q],(1/B_males[i,q]), log=FALSE )
                 }
    lppd <- lppd + log(mean(ll,na.rm=T))
}  }

# WAIC
WAIC<- -2*(lppd - pD)

# Return Results
waicM<-round(cbind(pD,lppd,WAIC),2)


# Set Range of Samples to Use for WAIC Computation
Q1 <-1
Q <- Samples

# Then for each data point, for each sample, compute log-likelihood
# pD is sum over data of variance of log-lik
pD <- 0

for ( i in 1:N_females ) {
 if(  A_females[i,1]>0){
    ll<-c()
# compute variance of log-lik for point i over posterior
 for(q in Q1:Q){
    ll[q] <- dgampois(RS_females[i], A_females[i,q],(1/B_females[i,q]), log=TRUE )
               }
    pD <- pD + var(ll,na.rm=T)
}}

# Then for each data point, for each sample, compute log of average likelihood
# lppd is log posterior likelihood
lppd <- 0
for ( i in 1:N_females ){
 if(  A_females[i,1]>0){
  ll<-c()
    # compute log of average likelihood
    for(q in Q1:Q){
    ll[q] <- dgampois(RS_females[i], A_females[i,q],(1/B_females[i,q]), log=FALSE )
                 }
    lppd <- lppd + log(mean(ll,na.rm=T))
}  }

# WAIC
WAIC<- -2*(lppd - pD)

# Return Results
waicF<-round(cbind(pD,lppd,WAIC),2)

ResTab <- rbind(c(parM,waicM),c(parF,waicF), c(parD,rep(NA,3)))
colnames(ResTab) <- c("Irs","IrsL","IrsH","Mrs","MrsL","MrsH","Ims","ImsL","ImsH","Mms","MmsL","MmsH",
                       "Incpt","IncptL","IncptH","Age","AgeL","AgeH","Spouses","SpousesL","SpousesH",
                       "SpouseYears","SpouseYearsL","SpouseYearsH","pD","lppd","WAIC" )
rownames(ResTab) <- c("Males","Females","Delta")

write.csv(ResTab,"Bateman-TWSY.csv")

