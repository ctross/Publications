################################# Load Data ####################################################################################
library(rstan)
library(rethinking)

 g<-read.csv(file.choose())    # MapFileData-WithCountyResultsAndCovariates.csv
 goog<-read.csv(file.choose()) # RacismData_Google-Stephens-Davidowitz.csv

################################################################ Extract Data
     Ym<- g$m.log.RR_Black_Unarmed_Versus_White_Unarmed     # First Outcome
     Ysd<- g$sd.log.RR_Black_Unarmed_Versus_White_Unarmed   #

############################################################### Define race-specific crime rates
     WhiteAssault <- (g$AssaultsWhite.sum/g$WA_TOT)
     BlackAssault <- (g$AssaultsBlack.sum/g$BAC_TOT)

     WhiteWeapons <- (g$WeaponsWhite.sum/g$WA_TOT)
     BlackWeapons <- (g$WeaponsBlack.sum/g$BAC_TOT)

############################################################## Extract other data
     Wealth <- (g$Median.Income)

     Pop<-g$TOT_POP

     BlackRatio<-(g$BAC_TOT+1)/Pop

     Gini<-g$Gini

     State<-g$State

     DMA<-g$DMA

################################################################### Compile Data
#### then reduce to only data with estimates of racial bias in police shooting
     g2<- data.frame(g$County.FIPS.Code,Ym,Ysd, Gini, Wealth,Pop,State,BlackRatio,Gini,DMA,WhiteAssault,BlackAssault,WhiteWeapons,BlackWeapons)
     g3 <- g2[complete.cases(g2$Ym),]

     Ym  <- g3$Ym
     Ysd <- g3$Ysd

     N<-length(Ym)

     DMA <- g3$DMA

     Pop <-g3$Pop/sd(g3$Pop,na.rm=T)
     BlackRatio<-g3$BlackRatio

     Wealth <-g3$Wealth/sd(g3$Wealth,na.rm=T)
     Gini<-g3$Gini

#### Convert measures to similar units (deviations), and prep missing data for
# Bayesian imputation if needed.

     WhiteAssault <-g3$WhiteAssault/sd(g3$WhiteAssault,na.rm=T)
     MaxWhiteAssault<-max(WhiteAssault,na.rm=T)
     MinWhiteAssault<-min(WhiteAssault,na.rm=T)

     BlackAssault <-g3$BlackAssault/sd(g3$BlackAssault,na.rm=T)
     MaxBlackAssault<-max(BlackAssault,na.rm=T)
     MinBlackAssault<-min(BlackAssault,na.rm=T)

     WhiteWeapons <-g3$WhiteWeapons/sd(g3$WhiteWeapons,na.rm=T)
     MaxWhiteWeapons<-max(WhiteWeapons,na.rm=T)
     MinWhiteWeapons<-min(WhiteWeapons,na.rm=T)

     BlackWeapons <-g3$BlackWeapons/sd(g3$BlackWeapons,na.rm=T)
     MaxBlackWeapons<-max(BlackWeapons,na.rm=T)
     MinBlackWeapons<-min(BlackWeapons,na.rm=T)

     WhiteAssault <- ifelse(WhiteAssault==0,NA,WhiteAssault)
     MissCumSumWhiteAssault <-cumsum(is.na(WhiteAssault))
     MissCumSumWhiteAssault <-ifelse(MissCumSumWhiteAssault ==0,1,MissCumSumWhiteAssault )
     NonMissWhiteAssault  <-ifelse(is.na(WhiteAssault ),0,1)
     NmissWhiteAssault <-sum(is.na(WhiteAssault ))
     WhiteAssault [is.na(WhiteAssault )]<-9999999

        BlackAssault <- ifelse(BlackAssault==0,NA,BlackAssault)
     MissCumSumBlackAssault <-cumsum(is.na(BlackAssault))
     MissCumSumBlackAssault <-ifelse(MissCumSumBlackAssault ==0,1,MissCumSumBlackAssault )
     NonMissBlackAssault  <-ifelse(is.na(BlackAssault ),0,1)
     NmissBlackAssault <-sum(is.na(BlackAssault ))
     BlackAssault [is.na(BlackAssault )]<-9999999

         WhiteWeapons <- ifelse(WhiteWeapons==0,NA,WhiteWeapons)
     MissCumSumWhiteWeapons <-cumsum(is.na(WhiteWeapons))
     MissCumSumWhiteWeapons <-ifelse(MissCumSumWhiteWeapons ==0,1,MissCumSumWhiteWeapons )
     NonMissWhiteWeapons  <-ifelse(is.na(WhiteWeapons ),0,1)
     NmissWhiteWeapons <-sum(is.na(WhiteWeapons ))
     WhiteWeapons [is.na(WhiteWeapons )]<-9999999

        BlackWeapons <- ifelse(BlackWeapons==0,NA,BlackWeapons)
     MissCumSumBlackWeapons <-cumsum(is.na(BlackWeapons))
     MissCumSumBlackWeapons <-ifelse(MissCumSumBlackWeapons ==0,1,MissCumSumBlackWeapons )
     NonMissBlackWeapons  <-ifelse(is.na(BlackWeapons ),0,1)
     NmissBlackWeapons <-sum(is.na(BlackWeapons ))
     BlackWeapons [is.na(BlackWeapons )]<-9999999

     Ones<-rep(1,N)

     DMAIndex<-goog$dmaindex
     GoogleRacism<-goog$raciallychargedsearch


model_dat  <-list(N=N,
   Ym=Ym,
   Ysd=Ysd,

   MissCumSumWhiteAssault=MissCumSumWhiteAssault,
   NonMissWhiteAssault=NonMissWhiteAssault,
   NmissWhiteAssault=NmissWhiteAssault,
   WhiteAssault=WhiteAssault,
   MaxWhiteAssault=MaxWhiteAssault,
   MinWhiteAssault=MinWhiteAssault,

   MissCumSumBlackAssault=MissCumSumBlackAssault,
   NonMissBlackAssault=NonMissBlackAssault,
   NmissBlackAssault=NmissBlackAssault,
   BlackAssault=BlackAssault,
   MaxBlackAssault=MaxBlackAssault,
   MinBlackAssault=MinBlackAssault,

   MissCumSumWhiteWeapons=MissCumSumWhiteWeapons,
   NonMissWhiteWeapons=NonMissWhiteWeapons,
   NmissWhiteWeapons=NmissWhiteWeapons,
   WhiteWeapons=WhiteWeapons,
   MaxWhiteWeapons=MaxWhiteWeapons,
   MinWhiteWeapons=MinWhiteWeapons,

   MissCumSumBlackWeapons=MissCumSumBlackWeapons,
   NonMissBlackWeapons=NonMissBlackWeapons,
   NmissBlackWeapons=NmissBlackWeapons,
   BlackWeapons=BlackWeapons,
   MaxBlackWeapons=MaxBlackWeapons,
   MinBlackWeapons=MinBlackWeapons,

   BlackRatio=BlackRatio,
   Pop=Pop,

   Wealth=Wealth,
   Gini=Gini,

   Ones=Ones,

   DMAIndex=DMAIndex,
   GoogleRacism=GoogleRacism,
   DMA=DMA
  )


##############################################################################################################STAN MODEL Code
model_code<-"
########################################################################################################## Data Block
data {
int<lower=0> N;

vector[N] Ym;
vector<lower=0>[N] Ysd;

vector<lower=0>[N] BlackRatio;
vector<lower=0>[N] Pop;

vector<lower=0>[N] Wealth;
vector<lower=0>[N] Gini;

int DMA[N];
vector[201] GoogleRacism;

int MissCumSumWhiteAssault[N];
int NonMissWhiteAssault[N];
int NmissWhiteAssault;
vector[N] WhiteAssault;
real MaxWhiteAssault;
real MinWhiteAssault;

int MissCumSumBlackAssault[N];
int NonMissBlackAssault[N];
int NmissBlackAssault;
vector[N] BlackAssault;
real MaxBlackAssault;
real MinBlackAssault;

int MissCumSumWhiteWeapons[N];
int NonMissWhiteWeapons[N];
int NmissWhiteWeapons;
vector[N] WhiteWeapons;
real MaxWhiteWeapons;
real MinWhiteWeapons;

int MissCumSumBlackWeapons[N];
int NonMissBlackWeapons[N];
int NmissBlackWeapons;
vector[N] BlackWeapons;
real MaxBlackWeapons;
real MinBlackWeapons;

vector[N] Ones;
}

parameters {
 vector[10] Theta;
 vector[N] log_Y;
 real<lower=0> Sigma;
 real<lower=25,upper=155> iHate;
 real<lower=MinWhiteAssault,upper=MaxWhiteAssault> iWhiteAssault[NmissWhiteAssault];
 real<lower=MinBlackAssault,upper=MaxBlackAssault> iBlackAssault[NmissBlackAssault];
 real<lower=MinWhiteWeapons,upper=MaxWhiteWeapons> iWhiteWeapons[NmissWhiteWeapons];
 real<lower=MinBlackWeapons,upper=MaxBlackWeapons> iBlackWeapons[NmissBlackWeapons];
 }

transformed parameters{
 vector<lower=0>[N] DataWhiteAssault;
 vector[N] DataBlackAssault;
 vector<lower=0>[N] DataWhiteWeapons;
 vector[N] DataBlackWeapons;
 vector<lower=0>[N] DataHate;

for(t in 1:N){
 DataHate[t] <-   if_else(DMA[t]==156,iHate,GoogleRacism[DMA[t]] );

 DataWhiteAssault[t] <- if_else(NonMissWhiteAssault[t], WhiteAssault[t], iWhiteAssault[MissCumSumWhiteAssault[t]]);
 DataBlackAssault[t] <- if_else(NonMissBlackAssault[t], BlackAssault[t], iBlackAssault[MissCumSumBlackAssault[t]]);

 DataWhiteWeapons[t] <- if_else(NonMissWhiteWeapons[t], WhiteWeapons[t], iWhiteWeapons[MissCumSumWhiteWeapons[t]]);
 DataBlackWeapons[t] <- if_else(NonMissBlackWeapons[t], BlackWeapons[t], iBlackWeapons[MissCumSumBlackWeapons[t]]);
             }
}

model {
vector[N] Mu;

log_Y ~ normal(Ym,Ysd);
Theta ~ cauchy(0,5);
Sigma ~ exponential(1);

#Mu <- ( Theta[1]*(Ones)  );

#Mu <- ( Theta[1] + Theta[2]*log(Pop)  );

#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio)  );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(Wealth) );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(Gini)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(DataHate)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(DataWhiteAssault) + Theta[4]*log(DataBlackAssault)  );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(DataWhiteWeapons) + Theta[4]*log(DataBlackWeapons)   );

#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(DataHate)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(DataWhiteAssault) + Theta[5]*log(DataBlackAssault)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(DataWhiteWeapons) + Theta[5]*log(DataBlackWeapons)   );

#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) + Theta[5]*log(Gini)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) + Theta[5]*log(DataHate)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) + Theta[5]*log(DataWhiteAssault) + Theta[6]*log(DataBlackAssault)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) + Theta[5]*log(DataWhiteWeapons) + Theta[6]*log(DataBlackWeapons)   );

#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini) + Theta[5]*log(DataHate)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini) + Theta[5]*log(DataWhiteAssault) + Theta[6]*log(DataBlackAssault)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini) + Theta[5]*log(DataWhiteWeapons) + Theta[6]*log(DataBlackWeapons)   );

#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) + Theta[5]*log(DataHate)  + Theta[6]*log(DataWhiteAssault) + Theta[7]*log(DataBlackAssault)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Wealth) + Theta[5]*log(DataHate)  + Theta[6]*log(DataWhiteWeapons) + Theta[7]*log(DataBlackWeapons)   );

#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini) + Theta[5]*log(DataHate)  + Theta[6]*log(DataWhiteAssault) + Theta[7]*log(DataBlackAssault)   );
#Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini) + Theta[5]*log(DataHate)  + Theta[6]*log(DataWhiteWeapons) + Theta[7]*log(DataBlackWeapons)   );

Mu <- ( Theta[1] + Theta[2]*log(Pop) + Theta[3]*log(BlackRatio) + Theta[4]*log(Gini) + Theta[5]*log(Wealth) + Theta[6]*log(DataHate)  + Theta[7]*log(DataWhiteWeapons) + Theta[8]*log(DataBlackWeapons)  + Theta[9]*log(DataWhiteAssault) + Theta[10]*log(DataBlackAssault)  );

log_Y ~  normal(Mu,Sigma);
}
"


################################################################################ Fit the Model IN STAN!
iter<-4000
warmup<-2000
fitKilling <- stan(model_code=model_code, data = model_dat, thin=1, iter = iter, warmup=warmup,chains = 1,refresh=10,pars=c("Theta","Sigma"))

print(fitKilling,digits_summary=4,pars=c("Theta","Sigma"))










# Below is for checking a few plots
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('antiquewhite','darkred'))

#This adds a column of color values
# based on the y values

Pop<-d$TOT_POP
BlackRatio<-(d$BAC_TOT+1)/Pop
COL <- rbPal(1000)[as.numeric(cut(d$TOT_POP,breaks = 1000))]

plot(d[,58]~d[,59], ylab="Black Armed-Unarmed Ratio",xlab="White Armed-Unarmed Ratio",pch=20,col = COL,cex=1.5)
abline(0,1)


 #Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('goldenrod','black'))

#This adds a column of color values
# based on the y values

Pop<-d$TOT_POP
HRatio<-(d$H_TOT+1)/Pop
COL <- rbPal(100)[as.numeric(cut(HRatio,breaks = 100))]

plot(d[,60]~d[,59], ylab="Black Armed-Unarmed Ratio",xlab="White Armed-Unarmed Ratio",pch=20,col = COL,cex=1.5)
abline(0,1)



####################### Below is for the tables in the text
      d1<-d[1:25,2:12]
      d2<-d[28:52,2:12]

      d1<-m
      d2<-sd
       XX <- matrix(NA, nrow=15,ncol=11)

       for(i in 1:15){
       for(j in 1:11){
       XX[i,j] <- pmin(pnorm(0,mean=as.numeric(as.character(d1[i,j])),sd=as.numeric(as.character(d2[i,j])), lower.tail = TRUE, log.p = FALSE),
               1-pnorm(0,mean=as.numeric(as.character(d1[i,j])),sd=as.numeric(as.character(d2[i,j])), lower.tail = TRUE, log.p = FALSE))
             }}














