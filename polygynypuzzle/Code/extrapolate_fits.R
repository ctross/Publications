################################################################################
# Age Adjustment
################################################################################
rm(list=ls())
load("ModelResults.RData")

all_stuff <- ls()
rm(list=setdiff(all_stuff, c("fitPoly", "P", "N", 
  "Rival_1", "Rival_2", "MaxExposure", "Exposure", "Wives"
)))

source('./Code/project_support.R')
Samp <-2000

############################################################## Make Full Data
ShadowPrice <- extract(fitPoly,pars="ShadowPrice")$ShadowPrice

Beta_Rival_1 <- extract(fitPoly,pars="Beta_Rival_1")$Beta_Rival_1
B_Rival_1 <- extract(fitPoly,pars="B_Rival_1")$B_Rival_1

Beta_Rival_2 <- extract(fitPoly,pars="Beta_Rival_2")$Beta_Rival_2
B_Rival_2 <- extract(fitPoly,pars="B_Rival_2")$B_Rival_2

Beta_Wives <- extract(fitPoly,pars="Beta_Wives")$Beta_Wives
B_Wives <- extract(fitPoly,pars="B_Wives")$B_Wives

gammaZ_rng <-function(beta,x) {
  SAMP <- rep(NA,length(beta));
  for (i in 1:length(beta)){
    if( beta[i]>0){
      SAMP[i] <- rgamma(1,beta[i],x);
    } else{
      SAMP[i] <- 0;
    }
  }
  SAMP;
}

poisZ_rng <-function(beta) {
  SAMP <- rep(NA,length(beta));
    for (i in 1:length(beta)){
      if( beta[i]>0){
        SAMP[i] <- rpois(1,beta[i]);
      }else{
        SAMP[i] <- 0;
      }
    }
  SAMP;
}

Rival<-vector('list',P)

for(i in 1:P){
  Gmat<-matrix(NA,nrow=Samp,ncol=N[i])
  for(j in 1:Samp){
    Gmat[j,]<- head(Rival_1[,i],N[i]) + exp(ShadowPrice[j,i])*head(Rival_2[,i],N[i]);
  }
  Rival[[i]]  <- Gmat
}
print("calculated current rival wealth composite based on model shadow prices")

AdjustedRival1<-vector('list',P)

for(i in 1:P){
  Gmat<-matrix(NA,nrow=Samp,ncol=N[i])
  for(j in 1:Samp){
    Mu <- exp(Beta_Rival_1[j,i,1])*((MaxExposure^Beta_Rival_1[j,i,2]) - (head(Exposure[,i],N[i])^Beta_Rival_1[j,i,2]))
    Gmat[j,] <- head(Rival_1[,i],N[i]) + gammaZ_rng(Mu*B_Rival_1[j,i],B_Rival_1[j,i]);
  }
  AdjustedRival1[[i]] <- Gmat
}


AdjustedRival2<-vector('list',P)

for(i in 1:P){
  Gmat<-matrix(NA,nrow=Samp,ncol=N[i])
  for(j in 1:Samp){
    Mu <- exp(Beta_Rival_2[j,i,1])*((MaxExposure^Beta_Rival_2[j,i,2]) - (head(Exposure[,i],N[i])^Beta_Rival_2[j,i,2]))
    Gmat[j,] <- head(Rival_2[,i],N[i]) + gammaZ_rng(Mu*B_Rival_2[j,i],B_Rival_2[j,i]);
  }
  AdjustedRival2[[i]] <- Gmat
}


AdjustedRival<-vector('list',P)

for(i in 1:P){
  Gmat<-matrix(NA,nrow=Samp,ncol=N[i])
  for(j in 1:Samp){
    Gmat[j,]<- AdjustedRival1[[i]][j,] + exp(ShadowPrice[j,i])*AdjustedRival2[[i]][j,];
  }
  AdjustedRival[[i]]  <- Gmat
}
print("extrapolated rival wealth at end of reproductive career")

# RealWives<-vector('list',P)
# for(i in 1:P){
#   RealWives[[i]]  <- head(Wives[,i],N[i])
# }
# print("cumulative number of marriages at observation")

AdjustedWives<-vector('list',P)
UnAdjustedWives<-vector('list',P) 
AgeVec<-vector('list',P)


for(i in 1:P){
  Gmat<-matrix(NA,nrow=Samp,ncol=N[i])
  Gmat2<-matrix(NA,nrow=Samp,ncol=N[i])
  Gmat3<-matrix(NA,nrow=Samp,ncol=N[i])
  for(j in 1:Samp){
    MuWives <- exp(Beta_Wives[j,i,1])*(((AdjustedRival[[i]][j,]^Beta_Wives[j,i,2]) * (MaxExposure^Beta_Wives[j,i,3])) -
  ((Rival[[i]][j,]^Beta_Wives[j,i,2]) * (head(Exposure[,i],N[i])^Beta_Wives[j,i,3])))
    Gmat[j,]<-  head(Wives[,i],N[i]) + poisZ_rng(gammaZ_rng(MuWives*B_Wives[j,i],B_Wives[j,i]));
    Gmat2[j,]<-  head(Wives[,i],N[i]) 
    Gmat3[j,]<-  head(Exposure[,i],N[i]) 
  }
  AdjustedWives[[i]]  <- Gmat
  UnAdjustedWives[[i]]  <- Gmat2
  AgeVec[[i]]  <- Gmat3
}

print("extrapolated number of wives at end of reproductive career")


save( "AdjustedWives", "AdjustedRival", "UnAdjustedWives", "AgeVec", file="Extrapolations.RData")
