################################# Load Data ####################################################################################
library(rstan)
setwd("insert path")
 d<-read.csv("Data-A.Lub.Full.csv")
  
 Y<-cbind(d$Democracy, d$DoesntMatter , d$Authoritarian  )

################################# Pre-Processing in R ###########################################################################
L<-51
D<-2
K<-3
Y<-Y
S<-8

########################################################################################################## STAN MODEL Data 
model_dat  <-list(S=S,L=L,D=D,K=K,Parasites=d$sumPositive, Tested=d$sumTested,Outcome=Y)
    
  
##############################################################################################################STAN MODEL Code  
model_code<-"
########################################################################################################## Data Block 
data {
int<lower=2> K;
int<lower=0> L;
int<lower=1> D;
int<lower=1> S;

int<lower=0> Outcome[L,K];  
int<lower=0> Parasites[L];
int<lower=0> Tested[L];

}

parameters {
#######
# Theta hieracrchy
#######

real muGrand;
real<lower=0> sdGrand;
real muGranddisp;
real<lower=0> sdGranddisp;

vector[S] muNationTheta;
vector<lower=0>[S] sdNationTheta;

real Theta[L];

#### Regression
matrix[(K-1),D] Beta0;
}

transformed parameters {
vector[D] X[L]; 
matrix[K,D] Beta;
 
 for(l in 1:L){
 X[l,1]<-1;
 X[l,2] <- inv_logit(Theta[l]);
        }   
        
# Set one coef to identity
for (d in 1:D){
Beta[1,d] <-Beta0[1,d];
Beta[2,d]<-1;
Beta[3,d] <-Beta0[2,d];
}                                                                    
}

model {
muGrand~normal(0,5);
sdGrand~cauchy(0,2.5);
muGranddisp~normal(0,5);
sdGranddisp~cauchy(0,2.5);

for(i in 1:8){
muNationTheta[i] ~ normal(muGrand,sdGrand);  
sdNationTheta[i] ~ normal(muGranddisp,sdGranddisp); }  

# Sample Nation Level Prior Theta for Ascarsis
for(i in 1:4){
Theta[i]~normal(muNationTheta[1],exp(sdNationTheta[1])); }
for(i in 5:9){
Theta[i]~normal(muNationTheta[2],exp(sdNationTheta[2])); }
for(i in 10:13){
Theta[i]~normal(muNationTheta[3],exp(sdNationTheta[3])); }
for(i in 14:24){
Theta[i]~normal(muNationTheta[4],exp(sdNationTheta[4])); }
for(i in 25:31){
Theta[i]~normal(muNationTheta[5],exp(sdNationTheta[5])); }
#for(i in 32:32){
Theta[32]~normal(muNationTheta[6],exp(sdNationTheta[6])); 
for(i in 33:38){
Theta[i]~normal(muNationTheta[7],exp(sdNationTheta[7])); }
for(i in 39:51){
Theta[i]~normal(muNationTheta[8],exp(sdNationTheta[8])); }

# Model District Level for Acars 
for(l in 1:L){
  Parasites[l]~binomial(Tested[l],inv_logit(Theta[l]));
}

###################### Priors for Regression 
for (k in 1:(K-1)){
for (d in 1:D){
Beta0[k,d] ~ normal(0,10);}}  

##################### Main Regression
for (l in 1:L){
Outcome[l] ~ multinomial(softmax((Beta * X[l]))); 
}


}

generated quantities{
real<lower=0, upper=1> Prevalence[L];
vector[K] Pred[L];

for (l in 1:L){
Prevalence[l]<-inv_logit(Theta[l]);
}

for (l in 1:L){
Pred[l]<-softmax((Beta * X[l]));
}      

############################################# End Model###########################################
}"


################################################################################ Fit the Model IN STAN!
fitParasitesDem <- stan(model_code=model_code, data = model_dat,inits=0, iter = 10000, warmup=2000,chains = 1)

  print(fitParasitesDem,digits_summary=4)
  
  Prevalence<-extract(fitParasitesDem ,pars="Prevalence")$Prevalence
  Pred<-extract(fitParasitesDem ,pars="Pred")$Pred
  
  mPrevalence <-c()
  mPred<-matrix(NA, nrow=L,ncol=K)
  
  for(i in 1:L){
  mPrevalence[i]<-mean(Prevalence[,i])
  for(j in 1:K){
  mPred[i,j]<-mean(Pred[,i,j])
  }}
  
  
  
line(smooth.spline(mPred[,1]~mPrevalence),type="l",col="blue",ylim=c(0,1))
lines((mPred[,3]~mPrevalence),type="l",col="red")
lines((mPred[,2]~mPrevalence),type="l",col="orange")



############ Plot Results
library(Cairo)
CairoPDF(width=8, height=8, "DemocracyModel.pdf")
plot((mPred[,1]~mPrevalence),type="n",ylim=c(0,1), xlab="Parasite Prevalence", ylab="Probability of Selecting Category")

for( i in 1:100){
z<-i*80
lines(smooth.spline(Pred[z,,1]~Prevalence[z,]),type="l",col=rgb(0, 0, 1,alpha=.05),ylim=c(0,1))
}

for( i in 1:100){
z<-i*80
lines(smooth.spline(Pred[z,,2]~Prevalence[z,]),type="l",col=rgb(1, 0, 1,alpha=.05),ylim=c(0,1))
}

for( i in 1:100){
z<-i*80
lines(smooth.spline(Pred[z,,3]~Prevalence[z,]),type="l",col=rgb(1, 0, 0,alpha=.05),ylim=c(0,1))
}

legend("topleft", c("Democracy", "Indifferent", "Authoritarian"), col = c(col=rgb(0, 0, 1,alpha=.75), col=rgb(1, 0, 1,alpha=.75), col=rgb(1, 0, 0,alpha=.75)),
       text.col = "black", lwd = c(2, 2, 2),  merge = TRUE, bg = "gray99")
       
       dev.off()
       


