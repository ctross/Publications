################################# Load Data ####################################################################################
library(rstan)
setwd("insert path")
 d<-read.csv("Data-A.Lub.Full.csv")
  
 Y<-d[1:51,c(78,77)]
 N<-Y[,1]+Y[,2]

################################# Pre-Processing in R ###########################################################################
L<-51
D<-2
K<-4
Y<-Y
S<-8

########################################################################################################## STAN MODEL Data 
model_dat  <-list(S=S,L=L,D=D,K=K,Parasites=d$sumPositive, Tested=d$sumTested,Outcome=Y[,1],N=N)
    
  
##############################################################################################################STAN MODEL Code  
model_code<-"
########################################################################################################## Data Block 
data {
int<lower=2> K;
int<lower=0> L;
int<lower=1> D;
int<lower=1> S;

int<lower=0> Outcome[L]; 
int<lower=0> N[L];  
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
vector[D] Beta;
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
for (d in 1:D){
Beta[d] ~ normal(0,10);} 

##################### Main Regression
# The code used in the model implies the commented out model, but is more effcient computationally
#for (n in 1:N){
#Outcome[n] ~ bernoulli((Beta[1] + Beta[2]*inv_logit(Theta[Cluster[n]])  ));
#}
for (n in 1:L) {
Outcome[n]~ binomial(N[n], inv_logit(Beta[1] + Beta[2]*inv_logit(Theta[n])));
}
}"


################################################################################ Fit the Model IN STAN!
fitParasitesSewer<- stan(model_code=model_code, data = model_dat,inits=0, thin=10,iter = 20000, warmup=2000,chains = 1)
 
  print(fitParasitesSewer,digits_summary=3)
  traceplot(fitParasitesSewer,ask=T,pars=c("Beta"))
  
  Beta<-extract(fitParasitesSewer, pars="Beta")$Beta
  Theta<-seq(0,1,length.out=21)
 
 library(rethinking)
  Pred<-array(NA, c(1800,21))
 for(i in 1:1800){ 
  for (n in 1:21) {
 Pred[i,n]<-logistic(Beta[i,1] + Beta[i,2]*logistic(Theta[n]))
 }} 
  
  mPred<-c()
  
  for(i in 1:21){
  mPred[i]<-mean(Pred[,i])
  }
  
############ Plot Results
library(Cairo)
CairoPDF(width=8, height=8, "Sewer.pdf")

  plot(mPred~Theta,ylim=c(0,1),type="n",col=rgb(1, 0, 0,alpha=.9),ylab="Probability of Having Sewage System", xlab="Parasite Prevalence" )
 
  for(z in 1:900){
  zz<-z*2
  lines(Pred[zz,]~Theta, col=rgb(1, 0, 0,alpha=.05))
  }
        
       dev.off()
       


