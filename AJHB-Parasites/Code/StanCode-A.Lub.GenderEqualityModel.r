################################# Load Data ####################################################################################
library(rstan)

setwd("insert path")
 d<-read.csv("Data-A.Lub.Full.csv")

 
  
 Y<-d[1:51,141:138]

################################# Pre-Processing in R ###########################################################################
L<-51
D<-2
K<-4
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
vector[D] Beta;
ordered[K-1] C;
}


model {
vector[K] scrap;

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
#Outcome[n] ~ ordered_logistic((Beta[1] + Beta[2]*inv_logit(Theta[Cluster[n]])  ), C);
#}
for (n in 1:L) {
real eta;
eta <- Beta[1] + Beta[2]*inv_logit(Theta[n]);
scrap[1] <- 1 - inv_logit(eta - C[1]);
for (k in 2:(K-1))
scrap[k] <- inv_logit(eta - C[k-1]) - inv_logit(eta - C[k]);
scrap[K] <- inv_logit(eta - C[K-1]);
Outcome[n] ~ multinomial(scrap);

}
}"


################################################################################ Fit the Model IN STAN!
fitParasitesGender<- stan(model_code=model_code, data = model_dat,inits=0, iter = 4000, warmup=2000,chains = 1)
 
  print(fitParasitesGender,digits_summary=3)
  traceplot(fitParasitesGender,ask=T,pars=c("C","Beta"))
  
  Beta<-extract(fitParasitesGender, pars="Beta")$Beta
  C<-extract(fitParasitesGender, pars="C")$C
  Theta<-seq(0,0.65,length.out=21)
 
 library(rethinking)
  Pred<-array(NA, c(2000,(K),21))
 for(i in 1:2000){ 
  for (n in 1:21) {
eta <- Beta[i,1] + Beta[i,2]*logistic(Theta[n]);
scrap<-c()
scrap[1] <- 1 - logistic(eta - C[i,1]);
for (k in 2:(K-1)){
scrap[k] <- logistic(eta - C[i,k-1]) - logistic(eta - C[i,k]);}
scrap[K] <- logistic(eta - C[i,K-1]);

 Pred[i,,n]<-cumsum(scrap)
 }} 
  
  mPred<-matrix(NA,nrow=21,ncol=K)
  
  for(i in 1:21){
  for(j in 1:K){
  mPred[i,j]<-mean(Pred[,j,i])
  }}
  
############ Plot Results
library(Cairo)
CairoPDF(width=8, height=8, "GenderEqualityModel.pdf")

  plot(mPred[,1]~Theta,ylim=c(0,1),type="n",col=rgb(1, 0, 0,alpha=.9),ylab="Cumulative Probability", xlab="Parasite Prevalence" )
 
  for(i in 0:(K-2)){
  for(z in 1:200){
  zz<-z*10
  lines(Pred[zz,i+1,]~Theta, col=rgb((1-(i*.06)), 0, ((i*.06)),alpha=.06))}
  }
  
 # for(i in 0:(K-2)){
 # lines(mPred[,i+1]~Theta, col=rgb((1-(i*.06)), 0, ((i*.06)),alpha=.5))
 # }
  
  abline(h=.2)
  abline(h=.4)
  abline(h=.6)
  abline(h=.8)
  abline(h=1)
  abline(h=0)
       
       dev.off()
       


