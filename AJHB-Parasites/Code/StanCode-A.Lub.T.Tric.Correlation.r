################################## Load Data ####################################################################################
setwd("insert path")
 d<-read.csv("Data-A.Lub.T.Tric.csv")

 ParasitesA <- d$sumPositiveA[1:51]
 TestedA <- d$sumTestedA[1:51]
 
 ParasitesT <- d$sumPositiveT[1:51]
 TestedT <- d$sumTestedT[1:51]
 
 ParasitesT[is.na(ParasitesT)]<-99999999
 TestedT[is.na(TestedT)]<-99999999


########################################################################################################## STAN MODEL Data 
  L<-51
  S<-8
 model_dat  <-list(S=S,L=L,ParasitesA=ParasitesA, TestedA=TestedA,ParasitesT=ParasitesT, TestedT=TestedT)
    
  
##############################################################################################################STAN MODEL Code  
model_code<-"
########################################################################################################## Data Block 
data {

int<lower=1> S;
int<lower=1> L;
 
int<lower=0> ParasitesA[L];
int<lower=0> TestedA[L];

int<lower=0> ParasitesT[L];
int<lower=0> TestedT[L];

}

parameters {
#######
# Theta hieracrchy
#######

real muGrandA;
real<lower=0> sdGrandA;
real muGrandAdisp;
real<lower=0> sdGrandAdisp;

real muGrandT;
real<lower=0> sdGrandT;
real muGrandTdisp;
real<lower=0> sdGrandTdisp;

vector[S] muNationThetaA;
vector<lower=0>[S] sdNationThetaA;

vector[S] muNationThetaT;
vector<lower=0>[S] sdNationThetaT;

real ThetaA[L];
real ThetaT[L];
}


model {
muGrandA~normal(0,5);
sdGrandA~cauchy(0,2.5);
muGrandAdisp~normal(0,5);
sdGrandAdisp~cauchy(0,2.5);

muGrandT~normal(0,5);
sdGrandT~cauchy(0,2.5);
muGrandTdisp~normal(0,5);
sdGrandTdisp~cauchy(0,2.5);

for(i in 1:8){
muNationThetaA[i] ~ normal(muGrandA,sdGrandA);  
sdNationThetaA[i] ~ normal(muGrandAdisp,sdGrandAdisp); }  

for(i in 1:8){
muNationThetaT[i] ~ normal(muGrandT,sdGrandT);  
sdNationThetaT[i] ~ normal(muGrandTdisp,sdGrandTdisp); }  

# Sample Nation Level Prior Theta for Ascarsis
for(i in 1:4){
ThetaA[i]~normal(muNationThetaA[1],exp(sdNationThetaA[1])); }
for(i in 5:9){
ThetaA[i]~normal(muNationThetaA[2],exp(sdNationThetaA[2])); }
for(i in 10:13){
ThetaA[i]~normal(muNationThetaA[3],exp(sdNationThetaA[3])); }
for(i in 14:24){
ThetaA[i]~normal(muNationThetaA[4],exp(sdNationThetaA[4])); }
for(i in 25:31){
ThetaA[i]~normal(muNationThetaA[5],exp(sdNationThetaA[5])); }
#for(i in 32:32){
ThetaA[32]~normal(muNationThetaA[6],exp(sdNationThetaA[6])); 
for(i in 33:38){
ThetaA[i]~normal(muNationThetaA[7],exp(sdNationThetaA[7])); }
for(i in 39:51){
ThetaA[i]~normal(muNationThetaA[8],exp(sdNationThetaA[8])); }

# Sample Nation Level Prior Theta for Trich
for(i in 1:4){
ThetaT[i]~normal(muNationThetaT[1],exp(sdNationThetaT[1])); }
for(i in 5:9){
ThetaT[i]~normal(muNationThetaT[2],exp(sdNationThetaT[2])); }
for(i in 10:13){
ThetaT[i]~normal(muNationThetaT[3],exp(sdNationThetaT[3])); }
for(i in 14:24){
ThetaT[i]~normal(muNationThetaT[4],exp(sdNationThetaT[4])); }
for(i in 25:31){
ThetaT[i]~normal(muNationThetaT[5],exp(sdNationThetaT[5])); }
#for(i in 32:32){
ThetaT[32]~normal(muNationThetaT[6],exp(sdNationThetaT[6])); 
for(i in 33:38){
ThetaT[i]~normal(muNationThetaT[7],exp(sdNationThetaT[7])); }
for(i in 39:51){
ThetaT[i]~normal(muNationThetaT[8],exp(sdNationThetaT[8])); }

# Model District Level for Acars 
for(l in 1:L){
  ParasitesA[l]~binomial(TestedA[l],inv_logit(ThetaA[l]));
}
  
# Model District Level for Trich - note missing data
for(l in 1:22){
  ParasitesT[l]~binomial(TestedT[l],inv_logit(ThetaT[l]));
}

ParasitesT[24]~binomial(TestedT[24],inv_logit(ThetaT[24]));

for(l in 26:30){
  ParasitesT[l]~binomial(TestedT[l],inv_logit(ThetaT[l]));
}

ParasitesT[33]~binomial(TestedT[33],inv_logit(ThetaT[33]));
ParasitesT[34]~binomial(TestedT[34],inv_logit(ThetaT[34]));

for(l in 36:L){
  ParasitesT[l]~binomial(TestedT[l],inv_logit(ThetaT[l]));
}

}

generated quantities{
real predThetaA[L];   
real predThetaT[L];  
real rho; 

for(i in 1:L){
predThetaA[i]<-inv_logit(ThetaA[i]);
predThetaT[i]<-inv_logit(ThetaT[i]);
}


rho <- (dot_product(predThetaA, predThetaT) * pow(L, -1) - (mean(predThetaA) * mean(predThetaT))) * pow((sd(predThetaA) * sd(predThetaT)), -1);

############################################# End Model###########################################
}"


################################################################################ Fit the Model IN STAN!
fitParasitesCor <- stan(model_code=model_code, data = model_dat,inits=0, iter = 10000, warmup=2000,chains = 1)

  print(fitParasites3,digits_summary=4)
  
  A<-get_posterior_mean(fitParasites3,pars="predThetaA")
  T<-get_posterior_mean(fitParasites3,pars="predThetaT")
  
  plot(A~T)
  cor(A,T)


    A<-extract(fitParasites3,pars="predThetaA")$predThetaA
    T<-extract(fitParasites3,pars="predThetaT")$predThetaT
    
   rho<-c() 
    for(i in 1:dim(A)[1]){
    rho[i] <- cor(A[i,],T[i,])
    }

        A<-extract(fitParasites3,pars="predThetaA")$predThetaA
    T<-extract(fitParasites3,pars="predThetaT")$predThetaT
    
   rho<-c() 
    for(i in 1:dim(A)[1]){
    rho[i] <- cor(A[i,],T[i,])
    }



