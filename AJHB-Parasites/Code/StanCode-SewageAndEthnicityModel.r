################################# Load Data ####################################################################################
library(rstan)
library(rethinking)
setwd("insert path")
 d<-read.csv("Data-A.Lub.Full.csv")
  
 Y<-cbind(d$White, d$Mestizo , d$Mulato, d$Indigenous, d$Black  )
 
 Sewer<-d[1:51,c(78,77)]
 N<-Sewer[,1]+Sewer[,2]


################################# Pre-Processing in R ###########################################################################
L<-51
D<-2
K<-5
Y<-Y
S<-8


########################################################################################################## STAN MODEL Data 
model_dat  <-list(S=S,L=L,D=D,K=K,Sewer=Sewer[,1], Tested=N,Outcome=Y)
    
  
##############################################################################################################STAN MODEL Code  
model_code<-"
########################################################################################################## Data Block 
data {
int<lower=2> K;
int<lower=0> L;
int<lower=1> D;
int<lower=1> S;

int<lower=0> Outcome[L,K];  
int<lower=0> Sewer[L];
int<lower=0> Tested[L];

}

parameters {
#######
# Theta 
#######

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
Beta[1,d] <- 1;
Beta[2,d] <- Beta0[1,d];
Beta[3,d] <- Beta0[2,d];
Beta[4,d] <- Beta0[3,d];
Beta[5,d] <- Beta0[4,d];
}                                                                    
}

model {

Theta~normal(0,10);

# Model District Level for Acars 
for(l in 1:L){
  Sewer[l]~binomial(Tested[l],inv_logit(Theta[l]));
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
fitParasitesSewerEth <- stan(model_code=model_code, data = model_dat,inits=0, iter = 10000, warmup=2000,chains = 1)

  print(fitParasitesSewerEth,digits_summary=3)
  traceplot(fitParasitesSewerEth,ask=T,pars = paste0("Beta[", c(1, 1, 2, 2,3,3,4,4), ",", c(1, 2, 1, 2,1,2,1,2), "]"))
  
  Prevalence<-extract(fitParasitesSewerEth ,pars="Prevalence")$Prevalence
  Pred<-extract(fitParasitesSewerEth,pars="Pred")$Pred
  
  mPrevalence <-c()
  mPred<-matrix(NA, nrow=L,ncol=K)
  
  for(i in 1:L){
  mPrevalence[i]<-mean(Prevalence[,i])
  for(j in 1:K){
  mPred[i,j]<-mean(Pred[,i,j])
  }}
  
  

############ Plot Results
library(Cairo)
CairoPDF(width=8, height=8, "SewageAndEthnicityModel.pdf")
plot((mPred[,1]~mPrevalence),type="n",ylim=c(0,1),xlim=c(0,1), xlab="Percent of Population with Sewage and Sanitation Services", ylab="Ethnicity as Proportion of Population")

for( i in 1:200){
z<-i*40
lines(smooth.spline(Pred[z,,1]~Prevalence[z,]),type="l",col=rgb(as.vector(col2rgb(brewer.pal(12,"Paired")[2])/255)[1],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[2])/255)[2],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[2])/255)[3],alpha=.05),ylim=c(0,1))}
for( i in 1:200){
z<-i*40
lines(smooth.spline(Pred[z,,2]~Prevalence[z,]),type="l",col=rgb(as.vector(col2rgb(brewer.pal(12,"Paired")[4])/255)[1],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[4])/255)[2],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[4])/255)[3],alpha=.05),ylim=c(0,1))}
for( i in 1:200){
z<-i*40
lines(smooth.spline(Pred[z,,3]~Prevalence[z,]),type="l",col=rgb(as.vector(col2rgb(brewer.pal(12,"Paired")[6])/255)[1],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[6])/255)[2],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[6])/255)[3],alpha=.05),ylim=c(0,1))}
for( i in 1:200){
z<-i*40
lines(smooth.spline(Pred[z,,4]~Prevalence[z,]),type="l",col=rgb(as.vector(col2rgb(brewer.pal(12,"Paired")[8])/255)[1],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[8])/255)[2],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[8])/255)[3],alpha=.05),ylim=c(0,1))}
for( i in 1:200){
z<-i*40
lines(smooth.spline(Pred[z,,5]~Prevalence[z,]),type="l",col=rgb(as.vector(col2rgb(brewer.pal(12,"Paired")[10])/255)[1],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[10])/255)[2],
                                                                as.vector(col2rgb(brewer.pal(12,"Paired")[10])/255)[3],alpha=.05),ylim=c(0,1))}


legend("top", c("White", "Mestizo" , "Mulato", "Indigenous", "Black"), col = c(brewer.pal(12,"Paired")[c(2,4,6,8,10)]),
       text.col = "black", lwd = c(4, 4, 4,4,4),  merge = TRUE, bg = "gray99")
       
       dev.off()
       


