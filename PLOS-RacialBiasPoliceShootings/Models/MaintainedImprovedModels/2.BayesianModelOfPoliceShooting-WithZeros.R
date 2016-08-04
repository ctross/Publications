
################################## Load Data ###################################
 d<-read.csv("MapFileData-WithCountyResultsAndCovariates.csv") # 

 library(rstanmulticore)

##############################################################################
# This code is used to model the relative risk across race-ethnicity
# and armed-unarmed status

UnarmedBlack<-d$BlackUnarmed	
ArmedBlack<-d$BlackArmed
UnarmedHispanic<-d$HispanicUnarmed	
ArmedHispanic<-d$HispanicArmed		
UnarmedWhite<-d$WhiteUnarmed		
ArmedWhite<-d$WhiteArmed

UnarmedBlack[which(is.na(UnarmedBlack))]<-0
ArmedBlack[which(is.na(ArmedBlack))]<-0
UnarmedHispanic[which(is.na(UnarmedHispanic))]<-0
ArmedHispanic[which(is.na(ArmedHispanic))]<-0
UnarmedWhite[which(is.na(UnarmedWhite))]<-0
ArmedWhite[which(is.na(ArmedWhite))]<-0

Nblack<-d$BAC_TOT
Nwhite<-d$WA_TOT
Nhispanic<-d$H_TOT

Nblack[which(is.na(Nblack))]<-0
Nwhite[which(is.na(Nwhite))]<-0
Nhispanic[which(is.na(Nhispanic))]<-0

G<-data.frame(UnarmedBlack, ArmedBlack, UnarmedHispanic, ArmedHispanic, UnarmedWhite, ArmedWhite, Nblack, Nwhite, Nhispanic)
G<-G[which(Nblack*Nwhite*Nhispanic>0),]

UnarmedBlack<-G$UnarmedBlack
ArmedBlack<-G$ArmedBlack
UnarmedHispanic<-G$UnarmedHispanic
ArmedHispanic<-G$ArmedHispanic	
UnarmedWhite<-G$UnarmedWhite	
ArmedWhite<-G$ArmedWhite

Nblack<-G$Nblack
Nwhite<-G$Nwhite
Nhispanic<-G$Nhispanic

N<-length(Nblack)

model_dat  <-list(N=N,
Nblack=Nblack,
Nwhite=Nwhite,
Nhispanic=Nhispanic,
UnarmedBlack=UnarmedBlack,	
ArmedBlack=ArmedBlack,		
UnarmedHispanic=UnarmedHispanic,		
ArmedHispanic=ArmedHispanic,		
UnarmedWhite=UnarmedWhite,		
ArmedWhite=ArmedWhite
  )   
  
##############################################################################################################STAN MODEL Code  
model_code<-"
########################################################################################################## Data Block 
data {
  int<lower=0> N;

  int<lower=1> Nwhite[N];
  int<lower=1> Nblack[N];
  int<lower=1> Nhispanic[N];
  int<lower=0> UnarmedWhite[N];
  int<lower=0> UnarmedBlack[N];
  int<lower=0> UnarmedHispanic[N];

  int<lower=0> ArmedWhite[N];
  int<lower=0> ArmedBlack[N];
  int<lower=0> ArmedHispanic[N];
}

parameters {
   vector[6] Mu;
   vector<lower=0>[6] Sigma;
   corr_matrix[6] Rho;
 
   vector[6] Theta[N];
}
   
model {
######################################################### Priors
   Mu ~ normal(-14,4);
   Sigma ~ cauchy(0,5);

   for(i in 1:N){
   Theta[i] ~ multi_normal_cholesky(Mu, (diag_matrix(Sigma) * cholesky_decompose(Rho)) );
   }


####################################################### Data Modeling
for(i in 1:N){
   ArmedBlack[i]~binomial(Nblack[i],inv_logit(Theta[i,1]));
   ArmedWhite[i]~binomial(Nwhite[i],inv_logit(Theta[i,3]));
   ArmedHispanic[i]~binomial(Nhispanic[i],inv_logit(Theta[i,5]));  

   UnarmedBlack[i]~binomial(Nblack[i],inv_logit(Theta[i,2])); 
   UnarmedWhite[i]~binomial(Nwhite[i],inv_logit(Theta[i,4]));
   UnarmedHispanic[i]~binomial(Nhispanic[i],inv_logit(Theta[i,6]));
   }

  }
  
generated quantities{
######################################################### Mean Quanitities
real Mu_Black_Armed; 
real Mu_White_Armed; 
real Mu_Hispanic_Armed;

real Mu_Black_Unarmed; 
real Mu_White_Unarmed; 
real Mu_Hispanic_Unarmed; 
   
real Mu_RR_Black_Armed_Versus_Unarmed; 
real Mu_RR_White_Armed_Versus_Unarmed; 
real Mu_RR_Hispanic_Armed_Versus_Unarmed;

real Mu_RR_Black_Armed_Versus_White_Armed; 
real Mu_RR_Hispanic_Armed_Versus_White_Armed; 
real Mu_RR_Hispanic_Armed_Versus_Black_Armed; 

real Mu_RR_Black_Unarmed_Versus_White_Unarmed; 
real Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed; 
real Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed; 

real Mu_RR_Black_Unarmed_Versus_White_Armed; 
real Mu_RR_Hispanic_Unarmed_Versus_White_Armed; 

######################################################### Quanitities By County 
vector[N]  Black_Armed; 
vector[N]  White_Armed; 
vector[N]  Hispanic_Armed;

vector[N]  Black_Unarmed; 
vector[N]  White_Unarmed; 
vector[N]  Hispanic_Unarmed; 
   
vector[N]  RR_Black_Armed_Versus_Unarmed; 
vector[N]  RR_White_Armed_Versus_Unarmed; 
vector[N]  RR_Hispanic_Armed_Versus_Unarmed;

vector[N]  RR_Black_Armed_Versus_White_Armed; 
vector[N]  RR_Hispanic_Armed_Versus_White_Armed; 
vector[N]  RR_Hispanic_Armed_Versus_Black_Armed; 

vector[N]  RR_Black_Unarmed_Versus_White_Unarmed; 
vector[N]  RR_Hispanic_Unarmed_Versus_White_Unarmed; 
vector[N]  RR_Hispanic_Unarmed_Versus_Black_Unarmed; 

vector[N]  RR_Black_Unarmed_Versus_White_Armed; 
vector[N]  RR_Hispanic_Unarmed_Versus_White_Armed; 

############################################################################################ Calc Means
  Mu_Black_Armed<-inv_logit(Mu[1]); 
  Mu_White_Armed<-inv_logit(Mu[3]); 
  Mu_Hispanic_Armed<-inv_logit(Mu[5]);

  Mu_Black_Unarmed<-inv_logit(Mu[2]); 
  Mu_White_Unarmed<-inv_logit(Mu[4]); 
  Mu_Hispanic_Unarmed<-inv_logit(Mu[6]);


  Mu_RR_Black_Armed_Versus_Unarmed            <-  inv_logit(Mu[1])/inv_logit(Mu[2]);
  Mu_RR_White_Armed_Versus_Unarmed            <-  inv_logit(Mu[3])/inv_logit(Mu[4]);
  Mu_RR_Hispanic_Armed_Versus_Unarmed         <-  inv_logit(Mu[5])/inv_logit(Mu[6]);

  Mu_RR_Black_Armed_Versus_White_Armed        <-  inv_logit(Mu[1])/inv_logit(Mu[3]);
  Mu_RR_Hispanic_Armed_Versus_White_Armed     <-  inv_logit(Mu[5])/inv_logit(Mu[3]);
  Mu_RR_Hispanic_Armed_Versus_Black_Armed     <-  inv_logit(Mu[5])/inv_logit(Mu[1]);

  Mu_RR_Black_Unarmed_Versus_White_Unarmed    <-  inv_logit(Mu[2])/inv_logit(Mu[4]);
  Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed <-  inv_logit(Mu[6])/inv_logit(Mu[4]);
  Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed <-  inv_logit(Mu[6])/inv_logit(Mu[2]);

  Mu_RR_Black_Unarmed_Versus_White_Armed      <-  inv_logit(Mu[2])/inv_logit(Mu[3]); 
  Mu_RR_Hispanic_Unarmed_Versus_White_Armed   <-  inv_logit(Mu[6])/inv_logit(Mu[3]);

  ############################################################################################ Calc Full Vectors
  for(i in 1:N){
  Black_Armed[i]   <-inv_logit(Theta[i,1]); 
  White_Armed[i]   <-inv_logit(Theta[i,3]); 
  Hispanic_Armed[i]<-inv_logit(Theta[i,5]);

  Black_Unarmed[i]   <-inv_logit(Theta[i,2]); 
  White_Unarmed[i]   <-inv_logit(Theta[i,4]); 
  Hispanic_Unarmed[i]<-inv_logit(Theta[i,6]);


  RR_Black_Armed_Versus_Unarmed[i]    <-  inv_logit(Theta[i,1])/inv_logit(Theta[i,2]);
  RR_White_Armed_Versus_Unarmed[i]    <-  inv_logit(Theta[i,3])/inv_logit(Theta[i,4]);
  RR_Hispanic_Armed_Versus_Unarmed[i] <-  inv_logit(Theta[i,5])/inv_logit(Theta[i,6]);

  RR_Black_Armed_Versus_White_Armed[i]    <-  inv_logit(Theta[i,1])/inv_logit(Theta[i,3]);
  RR_Hispanic_Armed_Versus_White_Armed[i] <-  inv_logit(Theta[i,5])/inv_logit(Theta[i,3]);
  RR_Hispanic_Armed_Versus_Black_Armed[i] <-  inv_logit(Theta[i,5])/inv_logit(Theta[i,1]);

  RR_Black_Unarmed_Versus_White_Unarmed[i]   <-  inv_logit(Theta[i,2])/inv_logit(Theta[i,4]);
  RR_Hispanic_Unarmed_Versus_White_Unarmed[i]<-  inv_logit(Theta[i,6])/inv_logit(Theta[i,4]);
  RR_Hispanic_Unarmed_Versus_Black_Unarmed[i]<-  inv_logit(Theta[i,6])/inv_logit(Theta[i,2]);

  RR_Black_Unarmed_Versus_White_Armed[i]     <-  inv_logit(Theta[i,2])/inv_logit(Theta[i,3]); 
  RR_Hispanic_Unarmed_Versus_White_Armed[i]  <-  inv_logit(Theta[i,6])/inv_logit(Theta[i,3]);
  }
}"


################################################################################ Fit the Model IN STAN!

fitKilling <- pstan(model_code=model_code, data = model_dat, thin=1, iter = 4000, warmup=2000,chains = 6,refresh=1, control=list(metric="dense_e", max_treedepth=16))
