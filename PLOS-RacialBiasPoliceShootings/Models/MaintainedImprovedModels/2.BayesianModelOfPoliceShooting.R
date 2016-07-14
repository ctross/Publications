
################################## Load Data ###################################
 d<-read.csv(file.choose()) # MapFileData-WithCountyResultsAndCovariates.csv
 library(rstan)
 library(rethinking)

##############################################################################
# This code is used to model the relative risk across race-ethnicity
# and armed-unarmed status

############### The following line of code just limits the data set to counties
# in which we have data on police shooting.
d<-d[which(complete.cases(d$BlackArmed)),]

############### Extract population size data for each county
Nblack<-d$BAC_TOT
Nwhite<-d$WA_TOT
Nhispanic<-d$H_TOT

############### Extract shooting data for each county
UnarmedBlack<-d$BlackUnarmed	
ArmedBlack<-d$BlackArmed
UnarmedHispanic<-d$HispanicUnarmed	
ArmedHispanic<-d$HispanicArmed		
UnarmedWhite<-d$WhiteUnarmed		
ArmedWhite<-d$WhiteArmed

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

fitKilling <- stan(model_code=model_code, data = model_dat,init=0, thin=1, iter = 2000, warmup=1000,chains = 1,refresh=1)

print(fitKilling,digits_summary=4,pars=c("Mu"))

print(fitKilling,digits_summary=4,pars=c(  "Mu_RR_Black_Armed_Versus_Unarmed",
                                           "Mu_RR_White_Armed_Versus_Unarmed",           
                                           "Mu_RR_Hispanic_Armed_Versus_Unarmed",        
                                                                                            
                                           "Mu_RR_Black_Armed_Versus_White_Armed",       
                                           "Mu_RR_Hispanic_Armed_Versus_White_Armed",    
                                           "Mu_RR_Hispanic_Armed_Versus_Black_Armed",    
                                                                                            
                                           "Mu_RR_Black_Unarmed_Versus_White_Unarmed",   
                                           "Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed",
                                           "Mu_RR_Hispanic_Unarmed_Versus_Black_Unarmed",
                                                                                            
                                           "Mu_RR_Black_Unarmed_Versus_White_Armed",     
                                           "Mu_RR_Hispanic_Unarmed_Versus_White_Armed"))  
    
################################################################################################################################################ Post-Processing
RR_Black_Armed_Versus_Unarmed                <-extract(fitKilling,pars="RR_Black_Armed_Versus_Unarmed"            )$RR_Black_Armed_Versus_Unarmed
RR_White_Armed_Versus_Unarmed                <-extract(fitKilling,pars="RR_White_Armed_Versus_Unarmed"            )$RR_White_Armed_Versus_Unarmed
RR_Hispanic_Armed_Versus_Unarmed             <-extract(fitKilling,pars="RR_Hispanic_Armed_Versus_Unarmed"         )$RR_Hispanic_Armed_Versus_Unarmed
                                                                                      
RR_Black_Armed_Versus_White_Armed            <-extract(fitKilling,pars="RR_Black_Armed_Versus_White_Armed"        )$RR_Black_Armed_Versus_White_Armed
RR_Hispanic_Armed_Versus_White_Armed         <-extract(fitKilling,pars="RR_Hispanic_Armed_Versus_White_Armed"     )$RR_Hispanic_Armed_Versus_White_Armed
                                                                                      
                                                                                      
RR_Black_Unarmed_Versus_White_Unarmed        <-extract(fitKilling,pars="RR_Black_Unarmed_Versus_White_Unarmed"    )$RR_Black_Unarmed_Versus_White_Unarmed
RR_Hispanic_Unarmed_Versus_White_Unarmed     <-extract(fitKilling,pars="RR_Hispanic_Unarmed_Versus_White_Unarmed" )$RR_Hispanic_Unarmed_Versus_White_Unarmed
                                                                                      
RR_Black_Unarmed_Versus_White_Armed          <-extract(fitKilling,pars="RR_Black_Unarmed_Versus_White_Armed"      )$RR_Black_Unarmed_Versus_White_Armed
RR_Hispanic_Unarmed_Versus_White_Armed       <-extract(fitKilling,pars="RR_Hispanic_Unarmed_Versus_White_Armed"   )$RR_Hispanic_Unarmed_Versus_White_Armed

Mu_RR_Black_Armed_Versus_Unarmed             <-extract(fitKilling,pars="Mu_RR_Black_Armed_Versus_Unarmed"            )$Mu_RR_Black_Armed_Versus_Unarmed            
Mu_RR_White_Armed_Versus_Unarmed             <-extract(fitKilling,pars="Mu_RR_White_Armed_Versus_Unarmed"            )$Mu_RR_White_Armed_Versus_Unarmed            
Mu_RR_Hispanic_Armed_Versus_Unarmed          <-extract(fitKilling,pars="Mu_RR_Hispanic_Armed_Versus_Unarmed"         )$Mu_RR_Hispanic_Armed_Versus_Unarmed         
                                                                                      
Mu_RR_Black_Armed_Versus_White_Armed         <-extract(fitKilling,pars="Mu_RR_Black_Armed_Versus_White_Armed"        )$Mu_RR_Black_Armed_Versus_White_Armed        
Mu_RR_Hispanic_Armed_Versus_White_Armed      <-extract(fitKilling,pars="Mu_RR_Hispanic_Armed_Versus_White_Armed"     )$Mu_RR_Hispanic_Armed_Versus_White_Armed     
                                                                                      
                                                                                      
Mu_RR_Black_Unarmed_Versus_White_Unarmed     <-extract(fitKilling,pars="Mu_RR_Black_Unarmed_Versus_White_Unarmed"    )$Mu_RR_Black_Unarmed_Versus_White_Unarmed    
Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed  <-extract(fitKilling,pars="Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed" )$Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed 
                                                                                      
Mu_RR_Black_Unarmed_Versus_White_Armed       <-extract(fitKilling,pars="Mu_RR_Black_Unarmed_Versus_White_Armed"      )$Mu_RR_Black_Unarmed_Versus_White_Armed      
Mu_RR_Hispanic_Unarmed_Versus_White_Armed    <-extract(fitKilling,pars="Mu_RR_Hispanic_Unarmed_Versus_White_Armed"   )$Mu_RR_Hispanic_Unarmed_Versus_White_Armed   


m.RR_Black_Armed_Versus_Unarmed               <-c()  
m.RR_White_Armed_Versus_Unarmed               <-c()  
m.RR_Hispanic_Armed_Versus_Unarmed            <-c()  
                                        
m.RR_Black_Armed_Versus_White_Armed           <-c()
m.RR_Hispanic_Armed_Versus_White_Armed        <-c()
                                        
                                        
m.RR_Black_Unarmed_Versus_White_Unarmed       <-c()
m.RR_Hispanic_Unarmed_Versus_White_Unarmed    <-c()
                                        
m.RR_Black_Unarmed_Versus_White_Armed         <-c()
m.RR_Hispanic_Unarmed_Versus_White_Armed      <-c()

pc95l.RR_Black_Armed_Versus_Unarmed              <-c()
pc95l.RR_White_Armed_Versus_Unarmed              <-c()
pc95l.RR_Hispanic_Armed_Versus_Unarmed           <-c()
                                              
pc95l.RR_Black_Armed_Versus_White_Armed          <-c()
pc95l.RR_Hispanic_Armed_Versus_White_Armed       <-c()
                                              
                                              
pc95l.RR_Black_Unarmed_Versus_White_Unarmed      <-c()
pc95l.RR_Hispanic_Unarmed_Versus_White_Unarmed   <-c()
                                              
pc95l.RR_Black_Unarmed_Versus_White_Armed        <-c()
pc95l.RR_Hispanic_Unarmed_Versus_White_Armed     <-c()

pc95h.RR_Black_Armed_Versus_Unarmed              <-c()
pc95h.RR_White_Armed_Versus_Unarmed              <-c()
pc95h.RR_Hispanic_Armed_Versus_Unarmed           <-c()
                                              
pc95h.RR_Black_Armed_Versus_White_Armed          <-c()
pc95h.RR_Hispanic_Armed_Versus_White_Armed       <-c()
                                              
                                              
pc95h.RR_Black_Unarmed_Versus_White_Unarmed      <-c()
pc95h.RR_Hispanic_Unarmed_Versus_White_Unarmed   <-c()
                                              
pc95h.RR_Black_Unarmed_Versus_White_Armed        <-c()
pc95h.RR_Hispanic_Unarmed_Versus_White_Armed     <-c()
 
 
 
for(i in 1:N){
m.RR_Black_Armed_Versus_Unarmed[i]<-median(RR_Black_Armed_Versus_Unarmed[,i])
m.RR_White_Armed_Versus_Unarmed[i]<-median(RR_White_Armed_Versus_Unarmed[,i])
m.RR_Hispanic_Armed_Versus_Unarmed[i]<-median(RR_Hispanic_Armed_Versus_Unarmed[,i])

m.RR_Black_Armed_Versus_White_Armed[i]<-median(RR_Black_Armed_Versus_White_Armed[,i])
m.RR_Hispanic_Armed_Versus_White_Armed[i]<-median(RR_Hispanic_Armed_Versus_White_Armed[,i])

m.RR_Black_Unarmed_Versus_White_Unarmed[i]<-median(RR_Black_Unarmed_Versus_White_Unarmed[,i])
m.RR_Hispanic_Unarmed_Versus_White_Unarmed[i]<-median(RR_Hispanic_Unarmed_Versus_White_Unarmed[,i])

m.RR_Black_Unarmed_Versus_White_Armed[i]<-median(RR_Black_Unarmed_Versus_White_Armed[,i])
m.RR_Hispanic_Unarmed_Versus_White_Armed[i]<-median(RR_Hispanic_Unarmed_Versus_White_Armed[,i])

m.log.RR_Black_Unarmed_Versus_White_Unarmed  <- mean(log(RR_Black_Unarmed_Versus_White_Unarmed[,i]))
m.log.RR_Hispanic_Unarmed_Versus_White_Unarmed <-mean(log(RR_Hispanic_Unarmed_Versus_White_Unarmed[,i]))
m.log.RR_Black_Unarmed_Versus_White_Armed     <- mean(log(RR_Black_Unarmed_Versus_White_Armed[,i]))
m.log.RR_Hispanic_Unarmed_Versus_White_Armed   <- mean(log(RR_Hispanic_Unarmed_Versus_White_Armed[,i]))

sd.log.RR_Black_Unarmed_Versus_White_Unarmed   <- sd(log(RR_Black_Unarmed_Versus_White_Unarmed[,i]))
sd.log.RR_Hispanic_Unarmed_Versus_White_Unarmed <-sd(log(RR_Hispanic_Unarmed_Versus_White_Unarmed[,i]))
sd.log.RR_Black_Unarmed_Versus_White_Armed     <- sd(log(RR_Black_Unarmed_Versus_White_Armed[,i]))
sd.log.RR_Hispanic_Unarmed_Versus_White_Armed  <- sd(log(RR_Hispanic_Unarmed_Versus_White_Armed[,i]))

pc95h.RR_Black_Armed_Versus_Unarmed[i]<-PCI(RR_Black_Armed_Versus_Unarmed[,i])[2]
pc95h.RR_White_Armed_Versus_Unarmed[i]<-PCI(RR_White_Armed_Versus_Unarmed[,i])[2]
pc95h.RR_Hispanic_Armed_Versus_Unarmed[i]<-PCI(RR_Hispanic_Armed_Versus_Unarmed[,i])[2]
                                                                                            
pc95h.RR_Black_Armed_Versus_White_Armed[i]<-PCI(RR_Black_Armed_Versus_White_Armed[,i])[2]
pc95h.RR_Hispanic_Armed_Versus_White_Armed[i]<-PCI(RR_Hispanic_Armed_Versus_White_Armed[,i])[2]
                                                                                            
                                                                                            
pc95h.RR_Black_Unarmed_Versus_White_Unarmed[i]<-PCI(RR_Black_Unarmed_Versus_White_Unarmed[,i])[2]
pc95h.RR_Hispanic_Unarmed_Versus_White_Unarmed[i]<-PCI(RR_Hispanic_Unarmed_Versus_White_Unarmed[,i])[2]
                                                                                            
pc95h.RR_Black_Unarmed_Versus_White_Armed[i]<-PCI(RR_Black_Unarmed_Versus_White_Armed[,i])[2]
pc95h.RR_Hispanic_Unarmed_Versus_White_Armed[i]<-PCI(RR_Hispanic_Unarmed_Versus_White_Armed[,i])[2]

pc95l.RR_Black_Armed_Versus_Unarmed[i]<-PCI(RR_Black_Armed_Versus_Unarmed[,i])[1]
pc95l.RR_White_Armed_Versus_Unarmed[i]<-PCI(RR_White_Armed_Versus_Unarmed[,i])[1]
pc95l.RR_Hispanic_Armed_Versus_Unarmed[i]<-PCI(RR_Hispanic_Armed_Versus_Unarmed[,i])[1]
                                                                                            
pc95l.RR_Black_Armed_Versus_White_Armed[i]<-PCI(RR_Black_Armed_Versus_White_Armed[,i])[1]
pc95l.RR_Hispanic_Armed_Versus_White_Armed[i]<-PCI(RR_Hispanic_Armed_Versus_White_Armed[,i])[1]
                                                                                            
                                                                                            
pc95l.RR_Black_Unarmed_Versus_White_Unarmed[i]<-PCI(RR_Black_Unarmed_Versus_White_Unarmed[,i])[1]
pc95l.RR_Hispanic_Unarmed_Versus_White_Unarmed[i]<-PCI(RR_Hispanic_Unarmed_Versus_White_Unarmed[,i])[1]
                                                                                            
pc95l.RR_Black_Unarmed_Versus_White_Armed[i]<-PCI(RR_Black_Unarmed_Versus_White_Armed[,i])[1]
pc95l.RR_Hispanic_Unarmed_Versus_White_Armed[i]<-PCI(RR_Hispanic_Unarmed_Versus_White_Armed[,i])[1]
} 


Results <- data.frame(
m.RR_Black_Armed_Versus_Unarmed               ,  
m.RR_White_Armed_Versus_Unarmed               ,  
m.RR_Hispanic_Armed_Versus_Unarmed            ,  
                                                     
m.RR_Black_Armed_Versus_White_Armed           ,  
m.RR_Hispanic_Armed_Versus_White_Armed        ,  
                                                     
                                                     
m.RR_Black_Unarmed_Versus_White_Unarmed       ,  
m.RR_Hispanic_Unarmed_Versus_White_Unarmed    ,  
                                                     
m.RR_Black_Unarmed_Versus_White_Armed         ,  
m.RR_Hispanic_Unarmed_Versus_White_Armed      ,  
                                                     
pc95l.RR_Black_Armed_Versus_Unarmed              ,
pc95l.RR_White_Armed_Versus_Unarmed              ,
pc95l.RR_Hispanic_Armed_Versus_Unarmed           ,
                                                     
pc95l.RR_Black_Armed_Versus_White_Armed          ,
pc95l.RR_Hispanic_Armed_Versus_White_Armed       ,
                                                     
                                                     
pc95l.RR_Black_Unarmed_Versus_White_Unarmed      ,
pc95l.RR_Hispanic_Unarmed_Versus_White_Unarmed   ,
                                                     
pc95l.RR_Black_Unarmed_Versus_White_Armed        ,
pc95l.RR_Hispanic_Unarmed_Versus_White_Armed     ,

pc95h.RR_Black_Armed_Versus_Unarmed              ,
pc95h.RR_White_Armed_Versus_Unarmed              ,
pc95h.RR_Hispanic_Armed_Versus_Unarmed           ,

pc95h.RR_Black_Armed_Versus_White_Armed          ,
pc95h.RR_Hispanic_Armed_Versus_White_Armed       ,


pc95h.RR_Black_Unarmed_Versus_White_Unarmed      ,
pc95h.RR_Hispanic_Unarmed_Versus_White_Unarmed   ,

pc95h.RR_Black_Unarmed_Versus_White_Armed        ,
pc95h.RR_Hispanic_Unarmed_Versus_White_Armed   

m.log.RR_Black_Unarmed_Versus_White_Unarmed      ,
m.log.RR_Hispanic_Unarmed_Versus_White_Unarmed	 ,
m.log.RR_Black_Unarmed_Versus_White_Armed	       ,
m.log.RR_Hispanic_Unarmed_Versus_White_Armed	   ,

sd.log.RR_Black_Unarmed_Versus_White_Unarmed	   ,
sd.log.RR_Hispanic_Unarmed_Versus_White_Unarmed	 ,
sd.log.RR_Black_Unarmed_Versus_White_Armed	     ,
sd.log.RR_Hispanic_Unarmed_Versus_White_Armed)



write.csv(cbind(d,Results),"USPSD-Results.csv")





############################################################## Caterpillar Plot
m.RR1   <- m.RR_Black_Unarmed_Versus_White_Unarmed
pch.RR1 <- pc95h.RR_Black_Unarmed_Versus_White_Unarmed
pcl.RR1 <- pc95l.RR_Black_Unarmed_Versus_White_Unarmed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Black_Unarmed_Versus_White_Unarmed)
PL <- PCI(Mu_RR_Black_Unarmed_Versus_White_Unarmed)[1]
PH <- PCI(Mu_RR_Black_Unarmed_Versus_White_Unarmed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Black-Unarmed to White-Unarmed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )
  
  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3) 
  abline(h=log(32),col="black",lty=3) 
  abline(h=log(16),col="black",lty=3) 
  abline(h=log(8),col="black",lty=3) 
  abline(h=log(4),col="black",lty=3) 
  abline(h=log(2),col="black",lty=3) 
  abline(h=log(1),col="black",lty=1) 
  abline(h=log(.5),col="black",lty=3) 
   abline(h=log(.25),col="black",lty=3) 
      abline(h=log(.125),col="black",lty=3)
  
  
  abline(h=0,col="black",lty=3) 
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )
  
  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )





############################################################## Caterpillar Plot
m.RR1   <- m.RR_Hispanic_Unarmed_Versus_White_Unarmed
pch.RR1 <- pc95h.RR_Hispanic_Unarmed_Versus_White_Unarmed
pcl.RR1 <- pc95l.RR_Hispanic_Unarmed_Versus_White_Unarmed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed)
PL <- PCI(Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed)[1]
PH <- PCI(Mu_RR_Hispanic_Unarmed_Versus_White_Unarmed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Hispanic-Unarmed to White-Unarmed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )




############################################################## Caterpillar Plot
m.RR1   <- m.RR_Hispanic_Armed_Versus_White_Armed
pch.RR1 <- pc95h.RR_Hispanic_Armed_Versus_White_Armed
pcl.RR1 <- pc95l.RR_Hispanic_Armed_Versus_White_Armed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Hispanic_Armed_Versus_White_Armed)
PL <- PCI(Mu_RR_Hispanic_Armed_Versus_White_Armed)[1]
PH <- PCI(Mu_RR_Hispanic_Armed_Versus_White_Armed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Hispanic-Armed to White-Armed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )



############################################################## Caterpillar Plot
m.RR1   <- m.RR_Black_Armed_Versus_White_Armed
pch.RR1 <- pc95h.RR_Black_Armed_Versus_White_Armed
pcl.RR1 <- pc95l.RR_Black_Armed_Versus_White_Armed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Black_Armed_Versus_White_Armed)
PL <- PCI(Mu_RR_Black_Armed_Versus_White_Armed)[1]
PH <- PCI(Mu_RR_Black_Armed_Versus_White_Armed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Black-Armed to White-Armed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )


############################################################## Caterpillar Plot
m.RR1   <- m.RR_Black_Unarmed_Versus_White_Armed
pch.RR1 <- pc95h.RR_Black_Unarmed_Versus_White_Armed
pcl.RR1 <- pc95l.RR_Black_Unarmed_Versus_White_Armed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Black_Unarmed_Versus_White_Armed)
PL <- PCI(Mu_RR_Black_Unarmed_Versus_White_Armed)[1]
PH <- PCI(Mu_RR_Black_Unarmed_Versus_White_Armed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Black-Unarmed to White-Armed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )


  ############################################################## Caterpillar Plot
m.RR1   <- m.RR_Hispanic_Unarmed_Versus_White_Armed
pch.RR1 <- pc95h.RR_Hispanic_Unarmed_Versus_White_Armed
pcl.RR1 <- pc95l.RR_Hispanic_Unarmed_Versus_White_Armed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Hispanic_Unarmed_Versus_White_Armed)
PL <- PCI(Mu_RR_Hispanic_Unarmed_Versus_White_Armed)[1]
PH <- PCI(Mu_RR_Hispanic_Unarmed_Versus_White_Armed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Hispanic-Unarmed to White-Armed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )


    ############################################################## Caterpillar Plot
m.RR1   <- m.RR_Hispanic_Armed_Versus_Unarmed
pch.RR1 <- pc95h.RR_Hispanic_Armed_Versus_Unarmed
pcl.RR1 <- pc95l.RR_Hispanic_Armed_Versus_Unarmed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Hispanic_Armed_Versus_Unarmed)
PL <- PCI(Mu_RR_Hispanic_Armed_Versus_Unarmed)[1]
PH <- PCI(Mu_RR_Hispanic_Armed_Versus_Unarmed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Hispanic-Armed to Hispanic-Unarmed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )

  ############################################################## Caterpillar Plot
m.RR1   <- m.RR_Black_Armed_Versus_Unarmed
pch.RR1 <- pc95h.RR_Black_Armed_Versus_Unarmed
pcl.RR1 <- pc95l.RR_Black_Armed_Versus_Unarmed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_Black_Armed_Versus_Unarmed)
PL <- PCI(Mu_RR_Black_Armed_Versus_Unarmed)[1]
PH <- PCI(Mu_RR_Black_Armed_Versus_Unarmed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: Black-Armed to Black-Unarmed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )

  ############################################################## Caterpillar Plot
m.RR1   <- m.RR_White_Armed_Versus_Unarmed
pch.RR1 <- pc95h.RR_White_Armed_Versus_Unarmed
pcl.RR1 <- pc95l.RR_White_Armed_Versus_Unarmed

m.RR   <- m.RR1[order(m.RR1)]
pch.RR <- pch.RR1[order(m.RR1)]
pcl.RR <- pcl.RR1[order(m.RR1)]

Mu <- median(Mu_RR_White_Armed_Versus_Unarmed)
PL <- PCI(Mu_RR_White_Armed_Versus_Unarmed)[1]
PH <- PCI(Mu_RR_White_Armed_Versus_Unarmed)[2]

plot(c(log(m.RR),log(Mu))~c(1:(length(m.RR)+1)), type="n",xlab="County-Level Ratio of Probabilities of Being Shot By Police: White-Armed to White-Unarmed",
 ylab="Probability Ratio (Natural Units, Log-scale)", xaxt="n" ,  yaxt="n" ,
  pch=20, col="indianred",ylim=c(min(-2),max(c(5))) )

  x<-c(log(128),log(64),log(32),log(16),log(8),log(4),log(2),log(1),log(0.5),log(0.25),log(0.125))
  y<-c(128,64,32,16,8,4,2,1,0.5,0.25,0.125)
  axis(2, at=x,labels=y)

   abline(h=log(128),col="black",lty=3)
  abline(h=log(64),col="black",lty=3)
  abline(h=log(32),col="black",lty=3)
  abline(h=log(16),col="black",lty=3)
  abline(h=log(8),col="black",lty=3)
  abline(h=log(4),col="black",lty=3)
  abline(h=log(2),col="black",lty=3)
  abline(h=log(1),col="black",lty=1)
  abline(h=log(.5),col="black",lty=3)
   abline(h=log(.25),col="black",lty=3)
      abline(h=log(.125),col="black",lty=3)


  abline(h=0,col="black",lty=3)
  segments(c(1:(length(pch.RR))),log(pch.RR),c(1:length(pch.RR)), log(pcl.RR),lwd=2,col="darkgrey",lend=1)
  points(log(m.RR)~c(1:length(pch.RR)),    pch=20, col="black",cex=1.5 )

  segments((length(pch.RR)+1),log(PL),(length(pch.RR)+1), log(PH),lwd=3,col="slateblue",lend=1)
  points((length(pch.RR)+1),log(Mu) ,   pch=20, col="indianred",cex=1.5 )




