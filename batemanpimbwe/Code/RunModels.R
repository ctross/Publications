
Iter <- 4000
Warmup <- 2000   
Thin <- 1
Init <- 0
Cores <- 2
Chains <- 2
Refresh <- 1
Seed <- 123
MTD <- 15
AD <- 0.96
Metric <- "diag_e"

library(rethinking)
library(parallel)

result1 <- mclapply(1:6,function(z){
 stan(model_code=model_code[[z]], data = model_dat,thin=Thin, iter = Iter, warmup=Warmup, cores=Cores, chains=Chains, refresh=Refresh,init=Init,control=list(max_treedepth=MTD, adapt_delta=AD,metric=Metric))
 },
 mc.cores = 6*Cores)    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 