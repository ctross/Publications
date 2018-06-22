
path1 <- "C:\\Users\\cody_ross\\Dropbox\\Completed and Published Projects\\1-papers\\28-Bateman Pimbwe\\RSOS Submission\\PublicWorkflow"
setwd(path1)

load('PimbweBase.RData') # Builds a database from several data files, exports and R session

source('./Code/CollectModels.R') # Reads in all Stan models
source('./Code/RunModels.R') # Runs all Stan models on a server

source('./Code/Process-SpouseNumber.R')
source('./Code/Process-SpouseYears.R') 
source('./Code/Process-TimingWeights.R') 
source('./Code/Process-QualityWeights.R')   
source('./Code/Process-BothWeights.R') 
source('./Code/Process-Full.R')  

source('./Code/ResultsMain.R') 
source('./Code/ResultsDiffs.R')


