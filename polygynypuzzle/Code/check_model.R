###############################################################################
# Check Model

rm(list=ls())

load("ModelResults.RData")
all_stuff <- ls()
rm( list=setdiff(all_stuff, "fitPoly") )
source('./Code/project_support.R')


print(fitPoly)
set.seed(714)
GG<-traceplot(fitPoly,pars=sample(names(fitPoly@sim$samples[[1]]),30))

ggsave("TracePPP.pdf",GG,width=11,height=8.5)
