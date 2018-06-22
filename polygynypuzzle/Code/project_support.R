
library(rgdal)
library(classInt)
library(RColorBrewer)
library(rethinking) # dont put rethinking after maps, ggplot2 tries to call maps::map
library(maps)
library(maptools)
library(scales)
library(sp)
library(plotrix)
library(ggmap)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(ineq)
library(directlabels)
library(scales)
library(MASS)
library(xtable)

set.seed(8675309)

cbbPalette <- c(
  "#000000", 
  "#E69F00", 
  "#009E73", 
  "#0072B2", 
  "indianred", 
  "#CC79A7")


Is_Rich <- function(Z, fracs=2){
  is_rich <- !cumsum(sort(Z,decreasing=TRUE))/(sum(Z)/fracs) > 1 
  is_rich
}

Perc_Rich <- function(Z, fracs=2){
  RES <- mean(Is_Rich(Z, fracs=fracs))
  return(RES)
}

Wealth_Ratio <- function(Z, fracs=2){
  wealth_vector <- sort(Z, decreasing=TRUE)
  numerator <- mean( wealth_vector[which( Is_Rich(Z, fracs=fracs))] )
  denominator <- mean( wealth_vector[which( !Is_Rich(Z, fracs=fracs))] )
  RES <- numerator / denominator
  return(RES)
}

Pct_Male_Poly <- function(Z){
  RES<-length(Z[which(Z>1)])/length(Z)
  RES
}

Pct_Female_Poly <- function(Z){
  RES<-sum(Z[which(Z>1)])/sum(Z)
  RES
}

Results_Wives <- function(M){
  RES <- matrix(NA,ncol=2,nrow=dim(M)[1])
  for(i in 1:dim(M)[1]){
    X <- M[i,]
    RES[i,] <-c(Pct_Female_Poly(X),Pct_Male_Poly(X))
  }
  RES
}


Results_Rival <- function(M,F){
  RES <- matrix(NA,ncol=3,nrow=dim(M)[1])
  for(i in 1:nrow(M)){
    X <- M[i,]
    RES[i,] <- c(Gini(X), Perc_Rich(X,fracs=F), Wealth_Ratio(X,fracs=F))
  }
  RES
}

Perc_Rich2 <- function(Z,W){
  dat<-data.frame(factor=round(W,0), value=Z)
  if(length(which(dat$factor==1))==0){
    RES <- NA
  }else{
    thresh1<-mean(  with(dat, tapply(value, factor==2, mean))[2] ,
         with(dat, tapply(value, factor==1, mean))[2] )

    RES <- length(which(Z>thresh1))/length(Z)
    RES
  } 
}

Wealth_Ratio2 <- function(Z,W){
  dat<-data.frame(factor=round(W,0), value=Z)
  if(length(which(dat$factor==1))==0){
    RES <- NA
  }else{
    dat<-data.frame(factor=round(W,0), value=Z)
    thresh1<-mean(  with(dat, tapply(value, factor==2, mean))[2] ,
         with(dat, tapply(value, factor==1, mean))[2] )
    RES <- mean(Z[which(Z>thresh1)])/mean(Z[which(Z<=thresh1)])
    RES
  }
}

Results_Rival2 <- function(M,Q){
  RES <- matrix(NA,ncol=3,nrow=dim(M)[1])
  for(i in 1:dim(M)[1]){
    X <- M[i,]
    Q2 <- Q[i,]
    RES[i,] <-c(Gini(X),Perc_Rich2(X,Q2),Wealth_Ratio2(X,Q2))
  }
  RES
}


prettyM2 <- function(X){
  RES <- sprintf("%.1f",round(median(na.omit(X)),2))
  return(RES)
}

prettyCI2 <- function(X,Y){
  RES <- sprintf("%.1f",round(PCI(na.omit(X),0.9)[Y],2))
  return(RES)
}


prettyM <- function(X){
  RES <- sprintf("%.2f",round(median(na.omit(X)),2))
  return(RES)
}

prettyCI <- function(X,Y){
  RES <- sprintf("%.2f",round(PCI(na.omit(X),0.9)[Y],2))
  return(RES)
}

RMC <- function(X,Y){
  RES <- rowMeans(na.omit(X[,which(Subsist==Y)]))
  return(RES)
}