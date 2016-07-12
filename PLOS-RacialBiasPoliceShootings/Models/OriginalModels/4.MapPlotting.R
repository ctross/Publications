########################################################### Load Libraries
# Begin by opening R using the MapFile.RData file provided
library(rgdal) # Loads SP package by default
library(classInt)
library(RColorBrewer)
library(maps)
library(maptools)#for shapefiles
library(scales)  #for transparency
library(sp)
library(plotrix)

g<-read.csv(file.choose()) # MapFileData-WithCountyResultsAndCovariates.csv
gadm@data<-g               # Merge Data File and Shapefile

############################################### Begin Plotting Results and Data

################################## White Armed- Unarmed
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_White_Armed_Versus_Unarmed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)
   

   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]
 
   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")

                  

################################## Black Armed- Unarmed
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Black_Armed_Versus_Unarmed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


################################## Hispanic Armed- Unarmed
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Hispanic_Armed_Versus_Unarmed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")



################################## Armed Black to Armed White
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Black_Armed_Versus_White_Armed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")




################################## Armed Hispanic to Armed White
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Hispanic_Armed_Versus_White_Armed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


################################## Unarmed Black to Unarmed White
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Black_Unarmed_Versus_White_Unarmed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")




################################## Unarmed Hispanic to Unarmed White
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Hispanic_Unarmed_Versus_White_Unarmed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


 ################################## Unarmed Black to Armed White
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Black_Unarmed_Versus_White_Armed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")




################################## Unarmed Hispanic to Armed White
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(gadm@data$m.RR_Hispanic_Unarmed_Versus_White_Armed )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


 ################################## White Assault Rate
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(ifelse(gadm@data$AssaultsWhite.sum>1,gadm@data$AssaultsWhite.sum/gadm@data$WA_TOT + min((gadm@data$AssaultsWhite.sum/gadm@data$WA_TOT)[which((gadm@data$AssaultsWhite.sum/gadm@data$WA_TOT)>0)],na.rm=T),NA ))
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G)*10000,1))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


 ################################## Black Assault Rate
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(ifelse(gadm@data$AssaultsBlack.sum>1,gadm@data$AssaultsBlack.sum/gadm@data$BAC_TOT + min((gadm@data$AssaultsBlack.sum/gadm@data$BAC_TOT)[which((gadm@data$AssaultsBlack.sum/gadm@data$BAC_TOT)>0)],na.rm=T),NA ) )
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G)*10000,1))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")



 ################################## White Weapons Rate
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(ifelse(gadm@data$WeaponsWhite.sum>1,gadm@data$WeaponsWhite.sum/gadm@data$WA_TOT + min((gadm@data$WeaponsWhite.sum/gadm@data$WA_TOT)[which((gadm@data$WeaponsWhite.sum/gadm@data$WA_TOT)>0)],na.rm=T),NA ))
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G)*10000,1))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


 ################################## Black Weapons Rate
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(ifelse(gadm@data$WeaponsBlack.sum>1,gadm@data$WeaponsBlack.sum/gadm@data$BAC_TOT + min((gadm@data$WeaponsBlack.sum/gadm@data$BAC_TOT)[which((gadm@data$WeaponsBlack.sum/gadm@data$BAC_TOT)>0)],na.rm=T),NA ))
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G)*10000,1))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


################################## Google Hate Data
goog<-read.csv(file.choose())
      GoogleRacism<-goog$raciallychargedsearch
     GoogleHate<-c()

     for(t in 1:length(g$DMA)){
 GoogleHate[t] <-   ifelse(g$DMA[t]==156,NA,GoogleRacism[g$DMA[t]] );
 }
  GoogleHate[which(GoogleHate==999999999)]<-NA
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-(GoogleHate)
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round((G),1))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")


################################## Wealth
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(g$Median.Income)
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G)/1000,1))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")



################################## Gini
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-(g$Gini)
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round((G),2))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")



################################## Population
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(g$TOT_POP)
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G)/10000,3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")





################################## Black Population
par(mar=c(1,5,1,1))
  palette(c(brewer.pal(9, "YlOrRd")))
  XXX<-log(((g$BAC_TOT+1)/g$TOT_POP))
 # plot(gadm.state, xlim=c(-125,-65),ylim=c(35,42))
  plot(gadm, col=plotrix:::rescale(XXX,c(1,9)),  xlim=c(-125,-65),ylim=c(35,42), border=0)
   maps::map("state", interior = FALSE, col="black", lwd=1, add=TRUE)
   maps::map("state", boundary = FALSE, col="black",lty = 2, add = TRUE)


   ############ Scaling

   X<- rescale(XXX,c(1,9))
   X<-c(c(1:9),X)
   G <- rescale(X,range( XXX,na.rm=T))
   G<-G[1:9]

   testcol<-c("white",brewer.pal(9, "YlOrRd"))
   col.labels<-as.character(round(exp(G),3))
   col.labels<-c("No Data",col.labels)
   col.labels[2]<-paste0("<" ,col.labels[2])
   col.labels[10]<-paste0(">" ,col.labels[10])
   color.legend(-129,25,-127,40,col.labels,testcol,gradient="y")
