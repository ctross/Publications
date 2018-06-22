rm(list=ls())
load("Extrapolations.RData")
source('./Code/project_support.R')

P<-29

set.seed(1)

# Produces figures 5A, 5B, 5C
fignames_list <- c("PPP-PR-33-Hist.pdf", "PPP-PR-50-Hist.pdf", "PPP-PR-66-Hist.pdf")
fracs_list <- c(3, 2, 1.5)
phi_values <- c("=0.33", "=0.50", "=0.66")

for(j in 1:length(fignames_list)){
  GiniWealth <-  PercentWithCowives <- PercentRich <- WealthRatio <- matrix(999999, nrow=2000,ncol=P)

  for(i in 1:P){
    Scrap <- Results_Rival(AdjustedRival[[i]],fracs_list[j])
    GiniWealth[1:dim(Scrap)[1],i] <- Scrap[,1]
    PercentRich[1:dim(Scrap)[1],i] <- Scrap[,2]
    WealthRatio[1:dim(Scrap)[1],i] <- Scrap[,3]
  }

  GR <- GiniWealth
  PR <- PercentRich
  WR <- WealthRatio
  
  df2 <- read.csv("./Data/ExtraData/EarlyAgriculturalData.csv")
 # df2<-df2[which(df2$Proxy != "Income"),]
  pr_hist <- df2[,c("Phi33", "Phi50", "Phi66")]

  d <- read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")
  Subsist <- d$Subsistence
  Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )
  Subsist <- as.character(Subsist)

  PRH <- c(colMeans(PR),pr_hist[,j])
  Subsist<-c(Subsist,rep("Agricultural",length(pr_hist[,j])))

  # Check t test for each frame
  qq1 <- t.test( PRH[which(Subsist=="Forager")], PRH[which(Subsist=="Agricultural")],conf.level = 0.9) #
  qq2 <- t.test( PRH[which(Subsist=="Horticultural")], PRH[which(Subsist=="Agricultural")],conf.level = 0.9)
  qq3 <- t.test( PRH[which(Subsist=="Agropastoral")], PRH[which(Subsist=="Agricultural")],conf.level = 0.9)

  Subsist <- ifelse(Subsist=="Forager",paste0("Forager: ",paste0(round(qq1$estimate[1]-qq1$estimate[2],2)," (", round(qq1$conf.int[1],2),", ",round(qq1$conf.int[2],2),")" )),Subsist)
  Subsist <- ifelse(Subsist=="Horticultural",paste0("Horticultural: ",paste0(round(qq2$estimate[1]-qq2$estimate[2],2)," (", round(qq2$conf.int[1],2),", ",round(qq2$conf.int[2],2),")" )),Subsist)
  Subsist <- ifelse(Subsist=="Agropastoral",paste0("Agropastoral: ",paste0(round(qq3$estimate[1]-qq3$estimate[2],2)," (", round(qq3$conf.int[1],2),", ",round(qq3$conf.int[2],2),")" )),Subsist)
  Subsist <- ifelse(Subsist=="Agricultural","Agricultural",Subsist)

  Subsist <- factor(Subsist, levels= c(paste0("Forager: ",paste0(round(qq1$estimate[1]-qq1$estimate[2],2)," (", round(qq1$conf.int[1],2),", ",round(qq1$conf.int[2],2),")" )),
                                    paste0("Horticultural: ",paste0(round(qq2$estimate[1]-qq2$estimate[2],2)," (", round(qq2$conf.int[1],2),", ",round(qq2$conf.int[2],2),")" )),
                                    paste0("Agropastoral: ",paste0(round(qq3$estimate[1]-qq3$estimate[2],2)," (", round(qq3$conf.int[1],2),", ",round(qq3$conf.int[2],2),")" )),
                                    "Agricultural")   )

  df<-data.frame(PercentRich=PRH, Subsistence=Subsist)
  
 

  p <- ggplot(df, aes(Subsistence, PercentRich)) + labs(y=bquote("Frequency of rich males, " ~ phi == .(phi_values[j]))) + theme(legend.title=element_text(size=14),legend.text=element_text(size=14)) +
  theme(strip.text.x = element_text(size=14),axis.text=element_text(size=12), axis.title=element_text(size=14)) +
  geom_point(aes(colour = Subsistence),size=3) + scale_colour_manual("Mean difference between focal\n subsistence type and agriculture:",values=cbbPalette[c(5,1,3,2)])  + theme( axis.text.x=element_blank(),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank())
  p

  ggsave(fignames_list[j] , p, width=5, height=4)

  print(fignames_list[j])
}

                       

