rm(list=ls())
load("Extrapolations.RData")
source('./Code/project_support.R')

P<-29

set.seed(1)

# Produces figures 5A, 5B, 5C
fignames_list <- c("PPP-PR-33.pdf", "PPP-PR-50.pdf", "PPP-PR-66.pdf")
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

  d <- read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")
  Subsist <- d$Subsistence
  Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )
  Subsist <- as.character(Subsist)

  # Check t test for each frame
  qq1 <- t.test( colMeans(PR)[which(Subsist=="Forager")], colMeans(PR)[which(Subsist=="Agricultural")],conf.level = 0.9) #
  qq2 <- t.test( colMeans(PR)[which(Subsist=="Horticultural")], colMeans(PR)[which(Subsist=="Agricultural")],conf.level = 0.9)
  qq3 <- t.test( colMeans(PR)[which(Subsist=="Agropastoral")], colMeans(PR)[which(Subsist=="Agricultural")],conf.level = 0.9)

  Subsist <- ifelse(Subsist=="Forager",paste0("Forager: ",paste0(round(qq1$estimate[1]-qq1$estimate[2],2)," (", round(qq1$conf.int[1],2),", ",round(qq1$conf.int[2],2),")" )),Subsist)
  Subsist <- ifelse(Subsist=="Horticultural",paste0("Horticultural: ",paste0(round(qq2$estimate[1]-qq2$estimate[2],2)," (", round(qq2$conf.int[1],2),", ",round(qq2$conf.int[2],2),")" )),Subsist)
  Subsist <- ifelse(Subsist=="Agropastoral",paste0("Agropastoral: ",paste0(round(qq3$estimate[1]-qq3$estimate[2],2)," (", round(qq3$conf.int[1],2),", ",round(qq3$conf.int[2],2),")" )),Subsist)
  Subsist <- ifelse(Subsist=="Agricultural","Agricultural",Subsist)

  Subsist <- factor(Subsist, levels= c(paste0("Forager: ",paste0(round(qq1$estimate[1]-qq1$estimate[2],2)," (", round(qq1$conf.int[1],2),", ",round(qq1$conf.int[2],2),")" )),
                                    paste0("Horticultural: ",paste0(round(qq2$estimate[1]-qq2$estimate[2],2)," (", round(qq2$conf.int[1],2),", ",round(qq2$conf.int[2],2),")" )),
                                    paste0("Agropastoral: ",paste0(round(qq3$estimate[1]-qq3$estimate[2],2)," (", round(qq3$conf.int[1],2),", ",round(qq3$conf.int[2],2),")" )),
                                    "Agricultural")   )

  df<-data.frame(PercentRich=colMeans(PR), Subsistence=Subsist)

  p <- ggplot(df, aes(Subsistence, PercentRich)) + labs(y=bquote("Frequency of rich males, " ~ phi == .(phi_values[j]))) + theme(legend.title=element_text(size=14),legend.text=element_text(size=14)) +
  theme(strip.text.x = element_text(size=14),axis.text=element_text(size=12), axis.title=element_text(size=14)) +
  geom_point(aes(colour = Subsistence),size=3) + scale_colour_manual("Mean difference between focal\n subsistence type and agriculture:",values=cbbPalette[c(5,1,3,2)])  + theme( axis.text.x=element_blank(),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank())
  p

  ggsave(fignames_list[j] , p, width=5, height=4)

  print(fignames_list[j])
}

                       


################################################################################
# Fig 1 c
################################################################################

fracs<-2

GiniWealth <-  PercentWithCowives <- PercentRich <- WealthRatio <- matrix(999999, nrow=2000,ncol=P)

for(i in 1:P){
  Scrap <- Results_Rival(AdjustedRival[[i]],fracs)
  GiniWealth[1:nrow(Scrap),i] <- Scrap[,1]
  PercentRich[1:nrow(Scrap),i] <- Scrap[,2]
  WealthRatio[1:nrow(Scrap),i] <- Scrap[,3]
}

GR <- GiniWealth
PR <- PercentRich
WR <- WealthRatio

d <- read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")
Subsist <- d$Subsistence
Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )

myData <- aggregate(colMeans(GR), by = list(Subsistence = Subsist), FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))

dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = myData$mean + 2*myData$se,
      ymin = myData$mean - 2*myData$se)

myData <- do.call(data.frame, myData)

myData$se <- myData$x.sd / sqrt(myData$x.n)

colnames(myData) <- c("Subsistence", "mean", "sd", "n", "se")

myData$names <- c("Foraging (8)", "Horticulture (9)", "Agropastoral (8)","Agriculture (4)"  )

p <- ggplot(data = myData, aes(x = factor(names,levels = myData$names), y = mean, fill = factor(names,levels = myData$names)))

p <- p + geom_bar(stat = "identity", position = dodge) +    theme_bw() +
geom_errorbar(limits, position = dodge, width = 0.25)  +  scale_fill_manual(values=c("grey58","grey58", "grey58","grey58"),name=element_blank()) +
theme(legend.position="none")+
ylab("Gini of rival wealth")  + theme(axis.text=element_text(size=14),
axis.title=element_text(size=14))+  theme(axis.title.x=element_blank())
ggsave("GiniType.pdf", p, height=4, width=6.5)

print("made GiniType figure")


################################################################################
# Fig 3
################################################################################

fracs<-2

GiniWealth <-  PercentWithCowives <- PercentRich <- WealthRatio <- matrix(999999, nrow=2000,ncol=P)

PercentWithCowives <- matrix(999999, nrow=2000,ncol=P)

for(i in 1:P){
  Scrap <- Results_Rival(AdjustedRival[[i]],fracs)
  GiniWealth[1:dim(Scrap)[1],i] <- Scrap[,1]
  PercentRich[1:dim(Scrap)[1],i] <- Scrap[,2]
  WealthRatio[1:dim(Scrap)[1],i] <- Scrap[,3]
}

GR <- GiniWealth
PR <- PercentRich
WR <- WealthRatio

for(i in 1:P){
  Scrap <- Results_Wives(AdjustedWives[[i]])
  PercentWithCowives[1:dim(Scrap)[1],i] <- Scrap[,1]
}

PFP <- PercentWithCowives

Type <- c(1,1,1,1,1,3,1,1,3,3,3,2,4,3,3,4,3,3,2,1,3,4,2,2,2,2,2,2,4,4,4,4,4,4)[1:P]
Subsist <- c("Agropastoral","Forager","Horticultural","Agricultural")[Type]
Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )

X <- 2000

df<-data.frame(Gini=colMeans(GR), Polygyny=colMeans(PFP),Subsistence=Subsist)
p <- ggplot(df, aes(Gini, Polygyny)) + labs(y="Completed percent female polygyny") + labs(x="Gini of completed rival wealth")+   theme(legend.title=element_text(size=14),legend.text=element_text(size=14))+
theme(strip.text.x = element_text(size=14),axis.text=element_text(size=14), axis.title=element_text(size=14)) +
geom_point(aes(colour = Subsistence),size=3) + scale_shape_manual(values=c(19,1)) + scale_colour_manual(values=cbbPalette[c(5,1,3,2)])
p

ggsave("PPP-Gini.pdf",p,width=8,height=6)

print("made PPP-Gini figure")

###################
fitlm<-lm(df$Polygyny[which(df$Subsistence=="Forager")]~df$Gini[df$Subsistence=="Forager"])
print(summary(fitlm))
print(confint(fitlm, level=0.9))

fitlm<-lm(df$Polygyny[which(df$Subsistence=="Horticultural")]~df$Gini[df$Subsistence=="Horticultural"])
print(summary(fitlm))
print(confint(fitlm, level=0.9))

fitlm<-lm(df$Polygyny[which(df$Subsistence=="Agropastoral")]~df$Gini[df$Subsistence=="Agropastoral"])
print(summary(fitlm))
print(confint(fitlm, level=0.9))

fitlm<-lm(df$Polygyny[which(df$Subsistence=="Agricultural")]~df$Gini[df$Subsistence=="Agricultural"])
print(summary(fitlm))
print(confint(fitlm, level=0.9))

fitlm<-lm(df$Polygyny[which(df$Subsistence=="Agropastoral" | df$Subsistence=="Horticultural")]~df$Gini[which(df$Subsistence=="Agropastoral" | df$Subsistence=="Horticultural")])
print(summary(fitlm))
print(confint(fitlm, level=0.9))
###################

################################################################################
# Fig 5 d
################################################################################

WealthRatio2 <- GiniWealth2 <- PercentRich2  <- matrix(NA, nrow=2000,ncol=P)

for(i in c(1:29)[c(-13,-14)]){         # Drop cases with no variation in wives
  Scrap <- Results_Rival2(AdjustedRival[[i]],AdjustedWives[[i]])
  GiniWealth2[1:dim(Scrap)[1],i] <- Scrap[,1]
  PercentRich2[1:dim(Scrap)[1],i] <- Scrap[,2]
  WealthRatio2[1:dim(Scrap)[1],i] <- Scrap[,3]
}

PR <- PercentRich2
d <- read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")
Subsist <- d$Subsistence
Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )
Subsist <- as.character(Subsist)

# Check t test for each frame
qq1 <- t.test( colMeans(PR,na.rm=TRUE)[which(Subsist=="Forager")], colMeans(PR,na.rm=TRUE)[which(Subsist=="Agricultural")],conf.level = 0.9) #
qq2 <- t.test( colMeans(PR,na.rm=TRUE)[which(Subsist=="Horticultural")], colMeans(PR,na.rm=TRUE)[which(Subsist=="Agricultural")],conf.level = 0.9)
qq3 <- t.test( colMeans(PR,na.rm=TRUE)[which(Subsist=="Agropastoral")], colMeans(PR,na.rm=TRUE)[which(Subsist=="Agricultural")],conf.level = 0.9)

Subsist <- ifelse(Subsist=="Forager",paste0("Forager: ",paste0(round(qq1$estimate[1]-qq1$estimate[2],2)," (", round(qq1$conf.int[1],2),", ",round(qq1$conf.int[2],2),")" )),Subsist)
Subsist <- ifelse(Subsist=="Horticultural",paste0("Horticultural: ",paste0(round(qq2$estimate[1]-qq2$estimate[2],2)," (", round(qq2$conf.int[1],2),", ",round(qq2$conf.int[2],2),")" )),Subsist)
Subsist <- ifelse(Subsist=="Agropastoral",paste0("Agropastoral: ",paste0(round(qq3$estimate[1]-qq3$estimate[2],2)," (", round(qq3$conf.int[1],2),", ",round(qq3$conf.int[2],2),")" )),Subsist)
Subsist <- ifelse(Subsist=="Agricultural","Agricultural",Subsist)

Subsist <- factor(Subsist, levels= c(paste0("Forager: ",paste0(round(qq1$estimate[1]-qq1$estimate[2],2)," (", round(qq1$conf.int[1],2),", ",round(qq1$conf.int[2],2),")" )),
                                paste0("Horticultural: ",paste0(round(qq2$estimate[1]-qq2$estimate[2],2)," (", round(qq2$conf.int[1],2),", ",round(qq2$conf.int[2],2),")" )),
                                paste0("Agropastoral: ",paste0(round(qq3$estimate[1]-qq3$estimate[2],2)," (", round(qq3$conf.int[1],2),", ",round(qq3$conf.int[2],2),")" )),
                                "Agricultural")   )

df<-data.frame(PercentRich=colMeans(PercentRich2,na.rm=TRUE), Subsistence=Subsist)
p <- ggplot(df, aes(Subsistence, PercentRich)) + labs(y="Frequency of rich males (> threshold) ") + theme(legend.title=element_text(size=14),legend.text=element_text(size=14)) +
theme(strip.text.x = element_text(size=14),axis.text=element_text(size=12), axis.title=element_text(size=14)) +
geom_point(aes(colour = Subsistence),size=3) + scale_colour_manual("Mean difference between focal\n subsistence type and agriculture:",values=cbbPalette[c(5,1,3,2)]) + theme( axis.text.x=element_blank(),
axis.title.x=element_blank(),
axis.ticks.x=element_blank())                                       
p

ggsave("PPP-PR-Emp.pdf",p,width=5,height=4)

print("made PPP-PR-Emp figure")

################################################################################
# Step 5: Calculate Estimates of Interest for Table 2
################################################################################

rm(list=ls())

load("Extrapolations.RData")

source('./Code/project_support.R')

Names <- c(
"Kipsigis",
"Sangu",
"Maasai",
"Maasai2",
"Koore",
"Chagga",
"Sidama",
"Sukuma",
"Pimbwe",
"Maya2",
"Tsimane",
"Dolgan",
"Polish",
"Maya",
"Matsigenka",
"Bangladesh",
"Tsimane2",
"Makushi",
"Lamalera",
"Himba",
"Mayangna",
"English",
"Ache",
"Agta",
"Aka",
"Hadza",
"Kung",
"Pume",
"Krummhorn"
)

fracs<-2
P <-29

GiniWealth <-  PercentWithCowives <- PercentRich <- WealthRatio <- matrix(999999, nrow=2000,ncol=P)

for(i in 1:P){
  Scrap <- Results_Rival(AdjustedRival[[i]],fracs)
  GiniWealth[1:dim(Scrap)[1],i] <- Scrap[,1]
  PercentRich[1:dim(Scrap)[1],i] <- Scrap[,2]
  WealthRatio[1:dim(Scrap)[1],i] <- Scrap[,3]
}

GR <- GiniWealth
PR <- PercentRich
WR <- WealthRatio

UPWC <- c()

for(i in 1:P){
  Scrap <- Results_Wives(AdjustedWives[[i]])
  PercentWithCowives[1:dim(Scrap)[1],i] <- Scrap[,1]
  
  UPWC[i] <- sum(UnAdjustedWives[[i]][1,][which(AgeVec[[i]][1,]>44 & UnAdjustedWives[[i]][1,]>1)])/
             sum(UnAdjustedWives[[i]][1,][which(AgeVec[[i]][1,]>44)])
}

PFP <- PercentWithCowives

d <- read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")
Subsist <- d$Subsistence
Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )

######################################################################## Table 2

RES2 <- matrix( NA, nrow=P, ncol=7)
for(i in 1:P){
  RES2[i,1] <-  Names[i]
  RES2[i,2] <-  as.character(Subsist[i])
  RES2[i,3] <-  paste0(prettyM(GR[,i]),  ", (", prettyCI(GR[,i],1),  ", ", prettyCI(GR[,i],2),  ")")
  RES2[i,4] <-  paste0(prettyM2(WR[,i]),  ", (", prettyCI2(WR[,i],1),  ", ", prettyCI2(WR[,i],2),  ")")
  RES2[i,5] <-  paste0(prettyM(PR[,i]),  ", (", prettyCI(PR[,i],1),  ", ", prettyCI(PR[,i],2),  ")")
  RES2[i,6] <-  paste0(prettyM(PFP[,i]), ", (", prettyCI(PFP[,i],1), ", ", prettyCI(PFP[,i],2), ")")
  RES2[i,7] <-  round(UPWC[i],2)
}

RES2<-RES2[order(Subsist),]
RESe <- matrix( NA, nrow=4, ncol=7)

RESe[1,2] <- "Forager"
RESe[1,3] <-  paste0(prettyM(RMC(GR,RESe[1,2])),  ", (", prettyCI(RMC(GR,RESe[1,2]),1)  ,", ", prettyCI(RMC(GR,RESe[1,2]),2),  ")")
RESe[1,4] <-  paste0(prettyM2(RMC(WR,RESe[1,2])),  ", (", prettyCI2(RMC(WR,RESe[1,2]),1)  ,", ", prettyCI2(RMC(WR,RESe[1,2]),2),  ")")
RESe[1,5] <-  paste0(prettyM(RMC(PR,RESe[1,2])),  ", (", prettyCI(RMC(PR,RESe[1,2]),1)  ,", ", prettyCI(RMC(PR,RESe[1,2]),2),  ")")
RESe[1,6] <-  paste0(prettyM(RMC(PFP,RESe[1,2])), ", (", prettyCI(RMC(PFP,RESe[1,2]),1) ,", ", prettyCI(RMC(PFP,RESe[1,2]),2), ")")
RESe[1,7] <-  paste0(prettyM(UPWC[which(Subsist=="Forager")]), ", (", prettyCI(UPWC[which(Subsist=="Forager")],1) ,", ", 
               prettyCI(UPWC[which(Subsist=="Forager")],2), ")")

RESe[2,2] <- "Horticultural"
RESe[2,3] <-  paste0(prettyM(RMC(GR,RESe[2,2])),  ", (", prettyCI(RMC(GR,RESe[2,2]),1),  ", ", prettyCI(RMC(GR,RESe[2,2]),2)  ,")")
RESe[2,4] <-  paste0(prettyM2(RMC(WR,RESe[2,2])),  ", (", prettyCI2(RMC(WR,RESe[2,2]),1),  ", ", prettyCI2(RMC(WR,RESe[2,2]),2)  ,")")
RESe[2,5] <-  paste0(prettyM(RMC(PR,RESe[2,2])),  ", (", prettyCI(RMC(PR,RESe[2,2]),1),  ", ", prettyCI(RMC(PR,RESe[2,2]),2)  ,")")
RESe[2,6] <-  paste0(prettyM(RMC(PFP,RESe[2,2])), ", (", prettyCI(RMC(PFP,RESe[2,2]),1), ", ", prettyCI(RMC(PFP,RESe[2,2]),2) ,")")
RESe[2,7] <-  paste0(prettyM(UPWC[which(Subsist=="Horticultural")]), ", (", prettyCI(UPWC[which(Subsist=="Horticultural")],1) ,", ", 
               prettyCI(UPWC[which(Subsist=="Horticultural")],2), ")")

RESe[3,2] <- "Agropastoral"
RESe[3,3] <-  paste0(prettyM(RMC(GR,RESe[3,2])),  ", (", prettyCI(RMC(GR,RESe[3,2]),1),  ", ", prettyCI(RMC(GR,RESe[3,2]),2)  ,")")
RESe[3,4] <-  paste0(prettyM2(RMC(WR,RESe[3,2])),  ", (", prettyCI2(RMC(WR,RESe[3,2]),1),  ", ", prettyCI2(RMC(WR,RESe[3,2]),2)  ,")")
RESe[3,5] <-  paste0(prettyM(RMC(PR,RESe[3,2])),  ", (", prettyCI(RMC(PR,RESe[3,2]),1),  ", ", prettyCI(RMC(PR,RESe[3,2]),2)  ,")")
RESe[3,6] <-  paste0(prettyM(RMC(PFP,RESe[3,2])), ", (", prettyCI(RMC(PFP,RESe[3,2]),1), ", ", prettyCI(RMC(PFP,RESe[3,2]),2) ,")")
RESe[3,7] <-  paste0(prettyM(UPWC[which(Subsist=="Agropastoral")]), ", (", prettyCI(UPWC[which(Subsist=="Agropastoral")],1) ,", ", 
               prettyCI(UPWC[which(Subsist=="Agropastoral")],2), ")")

RESe[4,2] <- "Agricultural"
RESe[4,3] <-  paste0(prettyM(RMC(GR,RESe[4,2])),  ", (", prettyCI(RMC(GR,RESe[4,2]),1),  ", ", prettyCI(RMC(GR,RESe[4,2]),2),  ")")
RESe[4,4] <-  paste0(prettyM2(RMC(WR,RESe[4,2])),  ", (", prettyCI2(RMC(WR,RESe[4,2]),1),  ", ", prettyCI2(RMC(WR,RESe[4,2]),2),  ")")
RESe[4,5] <-  paste0(prettyM(RMC(PR,RESe[4,2])),  ", (", prettyCI(RMC(PR,RESe[4,2]),1),  ", ", prettyCI(RMC(PR,RESe[4,2]),2),  ")")
RESe[4,6] <-  paste0(prettyM(RMC(PFP,RESe[4,2])), ", (", prettyCI(RMC(PFP,RESe[4,2]),1), ", ", prettyCI(RMC(PFP,RESe[4,2]),2), ")")
RESe[4,7] <-  paste0(prettyM(UPWC[which(Subsist=="Agricultural")]), ", (", prettyCI(UPWC[which(Subsist=="Agricultural")],1) ,", ", 
               prettyCI(UPWC[which(Subsist=="Agricultural")],2), ")")


################# Re-order
RES2[which(RES2[,2]=="Forager"),] <- RES2[which(RES2[,2]=="Forager"),][order(RES2[which(RES2[,2]=="Forager"),1]),]
RES2[which(RES2[,2]=="Horticultural"),] <- RES2[which(RES2[,2]=="Horticultural"),][order(RES2[which(RES2[,2]=="Horticultural"),1]),]
RES2[which(RES2[,2]=="Agropastoral"),] <- RES2[which(RES2[,2]=="Agropastoral"),][order(RES2[which(RES2[,2]=="Agropastoral"),1]),]
RES2[which(RES2[,2]=="Agricultural"),] <- RES2[which(RES2[,2]=="Agricultural"),][order(RES2[which(RES2[,2]=="Agricultural"),1]),]

RES<- rbind(RES2,RESe)

x <- print(xtable(RES), include.rownames=TRUE)

writeLines(x, "table2.txt")


print("made Table 2")
