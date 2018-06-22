###############################################################################
# Fig 6 a and b

rm(list=ls())

load("ModelResults.RData")
all_stuff <- ls()
rm( list=setdiff(all_stuff, "fitPoly") )
source('./Code/project_support.R')


d <- read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")
Subsist <- d$Subsistence
Subsist <- factor(Subsist, levels= c("Forager","Horticultural","Agropastoral","Agricultural")   )

Pops<-as.character(d$Label)
Pops[which(Pops=="!Kung")]<- "Kung"

################################################################### Make fig a
mu <- rstan::extract(fitPoly,pars="Beta_RS")$Beta_RS[,,2]
m1 <- (mu)
sample_eff <- apply(m1,2,quantile,probs=c(0.05,0.5,0.95))
df_all1 <- data.frame(Order=c(1:length(Subsist)), Population=Pops, Subsistence=Subsist, LI=sample_eff[1,],
                     Median=sample_eff[2,], HI=sample_eff[3,])


df_all1 <- df_all1[order(df_all1$Subsistence, df_all1$Population),]
df_all1$x <- 1:length(Subsist)

beta <- rstan::extract(fitPoly,pars="Beta_RS")$Beta_RS[,,3]
sample_eff <- apply(beta,2,quantile,probs=c(0.05,0.5,0.95))

sample_eff[,13]  <- c(NA,NA,NA) # These arent actualy estimates, just the prior, since wife number is fixed
sample_eff[,14]  <- c(NA,NA,NA) #

df_all2 <- data.frame(Order=c(1:length(Subsist)), Population=Pops,Subsistence=Subsist, LI=sample_eff[1,],
                     Median=sample_eff[2,], HI=sample_eff[3,])
df_all2 <- df_all2[order(df_all2$Subsistence, df_all2$Population),]
df_all2$x <- 1:length(Subsist) + 0.22

df_all <- rbind(df_all1,df_all2)
df_all$Estimate <- c(rep("\u03BC",length(Subsist)),rep("\u03B4 - \u03BC",length(Subsist)))

ssss <- (1.3)

p <- ggplot(df_all,aes(x=x,y=Median,colour=Estimate,shape=Estimate,size=Estimate,fill=Estimate))+geom_point(size=2)+
    geom_linerange(aes(ymin=LI,ymax=HI,colour=Estimate),size=1) +
    theme(legend.title=element_text(size=14),legend.text=element_text(size=14)) +
    geom_hline(aes(yintercept=1),linetype="dashed")+
    geom_hline(aes(yintercept=0),linetype="dashed")+
    geom_vline(aes(xintercept=8.75),size=1,linetype="dashed")+
    geom_vline(aes(xintercept=17.75),size=1,linetype="dashed")+
    geom_vline(aes(xintercept=25.75),size=1,linetype="dashed")+
    labs(y="Elasticity") + theme(strip.text.x = element_text(size=14,face="bold"),axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"))  +
    scale_x_discrete(name ="Population", 
                    limits=as.character(c(1:29))) +
    annotate("text", x = ((0+9)/2), y = ssss+ -0.08, label = "Foraging", size=5) +
    annotate("text", x = ((10+17)/2), y = ssss+ -0.08, label = "Horticultural", size=5) +
    annotate("text", x = ((18+25)/2), y = ssss+ -0.08, label = "Agropastoral", size=5) +
    annotate("text", x = ((24+32)/2), y = ssss+ -0.08, label = "Agricultural", size=5) + 
    expand_limits(x = 30) + scale_colour_manual(values = c("firebrick", "darkslateblue"))+
    scale_fill_manual(values = c("firebrick", "darkslateblue"))+ scale_shape_manual(values=c(21, 23)) 

ggsave("PPP-Elast1.pdf",p,width=8.5,height=3, device=cairo_pdf)

print("made PPP-Elast1 figure")

################################################################### Make fig b
mu <- rstan::extract(fitPoly,pars="Beta_RS")$Beta_RS[,,2]

sample_eff <- apply(mu,2,quantile,probs=c(0.05,0.5,0.95))
df_all1 <- data.frame(Order=c(1:length(Subsist)),Population=Pops, Subsistence=Subsist, LI=sample_eff[1,],
                     Median=sample_eff[2,], HI=sample_eff[3,])
df_all1 <- df_all1[order(df_all1$Subsistence, df_all1$Population),]
df_all1$x <- 1:length(Subsist)

mu <- rstan::extract(fitPoly,pars="Beta_RS")$Beta_RS[,,2]
beta <- rstan::extract(fitPoly,pars="Beta_RS")$Beta_RS[,,3]
mpb <-  beta + mu

sample_eff <- apply(mpb,2,quantile,probs=c(0.05,0.5,0.95))

sample_eff[,13]  <-c(NA,NA,NA) # These arent actualy estimates, just the prior, since wife number is fixed
sample_eff[,14]  <-c(NA,NA,NA) #

df_all2 <- data.frame(Order=c(1:length(Subsist)),Population=Pops, Subsistence=Subsist, LI=sample_eff[1,],
                     Median=sample_eff[2,], HI=sample_eff[3,])
df_all2 <- df_all2[order(df_all2$Subsistence, df_all2$Population),]
df_all2$x <- 1:length(Subsist) +0.22

df_all <- rbind(df_all2)
df_all$Estimate <- c(rep("\u03B4",29))

ssss<- 1.45

p <- ggplot(df_all,aes(x=x,y=Median,colour=Estimate,shape=Estimate,size=Estimate,fill=Estimate))+geom_point(size=2)+
   geom_linerange(aes(ymin=LI,ymax=HI,colour=Estimate),size=1) +
   theme(legend.title=element_text(size=14),legend.text=element_text(size=14)) +
   geom_hline(aes(yintercept=1),linetype="dashed")+
   geom_hline(aes(yintercept=0),linetype="dashed")+
   geom_vline(aes(xintercept=8.75),size=1,linetype="dashed")+
   geom_vline(aes(xintercept=17.75),size=1,linetype="dashed")+
   geom_vline(aes(xintercept=25.75),size=1,linetype="dashed")+
   labs(y="Elasticity") + theme(strip.text.x = element_text(size=14,face="bold"),axis.text=element_text(size=12),
   axis.title=element_text(size=14,face="bold"))+ scale_x_discrete(name ="Population", 
   limits=as.character(c(1:29))) +
   annotate("text", x = ((0+9)/2), y = ssss+ -0.08, label = "Foraging", size=5) +
   annotate("text", x = ((10+17)/2), y = ssss+ -0.08, label = "Horticultural", size=5) +
   annotate("text", x = ((18+25)/2), y = ssss+ -0.08, label = "Agropastoral", size=5) +
   annotate("text", x = ((24+32)/2), y = ssss+ -0.08, label = "Agricultural", size=5)+ 
   expand_limits(x = 30)  + scale_colour_manual(values = c("springgreen4"))+
   scale_fill_manual(values = c("springgreen4"))+ scale_shape_manual(values=c(21, 23)) 

ggsave("PPP-Elast2.pdf",p,width=8.5,height=3, device=cairo_pdf)

print("made PPP-Elast2 figure")


###############
QQ<-extract(fitPoly,pars="MU_RS")$MU_RS
print("Mean mu:")
print(paste0(prettyM(QQ[,2]),  ", (", prettyCI(QQ[,2],1),  ", ", prettyCI(QQ[,2],2)  ,")")   )

