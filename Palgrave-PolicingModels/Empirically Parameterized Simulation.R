################################################################################
# Resolution of Apparent Paradoxes in the Racially Biased  
# Use of Lethal and Non-Lethal Force by Police
#
# Empirically Parameterized Model

###################################################### Libraries
  library(bbmle)
  library(rethinking)
  library(Cairo)
  library(ggplot2)

# Set Control Parameters
  Iters <- 33       # Maximum percent of encounters arising through Escalated Non-Lethal force model
  Samples <- 500    # Times each model is sampled to build posterior

# Open Array Storage
Results <- array(NA,c(9,Iters,Iters,Samples))

# Begin Generative Model
for(j in 1:Samples){        # Sample each model to build up a posterior
for(k in 0:(Iters-1)){      # Loop through possible level of Escalated Non-Lethal force in black subpopulation
for(q in 0:(Iters-1)){      # Loop through possible level of Escalated Non-Lethal force in white subpopulation

################################################################################
######################################################### Set Parameters

###################################################### Population Data
# Make Population Data
 N_W <- 198000000
 N_B <-  37000000
 N   <- N_W + N_B

###################################################### Standard Policing Model
# Make Encouter Rates
 e_W <- 0.0503-(0.0503*q*0.01)
 e_B <- 0.1404-(0.1404*k*0.01)

# Make Non-Lethal Force Rates
 g_W <- 0.0015
 g_B <- 0.0015

# Make Lethal Force Rates
 s_W <- 0.000065
 s_B <- 0.000065

###################################################### Escalated Non-Lethal Force Model
# Make Encouter Rates
 e_hat_W <-(0.0503*q*0.01)
 e_hat_B <- (0.1404*k*0.01)

# Make Non-Lethal Force Rates
 g_hat_W <- 0.015
 g_hat_B <- 0.015

# Make Lethal Force Rates
 s_hat_W <- 0.0000065
 s_hat_B <- 0.0000065


################################################################################
######################################################### Generate Data
###################################### Standard Policing Model
# Make Encouters
 E_W <- rbinom(1,N_W,e_W)
 E_B <- rbinom(1,N_B,e_B)

# Make Violence Outcomes
 S_W <- rmultinom(1,E_W,c(s_W, g_W, (1-(g_W+s_W))))
 S_B <- rmultinom(1,E_B,c(s_B, g_B, (1-(g_B+s_B))))

###################################### Escalated Non-Lethal Force Model
# Make Encouters
 E_hat_W <- rbinom(1,N_W,e_hat_W)
 E_hat_B <- rbinom(1,N_B,e_hat_B)

# Make Violence Outcomes
 S_hat_W <- rmultinom(1,E_hat_W,c(s_hat_W, g_hat_W, (1-(g_hat_W+s_hat_W))))
 S_hat_B <- rmultinom(1,E_hat_B,c(s_hat_B, g_hat_B, (1-(g_hat_B+s_hat_B))))


################################################################################
######################################################### Sum Data
################################################################################
# Make Total Encounters
 E_sum_W <- E_W + E_hat_W
 E_sum_B <- E_B + E_hat_B

# Make Total Outcomes
 S_sum_W <- S_W + S_hat_W
 S_sum_B <- S_B + S_hat_B

################################################################################
######################################################### Get Relative Risks
################################################################################
# Make Total Encounter Rate
 RE_W <- E_sum_W/N_W
 RE_B <- E_sum_B/N_B

 Phi <- RE_B/RE_W
 Phi_hat  <- e_hat_B/e_hat_W

# Make Total Outcome Rates
 RS_W <- S_sum_W[1]/N_W
 RS_B <- S_sum_B[1]/N_B

 Psi <- RS_B/RS_W

 RG_W <- S_sum_W[2]/N_W
 RG_B <- S_sum_B[2]/N_B

 Delta <- RG_B/RG_W

# Make Conditional Outcome Rates
 CS_W <- S_sum_W[1]/E_sum_W
 CS_B <- S_sum_B[1]/E_sum_B

 Theta <- CS_B/CS_W

 CG_W <- S_sum_W[2]/E_sum_W
 CG_B <- S_sum_B[2]/E_sum_B

 Gamma <- CG_B/CG_W

################# Analize Simulated Data
 E <- c(E_sum_W, E_sum_B)
 S <- c(S_sum_W[1], S_sum_B[1])
 G <- c(S_sum_W[2], S_sum_B[2])
 B <- c(0,1)

 M_L  <- mle2(S~dbinom(logistic(a + b*B),size=E), start=list(a=-5,b=0),data=data.frame(S=S,B=B,E=E))
 M_NL <- mle2(G~dbinom(logistic(a + b*B),size=E), start=list(a=-3,b=0),data=data.frame(G=G,B=B,E=E))
 M_R  <- mle2(S~dbinom(logistic(a + b*B),size=(G+S)), start=list(a=-3,b=0),data=data.frame(G=G,B=B,E=E,S=S))

 S_M_L <- summary(M_L)
 S_M_NL <- summary(M_NL)
 S_M_R <- summary(M_R)

 Beta_L <- rnorm(1,S_M_L@coef[2,1], S_M_L@coef[2,2])
 Beta_NL <- rnorm(1,S_M_NL@coef[2,1], S_M_NL@coef[2,2])
 Beta_R <- rnorm(1,S_M_R@coef[2,1], S_M_R@coef[2,2])

Results[,k+1,q+1,j]<-c(Phi, Psi, Delta, Theta, Gamma, Beta_L, Beta_NL, Beta_R, Phi_hat)
  }}}

################################################################################
########################################################### Process Model Output
SumRes <- array(NA,c(Iters,Iters,9,3))

for(k in 1:Iters){      # Calculate Summary Stats
for(q in 1:Iters){
for(j in 1:9){
 SumRes[k,q,j,] <- c(median(c(na.omit(c(Results[j,k,q,])))),PCI(c(na.omit(c(Results[j,k,q,]))),0.99)[1],PCI(c(na.omit(c(Results[j,k,q,]))),0.99)[2])
      }}}

 RES<-as.matrix(SumRes[,,3,1])
 colnames(RES)<-seq(0,0.32, length.out=33)

## color brewing
library("RColorBrewer"); library("lattice");
brewer.div <- colorRampPalette(brewer.pal(11, "Spectral"), interpolate = "spline")

################################################################### Plot PsiU
 RES<-as.matrix(SumRes[,,2,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c(log(RES)))

 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso", main=expression(paste("log(", Psi[u],") - Population-Level Lethal Force Ratio")),
  xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

################################################################### Plot DeltaU
 RES<-as.matrix(SumRes[,,3,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c(log(RES)))
 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso",main=expression(paste("log(", Delta[u],") - Population-Level Non-Lethal Force Ratio")),
  xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

################################################################### Plot ThetaU
 RES<-as.matrix(SumRes[,,4,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c(log(RES)))
 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso",main=expression(paste("log(", Theta[u],") - Encounter-Conditional Lethal Force Ratio")),
  xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

################################################################### Plot GammaU
 RES<-as.matrix(SumRes[,,5,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c(log(RES)))
 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso",main=expression(paste("log(", Gamma[u],") - Encounter-Conditional Non-Lethal Force Ratio")),
  xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

################################################################### Plot BetaL
 RES<-as.matrix(SumRes[,,6,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c((RES)))
 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso",main=expression(paste(beta[L]," - Lethal Force From Encounters")),
 xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

################################################################### Plot BetaN
 RES<-as.matrix(SumRes[,,7,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c((RES)))
 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso",main=expression(paste(beta[N]," - Non-Lethal Force From Encounters")),
  xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

################################################################### Plot BetaR
 RES<-as.matrix(SumRes[,,8,1])
 d = data.frame(x=rep(seq(0, 0.32, length=nrow(RES)), ncol(RES)),
               y=rep(seq(0, 0.32, length=ncol(RES)), each=nrow(RES)),
               z=c((RES)))
 levelplot (z~x*y,data=d, cuts = 199, col.regions = brewer.div(200), aspect = "iso",main=expression(paste(beta[R]," - Lethal Force From Any Force")),
   xlab=expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model")),
   ylab=expression(paste(epsilon[W] - "Fraction of White Encounters from Escalated Non-Lethal Force Model")),
  panel = function(...) {
              panel.levelplot(...)
              panel.abline(c(0,1))
                })

############################################################## Plot CIs with Slices
 SumRes2<-log(SumRes)
        q<-5
 dat <- data.frame(M=c(SumRes2[,q,2,1],SumRes2[,q,3,1],SumRes2[,q,4,1],SumRes2[,q,5,1]),
                   L=c(SumRes2[,q,2,2],SumRes2[,q,3,2],SumRes2[,q,4,2],SumRes2[,q,5,2]),
                   H=c(SumRes2[,q,2,3],SumRes2[,q,3,3],SumRes2[,q,4,3],SumRes2[,q,5,3]),
                   Metric=c(rep("Psi",33),rep("Delta",33),rep("Theta",33),rep("Gamma",33)),
                   Ehat = c(seq(0,0.32,length.out=33),seq(0,0.32,length.out=33),seq(0,0.32,length.out=33),seq(0,0.32,length.out=33)))

 dat$Metric <- factor(dat$Metric, levels = c("Psi","Delta","Theta","Gamma" ))

 Pp <- ggplot(dat, aes(x=Ehat, y=M)) +
 geom_line(aes(colour=Metric),size=1) +
 geom_ribbon(aes(ymax=H, ymin=L, fill=Metric), alpha = 0.7) +
 geom_hline(aes(yintercept=0))+
 scale_fill_manual(values = c("#1f78b4","#a6cee3","#33a02c","#b2df8a"),labels=c(expression(Psi[ u]), expression(Delta[   u]), expression(Theta[ u]),expression(Gamma[ u]))) +
 scale_colour_manual(values = c("#1f78b4","#a6cee3","#33a02c","#b2df8a"),labels=c(expression(Psi[ u]), expression(Delta[  u]), expression(Theta[ u]),expression(Gamma[ u]))) +
 xlab(expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model"))) +
  ylab("log(Estimated Value)")+
   theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
 legend.text=element_text(size=18), legend.title=element_text(size=20,face="bold"))
       Pp
ggsave("DataEmpirical.pdf",height=8,width=10, device=cairo_pdf)


 dat <- data.frame(M=c(SumRes[,q,6,1],SumRes[,q,7,1],SumRes[,q,8,1]) ,
                   L=c(SumRes[,q,6,2],SumRes[,q,7,2],SumRes[,q,8,2]),
                   H=c(SumRes[,q,6,3],SumRes[,q,7,3],SumRes[,q,8,3]),
                   Metric=c(rep("Beta_L",33),rep("Beta_N",33),rep("Beta_R",33)),
                   Ehat = c(seq(0,0.32,length.out=33),seq(0,0.32,length.out=33),seq(0,0.32,length.out=33)))

 dat$Metric <- factor(dat$Metric, levels = c("Beta_L","Beta_N","Beta_R" ))

 Pp <- ggplot(dat, aes(x=Ehat, y=M)) +
 geom_line(aes(colour=Metric),size=1) +
 geom_ribbon(aes(ymax=H, ymin=L, fill=Metric), alpha = 0.7) +
 geom_hline(aes(yintercept=0))+
 scale_fill_manual(values = c("#e31a1c","#fb9a99","#cab2d6"),labels=c(expression(beta[L]),expression(beta[N]),expression(beta[R]))) +
 scale_colour_manual(values = c("#e31a1c","#fb9a99","#cab2d6"),labels=c(expression(beta[L]),expression(beta[N]),expression(beta[R]))) +
 xlab(expression(paste(epsilon[B] - "Fraction of Black Encounters from Escalated Non-Lethal Force Model"))) +
  ylab("Estimated Value")+
   theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
 legend.text=element_text(size=18), legend.title=element_text(size=20,face="bold"))
 Pp
ggsave("ParametersEmpricial.pdf",height=8,width=10, device=cairo_pdf)
