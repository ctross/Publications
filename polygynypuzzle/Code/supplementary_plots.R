library(ineq)

P <- function(s,t,m,u,c){
       X <- s*t*m*(d-u)*(1/(d*c))
       return(ifelse(X>1,NA,X))
       }

G <- function(t,m){
       X <- ((t*m )/(t*m + (1-t)))- t
       return(X)
       }

############################################  Plot 1
s <- 1
m <- 1.5
u <- 0.10
c <- 0.95
d <- 0.15

Size2<-1.75

t <- seq(0.001,0.999,length.out=1000)
pdf("SimulationResults-d15-u10.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(mar=c(5,5,5,5))
Size<-4
plot(P(s,t,m+1,u,c)~I(1-t),ylim=c(0,1),ylab="Percent female polygyny, P",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(P(s,t,m+2,u,c)~I(1-t),col="firebrick",lwd=Size)
lines(P(s,t,m+4,u,c)~I(1-t),col="peru",lwd=Size)
lines(P(s,t,m+8,u,c)~I(1-t),col="seagreen4",lwd=Size)
lines(P(s,t,m+16,u,c)~I(1-t),col="darkslateblue",lwd=Size)
lines(P(s,t,m+32,u,c)~I(1-t),col="mediumorchid4",lwd=Size)

plot(G(t,m+1)~I(1-t),ylim=c(0,1),ylab="Gini of wealth, G",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(G(t,m+2)~I(1-t),col="firebrick",lwd=Size)
lines(G(t,m+4)~I(1-t),col="peru",lwd=Size)
lines(G(t,m+8)~I(1-t),col="seagreen4",lwd=Size)
lines(G(t,m+16)~I(1-t),col="darkslateblue",lwd=Size)
lines(G(t,m+32)~I(1-t),col="mediumorchid4",lwd=Size)

plot(I(1/(G(t,m+1)/P(s,t,m+1,u,c)))~I(1-t),ylim=c(0,2),ylab="Ratio of P to G",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(I(1/(G(t,m+2)/P(s,t,m+2,u,c)))~I(1-t),col="firebrick",lwd=Size)
lines(I(1/(G(t,m+4)/P(s,t,m+4,u,c)))~I(1-t),col="peru",lwd=Size)
lines(I(1/(G(t,m+8)/P(s,t,m+8,u,c)))~I(1-t),col="seagreen4",lwd=Size)
lines(I(1/(G(t,m+16)/P(s,t,m+16,u,c)))~I(1-t),col="darkslateblue",lwd=Size)
lines(I(1/(G(t,m+32)/P(s,t,m+32,u,c)))~I(1-t),col="mediumorchid4",lwd=Size)
  dev.off()
  


############################################  Plot 2
s <- 1
m <- 1.5
u <- 0.10
c <- 0.95
d <- 0.50

t <- seq(0.001,0.999,length.out=1000)
pdf("SimulationResults-d50-u10.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(mar=c(5,5,5,5))
Size<-4
plot(P(s,t,m+1,u,c)~I(1-t),ylim=c(0,1),ylab="Percent female polygyny, P",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(P(s,t,m+2,u,c)~I(1-t),col="firebrick",lwd=Size)
lines(P(s,t,m+4,u,c)~I(1-t),col="peru",lwd=Size)
lines(P(s,t,m+8,u,c)~I(1-t),col="seagreen4",lwd=Size)
lines(P(s,t,m+16,u,c)~I(1-t),col="darkslateblue",lwd=Size)
lines(P(s,t,m+32,u,c)~I(1-t),col="mediumorchid4",lwd=Size)

plot(G(t,m+1)~I(1-t),ylim=c(0,1),ylab="Gini of wealth, G",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(G(t,m+2)~I(1-t),col="firebrick",lwd=Size)
lines(G(t,m+4)~I(1-t),col="peru",lwd=Size)
lines(G(t,m+8)~I(1-t),col="seagreen4",lwd=Size)
lines(G(t,m+16)~I(1-t),col="darkslateblue",lwd=Size)
lines(G(t,m+32)~I(1-t),col="mediumorchid4",lwd=Size)

plot(I(1/(G(t,m+1)/P(s,t,m+1,u,c)))~I(1-t),ylim=c(0,2),ylab="Ratio of P to G",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(I(1/(G(t,m+2)/P(s,t,m+2,u,c)))~I(1-t),col="firebrick",lwd=Size)
lines(I(1/(G(t,m+4)/P(s,t,m+4,u,c)))~I(1-t),col="peru",lwd=Size)
lines(I(1/(G(t,m+8)/P(s,t,m+8,u,c)))~I(1-t),col="seagreen4",lwd=Size)
lines(I(1/(G(t,m+16)/P(s,t,m+16,u,c)))~I(1-t),col="darkslateblue",lwd=Size)
lines(I(1/(G(t,m+32)/P(s,t,m+32,u,c)))~I(1-t),col="mediumorchid4",lwd=Size)
  dev.off()
                          



############################################  Plot 3
s <- 1
m <- 1.5
u <- 0.45
c <- 0.95
d <- 0.50

t <- seq(0.001,0.999,length.out=1000)
pdf("SimulationResults-d50-u45.pdf",width=12,height=4)
par(mfrow=c(1,3))
par(mar=c(5,5,5,5))
Size<-4

plot(P(s,t,m+1,u,c)~I(1-t),ylim=c(0,1),ylab="Percent female polygyny, P",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(P(s,t,m+2,u,c)~I(1-t),col="firebrick",lwd=Size)
lines(P(s,t,m+4,u,c)~I(1-t),col="peru",lwd=Size)
lines(P(s,t,m+8,u,c)~I(1-t),col="seagreen4",lwd=Size)
lines(P(s,t,m+16,u,c)~I(1-t),col="darkslateblue",lwd=Size)
lines(P(s,t,m+32,u,c)~I(1-t),col="mediumorchid4",lwd=Size)

plot(G(t,m+1)~I(1-t),ylim=c(0,1),ylab="Gini of wealth, G",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(G(t,m+2)~I(1-t),col="firebrick",lwd=Size)
lines(G(t,m+4)~I(1-t),col="peru",lwd=Size)
lines(G(t,m+8)~I(1-t),col="seagreen4",lwd=Size)
lines(G(t,m+16)~I(1-t),col="darkslateblue",lwd=Size)
lines(G(t,m+32)~I(1-t),col="mediumorchid4",lwd=Size)

plot(I(1/(G(t,m+1)/P(s,t,m+1,u,c)))~I(1-t),ylim=c(0,2),ylab="Ratio of P to G",xlab="Percent of males in poor wealth class",type="l",lwd=Size,cex.lab=Size2, cex.axis=Size2, cex.main=Size2, cex.sub=Size2)
lines(I(1/(G(t,m+2)/P(s,t,m+2,u,c)))~I(1-t),col="firebrick",lwd=Size)
lines(I(1/(G(t,m+4)/P(s,t,m+4,u,c)))~I(1-t),col="peru",lwd=Size)
lines(I(1/(G(t,m+8)/P(s,t,m+8,u,c)))~I(1-t),col="seagreen4",lwd=Size)
lines(I(1/(G(t,m+16)/P(s,t,m+16,u,c)))~I(1-t),col="darkslateblue",lwd=Size)
lines(I(1/(G(t,m+32)/P(s,t,m+32,u,c)))~I(1-t),col="mediumorchid4",lwd=Size)
  dev.off()
                   


