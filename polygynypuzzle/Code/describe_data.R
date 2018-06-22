################################################################################
######################### Describe data
################################################################################
############################################################### Figure 1
# Load Libraries
source('./Code/project_support.R')

# Load Data
d<-read.csv("./Data/ExtraData/SCCS_with_MBM_corrections.csv")
d<-d[!is.na(d$Percent_md_women_polyg_md_v872),]
d<-d[!is.na(d$subsistence_type_categorized),]

myData <- aggregate(d$Percent_md_women_polyg_md_v872, by = list(Subsistence = d$subsistence_type_categorized), FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))

dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = myData$mean + 2*myData$se, ymin = myData$mean - 2*myData$se)

myData <- do.call(data.frame, myData)

myData$se <- myData$x.sd / sqrt(myData$x.n)

colnames(myData) <- c("Subsistence", "mean", "sd", "n", "se")
x<-myData[4,]             # Swap order of rows
myData[4,]<-myData[3,]    #
myData[3,] <- x           #

myData$names <- c("Foraging (28)", "Horticulture (26)", "Pastoral (14)", "Agriculture (24)")

######################################################## Frame A
p <- ggplot(data = myData, aes(x = factor(names,levels = myData$names), y = mean, fill = factor(names,levels = myData$names)))

p <- p + geom_bar(stat = "identity", position = dodge) +    theme_bw() +
    geom_errorbar(limits, position = dodge, width = 0.25)  +  scale_fill_manual(values=c("grey58","grey58", "grey58","grey58"),name=element_blank()) +
    theme(legend.position="none")+
    ylab("Mean percent female polygyny")  + theme(axis.text=element_text(size=14),
    axis.title=element_text(size=14))+  theme(axis.title.x=element_blank())
ggsave("PercentFemalePolygyny.pdf", p, height=4, width=6.5)

print("made PercentFemalePolygyny figure")


######################################################## Frame B

# Load Data
d<-read.csv("./Data/ExtraData/SCCS_with_MBM_corrections.csv")
d<-d[!is.na(d$Rev_Standard_Polygamy_code_v861),]
d<-d[!is.na(d$subsistence_type_categorized),]

colnames(d)[8]<-"Marriage"
colnames(d)[9]<-"Subsistence"

df<-data.frame(Subsistence=d$Subsistence, Marriage=d$Marriage)

df2<-table(df$Subsistence,df$Marriage)/rowSums(table(df$Subsistence,df$Marriage))
x<-df2[4,]         # Swap row order
df2[4,]<-df2[3,]   #
df2[3,] <- x       #

rownames(df2)<- c("Foraging (29)", "Horticulture (30)", "Pastoral (18)","Agriculture (29)")
colnames(df2)<- c("Monogamy prescribed", "Monogamy prefered", "Limited polygyny", "Full polygyny")

df2<-as.data.frame(df2)

p <- ggplot(df2, aes(Var1, Freq, fill = Var2)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values=c("grey18","grey38", "grey58","grey78"),name=element_blank()) +
    theme_bw() +     ylab("Frequency")  + theme(axis.text=element_text(size=14),legend.text=element_text(size=14),
    axis.title=element_text(size=14))+  theme(axis.title.x=element_blank())

ggsave("StandardPolygynyCodes.pdf", p, height=4, width=9)

print("made StandardPolygynyCodes figure")

################################################################################
######################### Figure 2
################################################################################

# Load Data
d<-read.csv("./Data/ExtraData/PolygynyProjectPopulationsGPS.csv")

Latitude<-d$Latitude
Longitude<-d$Longitude
Label<-d$Label
Voffset <-d$Voffset
Hoffset <-d$Hoffset


# Using ggplot, plot the Base World Map
mp <- NULL

mapWorld <- borders("world", colour="burlywood3", fill="burlywood") 

## ggplot2 map_data references maps::map as just "map", so you cant do maps and then rethinking! 

mp <- ggplot() +  mapWorld    + theme(axis.text=element_text(size=14),
     axis.title=element_text(size=14))+
labs(x="Longitude",y="Latitude")  + coord_map(projection = "mercator",xlim=c(-180,180),ylim=c(-65, 75))

# Now Layer the populations on top
mp <- mp + geom_point(aes(x=Longitude, y=Latitude) ,color="black", size=2)
mp <- mp +  geom_text_repel(aes(x=Longitude, y=Latitude, label = Label),color="black",size=4)
    mp

ggsave( "PPPMap.pdf", mp, width = 10, height = 6.5)

print("made PPPMap figure")

################################################################################
######################### Fig 4 Frame A
################################################################################
# Load Libraries
mu <- 0.4
delta <- 0.6

# Set Data
x <- c(rep(3,95),rep(18,100-95))
y <- c(rep(1,24),rep(6,100-24))

x_rich <- mean(x[which(x>3)])
x_poor <- mean(x[which(x<=3)])
wrx<-x_rich/x_poor

y_rich <- mean(y[which(y>1)])
y_poor <- mean(y[which(y<=1)])
wry<-y_rich/y_poor

print(paste("Wealth Ratio x = ",wrx))
print(paste("Wealth Ratio = ",wry))

print(paste("Gini x = ",Gini(x)))
print(paste("Gini y = ",Gini(y)))

pct_x <- length(x[(x>3)])/length(x)
pct_y <- length(y[(y>1)])/length(y)

print(paste("Pct Rich x = ",pct_x))
print(paste("Pct Rich y = ",pct_y))

pctf_x <- min(1,pct_x*(x_rich/x_poor)*((delta-mu)/delta) )
pctf_y <- min(1,pct_y*(y_rich/y_poor)*((delta-mu)/delta) )

print(paste("Pct Female Polygyny x = ",pctf_x))
print(paste("Pct Polygyny y = ",pctf_y))

# Plot Lorenze Curves
gp<-data.frame(x=Lc(x)$L,y=Lc(y)$L,p=Lc(y)$p)
v <- ggplot(gp, aes(x=p, y=y))
v<-v+geom_ribbon(data=gp,aes(ymin=x,ymax=p),alpha=0.3,fill="firebrick")+
      geom_ribbon(aes(ymin=y,ymax=p),alpha=0.3,fill="darkslateblue" )
v<-v + geom_line(data=melt(gp),aes(x=rep(gp$p,3), y=value, colour=variable,linetype = variable),size=1.5) +
      scale_colour_manual(name= "Lines",values=c("firebrick","darkslateblue","black")) +
      scale_linetype_manual(name= "Lines",values=c("twodash", "solid","solid"))
v<-v +  ylab("Cumulative wealth")
v<-v +  xlab("Cumulative population")  #+ geom_abline(intercept = 1, slope = -1,linetype="dotted")
v<-v +  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
v<-v + theme(legend.title=element_text(size=14),legend.text=element_text(size=14))
v

ggsave("SimulatedLorenz1.pdf",v,width=6,height=5)

print("made SimulatedLorenz1 figure")


################################################################################
######################### Fig 4 Frame B
################################################################################
# Load Libraries
mu <- 0.4
delta <- 0.5

# Set Data
x <- c(rep(1,95),rep(35,100-95))
y <- c(rep(1,55),rep(5,100-55))

x_rich <- mean(x[which(x>1)])
x_poor <- mean(x[which(x<=1)])
wrx<-x_rich/x_poor

y_rich <- mean(y[which(y>1)])
y_poor <- mean(y[which(y<=1)])
wry<-y_rich/y_poor

print(paste("Wealth Ratio x = ",wrx))
print(paste("Wealth Ratio = ",wry))

print(paste("Gini x = ",Gini(x)))
print(paste("Gini y = ",Gini(y)))

pct_x <- length(x[(x>1)])/length(x)
pct_y <- length(y[(y>1)])/length(y)

print(paste("Pct Rich x = ",pct_x))
print(paste("Pct Rich y = ",pct_y))

pctf_x <- min(1,pct_x*(x_rich/x_poor)*((delta-mu)/delta) )
pctf_y <- min(1,pct_y*(y_rich/y_poor)*((delta-mu)/delta) )

print(paste("Pct Female Polygyny x = ",pctf_x))
print(paste("Pct Polygyny y = ",pctf_y))

# Plot Lorenze Curves
gp<-data.frame(x=Lc(x)$L,y=Lc(y)$L,p=Lc(y)$p)
v <- ggplot(gp, aes(x=p, y=y))
v<-v+geom_ribbon(data=gp,aes(ymin=x,ymax=p),alpha=0.3,fill="firebrick")+
      geom_ribbon(aes(ymin=y,ymax=p),alpha=0.3,fill="darkslateblue")
v<-v + geom_line(data=melt(gp),aes(x=rep(gp$p,3), y=value, colour=variable,linetype = variable),size=1.5) +
       scale_colour_manual(name= "Lines",values=c("firebrick","darkslateblue","black")) +
       scale_linetype_manual(name= "Lines",values=c("twodash", "solid","solid"))
v<-v +  ylab("Cumulative wealth")
v<-v +  xlab("Cumulative population")  #+ geom_abline(intercept = 1, slope = -1,linetype="dotted")
v<-v +  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
v<-v + theme(legend.title=element_text(size=14),legend.text=element_text(size=14))
v

ggsave("SimulatedLorenz2.pdf",v,width=6,height=5)

print("made SimulatedLorenz2 figure")

