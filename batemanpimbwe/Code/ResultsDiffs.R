################################################################# Load libraries
 library(rethinking)

###################################################################### Load data
 sn <- read.csv("Bateman-SN.csv")
 sy <- read.csv("Bateman-SY.csv")
 pq <- read.csv("Bateman-SQWSY.csv")
 sq <- read.csv("Bateman-TWSY.csv")
 sqpq <- read.csv("Bateman-TSQWSY.csv")
 snpqsq <- read.csv("Bateman-FM.csv")

################################################################# Calculate WAIC
 WAICm <- c(sn$WAIC[1],sy$WAIC[1],sq$WAIC[1],pq$WAIC[1],sqpq$WAIC[1],snpqsq$WAIC[1])
 WAICf <- c(sn$WAIC[2],sy$WAIC[2],sq$WAIC[2],pq$WAIC[2],sqpq$WAIC[2],snpqsq$WAIC[2])

 dWAICm <- WAICm - min(WAICm)
 dWAICf <- WAICf - min(WAICf)

 sWAICm <- exp(-0.5 *dWAICm)
 sWAICf <- exp(-0.5 *dWAICf)

 wWAICm <- sWAICm/sum(sWAICm)
 wWAICf <- sWAICf/sum(sWAICf)

 sn$dWAIC <- c(dWAICm[1],dWAICf[1],NA)
 sn$wWAIC <- c(wWAICm[1],wWAICf[1],NA)

 sy$dWAIC <- c(dWAICm[2],dWAICf[2],NA)
 sy$wWAIC <- c(wWAICm[2],wWAICf[2],NA)

 sq$dWAIC <- c(dWAICm[3],dWAICf[3],NA)
 sq$wWAIC <- c(wWAICm[3],wWAICf[3],NA)

 pq$dWAIC <- c(dWAICm[4],dWAICf[4],NA)
 pq$wWAIC <- c(wWAICm[4],wWAICf[4],NA)

 sqpq$dWAIC <- c(dWAICm[5],dWAICf[5],NA)
 sqpq$wWAIC <- c(wWAICm[5],wWAICf[5],NA)

 snpqsq$dWAIC <- c(dWAICm[6],dWAICf[6],NA)
 snpqsq$wWAIC <- c(wWAICm[6],wWAICf[6],NA)

######### Set RS I the same for all models, only difference is due to resampling
 sn[,1:7] <- sy[,1:7] <- pq[,1:7] <- sq[,1:7] <- sqpq[,1:7] <- snpqsq[,1:7]
 sink("ResultsDiffs.txt")
cat(
"/toprule \n",
"/rowcolor[HTML]{FFFFFF} \n",
"/multicolumn{1}{l}{/cellcolor[HTML]{FFFFFF}/textbf{}} & /multicolumn{1}{c}{/cellcolor[HTML]{FFFFFF}/textbf{Principle 2, MS}} & /multicolumn{2}{c}{/cellcolor[HTML]{FFFFFF}/textbf{Principle 3, Bateman Gradient}} // /midrule \n",
"/rowcolor[HTML]{FFFFFF} \n",
"/multicolumn{1}{l}{/cellcolor[HTML]{FFFFFF}/textbf{}} & /textbf{I$_{/text{s}}$}  & /textbf{N} & /textbf{Y}  // /midrule \n",
paste0("/rowcolor[HTML]{EFEFEF}/cellcolor[HTML]{FFFFFF}/textbf{Spouse Number} & ",     sn$Ims[3]," (", sn$ImsL[3],", ",sn$ImsH[3],") &",   sn$Spouses[3]," (", sn$SpousesL[3] ,", ", sn$SpousesH[3],") & // \n"   ),
paste0("/rowcolor[HTML]{C0C0C0}/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years} & ",   sy$Ims[3]," (", sy$ImsL[3],", ",sy$ImsH[3],") &", "","", "" ,"", ""," &", sy$SpouseYears[3]," (", sy$SpouseYearsL[3],", ",  sy$SpouseYearsH[3] ,") // \n"   ),
paste0("/rowcolor[HTML]{EFEFEF}/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years: Timing Weights} & ",    sq$Ims[3]," (", sq$ImsL[3],", ",sq$ImsH[3],") &", "","", "" ,"", ""," &", sq$SpouseYears[3]," (", sq$SpouseYearsL[3],", ",  sq$SpouseYearsH[3] ,") //  \n"   ),
paste0("/rowcolor[HTML]{C0C0C0}/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years: Spousal Quality Weights} & ",    pq$Ims[3]," (", pq$ImsL[3],", ",pq$ImsH[3],") &", "","", "" ,"", ""," &", pq$SpouseYears[3]," (", pq$SpouseYearsL[3],", ",  pq$SpouseYearsH[3] ,") // \n"   ),
paste0("/rowcolor[HTML]{EFEFEF}/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years: Both Weights} & ", sqpq$Ims[3]," (", sqpq$ImsL[3],", ",sqpq$ImsH[3],") &" ,"","", "" ,"", ""," &", sqpq$SpouseYears[3]," (", sqpq$SpouseYearsL[3],", ",  sqpq$SpouseYearsH[3] ,") //  \n"   ),
paste0("/rowcolor[HTML]{C0C0C0}/cellcolor[HTML]{FFFFFF}/textbf{Full Model} & ",    snpqsq$Ims[3]," (", snpqsq$ImsL[3],", ",snpqsq$ImsH[3],") &",  snpqsq$Spouses[3]," (", snpqsq$SpousesL[3] ,", ", snpqsq$SpousesH[3],") &", snpqsq$SpouseYears[3]," (", snpqsq$SpouseYearsL[3],", ",  snpqsq$SpouseYearsH[3] ,") // /bottomrule"   )
  )
sink()
# Copy and paste result into a text editor, batch replace the front slash with a backslash, then into LaTeX it goes!












