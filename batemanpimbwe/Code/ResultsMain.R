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

 sink("ResultsMain.txt")
cat(
"/toprule \n",
"/rowcolor[HTML]{FFFFFF} \n",
"/multicolumn{1}{l}{/cellcolor[HTML]{FFFFFF}/textbf{}} & /textbf{} & /multicolumn{1}{c}{/cellcolor[HTML]{FFFFFF}/textbf{Principle 2, MS}} & /multicolumn{4}{c}{/cellcolor[HTML]{FFFFFF}/textbf{Principle 3, Bateman Gradient}} & /multicolumn{3}{c}{/cellcolor[HTML]{FFFFFF}/textbf{WAIC}} // /midrule \n",
"/rowcolor[HTML]{FFFFFF} \n",
"/multicolumn{1}{l}{/cellcolor[HTML]{FFFFFF}/textbf{}} & /textbf{} & /textbf{I$_{/text{s}}$} & /textbf{Intercept} & /textbf{Age} & /textbf{N} & /textbf{Y}  & /textbf{WAIC} & /textbf{$/Delta$} & /textbf{$/omega$} // /midrule \n",
"/rowcolor[HTML]{EFEFEF} \n",
paste0("/cellcolor[HTML]{FFFFFF} & ","M &",     sn$Ims[1]," (", sn$ImsL[1],", ",sn$ImsH[1],") &", sn$Incpt[1]," (", sn$IncptL[1] ,", ", sn$IncptH[1],") &",  sn$Age[1]," (",  sn$AgeL[1],", ",  sn$AgeH[1],") &", sn$Spouses[1]," (", sn$SpousesL[1] ,", ", sn$SpousesH[1],") &", "","","","",  "" ," &",  round(sn$WAIC[1],2) ," &" ,  round(sn$dWAIC[1],2)," &" ,  round(sn$wWAIC[1],2), "// /cmidrule(l){2-15} \n"   ),
"/rowcolor[HTML]{EFEFEF} \n",
paste0("/multirow{-2}{*}{/cellcolor[HTML]{FFFFFF}/textbf{Spouse Number}} & ","F &",     sn$Ims[2]," (", sn$ImsL[2],", ",sn$ImsH[2],") &", sn$Incpt[2]," (", sn$IncptL[2] ,", ", sn$IncptH[2],") &",  sn$Age[2]," (",  sn$AgeL[2],", ",  sn$AgeH[2],") &", sn$Spouses[2]," (", sn$SpousesL[2] ,", ", sn$SpousesH[2],") &", ""," ", ""," ", ""," &",   round(sn$WAIC[2],2) ," &" ,  round(sn$dWAIC[2],2)," &" ,  round(sn$wWAIC[2],2) ,"// /midrule \n"  ),
"/rowcolor[HTML]{C0C0C0} \n",
paste0("/cellcolor[HTML]{FFFFFF} & ","M &",   sy$Ims[1]," (", sy$ImsL[1],", ",sy$ImsH[1],") &", sy$Incpt[1]," (", sy$IncptL[1] ,", ", sy$IncptH[1],") &",  sy$Age[1]," (",  sy$AgeL[1],", ",  sy$AgeH[1],") &", ""," ", "" ,"", ""," &", sy$SpouseYears[1]," (", sy$SpouseYearsL[1],", ",  sy$SpouseYearsH[1] ,") &",    round(sy$WAIC[1],2) ," &" ,  round(sy$dWAIC[1],2)," &" ,  round(sy$wWAIC[1],2), "// /cmidrule(l){2-15}\n"   ),
"/rowcolor[HTML]{C0C0C0} \n",
paste0("/multirow{-2}{*}{/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years}} & ","F& ",   sy$Ims[2]," (", sy$ImsL[2],", ",sy$ImsH[2],") &", sy$Incpt[2]," (", sy$IncptL[2] ,", ", sy$IncptH[2],") &",  sy$Age[2]," (",  sy$AgeL[2],", ",  sy$AgeH[2],") &", ""," ", "" ,"", ""," &", sy$SpouseYears[2]," (", sy$SpouseYearsL[2],", ",  sy$SpouseYearsH[2] ,") &",  round(sy$WAIC[2],2) ," &" ,  round(sy$dWAIC[2],2)," &" ,  round(sy$wWAIC[2],2),  "// /midrule \n"   ),
"/rowcolor[HTML]{EFEFEF} \n ",
paste0("/cellcolor[HTML]{FFFFFF} & ","M &",   sq$Ims[1]," (", sq$ImsL[1],", ",sq$ImsH[1],") &", sq$Incpt[1]," (", sq$IncptL[1] ,", ", sq$IncptH[1],") &",  sq$Age[1]," (",  sq$AgeL[1],", ",  sq$AgeH[1],") &", "","", "" ,"","","&", sq$SpouseYears[1]," (", sq$SpouseYearsL[1],", ",  sq$SpouseYearsH[1] ,") &",   round(sq$WAIC[1],2) ," &" ,  round(sq$dWAIC[1],2)," &" ,  round(sq$wWAIC[1],2), "// /cmidrule(l){2-15} \n"   ),
"/rowcolor[HTML]{EFEFEF} \n",
paste0("/multirow{-2}{*}{/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years: Timing Weights}} & ","F &",   sq$Ims[2]," (", sq$ImsL[2],", ",sq$ImsH[2],") &", sq$Incpt[2]," (", sq$IncptL[2] ,", ", sq$IncptH[2],") &",  sq$Age[2]," (",  sq$AgeL[2],", ",  sq$AgeH[2],") &", "","", "" ,"","","&", sq$SpouseYears[2]," (", sq$SpouseYearsL[2],", ",  sq$SpouseYearsH[2] ,") &",  round(sq$WAIC[2],2) ," &" ,  round(sq$dWAIC[2],2)," &" ,  round(sq$wWAIC[2],2) ,"// /midrule \n"  ),
"/rowcolor[HTML]{C0C0C0} \n",
paste0("/cellcolor[HTML]{FFFFFF} & ","M &",    pq$Ims[1]," (", pq$ImsL[1],", ",pq$ImsH[1],") &", pq$Incpt[1]," (", pq$IncptL[1] ,", ", pq$IncptH[1],") &",  pq$Age[1]," (",  pq$AgeL[1],", ",  pq$AgeH[1],") &", ""," ", "" ,"", ""," &", pq$SpouseYears[1]," (", pq$SpouseYearsL[1],", ",  pq$SpouseYearsH[1] ,") &",  round(pq$WAIC[1],2) ," &" ,  round(pq$dWAIC[1],2)," &" ,  round(pq$wWAIC[1],2), "// /cmidrule(l){2-15} \n"   ),
"/rowcolor[HTML]{C0C0C0} \n",
paste0("/multirow{-2}{*}{/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years: Spousal Quality Weights}} & ","F& ",     pq$Ims[2]," (", pq$ImsL[2],", ",pq$ImsH[2],") &", pq$Incpt[2]," (", pq$IncptL[2] ,", ", pq$IncptH[2],") &",  pq$Age[2]," (",  pq$AgeL[2],", ",  pq$AgeH[2],") &", ""," ", "" ,"", ""," &", pq$SpouseYears[2]," (", pq$SpouseYearsL[2],", ",  pq$SpouseYearsH[2] ,") &",    round(pq$WAIC[2],2) ," &" ,  round(pq$dWAIC[2],2)," &" ,  round(pq$wWAIC[2],2),  "// /midrule \n"   ),
"/rowcolor[HTML]{EFEFEF} \n ",
paste0("/cellcolor[HTML]{FFFFFF} & ","M &",    sqpq$Ims[1]," (", sqpq$ImsL[1],", ",sqpq$ImsH[1],") &", sqpq$Incpt[1]," (", sqpq$IncptL[1] ,", ", sqpq$IncptH[1],") &",  sqpq$Age[1]," (",  sqpq$AgeL[1],", ",  sqpq$AgeH[1],") &", ""," ", "" ,"", ""," &", sqpq$SpouseYears[1]," (", sqpq$SpouseYearsL[1],", ",  sqpq$SpouseYearsH[1] ,") &",   round(sqpq$WAIC[1],2) ," &" ,  round(sqpq$dWAIC[1],2)," &" ,  round(sqpq$wWAIC[1],2), "// /cmidrule(l){2-15} \n"   ),
"/rowcolor[HTML]{EFEFEF} \n",
paste0("/multirow{-2}{*}{/cellcolor[HTML]{FFFFFF}/textbf{Spouse Years: Both Weights}} & ","F &",    sqpq$Ims[2]," (", sqpq$ImsL[2],", ",sqpq$ImsH[2],") &", sqpq$Incpt[2]," (", sqpq$IncptL[2] ,", ", sqpq$IncptH[2],") &",  sqpq$Age[2]," (",  sqpq$AgeL[2],", ",  sqpq$AgeH[2],") &", "","", "" ,"","","&", sqpq$SpouseYears[2]," (", sqpq$SpouseYearsL[2],", ",  sqpq$SpouseYearsH[2] ,") &",  round(sqpq$WAIC[2],2) ," &" ,  round(sqpq$dWAIC[2],2)," &" ,  round(sqpq$wWAIC[2],2) ,"// /midrule \n"  ),
"/rowcolor[HTML]{C0C0C0} \n",
paste0("/cellcolor[HTML]{FFFFFF} & ","M &",   snpqsq$Ims[1]," (", snpqsq$ImsL[1],", ",snpqsq$ImsH[1],") &", snpqsq$Incpt[1]," (", snpqsq$IncptL[1] ,", ", snpqsq$IncptH[1],") &",  snpqsq$Age[1]," (",  snpqsq$AgeL[1],", ",  snpqsq$AgeH[1],") &", snpqsq$Spouses[1]," (", snpqsq$SpousesL[1] ,", ", snpqsq$SpousesH[1],") &", snpqsq$SpouseYears[1]," (", snpqsq$SpouseYearsL[1],", ",  snpqsq$SpouseYearsH[1] ,") &",    round(snpqsq$WAIC[1],2) ," &" ,  round(snpqsq$dWAIC[1],2)," &" ,  round(snpqsq$wWAIC[1],2), "// /cmidrule(l){2-15}"   ),
"/rowcolor[HTML]{C0C0C0} \n",
paste0("/multirow{-2}{*}{/cellcolor[HTML]{FFFFFF}/textbf{Full Model}} & ","F& ",   snpqsq$Ims[2]," (", snpqsq$ImsL[2],", ",snpqsq$ImsH[2],") &", snpqsq$Incpt[2]," (", snpqsq$IncptL[2] ,", ", snpqsq$IncptH[2],") &",  snpqsq$Age[2]," (",  snpqsq$AgeL[2],", ",  snpqsq$AgeH[2],") &", snpqsq$Spouses[2]," (", snpqsq$SpousesL[2] ,", ", snpqsq$SpousesH[2],") &", snpqsq$SpouseYears[2]," (", snpqsq$SpouseYearsL[2],", ",  snpqsq$SpouseYearsH[2] ,") &",    round(snpqsq$WAIC[2],2) ," &" ,  round(snpqsq$dWAIC[2],2)," &" ,  round(snpqsq$wWAIC[2],2),  "// /bottomrule "   )
)
sink()
# Copy and paste result into a text editor, batch replace the front slash with a backslash, then into LaTeX it goes!









