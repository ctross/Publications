###################################################################### Load Libraries
rm(list=ls())
rstan_options(auto_write = TRUE)

###################################################################### Read Data
d1<-read.csv("./Data/Populations/Kipsigis-BorgerhoffMulder.csv")
d2<-read.csv("./Data/Populations/Sangu-McElreath.csv")
d3<-read.csv("./Data/Populations/Makushi-Schacht.csv")
d4<-read.csv("./Data/Populations/Maasai-Caudell.csv")
d5<-read.csv("./Data/Populations/Maasai2-Caudell.csv")
d6<-read.csv("./Data/Populations/Ache-Hill.csv")
d7<-read.csv("./Data/Populations/Agta-Headland.csv")
d8<-read.csv("./Data/Populations/Maya-Cortez.csv")
d9<-read.csv("./Data/Populations/Matsigenka-Bunce.csv")
d10<-read.csv("./Data/Populations/Lamalera-Nolin.csv")
d11<-read.csv("./Data/Populations/Kung-HowellDraper.csv")
d12<-read.csv("./Data/Populations/Koore-Caudell.csv")
d13<-read.csv("./Data/Populations/Himba-Scelza.csv")
d14<-read.csv("./Data/Populations/Hadza-Marlow.csv")
d15<-read.csv("./Data/Populations/Chagga-Caudell.csv")
d16<-read.csv("./Data/Populations/Aka-Hewlett.csv")
d17<-read.csv("./Data/Populations/Maya-Kramer.csv")
d18<-read.csv("./Data/Populations/Mayangna-Koster.csv")
d19<-read.csv("./Data/Populations/Pume-Kramer.csv")
d20<-read.csv("./Data/Populations/Sidama-Caudell.csv")
d21<-read.csv("./Data/Populations/Sukuma-BorgerhoffMulder.csv")
d22<-read.csv("./Data/Populations/Tsimane-Godoy.csv")
d23<-read.csv("./Data/Populations/Dolgan-Ziker.csv")
d24<-read.csv("./Data/Populations/Bangladesh-Shenck.csv")
d25<-read.csv("./Data/Populations/English-Clark.csv")
d26<-read.csv("./Data/Populations/Pimbwe-BorgerhoffMulder.csv")
d27<-read.csv("./Data/Populations/Tsimane-Gurven.csv")
d28<-read.csv("./Data/Populations/Krummhorn-Willfuehr.csv")
d29<-read.csv("./Data/Populations/Poland-Colleran.csv")

###################################################################### Load Data
MaxExposure <- 47  # Maximum year of exposure to risk of RS
P<-29              # Populations
M<-13              # Age at first risk of RS

################################### Load the data for each population
# Note that data are shifted to remove zeros,
# normally by adding the smallest observed value
# Sometimes we just have to add a constant that is small on the measurement scale

Kipsigis_Exposure <- ifelse((d1$age-M)>MaxExposure,MaxExposure,(d1$age-M))
Kipsigis_Wives <- d1$wives
Kipsigis_RS <- d1$soff
Kipsigis_Cows <- (d1$cows +  min(d1$cows[which(d1$cows>0)],na.rm=T))/sd(d1$cows+ min(d1$cows[which(d1$cows>0)],na.rm=T))
Kipsigis_Land <- (d1$totalacres+ min(d1$totalacres[which(d1$totalacres>0)],na.rm=T))/sd(d1$totalacres +  min(d1$totalacres[which(d1$totalacres>0)],na.rm=T))

Sangu_Exposure <- ifelse((d2$Age-M)>MaxExposure,MaxExposure,(d2$Age-M))
Sangu_Wives <- d2$Wives
Sangu_RS <- d2$rs
Sangu_Cows <- (d2$Cattle+min(d2$Cattle[which(d2$Cattle>0)],na.rm=T))/sd(d2$Cattle+min(d2$Cattle[which(d2$Cattle>0)],na.rm=T))
LLsangu <-  d2$AcresCorn + d2$AcresRice
Sangu_Land <- (d2$AcresCorn + d2$AcresRice + min(LLsangu[which(LLsangu>0)],na.rm=T))/sd((d2$AcresCorn + d2$AcresRice + min(LLsangu[which(LLsangu>0)],na.rm=T)))

Maasai_Exposure <- ifelse((d4$age-M)>MaxExposure,MaxExposure,(d4$age-M))
Maasai_Wives <- d4$num_wives
Maasai_RS <- d4$rs
Maasai_Cows <- d4$cattle/sd(d4$cattle)
Maasai_Land <- (d4$farming_land + min(d4$farming_land[which(d4$farming_land>0)],na.rm=T))/sd(d4$farming_land + min(d4$farming_land[which(d4$farming_land>0)],na.rm=T))

Maasai2_Exposure <- ifelse((d5$age-M)>MaxExposure,MaxExposure,(d5$age-M))
Maasai2_Wives <- d5$timesmarried
Maasai2_RS <- d5$rs
Maasai2_Cows <- ((d5$livestock_cattle_total+min(d5$livestock_cattle_total[which(d5$livestock_cattle_total>0)],na.rm=T))/(1+d5$bro_living))/sd((d5$livestock_cattle_total+min(d5$livestock_cattle_total[which(d5$livestock_cattle_total>0)],na.rm=T))/(1+d5$bro_living))
Maasai2_Land <- (d5$farmed_land+min(d5$farmed_land[which(d5$farmed_land>0)],na.rm=T))/sd(d5$farmed_land+min(d5$farmed_land[which(d5$farmed_land>0)],na.rm=T))

Koore_Exposure <- ifelse((d12$age-M)>MaxExposure,MaxExposure,(d12$age-M))
Koore_Wives <- d12$num_wives
Koore_RS <- d12$rs
Koore_Cows <- (d12$cattle+min(d12$cattle[which(d12$cattle>0)],na.rm=T))/sd(d12$cattle+min(d12$cattle[which(d12$cattle>0)],na.rm=T))
Koore_Land <- (d12$farming_land+min(d12$farming_land[which(d12$farming_land>0)],na.rm=T))/sd(d12$farming_land+min(d12$farming_land[which(d12$farming_land>0)],na.rm=T))

Chagga_Exposure <- ifelse((d15$age-M)>MaxExposure,MaxExposure,(d15$age-M))
Chagga_Wives <- d15$timesmarried
Chagga_RS <- d15$rs
Chagga_Cows <- (d15$livestock_cattle_total+min(d15$livestock_cattle_total[which(d15$livestock_cattle_total>0)],na.rm=T))/sd(d15$livestock_cattle_total+min(d15$livestock_cattle_total[which(d15$livestock_cattle_total>0)],na.rm=T))
Chagga_Land <- (d15$farmed_land+min(d15$farmed_land[d15$farmed_land>0],na.rm=T))/sd(d15$farmed_land+min(d15$farmed_land[d15$farmed_land>0],na.rm=T))

Sidama_Exposure <- ifelse((d20$ego_age-M)>MaxExposure,MaxExposure,(d20$ego_age-M))
Sidama_Wives <- d20$ego_spouses
Sidama_RS <- d20$rs
Sidama_Cows <- (d20$ego_cattle+min(d20$ego_cattle[which(d20$ego_cattle>0)],na.rm=T))/sd(d20$ego_cattle+min(d20$ego_cattle[which(d20$ego_cattle>0)],na.rm=T))
Sidama_Land <- (d20$ego_land+min(d20$ego_land[which(d20$ego_land>0)],na.rm=T))/sd(d20$ego_land+min(d20$ego_land[which(d20$ego_land>0)],na.rm=T))

Sukuma_Exposure <- ifelse((d21$Age-M)>MaxExposure,MaxExposure,(d21$Age-M))
Sukuma_Wives <- d21$wives
Sukuma_RS <- d21$kidsgt5
Sukuma_Cows <- (d21$cattle+min(d21$cattle[which(d21$cattle>0)],na.rm=T))/sd(d21$cattle+min(d21$cattle[which(d21$cattle>0)],na.rm=T))
Sukuma_Land <- (d21$hec+min(d21$hec[which(d21$hec>0)],na.rm=T))/sd(d21$hec+min(d21$hec[which(d21$hec>0)],na.rm=T))

Pimbwe_Exposure <- ifelse((d26$Age-M)>MaxExposure,MaxExposure,(d26$Age-M))
Pimbwe_Wives <- d26$Wives
Pimbwe_RS <- d26$soff
Pimbwe_Cows <- (d26$MeanLivestock+min(d26$MeanLivestock[which(d26$MeanLivestock>0)],na.rm=T))/sd(d26$MeanLivestock+min(d26$MeanLivestock[which(d26$MeanLivestock>0)],na.rm=T))
Pimbwe_Land <- (d26$MeanLandHect+min(d26$MeanLandHect[which(d26$MeanLandHect>0)],na.rm=T))/sd(d26$MeanLandHect+min(d26$MeanLandHect[which(d26$MeanLandHect>0)],na.rm=T))

Makushi_Exposure <- ifelse((d3$AGE-M)>MaxExposure,MaxExposure,(d3$AGE-M))
Makushi_Wives <- d3$BabyMommas
Makushi_RS <- d3$RS
Makushi_Land <- (d3$FARMLAND+min(d3$FARMLAND[which(d3$FARMLAND>0)],na.rm=T))/sd(d3$FARMLAND+min(d3$FARMLAND[which(d3$FARMLAND>0)],na.rm=T))
#Makushi_Income <- (d3$MONTLY_INCOME+min(d3$MONTLY_INCOME[which(d3$MONTLY_INCOME>0)],na.rm=T))/sd(d3$MONTLY_INCOME+min(d3$MONTLY_INCOME[which(d3$MONTLY_INCOME>0)],na.rm=T))

Maya_Exposure <- ifelse((d8$Age-M)>MaxExposure,MaxExposure,(d8$Age-M))
Maya_Wives <- d8$Wives
Maya_RS <- d8$RS
Maya_Land <- d8$FieldAcres/sd(d8$FieldAcres)
# Maya_Harvest <- (d8$MaizeHarvest+min(d8$MaizeHarvest[which(d8$MaizeHarvest>0)],na.rm=T))/sd(d8$MaizeHarvest+min(d8$MaizeHarvest[which(d8$MaizeHarvest>0)],na.rm=T))

Maya2_Exposure <- ifelse((d17$Age-M)>MaxExposure,MaxExposure,(d17$Age-M))
Maya2_Wives <- d17$Wives
Maya2_RS <- d17$RS5
Maya2_Land <- (d17$Hectares+min(d17$Hectares[which(d17$Hectares>0)],na.rm=T))/sd(d17$Hectares+min(d17$Hectares[which(d17$Hectares>0)],na.rm=T))
Maya2_Transport <- (d17$TruckMoto+min(d17$TruckMoto[which(d17$TruckMoto>0)],na.rm=T) )/sd(d17$TruckMoto+min(d17$TruckMoto[which(d17$TruckMoto>0)],na.rm=T) )

Bangladesh_Exposure <- ifelse((d24$Age-M)>MaxExposure,MaxExposure,(d24$Age-M))
Bangladesh_Wives <-  d24$NUM_MAR
Bangladesh_RS <-  d24$rs
Bangladesh_Land <- (d24$LAND_DECIMAL+min(d24$LAND_DECIMAL[which(d24$LAND_DECIMAL>0)],na.rm=T) )/sd(d24$LAND_DECIMAL+min(d24$LAND_DECIMAL[which(d24$LAND_DECIMAL>0)],na.rm=T) )
# Bangladesh_Income <- (d24$AN_INCOME)/sd(d24$AN_INCOME)

Tsimane_Exposure <- ifelse((d22$age-M)>MaxExposure,MaxExposure,(d22$age-M))
Tsimane_Wives <- d22$spouses
Tsimane_RS <- d22$rsa
Tsimane_Wealth <- (d22$totWealth)/sd(d22$totWealth)
Tsimane_Land <- (d22$Hect+min(d22$Hect[which(d22$Hect>0)],na.rm=T) )/sd(d22$Hect+min(d22$Hect[which(d22$Hect>0)],na.rm=T) )

Tsimane2_Exposure <- ifelse((d27$Age-M)>MaxExposure,MaxExposure,(d27$Age-M))
Tsimane2_Wives <- d27$Wives
Tsimane2_RS <- d27$rs
Tsimane2_Wealth <- d27$Wealth/sd(d27$Wealth)
#LLtsimane <- (d27$WageIncome+d27$ProduceIncome)
#Tsimane2_Income <- (LLtsimane + min(LLtsimane[which(LLtsimane >0)],na.rm=T))/sd(LLtsimane +min(LLtsimane[which(LLtsimane >0)],na.rm=T))

Lamalera_Exposure <- ifelse(round(d10$age-M,0)>MaxExposure,MaxExposure,round(d10$age-M,0))
Lamalera_Wives  <- d10$num_wives
Lamalera_RS  <- d10$rs
Lamalera_BoatShares <- (d10$shares06+ min(d10$shares06[which(d10$shares06 >0)],na.rm=T))/sd(d10$shares06+ min(d10$shares06[which(d10$shares06 >0)],na.rm=T))
# Lamalera_WealthCategories <- (d10$wealth)/sd(d10$wealth)

Matsigenka_Exposure <- ifelse((d9$Age-M)>MaxExposure,MaxExposure,(d9$Age-M))
Matsigenka_Wives <- d9$BabyMommas
Matsigenka_RS <- d9$RS
Matsigenka_Boat <- d9$Boat  + 0.2      # Constant is arbitrary choice, check robustness to varying this
# Matsigenka_ComputerTV <- d9$ComputerTV

Dolgan_Exposure <- ifelse((d23$Age-M)>MaxExposure,MaxExposure,(d23$Age-M))
Dolgan_Wives <- d23$Wives
Dolgan_RS <- d23$Rsliving
Dolgan_Transport <- d23$Transport + 0.2  # Constant is arbitrary choice, check robustness to varying this
Dolgan_Territory <- d23$Territory + 0.2  #

Himba_Exposure <- ifelse((d13$Age-M)>MaxExposure,MaxExposure,(d13$Age-M))
Himba_Wives <- d13$Wives
Himba_RS <- d13$rs
Himba_Cows <- (d13$Cows+1)/sd(d13$Cows+1)

Mayangna_Exposure <- ifelse((d18$age-M)>MaxExposure,MaxExposure,(d18$age-M))
Mayangna_Wives <- d18$nwives
Mayangna_RS <- d18$rs
Mayangna_Wealth <- d18$household.wealth..rival./sd(d18$household.wealth..rival.)

English_Exposure <- ifelse((d25$Age-M)>MaxExposure,MaxExposure,(d25$Age-M))
English_Wives <- d25$Wives
English_RS <- d25$rs
English_Wealth <- d25$WealthRatio/sd(d25$WealthRatio)

Ache_Exposure <- ifelse((d6$Age-M)>MaxExposure,MaxExposure,(d6$Age-M))
Ache_Wives <- d6$Wives
Ache_RS <- d6$rs
Ache_Weight <- (d6$Weight-min(d6$Weight)+1)/sd((d6$Weight-min(d6$Weight)+1))

Agta_Exposure <- ifelse((d7$Age-M)>MaxExposure,MaxExposure,(d7$Age-M))
Agta_Wives <- d7$BabyMammas15
Agta_RS <- d7$SumRS15
Agta_Weight <- (d7$AvgOfWeight-min(d7$AvgOfWeight)+1)/sd((d7$AvgOfWeight-min(d7$AvgOfWeight)+1))

Pume_Exposure <- ifelse((d19$Age-M)>MaxExposure,MaxExposure,(d19$Age-M))
Pume_Wives <- d19$Wives
Pume_RS <- d19$rs
Pume_Weight <- (d19$Weight-min(d19$Weight)+1)/sd((d19$Weight-min(d19$Weight)+1))

Kung_Exposure <- ifelse((d11$age68-M)>MaxExposure,MaxExposure,(d11$age68-M))
Kung_Wives <- d11$spouses.polygynoushave2wivesbyfiat
Kung_RS <- d11$Rsalive
Kung_Weight <- (d11$wtkgs-min(d11$wtkgs)+1)/sd((d11$wtkgs-min(d11$wtkgs)+1))

Hadza_Exposure <- ifelse((d14$age-M)>MaxExposure,MaxExposure,(d14$age-M))
Hadza_Wives <- d14$BabyMommas
Hadza_RS <- d14$RSestimate
Hadza_Weight <- (d14$weight_kg-min(d14$weight_kg)+1)/sd((d14$weight_kg-min(d14$weight_kg)+1))

Aka_Exposure <- ifelse((d16$age-M)>MaxExposure,MaxExposure,(d16$age-M))
Aka_Wives <- d16$totalwives
Aka_RS <- d16$rs
Aka_Weight <- (d16$WtKg-min(d16$WtKg)+1)/sd((d16$WtKg-min(d16$WtKg)+1))

Krummhorn_Exposure <- ifelse((d28$AOD-M)>MaxExposure,MaxExposure,(d28$AOD-M))
Krummhorn_Wives <- d28$Nmarr
Krummhorn_RS <- d28$Nsur15
Krummhorn_Land <- (d28$grasen + 1)/sd(d28$grasen + 1)

Polish_Exposure <- ifelse((d29$age-M)>MaxExposure,MaxExposure,(d29$age-M))
Polish_Wives <- d29$Nwives
Polish_RS <- d29$noSurvived15
Polish_Land <- (d29$arableLand + 0.05)/sd(d29$arableLand + 0.05)
Polish_Cows <- (d29$cows + d29$bulls + 1)/sd(d29$cows + d29$bulls + 1)


N <- c(length(Kipsigis_Wives),  # 1
       length(Sangu_Wives),     # 1
       length(Maasai_Wives),    # 1
       length(Maasai2_Wives),   # 1
       length(Koore_Wives),     # 1
       length(Chagga_Wives),    # 4
       length(Sidama_Wives),    # 1
       length(Sukuma_Wives),    # 1
       length(Pimbwe_Wives),    # 3
       length(Maya2_Wives),     # 4
       length(Tsimane_Wives),   # 3
       length(Dolgan_Wives),    # 2
       length(Polish_Wives),    # 2
       length(Maya_Wives),      # 4
       length(Matsigenka_Wives),# 3
       length(Bangladesh_Wives),# 4
       length(Tsimane2_Wives),  # 3
       length(Makushi_Wives),   # 3
       length(Lamalera_Wives),  # 2
       length(Himba_Wives),     # 1
       length(Mayangna_Wives),  # 3
       length(English_Wives),   # 4
       length(Ache_Wives),      # 2
       length(Agta_Wives),      # 2
       length(Aka_Wives),       # 2
       length(Hadza_Wives),     # 2
       length(Kung_Wives),      # 2
       length(Pume_Wives),      # 2
       length(Krummhorn_Wives)  # 4
       )

X <- max(N)

Rival_1 <- Rival_2 <- RS <- Wives <- Exposure <- matrix(999999, nrow=X,ncol=P)

Rival_1[1:length(Kipsigis_Cows),1] <- Kipsigis_Cows
Rival_1[1:length(Sangu_Cows),2]    <- Sangu_Cows
Rival_1[1:length(Maasai_Cows),3]   <- Maasai_Cows
Rival_1[1:length(Maasai2_Cows),4]  <- Maasai2_Cows
Rival_1[1:length(Koore_Cows),5]    <- Koore_Cows
Rival_1[1:length(Chagga_Cows),6]   <- Chagga_Cows
Rival_1[1:length(Sidama_Cows),7]   <- Sidama_Cows
Rival_1[1:length(Sukuma_Cows),8]   <- Sukuma_Cows
Rival_1[1:length(Pimbwe_Cows),9]   <- Pimbwe_Cows
Rival_1[1:length(Maya2_Land ),10]   <- Maya2_Land
Rival_1[1:length(Tsimane_Wealth),11]   <- Tsimane_Wealth
Rival_1[1:length(Dolgan_Transport),12]   <-  Dolgan_Transport
Rival_1[1:length(Polish_Cows),13]   <-  Polish_Cows
Rival_1[1:length(Maya_Land),14]   <- Maya_Land
Rival_1[1:length(Matsigenka_Boat),15]   <-  Matsigenka_Boat
Rival_1[1:length(Bangladesh_Land),16]   <- Bangladesh_Land
Rival_1[1:length(Tsimane2_Wealth),17]   <- Tsimane2_Wealth
Rival_1[1:length(Makushi_Land),18]   <- Makushi_Land
Rival_1[1:length(Lamalera_BoatShares),19]   <- Lamalera_BoatShares
Rival_1[1:length(Himba_Cows),20]   <-  Himba_Cows
Rival_1[1:length(Mayangna_Wealth),21]   <-  Mayangna_Wealth
Rival_1[1:length(English_Wealth),22]   <-  English_Wealth
Rival_1[1:length(Ache_Weight),23]   <-  Ache_Weight
Rival_1[1:length(Agta_Weight),24]   <-  Agta_Weight
Rival_1[1:length(Aka_Weight),25]   <-  Aka_Weight
Rival_1[1:length(Hadza_Weight),26]   <-  Hadza_Weight
Rival_1[1:length(Kung_Weight),27]   <-  Kung_Weight
Rival_1[1:length(Pume_Weight),28]   <-  Pume_Weight
Rival_1[1:length(Krummhorn_Land),29]   <-  Krummhorn_Land

Rival_2[1:length(Kipsigis_Land),1] <- Kipsigis_Land
Rival_2[1:length(Sangu_Land),2]    <- Sangu_Land
Rival_2[1:length(Maasai_Land),3]   <- Maasai_Land
Rival_2[1:length(Maasai2_Land),4]  <- Maasai2_Land
Rival_2[1:length(Koore_Land),5]    <- Koore_Land
Rival_2[1:length(Chagga_Land),6]   <- Chagga_Land
Rival_2[1:length(Sidama_Land),7]   <- Sidama_Land
Rival_2[1:length(Sukuma_Land),8]   <- Sukuma_Land
Rival_2[1:length(Pimbwe_Land),9]   <- Pimbwe_Land
Rival_2[1:length(Maya2_Transport),10]   <- Maya2_Transport
Rival_2[1:length(Tsimane_Land),11]   <- Tsimane_Land
Rival_2[1:length(Dolgan_Territory),12]   <-  Dolgan_Territory
Rival_2[1:length(Polish_Land),13]   <-  Polish_Land
Rival_2[1:length(Maya_Land),14]      <- rep(0,length(Maya_Land))
Rival_2[1:length(Matsigenka_Boat),15]   <-  rep(0,length(Matsigenka_Boat))
Rival_2[1:length(Bangladesh_Land),16]   <- rep(0,length(Bangladesh_Land))
Rival_2[1:length(Tsimane2_Wealth),17]   <- rep(0,length(Tsimane2_Wealth))
Rival_2[1:length(Makushi_Land),18]     <- rep(0,length(Makushi_Land))
Rival_2[1:length(Lamalera_BoatShares),19]   <- rep(0,length(Lamalera_BoatShares))
Rival_2[1:length(Himba_Cows),20]   <-  rep(0,length(Himba_Cows))
Rival_2[1:length(Mayangna_Wealth),21]   <-  rep(0,length(Mayangna_Wealth))
Rival_2[1:length(English_Wealth),22]   <-  rep(0,length(English_Wealth))
Rival_2[1:length(Ache_Weight),23]   <-  rep(0,length(Ache_Weight))
Rival_2[1:length(Agta_Weight),24]   <-  rep(0,length(Agta_Weight))
Rival_2[1:length(Aka_Weight),25]   <-  rep(0,length(Aka_Weight))
Rival_2[1:length(Hadza_Weight),26]   <-  rep(0,length(Hadza_Weight))
Rival_2[1:length(Kung_Weight),27]   <-  rep(0,length(Kung_Weight))
Rival_2[1:length(Pume_Weight),28]   <-  rep(0,length(Pume_Weight))
Rival_2[1:length(Krummhorn_Land),29]   <-  rep(0,length(Krummhorn_Land))

RS[1:length(Kipsigis_RS),1] <- Kipsigis_RS
RS[1:length(Sangu_RS),2]    <- Sangu_RS
RS[1:length(Maasai_RS),3]   <- Maasai_RS
RS[1:length(Maasai2_RS),4]  <- Maasai2_RS
RS[1:length(Koore_RS),5]    <- Koore_RS
RS[1:length(Chagga_RS),6]   <- Chagga_RS
RS[1:length(Sidama_RS),7]   <- Sidama_RS
RS[1:length(Sukuma_RS),8]   <- Sukuma_RS
RS[1:length(Pimbwe_RS),9]   <- Pimbwe_RS
RS[1:length(Maya2_RS),10]   <- Maya2_RS
RS[1:length(Tsimane_RS),11]  <- Tsimane_RS
RS[1:length(Dolgan_RS),12]   <- Dolgan_RS
RS[1:length(Polish_RS),13]   <- Polish_RS
RS[1:length(Maya_RS),14]   <- Maya_RS
RS[1:length(Matsigenka_RS),15]   <- Matsigenka_RS
RS[1:length(Bangladesh_RS),16] <- Bangladesh_RS
RS[1:length(Tsimane2_RS),17]   <- Tsimane2_RS
RS[1:length(Makushi_RS),18]<-  Makushi_RS
RS[1:length(Lamalera_RS),19]   <- Lamalera_RS
RS[1:length(Himba_RS),20]   <- Himba_RS
RS[1:length(Mayangna_RS),21]   <- Mayangna_RS
RS[1:length(English_RS),22]   <- English_RS
RS[1:length(Ache_RS),23]   <- Ache_RS
RS[1:length(Agta_RS),24]   <- Agta_RS
RS[1:length(Aka_RS),25]   <- Aka_RS
RS[1:length(Hadza_RS),26]   <- Hadza_RS
RS[1:length(Kung_RS),27]   <- Kung_RS
RS[1:length(Pume_RS),28]   <- Pume_RS
RS[1:length(Krummhorn_RS),29]   <- Krummhorn_RS

Wives[1:length(Kipsigis_Wives),1] <- Kipsigis_Wives
Wives[1:length(Sangu_Wives),2]    <- Sangu_Wives
Wives[1:length(Maasai_Wives),3]   <- Maasai_Wives
Wives[1:length(Maasai2_Wives),4]  <- Maasai2_Wives
Wives[1:length(Koore_Wives),5]    <- Koore_Wives
Wives[1:length(Chagga_Wives),6]   <- Chagga_Wives
Wives[1:length(Sidama_Wives),7]   <- Sidama_Wives
Wives[1:length(Sukuma_Wives),8]   <- Sukuma_Wives
Wives[1:length(Pimbwe_Wives),9]   <- Pimbwe_Wives
Wives[1:length(Maya2_Wives),10]   <- Maya2_Wives
Wives[1:length(Tsimane_Wives),11]  <- Tsimane_Wives
Wives[1:length(Dolgan_Wives),12]   <- Dolgan_Wives
Wives[1:length(Polish_Wives),13]   <- Polish_Wives
Wives[1:length(Maya_Wives),14]   <- Maya_Wives
Wives[1:length(Matsigenka_Wives),15]   <- Matsigenka_Wives
Wives[1:length(Bangladesh_Wives),16] <- Bangladesh_Wives
Wives[1:length(Tsimane2_Wives),17]   <- Tsimane2_Wives
Wives[1:length(Makushi_Wives),18]<-  Makushi_Wives
Wives[1:length(Lamalera_Wives),19]   <- Lamalera_Wives
Wives[1:length(Himba_Wives),20]   <- Himba_Wives
Wives[1:length(Mayangna_Wives),21]   <- Mayangna_Wives
Wives[1:length(English_Wives),22]   <- English_Wives
Wives[1:length(Ache_Wives),23]   <- Ache_Wives
Wives[1:length(Agta_Wives),24]   <- Agta_Wives
Wives[1:length(Aka_Wives),25]   <- Aka_Wives
Wives[1:length(Hadza_Wives),26]   <- Hadza_Wives
Wives[1:length(Kung_Wives),27]   <- Kung_Wives
Wives[1:length(Pume_Wives),28]   <- Pume_Wives
Wives[1:length(Krummhorn_Wives),29]   <- Krummhorn_Wives

Exposure[1:length(Kipsigis_Exposure),1] <- Kipsigis_Exposure
Exposure[1:length(Sangu_Exposure),2]    <- Sangu_Exposure
Exposure[1:length(Maasai_Exposure),3]   <- Maasai_Exposure
Exposure[1:length(Maasai2_Exposure),4]  <- Maasai2_Exposure
Exposure[1:length(Koore_Exposure),5]    <- Koore_Exposure
Exposure[1:length(Chagga_Exposure),6]   <- Chagga_Exposure
Exposure[1:length(Sidama_Exposure),7]   <- Sidama_Exposure
Exposure[1:length(Sukuma_Exposure),8]   <- Sukuma_Exposure
Exposure[1:length(Pimbwe_Exposure),9]   <- Pimbwe_Exposure
Exposure[1:length(Maya2_Exposure),10]   <- Maya2_Exposure
Exposure[1:length(Tsimane_Exposure),11]  <- Tsimane_Exposure
Exposure[1:length(Dolgan_Exposure),12]   <- Dolgan_Exposure
Exposure[1:length(Polish_Exposure),13]   <- Polish_Exposure
Exposure[1:length(Maya_Exposure),14]   <- Maya_Exposure
Exposure[1:length(Matsigenka_Exposure),15]   <- Matsigenka_Exposure
Exposure[1:length(Bangladesh_Exposure),16] <- Bangladesh_Exposure
Exposure[1:length(Tsimane2_Exposure),17]   <- Tsimane2_Exposure
Exposure[1:length(Makushi_Exposure ),18]<-  Makushi_Exposure
Exposure[1:length(Lamalera_Exposure),19]   <- Lamalera_Exposure
Exposure[1:length(Himba_Exposure),20]   <- Himba_Exposure
Exposure[1:length(Mayangna_Exposure),21]   <- Mayangna_Exposure
Exposure[1:length(English_Exposure),22]   <- English_Exposure
Exposure[1:length(Ache_Exposure),23]   <- Ache_Exposure
Exposure[1:length(Agta_Exposure),24]   <- Agta_Exposure
Exposure[1:length(Aka_Exposure),25]   <- Aka_Exposure
Exposure[1:length(Hadza_Exposure),26]   <- Hadza_Exposure
Exposure[1:length(Kung_Exposure),27]   <- Kung_Exposure
Exposure[1:length(Pume_Exposure),28]   <- Pume_Exposure
Exposure[1:length(Krummhorn_Exposure),29]   <- Krummhorn_Exposure

Proxy <- c(rep(1,13),rep(0,P-13))

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

################################################################################
# Step 1: Fit Model
################################################################################
model_dat<-list(
  P=P,
  X=X,
  N=N,
  MaxExposure=MaxExposure,
  Proxy=Proxy,
  Exposure=Exposure,
  Wives=Wives,
  RS=RS,
  Rival_1=Rival_1,
  Rival_2=Rival_2
)

fitPoly <- stan(file="./Code/model.stan",
  data = model_dat,thin=1, iter = 2000, 
  warmup=1000, chains = 2, cores = 2,
  refresh=50, seed=8675309, control=list(max_treedepth=13))

save.image("ModelResults.RData")
















