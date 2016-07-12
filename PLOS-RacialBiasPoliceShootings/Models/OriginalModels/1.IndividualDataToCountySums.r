
############################# Read Data
d <- read.csv(file.choose()) # U.S._Police_Shootings_Data_(Cleaned).csv

############################# This code is used to translate the individual level data in the file named
# U.S._Police_Shootings_Data_(Cleaned).csv, into the county-level sums in the file named
# MapFileData-WithCountyResultsAndCovariates.csv

#############################
Results<-c()
Results <- ifelse(d$Armed.or.Unarmed.=="Armed" & d$Ethnicity=="Black", "BlackArmed",NA)
Results <- ifelse(d$Armed.or.Unarmed.=="Unarmed" & d$Ethnicity=="Black", "BlackUnarmed",Results)
Results <- ifelse(d$Armed.or.Unarmed.=="Armed" & d$Ethnicity=="Hispanic", "HispanicArmed",Results)
Results <- ifelse(d$Armed.or.Unarmed.=="Unarmed" & d$Ethnicity=="Hispanic", "HispanicUnarmed",Results)
Results <- ifelse(d$Armed.or.Unarmed.=="Armed" & d$Ethnicity=="White", "WhiteArmed",Results)
Results <- ifelse(d$Armed.or.Unarmed.=="Unarmed" & d$Ethnicity=="White", "WhiteUnarmed",Results)

#############################
table(d$ID,Results)

############################# This table is exported to the map data file