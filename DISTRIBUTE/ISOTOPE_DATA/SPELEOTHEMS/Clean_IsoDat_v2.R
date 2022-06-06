# *************************************************************** # 
# Script to clean / smooth the Amazon transect isotope data       #
# -------                                                         #
# T. Kukla (Stanford Univ. 2019)                                  #
#                                                                 #
# Date created: May 11, 2019                                      #
# Last modified: June 1, 2022                                     #
# *************************************************************** #
library(dplyr)
library(ggplot2)
library(ggmap)
library(ggpubr)
# library(gganimate)
library(tidyquant)
library(viridis)
library(data.table)

rm(list=ls())

# directory where DISTRIBUTE is located
user.dir <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data'
# shouldn't need to change
mainPath <- 'DISTRIBUTE/ISOTOPE_DATA/SPELEOTHEMS'
setwd(paste(user.dir, mainPath, sep='/'))


# -----------------------------------------------------------------
# Read in the data 
dfParaiso <- as.data.table(read.csv('Paraiso.csv')) ; colnames(dfParaiso) <- c("year_BP", "d18O", "d13C", "ref")
dfDiamante <- as.data.table(read.csv('Diamante.csv')) ; colnames(dfDiamante)[1] <- "Speleothem_ID"
dfRN <- as.data.table(read.csv('RioGrandeDoNorte.csv')) ; colnames(dfRN)[1] <- "Speleothem_ID"
dfTigre <- as.data.table(read.csv('TigrePerdido.csv'))

# -----------------------------------------------------------------
# clean up the data 
#... the central record
dfCent <- dfParaiso[ , ':='(ref = NULL)]
dfCent <- dfCent[order(dfCent$year_BP), ]
#... add NA's
gapThreshold <- -1e3
dfCent$gap <- c(NA, dfCent[1:(length(dfCent$year_BP)-1), ]$year_BP - dfCent[2:length(dfCent$year_BP), ]$year_BP)
dfCent$BeyondGap <- dfCent$gap < gapThreshold
gapDex <- which(dfCent$BeyondGap==T)
for(i in 1:length(gapDex)){
  #... collect the years 
  oldYear <- dfCent[gapDex[i], ]$year_BP
  youngYear <- dfCent[gapDex[i]-1, ]$year_BP
  theNAyears <- seq(youngYear+1, oldYear-1, by=10)
  theAddition <- as.data.table(cbind(theNAyears, -9999, -9999, -9999, -9999))
  colnames(theAddition) <- c('year_BP', 'd18O', 'd13C', 'gap', 'BeyondGap')
  
  #... add to the data
  if(i==1){
    dfCent_x <- rbind(dfCent, theAddition)
  } else{
    dfCent_x <- rbind(dfCent_x, theAddition)
  }
}

#... to deal with the NAs we first interpolate and then do the moving average
myAges <- seq(0,50e3, length=5e3)
CentRecord <- approxfun(x=dfCent_x$year_BP, y=dfCent_x$d18O, rule=1)
CentDF <- as.data.table(cbind(myAges, CentRecord(myAges))) ; colnames(CentDF) <- c('age', 'd18')
CentDF$d18[which(CentDF$d18 < -100)] <- NA
dfCent <- CentDF
dfCent$site <- 'central'
dfCent$age <- dfCent$age / 1e3
dfCent <- dfCent[order(dfCent$age, decreasing=TRUE), ]

#... the eastern record
dfEast <- dfRN[ , ':='(Ref = NULL)]
dfEast <- dfEast[order(dfEast$age_ky, decreasing=TRUE), ]
#... add NA's
gapThreshold <- 0.5
dfEast$gap <- c(NA, dfEast[1:(length(dfEast$age_ky)-1), ]$age_ky - dfEast[2:length(dfEast$age_ky), ]$age_ky)
dfEast$BeyondGap <- abs(dfEast$gap) > gapThreshold
gapDex <- which(dfEast$BeyondGap==T)
for(i in 1:length(gapDex)){
  #... collect the years 
  oldYear <- dfEast[gapDex[i], ]$age_ky
  youngYear <- dfEast[(gapDex[i]-1), ]$age_ky
  theNAyears <- seq(youngYear+.001, oldYear-.001, by=-.01)
  theAddition <- as.data.table(cbind(NA, theNAyears, -9999, -9999, -9999))
  colnames(theAddition) <- c('Speleothem_ID', 'age_ky', 'd18O', 'gap', 'BeyondGap')
  
  #... add to the data
  if(i==1){
    dfEast_x <- rbind(dfEast, theAddition)
  } else{
    dfEast_x <- rbind(dfEast_x, theAddition)
  }
}
myAges <- seq(0,50e3, length=5e3)
EastRecord <- approxfun(x=dfEast_x$age_ky*1e3, y=dfEast_x$d18O, rule=1)
EastDF <- as.data.table(cbind(myAges, EastRecord(myAges))) ; colnames(EastDF) <- c('age','d18')
EastDF$d18[which(EastDF$d18 < -100)] <- NA
dfEast <- EastDF
dfEast$age <- dfEast$age / 1e3
dfEast$site <- 'east'
dfEast <- dfEast[order(dfEast$age, decreasing=TRUE), ]

#... the western record
dfTigre <- dfTigre[order(dfTigre$age), ]   # sort the matrix
dfDiamante <- dfDiamante[order(dfDiamante$Age_kyBP), ]
westAge <- c((dfTigre$age/1e3), dfDiamante$Age_kyBP)
westd18 <- c(dfTigre$d18, dfDiamante$d18O)
westd13 <- c(dfTigre$d13, dfDiamante$d13C)
dfWest <- as.data.table(cbind(westAge, westd18, westd13)) ; colnames(dfWest) <- c('age', 'd18', 'd13')
dfWest$d18_corrected <- dfWest$d18 - 1.4     # correction of Wang et al., 2017
dfWest <- dfWest[order(dfWest$age, decreasing=TRUE), ]
#... add NA's
gapThreshold <- 1
dfWest$gap <- c(NA, dfWest[1:(length(dfWest$age)-1), ]$age - dfWest[2:length(dfWest$age), ]$age)
dfWest$BeyondGap <- abs(dfWest$gap) > gapThreshold
gapDex <- which(dfWest$BeyondGap==T)
for(i in 1:length(gapDex)){
  #... collect the years 
  oldYear <- dfWest[gapDex[i], ]$age
  youngYear <- dfWest[(gapDex[i]-1), ]$age
  theNAyears <- seq(youngYear+.001, oldYear-.001, by=-.01)
  theAddition <- as.data.table(cbind(theNAyears, -9999, -9999, -9999, NA, NA))
  colnames(theAddition) <- c('age', 'd18', 'd13', 'd18_corrected', 'gap', 'BeyondGap')
  
  #... add to the data
  if(i==1){
    dfWest_x <- rbind(dfWest, theAddition)
  } else{
    dfWest_x <- rbind(dfWest_x, theAddition)
  }
}
myAges <- seq(0,50e3, length=5e3)
WestRecord <- approxfun(x=dfWest_x$age*1e3, y=dfWest_x$d18_corrected, rule=1)
WestDF <- as.data.table(cbind(myAges, WestRecord(myAges))) ; colnames(WestDF) <- c('age','d18')
WestDF$d18[which(WestDF$d18 < -100)] <- NA
dfWest <- WestDF
dfWest$site <- 'west'
dfWest$age <- dfWest$age / 1e3
dfWest <- dfWest[order(dfWest$age, decreasing=TRUE), ]



# remove everything older than what we care about 
dfCent <- dfCent[which(dfCent$age <= 50), ]
dfEast <- dfEast[which(dfEast$age <= 50), ]
dfWest <- dfWest[which(dfWest$age <= 50), ]

# calculate differences
westDist <- 2.4                                # thousands of kilometers
eastDist <- 1.5                                # thousands of kilometers
DdEast <- (dfCent$d18 - dfEast$d18) / eastDist
DdWest <- (dfWest$d18 - dfCent$d18) / westDist

# bring together 
dfinterp <- as_tibble(cbind(dfEast$age, dfEast$d18, dfCent$d18, dfWest$d18, DdEast, DdWest)) 
colnames(dfinterp) <- c('age_ky', 'd18east', 'd18cent', 'd18west', 'Dd_East','Dd_West')

# smooth the isotope data 
mybdwth <- 0.5
dfinterp$d18east_s <- rev(ksmooth(x=dfinterp$age_ky, y=dfinterp$d18east, bandwidth=mybdwth)$y)
dfinterp$d18west_s <- rev(ksmooth(x=dfinterp$age_ky, y=dfinterp$d18west, bandwidth=mybdwth)$y)
dfinterp$d18cent_s <- rev(ksmooth(x=dfinterp$age_ky, y=dfinterp$d18cent, bandwidth=mybdwth)$y)

dfinterp$Dd_East_s <- (dfinterp$d18cent_s - dfinterp$d18east_s) / eastDist
dfinterp$Dd_West_s <- (dfinterp$d18west_s - dfinterp$d18cent_s) / westDist


#... save the resulting data file 
# saveRDS(dfinterp, paste(mainPath, 'IsoDF_cleaned.RDS', sep='/'))
