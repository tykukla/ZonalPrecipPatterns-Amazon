# --------------------------------- #
# Pull out the EFE / EFPM analysis  #
# from the Boos & Korty model       #
# ---                               #
# T Kukla (Stanford Univ. 2021)     #
#                                   #
# Date created: Sept 1, 2021        #
# Last modified: May 11, 2022       #
# --------------------------------- #
library(ncdf4)
library(ncdf4.helpers)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(data.table)
# library(metR)
library(dplyr)
library(plyr)
# library(RColorBrewer)
# library(scico)
library(raster)
library(cmocean)
# library(oce)
# library(paletteer)

rm(list=ls())
setwd('C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data')

# set the base directory (where output and figures are)
basedir <- 'DISTRIBUTE/ENERGY_BALANCE/MID-HOLOCENE'
fxndir <- 'DISTRIBUTE/ENERGY_BALANCE/R-analysis'

## ... and the functions 
source(paste(fxndir, 'ncdfExtract_fxns_multimodel.R', sep='/'))
source(paste(fxndir, 'Intersect_Calculator.R', sep='/'))

## ... set data file names
datDir <- paste(basedir, "Output/", sep='/')
filens <- list.files(path = datDir, pattern = "NDJFMAM_")
force_vec <- vector()
for(j in c(1:length(filens))){
  force_vec[j] <- as.numeric(strsplit(filens, split="-")[[j]][2])
}

myVars <- c('EFPot' = 'Chi', 'div_EFPot' = 'div_Chi', 'MSE_source' = 'MSE_source',
            'MSE_udiv' = 'MSE_udiv', 'MSE_vdiv' = 'MSE_vdiv', 'u_gradX' = 'u_gradX',
            'u_gradY' = 'u_gradY', 'v_gradX' = 'v_gradX', 'v_gradY' = 'v_gradY')
myPLev <- rep('N', length=length(myVars))

## ... read in files
#... set up
mytsl <- NA                # set up the timeslices you want
tslName <- 'NDJFMAM'
latname <- 'lat'             # character vector for latitude in ncdf files
lonname <- 'lon'             # character vector for longitude in ncdf files
modelname <- 'model'         # character vector for model in ncdf files
# what to average out
time_averaged <- T       # whether to take mean of months

#... define a bounding box if desired
maxlat <- NULL
minlat <- NULL
maxlon <- NULL
minlon <- NULL

#... list to hold output data
outlist <- list()
for(i in 1:length(filens)){
  #... track progress
  print(paste("Now extracting:", i, 'of', length(filens), sep=' '))
  #... set file path
  myNC <- paste(datDir, filens[i], sep='')
  thisExperiment <- paste("Force-Wm2", force_vec[i], sep='_')
  
  #... open the netCDF
  ncin <- nc_open(myNC)
  
  #... loop through variables
  for(j in 1:length(myVars)){
    print(paste("Pulling out var:", j, "of", length(myVars), sep=' '))
    thisVar <- as.character(myVars[j])
    
    dfx <- as.data.table(df_fun(ncin=ncin, dname=thisVar, latname=latname, lonname = lonname))
    
    
    if(j==1){
      outdf <- dfx
    } else{
      outdf[, thisVar] <- dfx[, ..thisVar]
    }
    
  }
  
  # shut down nc file
  nc_close(ncin) 
  
  # correct longitude
  outdf$lon <- ifelse(outdf$lon > 180, outdf$lon - 360, outdf$lon)
  # add column names and month and experiment code
  outdf$Experiment <- thisExperiment
  outdf$months_avgd <- tslName
  outdf$Force_Wm2 <- force_vec[i]
  # save in a list
  outlist[[i]] <- outdf
}

#... bring it all together
names(outlist) <- paste("Force-Wm2", force_vec, sep='_')

# ------------------------------------------------------------------------------- #
## NOW PULL OUT THE INTERSECTS
for(j in 1:length(outlist)){
  tempdf <- outlist[[j]]
  thisForce <- tempdf$Force_Wm2[1]
  thisOutput <- efe_efpm.fun(df=tempdf, ensemble=T, returnContours = T, Chi_Filter=T,
                             myLims = c(min.x = -75, max.x = -15, min.y = -30, max.y = 20))
  thisInt <- thisOutput[[1]]
  thisContour <- thisOutput[[2]]
  # add forcing
  thisInt$force_Wm2 <- thisForce
  thisContour$force_Wm2 <- thisForce
  
  if(j==1){
    dfx <- thisInt
    dfcontour <- thisContour
  } else{
    dfx <- as.data.table(rbind(dfx, thisInt))
    dfcontour <- as.data.table(rbind(dfcontour, thisContour))
  }
}
# order 
dfx <- dfx[order(dfx$force_Wm2),]
dfcontour <- dfcontour[order(dfcontour$force_Wm2), ]
dfx.p <- na.omit(dfx)

# put it together to one saveable list
saveList <- list(dfcontour, dfx.p)
names(saveList) <- c("contours", "intersects")



# -------------------------------------- #
# SAVE THE RESULT ---------------------- #
# saveRDS(saveList, file=paste(basedir, "MH-GrSaharaSmallBand-contours-pts.RDS", sep='/'))
# saveRDS(outlist, file=paste(basedir, 'MH-GrSahara_ALL.RDS', sep='/'))
# -------------------------------------- #




# ------------------------------------------------------------------------------- #
