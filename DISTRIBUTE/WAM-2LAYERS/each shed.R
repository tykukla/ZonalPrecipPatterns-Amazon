###### calculation ######

# tracking E
etrack_extract <- function(){
  
  # track for each month, and all months together
  etrack01 <- etrack02 <- etrack03 <- etrack04 <- etrack05 <- etrack06 <- 
    etrack07 <- etrack08 <- etrack09 <- etrack10 <- etrack11 <- etrack12 <- 
    etrack_year <- array(0,dim=c(107,240))
  
  # total track for each month
  for (i in 1:30){
    etrack01 <- etrack01 + Etrack_monthly[i,1,,]
    etrack02 <- etrack02 + Etrack_monthly[i,2,,]
    etrack03 <- etrack03 + Etrack_monthly[i,3,,]
    etrack04 <- etrack04 + Etrack_monthly[i,4,,]
    etrack05 <- etrack05 + Etrack_monthly[i,5,,]
    etrack06 <- etrack06 + Etrack_monthly[i,6,,]
    etrack07 <- etrack07 + Etrack_monthly[i,7,,]
    etrack08 <- etrack08 + Etrack_monthly[i,8,,]
    etrack09 <- etrack09 + Etrack_monthly[i,9,,]
    etrack10 <- etrack10 + Etrack_monthly[i,10,,]
    etrack11 <- etrack11 + Etrack_monthly[i,11,,]
    etrack12 <- etrack12 + Etrack_monthly[i,12,,]
  }
  etrack_year <- etrack01 + etrack02+ etrack03 + etrack04 + etrack05 + etrack06 + etrack07 + etrack08 + etrack09 + etrack10 + etrack11 + etrack12
  
  # modify longitudes
  etrack <- array(0,c(13,107,240)) # monthly data plus yearly data, lat (79.5 ~ -79.5), lon (0 ~ 358.5)
  etrack[1,,] <- cbind(etrack01[,121:240],etrack01[,1:120]) # lon (180 ~ 358.5 or -180 ~ -2.5, 0 ~ 178.5)
  etrack[2,,] <- cbind(etrack02[,121:240],etrack02[,1:120])
  etrack[3,,] <- cbind(etrack03[,121:240],etrack03[,1:120])
  etrack[4,,] <- cbind(etrack04[,121:240],etrack04[,1:120])
  etrack[5,,] <- cbind(etrack05[,121:240],etrack05[,1:120])
  etrack[6,,] <- cbind(etrack06[,121:240],etrack06[,1:120])
  etrack[7,,] <- cbind(etrack07[,121:240],etrack07[,1:120])
  etrack[8,,] <- cbind(etrack08[,121:240],etrack08[,1:120])
  etrack[9,,] <- cbind(etrack09[,121:240],etrack09[,1:120])
  etrack[10,,] <- cbind(etrack10[,121:240],etrack10[,1:120])
  etrack[11,,] <- cbind(etrack11[,121:240],etrack11[,1:120])
  etrack[12,,] <- cbind(etrack12[,121:240],etrack12[,1:120])
  etrack[13,,] <- cbind(etrack_year[,121:240],etrack_year[,1:120])

  return(etrack)
}

# tracking P
ptrack_extract <- function(){
  
  # track for each month, and all months together
  ptrack01 <- ptrack02 <- ptrack03 <- ptrack04 <- ptrack05 <- ptrack06 <- 
    ptrack07 <- ptrack08 <- ptrack09 <- ptrack10 <- ptrack11 <- ptrack12 <- 
    ptrack_year <- array(0,dim=c(107,240))
  
  # total track for each month
  for (i in 1:30){
    ptrack01 <- ptrack01 + Ptrack_monthly[i,1,,]
    ptrack02 <- ptrack02 + Ptrack_monthly[i,2,,]
    ptrack03 <- ptrack03 + Ptrack_monthly[i,3,,]
    ptrack04 <- ptrack04 + Ptrack_monthly[i,4,,]
    ptrack05 <- ptrack05 + Ptrack_monthly[i,5,,]
    ptrack06 <- ptrack06 + Ptrack_monthly[i,6,,]
    ptrack07 <- ptrack07 + Ptrack_monthly[i,7,,]
    ptrack08 <- ptrack08 + Ptrack_monthly[i,8,,]
    ptrack09 <- ptrack09 + Ptrack_monthly[i,9,,]
    ptrack10 <- ptrack10 + Ptrack_monthly[i,10,,]
    ptrack11 <- ptrack11 + Ptrack_monthly[i,11,,]
    ptrack12 <- ptrack12 + Ptrack_monthly[i,12,,]
  }
  ptrack_year <- ptrack01 + ptrack02+ ptrack03 + ptrack04 + ptrack05 + ptrack06 + ptrack07 + ptrack08 + ptrack09 + ptrack10 + ptrack11 + ptrack12
  
  # modify longitudes
  ptrack <- array(0,c(13,107,240)) # monthly data plus yearly data, lat (79.5 ~ -79.5), lon (0 ~ 358.5)
  ptrack[1,,] <- cbind(ptrack01[,121:240],ptrack01[,1:120]) # lon (180 ~ 358.5 or -180 ~ -2.5, 0 ~ 178.5)
  ptrack[2,,] <- cbind(ptrack02[,121:240],ptrack02[,1:120])
  ptrack[3,,] <- cbind(ptrack03[,121:240],ptrack03[,1:120])
  ptrack[4,,] <- cbind(ptrack04[,121:240],ptrack04[,1:120])
  ptrack[5,,] <- cbind(ptrack05[,121:240],ptrack05[,1:120])
  ptrack[6,,] <- cbind(ptrack06[,121:240],ptrack06[,1:120])
  ptrack[7,,] <- cbind(ptrack07[,121:240],ptrack07[,1:120])
  ptrack[8,,] <- cbind(ptrack08[,121:240],ptrack08[,1:120])
  ptrack[9,,] <- cbind(ptrack09[,121:240],ptrack09[,1:120])
  ptrack[10,,] <- cbind(ptrack10[,121:240],ptrack10[,1:120])
  ptrack[11,,] <- cbind(ptrack11[,121:240],ptrack11[,1:120])
  ptrack[12,,] <- cbind(ptrack12[,121:240],ptrack12[,1:120])
  ptrack[13,,] <- cbind(ptrack_year[,121:240],ptrack_year[,1:120])
  
  return(ptrack)
}

# tracked E and P for each site
load(paste(my.dir, "DTP_Etrack_monthly.RData", sep='/'))
etrack_DTP <- etrack_extract()
load(paste(my.dir, "PAR_Etrack_monthly.RData", sep='/'))
etrack_PAR <- etrack_extract()
load(paste(my.dir, "RGN_Etrack_monthly.RData", sep='/'))
etrack_RGN <- etrack_extract()
load(paste(my.dir, "DTP_Ptrack_monthly.RData", sep='/'))
ptrack_DTP <- ptrack_extract()
load(paste(my.dir, "PAR_Ptrack_monthly.RData", sep='/'))
ptrack_PAR <- ptrack_extract()
load(paste(my.dir, "RGN_Ptrack_monthly.RData", sep='/'))
ptrack_RGN <- ptrack_extract()

# precipitationshed based on tracked E
etrack_pshed <- function(etrack,cutoff){
  shed <- array(0,c(13,107,240))
  for (i in 1:13){ # 12 monthly + yearly
    for (n in 1:length(etrack[i,,])){
      if ((sum(etrack[i,,][order(etrack[i,,],decreasing=T)[1:n]])/sum(etrack[i,,])) > cutoff){
        shed[i,,][order(etrack[i,,],decreasing=T)[1:n]] <- 1 # for contours
        break
      }
    }
  }
  return(shed)
}

# evaporationshed based on tracked P
ptrack_eshed <- function(ptrack,cutoff){
  shed <- array(0,c(13,107,240))
  for (i in 1:13){ # 12 monthly + yearly
    for (n in 1:length(ptrack[i,,])){
      if ((sum(ptrack[i,,][order(ptrack[i,,],decreasing=T)[1:n]])/sum(ptrack[i,,])) > cutoff){
        shed[i,,][order(ptrack[i,,],decreasing=T)[1:n]] <- 1 # for contours
        break
      }
    }
  }
  return(shed)
}

# Esheds and Psheds, 0.7 threshold, calculation
pshed_DTP <- etrack_pshed(etrack_DTP,0.7)
eshed_DTP <- ptrack_eshed(ptrack_DTP,0.7)

pshed_PAR <- etrack_pshed(etrack_PAR,0.7)
eshed_PAR <- ptrack_eshed(ptrack_PAR,0.7)

pshed_RGN <- etrack_pshed(etrack_RGN,0.7)
eshed_RGN <- ptrack_eshed(ptrack_RGN,0.7)