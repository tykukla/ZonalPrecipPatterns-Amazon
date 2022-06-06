# ----------------------------------------------------- # 
# Amazon Hydroclimate GCM/EBM data extraction           #
# -- FXN script to extract relevant model parameters    #
# ------------------                                    #
# T. Kukla (Stanford Univ. 2018)                        #
# ----------------------------------------------------- #



df_fun_multimodel <- function(ncin, modelname, model.idx, dname, latname, 
                              lonname, tsl.idx=NULL, fillNA=F, tsl.mean = F,
                              model.mean = F){
  
  # extract lat, lon, variable
  lat <- ncvar_get(ncin, latname, verbose = F)   # lat vector
  lon <- ncvar_get(ncin, lonname, verbose = F)   # lon vector
  model <- ncvar_get(ncin, modelname, verbose=F)     # model index
  dat <- ncvar_get(ncin, dname, verbose = F)     # data array
  # find indices 
  lat_dex <- which(dim(dat)==length(lat))
  lon_dex <- which(dim(dat)==length(lon))
  model_dex <- which(dim(dat)==length(model))
  
  # create a vector for latitude and longitude
  lon_vec <- rep(lon, length(lat))
  lat_vec <- rep(lat, each=length(lon))
  
  for(m in model.idx){
    for(t in tsl.idx){
      thesedat <- dat[ , , t, m]
      vec.temp <- as.vector(thesedat)
      df.temp <- as.data.table(cbind(lon_vec, lat_vec, vec.temp)) 
      colnames(df.temp) <- c('lon', 'lat', dname)
      df.temp$month <- t
      df.temp$model <- model[m]
      
      if(m==model.idx[1] & t==tsl.idx[1]){
        dfout <- df.temp
      } else{
        dfout <- as.data.table(rbind(dfout, df.temp))
      }
      
    }
  }
  
  if(tsl.mean==T & model.mean==T){
    dfout %>%
      group_by(lon, lat, model) %>%
      summarise_if(is.numeric, mean) %>%
      group_by(lon, lat) %>%
      summarise_if(is.numeric, mean) -> dfout.final
      
    dfout.final <- as.data.table(dfout.final)
  } else if(tsl.mean==F & model.mean==T){
    dfout %>%
      group_by(lon, lat, month) %>%
      summarise_if(is.numeric, mean) -> dfout.final
    
    dfout.final <- as.data.table(dfout.final)
  } else if(tsl.mean==T & model.mean==F){
    dfout%>%
      group_by(lon, lat, model) %>%
      summarise_if(is.numeric, mean) -> dfout.final
    
    dfout.final <- as.data.table(dfout.final)
  } else{dfout.final <- as.data.table(dfout)}
  
  #... fix lon to -180 to 180
  dfout.final$lon <- ifelse(dfout.final$lon > 180, dfout.final$lon - 360, dfout.final$lon)
  #... return result
  return(dfout.final)
}







# ------------- FUNCTION 1 ------------- #
# ------ the dataframe function -------- #
#     (mono-layer data files only)       #
# -------------------------------------- # 
## Function to build a 3-D array with the following dimensions:
#.. [longitude, latitude, month] 
## OR build a 2-D matrix with the following dimensions:
#.. [grid cell, parameter]
# -------------------------------------- #
# FUNCTION VARIABLES:
# 1) ncin --> the input ncdf file (already opened using the nc_open fxn)
# 2) dname --> character of the data variable's name (i.e. 'tas' for surface temperature)
# 3) latname --> character for latitude in the ncdf
# 4) lonname --> character for longitude in the ncdf
# 5) tsl --> single integer or vector denoting which timeslices over which to take the average
#    (NULL if you want every timeslice output OR if the ncdf only contains one timeslice)
# 6) fillNA --> whether to fill NA values in the dataframe (default is FALSE)
# 7) returnfile --> options: 'dataframe', 'array', 'monthDF'
#    (lets you decide if you want the function to return a 2-D matrix or a 3-D array, or a 2-D matrix with monthly data)
#    (the array maintains spatial structure of the ncdf file by outputting a grid)
# 8) maxlat ; maxlon ; minlat ; minlon --> specify a bounding box using these terms if you want a smaller output file
#    (format should be -90 to 90 for latitude, and -180 to 180 for longitude)
# 9) multilayer_data --> options: 'average' or 'surface'. If it's a multilayer dataset should the output be the average
#    of all layers or just the surface value

df_fun <- function(ncin, dname, latname, lonname, tsl=NULL, fillNA=F, returnfile='dataframe', 
                   maxlat=NULL, minlat=NULL, maxlon=NULL, minlon=NULL, multilayer_data='average'){
  
  # extract lat, lon, variable
  lat <- ncvar_get(ncin, latname, verbose = F)   # lat vector
  lon <- ncvar_get(ncin, lonname, verbose = F)   # lon vector
  # model <- ncvar_get(ncin, modelname, verbose=F)     # model index
  dat <- ncvar_get(ncin, dname, verbose = F)     # data array
  
  # # find out if there's a pressure coordinate
  # pressure_test <- tryCatch({ncvar_get(ncin, 'plev', verbose = F) }, error=function(e){})
  # if(is.numeric(pressure_test)){ # average over all pressures or select surface pressure
  #   
  # }
  
  #... is there a bounding box specified? 
  # SWITCH: 1= YES, obey the bounding box
  #         2= NO, return all data
  if(is.numeric(maxlat)){
    myswitch=1
  } else{myswitch=2}
  #... run the switch
  switch(myswitch,
         
         {#if myswitch = 1 (crop the dataset)   
           # make sure longitude is correct coordinates
           if(max(lon) > 180){
             if(maxlon < 0){maxlon <- maxlon + 360}   # convert to correct coordinates
             if(minlon < 0){minlon <- minlon + 360
             } else{}
           }
           # ==================== DATAFRAME ====================== #
           ## 1) first here's how we return a dataframe
           if(returnfile=='dataframe'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             
             #... check if we need to average across timeslices
             # ------ WE TAKE MEAN OF ALL TIMESLICES ------- #
             if(length(dim(dat)) > 2 & is.null(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=dim(dat)[t_dex], ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               for(i in 1:dim(dat)[length(dim(dat))]){
                 outmat[i, ] <- as.vector(dat[,,i])
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               # cut out just the area of interest 
               outdf %>%
                 filter(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon) -> returnMe
               return(returnMe)
               
               # ------ WE TAKE MEAN OF SOME TIMESLICES ------- #
             } else if(length(dim(dat)) > 2 & is.numeric(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=length(tsl), ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               j <- 1 # track the timeslice index
               for(i in tsl){
                 outmat[j, ] <- as.vector(dat[,,i])
                 j <- j + 1
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               # cut out just the area of interest 
               outdf %>%
                 filter(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon) -> returnMe
               return(returnMe)
               
               # ------ THERE ARE NO TIMESLICES ------- #
             } else if(length(dim(dat))==2){
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # build the dataframe 
               outdat <- as.vector(dat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               # cut out just the area of interest 
               outdf %>%
                 filter(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon) -> returnMe
               return(returnMe)
             } else{stop("Check that data dimensions are compatible with function")}
             
             # ==================== ARRAY ====================== #
           } else if(returnfile == 'array'){
             latdims <- which(lat <= maxlat & lat >= minlat)
             londims <- which(lon <= maxlon & lon >= minlon)
              if(lon_dex == 1){
                outar <- dat[londims, latdims, ]
              } else{ outar <- dat[latdims, londims, ] }
             #.. return the cropped data
             return(outar)
         # ==================== MONTH DATA FRAME ====================== #
         } else if(returnfile == 'monthDF'){
           #... first find out which dimensions go to what
           lat_dex <- which(dim(dat)==length(lat))
           lon_dex <- which(dim(dat)==length(lon))
           
           #... the IPSL model is 96x96 in which case we assume that the first dimension is longitude
           if(length(lat_dex) > 1 | length(lon_dex) > 1){
             lat_dex <- 2; lon_dex <- 1
           }
           
           # create a vector for latitude and longitude
           if(lon_dex==1){
             lon_vec <- rep(lon, length(lat))
             lat_vec <- rep(lat, each=length(lon))
           } else{
             lon_vec <- rep(lon, each=length(lat))
             lat_vec <- rep(lat, length(lon))}
           # loop through all timeslices and create a matrix to take averages
           outmat <- matrix(ncol=dim(dat)[length(dim(dat))], nrow= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
           for(i in 1:dim(dat)[length(dim(dat))]){
             outmat[ , i] <- as.vector(dat[,,i])
           }
           outmat <- as_tibble(outmat)
           colnames(outmat) <- paste('month_', c(1:12), sep='')
           outmat$lat <- lat_vec ; outmat$lon <- lon_vec
           outdf <- as_tibble( reshape2::melt(outmat, id.vars = c('lat', 'lon')) )
           colnames(outdf) <- c('lat', 'lon', 'month', dname)
           #... return the dataframe
           return(outdf)
           
         } else{stop("returnfile must equal 'array' or 'dataframe' ")} 
           
           },
         
         # --------------------------------------------------------------------- #
         {#if myswitch = 2 (return all data)  
           
           # ==================== DATAFRAME ====================== #
           ## 1) first here's how we return a dataframe
           if(returnfile=='dataframe'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             # model_dex <- which(dim(dat)==length(model))
             
             
             #... check if we need to average across timeslices
             # ------ WE TAKE MEAN OF ALL TIMESLICES ------- #
             if(length(dim(dat)) > 2 & is.null(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=dim(dat)[length(dim(dat))], ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               for(i in 1:dim(dat)[length(dim(dat))]){
                 outmat[i, ] <- as.vector(dat[,,i])
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               return(outdf)
               
               # ------ WE TAKE MEAN OF SOME TIMESLICES ------- #
             } else if(length(dim(dat)) > 2 & is.numeric(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=length(tsl), ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               j <- 1 # track the timeslice index
               for(i in tsl){
                 outmat[j, ] <- as.vector(dat[,,i])
                 j <- j + 1
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               return(outdf)
               
               # ------ THERE ARE NO TIMESLICES ------- #
             } else if(length(dim(dat))==2){
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # build the dataframe 
               outdat <- as.vector(dat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               return(outdf)
             } else{stop("Check that data dimensions are compatible with function")}
             
             # ==================== ARRAY ====================== #
           } else if(returnfile == 'array'){
             return(dat)
             
             
             # ==================== MONTH DATA FRAME ====================== #
           } else if(returnfile == 'monthDF'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             
             #... the IPSL model is 96x96 in which case we assume that the first dimension is longitude
             if(length(lat_dex) > 1 | length(lon_dex) > 1){
               lat_dex <- 2; lon_dex <- 1
             }
             
             # create a vector for latitude and longitude
             if(lon_dex==1){
               lon_vec <- rep(lon, length(lat))
               lat_vec <- rep(lat, each=length(lon))
             } else{
               lon_vec <- rep(lon, each=length(lat))
               lat_vec <- rep(lat, length(lon))}
             # loop through all timeslices and create a matrix to take averages
             outmat <- matrix(ncol=dim(dat)[length(dim(dat))], nrow= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
             for(i in 1:dim(dat)[length(dim(dat))]){
               outmat[ , i] <- as.vector(dat[,,i])
             }
             outmat <- as_tibble(outmat)
             colnames(outmat) <- paste('month_', c(1:12), sep='')
             outmat$lat <- lat_vec ; outmat$lon <- lon_vec
             outdf <- as_tibble( reshape2::melt(outmat, id.vars = c('lat', 'lon')) )
             colnames(outdf) <- c('lat', 'lon', 'month', dname)
             #... return the dataframe
             return(outdf)
             
           } else{stop('returnfile must equal "array" or "dataframe" or "monthDF" ' )}
         }
  )
}




# ==================================================================================================================== #
# ------------- FUNCTION 2 ------------- #
# ------ the dataframe function -------- #
#   (pressure level data files only)     #
# -------------------------------------- # 
## Function to build a 3-D array with the following dimensions:
#.. [longitude, latitude, pressure_level, month] 
## OR build a 2-D matrix with the following dimensions:
#.. [grid cell, parameter]
# -------------------------------------- #
## New inputs:
## plevelname --> the variable name for the pressure level (generally 'plev')
## pressure_level --> the numeric value of the single pressure level to extract (the index, specifically, not the number with units)

df_fun_SinglePlevel <- function(ncin, dname, latname, lonname, plevelname, tsl=NULL, fillNA=F, returnfile='dataframe', 
                   pressure_level, maxlat=NULL, minlat=NULL, maxlon=NULL, minlon=NULL, multilayer_data='average'){
  
  # extract lat, lon, variable
  lat <- ncvar_get(ncin, latname, verbose = F)   # lat vector
  lon <- ncvar_get(ncin, lonname, verbose = F)   # lon vector
  plev <- ncvar_get(ncin, plevelname, verbose = F) # pressure level vector
  plev_dex <- pressure_level
  dat <- ncvar_get(ncin, dname, verbose = F)     # data array
  
  #... is there a bounding box specified? 
  # SWITCH: 1= YES, obey the bounding box
  #         2= NO, return all data
  if(is.numeric(maxlat)){
    myswitch=1
  } else{myswitch=2}
  #... run the switch
  switch(myswitch,
         
         {#if myswitch = 1 (crop the dataset)   
           # make sure longitude is correct coordinates
           if(max(lon) > 180){
             if(maxlon < 0){maxlon <- maxlon + 360}   # convert to correct coordinates
             if(minlon < 0){minlon <- minlon + 360
             } else{}
           }
           # ==================== DATAFRAME ====================== #
           ## 1) first here's how we return a dataframe
           if(returnfile=='dataframe'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             
             #... check if we need to average across timeslices
             # ------ WE TAKE MEAN OF ALL TIMESLICES ------- #
             if(length(dim(dat)) > 2 & is.null(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=dim(dat)[4], ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               for(i in 1:dim(dat)[length(dim(dat))]){
                 outmat[i, ] <- as.vector(dat[,,plev_dex ,i])
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               # cut out just the area of interest 
               outdf %>%
                 filter(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon) -> returnMe
               return(returnMe)
               
               # ------ WE TAKE MEAN OF SOME TIMESLICES ------- #
             } else if(length(dim(dat)) > 2 & is.numeric(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=length(tsl), ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               j <- 1 # track the timeslice index
               for(i in tsl){
                 outmat[j, ] <- as.vector(dat[,,plev_dex ,i])
                 j <- j + 1
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               # cut out just the area of interest 
               outdf %>%
                 filter(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon) -> returnMe
               return(returnMe)
               
               # ------ THERE ARE NO TIMESLICES ------- #
             } else if(length(dim(dat))==2){
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # build the dataframe 
               outdat <- as.vector(dat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               # cut out just the area of interest 
               outdf %>%
                 filter(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon) -> returnMe
               return(returnMe)
             } else{stop("Check that data dimensions are compatible with function")}
             
             # ==================== ARRAY ====================== #
           } else if(returnfile == 'array'){
             latdims <- which(lat <= maxlat & lat >= minlat)
             londims <- which(lon <= maxlon & lon >= minlon)
             if(lon_dex == 1){
               outar <- dat[londims, latdims, ]
             } else{ outar <- dat[latdims, londims, ] }
             #.. return the cropped data
             return(outar)
             
             # ==================== MONTH DATA FRAME ====================== #
           } else if(returnfile == 'monthDF'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             
             #... the IPSL model is 96x96 in which case we assume that the first dimension is longitude
             if(length(lat_dex) > 1 | length(lon_dex) > 1){
               lat_dex <- 2; lon_dex <- 1
             }
             
             # create a vector for latitude and longitude
             if(lon_dex==1){
               lon_vec <- rep(lon, length(lat))
               lat_vec <- rep(lat, each=length(lon))
             } else{
               lon_vec <- rep(lon, each=length(lat))
               lat_vec <- rep(lat, length(lon))}
             # loop through all timeslices and create a matrix to take averages
             outmat <- matrix(ncol=dim(dat)[length(dim(dat))], nrow= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
             for(i in 1:dim(dat)[length(dim(dat))]){
               outmat[ , i] <- as.vector(dat[,,plev_dex ,i])
             }
             outmat <- as_tibble(outmat)
             colnames(outmat) <- paste('month_', c(1:12), sep='')
             outmat$lat <- lat_vec ; outmat$lon <- lon_vec
             outdf <- as_tibble( reshape2::melt(outmat, id.vars = c('lat', 'lon')) )
             colnames(outdf) <- c('lat', 'lon', 'month', dname)
             #... return the dataframe
             return(outdf)
             
           } else{stop("returnfile must equal 'array' or 'dataframe' ")}
           
         },
         
         # --------------------------------------------------------------------- #
         {#if myswitch = 2 (return all data)  
           
           # ==================== DATAFRAME ====================== #
           ## 1) first here's how we return a dataframe
           if(returnfile=='dataframe'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             
             #... the IPSL model is 96x96 in which case we assume that the first dimension is longitude
             if(length(lat_dex) > 1 | length(lon_dex) > 1){
               lat_dex <- 2; lon_dex <- 1
             }
             
             #... check if we need to average across timeslices
             # ------ WE TAKE MEAN OF ALL TIMESLICES ------- #
             if(length(dim(dat)) > 2 & is.null(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=dim(dat)[length(dim(dat))], ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               for(i in 1:dim(dat)[length(dim(dat))]){
                 outmat[i, ] <- as.vector(dat[,,plev_dex ,i])
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               return(outdf)
               
               # ------ WE TAKE MEAN OF SOME TIMESLICES ------- #
             } else if(length(dim(dat)) > 2 & is.numeric(tsl)){
               # assume the time dimension is last
               if(lat_dex == length(dim(dat)) | lon_dex == length(dim(dat))){
                 stop("The time dimension is not last in the array... structure inconsistent with function!!")
               }
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # loop through all timeslices and create a matrix to take averages
               outmat <- matrix(nrow=length(tsl), ncol= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
               j <- 1 # track the timeslice index
               for(i in tsl){
                 outmat[j, ] <- as.vector(dat[,,plev_dex ,i])
                 j <- j + 1
               }
               # take the mean of timeslices and create dataframe
               outdat <- colMeans(outmat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               return(outdf)
               
               # ------ THERE ARE NO TIMESLICES ------- #
             } else if(length(dim(dat))==2){
               # create a vector for latitude and longitude
               if(lon_dex==1){
                 lon_vec <- rep(lon, length(lat))
                 lat_vec <- rep(lat, each=length(lon))
               } else{
                 lon_vec <- rep(lon, each=length(lat))
                 lat_vec <- rep(lat, length(lon))}
               # build the dataframe 
               outdat <- as.vector(dat)
               outdf <- as_tibble(cbind(lat_vec, lon_vec, outdat))
               colnames(outdf) <- c('lat', 'lon', dname)
               return(outdf)
             } else{stop("Check that data dimensions are compatible with function")}
             
             # ==================== ARRAY ====================== #
           } else if(returnfile == 'array'){
             return(dat)
             
             
             # ==================== MONTH DATA FRAME ====================== #
           } else if(returnfile == 'monthDF'){
             #... first find out which dimensions go to what
             lat_dex <- which(dim(dat)==length(lat))
             lon_dex <- which(dim(dat)==length(lon))
             
             #... the IPSL model is 96x96 in which case we assume that the first dimension is longitude
             if(length(lat_dex) > 1 | length(lon_dex) > 1){
               lat_dex <- 2; lon_dex <- 1
             }
             
             # create a vector for latitude and longitude
             if(lon_dex==1){
               lon_vec <- rep(lon, length(lat))
               lat_vec <- rep(lat, each=length(lon))
             } else{
               lon_vec <- rep(lon, each=length(lat))
               lat_vec <- rep(lat, length(lon))}
             # loop through all timeslices and create a matrix to take averages
             outmat <- matrix(ncol=dim(dat)[length(dim(dat))], nrow= (dim(dat)[lat_dex]*dim(dat)[lon_dex]) )
             for(i in 1:dim(dat)[length(dim(dat))]){
               outmat[ , i] <- as.vector(dat[,,plev_dex ,i])
             }
             outmat <- as_tibble(outmat)
             colnames(outmat) <- paste('month_', c(1:12), sep='')
             outmat$lat <- lat_vec ; outmat$lon <- lon_vec
             outdf <- as_tibble( reshape2::melt(outmat, id.vars = c('lat', 'lon')) )
             colnames(outdf) <- c('lat', 'lon', 'month', dname)
             #... return the dataframe
             return(outdf)
             
           } else{stop('returnfile must equal "array" or "dataframe" or "monthDF" ')}
         }
  )
}





# ==================================================================================================================== #

# ------------- FUNCTION 3 ------------- #
# ------ the wind calc function -------- #
#     (mono-layer data files only)       #
# -------------------------------------- # 
## Function to calculate wind direction and magnitude from xyz data
# -------------------------------------- #
# FUNCTION VARIABLES:
# 1) Uwind --> vector of eastward wind velocity data
# 2) Vwind --> vector of northward wind velocity data

Wind_calculator <- function(Uwind, Vwind){
  #... pull out the U and the V
  myU <- as_tibble(Uwind) ; colnames(myU) <- 'uas'
  myV <- as_tibble(Vwind) ; colnames(myV) <- 'vas'
  
  # ************* WIND DIRECTION CALCULATION ************* #
  #... first extract all of the possible directions for the calculation
  d1_dex <- which(myU > 0)                   # X > 0
  d2_dex <- which(myV >= 0 & myU < 0)        # Y >= 0, X < 0
  d3_dex <- which(myV < 0 & myU < 0)         # Y < 0 , X < 0
  d4_dex <- which(myV > 0 & myU == 0)        # Y > 0 , X = 0
  d5_dex <- which(myV < 0 & myU == 0)        # Y < 0 , X = 0
  d6_dex <- which(myV == 0 & myU == 0)       # Y = 0 , X = 0
  d7_dex <- which(is.na(myV) & is.na(myU))   # in case there is no data (still preserve dimensions)
  
  #... create an output vector and fill it up
  # ***************** #
  dir_vec <- vector() #
  # ***************** #
  dir_vec[d1_dex] <- rad2deg( as.numeric( atan(myV$vas[d1_dex] / myU$uas[d1_dex] )))
  dir_vec[d2_dex] <- rad2deg( as.numeric( atan(myV$vas[d2_dex] / myU$uas[d2_dex]) + pi) )
  dir_vec[d3_dex] <- rad2deg( as.numeric( atan(myV$vas[d3_dex] / myU$uas[d3_dex]) - pi ) )
  dir_vec[d4_dex] <- rad2deg( as.numeric( pi / 2 ) )                       
  dir_vec[d5_dex] <- rad2deg( as.numeric( -pi / 2 ) )
  dir_vec[d6_dex] <- rad2deg( as.numeric( -pi / 2 ) )
  dir_vec[d7_dex] <- NA
  
  
  # **************** WIND SPEED CALCULATION ************** #
  #... Simple application of pythagorean theorem to get magnitude of wind speed 
  speed_vec <- sqrt(myU$uas**2 + myV$vas**2)
  
  
  # **************** ADD NEW COLS TO INPUT AND RETURN ************** #
  outlist <- list()
  outlist[[1]] <- dir_vec
  outlist[[2]] <- speed_vec
  names(outlist) <- c('windDir', 'windVel')
  
  # ----
  return(outlist)
}




