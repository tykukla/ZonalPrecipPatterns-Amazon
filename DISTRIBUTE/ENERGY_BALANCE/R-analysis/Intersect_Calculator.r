

# ---- FUNCTION to find EFE / EFPM intersection
# Chi_Filter only looks for EFE/EFPMs when Chi is less than zero
# myLims refers to lat (y) and long (x)
efe_efpm.fun <- function(df, myLims = c('min.x'=-75, 'max.x'=-30, 'min.y'=-30, 'max.y'=20), ensemble = T, returnContours = F, Chi_Filter = T ){
  if(ensemble == T){
    if(Chi_Filter==T){
      # ... get the zero contour in the x direction
      m.x <- ggplot() + 
        geom_contour(data=df[which(Chi < 0)], aes(x=lon, y=lat, z=MSE_vdiv), color='blue', breaks=0) 
      m.y <- ggplot() + 
        geom_contour(data=df[which(Chi < 0)], aes(x=lon, y=lat, z=MSE_udiv), color='blue', breaks=0) 
    } else{
      # ... get the zero contour in the x direction
      m.x <- ggplot() + 
        geom_contour(data=df, aes(x=lon, y=lat, z=MSE_vdiv), color='blue', breaks=0) 
      m.y <- ggplot() + 
        geom_contour(data=df, aes(x=lon, y=lat, z=MSE_udiv), color='blue', breaks=0) 
    }
    
    # ... crop to window 
    mycontx <- ggplot_build(m.x)$data[[1]]
    myconty <- ggplot_build(m.y)$data[[1]]
    mycont2x <- mycontx[which(mycontx$x > myLims['min.x'] & mycontx$x < myLims['max.x'] & mycontx$y > myLims['min.y'] & mycontx$y < myLims['max.y']),]
    mycont2y <- myconty[which(myconty$x > myLims['min.x'] & myconty$x < myLims['max.x'] & myconty$y > myLims['min.y'] & myconty$y < myLims['max.y']),]
    # get just xy points
    mycont2x.s <- as.data.table(cbind(mycont2x$x, mycont2x$y))
    mycont2y.s <- as.data.table(cbind(mycont2y$x, mycont2y$y))
    # find the distance between all points 
    if(nrow(mycont2y.s) == 0 | nrow(mycont2x.s) == 0){
      outpoint <- as.data.table(cbind(NA, NA)) ; colnames(outpoint) <- c('x', 'y')
      if(returnContours == T){
        colnames(mycont2x.s) <- c('lon', 'lat')
        mycont2x.s$line <- 'EFE'
        colnames(mycont2y.s) <- c('lon', 'lat')
        mycont2y.s$line <- 'EFPM'
        out_contour <- as.data.table(rbind(mycont2x.s, mycont2y.s))
        return(list(outpoint, out_contour))
      } else{
        return(outpoint)
      }
      
    } else{
      tx <- pointDistance(mycont2x.s, mycont2y.s, lonlat=T, allpairs=T)
      my.min <- which(tx==min(tx), arr.ind=TRUE)
      this.lon <- mean(c(mycont2x.s[my.min[1,][1]]$V1, mycont2y.s[my.min[1,][2]]$V1))
      this.lat <- mean(c(mycont2x.s[my.min[1,][1]]$V2, mycont2y.s[my.min[1,][2]]$V2))
      
      # build output df
      outpoint <- as.data.table(cbind(this.lon, this.lat)) ; colnames(outpoint) <- c('x', 'y')
      
      if(returnContours == T){
        colnames(mycont2x.s) <- c('lon', 'lat')
        mycont2x.s$line <- 'EFE'
        colnames(mycont2y.s) <- c('lon', 'lat')
        mycont2y.s$line <- 'EFPM'
        out_contour <- as.data.table(rbind(mycont2x.s, mycont2y.s))
        return(list(outpoint, out_contour))
      } else{
        return(outpoint)
        }
    }
    
  } else{
    theseModels <- unique(df$model)
    # ... loop by model 
    for(i in 1:length(theseModels)){
      if(Chi_Filter==T){
        # ... get the zero contour in the x direction
        m.x <- ggplot() + 
          geom_contour(data=df[which(df$Chi < 0 & df$model==theseModels[i])], aes(x=lon, y=lat, z=MSE_vdiv), color='blue', breaks=0) 
        m.y <- ggplot() + 
          geom_contour(data=df[which(df$Chi < 0 & df$model==theseModels[i])], aes(x=lon, y=lat, z=MSE_udiv), color='blue', breaks=0) 
      } else{
        # ... get the zero contour in the x direction
        m.x <- ggplot() + 
          geom_contour(data=df[which(df$model==theseModels[i])], aes(x=lon, y=lat, z=MSE_vdiv), color='blue', breaks=0) 
        m.y <- ggplot() + 
          geom_contour(data=df[which(df$model==theseModels[i])], aes(x=lon, y=lat, z=MSE_udiv), color='blue', breaks=0) 
      }
      
      # ... crop to window 
      mycontx <- ggplot_build(m.x)$data[[1]]
      myconty <- ggplot_build(m.y)$data[[1]]
      mycont2x <- mycontx[which(mycontx$x > myLims['min.x'] & mycontx$x < myLims['max.x'] & mycontx$y > myLims['min.y'] & mycontx$y < myLims['max.y']),]
      mycont2y <- myconty[which(myconty$x > myLims['min.x'] & myconty$x < myLims['max.x'] & myconty$y > myLims['min.y'] & myconty$y < myLims['max.y']),]
      # get just xy points
      mycont2x.s <- as.data.table(cbind(mycont2x$x, mycont2x$y))
      mycont2y.s <- as.data.table(cbind(mycont2y$x, mycont2y$y))
      if(length(mycont2x.s)==0){
        this.lon <- NA
        this.lat <- NA
      } else{
        # find the distance between all points 
        tx <- pointDistance(mycont2x.s, mycont2y.s, lonlat=T, allpairs=T)
        my.min <- which(tx==min(tx), arr.ind=TRUE)
        this.lon <- mean(c(mycont2x.s[my.min[1,][1]]$V1, mycont2y.s[my.min[1,][2]]$V1))
        this.lat <- mean(c(mycont2x.s[my.min[1,][1]]$V2, mycont2y.s[my.min[1,][2]]$V2))
      }
      
      if(i==1){
        # build output df
        outpoint <- as.data.table(cbind(this.lon, this.lat)) ; colnames(outpoint) <- c('x', 'y')
        outpoint$model <- theseModels[i]
      } else{
        temppoint <- as.data.table(cbind(this.lon, this.lat)) ; colnames(temppoint) <- c('x', 'y')
        temppoint$model <- theseModels[i]
        outpoint <- as.data.table(rbind(outpoint, temppoint))
      }
      
    }
    
    return(outpoint)
  }
  
}

# -- TEST -- # 
# efe_efpm.fun(df=outlist$MH)
