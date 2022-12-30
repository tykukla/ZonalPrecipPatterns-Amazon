# ------------------------------------- # 
# PRECIPITATION AND EVAPORATION SHED    #
# PLOTS -- using Xia output from        # 
# WAM-2layers                           #
# ---                                   #
# T Kukla (UW; 2022)                    #
# ------------------------------------- # 
library(ggplot2)
library(ggthemes)
library(reshape2)
library(cmocean)
library(data.table)
library(mapdata)
rm(list=ls())

# working dir is file dir
my.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my.dir)   # location of WAM-2layers output

## FIRST CALCULATE RESULTS
# calculate Psheds and Esheds for use
source("each shed.R") # this files is for calculation

## NOW ORGANIZE IN A LIST
etracks <- list('dtp' = etrack_DTP, 'par' = etrack_PAR, 'rgn' = etrack_RGN) 
esheds <- list('dtp' = eshed_DTP, 'par' = eshed_PAR, 'rgn' = eshed_RGN)
ptracks <- list('dtp' = ptrack_DTP, 'par' = ptrack_PAR, 'rgn' = ptrack_RGN)
psheds <- list('dtp' = pshed_DTP, 'par' = pshed_PAR, 'rgn' = pshed_RGN)

## Where to save results
save.here <- paste(my.dir, 'panels', sep='/')  # PATH FOR STORING FIGURES 
## BY SITE
SITE.NAMES <- c('dtp', 'par', 'rgn')
TIME.NAMES <- list('all'=13, 'djf' = c(12,1,2), 'mam' = c(3,4,5), 
                   'jja' = c(6,7,8), 'son' = c(9,10,11), 'ndjfm' = c(11,12,1,2,3),
                   'ndjfmam' = c(11,12,1,2,3,4,5))

## COORDS for proxy sites
lat.prx <- c(-5.73, -4.06, -5.6) ; lon.prx <- c(-77.49, -55.45, -37.73); 
names.prx <- c('dtp', 'par', 'rgn')
prx.coords <- as.data.table(cbind(lat.prx, lon.prx)) ; colnames(prx.coords) <- c('lat', 'lon')
prx.coords$site <- names.prx


for(times in 1:length(TIME.NAMES)){
  THIS.TIME <- idx <- TIME.NAMES[[times]]
  THIS.TIME.NAME <- names(TIME.NAMES)[times]
  for(sites in 1:length(SITE.NAMES)){
    THIS.SITE <- SITE.NAMES[sites]
    
    # set up map
    worldmap <- map_data("world", regions=".", interior=F) #  map('world',fill=F,interior=F,plot=F)  # map_data("world") 
    lims.lat <- c(-40, 40) ; lims.lon <- c(-130, 110)
    # worldmap <- map2SpatialLines(worldmap,proj4string=CRS("+proj=longlat"))
    lakemap <- map("lakes",fill=F,interior=F,plot=F)
    # lakemap <- map2SpatialLines(lakemap,proj4string=CRS("+proj=longlat"))
    # define grid
    lats <- seq(79.5,-79.5,-1.5)
    lons <- seq(-180,180,1.5)
    
    # ... CONTOUR
    my.etrack <- etracks[[THIS.SITE]] ; my.eshed <- esheds[[THIS.SITE]]
    my.ptrack <- ptracks[[THIS.SITE]] ; my.pshed <- psheds[[THIS.SITE]]
    ## [PSHEDS]
    # colmeans takes average of idx (breaks if length(idx)==1, so we remove in this case)
    if(length(idx) > 1){
      e.1 <- cbind(colMeans(my.etrack[idx,,]), colMeans(my.etrack[idx,,1]))/sum(colMeans(my.etrack[idx,,]))*100 # plot raster grids require adding one more longitude
    } else{
      e.1 <- cbind(my.etrack[idx,,], my.etrack[idx,,1])/sum(my.etrack[idx,,])*100 # plot raster grids require adding one more longitude
    }
    e <- structure(as.vector(e.1), .Dim=c(107L, 241L), .Dimnames = list(NULL, NULL))
    # add lat / lon
    dimnames(e) <- list(x=lats, y=lons)
    e.df <- reshape2::melt(e) ; colnames(e.df) <- c('lat', 'lon', 'pshed') ; e.df <- as.data.table(e.df)
    e.df$pshed <- ifelse(e.df$pshed > 1, 1, e.df$pshed)  # max 1 for color purposes
    # ... THRESHOLD
    if(length(idx) > 1){
      p.thresh <- t(cbind(colMeans(my.pshed[idx,,]),colMeans(my.pshed[idx,,1])))[,107:1]
    } else{
      p.thresh <- t(cbind(my.pshed[idx,,],my.pshed[idx,,1]))[,107:1]
    }
    pt.mat <- structure(as.vector(p.thresh), .Dim=c(241L, 107L), .Dimnames = list(NULL, NULL))
    dimnames(pt.mat) <- list(x=lons, y=rev(lats))
    pt.df <- reshape2::melt(pt.mat) ; colnames(pt.df) <- c('lon', 'lat', 'pthreshold') ; pt.df <- as.data.table(pt.df)
    
    ## [ESHEDS]
    # colmeans takes average of idx (breaks if length(idx)==1, so we remove in this case)
    if(length(idx) > 1){
      p.1 <- cbind(colMeans(my.ptrack[idx,,]), colMeans(my.ptrack[idx,,1]))/sum(colMeans(my.ptrack[idx,,]))*100 # plot raster grids require adding one more longitude
    } else{
      p.1 <- cbind(my.ptrack[idx,,], my.ptrack[idx,,1])/sum(my.ptrack[idx,,])*100 # plot raster grids require adding one more longitude
    }
    p <- structure(as.vector(p.1), .Dim=c(107L, 241L), .Dimnames = list(NULL, NULL))
    # add lat / lon
    dimnames(p) <- list(x=lats, y=lons)
    p.df <- reshape2::melt(p) ; colnames(p.df) <- c('lat', 'lon', 'eshed') ; p.df <- as.data.table(p.df)
    p.df$eshed <- ifelse(p.df$eshed > 1, 1, p.df$eshed)  # max 1 for color purposes
    # ... THRESHOLD
    if(length(idx) > 1){
      e.thresh <- t(cbind(colMeans(my.eshed[idx,,]),colMeans(my.eshed[idx,,1])))[,107:1]
    } else{
      e.thresh <- t(cbind(my.eshed[idx,,],my.eshed[idx,,1]))[,107:1]
    }
    e.mat <- structure(as.vector(e.thresh), .Dim=c(241L, 107L), .Dimnames = list(NULL, NULL))
    dimnames(e.mat) <- list(x=lons, y=rev(lats))
    et.df <- reshape2::melt(e.mat) ; colnames(et.df) <- c('lon', 'lat', 'ethreshold') ; et.df <- as.data.table(et.df)
    
    
    
    ## [PLOT]
    ## [0] -- LEGEND IF NEEDED
    if(sites == 1 & times == 1){
      ## COLORBAR ON BOTTOM: https://stackoverflow.com/questions/54029017/ggplot-extend-legend-colorbar 
      p.legend <- ggplot(e.df[pshed > 0.02]) +
        coord_cartesian(xlim=lims.lon, ylim=lims.lat) + 
        geom_raster(aes(x=lon, y=lat, fill=pshed)) + 
        scale_fill_cmocean(name='ice', direction=-1, limits=c(0,1)) + 
        geom_polygon(data=worldmap, aes(long, lat, group = group), fill=NA, color='black') +
        geom_polygon(data=lakemap, aes(long, lat, group = group), fill=NA, color='black') +
        geom_point(data=prx.coords, aes(x=lon, y=lat), shape=24, fill='#BD2A2E', color='black', size=30, stroke=1) + 
        theme_few() +
        theme(axis.text = element_blank(), axis.title=element_blank(), axis.ticks = element_blank(),
              legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
      
      # p.legend
      save.name <- 'LEGEND.png'
      save.path <- paste(save.here, save.name, sep='/')
      ggsave(save.path, plot=p.legend, width=20, height=8, units='cm')
    }
    ## [1] -- P SHEDS
    p.1 <- ggplot(e.df[pshed > 0.02]) +
      coord_cartesian(xlim=lims.lon, ylim=lims.lat) + 
      geom_raster(aes(x=lon, y=lat, fill=pshed)) + 
      scale_fill_cmocean(name='ice', direction=-1, limits=c(0,1), guide='none') + 
      # geom_contour_filled(aes(x=lon, y=lat, z=dtp_pshed)) +
      geom_polygon(data=worldmap, aes(long, lat, group = group), fill=NA, color='black') +
      geom_polygon(data=lakemap, aes(long, lat, group = group), fill=NA, color='black') +
      geom_contour(data=pt.df, aes(x=lon, y=lat, z=pthreshold), breaks=1, color='magenta', size=1) +
      geom_point(data=prx.coords, aes(x=lon, y=lat), shape=24, fill='#BD2A2E', color='black', size=4, stroke=1) + 
      theme_few() +
      theme(axis.text = element_blank(), axis.title=element_blank(), axis.ticks = element_blank())
    
    # p.1
    save.name <- paste(paste('PSHED', THIS.SITE, THIS.TIME.NAME, sep='_'), '.png', sep='')
    save.path <- paste(save.here, save.name, sep='/')
    ggsave(save.path, plot=p.1, width=20, height=8, units='cm')
    
    ## [2] -- E SHEDS
    p.2 <- ggplot(p.df[eshed > 0.02]) +
      coord_cartesian(xlim=lims.lon, ylim=lims.lat) + 
      geom_tile(aes(x=lon, y=lat, fill=eshed)) + 
      scale_fill_cmocean(name='ice', direction=-1, limits=c(0,1), guide='none') + 
      # geom_contour_filled(aes(x=lon, y=lat, z=dtp_pshed)) +
      geom_polygon(data=worldmap, aes(long, lat, group = group), fill=NA, color='black') +
      geom_polygon(data=lakemap, aes(long, lat, group = group), fill=NA, color='black') +
      geom_contour(data=et.df, aes(x=lon, y=lat, z=ethreshold), breaks=1, color='magenta', size=1) + 
      geom_point(data=prx.coords, aes(x=lon, y=lat), shape=24, fill='#BD2A2E', color='black', size=4, stroke=1) +
      theme_few() +
      theme(axis.text = element_blank(), axis.title=element_blank(), axis.ticks = element_blank())
    
    p.2
    save.name <- paste(paste('ESHED', THIS.SITE, THIS.TIME.NAME, sep='_'), '.png', sep='')
    save.path <- paste(save.here, save.name, sep='/')
    ggsave(save.path, plot=p.2, width=20, height=8, units='cm')
    
    
  } # end site for loop
} # end times loop

