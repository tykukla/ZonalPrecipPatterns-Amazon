# --------------------------------- #
# Pull out the EFE / EFPM analysis  #
# from the Boos & Korty model       #
# ---                               #
# T Kukla (Stanford Univ. 2021)     #
#                                   #
# Date created: Sept 1, 2021        #
# Last modified: Oct 28, 2021       #
# --------------------------------- #
library(ncdf4)
library(ncdf4.helpers)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(data.table)
library(metR)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(scico)
library(raster)
library(cmocean)
library(oce)
library(paletteer)
library(readr)

rm(list=ls())
# where DISTRIBUTE folder is located
setwd('C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data')

# set directory to save everything 
saveDir <- 'DISTRIBUTE/ENERGY_BALANCE/MID-HOLOCENE'

# cave sites for longitudinal reference
df <- as_tibble(read.csv('DISTRIBUTE/PROXY_COMP/Amz_proxyComp_Rversion.csv'))
mypts <- df[which(df$CaveSite=='Y'), ]
forcingline <- seq(-5,100, length=500)
cave.pts <- as.data.table(cbind(rep(mypts$Lat, each=length(forcingline)),
                                rep(mypts$Lon, each=length(forcingline))))
colnames(cave.pts) <- c("Lat", "Lon")
cave.pts$Site <- rep(mypts$Site_name, each=length(forcingline))
cave.pts$forceline <- rep(forcingline, nrow(mypts))



# ----------------------------------------------------------------- #
# READ IN THE OUTPUT (from ExtractPlot... scripts)
# [MH GREEN SAHARA]
dfmh <- readRDS(paste(saveDir, 'MH-GrSaharaSmallBand-contours-pts.RDS', sep='/'))
dfcontour <- dfmh$contours
dfx.p <- dfmh$intersects

## PLOT POINTS ON MAP
pal <- c('algae', 'amp', 'balance', 'diff',   
         'gray', 'curl', 'deep', 'delta',
         'dense', 'haline', 'ice', 'matter',
         'oxy', 'phase', 'rain', 'solar',
         'speed', 'tarn', 'tempo', 'thermal',
         'topo', 'turbid')
pal.num <- 17
#... bring in map for plotting
globePoly <- map_data('world')
landColor <- '#F2F2F2'
myLims <- c('min.x'=-85, 'max.x'=-20, 'min.y'=-20, 'max.y'=20)
ptColor <- '#ca0020'  # cave site colors

# ... pick contour lines to add
myForce <- unique(dfcontour$force_Wm2)
myLines <- c(myForce[1], myForce[5], myForce[10], myForce[15], myForce[20])

# PLOT map
pmap <- ggplot() +
  geom_polygon(data=globePoly, aes(x=long, y=lat, group=group), fill=landColor, color='#404040') + 
  coord_fixed(1.3, xlim=c(myLims['min.x'], myLims['max.x']), 
              ylim=c(myLims['min.y'], myLims['max.y'])) +
  # speleothem sites
  geom_text(data=mypts, aes(x=Lon, y=Lat, label="★"), size=10, family = "HiraKakuPro-W3",
            color=ptColor) +
  # add contour lines
  # -- ALL
  # geom_path(data=dfcontour[which(dfcontour$line=="EFPM")], 
  #           aes(x=lon, y=lat, group=force_Wm2, color=force_Wm2), size=1) +
  # -- SELECTED
  geom_path(data=dfcontour[dfcontour$force_Wm2 %in% myLines & dfcontour$line=="EFPM"],
            aes(x=lon, y=lat, group=force_Wm2), color='#594839', size=0.7) +
  # add data points
  # geom_label_repel(data=df[which(df$Depth_Profile=='Y'), ], aes(x=Lon, y=Lat, label=Sample_ID)) +
  # geom_point(data=dfx.p, aes(x=x, y=y), fill=NA, 
  #            shape=21, color='black', size=6, stroke = 1) +
  geom_point(data=dfx.p[force_Wm2 %in% myLines], aes(x=x, y=y, size=force_Wm2), 
             shape=21, color='#594839', fill='#F2EFBD', stroke = 1.5) +
  # scale_fill_cmocean(name=pal[pal.num], 'MSE Forcing\n(Wm2)', limits=c(0,70), start = 0.1) +
  # scale_color_cmocean(name=pal[pal.num], 'MSE Forcing\n(Wm2)', limits=c(0,70), start=0.1) +
  scale_size_continuous("MSE Forcing\n(Wm2)", breaks=c(0, 15, 35, 50, 70), range=c(1.5, 9)) +
  scale_x_continuous(expand=c(0,0), name="longitude") +
  scale_y_continuous(expand=c(0,0), name="latitude") +
  theme_few() +
  theme(panel.background = element_rect(fill='lightblue'))

pmap

saveFile <- paste(saveDir, 'MH_GrSah_EBM_Map.png', sep='/')
# ggsave(saveFile, plot=pmap, width=14, height=13, units=c('cm'))



# plotted vs long
plong <- ggplot(dfx.p) + 
  coord_cartesian(xlim=c(myLims['min.x'], myLims['max.x'])) + 
  # add cave longitudes 
  # geom_line(data=cave.pts, aes(x=forceline, y=Lon, group=Site), linetype='dotdash') +
  # EF Intersection long
  geom_point(aes(y=force_Wm2, x=x), size=4, shape=21, stroke =1) + 
  geom_line(aes(y=force_Wm2, x=x)) +
  labs(y=expression("Forcing (W/m"^"2"*")"), x="EF intersection longitude") +
  scale_x_continuous(expand=c(0,0)) +
  theme_few() +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=17),
        plot.margin = margin(1, 1, 1, 1, "cm"))

plong 

saveFile <- paste(saveDir, 'MH_GrSah_EBM_EFint.png', sep='/')
# ggsave(saveFile, plot=plong, width=20, height=8, units=c('cm'))



## ----- MAKE A FORCING MAP ---- ## 
df.all <- readRDS(paste(saveDir, 'MH-GrSahara_ALL.RDS', sep='/'))
# pick any two to take anomaly
df.x1 <- df.all$`Force-Wm2_0` 
df.x2 <- df.all$`Force-Wm2_37`
this.anom <- df.x1$MSE_source - df.x2$MSE_source
df.anom <- as.data.table(cbind(df.x1$lat, df.x1$lon, this.anom))
colnames(df.anom) <- c('lat', 'lon', 'MSE_source_anom')

# map lims
# myLims.anom <- c('min.x'=-130, 'max.x'=56, 'min.y'=-15, 'max.y'=40)
myLims.anom <- c('min.x'=-150, 'max.x'=170, 'min.y'=-15, 'max.y'=60)

p.force <- ggplot() +
  coord_fixed(1.3, xlim=c(myLims.anom['min.x'], myLims.anom['max.x']), 
              ylim=c(myLims.anom['min.y'], myLims.anom['max.y'])) +
  geom_polygon(data=globePoly, aes(x=long, y=lat, group=group), fill=NA, color='#404040') + 
  # raster dat
  geom_tile(data=df.anom, aes(x=lon, y=lat, fill=MSE_source_anom, alpha=-MSE_source_anom)) +
  # speleothem sites
  geom_text(data=mypts, aes(x=Lon, y=Lat, label="★"), size=5, family = "HiraKakuPro-W3",
            color=ptColor)  +
  scale_fill_cmocean(name="gray", guide='none') +
  scale_alpha_binned(guide='none') +
  theme_few() +
  theme(axis.title = element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) 

p.force 

saveFile <- paste(saveDir, 'MH_GrSah_ForceMap.png', sep='/')
# ggsave(saveFile, plot=p.force, width=20, height=8, units=c('cm'))



