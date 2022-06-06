# ----------------------------------------- # 
# Script to plot SISAL range data and map   #
# For amazon paper                          #
# ---                                       # 
# T Kukla (Colostate Univ. 2022)            #
#                                           #
# Date created: May 2, 2022                 #
# Last modified: May 11, 2022               #
# ----------------------------------------- #
library(ggthemes)
library(data.table)
library(ggpubr)
library(ggplot2)
library(ggtext)

rm(list=ls())

# directory where DISTRIBUTE is located
user.dir <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data'
# shouldn't need to change
mainPath <- 'DISTRIBUTE/ISOTOPE_DATA/SISAL'
setwd(paste(user.dir, mainPath, sep='/'))

## READ IN DATA
df <- readRDS('SISAL_Smooth_10kyrMin-1kyrAvg.RDS')
df.sum.amz <- df$summary_all
df.sum <- df$`summary_no-amazon`
df.sum.noDTP <- df$`summary_no-added-DTP`
# add quantile data
df.sum$qtile <- ecdf(df.sum$d18_range)(df.sum$d18_range)
df.sum.amz$qtile <- ecdf(df.sum.amz$d18_range)(df.sum.amz$d18_range)


## WHERE TO SAVE
save.here <- mainPath

## GET STATS for amazon sites
# ----------------------------------------------------------------------------- 
# -----------------------------------------------------------------------------
# set limits
lat.limit <- 40 ; duration.limit <- 1e5
# add quantile data
df.sum$qtile <- ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18_range)(df.sum$d18_range)
df.sum.amz$qtile <- ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18_range)(df.sum.amz$d18_range)

# get range for amz sites
PAR.range <- df.sum.amz[site_id == 3]$d18_range
RGN.range <- df.sum.amz[site_id==111]$d18_range
DTP.range <- df.sum.amz[site_id == -9999]$d18_range
PAR.sd <- df.sum.amz[site_id == 3]$d18.sd
RGN.sd <- df.sum.amz[site_id==111]$d18.sd
DTP.sd <- df.sum.amz[site_id == -9999]$d18.sd
# quick summary statistics
ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18_range)(RGN.range)
ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18_range)(PAR.range)
ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18_range)(DTP.range)
ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18.sd)(RGN.sd)
ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18.sd)(PAR.sd)
ecdf(df.sum[abs(lat) < lat.limit & age_range < duration.limit]$d18.sd)(DTP.sd)
# age range restriction
ecdf(df.sum[age_range < 1e5]$d18_range)(RGN.range)
ecdf(df.sum[age_range < 1e5]$d18_range)(PAR.range)
ecdf(df.sum[age_range < 1e5]$d18_range)(DTP.range)

quantile(na.omit(df.sum$d18_range), c(0.1, 0.5, 0.9))
quantile(na.omit(df.sum$d18.sd), c(0.1, 0.5, 0.9))


# ----------------------------------------------------------------------------- 
## HISTOGRAM PLOT
p.hist <- ggplot(df.sum[amz.site==F & abs(lat) < lat.limit & age_range < duration.limit]) +
  coord_cartesian(xlim=c(0, 8), ylim=c(0,11)) +
  annotate('segment', x=RGN.range, xend=RGN.range, y=0, yend=15, color='darkred', size=1.5) +
  annotate('segment', x=PAR.range, xend=PAR.range, y=0, yend=15, color='darkred', size=1.5) +
  annotate('segment', x=DTP.range, xend=DTP.range, y=0, yend=15, color='darkred', size=1.5) +
  geom_histogram(aes(x=d18_range, fill=lat), binwidth = 0.5, color='white') + 
  scale_y_continuous(expand=c(0,0), name='count', breaks=c(0,2,4,6,8,10,12)) + 
  scale_x_continuous(expand=c(0,0), name="*&delta;*<sup>18</sup>O range (&permil;)") +
  theme_few() +
  theme(axis.title.x = element_markdown(size=20), axis.text = element_text(size=15),
        axis.title.y = element_markdown(size=20))
  
p.hist

save.name <- 'range_hist-1kyrAvg-1e5yrDurLim-sub40Lat.png'
# ggsave(paste(save.here, save.name, sep='/'), p.hist, width=15, height=10, units='cm')


# ----------------------------------------------------------------------------- 
## MAP PLOT
globePoly <- map_data('world')
landColor <- '#F2F2F2'
ptColor <- '#ca0020'

myLims <- c('min.x'=-170, 'max.x'=170, 'min.y'=-60, 'max.y'=70)

p.map <- ggplot() + 
  # the map and axes
  geom_polygon(data=globePoly, aes(x=long, y=lat, group=group), fill=landColor, color='#404040') + 
  coord_fixed(1.3, xlim=c(myLims['min.x'], myLims['max.x']), 
              ylim=c(myLims['min.y'], myLims['max.y'])) +
  geom_point(data=df.sum.amz[abs(lat) < lat.limit & age_range < duration.limit & qtile < 0.95], aes(x=lon, y=lat),
             fill='black', color='darkgray', shape=21, size=4, stroke = 0.7) +
  geom_point(data=df.sum.amz[abs(lat) < lat.limit & age_range < duration.limit & qtile >= 0.95], aes(x=lon, y=lat),
             fill=ptColor, color='black', shape=21, size=7, stroke = 0.9) +
  # geom_point(data=df.sum.amz[amz.site==T], aes(x=lon, y=lat, color=d18_range, size=d18_range)) +
  labs(x=NULL, y=NULL) +
  theme_few() + 
  theme(panel.background = element_rect(fill='lightblue'), axis.text = element_blank())

p.map 

save.name <- 'range_Map-95qtile-1kyrAvg-1e5yrDurLim-sub40Lat.png'
# ggsave(paste(save.here, save.name, sep='/'), p.map, width=20, height=14, units='cm')


