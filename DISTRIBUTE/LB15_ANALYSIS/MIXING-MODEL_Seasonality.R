# ------------------------------------- # 
# COMPUTE ANOMALY FROM PRECIP-WTD MEANS #
# OF LIU AND BATTISTI 2015 RESULTS      #
# ---                                   #
# T Kukla (UW; 2022)                    #
# ------------------------------------- # 
library(ggplot2)
library(ggthemes)
library(reshape2)
library(cmocean)
library(data.table)

rm(list=ls())

# working dir is file dir
my.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(my.dir)   # location of LB15 files

# read in graph-clicked data 
fn <- 'LB15_final.csv'
df <- as.data.table(read.csv(fn))
colnames(df) <- c('month', 'pr_mmday', 'd18O', 'case')

# calculate weighted means
cases <- c('low_insol', 'high_insol'); d18.wt <- vector()
for(i in 1:length(cases)){
  this.case <- cases[i]
  d18.wt[i] <- weighted.mean(df[case==this.case]$d18O, df[case==this.case]$pr_mmday)
}
# bring together
df.wt <- as.data.table(cbind(d18.wt))
df.wt$case <- cases
LB15.anom <- df.wt[case == 'low_insol']$d18.wt - df.wt[case == 'high_insol']$d18.wt

# [QUICK PLOT]
# where to save
save.here <- paste(my.dir, 'panels', sep='/')
# ... precip
p1 <- ggplot(df) + 
  geom_line(aes(x=month, y=pr_mmday, color=case), size=2) + 
  geom_point(aes(x=month, y=pr_mmday, fill=case), size=6, stroke=1, shape=21, color='black') + 
  scale_color_cmocean(discrete=T, name='tempo', start=0.2, end = 0.9) + 
  scale_fill_cmocean(discrete=T, name='tempo', start=0.2, end = 0.9) + 
  scale_x_continuous(name='month', breaks=c(1:12)) + 
  scale_y_continuous(name = "Precipitation (mm/month)") +
  theme_few() +
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18))
p1
save.name <- paste(save.here, 'pr_monthly_LB15.png', sep='/')
# ggsave(save.name, p1, width=20, height=12, units='cm')

# ... d18O
p2 <- ggplot(df) + 
  geom_line(aes(x=month, y=d18O, color=case), size=2) + 
  geom_point(aes(x=month, y=d18O, fill=case), size=6, stroke=1, shape=21, color='black') + 
  scale_color_cmocean(discrete=T, name='tempo', start=0.2, end = 0.9) + 
  scale_fill_cmocean(discrete=T, name='tempo', start=0.2, end = 0.9) + 
  scale_x_continuous(name='month', breaks=c(1:12)) + 
  scale_y_continuous(name=expression(delta^"18"*'O'['precipitation'])) +
  theme_few() +
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18))
p2 
save.name <- paste(save.here, 'd18_monthly_LB15.png', sep='/')
# ggsave(save.name, p2, width=20, height=12, units='cm')




## SIGNAL LIMITS
max.sig <- -8  # add this value to make sure we get the right color bar

# -------------------------------------------------------------------------------------------------
## [HOLD d18O CONSTANT, CHANGE PRECIP SEASONALITY] ---------------------------------------------
## [SET MODERN CONTROL TO OBS]
df.x <- df
df.x[case == 'high_insol']$pr_mmday <- df.x[case == 'obs']$pr_mmday
iters <- 1e3 # not used
# these are "additional anoms" -- the additional change relative to control
winter.P.anoms <- seq(-2, 2, length=50) ; summer.P.anoms <- seq(-2, 15, length = 50)
apply.anom.to <- 'low_insol' ; ctrl.case <- 'high_insol'
summer <- c(12,1,2,3,4) ; winter <- c(6,7,8)  # months 

this.wP.anom <- this.sP.anom <- this.d18.signal <- vector()
idx <- 1
## this may take a minute or so... apologies for lazy coding
for(j in 1:length(winter.P.anoms)){
  for(i in 1:length(summer.P.anoms)){
    # create temporary df
    tdf <- df.x
    # randomly sample anoms
    this.wP.anom[idx] <- winter.P.anoms[j]
    this.sP.anom[idx] <- summer.P.anoms[i]
    
    # modify the low insol case
    tdf[case==apply.anom.to & month %in% winter]$pr_mmday <- tdf[case==apply.anom.to & month %in% winter]$pr_mmday + this.wP.anom[idx]
    tdf[case==apply.anom.to & month %in% summer]$pr_mmday <- tdf[case==apply.anom.to & month %in% summer]$pr_mmday + this.sP.anom[idx]
    # remove negative rainfall
    tdf[case==apply.anom.to & month %in% winter & pr_mmday < 0]$pr_mmday <- 0
    tdf[case==apply.anom.to & month %in% summer & pr_mmday < 0]$pr_mmday <- 0
    
    # calculate weighted mean d18O 
    d18.ctrl <- weighted.mean(tdf[case == ctrl.case]$d18O, tdf[case == ctrl.case]$pr_mmday)
    d18.case <- weighted.mean(tdf[case == apply.anom.to]$d18O, tdf[case == apply.anom.to]$pr_mmday)
    
    # get signal
    this.d18.signal[idx] <- d18.case - d18.ctrl
    
    # set idx
    idx <- idx + 1 
  }
}

# bring together
df.anom <- as.data.table(cbind(this.wP.anom, this.sP.anom, this.d18.signal))
colnames(df.anom) <- c('winterP.add_anom', 'summerP.add_anom', 'd18.signal')
# calculate the unaltered anom
orig.sP.anom <- mean(df.x[case == apply.anom.to & month %in% summer]$pr_mmday - df.x[case == ctrl.case & month %in% summer]$pr_mmday)
orig.wP.anom <- mean(df.x[case == apply.anom.to & month %in% winter]$pr_mmday - df.x[case == ctrl.case & month %in% winter]$pr_mmday)
# get total anoms
df.anom$summerP.tot_anom <- df.anom$summerP.add_anom + orig.sP.anom
df.anom$winterP.tot_anom <- df.anom$winterP.add_anom + orig.wP.anom

# add value to get color bar working 
tdf.anom <- as.data.table(cbind(-1e3, -1e3, max.sig, -1e3, -1e3)) ; colnames(tdf.anom) <- colnames(df.anom)
df.anom <- as.data.table(rbind(df.anom, tdf.anom))
for.lims <- df.anom[winterP.add_anom > -100]

## [QUICK PLOT]
p4 <- ggplot(df.anom) + 
  coord_cartesian(xlim=c(min(for.lims$summerP.tot_anom), max(for.lims$summerP.tot_anom)),
                  ylim=c(min(for.lims$winterP.tot_anom), max(for.lims$winterP.tot_anom))) +
  geom_tile(aes(x=summerP.tot_anom, y=winterP.tot_anom, fill=d18.signal), alpha=0.7) +
  scale_fill_cmocean(name='thermal') +
  scale_x_continuous(expand=c(0,0), name='DJFMA precip. anomaly (mm/month)') + 
  scale_y_continuous(expand=c(0,0), name='JJA precip. anomaly (mm/month)') + 
  theme_few()  +
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18))
p4
save.name <- paste(save.here, 'LB15_Const-d18O.png', sep='/')
# ggsave(save.name, p4, width=20, height=12, units='cm')




# -------------------------------------------------------------------------------------------------
## [HOLD PRECIP CONSTANT, CHANGE d18O SEASONALITY] ---------------------------------------------
df.x <- df
df.x[case == 'high_insol']$pr_mmday <- df.x[case == 'obs']$pr_mmday
iters <- 1e3 # not used
# these are "additional anoms" -- the additional change relative to control
winter.d18.anoms <- seq(-8, 3, length=50) ; summer.d18.anoms <- seq(-8, 3, length = 50)
apply.anom.to <- 'low_insol' ; ctrl.case <- 'high_insol'
summer <- c(12,1,2,3,4) ; winter <- c(6,7,8)  # months 

this.wd18.anom <- this.sd18.anom <- this.d18.signal <- vector()
idx <- 1
## this may take a minute or so... apologies for lazy coding
for(j in 1:length(winter.d18.anoms)){
  for(i in 1:length(summer.d18.anoms)){
    # create temporary df
    tdf <- df.x
    # randomly sample anoms
    this.wd18.anom[idx] <- winter.d18.anoms[j]
    this.sd18.anom[idx] <- summer.d18.anoms[i]
    
    # modify the low insol case
    tdf[case==apply.anom.to & month %in% winter]$d18O <- tdf[case==apply.anom.to & month %in% winter]$d18O + this.wd18.anom[idx]
    tdf[case==apply.anom.to & month %in% summer]$d18O <- tdf[case==apply.anom.to & month %in% summer]$d18O + this.sd18.anom[idx]
    
    # calculate weighted mean d18O 
    d18.ctrl <- weighted.mean(tdf[case == ctrl.case]$d18O, tdf[case == ctrl.case]$pr_mmday)
    d18.case <- weighted.mean(tdf[case == apply.anom.to]$d18O, tdf[case == apply.anom.to]$pr_mmday)
    
    # get signal
    this.d18.signal[idx] <- d18.case - d18.ctrl
    
    # set idx
    idx <- idx + 1 
  }
}

# bring together
df.anom <- as.data.table(cbind(this.wd18.anom, this.sd18.anom, this.d18.signal))
colnames(df.anom) <- c('winter.d18.add_anom', 'summer.d18.add_anom', 'd18.signal')
# calculate the unaltered anom
orig.sd18.anom <- mean(df.x[case == apply.anom.to & month %in% summer]$d18O - df.x[case == ctrl.case & month %in% summer]$d18O)
orig.wd18.anom <- mean(df.x[case == apply.anom.to & month %in% winter]$d18O - df.x[case == ctrl.case & month %in% winter]$d18O)
# get total anoms
df.anom$summer.d18.tot_anom <- df.anom$summer.d18.add_anom + orig.sd18.anom
df.anom$winter.d18.tot_anom <- df.anom$winter.d18.add_anom + orig.wd18.anom

# add value to get color bar working 
tdf.anom <- as.data.table(cbind(-1e3, -1e3, max.sig, -1e3, -1e3)) ; colnames(tdf.anom) <- colnames(df.anom)
df.anom <- as.data.table(rbind(df.anom, tdf.anom))
for.lims <- df.anom[winter.d18.add_anom > -100]
# overwrite max
df.anom[d18.signal < max.sig]$d18.signal <- max.sig

## [QUICK PLOT]
p5 <- ggplot(df.anom) + 
  coord_cartesian(xlim=c(min(for.lims$summer.d18.tot_anom), max(for.lims$summer.d18.tot_anom)),
                  ylim=c(min(for.lims$winter.d18.tot_anom), max(for.lims$winter.d18.tot_anom))) +
  geom_tile(aes(x=summer.d18.tot_anom, y=winter.d18.tot_anom, fill=d18.signal), alpha=0.6) +
  geom_tile(data=df.anom[d18.signal > -7 & d18.signal < -5],
            aes(x=summer.d18.tot_anom, y=winter.d18.tot_anom, fill=d18.signal)) +
  geom_contour(aes(x=summer.d18.tot_anom, y=winter.d18.tot_anom, z=d18.signal), 
               breaks=c(-5, -7), color='white', size=3) +
  scale_fill_cmocean(name='thermal') +
  scale_x_continuous(expand=c(0,0), name=expression("DJFMA"~delta^"18"*"O anomaly (???)")) +
  scale_y_continuous(expand=c(0,0), name=expression("JJA"~delta^"18"*"O anomaly (???)")) +
  theme_few() +
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18))
p5
save.name <- paste(save.here, 'LB15_Const-Pr.png', sep='/')
ggsave(save.name, p5, width=20, height=12, units='cm')




# -------------------------------------------------------------------------------------------------
## [CHANGE DJFMA SUMMER P AND SUMMER d18O] ---------------------------------------------
## [SET MODERN CONTROL TO OBS]
df.x <- df
df.x[case == 'high_insol']$pr_mmday <- df.x[case == 'obs']$pr_mmday
iters <- 1e3
summer.d18.anoms <- seq(-8, 2, length=50) ; summer.p.anoms <- seq(-2, 10, length = 50)
apply.anom.to <- 'low_insol' ; ctrl.case <- 'high_insol'
summer <- c(12,1,2,3,4) ; winter <- c(6,7,8)

this.d18.anom <- this.p.anom <- this.d18.signal <- vector()
idx <- 1
## this may take a minute or so... apologies for lazy coding
for(j in 1:length(summer.d18.anoms)){
  for(i in 1:length(summer.p.anoms)){
    # create temporary df
    tdf <- df.x
    # randomly sample anoms
    this.d18.anom[idx] <- summer.d18.anoms[j]
    this.p.anom[idx] <- summer.p.anoms[i]
    
    # modify the low insol case
    tdf[case==apply.anom.to & month %in% summer]$d18O <- tdf[case==apply.anom.to & month %in% summer]$d18O + this.d18.anom[idx]
    tdf[case==apply.anom.to & month %in% summer]$pr_mmday <- tdf[case==apply.anom.to & month %in% summer]$pr_mmday + this.p.anom[idx]
    # remove negative rainfall
    tdf[case==apply.anom.to & month %in% summer & pr_mmday < 0]$pr_mmday <- 0
    
    # calculate weighted mean d18O 
    d18.ctrl <- weighted.mean(tdf[case == ctrl.case]$d18O, tdf[case == ctrl.case]$pr_mmday)
    d18.case <- weighted.mean(tdf[case == apply.anom.to]$d18O, tdf[case == apply.anom.to]$pr_mmday)
    # get signal
    this.d18.signal[idx] <- d18.case - d18.ctrl
    
    # set idx
    idx <- idx + 1 
  }
}

# bring together
df.anom <- as.data.table(cbind(this.d18.anom, this.p.anom, this.d18.signal))
colnames(df.anom) <- c('summer_d18.add_anom', 'summer_p.add_anom', 'd18.signal')
# calculate the unaltered anom
orig.d18.anom <- mean(df.x[case == apply.anom.to & month %in% summer]$d18O - df.x[case == ctrl.case & month %in% summer]$d18O)
orig.p.anom <- mean(df.x[case == apply.anom.to & month %in% summer]$pr_mmday - df.x[case == ctrl.case & month %in% summer]$pr_mmday)
# get total anoms
df.anom$summer_d18.tot_anom <- df.anom$summer_d18.add_anom + orig.d18.anom
df.anom$summer_p.tot_anom <- df.anom$summer_p.add_anom + orig.p.anom

# add value to get color bar working 
tdf.anom <- as.data.table(cbind(-1e3, -1e3, max.sig, -1e3, -1e3)) ; colnames(tdf.anom) <- colnames(df.anom)
df.anom <- as.data.table(rbind(df.anom, tdf.anom))
for.lims <- df.anom[summer_d18.add_anom > -100]
# overwrite max
df.anom[d18.signal < max.sig]$d18.signal <- max.sig

## [QUICK PLOT]
p5 <- ggplot(df.anom) + 
  coord_cartesian(xlim=c(min(for.lims$summer_d18.tot_anom), max(for.lims$summer_d18.tot_anom)),
                  ylim=c(min(for.lims$summer_p.tot_anom), max(for.lims$summer_p.tot_anom))) +
  geom_tile(aes(x=summer_d18.tot_anom, y=summer_p.tot_anom, fill=d18.signal), alpha=0.6) +
  geom_tile(data=df.anom[d18.signal > -7 & d18.signal < -5],
            aes(x=summer_d18.tot_anom, y=summer_p.tot_anom, fill=d18.signal)) +
  geom_contour(aes(x=summer_d18.tot_anom, y=summer_p.tot_anom, z=d18.signal), 
               breaks=c(-5, -7), color='white', size=3) +
  scale_fill_cmocean(name='thermal') +
  scale_x_continuous(expand=c(0,0), name=expression("DJFMA"~delta^"18"*"O anomaly (???)")) +
  scale_y_continuous(expand=c(0,0), name=expression("DJFMA precip. anomaly (mm/month)")) +
  theme_few()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18))
p5
save.name <- paste(save.here, 'LB15_Winter.png', sep='/')
ggsave(save.name, p5, width=20, height=12, units='cm')

## JUST FOR THE LEGEND
pleg <- ggplot(df.anom) + 
  coord_cartesian(xlim=c(min(for.lims$summer_d18.tot_anom), max(for.lims$summer_d18.tot_anom)),
                  ylim=c(min(for.lims$summer_p.tot_anom), max(for.lims$summer_p.tot_anom))) +
  geom_tile(aes(x=summer_d18.tot_anom, y=summer_p.tot_anom, fill=d18.signal), alpha=0.6) +
  geom_tile(data=df.anom[d18.signal > -7 & d18.signal < -5],
            aes(x=summer_d18.tot_anom, y=summer_p.tot_anom, fill=d18.signal)) +
  geom_contour(aes(x=summer_d18.tot_anom, y=summer_p.tot_anom, z=d18.signal), 
               breaks=c(-5, -7), color='white', size=3) +
  scale_fill_cmocean(name='thermal') +
  scale_x_continuous(expand=c(0,0), name=expression("DJFMA"~delta^"18"*"O anomaly (???)")) +
  scale_y_continuous(expand=c(0,0), name=expression("DJFMA precip. anomaly (mm/month)")) +
  theme_few()+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=18),
        legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
pleg
save.name <- paste(save.here, 'LB15_LEG.png', sep='/')
ggsave(save.name, pleg, width=20, height=12, units='cm')



