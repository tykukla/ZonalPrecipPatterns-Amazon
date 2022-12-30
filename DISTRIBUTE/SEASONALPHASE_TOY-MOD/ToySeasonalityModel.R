# ---------------------------------------------- # 
# Toy model for impact of spatially variable     #
# seasonality on isotope gradient across a storm # 
# track                                          #
# ---                                            #
# T Kukla (UW; Colostate Univ. 2022)             #
# ---------------------------------------------- #
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(cmocean)

rm(list=ls())

setwd('C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Fsupp_SeasonalityMod')

# Seasonal cycle parameters
tx <- seq(0,1,length=100)    # time axis (0-1 spans 1 year)
amp <- 1                     # d18O seasonal amplitude
phases <- seq(0, 2*pi, length=100)  # phase of downwind p (upwind P phase held constant)
d18.sensitivity <- seq((1/3), 3, length=10)  # differences in the sensitivity of d18 to P (basically changing relative magnitude of d18 change)

## [LOOP THROUGH PHASES]
idx <- 1 ; this.phase <- this.sensitivity <- Dd18.pwt <- Dd18.magnitude <- vector()
for(i in 1:length(phases)){  # loop through phases
  for(j in 1:length(d18.sensitivity)){   # loop through d18 sensitivity
    this.phase[idx] <- phases[i]
    this.sensitivity[idx] <- d18.sensitivity[j]
    
    ## [UPWIND SITE]
    # upwind P is a cosine fxn
    upwind.P <- amp * cos(2 * pi * tx) + amp   # add amp at end so P is never negative
    # upwind d18 anti-correlated w/ upwind P
    upwind.d18 <- 0 + ((upwind.P*-1)/this.sensitivity[idx])   # d18O anticorrelated with P and restricted to be negative (zero is for consistency with downwind eq'n)
    
    ## [DOWNWIND SITE]
    # downwind P is a cosine fxn offset from upwind by 'phases[i]'
    downwind.P <- amp * cos(2 * pi * tx + this.phase[idx]) + amp  # this is our 'rainout' across domain
    # downwind d18 depends on downwind P relative to the initial (upwind) signal
    downwind.d18 <- upwind.d18 + ((downwind.P*-1)/this.sensitivity[idx])
    
    ## [ISOTOPE GRADIENT]
    # seasonal iso grad
    Dd18 <- downwind.d18 - upwind.d18  # assumes sites separated by 1000 km (so units are consistent w/ text)
    # weighted mean iso grad
    Dd18.pwt[idx] <- weighted.mean(Dd18, downwind.P) # downwind.P, recall, is net rainout metric
    Dd18.magnitude[idx] <- max(Dd18) - min(Dd18)
    ## [OUTPUT DATA]
    if(i==1 & j==1){
      outdf <- as.data.table(cbind(tx, upwind.d18, upwind.P, downwind.d18, downwind.P, Dd18))
      outdf$phase <- this.phase[idx]
      outdf$d18.sensitivity <- 1/this.sensitivity[idx]
    } else{
      tdf <- as.data.table(cbind(tx, upwind.d18, upwind.P, downwind.d18, downwind.P, Dd18))
      tdf$phase <- this.phase[idx]
      tdf$d18.sensitivity <- 1/this.sensitivity[idx]
      outdf <- as.data.table(rbind(outdf, tdf))
    }
    
    # update idx
    idx <- idx + 1 
  }
  
}
# bring together summary Dd18 dat
df.Dd18 <- as.data.table(cbind(Dd18.pwt, Dd18.magnitude, this.phase, 1/this.sensitivity))
colnames(df.Dd18) <- c('Dd18.pwt', 'Dd18.magnitude', 'phase', 'd18.sens')
# get phase anomaly relative to no change
# (error calculated as precipitation weighted Dd18O when in-phase versus out of phase)
df.Dd18$Dd18.pwt_minusZeroPhase <- abs(rep(df.Dd18[phase == 0]$Dd18.pwt, length(phases)) - df.Dd18$Dd18.pwt)

# get Dd18.pwt phase sensitivity
phase.effect <- phase.effect.nrmlzd <- this.sens <- vector()
for(i in 1:length(unique(df.Dd18$d18.sens))){ # loop through d18 sensitivities to get phase sensitivity at each step
  this.sens[i] <- unique(df.Dd18$d18.sens)[i]
  tdf <- df.Dd18[d18.sens == this.sens[i]]
  # calculate phase effect on Dd18
  phase.effect[i] <- max(tdf$Dd18.pwt) - min(tdf$Dd18.pwt)
  phase.effect.nrmlzd[i] <- phase.effect[i] / mean(tdf$Dd18.magnitude, na.rm=T)
}
# bring together
sum.df <- as.data.table(cbind(this.sens, phase.effect, phase.effect.nrmlzd))
colnames(sum.df) <- c('d18.sensitivity', 'phase.effect', 'phase.effect_Dd18.magnitude')


## [PLOT RESULTS]
# identify 0, pi/2, pi, phases
p.0 <- phases[1]
p.pi_2 <- phases[which.min(abs(phases - (pi/2)))]
p.pi <- phases[which.min(abs(phases - (pi)))]
txt.size <- 14; titl.size <- 17

## [EXAMPLE]
## ZERO PHASE DIFFERENCE
THIS.DIFF <- 'zero' ; this.phase <- p.0
p.1 <- ggplot(outdf[phase == this.phase & d18.sensitivity == d18.sensitivity[1]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
  geom_line(aes(x=tx, y=upwind.P), color='lightblue', size=2) + 
  geom_line(aes(x=tx, y=downwind.P), color='darkblue', size=2, linetype='dashed') +
  scale_x_continuous(name='Time of year', expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5','1')) +
  scale_y_continuous(name='Precipitation\nrate') +
  theme_few() +
  theme(axis.text = element_text(size=14), axis.title=element_text(size=17))
p.1
save.name <- paste('Precip_', THIS.DIFF, '.png', sep = '')
# ggsave(save.name, p.1, width=20, height=5, units='cm')

p.2 <- ggplot(outdf[phase == this.phase & d18.sensitivity == d18.sensitivity[1]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-12.5, 0)) +
  geom_line(aes(x=tx, y=upwind.d18), color='lightblue', size=2) + 
  geom_line(aes(x=tx, y=downwind.d18), color='darkblue', size=2) +
  geom_line(aes(x=tx, y=Dd18), color='red', size=1, linetype='dashed') +
  scale_x_continuous(name='Time of year', expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5','1')) +
  scale_y_continuous(name=expression(delta^"18"*"O")) +
  theme_few() +
  theme(axis.text = element_text(size=txt.size), axis.title=element_text(size=titl.size))
p.2
save.name <- paste('d18_', THIS.DIFF, '.png', sep = '')
# ggsave(save.name, p.2, width=20, height=8, units='cm')

## PI OVER 2 PHASE DIFFERENCE
THIS.DIFF <- 'pi_2' ; this.phase <- p.pi_2
p.1 <- ggplot(outdf[phase == this.phase & d18.sensitivity == d18.sensitivity[1]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
  geom_line(aes(x=tx, y=upwind.P), color='lightblue', size=2) + 
  geom_line(aes(x=tx, y=downwind.P), color='darkblue', size=2) +
  scale_x_continuous(name='Time of year', expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5','1')) +
  scale_y_continuous(name='Precipitation\nrate') +
  theme_few() +
  theme(axis.text = element_text(size=14), axis.title=element_text(size=17))
p.1
save.name <- paste('Precip_', THIS.DIFF, '.png', sep = '')
# ggsave(save.name, p.1, width=20, height=5, units='cm')

p.2 <- ggplot(outdf[phase == this.phase & d18.sensitivity == d18.sensitivity[1]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-12.5,0)) +
  geom_line(aes(x=tx, y=upwind.d18), color='lightblue', size=2) + 
  geom_line(aes(x=tx, y=downwind.d18), color='darkblue', size=2) +
  geom_line(aes(x=tx, y=Dd18), color='red', size=1, linetype='dashed') +
  scale_x_continuous(name='Time of year', expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5','1')) +
  scale_y_continuous(name=expression(delta^"18"*"O")) +
  theme_few() +
  theme(axis.text = element_text(size=txt.size), axis.title=element_text(size=titl.size))
p.2
save.name <- paste('d18_', THIS.DIFF, '.png', sep = '')
# ggsave(save.name, p.2, width=20, height=8, units='cm')


## PI PHASE DIFFERENCE
THIS.DIFF <- 'pi' ; this.phase <- p.pi
p.1 <- ggplot(outdf[phase == this.phase & d18.sensitivity == d18.sensitivity[1]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
  geom_line(aes(x=tx, y=upwind.P), color='lightblue', size=2) + 
  geom_line(aes(x=tx, y=downwind.P), color='darkblue', size=2) +
  scale_x_continuous(name='Time of year', expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5','1')) +
  scale_y_continuous(name='Precipitation\nrate') +
  theme_few() +
  theme(axis.text = element_text(size=14), axis.title=element_text(size=17))
p.1
save.name <- paste('Precip_', THIS.DIFF, '.png', sep = '')
# ggsave(save.name, p.1, width=20, height=5, units='cm')

p.2 <- ggplot(outdf[phase == this.phase & d18.sensitivity == d18.sensitivity[1]]) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-12.5,0)) +
  geom_line(aes(x=tx, y=upwind.d18), color='lightblue', size=2) + 
  geom_line(aes(x=tx, y=downwind.d18), color='darkblue', size=2) +
  geom_line(aes(x=tx, y=Dd18), color='red', size=1, linetype='dashed') +
  scale_x_continuous(name='Time of year', expand=c(0,0), breaks=c(0,0.5,1), labels=c('0', '0.5','1')) +
  scale_y_continuous(name=expression(delta^"18"*"O")) +
  theme_few() +
  theme(axis.text = element_text(size=txt.size), axis.title=element_text(size=titl.size))
p.2
save.name <- paste('d18_', THIS.DIFF, '.png', sep = '')
# ggsave(save.name, p.2, width=20, height=8, units='cm')




## [SUMMARY]
p.a <- ggplot(df.Dd18) + 
  coord_cartesian(xlim = c(0, 2*pi), ylim = c(1e-3, 0.5e-1)) +
  geom_line(aes(x=phase, y=Dd18.pwt_minusZeroPhase, color=d18.sens, group=d18.sens), size=2) +
  scale_color_cmocean(name='rain', trans='log10') + 
  scale_x_continuous(name='Phase difference (radians)', expand=c(0,0),
                     breaks=c(0, pi/2, pi, 3*pi/2, 2*pi), labels=c(0, expression(pi/'2'),
                                                                   expression(pi), expression('3'*pi/'2'),
                                                                   expression('2'*pi))) + 
  scale_y_continuous(name=expression('| '*Delta*delta^"18"*'O error |  (???/1000 km)')) + 
  theme_few() +
  theme(axis.text = element_text(size=txt.size), axis.title=element_text(size=titl.size),
        legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
p.a

save.name <- 'phase_err.png'
# ggsave(save.name, p.a, width=20, height=13, units='cm')



p.b <- ggplot(sum.df) + 
  coord_cartesian(ylim = c(1e-3, 0.5e-1)) +
  geom_point(aes(x=d18.sensitivity, y=phase.effect, fill=d18.sensitivity), shape=21, size=4) +
  scale_fill_cmocean(name='rain', trans='log10', guide='none') + 
  scale_x_log10() +
  scale_y_continuous(name=expression('| '*Delta*delta^"18"*'O error |  (???/1000 km)')) + 
  theme_few() +
  theme(axis.text = element_text(size=txt.size), axis.title=element_text(size=titl.size))
p.b
save.name <- 'phase_errMax.png'
ggsave(save.name, p.b, width=8, height=13, units='cm')


p.c <- ggplot(sum.df) + 
  coord_cartesian(ylim = c(1e-3, 0.5e-1)) +
  geom_point(aes(x=d18.sensitivity, y=phase.effect_Dd18.magnitude, fill=d18.sensitivity), shape=21, size=4) +
  scale_fill_cmocean(name='rain', trans='log10', guide='none') + 
  scale_x_log10() +
  scale_y_continuous(name=expression('| '*Delta*delta^"18"*'O error |  (???/1000 km)')) + 
  theme_few()+
  theme(axis.text = element_text(size=txt.size), axis.title=element_text(size=titl.size))
p.c
save.name <- 'phase_errMax_overSzn.png'
ggsave(save.name, p.c, width=8, height=13, units='cm')

