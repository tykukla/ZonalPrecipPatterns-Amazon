# ************************************************************** # 
# Script to plot the MAP distribution from the monte carlo       # 
# results -- using the wettest 25% of solutions for the          #
# hydrostat scenario simulations                                 #
# -------------                                                  #
# T. Kukla (Stanford Univ. 2018)                                 #
# ************************************************************** # 
library(GGally)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(data.table)

rm(list=ls()) 
# directory where DISTRIBUTE is located
user.dir <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data'
# shouldn't need to change
mainPath <- 'DISTRIBUTE/REACTIVE_TRANSPORT'
setwd(paste(user.dir, mainPath, sep='/'))


# function to convert [kg m-2 s-1] to [mm yr-1]
flux_conv <- function(flux){
  secs_per_min <- 60
  min_per_hr <- 60
  hr_per_day <- 24
  day_per_yr <- 365

  secs_per_yr <- secs_per_min * min_per_hr * hr_per_day * day_per_yr

  flux * secs_per_yr
}

#... function to calculate slope
run_slope_CavePts <- function(df, cave1_kmloc=0.01e3, cave2_kmloc=3.5e3){
  cave1_dex <- which.min(abs(df$x - cave1_kmloc))
  cave2_dex <- which.min(abs(df$x - cave2_kmloc))
  cave_sep <- cave2_kmloc - cave1_kmloc
  (df$d18P[cave1_dex] - df$d18P[cave2_dex]) / (cave_sep / 1e3)
}

# ---- set cave locations
cave1 <- 10 ; cave2 <- 2400

 
# #... Read in the data 
# ## --------------- LGM --------------- ##
# LGMdir <- 'output'
# dfoutLGM <- readRDS(paste(LGMdir,'MCdf_7par_LGM_Mw2.RDS', sep='/'))
# outlist <- readRDS(paste(LGMdir, 'MClist_7par_LGM_Mw2.RDS', sep='/'))
# # re-calculate the gradient
# thisSlope <- thisP <- vector()
# for(i in 1:length(outlist)){
#   # calculate n=1 to n=2.4e3 gradient
#   thisSlope[i] <- run_slope_CavePts(df=outlist[[i]], cave1_kmloc = cave1, cave2_kmloc = cave2)
#   thisP[i] <- mean(outlist[[i]]$P[which(outlist[[i]]$x >= cave1 & outlist[[i]]$x <= cave2)])
#   print(i)
# }
# dfoutLGM$iso_Dd2 <- thisSlope
# dfoutLGM$pr_mean2 <- thisP
# # get rid of the list because it's a huge file
# rm(outlist)
# 
# # --- SAVE --- #
# # saveRDS(dfoutLGM, file='LGM_isoPr.RDS')

# ## --------------- Mid-Hol --------------- ##
# midHoldir <- 'output'
# dfoutMH <- readRDS(paste(midHoldir,'MCdf_7par_midHol_Mw2.RDS', sep='/'))
# outlist <- readRDS(paste(midHoldir, 'MClist_7par_midHol_Mw2.RDS', sep='/'))
# # re-calculate the gradient
# thisSlope <- thisP <- thisP_median <- vector()
# for(i in 1:length(outlist)){
#   # calculate n=1 to n=2.4e3 gradient 
#   thisSlope[i] <- run_slope_CavePts(df=outlist[[i]], cave1_kmloc = cave1, cave2_kmloc = cave2)
#   thisP[i] <- mean(outlist[[i]]$P[which(outlist[[i]]$x >= cave1 & outlist[[i]]$x <= cave2)])
#   thisP_median[i] <- median(outlist[[i]]$P[which(outlist[[i]]$x >= cave1 & outlist[[i]]$x <= cave2)])
#   print(i)
# }
# dfoutMH$iso_Dd2 <- thisSlope
# dfoutMH$pr_mean2 <- thisP
# dfoutMH$pr_median2 <- thisP_median
# # get rid of the list because it's a huge file 
# rm(outlist)
# 
# # --- SAVE --- #
# # saveRDS(dfoutMH, file='MH_isoPr.RDS')
# 
# 
# 
# ## --------------- Pre-Industrial --------------- ##
# PICdir <- 'output'
# dfoutPIC <- readRDS(paste(PICdir,'MCdf_7par_PIC_Mw2.RDS', sep='/'))
# outlist <- readRDS(paste(PICdir, 'MClist_7par_PIC_Mw2.RDS', sep='/'))
# # re-calculate the gradient
# thisSlope <- thisP <- vector()
# for(i in 1:length(outlist)){
#   # calculate n=1 to n=2.4e3 gradient 
#   thisSlope[i] <- run_slope_CavePts(df=outlist[[i]], cave1_kmloc = cave1, cave2_kmloc = cave2)
#   thisP[i] <- mean(outlist[[i]]$P[which(outlist[[i]]$x >= cave1 & outlist[[i]]$x <= cave2)])
#   print(i)
# }
# dfoutPIC$iso_Dd2 <- thisSlope
# dfoutPIC$pr_mean2 <- thisP
# # get rid of the list because it's a huge file 
# rm(outlist)
# 
# # --- SAVE --- #
# # saveRDS(dfoutPIC, file='PIC_isoPr.RDS')




# ------ READ IN THE ALREADY PROCESSED DATA (processing code is above) 
# [MID HOLOCENE]
dfoutMH <- readRDS('MH_isoPr.RDS')
# [PRE INDUSTRIAL]
dfoutPIC <- readRDS('PIC_isoPr.RDS')
# [LAST GLACIAL MAXIMUM]
dfoutLGM <- readRDS('LGM_isoPr.RDS')


# ****************************************************** #
## RAW MONTE CARLO DATA FILES ------------------------- ##
# ****************************************************** #
# Extract the simulations that fit the isotope gradient
## LGM
grad_min <- 2.0 ; grad_max <- 2.3
dfxLGM <- dfoutLGM[which(dfoutLGM$iso_Dd2 > grad_min & dfoutLGM$iso_Dd2 < grad_max),]

## MidHol
grad_minMH <- 0 ; grad_maxMH <- 0.3
dfxMH <- dfoutMH[which(dfoutMH$iso_Dd2 > grad_minMH & dfoutMH$iso_Dd2 < grad_maxMH),]
#... for the mid-holocene extract the wettest 25% of results 
myPercentile <- 0.75        # value 0-to-1 ; we will extract values above this percentile 
percVal <- quantile(dfxMH$pr_median2, probs=myPercentile)
dfxMH_perc <- dfxMH[which(dfxMH$pr_median2 >= percVal), ]

## PIC
grad_minPIC <- 0.83 ; grad_maxPIC <- 1.3
dfxPIC <- dfoutPIC[which(dfoutPIC$iso_Dd2 > grad_minPIC & dfoutPIC$iso_Dd2 < grad_maxPIC),]


# ****************************************************** #
## MODIFIED MONTE CARLO DATA FILE (LGM only) ---------- ##
# ****************************************************** #
#... set the min vals for tfrac and u-intercept 
min_tfrac <- 0.64         # modern tfrac
min_uInt <- 8.5            # modern wind int

dfxLGM_modify <- dfxLGM[which(dfxLGM$tfrac >= min_tfrac & dfxLGM$u_intercept > min_uInt), ]



# *************************************************************************************************************** #
# ---------------------------------------------- PDF FIGURE(S) -------------------------------------------------- #
# *************************************************************************************************************** #
# --- MANUSCRIPT SIZE 
#... arrow annotation info
LGMmean1 <- median(flux_conv(dfxLGM$pr_mean2))
LGMmean2 <- median(flux_conv(dfxLGM_modify$pr_mean2))
PImean <- median(flux_conv(dfxPIC$pr_mean2))
MHmean <- median(flux_conv(dfxMH_perc$pr_mean2))

yDens <- 6.5e-4

#... fill colors 
LGMcol <- '#6DBCDB' 
wetLGMcol <- '#3485CC'
MHcol <- '#FC4349' 
PICcol <- '#7F8282'
arrowCol <- '#1B466B'
#... segment colors
PICsegcol <- '#787A7A'
LGMsegcol <- '#528EA6'
wetLGMsegcol <- '#28659C'
MHsegcol <- '#CC363B'
#... y limits of figures
LGM_ylim <- 1.2
MH_ylim <- 2.4


# re-organize the data
# .. remove last column of MH which is not shared by LGM or PI
dfxMHx <- dfxMH[,c(1:(ncol(dfxMH)-1))]
dfxMH_percx <- dfxMH_perc[,c(1:(ncol(dfxMH_perc)-1))]
# .. add idx labels to PI, LGM, MH
dfxMHx$timeslice <- dfxMH_percx$timeslice <-  'MH'
dfxMH_percx$subtimeslice <- 'MH'
dfxLGM$timeslice <- dfxLGM_modify$timeslice <- 'LGM'
dfxLGM$subtimeslice <- 'LGM1' ; dfxLGM_modify$subtimeslice <- 'LGM2'
dfxPIC$timeslice <- 'PIC'
dfxPIC$subtimeslice <- 'PIC'

alldf <- as.data.table(rbind(dfxPIC, dfxMH_percx, dfxLGM, dfxLGM_modify))


# --- 
# where to save
saveDir <- '/Users/tylerkukla/Box Sync/Amazon_ResilienceTEMPORARY/figPr_iso'

# -------- BOXPLOT -------------
pr.px <- ggplot(alldf) + 
  coord_cartesian(ylim=c(0.5,4.8)) +
  # [ -- OPTION 1 -- ] Violin plot with quantiles
  # geom_violin(aes(x=timeslice, y=flux_conv(pr_mean2)/1e3, group=subtimeslice, fill=subtimeslice),
  #             draw_quantiles = c(0.25, 0.5, 0.75)) +
  # [ -- OPTION 2 -- ] Boxplot 
  # geom_boxplot(aes(x=timeslice, y=flux_conv(pr_mean2)/1e3, group=subtimeslice, fill=subtimeslice), width=0.5) + 
  # [ -- OPTION 3 -- ] Boxplot + Violin plot
  geom_violin(aes(x=timeslice, y=flux_conv(pr_mean2)/1e3, group=subtimeslice, fill=subtimeslice),
              position=position_dodge(1), width=1, alpha=0.4) +
  geom_boxplot(aes(x=timeslice, y=flux_conv(pr_mean2)/1e3, group=subtimeslice, fill=subtimeslice), width=0.2,
               position=position_dodge(1)) +
  # other plot stuffs
  scale_fill_manual(values=c('LGM1'=LGMcol, 'LGM2'=wetLGMcol, 'MH'=MHcol, 'PIC'=PICcol), guide='none') +
  scale_y_continuous(name="Annual precipitation (m)") +
  scale_x_discrete(labels=c("LGM" = "LGM", "MH" = "MH", "PIC" = "LH")) +
  theme_few() +
  theme(plot.background=element_rect(fill="transparent", color="transparent"),
        panel.background=element_rect(fill="transparent", color="transparent"),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title.y = element_text(size=15))

pr.px

saveFile <- paste(saveDir, 'PrBoxplot+Violin_iso.png', sep='/')
# ggsave(filename=saveFile, plot=pr.px, width=11, height=9, units=c('cm'), bg='transparent')



# 

