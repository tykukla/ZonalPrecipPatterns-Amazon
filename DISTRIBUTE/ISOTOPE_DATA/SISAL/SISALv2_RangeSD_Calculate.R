# ------------------------------------------------------- # 
# Script to pull out the SISAL database range/sd data     #
# specifically for the Amazon basin LGM to modern         #
# (and comparison to global)                              #
# ---                                                     #
# T Kukla (Stanford Univ. 2018)                           # 
# ------------------------------------------------------- #
library(ggplot2)
library(ggmap)
library(dplyr)
# library(gifski)
# library(ggrepel)
library(data.table)
library(maps)


rm(list=ls())

# directory where DISTRIBUTE is located
user.dir <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data'
# shouldn't need to change
mainPath <- 'DISTRIBUTE/ISOTOPE_DATA/SISAL/sisalv2_csv'
setwd(paste(user.dir, mainPath, sep='/'))

# ----------------------------------------------------------------------------- 
# read in data (one at a time)
# ----------------------------------------------------------------------------- 

#... Isotopes
d13 <- as_tibble(read.csv('d13C.csv'))
d18 <- as_tibble(read.csv('d18O.csv'))

#... Sites 
site <- as_tibble(read.csv('site.csv'))

#... Entity
entity <- as_tibble(read.csv('entity.csv'))

#... Chronology
ages <- as_tibble(read.csv('original_chronology.csv'))

#... Samples
sampledf <- as_tibble(read.csv('sample.csv'))


# ----------------------------------------------------------------------------- 
# bring it all together
# ----------------------------------------------------------------------------- 

d13d18 <- merge(d13, d18, by=c('sample_id'), all=TRUE)

sample_ent <- merge(d13d18, sampledf, by=c('sample_id'), all=TRUE)

ent <- merge(sample_ent, entity, by=c('entity_id'), all=TRUE)

sitedf <- merge(ent, site, by=c('site_id'), all=TRUE)

df <- as_tibble(merge(sitedf, ages, by=c('sample_id'), all=TRUE))
df <- as.data.table(df)
df <- df[d18O == "yes"]

# ... also bring in composite DTP curve to repeat analysis on that
dtp.path <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/IsotopeData'
dtp.name <- 'd18_Western_notSmoothed_cleaned.RDS'
df.dtp <- as.data.table(readRDS(paste(dtp.path, dtp.name, sep='/')))
df.dtp$site_id <- df.dtp$latitude <- df.dtp$longitude <- -9999
colnames(df.dtp) <- c("interp_age", "d18O_measurement", "longitude", "latitude", "site_id")
df.dtp.in <- df.dtp[ , c("d18O_measurement", "interp_age", "site_id", "latitude", "longitude")]
# get correct age units
df.dtp.in$interp_age <- df.dtp.in$interp_age * 1e3
# only look at last 25 kyr
df.dtp.in <- df.dtp.in[interp_age < 25e3]

# ----------------------------------------------------------------------------- 
# SMOOTH AND GET STATISTICS
# ----------------------------------------------------------------------------- 
df.in <- df[ , c("d18O_measurement", "interp_age", "site_id", "latitude", "longitude")]
# add dtp record
df.in <- as.data.table(rbind(df.in, df.dtp.in))
these.ids <- unique(df.in$site_id)
avg.bin <- 1e3 # [years to avg in each bin]
avg.sample.threshold <- 5  # need at least this many samples, otherwise bin is NA
record.dur.min <- 1e4  # [years] minimum duration of record to calculate its long-term d18O mean

for(j in 1:length(these.ids)){
  tmp.df <- df.in[site_id == these.ids[j]]
  min.age <- min(tmp.df$interp_age, na.rm=T)
  max.age <- max(tmp.df$interp_age, na.rm=T)
  if(is.infinite(min.age) | is.infinite(max.age)){
    avg.out[i] <- age.out[i] <- NA
    next
  } else{
    # get running avg
    these.avgs <- seq(min.age, max.age, by=avg.bin)
    avg.out <- age.out <- sd.out <- range.out <- vector()
    for(i in 1:length(these.avgs)){
      avg.df <- tmp.df[interp_age > these.avgs[i] & interp_age < these.avgs[i+1]]
      if(nrow(avg.df[!is.na(d18O_measurement)]) < avg.sample.threshold){ # if below threshold, assign NA
        avg.out[i] <- age.out[i] <- sd.out[i] <- range.out[i] <- NA
      } else{ # otherwise do calculations
        avg.out[i] <- mean(avg.df$d18O_measurement, na.rm = T)
        age.out[i] <- mean(c(these.avgs[i], these.avgs[i+1]), na.rm=T)
        sd.out[i] <- sd(avg.df$d18O_measurement, na.rm = T)
        range.out[i] <- range(avg.df$d18O_measurement, na.rm = T)[2]-range(avg.df$d18O_measurement, na.rm = T)[1]
      }
    }
  }
  
  if(j==1){
    outdf <- as.data.table(cbind(tmp.df$site_id[1], tmp.df$longitude[1], tmp.df$latitude[1], avg.out, age.out, sd.out, range.out))
    colnames(outdf) <- c("site_id", "longitude", "latitude", "d18_mean", "age_ka", "kyr_sd", "kyr_range")
  } else{
    outdf.tmp <- as.data.table(cbind(tmp.df$site_id[1], tmp.df$longitude[1], tmp.df$latitude[1], avg.out, age.out, sd.out, range.out))
    colnames(outdf.tmp) <- c("site_id", "longitude", "latitude", "d18_mean", "age_ka", "kyr_sd", "kyr_range")
    outdf <- as.data.table(rbind(outdf, outdf.tmp))
  }
  # ... print progress
  print(paste("solving", j+1, "of", length(these.ids), sep=' '))
}

# remove records that are too short, calculate st.dev and range
# first remove d18O and ages that are na
outdf <- outdf[!is.na(d18_mean) & !is.na(age_ka)]
stdev_out <- min.out <- max.out <- range_out <- vector()
this.lat <- this.lon <- this.site <- this.age.mean <- vector()
age.range <- kyr.range.mean <- kyr.sd.mean <- vector()
these.ids2 <- unique(outdf$site_id)
idx <- 1
for(k in 1:length(these.ids2)){
  tmp.df <- outdf[site_id == these.ids2[k]]
  # check if duration is long enough 
  if(nrow(tmp.df) >= (record.dur.min/avg.bin)){ # get summary statistics
    stdev_out[idx] <- sd(tmp.df$d18_mean, na.rm=T)
    min.out[idx] <- range(tmp.df$d18_mean, na.rm=T)[1]
    max.out[idx] <- range(tmp.df$d18_mean, na.rm=T)[2]
    range_out[idx] <- max.out[idx] - min.out[idx]
    this.lat[idx] <- tmp.df$latitude[1]
    this.lon[idx] <- tmp.df$longitude[1]
    this.site[idx] <- tmp.df$site_id[1]
    this.age.mean[idx] <- mean(tmp.df$age_ka, na.rm = T)
    age.range[idx] <- range(tmp.df$age_ka, na.rm=T)[2] - range(tmp.df$age_ka, na.rm=T)[1]
    kyr.range.mean[idx] <- range(tmp.df$kyr_range, na.rm=T)[2] - range(tmp.df$kyr_range, na.rm=T)[1]
    kyr.sd.mean[idx] <- sd(tmp.df$kyr_sd, na.rm=T)
    # up the counter
    idx <- idx + 1 
  } else{
    next
  }
  
}

# bring together
df.sum <- as.data.table(cbind(this.lat, this.lon, this.site, this.age.mean, age.range,
                              range_out, min.out, max.out, stdev_out, kyr.range.mean, kyr.sd.mean))
colnames(df.sum) <- c("lat", "lon", "site_id", "age_mean", "age_range", "d18_range", "d18.min", "d18.max", "d18.sd",
                      "kyr_range_mean", "kyr_sd_mean")


# find the amazon data points
myLims <- c(minlat = -40, maxlat = 40, maxlon = -30, minlon = -108) 
df.sum.amz <- df.sum[lat > myLims['minlat'] & lat < myLims['maxlat'] & lon > myLims['minlon'] & lon < myLims['maxlon']]
# RGN
ggplot(df[site_id == 111]) + 
  geom_line(aes(x=interp_age, y=d18O_measurement)) 
# PAR 
ggplot(df[site_id == 3]) + 
  geom_line(aes(x=interp_age, y=d18O_measurement)) 
# DTP 
ggplot() + 
  geom_line(data = df[site_id %in% c(235, 226)], 
            aes(x=interp_age, y=d18O_measurement, group=site_id, color=site_id)) 


amz.sites <- c(111, 3, 235, 226, -9999)
df.sum$amz.site <- ifelse(df.sum$site_id %in% amz.sites, TRUE, FALSE)
df.sum.amz <- df.sum
df.sum <- df.sum[amz.site == FALSE]
df.sum.noDTP <- df.sum.amz[!(site_id %in% amz.sites[5])]


# ------------------------------------------------------------ # 
# SAVE RESULTS ----------------------------------------------- # 
OUTLIST <- list(df.sum, df.sum.noDTP, df.sum.amz) 
names(OUTLIST) <- c("summary_no-amazon", "summary_no-added-DTP", "summary_all")
save.here <- mainPath
save.name <- 'SISAL_Smooth_10kyrMin-1kyrAvg.RDS'
# saveRDS(OUTLIST, paste(save.here, save.name, sep='/'))
# ------------------------------------------------------------ # 
# ------------------------------------------------------------ # 




