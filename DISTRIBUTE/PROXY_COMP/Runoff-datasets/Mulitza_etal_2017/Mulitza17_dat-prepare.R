# ------------------------------------ # 
# bring together Mulitza et al 2017    #
# data                                 #
# ------------------------------------ #

rm(list=ls())

# read in data
dat.path <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data/Mulitza_etal_2017'
fn.age <- 'GeoB16202-2_age_model-Rready.txt'
fn.dat <- 'GeoB16202-2_Fe_Ca_ratio-Rready.txt'

# bring in 
df.age <- as.data.table(read.delim(paste(dat.path, fn.age, sep='/')))
df.dat <- as.data.table(read.delim(paste(dat.path, fn.dat, sep='/')))
# (age is already in data)

colnames(df.dat) <- c('depth_m', 'age_kabp', 'Fe.Ca')

ggplot(df.dat) + 
  geom_line(aes(x=age_kabp, y=log(Fe.Ca))) + 
  theme_few()


# -------------------------------------- # 
# SAVE RESULT
save.here <- dat.path
save.name <- 'Mulitza17_FeCa.RDS'
# saveRDS(df.dat, paste(save.here, save.name, sep='/'))
# -------------------------------------- # 