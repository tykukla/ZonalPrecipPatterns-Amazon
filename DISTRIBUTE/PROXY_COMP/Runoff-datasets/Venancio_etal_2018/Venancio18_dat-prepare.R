# ------------------------------------ # 
# bring together Venancio et al 2018   #
# data                                 #
# ------------------------------------ #

rm(list=ls())

# read in data
dat.path <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data/Venancio_etal_2018'
fn.dat <- 'Venancio18_TiCa-dat.csv'

# bring in 
df.dat <- as.data.table(read.csv(paste(dat.path, fn.dat, sep='/')))
# (age is already in data)

colnames(df.dat) <- c('depth_m', 'age_kabp', 'Ti.Ca')

ggplot(df.dat) + 
  geom_line(aes(x=age_kabp, y=log(Ti.Ca))) +
  scale_x_continuous(limits=c(0, 25)) +
  theme_few()


# -------------------------------------- # 
# SAVE RESULT
save.here <- dat.path
save.name <- 'Venancio17_TiCa.RDS'
# saveRDS(df.dat, paste(save.here, save.name, sep='/'))
# -------------------------------------- # 