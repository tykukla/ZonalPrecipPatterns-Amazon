# ----------------------------------------------------- # 
# Script to bring together Nace 2012 data (published in #
# Nace et al 2014) for Ti/Ca                            #
# ----------------------------------------------------- #

rm(list=ls())

# ... read in data
dat.path <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data/Nace_etal_2014'
fn <- 'CDH86-dist-R_ready.csv'
#fn1 <- 'age_model_for-R.csv'
#fn2 <- 'Nace_thesis_data-only-R.csv'
#df.ages <- as.data.table(read.csv(paste(dat.path, fn1, sep='/')))
df <- as.data.table(read.csv(paste(dat.path, fn, sep='/')))

# save result
# saveRDS(df, file=paste(dat.path, 'Nace14_dat.RDS', sep='/'))

# ... see how it looks
ggplot(df) + 
  geom_point(aes(x=Age..Cal.Yrs., y=Ti.Ca)) +
  scale_y_log10() +
  scale_x_continuous(limits=c(0, 25e3)) + 
  theme_few()



# ... apply age model
age.fun <- approxfun(df.ages$Depth_cm, df.ages$Foram_Cal_age_BP)
age.fun2 <- approxfun(df.ages$Depth_cm, df.ages$OrgC_Cal_age_BP)
df.x$foram_age <- age.fun(df.x$depth_cm)
df.x$orgC_age <- age.fun2(df.x$depth_cm)
df.x$age <- ifelse(is.na(df.x$orgC_age), df.x$foram_age, df.x$orgC_age)
# add NA where there's a gap
View(df.x[age > 17e3 & age < 20e3])

# ... see how it looks
ggplot(df.x) + 
  geom_line(aes(x=depth_cm, y=TiCa)) + 
  scale_y_log10() +
  theme_few()


ggplot(df.x) + 
  coord_cartesian(xlim=c(0, 75000)) + 
  geom_line(aes(x=orgC_age, y=log(TiCa)), color='orange') + 
  # scale_y_continuous(trans="log2") +
  theme_few()


ggplot(df.x) + 
  coord_cartesian(xlim=c(0, 75000)) + 
  geom_line(aes(x=age, y=TiCa), color='orange') + 
  scale_y_continuous(trans="log2") +
  theme_few()

ggplot(df.x) + 
  coord_cartesian(xlim=c(0, 75000)) + 
  geom_line(aes(x=age, y=log(TiCa)), color='orange') + 
  theme_few()


