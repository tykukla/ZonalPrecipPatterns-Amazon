# ----------------------------------------------------- # 
# Script to bring together Nace 2012 data (published in #
# Nace et al 2014) for Ti/Ca                            #
# ----------------------------------------------------- #

rm(list=ls())

# ... read in data
dat.path <- 'C:/Users/tkukl/OneDrive/Documents/Amazon_ResilienceTEMPORARY/Data/Nace_etal_2014'
fn1 <- 'age_model_for-R.csv'
fn2 <- 'Nace_thesis_data-only-R.csv'
df.ages <- as.data.table(read.csv(paste(dat.path, fn1, sep='/')))
df <- as.data.table(read.csv(paste(dat.path, fn2, sep='/')))


# ... get df data in one set of columns
# split into diff dfs
df.1 <- df[ ,c('ï..depth_cm', 'TiCa_1', 'FeK_1')]
colnames(df.1) <- c('depth_cm', 'TiCa', 'FeK')
df.2 <- df[ ,c('depth_cm2', 'TiCa_2', 'FeK_2')]
colnames(df.2) <- c('depth_cm', 'TiCa', 'FeK')
df.3 <- df[ ,c('depth_cm3', 'TiCa_3', 'FeK_3')]
colnames(df.3) <-  c('depth_cm', 'TiCa', 'FeK')
df.4 <- df[ ,c('Depth_cm_4', 'TiCa4', 'FeK4')]
colnames(df.4) <-  c('depth_cm', 'TiCa', 'FeK')
df.5 <- df[ ,c('depth_cm5', 'TiCa5', 'FeK5')]
colnames(df.5) <- c('depth_cm', 'TiCa', 'FeK')

df.x <- as.data.table(rbind(df.1, df.2, df.3, df.4, df.5))


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


