
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
library("pacman")
p_load(data.table,tidyr,dplyr,vroom) 

library(TwoSampleMR)
biofsci <- fread("LDexposure12.csv") 
biofsci$samplesize = 1381+409237
head(biofsci)
biofsci$R2<-(2*biofsci$beta*biofsci$beta*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)/(2*biofsci$beta*biofsci$beta*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)+2*biofsci$se.exposure*biofsci$se.exposure*biofsci$samplesize*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)))     #计算R2

biofsci$F<-biofsci$R2*(biofsci$samplesize-2)/(1-biofsci$R2)    

biofsci=biofsci[as.numeric(biofsci$"F")>10,]

write.csv(biofsci, file="FLDexposure12.csv", row.names=F)


