 
rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)

 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 

 library("pacman")

 p_load(forestploter,grid,tidyr,ggplot2,data.table)

filename <- "allOR13.csv"  
inputfile <- "exposureInformation.csv"

biofsci <- fread(filename)
 mydata <- fread(inputfile)
 biof_imm <- merge(biofsci,mydata,by="id")

fwrite(biof_imm,"allOR_Inform14.csv")
ivw_data <- subset(biof_imm, method == "Inverse variance weighted")

fwrite(ivw_data,"ivwOR_Inform14.csv")

 