

library("pacman")

p_load(forestploter,grid,tidyr,ggplot2,data.table)

filename <- "filter_MRresult09.csv"  
inputfile <- "exposureInformation.csv"

biofsci <- fread(filename)
mydata <- fread(inputfile)
biof_imm <- merge(biofsci,mydata,by="id")

fwrite(biof_imm,"allOR_Inform10.csv")
ivw_data <- subset(biof_imm, method == "Inverse variance weighted")

fwrite(ivw_data,"ivwOR_Inform10.csv")

