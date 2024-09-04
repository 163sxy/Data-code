 
rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)

 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 

 library("pacman")

 p_load(data.table,ggplot2,dplyr) 

 mr_data <- fread("allOR_Inform14.csv")
ple_data <- fread("allpleiotropy13.csv")

 
ivwpfilter = 1

ivw <- data.frame()

 for (biofsci in unique(mr_data$id)) {
   biofdata <- mr_data[mr_data$id == biofsci,]
    if (biofdata[biofdata$method == "Inverse variance weighted", "pval"] < ivwpfilter) {
       if (T) {
       ivw <- rbind(ivw, biofdata)
    }
  }
}

 ple_data <- ple_data[ple_data$pval > 0.05,]
 immuneLists <- as.vector(ple_data$id)
 immuneLists <- sub("_pleiotropy.csv", "", immuneLists)
ivw$id <- sub("_OR.csv", "", ivw$id)
 pleout <- ivw[ivw$id %in% immuneLists,]
 fwrite(pleout,"filter_MRresult15.csv")

 ivw <- ivw[ivw$method == "Inverse variance weighted" & ivw$pval < ivwpfilter,]
 fwrite(ivw,"filter_ivw_MRresult15.csv")
 