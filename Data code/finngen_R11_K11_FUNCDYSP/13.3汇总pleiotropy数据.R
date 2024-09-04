  

rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 library(pacman)
 p_load(data.table,dplyr,tidyr,ggplot2)

 cwd0 <- getwd()
cwd <- paste0(cwd0,"/Result13")

 subdirs <- list.dirs(cwd, recursive = FALSE, full.names = TRUE)

 pleiotropy_files <- list()

 for (subdir in subdirs) {
   pleiotropy_file <- file.path(subdir, "pleiotropy.csv")
  
   if (file.exists(pleiotropy_file)) {
     pleiotropy_data <- read.csv(pleiotropy_file)
    
     pleiotropy_data$id <- basename(subdir)
    
    pleiotropy_data$id <- gsub("_buildGRCh37.tsv","",pleiotropy_data$id)
     pleiotropy_files[[length(pleiotropy_files) + 1]] <- pleiotropy_data
  }
}

 pleiotropy_data_combined <- do.call(rbind, pleiotropy_files)

 print(pleiotropy_data_combined)
 library(data.table)
fwrite(pleiotropy_data_combined, file = paste0("allpleiotropy13.csv"),sep = ",", quote = F)

 