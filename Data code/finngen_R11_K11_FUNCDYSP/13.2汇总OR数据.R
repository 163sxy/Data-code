 
rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 library(pacman)
 p_load(data.table,dplyr,tidyr,ggplot2)


 cwd0 <- getwd()
cwd <- paste0(cwd0,"/Result13")

 subdirs <- list.dirs(cwd, recursive = FALSE, full.names = TRUE)

 or_files <- list()

 for (subdir in subdirs) {
   or_file <- file.path(subdir, "OR.csv")
  
   if (file.exists(or_file)) {
     or_data <- read.csv(or_file)
    
     or_data$id <- basename(subdir)
    or_data$id <- gsub("_buildGRCh37.tsv","",or_data$id)
    
     or_files[[length(or_files) + 1]] <- or_data
  }
}

 or_data_combined <- do.call(rbind, or_files)


fwrite(or_data_combined , file = paste0("allOR13.csv"),sep = ",", quote = F)


 

 