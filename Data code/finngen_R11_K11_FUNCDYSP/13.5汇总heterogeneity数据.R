 

rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 library(pacman)
 p_load(data.table,dplyr,tidyr,ggplot2)

 cwd0 <- getwd()
cwd <- paste0(cwd0,"/Result13")

 subdirs <- list.dirs(cwd, recursive = FALSE, full.names = TRUE)

 heterogeneity_files <- list()

 for (subdir in subdirs) {
   heterogeneity_file <- file.path(subdir, "heterogeneity.csv")
  
   if (file.exists(heterogeneity_file)) {
     heterogeneity_data <- read.csv(heterogeneity_file)
    
     heterogeneity_data$id <- basename(subdir)
    heterogeneity_data$id <- gsub("_buildGRCh37.tsv","",heterogeneity_data$id)
    
     heterogeneity_files[[length(heterogeneity_files) + 1]] <- heterogeneity_data
  }
}

 heterogeneity_data_combined <- do.call(rbind, heterogeneity_files)

 print(heterogeneity_data_combined)
 library(data.table)
fwrite(heterogeneity_data_combined, file = paste0("heterogeneity13.csv"),sep = ",", quote = F)
 