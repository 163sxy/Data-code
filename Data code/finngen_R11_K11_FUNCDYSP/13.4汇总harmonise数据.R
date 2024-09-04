 rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 library(pacman)
 p_load(data.table,dplyr,tidyr,ggplot2)

  cwd <- getwd()
cwd <- paste0(cwd,"/Result13")

 subdirs <- list.dirs(cwd, recursive = FALSE, full.names = TRUE)

 harmonise_files <- list()

 for (subdir in subdirs) {
   harmonise_file <- file.path(subdir, "harmonise.csv")
  
   if (file.exists(harmonise_file)) {
     harmonise_data <- read.csv(harmonise_file)
    
     harmonise_data$id <- basename(subdir)
    
    
     harmonise_files[[length(harmonise_files) + 1]] <- harmonise_data
  }
}

 harmonise_data_combined <- do.call(rbind, harmonise_files)

 print(harmonise_data_combined)
 library(data.table)
fwrite(harmonise_data_combined, file = paste0("allharmonise13.csv"),sep = ",", quote = F)

 