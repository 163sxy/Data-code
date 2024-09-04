
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
if(!require("pacman")) install.packages("pacman",update = F,ask = F)

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 

library("pacman")
p_load(data.table, dplyr) 

df_list <- list()

files_path <- "heterogeneity07/" 
files <- list.files(path = files_path, pattern = "*heterogeneity.csv") 
full_files <- paste0(files_path, files) 
full_files

for(f in full_files){
  df <- fread(f)
  df$id <- f
  df$id <- gsub(".heterogeneity.csv", "", df$id)
  df$id <- gsub("heterogeneity07/", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
fwrite(combined_df,"allheterogeneity07.csv")


