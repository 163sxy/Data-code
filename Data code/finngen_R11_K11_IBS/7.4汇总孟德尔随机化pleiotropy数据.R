
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
if(!require("pacman")) install.packages("pacman",update = F,ask = F)

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 

library("pacman")
p_load(data.table, dplyr) 

df_list <- list()

files_path <- "Pleiotropydata07/" 
files <- list.files(path = files_path, pattern = "*_pleiotropy.csv") 
full_files <- paste0(files_path, files) 
full_files

for(f in full_files){
  df <- fread(f)
  df$id <- f
  df$id <- gsub("_pleiotropy.csv", "", df$id)
  df$id <- gsub("Pleiotropydata07/", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
fwrite(combined_df,"allpleiotropy07.csv")


