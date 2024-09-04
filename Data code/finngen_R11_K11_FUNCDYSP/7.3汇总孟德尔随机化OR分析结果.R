rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
if(!require("pacman")) install.packages("pacman",update = F,ask = F)

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 

library("pacman")
p_load(data.table, dplyr) 

df_list <- list()

files_path <- "ORdata07/" 
files <- list.files(path = files_path, pattern = "*_OR.csv") 
full_files <- paste0(files_path, files) 
full_files

for(f in full_files){
  df <- fread(f)
  df$id <- f
  df$id <- gsub("_OR.csv", "", df$id)
  df$id <- gsub("ORdata07/", "", df$id)
  df_list <- c(df_list, list(df))
}
combined_df <- Reduce(function(x,y) rbind(x,y), df_list)
fwrite(combined_df,"allOR07.csv")

ivw_data <- combined_df %>% 
  filter(method=="Inverse variance weighted")
fwrite(ivw_data,"ivw_allOR07.csv")

