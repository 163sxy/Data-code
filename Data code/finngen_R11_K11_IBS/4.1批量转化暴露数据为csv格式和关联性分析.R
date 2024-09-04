

rm(list=ls())

if(!require("pacman")) install.packages("pacman",update = F,ask = F)

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
library("pacman")
p_load(VariantAnnotation,gwasglue,dplyr,tidyr,CMplot)
p_load_gh("mrcieu/gwasglue")

library(data.table)
library(dplyr)

data_path <- "data04"

output_path <- "newdata04/"

vcf_files <- list.files(path = data_path, pattern = "*.tsv.gz$", full.names = TRUE)
vcf_files


if (!dir.exists(output_path)) {
  dir.create(output_path)
}

for (infile in vcf_files) {
  outputname0 <- basename(gsub("_buildGRCh37.tsv.gz$", "", infile))
  
  outputname <- paste0(output_path, outputname0, ".csv")
  
  expo_data_MR <- vroom::vroom(infile) 
  
  expo_data_MR1 =expo_data_MR %>% 
    dplyr::filter(p_value < 1e-5)
  
  fwrite(expo_data_MR1, outputname)
}

