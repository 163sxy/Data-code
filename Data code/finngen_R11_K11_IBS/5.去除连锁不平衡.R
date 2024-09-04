
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
library("pacman")
p_load(VariantAnnotation,gwasglue,dplyr,tidyr,CMplot)
p_load_gh("mrcieu/gwasglue")

library(data.table)
library(dplyr)
library(TwoSampleMR)


kbfilter=10000
r2filter=0.001

data_path <- "data05"

output_path <- "newdata05/"

vcf_files <- list.files(path = data_path, pattern = "*.csv$", full.names = TRUE)
vcf_files
if (!dir.exists(output_path)) {
  dir.create(output_path)
}


for (infile in vcf_files) {
  outputname0 <- basename(gsub(".csv$", "", infile))
  outputname <- paste0(output_path, outputname0, ".csv")
  biofsci <- read_exposure_data(infile,
                                sep = ",",       
                                #phenotype_col = "phenotype", 
                                snp_col = "SNP",     
                                beta_col = "beta", 
                                se_col = "standard_error",   
                                effect_allele_col = "effect_allele", 
                                other_allele_col = "other_allele", 
                                pval_col =  "p_value",
                                eaf_col = "effect_allele_frequency", 
                                samplesize_col = "total_sample_size", 
                                chr_col = "chromosome",
                                pos_col = "base_pair_location",
                                clump=F 
  )  
  



  
  expo_data <- biofsci %>%
    distinct(SNP, .keep_all=TRUE)  
  expo_data <- expo_data[expo_data$SNP != "", ]  
  colnames(expo_data)
  biof_iv <- expo_data[,c("SNP","pval.exposure")]  
  colnames(biof_iv) <- c("rsid","pval")
  
  clump_dat <- ld_clump_local(dat = biof_iv,  
                              clump_kb = kbfilter,  
                              clump_r2 = r2filter,  
                              clump_p = 1,  
                              bfile = "./data_maf0.01_rs_ref/data_maf0.01_rs_ref",  
                              plink_bin = "./plink_mac/plink"  
                              
                              )
  
  expo_data1 <- expo_data[which(expo_data$SNP %in% clump_dat$rsid),] 
  
  
  fwrite(expo_data1, outputname)
}

