
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
library("pacman")
p_load(data.table,tidyr,dplyr,vroom) 
p_load_gh("mrcieu/gwasglue")


library(TwoSampleMR)

MRFile="exposure12.csv"     
kbfilter=10000
r2filter=0.001
biofsci=fread(MRFile)
head(biofsci)

biofsci <- read_exposure_data(MRFile, 
                              sep = ",",      
                              snp_col = "SNP",    
                              beta_col = "BETA", 
                              se_col = "SE",   
                              effect_allele_col = "effect_allele", 
                              other_allele_col = "other_allele", 
                              pval_col = "P",
                              eaf_col = "EAF", 
                              #samplesize_col = "samplesize", 
                              chr_col = "CHR",
                              pos_col = "BP",
                              clump=F, 
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
                            # plink_bin = "./plink_win64_20231211/plink.exe"  
)

expo_data <- expo_data[which(expo_data$SNP %in% clump_dat$rsid),] 

write.table(expo_data, "LDexposure12.csv", row.names = FALSE, sep = ",", quote = FALSE)


