rm(list=ls())
library(vroom)
library(tidyr)
library(dplyr)
library(data.table)
biofsci <- vroom('finngen_R11_K11_FUNCDYSP', col_names = TRUE)
head(biofsci)
biofsci1 <- biofsci %>%
  rename(
    SNP = rsids,
    CHR = "#chrom",
    BP = pos,
    effect_allele = alt,
    other_allele = ref,
    P = pval,
    EAF = af_alt,
    BETA = beta,
    SE = sebeta
  ) %>%
  select(SNP, CHR, BP, effect_allele, other_allele, P, EAF, BETA, SE) %>%
  mutate(P = as.numeric(P)) 
biofsci2 <- biofsci1 %>%
  filter(P < 5e-8) 
write.csv(biofsci2, "exposure12.csv", row.names = FALSE)
library("pacman")
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
biofsci2 <- biofsci1 %>%
  filter(P < 5e-6) 
write.csv(biofsci2, "exposure12.csv", row.names = FALSE)
MRFile="exposure12.csv"     
kbfilter=100
r2filter=0.01
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
                            # plink_bin = "./plink_mac/plink"  
                            plink_bin = "./plink_win64_20231211/plink.exe"  
)
clump_dat <- ld_clump_local(dat = biof_iv,  
                            clump_kb = kbfilter,  
                            clump_r2 = r2filter, 
                            clump_p = 1,  
                            bfile = "./data_maf0.01_rs_ref/data_maf0.01_rs_ref",  
                            # plink_bin = "./plink_mac/plink"  
                            plink_bin = "./plink_win64_20231211/plink.exe"  
)
kbfilter=1000
r2filter=0.01
clump_dat <- ld_clump_local(dat = biof_iv,  
                            clump_kb = kbfilter,  
                            clump_r2 = r2filter, 
                            clump_p = 1,  
                            bfile = "./data_maf0.01_rs_ref/data_maf0.01_rs_ref",  
                            # plink_bin = "./plink_mac/plink"  
                            plink_bin = "./plink_win64_20231211/plink.exe"  
)
expo_data <- expo_data[which(expo_data$SNP %in% clump_dat$rsid),] 
write.table(expo_data, "LDexposure12.csv", row.names = FALSE, sep = ",", quote = FALSE)
biofsci <- fread("LDexposure12.csv")
View(biofsci)
biofsci$samplesize = 395933
head(biofsci)
biofsci$R2<-(2*biofsci$beta*biofsci$beta*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)/(2*biofsci$beta*biofsci$beta*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)+2*biofsci$se.exposure*biofsci$se.exposure*biofsci$samplesize*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)))     #计算R2
biofsci$F<-biofsci$R2*(biofsci$samplesize-2)/(1-biofsci$R2)     #计算F检验值
biofsci=biofsci[as.numeric(biofsci$"F")>10,]
write.csv(biofsci, file="FLDexposure12.csv", row.names=F)
library(FastDownloader)
res=look_trait(file_name="FLDexposure12.csv",pval=1e-5)
library(FastTraitR)
FastDownloader::install_pkg("FastTraitR")
res=look_trait(file_name="FLDexposure12.csv",pval=1e-5)
library(FastTraitR)
res=look_trait(file_name="FLDexposure12.csv",pval=1e-5)
library(vroom)
library(tidyr)
library(dplyr)
library(data.table)
biofsci <- vroom('finngen_R11_K11_IBS', col_names = TRUE)
head(biofsci)
colnames(biofsci)
biofsci1 <- biofsci %>%
  rename(
    SNP = rsids,
    CHR = "#chrom",
    BP = pos,
    effect_allele = alt,
    other_allele = ref,
    P = pval,
    EAF = af_alt,
    BETA = beta,
    SE = sebeta
  ) %>%
  select(SNP, CHR, BP, effect_allele, other_allele, P, EAF, BETA, SE) %>%
  mutate(P = as.numeric(P)) 
biofsci2 <- biofsci1 %>%
  filter(P < 5e-6) 
write.csv(biofsci2, "exposure12.csv", row.names = FALSE)
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
                            # plink_bin = "./plink_mac/plink"  
                            plink_bin = "./plink_win64_20231211/plink.exe" 
)
expo_data <- expo_data[which(expo_data$SNP %in% clump_dat$rsid),] 
write.table(expo_data, "LDexposure12.csv", row.names = FALSE, sep = ",", quote = FALSE)
biofsci <- fread("LDexposure12.csv")
biofsci$samplesize = 372135
head(biofsci)
biofsci$R2<-(2*biofsci$beta*biofsci$beta*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)/(2*biofsci$beta*biofsci$beta*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)+2*biofsci$se.exposure*biofsci$se.exposure*biofsci$samplesize*biofsci$eaf.exposure*(1-biofsci$eaf.exposure)))     #计算R2
biofsci$F<-biofsci$R2*(biofsci$samplesize-2)/(1-biofsci$R2)     
biofsci=biofsci[as.numeric(biofsci$"F")>10,]
write.csv(biofsci, file="FLDexposure12.csv", row.names=F)
library(FastTraitR)
res=look_trait(file_name="FLDexposure12.csv",pval=1e-5)