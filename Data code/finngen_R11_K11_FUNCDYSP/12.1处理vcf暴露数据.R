
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
library("pacman")
p_load(VariantAnnotation,gwasglue,dplyr,tidyr,CMplot)
p_load_gh("mrcieu/gwasglue")


infile="ieu-b-42.vcf.gz" 

expo_data_MR <- readVcf(infile) %>% 
  gwasvcf_to_TwoSampleMR(type = "exposrue") 
head(expo_data_MR)
colnames(expo_data_MR)
expo_data_MR1 <- expo_data_MR %>%
  rename(
    SNP = SNP,
    CHR = chr.exposrue,
    BP = pos.exposrue,
    effect_allele = effect_allele.exposrue,
    other_allele = other_allele.exposrue,
    P = pval.exposrue,
    EAF = eaf.exposrue,
    BETA = beta.exposrue,
    SE = se.exposrue,
    samplesize = samplesize.exposrue
  ) 

expo_data_MR2 <- expo_data_MR1 %>%
  select(SNP, CHR, BP, effect_allele, other_allele, P, EAF, BETA, SE,samplesize) %>%
  mutate(P = as.numeric(P))%>% 
  filter(P<5e-8) 
write.csv(expo_data_MR2, "exposure.csv", row.names = FALSE)
