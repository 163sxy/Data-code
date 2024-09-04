

rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
library("pacman")
p_load(data.table,tidyr,dplyr,vroom) 

infile="finngen_R11_K11_IBS.gz"                     
diseaseName="IBS"        

outputname = "finngen_R11_K11_IBS.csv"

samplesize = 372135

out_data_MR = vroom(infile) 
head(out_data_MR)
out_data_MR$PHENO <- diseaseName 
out_data_MR$samplesize=samplesize

out_data_MR1 <- out_data_MR %>%
  dplyr::rename(
    PHENO=PHENO,
    SNP = rsids,
    CHR = "#chrom",
    BP = pos,
    effect_allele = alt,
    other_allele = ref,
    P = pval,
    EAF = af_alt,
    BETA = beta,
    SE = sebeta,
    samplesize = samplesize
  ) 

out_data_MR2 <- out_data_MR1 %>%
  select(PHENO,SNP, CHR, BP, effect_allele, other_allele, P, EAF, BETA, SE,samplesize) %>%
  mutate(P = as.numeric(P)) 


fwrite(out_data_MR2,outputname)

