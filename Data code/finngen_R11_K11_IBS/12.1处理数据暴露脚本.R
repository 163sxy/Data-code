rm(list=ls())
library(vroom)
library(tidyr)
library(dplyr)
library(data.table)

biofsci <- vroom('finngen_R10_M13_ANKYLOSPON_STRICT.gz', col_names = TRUE)
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
  filter(P < 5e-8) 
write.csv(biofsci2, "exposure12.csv", row.names = FALSE)



