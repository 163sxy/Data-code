getwd()
setwd("D:/Users/ASUS/Desktop/finnibs/functional/finngen_R11_K11_FUNCDYSP")
library("pacman")


p_load(VariantAnnotation,gwasglue,dplyr,tidyr,CMplot)

p_load_gh("mrcieu/gwasglue")
library("VariantAnnotation")

library(data.table)
library("GenomeInfoDb")

library(vcfR)
infile="ieu-a-30.vcf.gz"                  
diseaseName="ieucd"        




expo_data_MR <- readVcf(infile) %>% 
  gwasvcf_to_TwoSampleMR(type = "outcome")

expo_data_MR$PHENO <- diseaseName 
colnames(expo_data_MR)
expo_data_MR <- expo_data_MR %>%
  rename(
    PHENO=PHENO,
    SNP = SNP,
    CHR = chr.outcome,
    BP = pos.outcome,
    effect_allele = effect_allele.outcome,
    other_allele = other_allele.outcome,
    P = pval.outcome,
    EAF = eaf.outcome,
    BETA = beta.outcome,
    SE = se.outcome,
    samplesize = samplesize.outcome
  ) 
colnames(expo_data_MR)
expo_data_MR <- expo_data_MR %>%
  select(PHENO,SNP, CHR, BP, effect_allele, other_allele, P, EAF, BETA, SE,samplesize) %>%
  mutate(P = as.numeric(P)) 


write.csv(expo_data_MR, "ieucdout.csv", row.names = FALSE)


library(pacman)

p_load(data.table,dplyr,tidyr,ggplot2)

library(ieugwasr)
library(MRInstruments)
library(plyr)
library(dplyr)
library(data.table)
library(TwoSampleMR)



if(!dir.exists("ORdata07")){
  dir.create("ORdata07")
}


pfilter=0.05

if(!dir.exists("ORdata07")){
  dir.create("ORdata07")
}

if(!dir.exists("Pleiotropydata07")){
  dir.create("Pleiotropydata07")
}

if(!dir.exists("Result07")){
  dir.create("Result07")
}


if(!dir.exists("heterogeneity07")){
  dir.create("heterogeneity07")
}


if(!dir.exists("PRESSO07")){
  dir.create("PRESSO07")
}

if(!dir.exists("harmonise07")){
  dir.create("harmonise07")
}



FileNames <- list.files(path="C:/Users/ASUS/Desktop/finnibs/7finnibs/7/exposure/", pattern = "*.csv")
FileNames
outcomefile="finngen_R11_K11_FUNCDYSP"

exp_dat <- list() 
ex_pore <- c()
for(i in c(1:length(FileNames))){
  IV <- fread(paste0("C:/Users/ASUS/Desktop/finnibs/7finnibs/7/exposure","/",FileNames[i]))
  IV$PHENO <- FileNames[i] 
  IV=as.data.frame(IV)
  IV1<-format_data(IV,
                   type="exposure",
                   phenotype_col = "PHENO", 
                   snp_col = "rsids2",
                   beta_col = "beta",
                   se_col = "standard_error",
                   eaf_col = "effect_allele_frequency",
                   effect_allele_col = "effect_allele",
                   other_allele_col = "other_allele",
                   pval_col = "p_value",
                   samplesize_col = "total_sample_size",
                   chr_col = "chromosome",
                   pos_col = "base_pair_location")
  exp_dat[[i]] <- IV1 
  
  ex_pore<-c(ex_pore,FileNames[i])
}

GWAS_1 <- fread(outcomefile)
allSNP <- do.call(rbind, exp_dat)
GWAS_2 <- subset(GWAS_1,GWAS_1$rsids %in% allSNP$SNP)  #GWAS_1的SNP肯可能要修改
rm(GWAS_1)

GWAS_2$PHENO<-"fd" 
head(GWAS_2)
GWAS_2 <- as.data.frame(GWAS_2)
out_data <- format_data(
  GWAS_2,
  type="outcome",
  phenotype_col = "PHENO",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  # eaf_col = "EAF",
  samplesize_col = "samplesize",
  pval_col = "pval",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  chr_col = "#chrom",
  pos_col = "pos"
  )

out_dat <- list() 
out_dat[[1]] <- out_data
out_come<-c("ieucd") 

results <- list()
for (i in seq_along(ex_pore)) {
  for (j in seq_along(out_come)) {
    tryCatch({
      dat <- harmonise_data(
        exposure_dat = exp_dat[[i]],
        outcome_dat = out_dat[[j]],
        action = 1
      )
      dat<-subset(dat,mr_keep==TRUE) 
      
      dat$R2 <- (2 * (dat$beta.exposure^2)) /
        (2 * (dat$beta.exposure^2) +
           2 * dat$samplesize.exposure * dat$se.exposure^2)
      dat$f <- dat$R2 * (dat$samplesize.exposure - 2) / (1 - dat$R2)
      dat$meanf<- mean( dat$f)
      dat<-dat[dat$f>10,]
      res <- mr(dat)
      result_or <- generate_odds_ratios(res)

      if (!is.na(result_or$pval[3]) && result_or$pval[3] < pfilter) {
        filename <- basename(sub("\\.txt$","",ex_pore[i])) 
        filename2 <- sub("\\.csv$", "", paste0(filename))
        dir.create(paste0("./Result07/",filename2))
        write.table(dat, 
                    file = paste0("./Result07/",filename2,"/harmonise.csv"),
                    row.names = F, sep = ",", quote = F)
        
        write.table(dat, 
                    file = paste0("./harmonise07/",filename2, "_harmonise.csv"), 
                    sep = ",", quote = F, row.names=F)
        
        
        write.table(result_or, 
                    file = paste0("./Result07/",filename2,"/OR.csv"),
                    row.names = F, sep = ",", quote = F)
        
        write.table(result_or, 
                    file = paste0("./ORdata07/",filename2, "_OR.csv"), 
                    row.names = FALSE, sep = ",", quote = F)
        p1 <- mr_scatter_plot(res, dat)
        ggsave(p1[[1]], file=paste0("./Result07/",filename2,"/scatter.pdf"), 
               width=8, height=8)
        
        pleiotropy <- mr_pleiotropy_test(dat)
        write.table(pleiotropy, file = paste0("./Result07/",filename2,"/pleiotropy.csv"),
                    sep = ",", quote = F, row.names=F)
        
        write.table(pleiotropy, 
                    file = paste0("./Pleiotropydata07/",filename2, "_pleiotropy.csv"), 
                    sep = ",", quote = F, row.names=F)
        
        
        heterogeneity <- mr_heterogeneity(dat)
        write.table(heterogeneity, file = paste0("./Result07/",filename2,"/heterogeneity.csv"),
                    sep = ",", quote = F, row.names=F)
        
        write.table(heterogeneity, 
                    file = paste0("./heterogeneity07/",filename2, "_heterogeneity.csv"),
                    sep = ",", quote = F, row.names=F)
        
        
        
        presso <- run_mr_presso(dat, NbDistribution = 1000)
        capture.output(presso, file = paste0("./Result07/",filename2,"/presso.csv"))
        

        write.table(presso[[1]]$`Main MR results`, 
                    file = paste0("./PRESSO07/",filename2, "_mrPRESSO_Main_MR_results.csv"),
                    sep = ",", quote = F, row.names=F)
        
        write.table(presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, 
                    file = paste0("./PRESSO07/",filename2, "_mrPRESSO_Global_Test_Pvalue.csv"),
                    sep = ",", quote = F, row.names=F)
        
        singlesnp_res <- mr_singlesnp(dat)
        singlesnpOR <- generate_odds_ratios(singlesnp_res)
        write.table(singlesnpOR, file=paste0("./Result07/",filename2,"/singlesnpOR.csv"),
                    row.names = F, sep = ",", quote = F)
        
        p2 <- mr_forest_plot(singlesnp_res)
        ggsave(p2[[1]], file=paste0("./Result07/",filename2,"/forest.pdf"), width=8, height=8)
        
        sen_res <- mr_leaveoneout(dat)
        p3 <- mr_leaveoneout_plot(sen_res)
        ggsave(p3[[1]], file=paste0("./Result07/",filename2,"/sensitivity-analysis.pdf"), 
               width=8, height=8)
        
        res_single <- mr_singlesnp(dat)
        p4 <- mr_funnel_plot(singlesnp_res)
        ggsave(p4[[1]], file=paste0("./Result07/",filename2,"/funnelplot.pdf"), width=8, height=8)
        res$exposure=ex_pore[i]
        res$outcome=out_come[j]
        

        results[[length(out_come)*(i-1)+j]] <- generate_odds_ratios(res)
      }
    }, error = function(e) {
      cat("Error occurred for file:", i, "\n")
      cat("Error message:", conditionMessage(e), "\n")
    })
  }
}
results_allIV <- do.call(rbind, results) 
fwrite(results_allIV,"result07.csv")


