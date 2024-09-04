
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/

if (!require("devtools")) {
  install.packages("devtools")
} else {}

if (!require("data.table")) {
  install.packages("data.table")
} else {}

if (!require("TwoSampleMR")) {
  devtools::install_github("MRCIEU/TwoSampleMR") 
} else {}

library(pacman)
p_load(data.table,dplyr,tidyr,ggplot2)

library(ieugwasr)
library(MRInstruments)
library(plyr)
library(dplyr)
library(data.table)
library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(purrr)


if(!dir.exists("Result13")){ 
  dir.create("Result13") 
}


if(!dir.exists("ORdata13")){
  dir.create("ORdata13")
}

if(!dir.exists("Pleiotropydata13")){
  dir.create("Pleiotropydata13")
}


if(!dir.exists("heterogeneity13")){
  dir.create("heterogeneity13")
}

if(!dir.exists("PRESSO13")){
  dir.create("PRESSO13")
}

if(!dir.exists("harmonise13")){
  dir.create("harmonise13")
}


exposure_dat=fread("FLDexposure12.csv")
workdir <- getwd() 
workdir <- paste0(getwd(),"/out")
workdir



filterdata=fread("filter_ivw_MRresult09.csv")
filterdataid=filterdata$id
outcomefile=paste0("./out/",filterdataid,".csv.gz")


print(outcomefile)

result=c()
for(i in outcomefile){
  tryCatch({  
    outdata=fread(i)
    outdata$PHENO <- gsub("\\.csv.gz$", "", basename(i))
    colnames(outdata)
    outdata=as.data.frame(outdata)
    outc_data <- format_data(
      outdata,
      type="outcome",
      snps=exposure_dat$SNP,
      header=T,
      phenotype_col="PHENO",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "standard_error", 
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "p_value",
      eaf_col = "effect_allele_frequency",
      chr_col = "chromosome",
      pos_col = "base_pair_location"
    )

    dat <- harmonise_data(exposure_dat, outc_data, action = 1)
  
    dat<-subset(dat,mr_keep==TRUE) 
     res <- mr(dat)
     result_or <- generate_odds_ratios(res)
       if (T) {
         
        
      filename <- gsub(".csv.gz$", "", basename(i))
        
       filename2 <- sub("\\.csv$", "", paste0(filename))
      dir.create(paste0("./Result13/",filename2))
      
      
       write.table(dat, 
                  file = paste0("./Result13/",filename2,"/harmonise.csv"),
                  row.names = F, sep = ",", quote = F)
      
       write.table(dat, 
                  file = paste0("./harmonise13/",filename, "_harmonise.csv"), 
                  sep = ",", quote = F, row.names=F)
      
      
      
      write.table(result_or, 
                  file = paste0("./Result13/",filename2,"/OR.csv"),
                  row.names = F, sep = ",", quote = F)
      
       write.table(result_or, 
                  file = paste0("./ORdata13/",filename, "_OR.csv"), 
                  row.names = FALSE, sep = ",", quote = F)
       p1 <- mr_scatter_plot(res, dat)
      ggsave(p1[[1]], file=paste0("./Result13/",filename2,"/scatter.pdf"), 
             width=8, height=8)
      
       pleiotropy <- mr_pleiotropy_test(dat)
      write.table(pleiotropy, file = paste0("./Result13/",filename2,"/pleiotropy.csv"),
                  sep = ",", quote = F, row.names=F)
      
       write.table(pleiotropy, 
                  file = paste0("./Pleiotropydata13/",filename, "_pleiotropy.csv"), 
                  sep = ",", quote = F, row.names=F)
      
       heterogeneity <- mr_heterogeneity(dat)
      write.table(heterogeneity, file = paste0("./Result13/",filename2,"/heterogeneity.csv"),
                  sep = ",", quote = F, row.names=F)
      
       write.table(heterogeneity, 
                  file = paste0("./heterogeneity13/",filename2, "_heterogeneity.csv"),
                  sep = ",", quote = F, row.names=F)
      
      
      
       presso <- run_mr_presso(dat, NbDistribution = 1000)
      
       write.table(presso[[1]]$`Main MR results`, 
                  file = paste0("./PRESSO13/",filename2, "_mrPRESSO_Main_MR_results.csv"),
                  sep = ",", quote = F, row.names=F)
      
      write.table(presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue, 
                  file = paste0("./PRESSO13/",filename2, "_mrPRESSO_Global_Test_Pvalue.csv"),
                  sep = ",", quote = F, row.names=F)
      
       
      
      capture.output(presso, file = paste0("./Result13/",filename2,"/presso.csv"))
      
      
  
       
       singlesnp_res <- mr_singlesnp(dat)
      singlesnpOR <- generate_odds_ratios(singlesnp_res)
      write.table(singlesnpOR, file=paste0("./Result13/",filename2,"/singlesnpOR.csv"),
                  row.names = F, sep = ",", quote = F)
      
       p2 <- mr_forest_plot(singlesnp_res)
      ggsave(p2[[1]], file=paste0("./Result13/",filename2,"/forest.pdf"), width=8, height=8)
      
       sen_res <- mr_leaveoneout(dat)
      p3 <- mr_leaveoneout_plot(sen_res)
      ggsave(p3[[1]], file=paste0("./Result13/",filename2,"/sensitivity-analysis.pdf"), 
             width=8, height=8)
      
       res_single <- mr_singlesnp(dat)
      p4 <- mr_funnel_plot(singlesnp_res)
      ggsave(p4[[1]], file=paste0("./Result13/",filename2,"/funnelplot.pdf"), width=8, height=8)
     }
  }, silent = TRUE)  
}

 