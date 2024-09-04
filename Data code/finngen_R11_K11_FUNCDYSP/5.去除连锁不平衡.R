#######################################################
######课程来源：生信私学                           ####
######零基础学生信，就到生信私学                   ####
######保姆式教学，一对一指导                       ####
######生信私学，让您出类拔萃                       ####
######关注微信公众号：回复 M40 获得脚本            ####
######脚本经常更新，请和我们联系更换               ####
#######################################################

rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
# 加载pacman包
library("pacman")
# 在加载这些包之前，“pacman”会先检查是否已经安装，如果没有会自动安装
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("VariantAnnotation",force = TRUE)
p_load(VariantAnnotation,gwasglue,dplyr,tidyr,CMplot)
p_load_gh("mrcieu/gwasglue")

library(data.table)
library(dplyr)
library(TwoSampleMR)

#关联性分析 1e-5 5e-8
#pfilter <- 1e-05

#去除连锁不平衡的SNP参数
kbfilter=10000
r2filter=0.001

# 设置路径到包含VCF文件的data文件夹
data_path <- "data05"

# 设置目标路径到输出CSV文件的newdata文件夹
output_path <- "newdata05/"

#data=fread("data05/ieu-a-1073.csv" )
# 查找data文件夹内所有以.csv结尾的文件
vcf_files <- list.files(path = data_path, pattern = "*.csv$", full.names = TRUE)
#vcf_files=vcf_files[1]
#vcf_files=vcf_files[16:46]
vcf_files
# 确保新文件夹newdata存在，如果不存在则创建
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

#data =fread("data05/GCST90026011.csv")
#colnames(data)

# 读取每个csv文件，转换它们，并写入CSV文件到newdata文件夹
for (infile in vcf_files) {
  # 提取文件名，不包括路径和扩展名
  outputname0 <- basename(gsub(".csv$", "", infile))
  # 创建CSV文件名
  outputname <- paste0(output_path, outputname0, ".csv")
  # 读取输入文件
  biofsci <- read_exposure_data(infile, # 读取文件名存储在MRFile中的CSV文件
                                sep = ",",       # 列分隔符为逗号
                                #phenotype_col = "phenotype", # 暴露表型列名
                                snp_col = "SNP",     # SNP ID列名
                                beta_col = "beta", # 暴露效应大小列名
                                se_col = "standard_error",   # 暴露效应标准误列名
                                effect_allele_col = "effect_allele", # 暴露效应等位基因列名
                                other_allele_col = "other_allele", # 暴露其他等位基因列名
                                pval_col =  "p_value",
                                eaf_col = "effect_allele_frequency", # 暴露等位基因频率列名
                                samplesize_col = "total_sample_size", # 暴露样本大小列名
                                chr_col = "chromosome",
                                pos_col = "base_pair_location",
                                clump=F # 进行连锁不平衡剔除
  )  
  
  #biofsci1 =biofsci %>%  dplyr::filter(pval.exposure< pfilter)
  
  # 在线clump方法
  # expo_data1 <- clump_data(biofsci, clump_kb=kbfilter,clump_r2 = r2filter,clump_p1 = 1,clump_p2 = 1,pop = "EUR")
  # 将结果写入CSV文件
  # fwrite(expo_data1, outputname)}


  
  # 备用本地clump方法
  expo_data <- biofsci %>%
    distinct(SNP, .keep_all=TRUE)  # 根据SNP列进行去重，保留所有列
  expo_data <- expo_data[expo_data$SNP != "", ]  # 过滤掉SNP为空的行
  #expo_data1 <- expo_data[expo_data$pval.exposure< 5e-8, ]  # 根据p值筛选出p值小于1e-6的行，这里可以根据需求进行修改
  colnames(expo_data)
  #本地clump
  biof_iv <- expo_data[,c("SNP","pval.exposure")]  #修改
  # 将biof_iv的数据框中的列名修改为rsid和pval
  colnames(biof_iv) <- c("rsid","pval")
  
  # 使用ld_clump函数进行SNP筛选,得到独立的SNP
  clump_dat <- ld_clump_local(dat = biof_iv,  # 进行局部LD聚类的数据集
                              clump_kb = kbfilter,  # LD聚类时使用的窗口大小（以基对为单位）
                              clump_r2 = r2filter,  # LD聚类时使用的相关系数阈值
                              clump_p = 1,  # LD聚类时使用的p值阈值
                              bfile = "./data_maf0.01_rs_ref/data_maf0.01_rs_ref",  # PLINK格式的参考数据集文件路径
                              plink_bin = "./plink_mac/plink"  # mac电脑用这句，windows电脑此句前加#
                              #plink_bin = "./plink_win64_20231211/plink.exe"  # windows电脑用这句，mac电脑此句前加#
                              )
  
  # 保留expo_data数据框中在clump_dat$rsid中的SNP                      
  expo_data1 <- expo_data[which(expo_data$SNP %in% clump_dat$rsid),] #注意不要少了逗号
  
  
  # 将结果写入CSV文件
  fwrite(expo_data1, outputname)
}

#######################################################
######课程来源：生信私学                           ####
######零基础学生信，就到生信私学                   ####
######保姆式教学，一对一指导                       ####
######生信私学，让您出类拔萃                       ####
######关注微信公众号：回复 M40 获得脚本            ####
######脚本经常更新，请和我们联系更换               ####
#######################################################