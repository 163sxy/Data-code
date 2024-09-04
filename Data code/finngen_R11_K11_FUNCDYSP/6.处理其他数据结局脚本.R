#######################################################
######课程来源：生信私学                           ####
######零基础学生信，就到生信私学                   ####
######保姆式教学，一对一指导                       ####
######生信私学，让您出类拔萃                       ####
######关注微信公众号：biofsci 回复 M40 获得脚本    ####
######脚本经常更新，请和我们联系更换               ####
#######################################################

# 其他格式的数据如何处理
# 《34.M8孟德尔随机化暴露和结局数据处理与读取》

rm(list=ls())
# 判断是否已经安装了“pacman”包，如果没有就安装它
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
# 设置Bioconductor镜像地址为中国科技大学的镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
# 加载“pacman”包，用于方便加载其他的R包
# 加载pacman包
library("pacman")
# 批量安装加载需要的包
p_load(data.table,tidyr,dplyr,vroom) 
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("VariantAnnotation",force = TRUE)

# 输入文件
infile="finngen_R11_K11_FUNCDYSP.gz"                     # 这里定义了输入的VCF格式文件
diseaseName="FD"        #疾病名称

outputname = "finngen_R11_K11_FUNCDYSP.csv"

#样本量
samplesize = 395933

out_data_MR = vroom(infile) #读取数据 
head(out_data_MR)
out_data_MR$PHENO <- diseaseName #要修改
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
  mutate(P = as.numeric(P)) # 转换P值为数值型


fwrite(out_data_MR2,outputname)

#######################################################
######课程来源：生信私学                           ####
######零基础学生信，就到生信私学                   ####
######保姆式教学，一对一指导                       ####
######生信私学，让您出类拔萃                       ####
######关注微信公众号：biofsci 回复 M40 获得脚本    ####
######脚本经常更新，请和我们联系更换               ####
#######################################################
