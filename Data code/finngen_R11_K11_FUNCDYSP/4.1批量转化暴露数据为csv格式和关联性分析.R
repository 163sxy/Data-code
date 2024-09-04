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

# 设置路径到包含VCF文件的data文件夹
data_path <- "data04"

# 设置目标路径到输出CSV文件的newdata文件夹
output_path <- "newdata04/"

# 查找data文件夹内所有以.vcf.gz结尾的文件
vcf_files <- list.files(path = data_path, pattern = "*.tsv.gz$", full.names = TRUE)
#vcf_files=vcf_files[1]
vcf_files

#data = vroom::vroom("data04/GCST90002426_buildGRCh37.tsv")


# 确保新文件夹newdata存在，如果不存在则创建
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

# 读取每个VCF文件，转换它们，并写入CSV文件到newdata文件夹
for (infile in vcf_files) {
  # 提取文件名，不包括路径和扩展名
  outputname0 <- basename(gsub("_buildGRCh37.tsv.gz$", "", infile))
  
  # 创建CSV文件名
  outputname <- paste0(output_path, outputname0, ".csv")
  
  # 读取VCF文件并进行转换
  expo_data_MR <- vroom::vroom(infile) 
  
  expo_data_MR1 =expo_data_MR %>% 
    dplyr::filter(p_value < 1e-5)
  
  # 将结果写入CSV文件
  fwrite(expo_data_MR1, outputname)
}

#######################################################
######课程来源：生信私学                           ####
######零基础学生信，就到生信私学                   ####
######保姆式教学，一对一指导                       ####
######生信私学，让您出类拔萃                       ####
######关注微信公众号：回复 M40 获得脚本            ####
######脚本经常更新，请和我们联系更换               ####
#######################################################