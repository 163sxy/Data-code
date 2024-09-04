
rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/



library(reshape2)
library(circlize)
library(ComplexHeatmap)

mrFile="allOR07.csv"      


biofsci=read.csv(mrFile, header=T, sep=",", check.names=F)
biofsci$method[biofsci$method=="Inverse variance weighted"]="IVW"

mat=acast(biofsci, exposure ~ method, value.var="pval")
mat
mat[is.na(mat)] <- 0
mat

col_color = colorRamp2(c(0, 0.5, 1), c("midnightblue", "gold", "firebrick"))

circos.clear()
circos.par(gap.after=c(30))
pdf(file="circos.pdf", width=10, height=10)

circos.heatmap(mat, 
               col = col_color,
               dend.side = "inside",        
               rownames.side = "outside",   
               bg.border = "black")          

cn = colnames(mat)     
n = length(cn)
circos.text(rep(CELL_META$cell.xlim[2], n) + 
              convert_x(1, "mm"), 1.2+(n:1)*0.8,     
             cn, cex = 1,      
             adj = c(0, 0.5),
             facing = "inside")
circos.clear()

lgd = Legend(title="Pvalue", col_fun=col_color)
grid.draw(lgd)
dev.off()

