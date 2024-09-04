 
rm(list=ls())
 if(!require("pacman")) install.packages("pacman",update = F,ask = F)
 options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #https://mirrors.pku.edu.cn/CRAN/
 
 

 library(reshape2)
library(circlize)
library(ComplexHeatmap)

mrFile="filter_MRresult15.csv"       

 biofsci=read.csv(mrFile, header=T, sep=",", check.names=F)
biofsci$method[biofsci$method=="Inverse variance weighted"]="IVW"

 mat=acast(biofsci, id ~ method, value.var="pval")
mat
mat[is.na(mat)] <- 0
mat

 col_color = colorRamp2(c(0, 0.5, 1), c("midnightblue", "gold", "firebrick"))#绘制圈图

circos.clear()
circos.par(gap.after=c(30))
pdf(file="circos.pdf", width=8, height=8)

circos.heatmap(mat, 
               col = col_color,
               dend.side = "inside",          
               rownames.side = "outside",     
               bg.border = "black")           

 cn = colnames(mat)      
n = length(cn)
circos.text(rep(CELL_META$cell.xlim[2], n) + 
              convert_x(1, "mm"), 1.2+(n:1)*0.7,      
            cn, cex = 0.6,       
            adj = c(0, 0.5),
            facing = "inside")
circos.clear()

 lgd = Legend(title="Pvalue", col_fun=col_color)
grid.draw(lgd)
dev.off()

 