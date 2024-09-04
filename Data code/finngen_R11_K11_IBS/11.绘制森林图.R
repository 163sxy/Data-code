

rm(list=ls())
if(!require("pacman")) install.packages("pacman",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
library("pacman")
p_load(forestploter,grid,ggplot2,data.table)
mrresultsfile <- "allOR_Inform10.csv"
biofsci = fread(mrresultsfile, header = T)
colnames(biofsci)
biofsci$pval = as.numeric(biofsci$pval)
biofsci$pval <- ifelse(biofsci$pval<0.001, "<0.001", sprintf("%.4f", biofsci$pval))

biofsci$estimate <- paste0(format(round(biofsci$or, 4), nsmall = 4), " (",  
                           format(round(biofsci$or_lci95, 4), nsmall = 4), "-",
                           format(round(biofsci$or_uci95, 4), nsmall = 4), ")")

biofsci$Trails <-  biofsci$Trait


biofsci$Trails = ifelse(is.na(biofsci$Trails), "", biofsci$Trails)
biofsci$method = ifelse(is.na(biofsci$method), "", biofsci$method)
biofsci$nsnp = ifelse(is.na(biofsci$nsnp), "", biofsci$nsnp)  
biofsci$pval = ifelse(is.na(biofsci$pval), "", biofsci$pval)

colnames(biofsci)
biofsci$` ` <- paste(rep(" ", 15), collapse = " ")   
biofsci$`OR(95%CI)` <- biofsci$estimate

colnames(biofsci)
biofsci$Trails[duplicated(biofsci$Trails)] <- " "
biofsci1=biofsci[,c("Trails","method", "nsnp","pval","OR(95%CI)"," ")]
colnames(biofsci1)

tm <- forest_theme(base_size = 10,
                   refline_col = "red",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

colnames(biofsci1)
pdf("allforestall11.pdf", width=10, height=20)
p=forest(biofsci1,  
         est = biofsci$or,
         lower = biofsci$or_lci95,
         upper = biofsci$or_uci95,
         #arrow_lab = c("Placebo Better", "Treatment Better"),
         sizes = 0.4,
         ci_column = 6,
         ref_line = 1,
         xlim = c(0.05, 3),
         footnote = "",
         theme = tm)
p
dev.off()


