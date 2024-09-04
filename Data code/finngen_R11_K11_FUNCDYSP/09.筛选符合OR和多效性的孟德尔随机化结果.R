
library("pacman")

p_load(data.table,ggplot2,dplyr) 

mr_data <- fread("allOR07.csv")
ple_data <- fread("allpleiotropy07.csv")

ivwpfilter = 0.05

ivw <- data.frame()

for (biofsci in unique(mr_data$id)) {
  biofdata <- mr_data[mr_data$id == biofsci,]
  if (biofdata[biofdata$method == "Inverse variance weighted", "pval"] < ivwpfilter) {
    if (all(biofdata$or > 1) || all(biofdata$or < 1)) {
      ivw <- rbind(ivw, biofdata)
    }
  }
}

ple_data <- ple_data[ple_data$pval > 0.05,]
immuneLists <- as.vector(ple_data$id)
immuneLists <- sub("_pleiotropy.csv", "", immuneLists)
ivw$id <- sub("_OR.csv", "", ivw$id)
pleout <- ivw[ivw$id %in% immuneLists,]
fwrite(pleout,"filter_MRresult09.csv")
pleout$method
ivw <- ivw[ivw$method == "Inverse variance weighted" & ivw$pval < ivwpfilter,]
fwrite(ivw,"filter_ivw_MRresult09.csv")


