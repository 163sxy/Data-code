
library(dplyr)

directory <- "PRESSO13"
suffix <- "_mrPRESSO_Main_MR_results.csv"

files <- list.files(directory, pattern = paste0(suffix, "$"))
files 
combined_df <- data.frame()

for (file in files) {
  file_path <- file.path(directory, file)
  
  df <- read.csv(file_path)
  df <- df[1,]
  
  id_value <- sub(suffix, "", file)


  df$id <- id_value
  df$id <- gsub("_buildGRCh37.tsv","",df$id)
  combined_df <- bind_rows(combined_df, df)
}

write.csv(combined_df, "allmrPRESSO_Main_MR_results13.csv", row.names = FALSE)

print(combined_df)

