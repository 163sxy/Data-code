
library(dplyr)

directory <- "PRESSO07"
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
  
  combined_df <- bind_rows(combined_df, df)
}

write.csv(combined_df, "allmrPRESSO_Main_MR_results07.csv", row.names = FALSE)

print(combined_df)



