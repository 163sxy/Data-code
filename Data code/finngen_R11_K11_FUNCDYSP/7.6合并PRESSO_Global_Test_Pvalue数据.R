library(dplyr)

directory <- "PRESSO07"
suffix <- "_mrPRESSO_Global_Test_Pvalue.csv"

files <- list.files(directory, pattern = paste0(suffix, "$"))
files 
combined_df <- data.frame()

for (file in files) {
  file_path <- file.path(directory, file)
  
  df <- read.csv(file_path)
  
  id_value <- sub(suffix, "", file)
  colnames(df)[1]="Global_Test_Pvalue"
  
  df$Global_Test_Pvalue <- as.character(df$Global_Test_Pvalue)
  
  df$id <- id_value
  
  combined_df <- bind_rows(combined_df, df)
}

write.csv(combined_df, "allmrPRESSO_Global_Test_Pvalue07.csv", row.names = FALSE)

print(combined_df)
