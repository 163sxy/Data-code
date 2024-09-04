 library(dplyr)

 directory <- "PRESSO13"
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
  df$id <- gsub("_buildGRCh37.tsv","",df$id)
  
   combined_df <- bind_rows(combined_df, df)
}

 write.csv(combined_df, "allmrPRESSO_Global_Test_Pvalue13.csv", row.names = FALSE)

 print(combined_df)

 