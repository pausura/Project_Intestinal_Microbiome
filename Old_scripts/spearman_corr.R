
  ## Spearman Correlation ##

getwd()
setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

  ## Open species table

  species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
 
  library(corrplot)

      my_species_corr <- cor(species_table)
      corrplot(my_species_corr, method="color", order = "hclust", tl.pos = "n")
      
      
      
      
      
      
      
      
      
      
      
      
      