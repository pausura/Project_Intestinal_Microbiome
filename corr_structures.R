
  ##Correlation test - correlation between phenotypes

library(psych)

metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep ="\t", header = T, row.names = 1, check.names = F)

  numeric_metadata <- metadata

  ## Convert the categorical factors in numeric
for (i in 1:ncol(numeric_metadata)){
  
  if (is.numeric(numeric_metadata[,i])){
      numeric_metadata[,i] <- numeric_metadata[,i]
     }

  else {
    numeric_metadata[,i] <- as.numeric(factor(numeric_metadata[,i]))
  }
}

  ##Calculate the correlation
ct_phenotypes <- corr.test(numeric_metadata, method = "spearman", adjust = "fdr")

    #Correlation table
corr_metadata <- as.matrix(ct_phenotypes$r)
      corr_metadata[is.na(corr_metadata)] <- 0
    #P_values table
fdr_corr_metadata <- as.matrix(ct_phenotypes$p)
      fdr_corr_metadata[is.na(fdr_corr_metadata)] <- 1

    
      ##Subset factors with correlations values > 0.8
library(dplyr)
library(reshape2)

  corr_metadata_melt <- arrange(melt(corr_metadata), -abs(value))
  fdr_metadata_melt <- arrange(melt(fdr_corr_metadata))
     
      corr_valid_values <- subset(corr_metadata_melt, corr_metadata_melt$value >= 0.70 | corr_metadata_melt$value <= -0.70)
      fdr_valid_values <- subset(fdr_metadata_melt, fdr_metadata_melt$value < 0.1)
      
    
  write.table(corr_valid_values, file="~/corr_valid_values.txt", sep = "\t", quote = F)
  write.table(fdr_valid_values, file="~/fdr_valid_values.txt", sep = "\t", quote= F)
  

    ## Create a table with correlations > 0.7 and p_value < 0.1
  correlation_factors = matrix(ncol=3, nrow=nrow(corr_metadata))
  
  for (x in 1:nrow(fdr_corr_metadata)) {
    
    for (y in 1:ncol(fdr_corr_metadata)) {
      
      if (fdr_corr_metadata[x,y]<0.1) {
        
        if (corr_metadata[x,y]>0.7 | corr_metadata[x,y]< -0.7)  {
       
           correlation_factors[x,1] = colnames(corr_metadata)[y]
           correlation_factors[x,2] = rownames(corr_metadata)[x]
           correlation_factors[x,3] = corr_metadata[x,y]
    
        }
      }
    }
  }
  
  correlation_factors <- as.data.frame(correlation_factors)

    ## Make sure that the columns are vectors (and not factors)
  for (i in 1:ncol(correlation_factors)){
    correlation_factors[,i]=as.vector(correlation_factors[,i])
  }
      ## Keep in a variable if the rows are equal or not
  same_values <- as.data.frame(ifelse (correlation_factors$V1==correlation_factors$V2, 1, 0))
      ## Create a new table with the new variable
  new_correlation_factors<-as.data.frame(cbind(correlation_factors,same_values))
      colnames(new_correlation_factors)[4] <- "V4"
  
      ## Remove that rows that have the same value in the two columns (=1 condition)
  new_correlation_factors <- new_correlation_factors[!(new_correlation_factors$V4=="1"),]
  
  write.table(new_correlation_factors, file = "~/final_correlation_factors.txt", sep = "\t", quote= F)
  