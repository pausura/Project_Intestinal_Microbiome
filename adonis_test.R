
    ## Adonis test ##

getwd()
setwd("/Users/paulasuredahorrach/Documents/Universitat/Holanda/Projecte/")

library(vegan)

#Open taxa data
  tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
      
#Open phenotypes data
  factor_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
 
  my_results <- matrix(ncol = 3, nrow=ncol(factor_table))       
        
  library(vegan)
  
#For each column in factor_table (factor)  
for (i in 1:ncol(factor_table)) {
  
  #Create a table for complete cases
   final_factor_table <- factor_table[complete.cases(factor_table[,i]),]
   
   filter_tax_table <- tax_table[rownames(final_factor_table),]
    
        ad <- adonis(formula = filter_tax_table ~ final_factor_table[,i] , data = final_factor_table, permutations = 1000, method = "bray")
        aov_table <- ad$aov.tab
        
        my_results[i,1]=aov_table[1,1]
        my_results[i,2]=aov_table[1,5]
        my_results[i,3]=aov_table[1,6]
        
    }

  
  rownames(my_results) = colnames(factor_table)
  colnames(my_results) = c("Df", "R2", "Pr(>F)")
  
  write.table(my_results, file="~/adonis_results.txt", quote=F, sep = "\t")
  
  
  ### Count nomber of NA (per factor)

#NA sum function
  nasum <- function(x){
    sum (is.na(x))
  }

  na_results <- matrix(ncol = 1, nrow=ncol(factor_table)) 
  
  for (i in 1:ncol(factor_table)) {
    
    na_sum = nasum(factor_table[,i])
    na_results[i,1] = na_sum
    
  }
  
  rownames(na_results) = colnames(factor_table) 
  colnames(na_results) = "N_NA"
   
  na_results <- as.data.frame(na_results)
  

    #Merge with adonis_results table
 
#Open adonis_results
  adonis_test <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/results_adonis.txt", sep = "\t", header = T, row.names = 1, check.names = F)
      colnames(adonis_test) = c("Df", "R2", "Pr(>F)")
  
  adonis_results <- merge(adonis_test, na_results, by="row.names")
    rownames(adonis_results) <- adonis_results[,1]
    adonis_results <- adonis_results[,-1]
  

    ### P-value correction
  
  library(stats)
    
  p_correction_fdr <- as.data.frame(p.adjust(adonis_results$`Pr(>F)`, method = "fdr"))
    rownames(p_correction_fdr) <- rownames(adonis_results)
  
  p_correction_bonferroni <- as.data.frame(p.adjust(adonis_results$`Pr(>F)`, method = "bonferroni"))
    rownames(p_correction_bonferroni) <- rownames(adonis_results)
    
    
  #Merge with adonis_results table
    
  final_adonis_results <- merge(adonis_results, p_correction_fdr, by="row.names")
    rownames(final_adonis_results) <- final_adonis_results[,1]
    final_adonis_results <- final_adonis_results[,-1]
  
  final_adonis_results <- merge(final_adonis_results, p_correction_bonferroni, by="row.names")
    rownames(final_adonis_results) <- final_adonis_results[,1]
    final_adonis_results <- final_adonis_results[,-1]
    
  #Change colnames
    colnames(final_adonis_results)[5] <- "FDR_p_value"
    colnames(final_adonis_results)[6] <- "Bonferroni_p_value"
    
  
  
    write.table(final_adonis_results, file="~/final_adonis_results.txt", sep= "\t", quote = F)
    
    
    