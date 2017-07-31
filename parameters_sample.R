### Calculate parameters per sample

getwd()

setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

#Reads table
  reads_table <- read.table("~/Documents/Universitat/Holanda/Projecte/reads_over_10M.txt", sep = "\t", header = T)
     
#Species table
  species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", sep = "\t", header = T, row.names = 1)

        ## Calculate number of species per sample ##

  species_results=as.data.frame(matrix(ncol = 1, nrow=ncol(species_table)))
 
   ##Function to calculate nº of non-0
  nsum <- function(x){
    sum (x!=0)  }
  
  ## Loop
  
  for (i in 1:ncol(species_table)) {
  
    a = nsum(species_table[,i])
    
    species_results[i,1] = a
    
  }

colnames(species_results) = c("Nº of species")
rownames(species_results) = colnames(species_table)
  
sp_results <- write.table(species_results, file = "~/filtered_species_results_DUDes.txt", quote = F, sep = "\t")
  
  #Introduce rownames to merge
    species_results$IBDFEC <- rownames(species_results)
  
#Shannon diversity table
  diversity_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/alpha_diversity_DUDes.txt", sep = "\t", header = T, row.names = 1)

  #Introduce rownames to merge
    diversity_table$IBDFEC <- rownames(diversity_table)
 

  #### MERGE: nº of reads / nº of species / shannon diversity
  MyMerge <- function(x, y){
    df <- merge(x, y, by= "IBDFEC", all.x= TRUE, all.y= TRUE)
    return(df)
  }
 
 sample_table <- Reduce(MyMerge, list(diversity_table, species_results, reads_table))
  rownames(sample_table) <- sample_table[,1]
  sample_table <- sample_table[,-1]
 
 merge_results <- write.table(sample_table, file = "~/filtered_sample_parameters_DUDes.txt", quote = F, sep = "\t")
