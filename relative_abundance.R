## Calculate relative abundance

getwd()

setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

#Open data
tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1)
  
## Merge metadata with intestinal categories file     

intestinal_categories <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep="\t", header = T, row.names = 1)

total_data <- merge(intestinal_categories, tax_table, by="row.names")  
    rownames(total_data) <- total_data[,1]
    total_data <- total_data[,-1]
    
    
## Divide data in groups
    #Normal Group
        table_normal <- subset(total_data,total_data$Group=="normal")
            table_normal <- table_normal[,-1]
    #Intermediate Group
        table_intermediate <- subset(total_data, total_data$Group=="intermediate")
            table_intermediate <- table_intermediate[,-1]
    #Small Intestine Group
        table_smallintestine <- subset(total_data, total_data$Group=="small intestine")
            table_smallintestine <- table_smallintestine[,-1]

  total_data$Group <- NULL

             
##Function to calculate mean excluding 0 values
 nzmean <- function(x){
   mean(x[x!=0])
 }

##Function to calculate nº of 0
 nzsum <- function(x){
   sum (x==0)
 }

##Function to calculate nº of non-0
 nsum <- function(x){
   sum (x!=0)
 }

my_results=matrix(ncol = 4, nrow=ncol(table_smallintestine)) 
  
## Loop

for (i in 1:ncol(table_smallintestine)) {
  
  x = nzmean(table_smallintestine[,i])
  y = mean(table_smallintestine[,i])
  z = nzsum(table_smallintestine[,i])
  a = nsum(table_smallintestine[,i])
  
  my_results[i,1] = x
  my_results[i,2] = y
  my_results[i,3] = z
  my_results[i,4] = a

}

# The column names from the original table = row names from the new table
rownames(my_results) = colnames(table_smallintestine)

# Give names to the columns of the new table
colnames(my_results) = c("Non-zero mean (I)", "Mean (I)", "Nº of 0 (I)", "Nº of non-0 (I)") 

results_normal <- as.data.frame(my_results)
  results_normal[is.na(results_normal)] <- 0
  
results_intermediate <- as.data.frame(my_results)
  results_intermediate[is.na(results_intermediate)] <- 0

results_smallintestine <- as.data.frame(my_results)
  results_smallintestine[is.na(results_smallintestine)] <- 0

results_total <- as.data.frame(my_results)
  results_total[is.na(results_total)] <- 0
    
  ##Merge all results in one table
  abundance_normal_intermediate <- merge(results_normal, results_intermediate, by="row.names")
      rownames(abundance_normal_intermediate) <- abundance_normal_intermediate[,1]
      abundance_normal_intermediate <- abundance_normal_intermediate[,-1]
 
  abundance_groups <- merge(abundance_normal_intermediate, results_smallintestine, by="row.names")   
      rownames(abundance_groups) <- abundance_groups[,1]
      abundance_groups <- abundance_groups[,-1]   
      
  
abundance_results <- write.table(abundance_groups, file = "~/abundance_per_groups.txt", quote = F, sep = "\t")
total_abundance_results <- write.table(results_total, file = "~/abundance_total.txt", quote = F, sep = "\t")
