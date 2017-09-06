
##########
#Metadata#
##########
metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/logistic_metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
group <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1, check.names = F)

metadata <- merge(metadata, group, by = "row.names")
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]

metadata_table <- metadata

 ##Remove NA values
        #Numeric: use median value
library (psych) 

for (i in 1:ncol(metadata_table)){
  
  for (j in 1:nrow(metadata_table)) {
    
    if (is.numeric(metadata_table[,i])){
      if (is.na(metadata_table[j,i])){
       
         x = describe(metadata_table[,i])
         a = x$median
         metadata_table[j,i] = a
        
} } } }

      #Categoric: remove row
metadata_table <- na.omit(metadata_table)


##########
#Taxonomy#
##########
taxonomy_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_taxonomy.txt", sep = "\t", header = T, row.names = 1, check.names = F)
p_a_table <- taxonomy_table

for (i in 1:ncol(taxonomy_table)) {
  
  for (j in 1:nrow(taxonomy_table)) {
    
    if (taxonomy_table[j,i]>0) {
      
      p_a_table[j,i] = 1
      
    } } }


##Function to calculate nº of 0
nzsum <- function(x){
  sum (x==0)
}
##Function to calculate nº of non-0
nsum <- function(x){
  sum (x!=0)
}

##Merge matrix list
MyMerge <- function(x, y){
  df <- rbind(x, y[, colnames(x)])
  return(df)
}


# NORMAL - SMALL INTESTINE #
normal_small_metadata <- metadata_table[!(metadata_table$Group=="intermediate"),]

# NORMAL - INTERMEDIATE #
normal_intermediate_metadata <- metadata_table[!(metadata_table$Group=="small intestine"),]

# SMALL - INTERMEDIATE #
small_intermediate_metadata <- metadata_table[!(metadata_table$Group=="normal"),]


matrix_list <- list()

for (x in 1:ncol(p_a_table))  {

  name_column <- colnames(p_a_table)[x]
  taxonomy <- subset(p_a_table, select = name_column)
  
  model_table <- merge(taxonomy, small_intermediate_metadata, by = "row.names" ) 
    row.names(model_table) <- model_table[,1]
    model_table <- model_table[,-1]
  
  colnames(model_table)[1] <- "Taxonomy"  

  model <- glm(Taxonomy ~ . , family = binomial(link = "logit"), data = model_table)
  summary_table <- summary(model)
  coefficients_table <- as.data.frame(summary_table$coefficients)
  
  normal_samples <- subset(model_table, model_table$Group=="intermediate")
    nsum1 = nsum(normal_samples[,1])
    nzsum1 = nzsum(normal_samples[,1])
  
  small_samples <- subset(model_table, model_table$Group=="small intestine")
    nsum2 = nsum(small_samples[,1])
    nzsum2 = nzsum(small_samples[,1])
  
  loop_matrix <- matrix(ncol = 8, nrow = nrow(coefficients_table))
      colnames(loop_matrix) <- c("Taxonomy","presence_intermediate", "absence_intermediate", "presence_small", "absence_small", "Variable","effect","p_value")
  
  a = 1
  
    for (jj in 1:nrow(coefficients_table)) {
      
        loop_matrix[a,1] = name_column
        loop_matrix[a,2] = nsum1
        loop_matrix[a,3] = nzsum1
        loop_matrix[a,4] = nsum2
        loop_matrix[a,5] = nzsum2
        loop_matrix[a,6] = rownames(coefficients_table)[jj]
        loop_matrix[a,7] = coefficients_table[jj,1]
        loop_matrix[a,8] = coefficients_table[jj,4]
      a <- a + 1  
    }
  
  loop_matrix <- na.omit(loop_matrix)
  loop_matrix1 <- loop_matrix
  
  for (kk in 1:nrow(loop_matrix)) {
    
    if (loop_matrix[kk,6]=="(Intercept)") {
      
      loop_matrix1 <- loop_matrix[-kk,]
    }
  }

  matrix_list[[x]] <- loop_matrix1

}

all_matrix <- as.data.frame(Reduce(MyMerge, matrix_list))
split_matrix <- split(all_matrix, all_matrix$Variable)

for (bb in 1:length(split_matrix)){
  
  matrix <- as.data.frame(split_matrix[bb])
  p_value <- as.vector(matrix[,8])
  
  corrected_pvalues <- p.adjust(p_value, method = "fdr")
  
  matrix <- cbind(matrix, corrected_pvalues)

      name <- colnames(matrix)[6]
      nc <- paste(name, ".txt", sep ="")
      assign(nc, matrix)
      final_name_matrix <- paste("~/", nc, sep = "")
  
  write.table(matrix, file = final_name_matrix, quote = F, sep = "\t")
 
  ## Filtering significant results
   
  filtered_matrix <- subset(matrix, matrix[,9]<0.05)

      name2 <- colnames(filtered_matrix)[6]
      nc2 <- paste(name2, "_filtered.txt", sep ="")
      
      assign(nc2, filtered_matrix)
      final_name_matrix2 <- paste("~/", nc2, sep = "")
  
  
 write.table(filtered_matrix, file = final_name_matrix2, quote = F, sep = "\t" ) 

}



## View_overlap

tax_small_intestine <- as.data.frame(Groupsmall.intestine.Variable_filtered.txt[,1])
tax_age <- as.data.frame(AgeAtFecalSampling.Variable_filtered.txt[,1])
tax_PPI <- as.data.frame(PPIyes.Variable_filtered.txt[,1])

common_small_age <- as.data.frame(intersect(tax_small_intestine[,1], tax_age[,1]))
common_small_PPI <- as.data.frame(intersect(tax_small_intestine[,1], tax_PPI[,1]))
common_age_PPI <- as.data.frame(intersect(tax_age[,1], tax_PPI[,1])) # No data
common_total <- as.data.frame(intersect(common_small_age[,1], tax_PPI[,1])) # No data

plot_table <- read.table("~/Documents/Universitat/Holanda/Projecte/p_a_plot_table_intermediate_small.txt", sep = "\t", header = T, row.names = 1, check.names = F)
group_plot_table <- merge(plot_table, group, by = "row.names")
  rownames(group_plot_table) <- group_plot_table[,1]
  group_plot_table <- group_plot_table[,-1]

plot_matrix = matrix(ncol = 3, nrow = 15)

normal_plot_table <- subset(group_plot_table, group_plot_table$Group=="normal")
  normal_plot_table$Group <- NULL
  
small_plot_table <- subset(group_plot_table, group_plot_table$Group=="small intestine")
  small_plot_table$Group <- NULL
  
intermediate_plot_table <- subset(group_plot_table, group_plot_table$Group=="intermediate")
  intermediate_plot_table$Group <- NULL

group_plot_table$Group <- NULL  
  
  for (k in 1:ncol(group_plot_table)) {
    
    percentage_normal = (nsum(intermediate_plot_table[,k])/nrow(intermediate_plot_table))
    percentage_small = (nsum(small_plot_table[,k])/nrow(small_plot_table))
    
    name_column <- colnames(group_plot_table)[k]
    name_column <- unlist(strsplit(name_column, split='|', fixed=TRUE))
    new_name_column <- tail(name_column, n = 1)
    
    plot_matrix[k,1] = new_name_column
    plot_matrix[k,2] = percentage_normal
    plot_matrix[k,3] = percentage_small
   
    
  }

  colnames(plot_matrix) <- c("Taxonomy", "%_intermediate", "%_small")
  
  library(ggplot2)
  library(reshape2)
  plot_matrix=as.data.frame(plot_matrix)
  test_plot1 = melt(plot_matrix, id=1)
  
  ggplot(test_plot1, aes(x=Taxonomy, y=as.numeric(value), group=variable, fill=variable)) + geom_bar(stat="identity",position="dodge") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
  
  
  

  
  
   