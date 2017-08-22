
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

ns_matrix_list <- list()

for (x in 1:ncol(p_a_table))  {

  name_column <- colnames(p_a_table)[x]
  taxonomy <- subset(p_a_table, select = name_column)
  
  model_table <- merge(taxonomy, normal_small_metadata, by = "row.names" ) 
    row.names(model_table) <- model_table[,1]
    model_table <- model_table[,-1]
  
  colnames(model_table)[1] <- "Taxonomy"  

  model <- glm(Taxonomy ~ . , family = binomial(link = "logit"), data = model_table)
  summary_table <- summary(model)
  coefficients_table <- as.data.frame(summary_table$coefficients)
  
  normal_samples <- subset(model_table, model_table$Group=="normal")
    n_nsum = nsum(normal_samples[,1])
    n_nzsum = nzsum(normal_samples[,1])
  
  small_samples <- subset(model_table, model_table$Group=="small intestine")
    s_nsum = nsum(small_samples[,1])
    s_nzsum = nzsum(small_samples[,1])
  
  loop_matrix <- matrix(ncol = 8, nrow = nrow(coefficients_table))
      colnames(loop_matrix) <- c("Taxonomy","Variable","effect","p_value", "presence_normal", "absence_normal", "presence_small", "absence_small")
  
  a = 1
  
    for (jj in 1:nrow(coefficients_table)) {
      
      if (coefficients_table[jj,4]<0.05) {
        
        loop_matrix[a,1] = name_column
        loop_matrix[a,2] = rownames(coefficients_table)[jj]
        loop_matrix[a,3] = coefficients_table[jj,1]
        loop_matrix[a,4] = coefficients_table[jj,4]
        loop_matrix[a,5] = n_nsum
        loop_matrix[a,6] = n_nzsum
        loop_matrix[a,7] = s_nsum
        loop_matrix[a,8] = s_nzsum
      a <- a + 1  
      }
    }
  loop_matrix <- na.omit(loop_matrix)
  ns_matrix_list[[x]] <- loop_matrix

}

all_matrix_normal_small <- as.data.frame(Reduce(MyMerge, ns_matrix_list))


# NORMAL - INTERMEDIATE #

normal_intermediate_metadata <- metadata_table[!(metadata_table$Group=="small intestine"),]

ni_matrix_list <- list()

for (x in 1:ncol(p_a_table))  {
  
  name_column <- colnames(p_a_table)[x]
  taxonomy <- subset(p_a_table, select = name_column)
  
  model_table <- merge(taxonomy, normal_intermediate_metadata, by = "row.names" ) 
  row.names(model_table) <- model_table[,1]
  model_table <- model_table[,-1]
  
  colnames(model_table)[1] <- "Taxonomy"  
  
  model <- glm(Taxonomy ~ . , family = binomial(link = "logit"), data = model_table)
  summary_table <- summary(model)
  coefficients_table <- as.data.frame(summary_table$coefficients)
  
  normal_samples <- subset(model_table, model_table$Group=="normal")
  n_nsum = nsum(normal_samples[,1])
  n_nzsum = nzsum(normal_samples[,1])
  
  intermediate_samples <- subset(model_table, model_table$Group=="intermediate")
  i_nsum = nsum(intermediate_samples[,1])
  i_nzsum = nzsum(intermediate_samples[,1])
  
  loop_matrix <- matrix(ncol = 8, nrow = nrow(coefficients_table))
  colnames(loop_matrix) <- c("Taxonomy","Variable","effect","p_value", "presence_normal", "absence_normal", "presence_intermediate", "absence_intermediate")
  
  a = 1
  
  for (jj in 1:nrow(coefficients_table)) {
    
    if (coefficients_table[jj,4]<0.05) {
      
      loop_matrix[a,1] = name_column
      loop_matrix[a,2] = rownames(coefficients_table)[jj]
      loop_matrix[a,3] = coefficients_table[jj,1]
      loop_matrix[a,4] = coefficients_table[jj,4]
      loop_matrix[a,5] = n_nsum
      loop_matrix[a,6] = n_nzsum
      loop_matrix[a,7] = i_nsum
      loop_matrix[a,8] = i_nzsum
      a <- a + 1  
    }
  }
  loop_matrix <- na.omit(loop_matrix)
  ni_matrix_list[[x]] <- loop_matrix
  
}

all_matrix_normal_intermediate <- as.data.frame(Reduce(MyMerge, ni_matrix_list))



# SMALL - INTERMEDIATE #

small_intermediate_metadata <- metadata_table[!(metadata_table$Group=="normal"),]

si_matrix_list <- list()

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
  
  small_samples <- subset(model_table, model_table$Group=="small")
  s_nsum = nsum(small_samples[,1])
  s_nzsum = nzsum(small_samples[,1])
  
  intermediate_samples <- subset(model_table, model_table$Group=="intermediate")
  i_nsum = nsum(intermediate_samples[,1])
  i_nzsum = nzsum(intermediate_samples[,1])
  
  loop_matrix <- matrix(ncol = 8, nrow = nrow(coefficients_table))
  colnames(loop_matrix) <- c("Taxonomy","Variable","effect","p_value", "presence_small", "absence_small", "presence_intermediate", "absence_intermediate")
  
  a = 1
  
  for (jj in 1:nrow(coefficients_table)) {
    
    if (coefficients_table[jj,4]<0.05) {
      
      loop_matrix[a,1] = name_column
      loop_matrix[a,2] = rownames(coefficients_table)[jj]
      loop_matrix[a,3] = coefficients_table[jj,1]
      loop_matrix[a,4] = coefficients_table[jj,4]
      loop_matrix[a,5] = s_nsum
      loop_matrix[a,6] = s_nzsum
      loop_matrix[a,7] = i_nsum
      loop_matrix[a,8] = i_nzsum
      a <- a + 1  
    }
  }
  loop_matrix <- na.omit(loop_matrix)
  si_matrix_list[[x]] <- loop_matrix
  
}

all_matrix_small_intermediate <- as.data.frame(Reduce(MyMerge, si_matrix_list))










  
  
   