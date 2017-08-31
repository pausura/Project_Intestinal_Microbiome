### Logistic regresion - Adonis data ###

metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/logistic_regression_metadata_adonis.txt", sep = "\t", header = T, row.names = 1, check.names = F)
metadata_table <- metadata

##Remove NA values
# Convert categorical values to numeric
for (i in 1:ncol(metadata_table)) {
  if (is.factor(metadata_table[,i]) & any(is.na(metadata_table[,i]))) {
    metadata_table[,i] <- as.integer(metadata_table[,i])
  }
}

# Replace NA values: median value
library (psych) 

for (ii in 1:ncol(metadata_table)){
  for (jj in 1:nrow(metadata_table)) {
    
    if (is.na(metadata_table[jj,ii])){
      
      x = describe(metadata_table[,ii])
      a = x$median
      metadata_table[jj,ii] = a
      
    }
  }
}

##Remove columns with only one level
metadata_table <- metadata_table[, sapply(metadata_table, function(col) length(unique(col))) > 1]

taxonomy_table <- read.table("~/Documents/Universitat/Holanda/Projecte/abundance_filtered_taxonomy.txt", sep = "\t", header = T, row.names = 1, check.names = F)

##Create presence/absence table
p_a_table <- taxonomy_table

for (i in 1:ncol(taxonomy_table)) {
  
  for (j in 1:nrow(taxonomy_table)) {
    
    if (taxonomy_table[j,i]>0) {
      
      p_a_table[j,i] = 1
      
    } 
  }
}

##Function to calculate nº of 0
nzsum <- function(x){
  sum (x==0)
}
##Function to calculate nº of non-0
nsum <- function(x){
  sum (x!=0)
}


# NORMAL - SMALL INTESTINE #
normal_small_metadata <- metadata_table[!(metadata_table$Group=="intermediate"),]

# NORMAL - INTERMEDIATE #
normal_intermediate_metadata <- metadata_table[!(metadata_table$Group=="small intestine"),]

# INTERMEDIATE - SMALL INTESTINE #
intermediate_small_metadata <- metadata_table[!(metadata_table$Group=="normal"),]

#Create a list of matrix to save the different results
matrix_list <- list()
table_variables <- matrix(ncol = 2, nrow = ncol(p_a_table))

# For each taxonomy
for (x in 1:ncol(p_a_table))  {
  
  #Get column
  name_column <- colnames(p_a_table)[x]
  taxonomy <- subset(p_a_table, select = name_column)
  
  #Create a table for the model. Merge taxonomy column with metadata  
  model_table <- merge(taxonomy, normal_small_metadata, by = "row.names" ) 
  row.names(model_table) <- model_table[,1]
  model_table <- model_table[,-1]
  #Change taxonomy name
  colnames(model_table)[1] <- "Taxonomy"  
  
  #Sort
  model_table <- model_table[ , order(names(model_table))]
  
  #Calculate model
  model <- glm(Taxonomy ~ . , family = binomial(link = "logit"), data = model_table)
  
  ##Calculate Anova
  anova_test <- anova(model, test = "Chisq")
  
  ##Keep significative variables for model2
  variables_model2 <- subset(anova_test, anova_test[,5]<0.05)
  
  list_variables_model1 <- as.vector(c(colnames(model_table)))
  list_variables_model1 <- paste(c(list_variables_model1), collapse=',' )
  
  
  #If there are significative variables
  if (nrow(variables_model2)>0) {
    
    # Save names of the significative variables: create a new metadata table with these variables (model_table2)
    matrix_variables_model2 <- as.data.frame(rownames(variables_model2))
    rownames(matrix_variables_model2) <- matrix_variables_model2[,1]
    t_model_table <- t(model_table)
    
    list_variables_model2 <- as.vector(matrix_variables_model2$`rownames(variables_model2)`)
    list_variables_model2 <- paste(c(list_variables_model2), collapse=', ' )
    
    t_model_table2 <- merge(matrix_variables_model2, t_model_table, by = "row.names")
    rownames(t_model_table2) <- t_model_table2[,1]
    t_model_table2[1:2] <- NULL
    
    model_table2 <- t(t_model_table2) 
    # Merge the new table with the taxonomy again (lost in the last merge)
    model_table2 <- merge(taxonomy, model_table2, by = "row.names" ) 
    row.names(model_table2) <- model_table2[,1]
    model_table2 <- model_table2[,-1]
    
    #Change taxonomy name
    colnames(model_table2)[1] <- "Taxonomy"
    
    # Change character to numeric
    for (ii in 1:ncol(model_table2))  {
      column_name2 <- colnames(model_table2)[ii]
      
      for (jjj in 1:ncol(model_table)) {
        column_name <- colnames(model_table)[jjj]
        
        if (column_name == column_name2 & is.numeric(model_table[,jjj])) {
          model_table2[,ii] <- as.numeric(as.character(model_table2[,ii]))
        }
      }
    }
    
    # Sort
    model_table2 <- model_table2[ , order(names(model_table2))]
    
    ## Calculate model2 with the new table: contains only the significative variables in anova test
    model2 <- glm(Taxonomy ~ . , family = binomial(link="logit"), data = model_table2)
    
    # Test the two models
    anova_test <- anova(model, model2, test = "Chisq")
    
    # Model 2 and model 1 are equal: save model 1 results
    if (is.na(anova_test[2,5])) {
      
      summary_table <- summary(model)  
      
      coefficients_table <- as.data.frame(summary_table$coefficients)
      
      category1_samples <- subset(model_table, model_table$Group=="normal")
      nsum1 = nsum(category1_samples$Taxonomy)
      nzsum1 = nzsum(category1_samples$Taxonomy)
      
      category2_samples <- subset(model_table, model_table$Group=="small intestine")
      nsum2 = nsum(category2_samples$Taxonomy)
      nzsum2 = nzsum(category2_samples$Taxonomy)
      
      loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
      colnames(loop_matrix) <- c("Taxonomy","presence_normal", "absence_normal", "presence_small", "absence_small", "Variable","effect","p_value", "Model")
      
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
        loop_matrix[a,9] = "model_1"
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
      
      table_variables[x,1] = name_column
      table_variables[x,2] = list_variables_model1
      
    }
    
    #Model 2 is not better than model 1: save model 1 results
    else if (anova_test[2,5]<0.05) {
      
      summary_table <- summary(model)  
      
      coefficients_table <- as.data.frame(summary_table$coefficients)
      
      category1_samples <- subset(model_table, model_table$Group=="normal")
      nsum1 = nsum(category1_samples$Taxonomy)
      nzsum1 = nzsum(category1_samples$Taxonomy)
      
      category2_samples <- subset(model_table, model_table$Group=="small intestine")
      nsum2 = nsum(category2_samples$Taxonomy)
      nzsum2 = nzsum(category2_samples$Taxonomy)
      
      loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
      colnames(loop_matrix) <- c("Taxonomy","presence_normal", "absence_normal", "presence_small", "absence_small", "Variable","effect","p_value", "Model")
      
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
        loop_matrix[a,9] = "model_1"
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
      
      table_variables[x,1] = name_column
      table_variables[x,2] = list_variables_model1
      
    }
    
    ##Model 2 is better than model 1: save model 2 results
    else {
      
      summary_table2 <- summary(model2)  
      
      #Save coefficients results
      coefficients_table <- as.data.frame(summary_table2$coefficients)
      
      #Get number presence/absence of the taxonomy for each category
      category1_samples <- subset(model_table2, model_table2$Group=="normal")
      #Presence
      nsum1 = nsum(category1_samples$Taxonomy)
      #Absence
      nzsum1 = nzsum(category1_samples$Taxonomy)
      
      category2_samples <- subset(model_table2, model_table2$Group=="small intestine")
      #Presence      
      nsum2 = nsum(category2_samples$Taxonomy)
      #Absence
      nzsum2 = nzsum(category2_samples$Taxonomy)
      
      loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
      colnames(loop_matrix) <- c("Taxonomy","presence_normal", "absence_normal", "presence_small", "absence_small", "Variable","effect","p_value", "Model")
      
      a = 1
      #Save in a new matrix:
      for (jj in 1:nrow(coefficients_table)) {
        
        loop_matrix[a,1] = name_column    #name taxonomy
        loop_matrix[a,2] = nsum1          #presence taxonomy in category1
        loop_matrix[a,3] = nzsum1         #absence taxonomy in category1
        loop_matrix[a,4] = nsum2          #presence taxonomy in category2
        loop_matrix[a,5] = nzsum2         #absence taxonomy in category2
        loop_matrix[a,6] = rownames(coefficients_table)[jj]   #name of the variable
        loop_matrix[a,7] = coefficients_table[jj,1]    #effect value
        loop_matrix[a,8] = coefficients_table[jj,4]    #p_value
        loop_matrix[a,9] = "model_2"
        a <- a + 1  
      }
      
      loop_matrix <- na.omit(loop_matrix)  #remove empty rows (NA values)
      loop_matrix1 <- loop_matrix
      
      for (kk in 1:nrow(loop_matrix)) {   #remove (Intercept) results
        
        if (loop_matrix[kk,6]=="(Intercept)") {
          
          loop_matrix1 <- loop_matrix[-kk,]
        }
      }
      
      matrix_list[[x]] <- loop_matrix1   #Save the new matrix in a list of matrix
      
      table_variables[x,1] = name_column
      table_variables[x,2] = list_variables_model2
    }
    
  }
  
  ## If are not significative variables in anova test: keep model1 results
  else {
    
    summary_table <- summary(model)  
    
    coefficients_table <- as.data.frame(summary_table$coefficients)
    
    category1_samples <- subset(model_table, model_table$Group=="normal")
    nsum1 = nsum(category1_samples$Taxonomy)
    nzsum1 = nzsum(category1_samples$Taxonomy)
    
    category2_samples <- subset(model_table, model_table$Group=="small intestine")
    nsum2 = nsum(category2_samples$Taxonomy)
    nzsum2 = nzsum(category2_samples$Taxonomy)
    
    loop_matrix <- matrix(ncol = 9, nrow = nrow(coefficients_table))
    colnames(loop_matrix) <- c("Taxonomy","presence_normal", "absence_normal", "presence_small", "absence_small", "Variable","effect","p_value", "Model")
    
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
      loop_matrix[a,9] = "model_1"
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
    
    table_variables[x,1] = name_column
    table_variables[x,2] = list_variables_model1
  }
}


#Save in the same matrix all the results
all_matrix <- as.data.frame(do.call(rbind, matrix_list))

colnames(table_variables) <- c("Taxonomy", "Variables")
write.table(table_variables, file = "~/table_variables_logistic_regression.txt", sep = "\t", quote = F)

#Split by Variable   
split_matrix <- split(all_matrix, all_matrix$Variable)

for (bb in 1:length(split_matrix)){
  
  #Correct by p_values
  matrix <- as.data.frame(split_matrix[bb])
  p_value <- as.vector(matrix[,8])
  
  corrected_pvalues <- p.adjust(p_value, method = "fdr")
  
  #Add a new column with the new p_values 
  matrix <- cbind(matrix, corrected_pvalues)
  
  name <- colnames(matrix)[6]
  nc <- paste(name, ".txt", sep ="")
  assign(nc, matrix)
  final_name_matrix <- paste("~/", nc, sep = "")
  
  write.table(matrix, file = final_name_matrix, quote = F, sep = "\t")
  
  ## Filtering significant results
  
  #Filter by the new p_values
  filtered_matrix <- subset(matrix, matrix[,10]<0.05)
  
  name2 <- colnames(filtered_matrix)[6]
  nc2 <- paste(name2, "_filtered.txt", sep ="")
  
  assign(nc2, filtered_matrix)
  final_name_matrix2 <- paste("~/", nc2, sep = "")
  
  #Not print empty tables
  if (nrow(filtered_matrix)>0) {
    
    write.table(filtered_matrix, file = final_name_matrix2, quote = F, sep = "\t" ) 
  }
  
}
