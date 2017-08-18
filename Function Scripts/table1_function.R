### Calculate table1 ###

##Open data
total_metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1, check.names = F)


calculate_table1(total_metadata, intestinal_groups)

calculate_table1 <- function (metadata_input, category_table) {
  
  # Create other functions to calculate the different parameters
  
      ## Categorical values - create function to calculate the counts and the percentage for categorical variables
         tblFun <- function(x) {
           # Create a table
             tbl <- table(x)
           # Combine columnes/rows to get the counts and percentage (creates new table -> res)
             res <- cbind(tbl,round(prop.table(tbl)*100,2))
           # Give names to the columns
             colnames(res) <- c('Count','Percentage')
             res
         }
         
      ## NA sum function - counts the number of NA
         nzsum <- function(x) {
            sum (is.na(x))
         }
         
  # Packages needed       
    library (psych)  #describe r function

  # Merge metadata_table with category_table
    function_metadata <- merge(category_table, metadata_input, by ="row.names")  
    # First column as rownames
      rownames(function_metadata) <- function_metadata[,1]
      function_metadata <- function_metadata[,-1]
  # Create a new column to assign a number to each category
    function_metadata$category <- as.integer(function_metadata[,1])      
    category_number <- nlevels(function_metadata[,1])
    
  # Save the different categories files in the global environment  
    matrix_list <- list()
    j = 1
  
  # Loop for to subset the different variables by category
    for (i in 1:category_number) {
     
      category_matrix <- subset(function_metadata, function_metadata$category== i )
      nam <- paste("Category", j , sep = "")
      matrix_list[[j]] <-  assign(nam, category_matrix)
      
      j <- j + 1  
    }
    
  # Save in different matrix the variables for each category 
    for (ii in 1:category_number) {
      
      new_matrix <- as.data.frame(matrix_list[[ii]])
      my_results = matrix(ncol = 6, nrow = ncol(new_matrix))
      
      # For each loop goes to the next column (numerical way)
        for (iii in 1:ncol(new_matrix)) { 
        
        # Condition: if the column values are numerical/continuous
          if (is.numeric(new_matrix[,iii])) {
          # Keep in "x" the result from describe function (done in the columns) - for each factor  
            x = describe(new_matrix[,iii])
            z = nzsum(new_matrix[,iii])
          # In the new table ("x"): keep different values in the different columns
            my_results[iii,1] = "numerical"
            my_results[iii,2] = x$median
            my_results[iii,3] = x$mean
            my_results[iii,4] = x$sd
            my_results[iii,5] = x$n
            my_results[iii,6] = z
          }
        # Condition: if the column values are categorical  
          else {
          # Keep in "x" the result from tblFun function (done in the columns) - for each factor
            x = tblFun(new_matrix[,iii])
            z = nzsum(new_matrix[,iii])
          # In the new table ("x"): keep different values in the different columns 
            my_results[iii,1]="categorical"
          # toString to keep the possible different values/categories in the same vector/column
            my_results[iii,2]=toString(rownames(x))
          # First column table x = 'Count'
            my_results[iii,3]=toString(x[,1]) 
          # Second column table x = 'Percentage'
            my_results[iii,4]=toString(x[,2])
          # Sum of the values on column1 ("x")
            my_results[iii,5]=sum(x[,1])
            my_results[iii,6]= z
          }
        }
      
      # The column names from the original table = row names from the new table 
        rownames(my_results) = colnames(new_matrix)
      # Give names to the columns of the new table 
        colnames(my_results) = c("Type", "Categories/Median", "Counts/Mean", "SD/%", "Number_non-zeros(n)", "Number_NA") 
      
      # Save the name of the variable to title the data.frame (table)
        name_category <- new_matrix[1,1] 
        name_matrix <- paste(name_category, "_metadata_table1.txt", sep = "")
        final_name_matrix <- paste("~/", name_matrix, sep = "")
        
      # Export the new table
        write.table (my_results, file = final_name_matrix , quote = F, sep = "\t") 
        
    }
    
    
    ## Calculate table1 with the whole data:
    
    my_results = matrix(ncol = 6, nrow = ncol(function_metadata))
    
    for (k in 1:ncol(function_metadata)){
      
      if (is.numeric(function_metadata[,k])) {
        # Keep in "x" the result from describe function (done in the columns) - for each factor  
          x = describe(function_metadata[,k])
          z = nzsum(function_metadata[,k])
        # In the new table ("x"): keep different values in the different columns
          my_results[k,1] = "numerical"
          my_results[k,2] = x$median
          my_results[k,3] = x$mean
          my_results[k,4] = x$sd
          my_results[k,5] = x$n
          my_results[k,6] = z
      }
      # Condition: if the column values are categorical  
      else {
        # Keep in "x" the result from tblFun function (done in the columns) - for each factor
          x = tblFun(function_metadata[,k])
          z = nzsum(function_metadata[,k])
        # In the new table ("x"): keep different values in the different columns 
          my_results[k,1]="categorical"
        # toString to keep the possible different values/categories in the same vector/column
          my_results[k,2]=toString(rownames(x))
        # First column table x = 'Count'
          my_results[k,3]=toString(x[,1]) 
        # Second column table x = 'Percentage'
          my_results[k,4]=toString(x[,2])
        # Sum of the values on column1 ("x")
          my_results[k,5]=sum(x[,1])
          my_results[k,6]= z
      }
    }
    
    # The column names from the original table = row names from the new table 
      rownames(my_results) = colnames(function_metadata)
    # Give names to the columns of the new table 
      colnames(my_results) = c("Type", "Categories/Median", "Counts/Mean", "SD/%", "Number_non-zeros(n)", "Number_NA") 
    
    # Export the new table
      write.table (my_results, file = "~/total_metadata_table1.txt" , quote = F, sep = "\t")  
    
}

    

  
  
  





