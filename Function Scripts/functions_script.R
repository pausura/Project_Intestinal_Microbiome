#### FUNCTIONS FILE ####

######################
#Filter data function#
######################
filter_data <- function(tax_input, metadata_input, parameter_file, parameter_threshold, metadata_threshold) {
  
  ##Filter the parameter_file by the given threshold
  parameter_threshold_table <- as.data.frame(parameter_file[with(parameter_file, parameter_file[,1]>parameter_threshold),])
  #Write table
  write.table(parameter_threshold_table, file = "~/parameter_threshold_table.txt", quote = F, sep ="\t")
  
  #Merge the parameter_threshold_table with tax_input file  
  filtered_tax <- merge(parameter_threshold_table, tax_input, by="row.names")
  #First column as Rownames (sample identification)
  rownames(filtered_tax) <- filtered_tax[,1]
  filtered_tax <- filtered_tax[,-1]
  #Remove parameter_file column and sample_names (repeated)  
  filtered_tax[1:2] <- NULL
  
  #Write table  
  write.table(filtered_tax, file = "~/filtered_tax.txt", quote= F, sep="\t")
  
  ##Filter the metadata_file by the given parameter_threshold
  filtered_metadata <- merge(parameter_threshold_table, metadata_input, by="row.names")
  #First column as Rownames (sample identification)
  rownames(filtered_metadata) <- filtered_metadata[,1]
  filtered_metadata <- filtered_metadata[,-1]
  #Remove parameter_file column and sample_names (repeated)  
  filtered_metadata[1:2] <- NULL
  
  ##Filter the filtered_metadata by the given metadata_threshold
  filtered_total_metadata <- filtered_metadata
  
  ##Remove factors/variables with >metadata_threshold (number max of tolerated NA)
  cond_na <- sapply(filtered_total_metadata, function(col) sum(is.na(col)) < metadata_threshold)
  filtered_total_metadata <- filtered_total_metadata[, cond_na, drop = FALSE]
  
  ##Remove columns with only one level
  filtered_total_metadata <- filtered_total_metadata[, sapply(filtered_total_metadata, function(col) length(unique(col))) > 1]
  
  #Write table
  write.table(filtered_total_metadata, file = "~/filtered_metadata.txt", quote = F, sep="\t")
  
}

#########################
#Taxonomy level function#
#########################
level_function <- function(taxonomy_table,level_taxonomy) {
  
  #Transpose taxonomy table
  t_tax_table <- as.data.frame(t(taxonomy_table))
  #Create a table to save the species taxonomy
  loop_table <- as.data.frame(matrix(nrow = nrow(t_tax_table) , ncol = ncol(t_tax_table)))
  
  
  ## Count only the desired levels
  for (i in 1:nrow(t_tax_table)) {
    
    if (count.fields(textConnection(row.names(t_tax_table[i,])), sep="|") == level_taxonomy){
      #print (paste0("Species found: ", row.names(tax_table[i,]))) ##Loop check
      
      loop_table[i,] = t_tax_table[i,]
      
    }
    
  }
  
  # Give row names to the new table
  row.names(loop_table) = row.names(t_tax_table)
  # Give column names to  new table
  colnames(loop_table) = colnames(t_tax_table)
  
  ##Remove all rows with NA values  
  level_table <- na.omit(loop_table)
  
  ##Transpose level table to get the taxonomy in columns and samples in rows
  t_level_table <- as.data.frame(t(level_table))
  
  for (i in 1:ncol(t_level_table)) {
    name_column <- colnames(t_level_table)[i]
    new_name_column = unlist(strsplit(name_column, split='|', fixed=TRUE))[level_taxonomy]
    colnames(t_level_table)[i] <- new_name_column
    
  }
  
  #Export table
  write.table(t_level_table, file = "~/level_table.txt", quote = F, sep = "\t")
  
}

########################
#Shannon Index function#
########################
shannon_function <- function(level_table, group_table) {
  
  ##Required packages - shannon diversity/ggplot
  library(vegan)
  library(ggplot2)
  
  ## Calculate shannon (alpha) for each sample
  alpha <- as.data.frame(diversity(level_table, index="shannon"))
  colnames(alpha)[1] <- "alpha_diversity"
  ## Divide alpha results by intestinal groups
  #Merge alpha with intestinal groups table
  group_taxa <- merge(group_table, alpha, by="row.names")
  rownames(group_taxa) <- group_taxa[,1]
  group_taxa <- group_taxa[,-1]
  #colnames(group_taxa)[2] <- "diversity"
  # Calculate the number of categories
  category_number <- nlevels(group_taxa[,1])
  # Create a new column to assign a number to each category
  group_taxa$category <- as.integer(group_taxa[,1])    
  # Create a new column to colour the plot depending on the category (level)        
  group_taxa$color = "none" 
  # Create a palette of colors depending on the number of categories
  my_palette <- matrix(brewer.pal(category_number,"Set1"))
  # For loop to assign one different color to each different category in a new column
  for (i in 1:category_number){
    group_taxa[group_taxa$category == i,]$color = my_palette[i,1]
  }
  
  ### Create a violin plot
  shannon_plot <- ggplot(group_taxa, aes(x=category, y=alpha_diversity, fill = group_taxa$color)) + labs (y="Shannon Diversity Index", x="Category") + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1) + theme_classic() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 1, size=16,color="black")) + scale_color_identity("All_categories", breaks = group_taxa$color, labels= group_taxa$category, guide = "legend")
  
  ##Save the plots as pdf.file
  pdf("shannon_plot.pdf")
  print(shannon_plot)
  
  dev.off()
  
  
}

###############################
#Taxonomy_composition function#
###############################
tax_composition_differences <- function(tax_level_table, category_table, top_tax_value) {
  
  # Merge tax_level_table with category_table file  
  filum_groups <- merge(category_table, tax_level_table, by = "row.names")
  # First column as Rownames (sample identification)  
  rownames(filum_groups) <- filum_groups[,1]
  filum_groups <- filum_groups[,-1]
  
  # Split filum_groups table by categories      
  split_categories <- split(filum_groups, filum_groups[,1])
  
  # Packages needed
  library(psych)
  library(reshape2) 
  library(ggplot2) 
  
  # Save the different categories files in the global environment
  j = 1
  others_position = top_tax_value + 1
  results_table <- matrix(nrow = others_position, ncol= 2 )
  
  matrix_list <- list()   
  
  # Loop for to calculate the filum mean by category
  for (i in split_categories) {
    # Calculate mean
    summary_table <- describe(i)
    # Find the rows with the highest values by the given value
    top_tax <- summary_table[order(summary_table$mean, decreasing = T)[1:top_tax_value],]
    # Sum values and create another category: others
    sum_top_tax <- sum(top_tax$mean)
    others_tax <- as.numeric((100-sum_top_tax))
    
    total_mean <- cbind(top_tax$mean)
    total_mean <- rbind(total_mean, others_tax)
    
    results_table[,1] <- total_mean
    results_table[,2] <- c(rownames(top_tax), "others")
    
    nam <- paste("Category", j , sep = "")
    matrix_list[[j]] <-  assign(nam, results_table)
    
    j <- j + 1        
  }
  
  # Change colnames from all matrix in the list
  for (k in seq_along(matrix_list)) {
    colnames(matrix_list[[k]]) <- c("mean","bacteria")
  }
  
  # Merge function 
  MyMerge <- function(x, y){
    df <- merge(x, y, by="bacteria" , all.x= TRUE, all.y= TRUE)
    return(df)
  }
  
  # Merge the matrix saved in the matrix_list in a same table
  composition_table <- Reduce(MyMerge, matrix_list)
  # First column as Rownames (sample identification)  
  rownames(composition_table) <- composition_table[,1]
  composition_table <- composition_table[,-1]
  
  # For loop to name the columns depending on the groyup/category in the composition_table  
  for (bb in 1:ncol(composition_table)) {
    nc <- paste("Category", bb , sep = "")
    colnames(composition_table)[bb] = nc
    assign(nc,composition_table)
  }
  
  # Export the table with the top4 abundances per each category          
  write.table(composition_table, file = "~/bacteria_composition.txt", sep = "\t", quote = F)
  
  ## Stacked barplot ##   
  composition_table$bacteria = row.names(composition_table) 
  row.names(composition_table) <- NULL
  my_plot_table <- melt(composition_table, id.vars = "bacteria")
  
  filum_plot <- ggplot (my_plot_table, aes(x=variable, y=as.numeric(value))) + geom_bar (aes(fill = bacteria), stat = "identity") + theme_classic() + xlab("Group") + ylab("relative_abundance")
  
  # Save the plots as pdf.file
  
  pdf("taxonomy_composition_plot.pdf")
  print(filum_plot)
  
  dev.off()
  
}

##################
#Table1 function#
##################
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

########################
#PCoA analysis function#
########################
pcoa_function <- function(tax_level_table, pcoa_elements, variable_table, top_value) {
  
  ##Required packages
  library(vegan)
  library(ggplot2)
  library(psych)
  library(RColorBrewer)
  
  ## Generate a distance matrix - Bray method
  beta <- vegdist(tax_level_table, method="bray")
  ## cmdscale -> multidimensional scaling of a data matrix (distance matrix) 
  my_pcoa <- as.data.frame(cmdscale(beta, k = pcoa_elements))
  
  # Variable that contains the ".pdf" string to save plots in pdf files 
  a <- ".pdf"  
  
  ## If the variable_table is numeric:
  
  if (is.numeric(variable_table[,1])) {
    # Variable that contains colors codes to color the plot 
    my_col=c("#0000FF","#62A4D1","#5BE55B","#FFF000", "#FF0000")
    
    # Merge the variable_table with my_pcoa table --> the tax_level_table and the variable_table need to have the same number of rows/samples
    numeric_new_table <- merge(variable_table, my_pcoa, by="row.names") 
    # First column as Rownames
    rownames(numeric_new_table) <- numeric_new_table[,1]
    numeric_new_table <- numeric_new_table[,-1]
    
    # Create new tables with new row number
    new_variable_table <- subset(numeric_new_table, select = 1:ncol(variable_table))
    new_pcoa <- numeric_new_table
    new_pcoa[1:ncol(variable_table)] <- NULL
    
    # If the category table has more than one column:
    if (ncol(variable_table)>1) {
      # Calculate the mean for each column of the variable_table
      summary_table <- describe(new_variable_table)
      # Find the rows with the highest mean values (more abundant variables): top_value
      top_tax <- summary_table[order(summary_table$mean, decreasing = T)[1:top_value],]
      # Save the names of the highest mean values
      total_tax_names <- cbind(rownames(top_tax))
      # First column as Rownames
      rownames(total_tax_names) <- total_tax_names[,1]
      
      # Transpose the category table to merge  
      t_variable_table <- as.data.frame(t(new_variable_table))    
      # Merge the t_variable_table with the total_tax_names (table with the highest mean values) 
      numeric_plot_table <- merge(total_tax_names, t_variable_table, by = "row.names")    
      # First column as Rownames
      rownames(numeric_plot_table) <- numeric_plot_table[,1]
      numeric_plot_table <- numeric_plot_table[,-1]
      # Remove repeated column (rownames)
      numeric_plot_table <- numeric_plot_table[,-1]
      
      # Transpose the new table to have the variables as columns  
      t_plot_table <- as.data.frame(t(numeric_plot_table))
      
      # For loop to get plots as number of variables in t_plot_table   
      for (jj in 1:ncol(t_plot_table)) {
        # Save the name of the variable to title the plot
        name_category <- colnames(t_plot_table)[jj]
        # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file
        name_pdf <- paste(name_category,a, sep = "")
        # Create the plot
        bla = ggplot(new_pcoa, aes(x=V1, y=V2, geom = "blank", colour = t_plot_table[,jj])) + geom_point() + scale_color_gradientn(colours = my_col, colnames(t_plot_table[,jj])) + theme_classic() + labs(x = "PCoA1", y = "PCoA2") + ggtitle(name_category)
        #x_axis: V1 from my_pcoa table (first PCoA element)
        #y_axis: V2 from my_pcoa table (second PCoA element)
        #colour: a different plot is generated for each column in t_plot_table 
        # Create the pdf file
        pdf(name_pdf)
        # Print the plot
        print(bla)
        # Empty the current device to create the next plot in the next loop
        dev.off()
      }
    }
    
    # The variable_table has only one column:
    else { 
      # Save the name of the variable to title the plot
      name_category <- colnames(new_variable_table)[1]
      # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file   
      name_pdf <- paste(name_category, a, sep = "")
      # Create the plot
      bla = ggplot(new_pcoa, aes(x=V1, y= V2, geom = "blank", colour = new_variable_table[,1])) + geom_point() + scale_color_gradientn(colours = my_col, colnames(variable_table[,1])) + theme_classic() + labs(x="PCoA1", y = "PCoA2") + ggtitle(name_category)
      #x_axis: V1 from my_pcoa table (first PCoA element)
      #y_axis: V2 from my_pcoa table (second PCoA element)
      #colour: colored by the column in the variable_table
      # Create the pdf file
      pdf(name_pdf)
      # Print the plot   
      print(bla)
      # Empty the current device  
      dev.off()
    }
  }
  
  ## The variable_table is categoric:          
  else { 
    
    # Merge the variable_table with my_pcoa table --> needed the same num of rows in both tables
    categoric_plot_table <- merge(variable_table, my_pcoa, by="row.names") 
    # First column as Rownames
    rownames(categoric_plot_table) <- categoric_plot_table[,1]
    categoric_plot_table <- categoric_plot_table[,-1]
    # Calculate the number of categories
    category_number <- nlevels(categoric_plot_table[,1])
    # Create a new column to assign a number to each category
    categoric_plot_table$category <- as.integer(categoric_plot_table[,1])
    # Create a new column to colour the plot depending on the category (level)        
    categoric_plot_table$color = "none" 
    # Create a palette of colors depending on the number of categories
    my_palette <- matrix(brewer.pal(category_number,"Set1"))
    # For loop to assign one different color to each different category in a new column
    for (i in 1:category_number){
      categoric_plot_table[categoric_plot_table$category == i,]$color = my_palette[i,1]
    }
    
    # Save the name of the variable to title the plot
    name_category <- colnames(variable_table)[1]
    # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file  
    name_pdf <- paste(name_category, a, sep = "")
    # Create the plot
    bla = ggplot (categoric_plot_table, aes(x=V1, y=V2, geom="blank", colour=color)) + geom_point () + scale_color_identity("All_categories", breaks=categoric_plot_table$color, labels=categoric_plot_table$category, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
    #x_axis: V1 from my_pcoa table (first PCoA element)
    #y_axis: V2 from my_pcoa table (second PCoA element)
    #colour: colored by the new color column 
    # Create the pdf file
    pdf(name_pdf)
    # Print the plot
    print(bla)
    # Empty the current device
    dev.off()
  }
  
}

################################
#Clustering dendrogram function#
################################
clustering_dendrogram <- function(tax_level_table, category_table, category_number) {
  
  #Packages needed
  library(vegan)
  library(gclus)
  library(data.table)
  
  ## Calculate distance matrix - Bray method
  beta <- vegdist(tax_level_table, method="bray")
  # Calculate AVERAGE distance in beta
  caver <- hclust(beta, method="aver")
  # Reorder the matrix generated    
  caver1 <- reorder.hclust(caver, beta)
  # Get the labels(sample identification) in the same order as are in caver1    
  order_caver1 = as.data.frame(caver1$labels[c(caver1$order)])  
  setDT(order_caver1, keep.rownames = TRUE)[]
  # Give names to the columns and use ID column to put rownames
  colnames(order_caver1)[1] <- "num"
  colnames(order_caver1)[2] <- "ID"
  
  ## Before merge: create a new colum in category_table with the ID
  category_table$ID <- rownames(category_table)    
  
  ## Merge order_caver1 with the category_table        
  order_cat <- merge(category_table, order_caver1, by="ID")    
  # First column as Rownames (sample identification)    
  rownames(order_cat) <- order_cat[,1]
  order_cat <- order_cat[,-1]
  # Create a new column to assign a number to each category
  order_cat$category <- as.integer(order_cat[,1])
  # Create a new column to colour the plot depending on the category (level)        
  order_cat$color = "none"     
  # Create a palette of colors depending on the number of categories
  my_palette <- matrix(brewer.pal(category_number,"Paired"))
  #Assign one different color to each different category
  for (i in 1:category_number){
    order_cat[order_cat$category == i,]$color = my_palette[i,1]
  }
  # Create a new variable to get the same value to plot as barplot        
  order_cat$values <- "5"
  order_cat$values <- as.numeric(as.character(order_cat$values))
  # Sort the table by "num" column (same order as caver1)    
  order_cat1 <- order_cat[order(as.numeric(order_cat$num)),] 
  
  ### Cluster dendrogram  
  
  ## Instruction to combine plots: get two plots in the same view    
  par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
  # Create the dendrogram clustering plot (based on caver1 distances) 
  plot(caver1, hang=-1, labels = FALSE, axes = FALSE, ylab = "", xlab="", sub="") 
  # Create a barplot colored by category (the samples are in the same order as in caver1)    
  cluster_dendrogram <- barplot(order_cat1$values, col=order_cat1$color, border = NA, yaxt="n")
  
  # Save the plots as pdf.file
  pdf("cluster_dendrogram_plot.pdf")
  print(cluster_dendrogram)
  
  dev.off()
}

################################
#Taxonomy abundance + filtering#
################################
taxonomy_abundance <- function(taxonomy_table,abundance_value,individuals_value) {
  
  ##Function to calculate mean excluding 0 values
  nzmean <- function(a){
    mean(a[a!=0])
  }
  
  ##Function to calculate nº of 0
  zsum <- function(a){
    sum (a==0)
  }
  
  ##Function to calculate nº of non-0
  nsum <- function(a){
    sum (a!=0)
  }
  
  my_results=matrix(ncol = 4, nrow=ncol(taxonomy_table)) 
  
  ## Loop for each column (taxonomy) in the taxonomy table
  for (i in 1:ncol(taxonomy_table)) {
    
    #Calculate mean for each column
    aa = mean(taxonomy_table[,i])
    #Calculate number of non-zeros (individuals)
    bb = nsum(taxonomy_table[,i])
    #Calculate mean without taking into account the 0
    cc = nzmean(taxonomy_table[,i])
    #Calculate number of zeros 
    dd = zsum(taxonomy_table[,i])
    
    my_results[i,1] = aa
    my_results[i,2] = bb
    my_results[i,3] = cc
    my_results[i,4] = dd
  }
  
  
  # The column names from the original table = row names from the new table
  rownames(my_results) = colnames(taxonomy_table)
  # Give names to the columns of the new table
  colnames(my_results) = c("Mean","N_of_non-0", "Non-0_Mean", "N_of_0") 
  
  tax_parameters <- as.data.frame(my_results)
  tax_parameters[is.na(tax_parameters)] <- 0
  
  write.table(tax_parameters, file="~/summary_taxonomy_filtering.txt", sep = "\t", quote=F)
  
  ##Filtering by ABUNDANCE or N_OF_0
  ##Remove mean rows with less than 0.01% values OR less than 15 individuals (non_0)
  filtered_taxonomy <- tax_parameters[!(tax_parameters$Mean<abundance_value) & !(tax_parameters$`N_of_non-0`<individuals_value), , FALSE]
  #Remove non-necessary columns
  filtered_taxonomy[2:4] <- NULL
  
  ##Transpose the taxonomy_table to merge
  t_taxonomy_table <- as.data.frame(t(taxonomy_table))
  
  #Merge both tables <- only the rows in common are in the new table
  t_filtered_taxonomy_table <- merge(filtered_taxonomy, t_taxonomy_table, by="row.names")
  rownames(t_filtered_taxonomy_table) <- t_filtered_taxonomy_table[,1]
  t_filtered_taxonomy_table <- t_filtered_taxonomy_table[,-1]
  ##Remove Mean column
  t_filtered_taxonomy_table <- t_filtered_taxonomy_table[,-1]
  
  #Remove duplicated rows and keep the last one
  t_filtered_taxonomy_table = as.data.frame(t_filtered_taxonomy_table[!duplicated(t_filtered_taxonomy_table, fromLast = T), ])
  
  ##Transpose the matrix to get the filtered taxonomy table
  filtered_taxonomy_table <- as.data.frame(t(t_filtered_taxonomy_table))
  
  write.table(filtered_taxonomy_table, file="~/filtered_taxonomy.txt", sep = "\t", quote = F)
  
  ##Insterection plot    
  filtered_plot_data <- as.data.frame(setNames(replicate(2, numeric(0), simplify = F), letters[1:2])) 
  colnames(filtered_plot_data)[1:2] <- c("abundance","individuals")
  
  ##Add filter columns by abundance and individuals conditions
  # 1 value: removed taxa
  # 0 value: not removed taxa
  
  for (i in 1:nrow(tax_parameters))  {
    
    if (tax_parameters[i,"Mean"]<abundance_value) {
      filtered_plot_data[i,"abundance"] <- 1
    }
    else {
      # 1 value
      filtered_plot_data[i,"abundance"] <- 0
    }
    
    if (tax_parameters[i,"N_of_non-0"]<individuals_value) {
      filtered_plot_data[i,"individuals"] <- 1
    }
    else {
      filtered_plot_data[i,"individuals"] <- 0
    }
  }
  
  rownames(filtered_plot_data) <- rownames(tax_parameters)   
  
  ##Add column: overlap between NOT removed taxa by both methods
  filtered_plot_data$overlap <- "none"
  
  ##Add overlap column to color the plot with the intersection
  for (i in 1:nrow(filtered_plot_data)) {
    #If the taxa is NOT removed by the two methods
    if (filtered_plot_data[i,"abundance"]==1 & filtered_plot_data[i,"individuals"]==1) {
      #overlap column gets a random number between 20 and 40
      filtered_plot_data[i, "overlap"] <- sample(20:40, 1)
    }
    else {
      filtered_plot_data[i,"overlap"] <- 0
    }
  }
  ##Get rownames as a column to plot
  filtered_plot_data <- cbind(Row.names = rownames(filtered_plot_data), filtered_plot_data)
  rownames(filtered_plot_data) <- NULL
  colnames(filtered_plot_data)[1] <- "Name"
  
  library(UpSetR)        
  
  #Function to color the intersection in NOT removed taxonomies by both methods
  yes <- function(row, min, max){
    #Only count/color the taxa with overlap values
    newData <- (row["overlap"] <= max) & (row["overlap"] >= min)
  } 
  
  ## Save the plot as pdf.file
  pdf("intersection_filtering_plot.pdf")
  
  ##Intersection plot: marked in red the number of taxonomies removed by both methods
  intersection_filtering_plot <- upset(filtered_plot_data, sets = c("abundance","individuals"), main.bar.color = "black", queries = list(list(query =intersects, params = list("abundance","individuals")), list( query=yes, params = list(20,40), color = "#528EE7", active = T)))
  
  print(intersection_filtering_plot)
  
  dev.off()
}
