  ## Taxonomy Analysis Function ##

  #Open the taxonomy table
tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
    
taxonomy_abundance(tax_table,0.01,15)

  #taxonomy_table
  #abundance_value: value for filtering by relative abundance mean
  #individuals_value: value for filtering by number of minimum individuals

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











