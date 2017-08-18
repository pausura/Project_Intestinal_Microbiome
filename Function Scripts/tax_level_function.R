### Create taxonomy level table ###

  #Open the taxonomy table
tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)

level_function(tax_table,2)

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
    
    
    
  write.table(t_level_table, file = "~/level_table.txt", quote = F, sep = "\t")

}
