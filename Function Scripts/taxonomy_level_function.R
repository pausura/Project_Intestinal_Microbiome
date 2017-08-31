### Create taxonomy level table ###

# For a given txt file with taxonomies this function creates a new txt file with taxonomies at specific level. For example,
# returns a table only with species/orders/family... Needs a file with taxonomies in columns and taxonomy levels separated by "|".
# Taxonomy levels are numeric (for example: in taxonomy name k_|p_|c_|o_|f_|g_|s_|t_, species is level 7)

taxonomy_level_function <- function(taxonomy_table,taxonomy_level) {
  
  #Transpose taxonomy table (taxonomies in rows)
  t_tax_table <- as.data.frame(t(taxonomy_table))
  #Create a table to save the taxonomies at specific level. Same size as the original one
  loop_table <- as.data.frame(matrix(nrow = nrow(t_tax_table) , ncol = ncol(t_tax_table)))
  
  ## For each row/taxonomy in the table:
  for (i in 1:nrow(t_tax_table)) {
    # If the taxonomy has the desired number of levels: 
    if (count.fields(textConnection(row.names(t_tax_table[i,])), sep="|") == taxonomy_level){
      #print (paste0("Species found: ", row.names(tax_table[i,]))) ##Loop check
      
      # Only the rows that meet the condition are filled, the rest get NA values
      loop_table[i,] = t_tax_table[i,]
      
    }
  }
    #Taxonomy names as rownames (same as in the original transposed table)
    row.names(loop_table) = row.names(t_tax_table)
    #Sample names as colnames (same as in the original table)
    colnames(loop_table) = colnames(t_tax_table)
          
  ##Remove all rows with NA values  
    level_table <- na.omit(loop_table)

  ##For each row/taxonomy in the new level table
    for (i in 1:nrow(level_table)){
      #Save only the last part of the taxonomy name = taxonomy level name  
      name_taxonomy <- rownames(level_table)[i]
      new_taxonomy_name <- unlist(strsplit(name_taxonomy, split = "|", fixed= TRUE))[taxonomy_level]
      #Save the new name as rowname
      rownames(level_table)[i] <- new_taxonomy_name
    }
    
  ##Transpose the new level table to get taxonomies in columns and samples in rows
    t_level_table <- as.data.frame(t(level_table))
  
  #Export the new table
  write.table(t_level_table, file = "~/level_table.txt", quote = F, sep = "\t")

}
