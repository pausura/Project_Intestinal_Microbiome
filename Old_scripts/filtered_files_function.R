 
  ## Create filtered files by SAMPLE quality control ##

#Input files
#1.Taxonomy_table
#2.Metadata_table
#3.Parameter_filter file: file with the parameter to filter
#4.Taxonomy_threshold: threshold (value) to filter taxonomy_table based on parameter_filter file
#5.Metadata_threshold: threshold (value) to filter metadata_table based on parameter_filter file



filter_files_by_sample_qc <- function(tax_input, metadata_input, parameter_file, taxonomy_threshold, metadata_threshold) {
   
  ##Filter the parameter_file by the given threshold
  taxonomy_threshold_table <- as.data.frame(parameter_file[with(parameter_file, parameter_file[,1]>taxonomy_threshold),])
  #Write table
  write.table(taxonomy_threshold_table, file = "./taxonomy_threshold_table.txt", quote = F, sep ="\t")
    
  #Merge the taxonomy_threshold_table with tax_input file  
  filtered_tax <- merge(taxonomy_threshold_table, tax_input, by="row.names")
  rownames(filtered_tax) <- filtered_tax[,1]
  filtered_tax <- filtered_tax[,-1]
  #Remove parameter_file column and sample_names (repeated)  
  filtered_tax[1:2] <- NULL
        
  #Write filtered_taxonomy  
  write.table(filtered_tax, file = "./filtered_tax.txt", quote= F, sep="\t")
      
  ##Filter the metadata_file by the given taxonomy_threshold
  filtered_metadata <- merge(taxonomy_threshold_table, metadata_input, by="row.names")
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
      
  #Write filtered_metadata
  write.table(filtered_total_metadata, file = "./filtered_metadata.txt", quote = F, sep="\t")
    
}