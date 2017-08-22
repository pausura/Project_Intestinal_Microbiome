 
  ## Create filtered files by SAMPLE quality ##

  # input_file: file to filtrate 
        # tax_input: taxa 
        # metadata_input: metadata
  # parameter_file: file with which to filter
  # parameter_threshold: threshold to filter tax_data based on parameter_file
  # metadata_threshold: number max of NA to tolerate in metadata file

    #Open parameter_file
reads_table <- read.table("~/Documents/Universitat/Holanda/Projecte/read_depth_per_sample_ibd.txt", sep = "\t", header = T, row.names= 1)
        #Add rownames as column
          reads_table$sample_names <- "none"
          reads_table$sample_names <- row.names(reads_table)

    #Open tax_input file
tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/taxa_data/IBD_taxonomy_DUDes.txt", sep = "\t", header = T, row.names = 1)
      ## Need taxonomy in columns and samples in rows
         tax_table <- t(tax_table)
                  
    #Open metadata_input file                
metadata_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/new_metadata_txt.txt", sep = "\t", header = T, row.names = 1)
                  
       filter_data(tax_table, metadata_table, reads_table, 10000000, 66)           

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


