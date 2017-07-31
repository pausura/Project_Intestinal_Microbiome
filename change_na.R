
    ## Change NA from metadata

  getwd()
  setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

    new_metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep = "\t", header = T, row.names = 1)
    intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)
    
    
    total_metadata <- merge(new_metadata, intestinal_groups, by="row.names")
        rownames(total_metadata) <- total_metadata[,1]
        total_metadata <- total_metadata[,-1]
    
  for(i in 1:ncol(total_metadata)){
    
    if (is.numeric(total_metadata[,i])){
  
      total_metadata[is.na(total_metadata[,i]), i] <- median(total_metadata[,i], na.rm = TRUE)
      
    }
    
    
    else {
      
      total_metadata[,i] <- as.numeric(factor(total_metadata[,i]))
      
      
    }
           }
 
 write.table(total_metadata, file="~/replaced_na_filtered_metadata.txt", quote = F, sep="\t")
 
 