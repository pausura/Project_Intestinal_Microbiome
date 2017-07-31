
    ## Create filtered files (reads over 10M)

getwd()
setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

  ## Filter samples <10M reads
reads_table <- read.table("~/Documents/Universitat/Holanda/Projecte/read_depth_per_sample_ibd.txt", sep = "\t", header = T, row.names= 1)
      reads_table$IBDFEC <- "none"
      reads_table$IBDFEC <- row.names(reads_table)

reads_over_10M <- as.data.frame(reads_table[with(reads_table, reads_table$Reads>10000000),])
      over_10M_reads <- write.table(reads_over_10M, file = "~/reads_over_10M.txt", quote = F, sep ="\t")
      reads_below_10M <- as.data.frame(reads_table[with(reads_table, reads_table$Reads<10000000),])
          #reads_below_10M <- write.table(reads_below_10M, file = "~/reads_below_10M.txt", quote = F, sep = "\t")

tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/Taxa data/IBD_taxonomy_DUDes.txt", sep = "\t", header = T, row.names = 1)
    t_tax_table <- as.data.frame(t(tax_table))
    
  #Filtered DUDes file
filtered_DUDes <- merge(t_tax_table, reads_over10M_table, by="row.names")
  rownames(filtered_DUDes) <- filtered_DUDes[,1]
  filtered_DUDes <- filtered_DUDes[,-1]
  filtered_DUDes$Reads <- NULL
      
      filtered_DUDes_tax <- write.table(filtered_DUDes, file = "~/filtered_tax_DUDes.txt", quote= F, sep="\t")
    
  #Filtered metaphlan file
filtered_metaphlan <- merge(t_tax_table, reads_over10M_table, by="row.names")
  rownames(filtered_metaphlan) <- filtered_metaphlan[,1]
  filtered_metaphlan <- filtered_metaphlan[,-1]
  filtered_metaphlan$Reads <- NULL
        
      filtered_metaphlan_tax <- write.table(filtered_metaphlan, file="~/filtered_tax_metaphlan.txt", quote = F, sep = "\t")
      
 
  ## Filtered metadata file (factors/phenotypes)  
new_metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/new_metadata_txt.txt", sep = "\t", header = T, row.names = 1)
  rownames(new_metadata) <- new_metadata[,3]
  new_metadata <- new_metadata[,-3]

  
  filtered_reads_metadata <- merge(new_metadata, reads_over_10M, by="row.names")
    rownames(filtered_reads_metadata) <- filtered_reads_metadata[,1]
    filtered_reads_metadata <- filtered_reads_metadata[,-1]
        filtered_reads_metadata$TotalReads <- NULL
        filtered_reads_metadata$IBDFEC <- NULL
    
filtered_total_metadata <- filtered_reads_metadata
  
  #Remove factors >66 NA
    cond_na <- sapply(filtered_total_metadata, function(col) sum(is.na(col)) < 66)
    filtered_total_metadata <- filtered_total_metadata[, cond_na, drop = FALSE]
  
  #Remove columns with only one level
    filtered_total_metadata <- filtered_total_metadata[, sapply(filtered_total_metadata, function(col) length(unique(col))) > 1]
    
  ##Add Group column
    
    intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names=1)
    
  filtered_total_metadata <- merge(intestinal_groups, filtered_total_metadata, by="row.names")
    rownames(filtered_total_metadata) <- filtered_total_metadata[,1]
    filtered_total_metadata <- filtered_total_metadata[,-1]
          filtered_total_metadata$Sex2 <- NULL
          filtered_total_metadata$WGSProjectID2 <- NULL
          filtered_total_metadata$MedicationMesalazines <- NULL
          filtered_total_metadata$MedicationPPI <- NULL
    
    
  filtered_metadata <- write.table(filtered_total_metadata, file = "~/filtered_metadata.txt", quote = F, sep="\t")
    
    
    
    
    