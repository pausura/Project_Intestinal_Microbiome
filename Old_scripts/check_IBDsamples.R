
  ## Check IBDFEC numbers in PCoA analysis

pcoa_table <- read.table("~/Documents/Universitat/Holanda/Projecte/DUDes_results/final_plot_table_DUDes.txt", sep = "\t", header = T, row.names = 1)

  # Check small intestine samples
small_intestine_table <- subset(pcoa_table, pcoa_table$Group=="small intestine")
      #Small intestine with PCoA1 values < 0 (normal zone)
      pcoa_num <- subset(small_intestine_table, small_intestine_table$PCoA1<0)
 
  # Check normal/intermediate samples
normal_table <- subset(pcoa_table, pcoa_table$Group=="normal"|pcoa_table$Group=="intermediate")
      #Normal/Intermediate with PCoA1 values >0.49 (small intestine zone)
      pcoa_num1 <- subset(normal_table, normal_table$PCoA1>0.49)

  ## Check samples with less than 10M reads
reads_table <- subset(pcoa_table, pcoa_table$Reads<10000000)
    summary(reads_table$Group)
