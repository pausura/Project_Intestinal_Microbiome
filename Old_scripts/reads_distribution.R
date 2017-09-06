  
  ### Plot reads distribution per sample ###

getwd()

setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

library(ggplot2)

#Num reads
reads_table <- read.table("~/Documents/Universitat/Holanda/Projecte/reads_over_10M.txt", sep = "\t", header = T, row.names = 1)
  reads_table$IBDFEC <- NULL
intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)

  ##Merge reads table with intestinal groups table to colour by groups
  
reads_table2 <- merge(reads_table, intestinal_groups, by="row.names")
  colnames(reads_table2)[1] <- "IBDFEC"  
  
  ## Plot colour by intestinal group
  reads_table2$color_merged="none"
  reads_table2[reads_table2$Group=="normal",]$color_merged="#9D9DAB"
  reads_table2[reads_table2$Group=="intermediate",]$color_merged="#2F2BFF"
  reads_table2[reads_table2$Group=="small intestine",]$color_merged="red"

reads_distribution_plot <- ggplot(reads_table2, aes(x=reorder(as.factor(IBDFEC), -Reads), y=Reads, colour=reads_table2$color_merged)) + geom_bar(stat="identity") + scale_color_identity("Nº of Reads", breaks=reads_table2$color_merged, labels=reads_table2$Group, guide = "legend") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
                                                                                                                                                                                                                                                                                       
reads_distribution_plot
