
### Principal coordinate analysis (PCoA) ###

getwd()

setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")


## Open species table (created in shannon_index script)
species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", sep = "\t", header = T, row.names = 1)
      ##Transpose the filum table for the diversity function
        t_species_table <- as.data.frame(t(species_table))

# Open reads table
reads_table <- read.table("~/Documents/Universitat/Holanda/Projecte/reads_over_10M.txt", sep = "\t", header = T, row.names = 1)
    reads_table$IBDFEC <- NULL
    
  library(vegan)
alpha <- as.data.frame(diversity(t_species_table,index="shannon"))

 #Open intestinal groups table
intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep="\t", header = T, row.names = 1)


## Open tax table to create filum table

tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
  t_tax_table <- as.data.frame(t(tax_table))
loop_table <- as.data.frame(matrix(nrow = nrow(t_tax_table) , ncol = ncol(t_tax_table)))

## Count only "|p_" --> 2 fields

for (i in 1:nrow(t_tax_table)) {
  
  if (count.fields(textConnection(row.names(t_tax_table[i,])), sep="|") == 2){
    #print (paste0("Species found: ", row.names(tax_table[i,]))) ##Loop check
    
    loop_table[i,] = t_tax_table[i,]
    
  }
  
}

# Give row names to the new table
row.names(loop_table) = row.names(t_tax_table)

# Give column names to  new table
colnames(loop_table) = colnames(t_tax_table)

##Remove all rows with NA values  
filum_table <- na.omit(loop_table)

fl_table <- write.table(filum_table, file = "~/filum_table_DUDes.txt", quote = F, sep = "\t")

  
library(vegan)
library(ggplot2)

      #Generate a distance matrix - Bray method
  beta <- vegdist(t_species_table, method="bray")
      #cmdscale -> multidimensional scaling of a data matrix (distance matrix)
  my_pcoa <- as.data.frame(cmdscale(beta, k = 4))
    colnames(my_pcoa)[1:4] <- c("PCoA1","PCoA2","PCoA3","PCoA4")

pcoa_table <- write.table(my_pcoa, file="~/filtered_pcoa_table_DUDes.txt", quote = F, sep = "\t")
  
    # Calculate the mean for each phyla to see the most abundants
t_filum_table <- as.data.frame(t(filum_table))
  summary(t_filum_table)
        ##Actinobacteria - Bacteroidetes - Firmicutes

## Create a full table
plot_table1 <- merge(t_filum_table, my_pcoa, by="row.names")
    rownames(plot_table1) <- plot_table1[,1]
    plot_table1 <- plot_table1[,-1]

plot_table2 <- merge(plot_table1, intestinal_groups, by="row.names")
    rownames(plot_table2) <- plot_table2[,1]
    plot_table2 <- plot_table2[,-1]

plot_table3 <- merge(plot_table2, alpha, by="row.names")
    rownames(plot_table3) <- plot_table3[,1]
    plot_table3 <- plot_table3[,-1]

final_plot_table <- merge(plot_table3, reads_table, by = "row.names")
   rownames(final_plot_table) <- final_plot_table[,1]
   final_plot_table <- final_plot_table[,-1]

final_table <- write.table(final_plot_table, file = "~/filtered_final_pcoa_table.txt", quote = F, sep= "\t")   
   

  ##### PCoA plot 

## Per filum

  # Variable that contains the colors used in the plot
  my_col=c("#0000FF","#62A4D1","#5BE55B","#FFF000", "#FF0000")
   
     ## Actinobacteria
a_plot <- ggplot (final_plot_table, aes(x=PCoA1, y=PCoA2, geom="blank", colour=final_plot_table$`k__Bacteria|p__Actinobacteria`)) + geom_point () + scale_color_gradientn(colours = my_col, "Actinobacteria") + theme_classic() + labs(x="PCoA1", y="PCoA2")
a_plot

    ## Bacteroidetes
b_plot <- ggplot (final_plot_table, aes(x=PCoA1, y= PCoA2, geom = "blank", colour = final_plot_table$`k__Bacteria|p__Bacteroidetes`)) + geom_point() + scale_color_gradientn(colours = my_col, "Bacteroidetes") + theme_classic() + labs(x = "PCoA1", y="PCoA2")
b_plot

    ## Firmicutes
f_plot <- ggplot (final_plot_table, aes(x=PCoA1, y= PCoA2, geom = "blank", colour = final_plot_table$`k__Bacteria|p__Firmicutes`)) + geom_point() + scale_color_gradientn(colours = my_col, "Firmicutes") + theme_classic() + labs(x = "PCoA1", y = "PCoA2")
f_plot


### Shannon Index
my_col2 <- c("#000080","#00ffff")
si_plot <- ggplot(final_plot_table, aes(x=PCoA1, y=PCoA2, geom = "blank", colour = final_plot_table$`diversity(t_species_table, index = "shannon")`)) + geom_point() + scale_color_gradientn(colours = my_col2, "Shannon Index") + theme_classic() + labs(x = "PCoA1", y = "PCoA2")
si_plot


### Per number of reads
final_plot_table$color="black"
    #Only the samples with <10000000 reads are in red
final_plot_table[final_plot_table$Reads<10000000,]$color="red"
reads_plot <- ggplot (final_plot_table, aes(x=PCoA1, y=PCoA2, geom="blank", colour=final_plot_table$color)) + geom_point () + scale_color_identity("Reads", breaks=final_plot_table$color, labels=final_plot_table$color, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
reads_plot

### Per group

  ## Total
final_plot_table$color="none"
final_plot_table[final_plot_table$Group=="normal",]$color="black"
final_plot_table[final_plot_table$Group=="intermediate",]$color="#2F2BFF"
final_plot_table[final_plot_table$Group=="small intestine",]$color="red"

total_groups_plot <- ggplot (final_plot_table, aes(x=PCoA1, y=PCoA2, geom="blank", colour=color)) + geom_point () + scale_color_identity("Microbiome", breaks=final_plot_table$color, labels=final_plot_table$Group, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
total_groups_plot

## Normal
final_plot_table[final_plot_table$Group=="normal",]$color="#ffa500"
final_plot_table[final_plot_table$Group=="intermediate",]$color="black"
final_plot_table[final_plot_table$Group=="small intestine",]$color="black"

normal_plot <- ggplot (final_plot_table, aes(x=PCoA1, y=PCoA2, geom="blank", colour=color)) + geom_point () + scale_color_identity("Microbiome", breaks=final_plot_table$color, labels=final_plot_table$Group, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
normal_plot

## Intermediate
final_plot_table[final_plot_table$Group=="normal",]$color="black"
final_plot_table[final_plot_table$Group=="intermediate",]$color="#ffa500"
final_plot_table[final_plot_table$Group=="small intestine",]$color="black"

intermediate_plot <- ggplot (final_plot_table, aes(x=PCoA1, y=PCoA2, geom="blank", colour=color)) + geom_point () + scale_color_identity("Microbiome", breaks=final_plot_table$color, labels=final_plot_table$Group, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
intermediate_plot

## Small Intestine
final_plot_table[final_plot_table$Group=="normal",]$color="black"
final_plot_table[final_plot_table$Group=="intermediate",]$color="black"
final_plot_table[final_plot_table$Group=="small intestine",]$color="#ffa500"

smallintestine_plot <- ggplot (final_plot_table, aes(x=PCoA1, y=PCoA2, geom="blank", colour=color)) + geom_point () + scale_color_identity("Microbiome", breaks=final_plot_table$color, labels=final_plot_table$Group, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
smallintestine_plot
  