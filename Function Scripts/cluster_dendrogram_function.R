### Clustering and Heatmap PCoA ###

## Open tax_level table
species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", header = T, row.names = 1, check.names = F)
      ##Need taxonomy in columns, transpose if it's necessary      
## Open category file
intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)
      
clustering_dendrogram(species_table, intestinal_groups, 3)

clustering_dendrogram <- function(tax_level_table, category_table, category_number) {
  
  #Packages needed
    library(vegan)
    library(gclus)
    library(data.table)
  
   ## Calculate distance matrix - Bray method
        beta <- vegdist(tax_level_table, method="bray")
    # Calculate AVERAGE distance in beta
        caver <- hclust(beta, method="aver")
    # Reorder the matrix generated    
        caver1 <- reorder.hclust(caver, beta)
    # Get the labels(sample identification) in the same order as are in caver1    
        order_caver1 = as.data.frame(caver1$labels[c(caver1$order)])  
        setDT(order_caver1, keep.rownames = TRUE)[]
          # Give names to the columns and use ID column to put rownames
            colnames(order_caver1)[1] <- "num"
            colnames(order_caver1)[2] <- "ID"
   
   ## Before merge: create a new colum in category_table with the ID
        category_table$ID <- rownames(category_table)    
            
   ## Merge order_caver1 with the category_table        
        order_cat <- merge(category_table, order_caver1, by="ID")    
          # First column as Rownames (sample identification)    
            rownames(order_cat) <- order_cat[,1]
            order_cat <- order_cat[,-1]
    # Create a new column to assign a number to each category
        order_cat$category <- as.integer(order_cat[,1])
    # Create a new column to colour the plot depending on the category (level)        
        order_cat$color = "none"     
    # Create a palette of colors depending on the number of categories
        my_palette <- matrix(brewer.pal(category_number,"Paired"))
    #Assign one different color to each different category
         for (i in 1:category_number){
           order_cat[order_cat$category == i,]$color = my_palette[i,1]
         }
    # Create a new variable to get the same value to plot as barplot        
        order_cat$values <- "5"
        order_cat$values <- as.numeric(as.character(order_cat$values))
    # Sort the table by "num" column (same order as caver1)    
        order_cat1 <- order_cat[order(as.numeric(order_cat$num)),] 

### Cluster dendrogram  
        
   ## Instruction to combine plots: get two plots in the same view    
        par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
    # Create the dendrogram clustering plot (based on caver1 distances) 
        plot(caver1, hang=-1, labels = FALSE, axes = FALSE, ylab = "", xlab="", sub="") 
    # Create a barplot colored by category (the samples are in the same order as in caver1)    
        cluster_dendrogram <- barplot(order_cat1$values, col=order_cat1$color, border = NA, yaxt="n")
    
    # Save the plots as pdf.file
        pdf("cluster_dendrogram_plot.pdf")
        print(cluster_dendrogram)
        
        dev.off()
}

  











