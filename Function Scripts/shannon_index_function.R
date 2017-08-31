### Calculate Shannon Index ###

# For a given taxonomy file at specific level (for example: species) this function calculates alpha diversity for each sample
# and creates a violin plot to represent differences between sample categories. Need a taxonomy file with taxonomies in columns
# and samples in rows and a sample category file. 

shannon_index_function <- function(taxonomy_level_table, category_table) {
  
  ##Required packages - shannon diversity/ggplot
  library(vegan)
  library(ggplot2)
  library(RColorBrewer)
  
  ## Calculate shannon index (alpha) for each sample in the taxonomy file: diversity function from vegan package
     alpha <- as.data.frame(diversity(taxonomy_level_table, index="shannon"))
      colnames(alpha)[1] <- "alpha_diversity"
  ## Divide alpha results by categories
      # Merge alpha with category table
          my_table <- merge(category_table, alpha, by="row.names")
              rownames(my_table) <- my_table[,1]
              my_table <- my_table[,-1]
      # Calculate the number of categories
        category_number <- nlevels(my_table[,1]) #First column because in merge the category table is the first one
      # Create a new column to assign a number to each category
        my_table$category <- as.integer(my_table[,1])    
      # Create a new column to colour the plot depending on the category (level)        
        my_table$color = "none" 
      # Create a palette of colors depending on the number of categories
        my_palette <- matrix(brewer.pal(category_number,"Set1"))
      # For loop to assign one different color to each different category in the new color column
        for (i in 1:category_number){
            my_table[my_table$category == i,]$color = my_palette[i,1]
          }
  
  ### Create a violin plot
shannon_index_plot <- ggplot(my_table, aes(x=category, y=alpha_diversity, fill = my_table$color)) + labs (y="Shannon Diversity Index", x="Category") + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1) + theme_classic() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 1, size=16,color="black")) + scale_color_identity("All_categories", breaks = my_table$color, labels= my_table$category, guide = "legend")

  ##Save the plot as pdf.file
  pdf("shannon_plot.pdf")
  print(shannon_index_plot)

dev.off()


}


