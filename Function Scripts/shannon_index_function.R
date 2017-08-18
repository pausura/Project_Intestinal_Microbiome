### Calculate Shannon Index ###

# Open species table
  species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
# Open group table
  intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)

shannon_function(species_table, intestinal_groups)

shannon_function <- function(level_table, group_table) {
  
  ##Required packages - shannon diversity/ggplot
  library(vegan)
  library(ggplot2)
  
  ## Calculate shannon (alpha) for each sample
   alpha <- as.data.frame(diversity(level_table, index="shannon"))
      colnames(alpha)[1] <- "alpha_diversity"
  ## Divide alpha results by intestinal groups
      #Merge alpha with intestinal groups table
          group_taxa <- merge(group_table, alpha, by="row.names")
              rownames(group_taxa) <- group_taxa[,1]
              group_taxa <- group_taxa[,-1]
          #colnames(group_taxa)[2] <- "diversity"
      # Calculate the number of categories
        category_number <- nlevels(group_taxa[,1])
      # Create a new column to assign a number to each category
        group_taxa$category <- as.integer(group_taxa[,1])    
      # Create a new column to colour the plot depending on the category (level)        
        group_taxa$color = "none" 
      # Create a palette of colors depending on the number of categories
        my_palette <- matrix(brewer.pal(category_number,"Set1"))
      # For loop to assign one different color to each different category in a new column
        for (i in 1:category_number){
            group_taxa[group_taxa$category == i,]$color = my_palette[i,1]
          }
  
  ### Create a violin plot
shannon_plot <- ggplot(group_taxa, aes(x=category, y=alpha_diversity, fill = group_taxa$color)) + labs (y="Shannon Diversity Index", x="Category") + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1) + theme_classic() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 1, size=16,color="black")) + scale_color_identity("All_categories", breaks = group_taxa$color, labels= group_taxa$category, guide = "legend")

  ##Save the plot as pdf.file
  pdf("shannon_plot.pdf")
  print(shannon_plot)

dev.off()


}


