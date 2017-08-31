
##Overlap##

#NORMAL SMALL
unique_group_normal_small <- read.table("~/Desktop/unique_group_normal_small.txt", header= F)
unique_group_normal_small <- as.data.frame(t(unique_group_normal_small))
colnames(unique_group_normal_small)[1] <- "Taxonomy"

Groupsmall.intestine.Variable_filtered <- read.table("~/Documents/Universitat/Holanda/logistic_regression_adonis_data/normal_small_results/filtered_normal_small/Groupsmall.intestine.Variable_filtered.txt", header = T, sep = "\t")
colnames(Groupsmall.intestine.Variable_filtered)[1] <- "Taxonomy"

unique_group_normal_small <- merge(unique_group_normal_small, Groupsmall.intestine.Variable_filtered, by = "Taxonomy")

write.table(unique_group, file = "~/effect_group_normal_small.txt", sep = "\t", quote = F)

#NORMAL INTERMEDIATE
unique_group_normal_intermediate <- read.table("~/Desktop/unique_group_normal_intermediate.txt", header= F)
colnames(unique_group_normal_intermediate)[1] <- "Taxonomy" 

Groupnormal.Variable_filtered <- read.table("~/Documents/Universitat/Holanda/logistic_regression_adonis_data/normal_intermediate_results/filtered_normal_intermediate/Groupnormal.Variable_filtered.txt", header = T, sep = "\t")
colnames(Groupnormal.Variable_filtered)[1] <- "Taxonomy"

unique_group_normal_intermediate <- merge(unique_group_normal_intermediate, Groupnormal.Variable_filtered, by = "Taxonomy")

write.table(unique_group_normal_intermediate, file = "~/effect_group_normal_intermediate.txt", sep = "\t", quote = F)


#SMALL INTERMEDIATE
unique_group_intermediate_small <- read.table("~/Desktop/unique_group_intermediate_small.txt", header = F)
colnames(unique_group_intermediate_small)[1] <- "Taxonomy"

Groupsmall.intestine.Variable_filtered <- read.table("~/Documents/Universitat/Holanda/logistic_regression_adonis_data/intermediate_small_results/filtered_intermediate_small/Groupsmall.intestine.Variable_filtered.txt", header = T, sep = "\t")
  colnames(Groupsmall.intestine.Variable_filtered)[1] <- "Taxonomy"

unique_group_intermediate_small <- merge(unique_group_intermediate_small, Groupsmall.intestine.Variable_filtered, by = "Taxonomy")

write.table(unique_group_intermediate_small, file = "~/effect_group_intermediate_small.txt", sep= "\t", quote = F)


#Plot

plot_table <- read.table("~/Documents/Universitat/Holanda/Projecte/plot_table_species.txt", sep = "\t", header = T, row.names = 1, check.names = F)
group <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_content_group.txt", sep ="\t", header = T, row.names = 1, check.names = F)

group_plot_table <- merge(plot_table, group, by = "row.names")
rownames(group_plot_table) <- group_plot_table[,1]
group_plot_table <- group_plot_table[,-1]

plot_matrix = matrix(ncol = 4, nrow = 15)

normal_plot_table <- subset(group_plot_table, group_plot_table$Group=="normal")
normal_plot_table$Group <- NULL

small_plot_table <- subset(group_plot_table, group_plot_table$Group=="small intestine")
small_plot_table$Group <- NULL

intermediate_plot_table <- subset(group_plot_table, group_plot_table$Group=="intermediate")
intermediate_plot_table$Group <- NULL

group_plot_table$Group <- NULL  

for (k in 1:ncol(group_plot_table)) {
  
  percentage_normal = (nsum(normal_plot_table[,k])/nrow(normal_plot_table))
  percentage_intermediate = (nsum(intermediate_plot_table[,k])/nrow(intermediate_plot_table))
  percentage_small = (nsum(small_plot_table[,k])/nrow(small_plot_table))
  
  name_column <- colnames(group_plot_table)[k]
  name_column <- unlist(strsplit(name_column, split='|', fixed=TRUE))
  new_name_column <- tail(name_column, n = 1)
  
  plot_matrix[k,1] = new_name_column
  plot_matrix[k,2] = percentage_normal
  plot_matrix[k,3] = percentage_intermediate
  plot_matrix[k,4] = percentage_small
  
  
}

colnames(plot_matrix) <- c("Taxonomy", "%_normal", "%_intermediate", "%_small")

library(ggplot2)
library(reshape2)
plot_matrix=as.data.frame(plot_matrix)
test_plot1 = melt(plot_matrix, id=1)

ggplot(test_plot1, aes(x=Taxonomy, y=as.numeric(value), group=variable, fill=variable)) + geom_bar(stat="identity",position="dodge") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))



  
  
  
  
  
  
  

