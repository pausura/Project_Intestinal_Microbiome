
  ## Statistical analysis - alfa-diversity/Shannon Index ##

getwd()
setwd("/Users/paulasuredahorrach/Documents/Universitat/Holanda/Projecte")

# Open diversity_table
shannon_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/alpha_diversity_DUDes.txt", sep = "\t", header = T, row.names = 1)
# Open intestinal_groups table
intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)

## Merge both tables
shannon_groups_table <- merge(shannon_table, intestinal_groups, by = "row.names")
    rownames(shannon_groups_table) <- shannon_groups_table[,1]
    shannon_groups_table <- shannon_groups_table[,-1]
    colnames(shannon_groups_table)[1] <- "alfa_diversity"

## Specify the order of factor levels
   library(dplyr)
   library(FSA)

shannon_groups_table1 <- mutate(shannon_groups_table, Group = factor(Group, levels = unique(Group)))
summ <- Summarize(alfa_diversity~Group, data = shannon_groups_table1)
summ

## Histograms and values across groups
   library(lattice)

hist_plot <- histogram (~ alfa_diversity | Group, data=shannon_groups_table1, layout=c(1,3)) #layout=c(1,3) -> columns and rows of individual plots
hist_plot

box_plot <- boxplot(alfa_diversity ~ Group, data = shannon_groups_table1, ylab="alfa_diversity", xlab="Group")
box_plot

### Kruskal-Wallis test ###

# For all groups
test1 <- kruskal.test(alfa_diversity ~ Group, data = shannon_groups_table1)
test1
            # DUDes: p-value <2.2e-16 --> significative differences between the three groups

### Mann-Whitney-Wilcoxon test ###
  # Multiple comparisons (per parelles de grups)

    ## Intermediate - Normal (in)
intermediate_normal <- shannon_groups_table[!(shannon_groups_table$Group=="small intestine"),]
in_test <- wilcox.test(alfa_diversity ~ Group,  data = intermediate_normal, paired = FALSE) 
in_test       ## DUDes: p-value = 1.836e-10 --> significative differences - intermediate/normal are non-identical populations
  
    ## Intermediate - Small Intestine (is)
intermediate_small <- shannon_groups_table[!(shannon_groups_table$Group=="normal"),]
is_test <- wilcox.test(alfa_diversity ~ Group,  data = intermediate_small, paired = FALSE) 
is_test       ## DUDes: p-value = 7.157e-13 --> significative differences - intermediate/small are non-identical populations

    ## Normal - Small Intestine (ns)
normal_small <- shannon_groups_table[!(shannon_groups_table$Group=="intermediate"),]
ns_test <- wilcox.test(alfa_diversity ~ Group, data = normal_small, paired = FALSE)
ns_test      ## DUDes: p-value = 2.2e-16 --> significative differences - normal/small are non-identical populations  

  ### The three groups are statistical different between them


