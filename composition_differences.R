
  ### Composition differences ### (Which differences?)

getwd()
setwd("/Users/paulasuredahorrach/Documents/Universitat/Holanda/Projecte")

intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)


# Filum level
    # Open filum_table
filum_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/filum_table_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
t_filum_table <- as.data.frame(t(filum_table))

filum_groups <- merge(t_filum_table, intestinal_groups, by = "row.names")
      rownames(filum_groups) <- filum_groups[,1]
      filum_groups <- filum_groups[,-1]

  ## Subset the table by categories    
library(psych)
       
##Intermediate
  intermediate_filum <- subset(filum_groups, filum_groups$Group=="intermediate")
  intermediate_filum <- intermediate_filum[,-13]
      ## Calculate mean
            resum_intermediate <- describe(intermediate_filum)
          ## Find the 4 rows with the highest values
            top4_table_intermediate <- resum_intermediate[order(resum_intermediate$mean, decreasing = T)[1:4],]
          #Sum values and create an other category: others
            sum_top4_intermediate <- sum(top4_table_intermediate$mean)
            others_intermediate <- (100-sum_top4_intermediate)
    
      ## Create a new table and add a row for "others" category
  intermediate_results <- subset(top4_table_intermediate, select = "mean")
      colnames(intermediate_results) <- "Intermediate"
  intermediate_results <- rbind(intermediate_results, others_intermediate)
      rownames(intermediate_results)[5] <- "others"    
      
##Normal
  normal_filum <- subset(filum_groups, filum_groups$Group=="normal")
  normal_filum <- normal_filum[,-13]
       ## Calculate mean
            resum_normal <- describe(normal_filum)
          ## Find the 4 rows with the highest values
            top4_table_normal <- resum_normal[order(resum_normal$mean, decreasing = T)[1:4],]
          #Sum values and create an other category: others
            sum_top4_normal <- sum(top4_table_normal$mean)
            others_normal <- (100-sum_top4_normal)
  
        ## Create a new table and add a row for "others" category
  normal_results <- subset(top4_table_normal, select = "mean")
      colnames(normal_results) <- "Normal"
  normal_results <- rbind(normal_results, others_normal)
      rownames(normal_results)[5] <- "others" 

##Small Intestine
  small_filum <- subset(filum_groups, filum_groups$Group=="small intestine")
  small_filum <- small_filum[,-13]
        ## Calculate mean
            resum_small <- describe(small_filum)
          ## Find the 4 rows with the highest values
            top4_table_small <- resum_small[order(resum_small$mean, decreasing = T)[1:4],]
          #Sum values and create an other category: others
            sum_top4_small <- sum(top4_table_small$mean)
            others_small <- (100-sum_top4_small)
  
        ## Create a new table and add a row for "others" category
  small_results <- subset(top4_table_small, select = "mean")
      colnames(small_results) <- "Small_Intestine"
  small_results <- rbind(small_results, others_small)
      rownames(small_results)[5] <- "others"       
      
 ## Merge the results for each categories in the same table
       
    normal_intermediate <- merge(normal_results, intermediate_results, by="row.names")
          rownames(normal_intermediate) <- normal_intermediate[,1]
          normal_intermediate <- normal_intermediate[,-1]
    
    filum_group_table <- merge(normal_intermediate, small_results, by = "row.names")  
          rownames(filum_group_table) <- filum_group_table[,1]
          filum_group_table <- filum_group_table[,-1]
   
  #### Stacked barplot ####

  library(reshape2) 
  library(ggplot2)
          
    filum_group_table$bacteria=row.names(filum_group_table)
          row.names(filum_group_table)=NULL
    my_table=melt(filum_group_table)

filum_plot <- ggplot (my_table, aes(x=variable, y=value)) + geom_bar (aes(fill = bacteria), stat = "identity") + theme_classic() + xlab("Group") + ylab("relative_abundance")
filum_plot


# Order level
    # Create a order_table
tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
  t_tax_table <- as.data.frame(t(tax_table))
loop_table <- as.data.frame(matrix(nrow = nrow(t_tax_table) , ncol = ncol(t_tax_table)))

        ## Count only "|o_" --> 4 fields
        for (i in 1:nrow(t_tax_table)) {
          if (count.fields(textConnection(row.names(t_tax_table[i,])), sep="|") == 4){
            #print (paste0("Species found: ", row.names(tax_table[i,]))) ##Loop check
            loop_table[i,] = t_tax_table[i,]
          }
        }
        
        # Give row names to the new table
        row.names(loop_table) = row.names(t_tax_table)
        # Give column names to  new table
        colnames(loop_table) = colnames(t_tax_table)
        # Remove all rows with NA values  
        order_table <- na.omit(loop_table)
        
  or_table <- write.table(order_table, file = "~/order_table_DUDes.txt", quote = F, sep = "\t")
        
  t_order_table <- as.data.frame(t(order_table))
  
  order_groups <- merge(t_order_table, intestinal_groups, by = "row.names")
      rownames(order_groups) <- order_groups[,1]
      order_groups <- order_groups[,-1]
  
  ## Subset the table by categories    
  library(psych)
  
  ##Intermediate
    intermediate_order <- subset(order_groups, order_groups$Group=="intermediate")
    intermediate_order <- intermediate_order[,-48]
        ## Calculate mean
           resum_intermediate_order <- describe(intermediate_order)
  
  ## Create a new table for mean column
    intermediate_results_order <- subset(resum_intermediate_order, select = "mean")
    colnames(intermediate_results_order) <- "Intermediate"
  
  ##Normal
    normal_order <- subset(order_groups, order_groups$Group=="normal")
    normal_order <- normal_order[,-48]
      ## Calculate mean
        resum_normal_order <- describe(normal_order)
  
  ## Create a new table for mean column
    normal_results_order <- subset(resum_normal_order, select = "mean")
    colnames(normal_results_order) <- "Normal"
  
  ##Small Intestine
    small_order <- subset(order_groups, order_groups$Group=="small intestine")
    small_order <- small_order[,-48]
      ## Calculate mean
          resum_small_order <- describe(small_order)
  
  ## Create a new table and add a row for "others" category
    small_results_order <- subset(resum_small_order, select = "mean")
    colnames(small_results_order) <- "Small_Intestine"
  
    
  ## Merge the results for each categories in the same table
  
  normal_intermediate_order <- merge(normal_results_order, intermediate_results_order, by="row.names")
    rownames(normal_intermediate_order) <- normal_intermediate_order[,1]
    normal_intermediate_order <- normal_intermediate_order[,-1]
  
  order_group_table <- merge(normal_intermediate_order, small_results_order, by = "row.names")  
    rownames(order_group_table) <- order_group_table[,1]
    order_group_table <- order_group_table[,-1]
  
  #### Stacked barplot ####
  
  order_group_table$bacteria=row.names(order_group_table)
      row.names(order_group_table)=NULL
  my_table_order=melt(order_group_table)
  
  order_plot <- ggplot (my_table_order, aes(x=variable, y=value)) + geom_bar (aes(fill = bacteria), stat = "identity") + theme(legend.position="top", legend.text = element_text(size=4)) + xlab("Group") + ylab("relative_abundance")  
  order_plot
  
