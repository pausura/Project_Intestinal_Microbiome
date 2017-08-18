
#Filum level

filum_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/filum_table_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)
groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1, check.names = F)

abundance_taxonomy <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/relative_abundance_per_groups_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)

p_a_table <- filum_table

for (i in 1:ncol(filum_table)) {
  
  for (j in 1:nrow(filum_table)) {
    
    if (filum_table[j,i]>0) {
      
      p_a_table[j,i] = 1
      
    }

  }
}

group_pa_table <- merge(p_a_table, groups, by = "row.names")
  rownames(group_pa_table) <- group_pa_table[,1]
  group_pa_table <- group_pa_table[,-1]
  
  
  normal_pa <- subset(group_pa_table, group_pa_table$Group =="normal")
  intermediate_pa <- subset(group_pa_table, group_pa_table$Group=="intermediate")
  small_pa <- subset(group_pa_table, group_pa_table$Group=="small intestine")
  
    normal_pa$Group <- NULL
    intermediate_pa$Group <- NULL
    small_pa$Group <- NULL
  
    sum_table <- matrix(ncol = ncol(normal_pa), nrow=3)
    
  for (i in 1:ncol(intermediate_pa)){
    
    a = sum(intermediate_pa[,i])
    
    sum_table[2,i] <- a
    
  }  
    
    rownames(sum_table) <- c("normal", "intermediate", "small_intestine")
    colnames(sum_table) <- colnames(normal_pa)
    
    
    
    