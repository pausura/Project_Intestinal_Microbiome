##### Merge table script #####

getwd()

setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

# Open data
##(No consider the header - for next merge)
new_ID <- read.table("~/Documents/Universitat/Holanda/Projecte/Dades/new_IBD_samples.txt", sep = "\t", row.names = 1)
##(Consider the header - for last merge)
old_ID <- read.table("~/Documents/Universitat/Holanda/Projecte/Dades/old_IBD_samples.txt", sep = "\t", header = TRUE)

# Change the name of the first row (for next merge)
rownames(new_ID) [1] <- "WGSProjectID"
colnames(old_ID) [1] <- "Taxa"

g_number <- read.table("~/Documents/Universitat/Holanda/Projecte/Dades/g_number.txt", sep = "\t", header = T)

# Transpose table new_ID
t_new_ID <- t(new_ID)
t_newID <- as.data.frame(t_new_ID)

## Merge t_new_ID table with g_number table by G_number
## all = TRUE condition -> keep all the columns (include this ones that don't match) - NA values
merged_table <- merge (t_newID, g_number, by="WGSProjectID")

## Change NA values for 0 - is.na function
merged_table[is.na(merged_table)] <- 0


##Export merged_table to see if it's fine merged
write.table (merged_table,file = "~/merged_table.txt", quote = F, sep = "\t")

## Delete repeated columns and re-organize the merged_table in excel (saved as merged_table2.txt)

## Transpose again the new merged table
merged_table2 <- read.table("~/Documents/Universitat/Holanda/Projecte/merged_table2.txt", sep = "\t", row.names = 1)

t_merged_table <- t(merged_table2)
## Change column1 name
colnames(t_merged_table)[1] <- "Taxa"

#Convert tables as dataframes
df_new = as.data.frame(t_merged_table)
df_old = as.data.frame(old_ID)


## Merge old_IBD + new_IBD
total_table <- merge (df_new, df_old, by = "Taxa", all = TRUE)
  total_table[is.na(total_table)] <- 0

  ##Export merged_table to see if it's fine merged
write.table (total_table,file = "~/total_table.txt", quote = F, sep = "\t")
