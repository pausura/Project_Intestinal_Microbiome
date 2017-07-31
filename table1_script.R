##### Table 1 script #####

##Indicates working directory
getwd()

##Changes workind directory
setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

##Open data
total_metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep = "\t", header = T, row.names = 1)
       
#Categorical values - create function
      tblFun <- function(x) {
          #Create a table
        tbl <- table(x)
          #Combine columnes/rows to get the counts and percentage (creates new table -> res)
        res <- cbind(tbl,round(prop.table(tbl)*100,2))
          #Give names to the columns
        colnames(res) <- c('Count','Percentage')
        res
      }

  # NA sum function
      nzsum <- function(x){
        sum (is.na(x))
      }
      
#Continuous variables - two packages
      library (psych)
        #describe()
      library (pastecs)
        #stat.desc()
        
##Divide data in groups
   #CD
      table_cd <- subset(total_metadata,total_metadata$DiagnosisCurrent=="CD")
   #UC
      table_uc <- subset(total_metadata,total_metadata$DiagnosisCurrent=="UC")        
   
   #Normal Group
      table_normal <- subset(total_metadata,total_metadata$Group=="normal")
   #Intermediate Group
      table_intermediate <- subset(total_metadata, total_metadata$Group=="intermediate")
   #Small Intestine Group
      table_smallintestine <- subset(total_metadata, total_metadata$Group=="small intestine")
   

### Loop for ###

#Create an empty table: nrows = ncolumns original table
my_results=matrix(ncol = 6, nrow=ncol(table_smallintestine))

for (i in 1:ncol(table_smallintestine)) { # For each loop goes to the next column (numerical way)
  #print (my_table[,i]) # Right loop verification
  
# Condition: if the column values are numerical/continuous
  if ( is.numeric(table_smallintestine[,i])){
    #print (paste0("Column numeric: ", i)) # Right loop verification
    
  # Keep in "x" the result from describe function (done in the columns) - for each factor
    x=describe(table_smallintestine[,i])
    z = nzsum(table_smallintestine[,i])
  # In the new table ("x"): keep different values in the different columns (in crhonological order)
    my_results[i,1]="numerical"
    my_results[i,2]=x$median
    my_results[i,3]=x$mean
    my_results[i,4]=x$sd
    my_results[i,5]=x$n
    my_results[i,6]=z
    
  } 
  
  # Condition: if the column values are categorical
  else {
    #print (paste0("Column categorical ", i)) # Right loop verification
    
  # Keep in "x" the result from tblFun function (done in the columns) - for each factor
    x=tblFun(table_smallintestine[,i])
    z = nzsum(table_smallintestine[,i])
  
  # In the new table ("x"): keep different values in the different columns (in crhonologial order)
    my_results[i,1]="categorical"
      # toString to keep the possible different values/categories in the same vector/column
    my_results[i,2]=toString(rownames(x))
      # First column table x = 'Count'
    my_results[i,3]=toString(x[,1]) 
      # Second column table x = 'Percentage'
    my_results[i,4]=toString(x[,2])
      # Sum of the values on column1 ("x")
    my_results[i,5]=sum(x[,1])
    my_results[i,6]= z
    
  }
}

# The column names from the original table = row names from the new table
rownames(my_results) = colnames(table_smallintestine)

# Give names to the columns of the new table
colnames(my_results) = c("Type", "Categories/Median", "Counts / Mean", "SD / %", "Number of non-zeros (n)", "Number of NA") 


## Export the new table
write.table (my_results, file = "~/small_intestine_summary_table.txt", quote = F, sep = "\t")







    



