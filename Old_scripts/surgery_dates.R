
## Surgery dates ##

  #Open date table
date_table <- read.table("~/Documents/Universitat/Holanda/Projecte/time_differences.txt", sep = "\t", header = T, row.names = 1, check.names = F)
date_table1 <- date_table

  #Remove all samples that don't have resection and pouch/stoma dates
  date_table1 <- date_table[!(is.na(date_table$`DayDifferencePouch/Stoma`) & is.na(date_table$DayDifferenceResection)) ,]
      
  #Remove unnecessary columns
 date_table1$FecalsampleProductionDate <- NULL
 date_table1$LastPouchStoma <- NULL
 date_table1$LastResection <- NULL
  
 #Add columns for the most recent and the farther difference in days
      date_table1$RecentDifference <- "none"
      date_table1$FartherDifference <- "none"

 
  for (i in 1:nrow(date_table1)) {
    
      # If PouchStoma is NA
    if(is.na(date_table1$`DayDifferencePouch/Stoma`)[i]) {
          # Further days = Resection days
        date_table1$FartherDifference[i]=date_table1$DayDifferenceResection[i]
          # Replace NA value
        date_table1$`DayDifferencePouch/Stoma`[i] <- 0
        
    }
      # If Resection is NA
    if(is.na(date_table1$DayDifferenceResection)[i]) {
        # Further days = PouchStoma days
      date_table1$FartherDifference[i] = date_table1$`DayDifferencePouch/Stoma`[i]
        # Replace NA value
      date_table1$DayDifferenceResection[i] <- 0
      
    }
    
      # If Pouch/Stoma days > Resection days
    if (date_table1$`DayDifferencePouch/Stoma`[i] > date_table1$DayDifferenceResection[i]) {
        # Further days = PouchStoma days
      date_table1$FartherDifference[i] = date_table1$`DayDifferencePouch/Stoma`[i]
        # Recent days = Resection days
      date_table1$RecentDifference[i] = date_table1$DayDifferenceResection[i]
    }
    
    else { #else: PouchStoma days < Resection days
        # Further days = Resection days
     date_table1$FartherDifference[i]=date_table1$DayDifferenceResection[i]
        # Recent days = Pouch/Stoma days
     date_table1$RecentDifference[i]=date_table1$`DayDifferencePouch/Stoma`[i]
    }
    
  }
 
      #Replace 0 values for NA in Recent and Further columns
  #final_date_table <- date_table1
  #final_date_table[final_date_table == 0] <- NA  
      
  write.table(date_table1, file = "~/final_date_table.txt", sep = "\t", quote = F)
  
  
  
  
  
      