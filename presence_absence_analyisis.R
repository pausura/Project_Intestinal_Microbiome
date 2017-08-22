
metadata_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
  sex <- subset(metadata_table, select= Sex)
  ppi <- subset(metadata_table, select = PPI)
  bmi <- subset(metadata_table, select = BMI)
  reads <- subset(metadata_table, select = TotalReads)
  antibiotics <- subset(metadata_table, select = MedicationAntibiotics)
  age <- subset(metadata_table, select = AgeAtFecalSampling)
  
  
metadata_analysis <- merge(sex, ppi, by = "row.names")
  rownames(metadata_analysis) <- metadata_analysis[,1]
  metadata_analysis <- metadata_analysis[,-1]
metadata_analysis <- merge(metadata_analysis, bmi, by = "row.names")
  rownames(metadata_analysis) <- metadata_analysis[,1]
  metadata_analysis <- metadata_analysis[,-1]
metadata_analysis <- merge(metadata_analysis, reads, by = "row.names")
  rownames(metadata_analysis) <- metadata_analysis[,1]
  metadata_analysis <- metadata_analysis[,-1]
metadata_analysis <- merge(metadata_analysis, antibiotics, by = "row.names")
  rownames(metadata_analysis) <- metadata_analysis[,1]
  metadata_analysis <- metadata_analysis[,-1]
metadata_analysis <- merge(metadata_analysis, age, by = "row.names")
  rownames(metadata_analysis) <- metadata_analysis[,1]
  metadata_analysis <- metadata_analysis[,-1]
  
  
groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1, check.names = F)
  groups$name <- rownames(groups)
tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_taxonomy.txt", sep = "\t", header = T, row.names = 1, check.names = F)

intermediate_normal <- shannon_groups_table[!(shannon_groups_table$Group=="small intestine"),]

normal_intermediate <- groups[!(groups$Group == "small intestine"),]
  normal_intermediate$name <- NULL
normal_small <- groups[!(groups$Group == "intermediate"),]
  normal_small$name <- NULL
intermediate_small <- groups[!(groups$Group == "normal"),]
  intermediate_small$name <- NULL

p_a_table <- tax_table
for (i in 1:ncol(tax_table)) {
  
  for (j in 1:nrow(tax_table)) {
    
    if (tax_table[j,i]>0) {
      
      p_a_table[j,i] = 1
      
    }
    
  }
}

archaea <- subset(p_a_table, select = k__Archaea)

###Merge metadata_analysis | precence_archaea | intestinal_group

##############
#Normal-Small#
##############

normal_small_data <- merge(normal_small, metadata_analysis, by = "row.names")
  rownames(normal_small_data) <- normal_small_data[,1]
  normal_small_data <- normal_small_data[,-1]
normal_small_data <- merge(normal_small_data, archaea, by = "row.names")  
  rownames(normal_small_data) <- normal_small_data[,1]
  normal_small_data <- normal_small_data[,-1]

normal_small_data1 <- normal_small_data  

## Convert integer variables to numeric variables

for (i in 1:ncol(normal_small_data)) {
  
  if (is.numeric(normal_small_data[,i])) {
    
    normal_small_data1[,i] <- as.numeric(normal_small_data[,i])
    
  }
  
  else {
    
    normal_small_data1[,i] <- as.numeric(normal_small_data[,i])
    
  }
  
  
}

normal_small_data2 <- normal_small_data1
### Change 2 and 3 values in Group column
for (j in 1:nrow(normal_small_data2)) {
  
  if (normal_small_data2[j,1] = 2) {
    
    normal_small_data2[j,1] = 1
    
  }

} 
  


##Examining correlations among variables
library(PerformanceAnalytics)
chart.Correlation(normal_small_data1, method = "spearman", histogram = T, pch = 16)

library(psych)
corr.test(normal_small_data1, use = "pairwise", method = "spearman", adjust = "fdr", alpha = 0.05)



### Determining model with step procedure
ns_data_omit <- na.omit(normal_small_data)

# Define full and null models and do step procedure

  # Null = no independent variables
ns_model_null <- glm(k__Archaea ~ 1, data = ns_data_omit, family = binomial(link = "logit"))
  # Full = independent variables
ns_model_full <- glm(k__Archaea ~ PPI + Sex + AgeAtFecalSampling + TotalReads + BMI + MedicationAntibiotics + Group, data = ns_data_omit, family = binomial(link = "logit"))

step(ns_model_null, scope = list(upper = ns_model_full), direction = "both", test = "Chisq", data = normal_small_data1)

# Final model (with the independent variables with significant Pr in last step procedure)
ns_model_final <- glm(k__Archaea ~ AgeAtFecalSampling + Group + MedicationAntibiotics, data = normal_small_data, family = binomial(link = "logit"), na.action(na.omit))
summary(ns_model_final)

Anova(ns_model_final, type = "II", test = "Wald")

### Overall p-value for model

# Create data frame with variables in final model and NA’s omitted
library(dplyr)
ns_data_final <- select(normal_small_data1, k__Archaea, AgeAtFecalSampling, Group, MedicationAntibiotics)
ns_data_final <- na.omit(ns_data_final)
### Define null models and compare to final model
ns_final_model_null <- glm(k__Archaea ~1, data = ns_data_final, family = binomial(link="logit"))
anova(ns_model_final, ns_final_model_null, test = "Chisq" )

library(lmtest)
lrtest(ns_model_final)

plot(fitted(ns_model_final), rstandard(ns_model_final))


# Simple plot of predicted values
library(dyplr)

ns_data_final$predy = predict(ns_model_final, type = "response")
plot(k__Archaea ~ predy, data = ns_data_final, pch = 16, xlab = "Predicted probability of 1 response", ylab = "Actual response")

