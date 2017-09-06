
setwd("~/Documents/Universitat/Holanda/logistic_regression_adonis_data/normal_small_results/filtered_normal_small/")

list.files(pattern=".txt")

list.filenames<-list.files(pattern=".txt")
list.filenames

list.data<-list()

for (i in 1:length(list.filenames)){
  list.data[[i]]<-read.csv(list.filenames[i])
}

for (ii in 1:length(list.data))

  new_matrix <- as.data.frame(list.data[[ii]])