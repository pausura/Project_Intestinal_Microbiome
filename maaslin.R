
    ##MaAsLin files

#Open taxa data
  tax_table <- read.table("~/Documents/Universitat/Holanda/Projecte/filtered_tax_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Open metadata
  metadata <- read.table("~/Documents/Universitat/Holanda/Projecte/Metadata/filtered_metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
  
maaslin_data <- merge(metadata, tax_table, by = "row.names")
  rownames(maaslin_data) <- maaslin_data[,1]
  maaslin_data <- maaslin_data[,-1]

  write.table(maaslin_data, file = "~/maaslin_data.tsv", sep ="\t", quote = F)

  
##Install older version of gam package
  require(devtools)
  install_version("gam", version = "1.14", repos = "http://cran.us.r-project.org")
  
  
  library(Maaslin)
  
  #Example Maaslin
  example(Maaslin)
      InputTSV <- system.file('extdata', 'maaslin_demo2.tsv', package = "Maaslin")
      InputConfig <- system.file('extdata', 'maaslin_demo2.read.config', package = "Maaslin")
  Maaslin(InputTSV, "maaslin_example_output", strInputConfig = InputConfig)

  
  ## MY DATA
  
  Maaslin("~/Documents/Universitat/Holanda/Projecte/maaslin_data.tsv", "maaslin_data_output", strInputConfig = "~/Documents/Universitat/Holanda/Projecte/maaslin_data.read.config")
  
  
  
  