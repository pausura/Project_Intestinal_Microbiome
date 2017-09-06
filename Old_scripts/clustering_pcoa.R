
 ### CLUSTERS ###

getwd()
setwd("/Users/paulasuredahorrach/documents/universitat/holanda/projecte")

  ## Open species table

species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", header = T, row.names = 1, check.names = F)
        
        ## Transpose species table
           t_species <- as.data.frame(t(species_table))

library(vegan)

      ## Calculate distance matrix - Bray method
           beta <- vegdist(t_species, method="bray")
          
   # Average distance 
   
   library(gclus)   
   caver <- hclust(beta, method="aver")
   caver1 <- reorder.hclust(caver, beta)
     
   library(data.table) 
      order_caver1= as.data.frame(caver1$labels[c(caver1$order)])  
      setDT(order_caver1, keep.rownames = TRUE)[]
          colnames(order_caver1)[1] <- "num"
          colnames(order_caver1)[2] <- "IBDFEC"
      
   intestinal_cat <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T)
      colnames(intestinal_cat)[1] <- "IBDFEC"
   
   order_cat <- merge(intestinal_cat, order_caver1, by="IBDFEC")
      rownames(order_cat) <- order_cat[,1]
      order_cat <- order_cat[,-1]
      
      order_cat$color ="none"
      order_cat[order_cat$Group=="normal",]$color = "#B0B2B3"
      order_cat[order_cat$Group=="small intestine",]$color = "red"
      order_cat[order_cat$Group=="intermediate",]$color = "#2F2BFF"
   
      order_cat$values <- "5"
      order_cat$values <- as.numeric(as.character(order_cat$values))
      
    ## Sort by "num" column    
  order_cat1 <- order_cat[order(as.numeric(order_cat$num)),] 
  
  par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
      plot(caver1, hang=-1, labels = FALSE, axes = FALSE, ylab = "", xlab="", sub="") 
      barplot(order_cat1$values, col=order_cat1$color, border = NA, yaxt="n")
  
  

      ## Linkage clusters -- 3 clusters assigned       
      # Single distance
      csin <- hclust(beta, method="single")
      plot(csin, hang=-1)
      rect.hclust(csin, 3)
      # Complete distance
      ccom <- hclust(beta, method="complete")
      plot(ccom, hang=-1)
      rect.hclust(ccom, 3)
      # Average distance 
      caver <- hclust(beta, method="aver")
      plot(caver, hang=-1, labels = FALSE)
      rect.hclust(caver, 3)
      
      cl <- cutree(caver, 3)  
      
      ## Clustering and ordination - PCoA
      
      genus_pcoa <- as.data.frame(cmdscale(beta, k = 2))
      ordiplot(genus_pcoa)
      
      # Optimized clustering
      ckm <- kmeans(decostand(t_genus, "hell"), 3)
      ckm$cluster
      
      ordiplot(genus_pcoa, dis="si")
      ordihull(genus_pcoa, ckm$cluster, col="red")
      
      genus_calinski <- cascadeKM(scale(genus_pcoa, center = TRUE,  scale = TRUE), 1, 5, iter = 1000)
      plot(genus_calinski, sortg = TRUE, grpmts.plot = TRUE)
      calinski_best <- as.numeric(which.max(genus_calinski$results[2,]))
      
      
      ## Mclust optimization
      library(mclust)
      
      genus_pcoa_clust <- Mclust(as.matrix(genus_pcoa), G=1:20)
      m_best <- dim(genus_pcoa_clust$z)[2]
      cat("Model-based optimal number of clusters:", m_best, "\n")
      
      plot(genus_pcoa_clust)
      
      
      ## pamk optimization
      library(fpc)
      library(cluster)
      
      pamk_best <- pamk(genus_pcoa)
      cat("Number of clusters estimated by optimum average silhouette width:", pamk_best$nc, "\n")
      
      plot(pam(genus_pcoa, pamk_best$nc))  
      
      
      ## Affinity propagation - Heat map
      library(apcluster)
      
      genus_pcoa1 <- genus_pcoa
      row.names(genus_pcoa1) <- NULL
      
      d_apclus <- apcluster(negDistMat(r=2), genus_pcoa1)
      cat("Affinity propogation optimal number of clusters:", length(d_apclus@clusters), "\n")
      
      heatmap(d_apclus)
      plot(d_apclus, genus_pcoa)
      
      
      ### Generate HeatMap
      library(gplots)
      source("http://bioconductor.org/biocLite.R")
      biocLite("Heatplus") 
      library(Heatplus)
      library(vegan)
      library(RColorBrewer)
      
      ## Create a matrix with 541 samples (and not 544)
      
      intestinal_groups <- read.table("~/Documents/Universitat/Holanda/Projecte/intestinal_groups.txt", sep = "\t", header = T, row.names = 1)
      
      new_species <- merge(intestinal_groups, t_species, by = "row.names")
      rownames(new_species) <- new_species[,1]
      new_species <- new_species[,-1]
      
      ## Remove group column
      new_species$Group <- NULL
      
      ## Calculate the new distance matrix
      new_beta <- vegdist(new_species, method="bray")
      new_pcoa <- as.data.frame(cmdscale(new_beta, k = 2))
      
      library(apcluster)
      
      d_apclus_new <- apcluster(negDistMat(r=2), new_pcoa)
      
      row_clus <- hclust(new_beta, "aver")
      col_clus <- hclust(new_beta, "aver")
      
      hc_dist= dist(new_genus)
      hc_clust= hclust(hc_dist)
      hr_dist= dist(t(new_genus))
      hr_clust= hclust(hr_dist)
      
      ## Generate a vector for the colors depending on the intestinal category     
      intestinal_groups$color ="none"
      intestinal_groups[intestinal_groups$Group=="normal",]$color = "#1781CE"
      intestinal_groups[intestinal_groups$Group=="small intestine",]$color = "#17CE63"
      intestinal_groups[intestinal_groups$Group=="intermediate",]$color = "#A41672"
      
      order_def <- order_cat[order(order_cat$num),]
      rownames(order_def) <- order_def[,2]
      order_def <- order_def[,-2]
      
      heatmap1 <- heatmap.2(as.matrix(d_apclus_new@sim), Rowv = as.dendrogram(hc_clust), Colv = as.dendrogram(hc_clust),  trace = "none", density.info = "none", lhei = c(2, 8), labRow = FALSE, labCol = FALSE, RowSideColors = order_def[,3]) 
      
  
  
  