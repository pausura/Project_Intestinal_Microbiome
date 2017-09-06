

species_table <- read.table("~/Documents/Universitat/Holanda/Projecte/Filtered_DUDes/species_table_DUDes.txt", sep = "\t", header = T, row.names = 1, check.names = F)

library(vegan)
library(ggplot2)

# Open dates table
dates_table <- read.table("~/Documents/Universitat/Holanda/Projecte/final_date_table.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Generate a distance matrix - Bray method
beta <- vegdist(species_table, method="bray")
# cmdscale -> multidimensional scaling of a data matrix (distance matrix) 
my_pcoa <- as.data.frame(cmdscale(beta, k = 4))
#colnames(my_pcoa)[1:4] <- c("PCoA1","PCoA2","PCoA3","PCoA4")

## Create a full table to color the PCoA plot (by phyla, reads, shannon index and Group)
plot_table1 <- merge(dates_table, my_pcoa, by="row.names")
    rownames(plot_table1) <- plot_table1[,1]
    plot_table1 <- plot_table1[,-1]

  
plot_table1$recent_scale = "none"
    plot_table1[plot_table1$RecentDifference == 0, ]$recent_scale = "#B0A7A7"
    plot_table1[plot_table1$RecentDifference > 1 & plot_table1$RecentDifference < 500 , ]$recent_scale = "#0000FF"
    plot_table1[plot_table1$RecentDifference >= 500 & plot_table1$RecentDifference < 1000 , ]$recent_scale = "#62A4D1"
    plot_table1[plot_table1$RecentDifference >= 1000 & plot_table1$RecentDifference < 3000 , ]$recent_scale = "#5BE55B"
    plot_table1[plot_table1$RecentDifference >= 3000, ]$recent_scale = "#FFF000"
       
recent_dates_plot <- ggplot (plot_table1, aes(x=V1, y=V2, geom="blank", colour=recent_scale)) + geom_point () + scale_color_identity("Recent_dates", breaks = plot_table1$recent_scale, labels = plot_table1$RecentDifference, guide = "legend") +  theme_classic() + labs(x="PCoA1", y="PCoA2")
recent_dates_plot

plot_table1$farther_scale = "none"
    plot_table1[plot_table1$FartherDifference >= 0 & plot_table1$FartherDifference < 500 , ]$farther_scale = "#20aaf9"
    plot_table1[plot_table1$FartherDifference >= 500 & plot_table1$FartherDifference < 1000 , ]$farther_scale = "#0000FF"
    plot_table1[plot_table1$FartherDifference >= 1000 & plot_table1$FartherDifference < 3000 , ]$farther_scale = "#62A4D1"
    plot_table1[plot_table1$FartherDifference >= 3000 & plot_table1$FartherDifference < 5000 , ]$farther_scale = "#5BE55B"
    plot_table1[plot_table1$FartherDifference >= 5000 & plot_table1$FartherDifference < 8000 , ]$farther_scale = "#FFF000"
    plot_table1[plot_table1$FartherDifference >= 8000, ]$farther_scale = "#FF0000"
    

farther_dates_plot <- ggplot (plot_table1, aes(x=V1, y=V2, geom="blank", colour=farther_scale)) + geom_point () + scale_color_identity("Farther_dates", breaks = plot_table1$farther_scale, labels = plot_table1$FartherDifference, guide = "legend") +  theme_classic() + labs(x="PCoA1", y="PCoA2")
farther_dates_plot




