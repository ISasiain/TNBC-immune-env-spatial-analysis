rm(list = ls())
setwd("/home/isc/Spatial_immune_env")


#
# LOADING THE REQUIRED PACKAGES
#

if (!require("SpatialExperiment", quietly = TRUE))
  BiocManager::install("SpatialExperiment")
library(SpatialExperiment)

if (!require("SPIAT", quietly = TRUE))
  BiocManager::install("SPIAT")
library(SPIAT)

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!require("stats", quietly = TRUE))
  install.packages("stats")
library(stats)

if (!require("fpc", quietly = TRUE))
  install.packages("fpc")
library(fpc)

if (!require("NbClust", quietly = TRUE))
  install.packages("NbClust")
library(NbClust)

if (!require("alphahull", quietly = TRUE))
  install.packages("alphahull")
library(alphahull)

#
# READING THE SPE OBJECT TO ANALYSE
#


#Reading spe
spe_path <- "./data_from_suze/data/spiat/PD31028a_BLOCK_3|L4_spiat.rds"
my_spe <- readRDS(spe_path)



#Generating matrix to store the values

# Spatial coors + densities
number_of_columns <- 2 + length(unique(my_spe$"Cell.Type"))

#Cell numbers
number_of_rows <- length(my_spe$"Cell.ID")

#Generating empty matrix and assigning column and row names
values_for_clustering <- matrix(NA, nrow=number_of_rows, ncol=number_of_columns)

rownames(values_for_clustering) <- my_spe$"Cell.ID"
colnames(values_for_clustering) <- c("X_coor", "Y_coor", unique(my_spe$"Cell.Type"))

#Generating data frame with all the cellIDs, Cell.Type and coordinates
spatial_data <- cbind(my_spe@colData[, c("Cell.ID", "Cell.Type")], spatialCoords(my_spe))


#
# FILLING THE MATRIX
#

#Getting spatial coordinates
values_for_clustering[,c("X_coor",  "Y_coor")] <- spatialCoords(my_spe)

#Getting DIN values

for (cell_type in unique(my_spe$"Cell.Type")) {
  
  #Use sapply to generate a vector with the DINs of each cell type
  values_for_clustering[,cell_type] <- sapply(
    
    #Calculating DIN for each cell in the sample
    rownames(values_for_clustering), 
    
    function(cell) {
      
      #Defining parameters
      radius <- 100
      X_coor <- values_for_clustering[cell, "X_coor"]
      Y_coor <- values_for_clustering[cell, "Y_coor"]
      
      #Creating or reassigning a column in spatial_data dataframe to calculate
      #the euclidean distance between cells
      spatial_data$"Distance" <-
        (spatial_data$Cell.X.Position - X_coor)^2 + 
        (spatial_data$Cell.Y.Position - Y_coor)^2
      
      #Filtering the cells not included in the area of interest and that do not
      #belong to the category of interest. Counting the cells that match the conditions.
      return(nrow(spatial_data[spatial_data$"Distance" <= radius^2 & spatial_data$"Cell.Type" == cell_type,])/(pi*radius^2))
      
    }
  )
}



#
# DETERMINING CELL CLUSTERS BASED ON THE VADENSITY VALUES
#

# Using K means clustering

#NbClust only accepts a minimum number of two clusters
optimal_clusters <- NbClust(values_for_clustering[,-c(1,2)], method = "kmeans")
optimal_clusters2 <- kmeans(values_for_clustering[,-c(1,2)], centers = 3)
mean(as.vector(dist(optimal_clusters2$centers)))


#If the number of clusters is set to be 2, check the number of cluster

#Determining optimal number of clusters
#optimal_clusters <- clusGap(values_for_clustering[,-c(1,2)],
#                            FUN = kmeans,
#                            K.max = 15)

#fviz_gap_stat(optimal_clusters)

#get_number_of_clusters <- function(clusGap_results) {
  
#  for (row in 1:(nrow(clusGap_results$Tab)-1)) {
    
#    if (clusGap_results$Tab[row,"gap"] >= (clusGap_results$Tab[row+1,"gap"]-clusGap_results$Tab[row+1,"SE.sim"])) {
    
#      return(row)
#      break
#    }}
    
#    print("No optimal number of clusters has been  detected. 1 cluster is assumed.")
#    return(1)
  
#}

#Determine number of clusters
#cl_num <- get_number_of_clusters(clusGap_results = optimal_clusters)

#Generate clusters
#clustering <- kmeans(values_for_clustering[,-c(1,2)], centers = cl_num)

#Plotting identified clusters
clusters_plot <- ggplot(
  as.data.frame(values_for_clustering), 
  aes(x = X_coor, y = Y_coor, colour=as.factor(optimal_clusters2$cluster))) +
  geom_point() +
  labs(color = "Clusters") +
  theme_classic()

clusters_plot
plot_cell_categories(my_spe, categories_of_interest = c("CD3", "p53", "CD20", "CD8"), colour=c("red", "grey", "green", "blue"))


par(mfrow = c(2, 1))
boxplot(as.data.frame(values_for_clustering[optimal_clusters2$cluster==1, c(3,4,5,6)]), ylim=c(0,0.0014))
boxplot(as.data.frame(values_for_clustering[optimal_clusters2$cluster==2, c(3,4,5,6)]), ylim=c(0,0.0014))
boxplot(as.data.frame(values_for_clustering[optimal_clusters2$cluster==3, c(3,4,5,6)]), ylim=c(0,0.0014))

#Generate stacked box plot of cell counts of each cluster
clusters <- sort(unique(unname(optimal_clusters$Best.partition)))
cells <- unique(my_spe$Cell.Type)

df_to_plot <- data.frame()

for (cl in clusters) {
  for (cell in cells) {
    
    new_row <- c(cl, cell, table(my_spe$Phenotype[optimal_clusters$Best.partition==cl])[cell])
    df_to_plot <- rbind(df_to_plot, new_row)
    
  }
}

colnames(df_to_plot) <- c("Clusters", "Cells", "Counts")
df_to_plot$Clusters <- as.factor(df_to_plot$Clusters)
df_to_plot$Counts <- as.numeric(df_to_plot$Counts)

#Stacked barplot of counts
ggplot(df_to_plot, aes(fill=Cells, y=Counts, x=Clusters)) + 
  geom_bar(position="stack", stat="identity")

#Stacked barplot of proportions
ggplot(df_to_plot, aes(fill=Cells, y=Counts, x=Clusters)) + 
  geom_bar(position="fill", stat="identity") +
  ylab("Cell proportions")

#Calculating attraction of each cell type.

#Defining general cell denisty
cell_density <- as.list(table(my_spe$Cell.Type) / (pi * 1500^2))


#Calculating attraction
attraction_values <- values_for_clustering

attraction_values[,"CD20"] <- attraction_values[,"CD20"] / cell_density$CD20
attraction_values[,"CD3"] <- attraction_values[,"CD3"] / cell_density$CD3
attraction_values[,"CD8"] <- attraction_values[,"CD8"] / cell_density$CD8
attraction_values[,"p53"] <- attraction_values[,"p53"] / cell_density$p53


#Define colours for plotting
breaks_values <- c(0,1,7)
colours <- c("blue", "grey", "red")

#Plotting attraction values of each cell type
p53_attraction_plot <- ggplot(
  as.data.frame(attraction_values), 
  aes(x = X_coor, y = Y_coor, colour=attraction_values[,"p53"])) +
  geom_point() +
  labs(color = "Attraction value") +
  theme_classic() +
  scale_color_gradient(low = "grey", high = "red", limits = c(0,25))

CD3_attraction_plot <- ggplot(
  as.data.frame(attraction_values), 
  aes(x = X_coor, y = Y_coor, colour=attraction_values[,"CD3"])) +
  geom_point() +
  labs(color = "Attraction value") +
  theme_classic() +
  scale_color_gradient(low = "grey", high = "red", limits = c(0,25))

CD8_attraction_plot <- ggplot(
  as.data.frame(attraction_values), 
  aes(x = X_coor, y = Y_coor, colour=attraction_values[,"CD8"])) +
  geom_point() +
  labs(color = "Attraction value") +
  theme_classic()  +
  scale_color_gradient(low = "grey", high = "red", limits = c(0,25))

CD20_attraction_plot <- ggplot(
  as.data.frame(attraction_values), 
  aes(x = X_coor, y = Y_coor, colour=attraction_values[,"CD20"])) +
  geom_point() +
  labs(color = "Attraction value") +
  theme_classic() +
  scale_color_gradient(low = "grey", high = "red", limits = c(0,25))

CD3_attraction_plot
CD8_attraction_plot
CD20_attraction_plot
p53_attraction_plot

#Generating boxplots for the atraction values of each detected cluster
boxplot(as.data.frame(attraction_values[optimal_clusters$Best.partition==1, c(3,4,5,6)]), ylim=c(0,29))
boxplot(as.data.frame(attraction_values[optimal_clusters$Best.partition==2, c(3,4,5,6)]), ylim=c(0,29))
boxplot(as.data.frame(attraction_values[optimal_clusters$Best.partition==3, c(3,4,5,6)]), ylim=c(0,29))


a <- ashape(values_for_clustering[optimal_clusters$Best.partition==1,c(1, 2)], alpha=40)
b <- ashape(values_for_clustering[optimal_clusters$Best.partition==2,c(1, 2)], alpha=75)
c <- ashape(values_for_clustering[optimal_clusters$Best.partition==3,c(1, 2)], alpha=75)

plot(values_for_clustering[, c(1,2)], pch=20)
plot(a, do.shape=TRUE, wpoints=FALSE)

test <- dbscan(values_for_clustering[optimal_clusters$Best.partition==1,c(1, 2)], MinPts = 25, eps = 125)



ggplot(
  as.data.frame(values_for_clustering[optimal_clusters$Best.partition==1,c(1, 2)]), 
  aes(x = X_coor, y = Y_coor, colour=as.factor(test$cluster))) +
  geom_point() +
  labs(color = "Clusters") +
  theme_classic()



#Test. Trying to estimate areas of clusters based on areas arond points

my_cells <- values_for_clustering[optimal_clusters$Best.partition==1,c(1, 2)]




X <- sample(0:100, 40, replace = FALSE)
Y <- sample(0:100, 40, replace = FALSE)

AH <-ahull(cbind(X, Y), alpha=100)

plot(cbind(X, Y)[c(31,35,6,37,39,14,5,13),])

par(mfrow = c(1, 3))
plot(a, col = c(4), xlim=c(0,3000), ylim=c(0,3000))
plot(b, col = c(1), xlim=c(0,3000), ylim=c(0,3000))
plot(c, col = c(2), xlim=c(0,3000), ylim=c(0,3000))

plot(a)

library(factoextra)
fviz_nbclust(values_for_clustering[,-c(1,2)], FUNcluster = kmeans)


AH$ashape.obj
